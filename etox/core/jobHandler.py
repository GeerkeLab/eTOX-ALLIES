# -*- coding: utf-8 -*-

import os
import tempfile
import json
import time
import logging
import operator
import glob
import signal
import shutil
from psutil import Process

from .. import settings


def killProcs(PID):
    #try:
    p=Process(PID)
    #try:
    while len(p.children())>0:
        try:
            child=p.children()[0]
            while len(child.children())>0:
                print len(child.children()), child
                child=child.children()[0]
            print "kill"
            child.terminate()
            print "Restart"
        except:
            pass
    try:
        p.terminate()
    except:
        pass


class WarmKiller:
  kill_now = False
  def __init__(self):
    logging.info('SIGNAL NOT ALLOWED IN THREAD')
    #signal.signal(signal.SIGINT, self.exitProg)
    #signal.signal(signal.SIGTERM, self.exitProg)

  def exitProg(self,signum, frame):
    self.kill_now = True


class jobCpd:
    ### jobs for cpds have 'job' extension
    ### jobs for sdf have 'dsr' extension (DataSetRun)
    ### jobs for cpd has prediction=None
    ### jobs for sdf has prediction=True or False
    def __init__(self):
        self.filename=''
        self.dirTemp=None
        self.sdf=''
        self.modelProt=''
        self.modelProtVer=''
        self.dateSub=''
        self.jobid='<not defined>'  # Job ID assigned by etoxwsapi
        self.prediction=None
        self.experiment=None
        self.status=None
        self.results={}

    def create(self,sdf,modelProt,modelProtVer,prediction=None,jobid='<not defined>'):
        
        if prediction is None:
            extFile='job'
        else:
            extFile='dsr'
        
        newJob=tempfile.NamedTemporaryFile(
            prefix=settings.get('etoxlie_prefixFn', 'iLIEtmp'),
            suffix='.%s'%extFile,
            dir=settings.get('etoxlie_folder', ''),
            delete=False
        )      
        self.filename=newJob.name
        
        if prediction is None:
            self.dirTemp=os.path.splitext(newJob.name)[0]
        else:
            self.dirTemp=None

        self.sdf=sdf
        self.modelProt=modelProt
        self.modelProtVer=modelProtVer
        self.prediction=prediction
        self.jobid=jobid
        self.experiment=None
        self.dateSub=time.time()
        self.status="SUBMITTED"
        self.results=None
        
        infoJob={"FileName":self.filename,
                 "TempDir":self.dirTemp,
                 "Model":self.modelProt,
                 "Version":self.modelProtVer,
                 "Prediction":self.prediction,
                 "Experiment":self.experiment,
                 "Date":time.strftime("%Y-%m-%d %H:%M:%S",time.localtime(self.dateSub)),
                 "EpochDate":self.dateSub,
                 "ApiJobId":self.jobid,
                 "Status":self.status,
                 "Results":self.results,
                 "SDF":self.sdf
                 }
        
        json.dump(infoJob,newJob,indent=4)
        newJob.close
        
        ## for compatibility with apache user in etoxsys
        os.chmod(self.filename, 0775)
        logging.debug("Created new job %s"%self.filename)

    def load(self,filename):
        with open(filename,'r') as infile:
            infoJob=json.load(infile)
    
        self.filename=infoJob['FileName']
        self.dirTemp=infoJob['TempDir']
        self.sdf=infoJob['SDF']
        self.modelProt=infoJob['Model']
        self.modelProtVer=infoJob['Version']
        self.prediction=infoJob['Prediction']
        self.experiment=infoJob['Experiment']
        self.dateSub=infoJob["EpochDate"]
        self.jobid=infoJob["ApiJobId"]
        self.status=infoJob["Status"]
        self.results=infoJob["Results"]       
        logging.debug("Loaded job %s"%self.filename)
        
    def update(self):
        infoJob={"FileName":self.filename,
                    "TempDir":self.dirTemp,
                    "Model":self.modelProt,
                    "Version":self.modelProtVer,
                    "Prediction":self.prediction,
                    "Experiment":self.experiment,
                    "EpochDate":self.dateSub,
                    "ApiJobId":self.jobid,
                    "Status":self.status,
                    "Results":self.results,
                    "SDF":self.sdf
                    }
                    
        with open(self.filename,'w') as infile:
            json.dump(infoJob,infile,indent=4)
        ## for compatibility with apache user in etoxsys
        os.chmod(self.filename, 0775)
        
    def remove(self):
       
        logging.info('Removing job with identifier: %s'%self.filename)
        if self.dirTemp is not None:
            if os.path.exists(self.dirTemp):
                shutil.rmtree(self.dirTemp, ignore_errors=True)
        
        os.remove(self.filename)
        return True
            
    def open2Group(self):
        # For eTOXsys make files accessible to 
        for site in os.walk(self.dirTemp, topdown=True):
            m=oct(os.stat(site[0]).st_mode)[-3:]
            nm=int(m[-3]+m[-3]+m[-1],8)
            os.chmod(site[0], 0775)
    
            for filenm in site[2]:
                pathfn=os.path.join(site[0],filenm)
                m=oct(os.stat(pathfn).st_mode)[-3:]
                nm=int(m[-3]+m[-3]+m[-1],8)
                os.chmod(pathfn, 0775)

def collectJobs(etoxlie_folder=None, cpd=True, prefName=None, sort=True):
    ### Collect all running job. extension:
    ### extFile=job for single compound run
    ###        =dsr for dataset run

    if prefName==None:
        prefName=settings.get('etoxlie_prefixFn')
    if etoxlie_folder==None:
        etoxlie_folder=settings.get('etoxlie_folder', '')

    os.chdir(etoxlie_folder)
    if cpd:
        extFile='job'
    else:
        extFile='dsr'

    listJobFile=glob.glob("%s*.%s"%(prefName,extFile))
    
    listJob=[]
    for jobFn in listJobFile:
        job=jobCpd()
        ### Insert exception in case no loading possible
        try:
            job.load(jobFn)
            listJob.append(job)
        except ValueError,e :
            logging.error("There were error in loading Job %s: %s"%(jobFn,e))
            pass

    
    logging.debug("%d jobs found"%len(listJob))    
    listJob.sort(key=operator.attrgetter('dateSub'))
       
    return listJob


def jobsWithCpd(jobid,etoxlie_folder=None,prefName=None):
    ### Collect job Ids of screenings containing the specified compound. extension:
    ### extFile=job for single compound run
    ###        =dsr for dataset run
    if prefName==None:
        prefName=settings.get('etoxlie_prefixFn')
    if etoxlie_folder==None:
        etoxlie_folder=settings.get('etoxlie_folder', '')
    
    extFile='dsr'
    listJobFile=glob.glob("%s*.%s"%(prefName,extFile))
    listJob=[]

    for jobFn in listJobFile:
        job=jobCpd()
        ### Insert exception in case no loading possible
        try:
            job.load(os.path.join(etoxlie_folder,jobFn))
            for cpd in job.results:
                if cpd['JobName']==jobid:
                    logging.debug("Job %s found in screen %s"%(jobid,job.filename))
                    listJob.append(job.filename)                    
                    break
        except ValueError,e :
            logging.error("There were error in loading Job %s: %s"%(jobFn,e))
            pass

    return listJob

def cpdsInJob():
    pass
