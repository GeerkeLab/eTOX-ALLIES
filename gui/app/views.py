import os
import sys
import logging
import time
import pybel
import json

from flask import render_template, url_for, request, Response, send_file
from numpy import corrcoef
from gui.app import app

from eTOX_ALLIES import __rootpath__
from eTOX_ALLIES.etox.core import jobHandler
from eTOX_ALLIES.etox.core import modelHandler
from eTOX_ALLIES.etox.core import plotting
from eTOX_ALLIES.etox.core.main import submitScreen, updateJobsStatus
import eTOX_ALLIES.etox.core.modelTemplates as modTemp
import eTOX_ALLIES.etox.core.eTOXlab_integration as eLab

json_settings = os.path.join(__rootpath__, 'data/settings.json')
with open(json_settings) as st:
    settings = json.load(st)


def jobstoDict(listJobs):
    dictList=[]
    for job in listJobs:
        
        id=os.path.split(job.filename)[1]
        jobDict={
                'name': job.filename,
                'id':id,
                'model': job.modelProt,
                'version': job.modelProtVer,
                'status': job.status,
                'cpds': job.results,
                'ncpd':len(job.results),
                'date': time.strftime("%Y/%m/%d\n%H:%M:%S", time.localtime(job.dateSub))
                }
        if job.prediction:
            jobDict['type']='PRED'
        else:
            jobDict['type']='CAL'
        
        dictList.append(jobDict)
    
    return dictList


def listModVers(modelDir):
    listModels=modelHandler.listModels(modelDir)
    listVers=[]
    for model in listModels:
        listVers=listVers+[ver for ver in model['vers']]
    listVers=list(set(listVers))
    return (listModels,listVers)

@app.route('/')
@app.route('/index')
@app.route('/home')
def home():
    return render_template('index.html')


@app.route('/submit')
def submit():
    listModels,listVers=listModVers(settings.get('modelDir'))
    return render_template('submit.html',listModels=listModels,listVers=listVers)


@app.route('/createJob', methods =['POST'])
def createJob():
    
    try:
        modelId={}
        modelId['modelProt']=request.form['ProtMod']
        
        #check model
        listModels=modelHandler.listModels(settings.get('modelDir'))
        model=next((item for item in listModels if item["name"] == modelId['modelProt']),None)
        if model is None:
            raise Exception('Model protein not found')
        
        # check file
        sdffile=request.files['sdfFile']
        if sdffile.filename !='':
            newSdfFn=os.path.split(sdffile.filename)[1]
            newSdfFn=os.path.join(settings.get('etoxlie_folder'),newSdfFn)
            sdffile.save(newSdfFn)
        else:
            raise Exception('You forgot to add the file!')
        if request.form['prediction'] == 'true':
            prediction=True
            fieldExp=''
            modelId['modelProtVer']=int(request.form['ProtModVer'])
            if modelId['modelProtVer'] not in model['vers']:
                raise Exception('Model version not found')
            misc=None
        else:
            prediction=False
            fieldExp=request.form['fieldActivity']
            
            isGamma=False
            if 'isGamma' in request.form:
                isGamma=True
            fixBeta=False
            if 'fixBeta' in request.form:
                fixBeta=True
            misc={'isGamma':isGamma,'fixBeta':fixBeta}
            
            if len(model['vers'])>0:
                vers=sorted(model['vers'])[-1]+1
            else:
                vers=1
            modelId['modelProtVer']=vers
        
        success, results=submitScreen(newSdfFn,modelId,prediction,etoxlie_folder=settings.get('etoxlie_folder'),fieldExp=str(fieldExp),misc=misc)
    except Exception,e:
        success=False
        results='Error in request processing: check your input! %s'%e
    
    if 'newSdfFn' in locals():
        if os.path.exists(newSdfFn):
            os.remove(newSdfFn)
    
    if success:
        msg='Job successfully submitted (ID: %s) <br> <a href="/">Go to the initial page</a>'%results
        logging.info('Job successfully submitted (ID: {0})'.format(results))
    else:
        msg='Job submission failed: %s <br> <a href="/">Go to the initial page</a>'%results
        logging.info('Job submission failed: {0}'.format(results))
    
    return msg
    

@app.route('/models')
def model():
    listModels,listVers=listModVers(settings.get('modelDir'))
    return render_template('models.html',listModels=listModels,listVers=listVers)


@app.route('/models/descModel/<prot>')
def modelDesc(prot):
    protDir=os.path.join(settings.get('modelDir'),prot)
    success, model=modelHandler.loadModel(protDir, prediction=False, verModel=0, modelFN='model.dat')
    if not success:
        return "Model not found :'-("
    
    # Source model parameters
    return render_template('modelDesc.html',model=model)


@app.route('/models/descModVer/<prot>/<ver>')
def modelDescVer(prot,ver):
    protDir=os.path.join(settings.get('modelDir'),prot)
    success, (model,params)=modelHandler.loadModel(protDir, prediction=True, verModel=int(ver), modelFN='model.dat')
    
    # compute pearson r, spearman s
    Gcalc=[float(x['Gcalc']) for x in params['trainSet']]
    Gexp=[float(x['Gexp']) for x in params['trainSet']]
    
    idxs=sorted(list(range(len(Gcalc))), key=lambda x: Gcalc[x])
    rankGcalc=[0]*len(idxs)
    for i, x in enumerate(idxs):
        rankGcalc[x] = i
    
    idxs=sorted(list(range(len(Gexp))), key=lambda x: Gexp[x])
    rankGexp=[0]*len(idxs)
    for i, x in enumerate(idxs):
        rankGexp[x] = i
    
    s=corrcoef(rankGexp,rankGcalc)[0,1]
    r=corrcoef(Gexp,Gcalc)[0, 1]
    
    params['r']=r
    params['s']=s
    params['date']=time.strftime("%Y/%m/%d\n%H:%M:%S", time.localtime(params['creation_date']))
    
    img=plotting.plotLIE(params['trainSet'])
    svg=plotting.fixSVG(img,h=400)
    
    return render_template('modelVersDesc.html',param=params,plot= svg)#send_file(a, mimetype='image/svg+xml'))


@app.route('/models/descModPlot/<prot>/<ver>')
def modPlot():
    a=plotting.plotLIE()
    return send_file(a, mimetype='image/svg+xml')


@app.route('/models/newModel')
def newModel():
    return render_template('newModel.html')


@app.route('/models/newModel/create', methods =['POST'])
def createModel():
    modelData={}
    try:
        modelData['name']= request.form['modName']
        modelData['forceField'] = 'amber'
        modelData['timeSim'] = float(request.form['time'])
        modelData['dockSoftware'] = 'paradocks'
        modelData['pH'] = 7.4
        modelData['pHCorr'] = True
        if 'isIon' in request.form:
            if request.form['isIon'] == 'true':
                modelData['pHCorr'] = False
        getResList=False
        if 'guessRes' in request.form:
            if request.form['guessRes'] == 'true':
                getResList=True
        if getResList:
            modelData['resSite'] = None
        else:
            modelData['resSite'] = request.form['resList']
        protList=[]
        nprot=0
        for pdb in request.files:
            idProt=pdb.strip('pdbFile')
            guessCoords=False
            if 'guessCoords%s'%idProt in request.form:
                if request.form['guessCoords%s'%idProt] == 'true':
                    guessCoords=True
            if guessCoords:
                center=None
            else:
                X=float(request.form['Xcoord%s'%idProt])
                Y=float(request.form['Ycoord%s'%idProt])
                Z=float(request.form['Zcoord%s'%idProt])
                center=[X,Y,Z]
            radius=float(request.form['radius%s'%idProt])
            content=request.files[pdb].read()
            
            protList.append({
                             'pdb' : content,
                             'center' : center,
                             'radius' : radius
                             })
        
        modelData['protConfs']=protList
        success,results=modTemp.prepareModel(settings.get('modelDir'),modelData,settings.get('etoxlie_folder'),radiusRes=16)
        if not success:
            raise Exception(results)
        
        msg='New model %s has been successfully created <br> <a href="/home">Go to the initial page</a>'%results
    except Exception, e:
        msg='%s<br><a href="/models/newModel">Go back</a>'%e
    
    return msg


@app.route('/models/toElab')
def connectToELAB():
    success,msg=eLab.linkElab(keyElab='ETOXLABHOME', overWrite=True)
    if success:
        return '%s <br> <a href="/models">Go to the model page</a>'%msg
    else:
        return 'Error occurred during the linking: %s <br> <a href="/models">Go to the model page</a>'%msg


@app.route('/jobs')
def jobs():
    updateJobsStatus()
    listJobs=jobstoDict(jobHandler.collectJobs(cpd=False, sort=True))
    return render_template('jobs.html',listJobs=listJobs)


@app.route('/showMol/<smi>')
def showMol(smi):
    
    smi=smi.replace('&equal&','=')
    smi=smi.replace('&hash&','#')
    smi=smi.replace('&slash&','/')
    
    mol=pybel.readstring('smi',str(smi))
    svgString=mol.write('svg',opt={"a":True,"P":300})
    
    return svgString


@app.route('/theory')
def theory():
    return 'In preparation. <br> <a href="/">Go to the initial page</a>'


@app.route('/about')
def about():
    return render_template('about.html')


@app.route('/Download/<typeReq>/<id>/<key>')
def downloadFile(typeReq,id,key):
    if typeReq=='model':
        sdir=settings.get('modelDir')
    elif typeReq=='job':
        sdir=settings.get('etoxlie_folder')
    
    filenm=os.path.join(sdir,id,key)
    with open(filenm,'r') as infile:
        out=infile.read()
    
    return Response(out, mimetype='text/csv')


@app.route('/cleanJob', methods =['POST'])
def cleanJob():
    jobid=request.form['jobName']
    jobDir,jobFn=os.path.split(jobid)
    job=jobHandler.jobCpd()
    job.load(jobFn)
    if job.dirTemp is not None:
        # if cpd: ceck in all mains. if only one, clean from main
        #1 get list of screenings with cpd
        listScreens=jobHandler.jobsWithCpd(job.filename)
        # if more than one: stop
        if len(listScreens) > 1:
            return 'Multiple screening are using this job. It won\'t be removed. <br> <a href="/">Go to the initial page</a>'
        elif len(listScreens) == 1:
            # Update screening job
            screening=jobHandler.jobCpd()
            screening.load(listScreens[0])
            screening.results=[x for x in screening.results if x['JobName'] != os.path.join(jobDir,jobFn) ]
            # clean from main job
            screening.update()
            #clean mol folder and index
            job.remove()
    else:
        for item in job.results:
            # check if used by other jobs, in case rmove job cpd.
            jobCpd=jobHandler.jobCpd()
            jobCpd.load(item['JobName'])
            listScreens=jobHandler.jobsWithCpd(jobCpd.filename)
            # if more than one: do not remove
            if len(listScreens) == 1:
                jobCpd.remove()
        job.remove()
    
    return '%s, %s REMOVED <br> <a href="/jobs">Go to the list jobs page</a>'%(jobDir, jobFn)
