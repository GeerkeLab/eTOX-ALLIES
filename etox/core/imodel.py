# -*- coding: utf-8 -*-

##    Description    eTOXlab model template
##                   
##    Authors:       Manuel Pastor (manuel.pastor@upf.edu) 
##
##    Copyright 2013 Manuel Pastor
##
##    This file is part of eTOXlab.
##
##    eTOXlab is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation version 3.
##
##    eTOXlab is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with eTOXlab.  If not, see <http://www.gnu.org/licenses/>.
import sys
import os
import logging
from time import sleep 

from eTOX_ALLIES.etox import settings
from eTOX_ALLIES.etox.core import jobHandler
from eTOX_ALLIES.etox.core.main import submitScreen, screenSDF

etoxlie_folder=settings.get('etoxlie_folder')

class imodel():
    def __init__ (self, vpath, jobid='<not defined>'):
        """
        LIE settings
        
        The jobid argument was added by Marc van Dijk, 21-10-2016.
        It passes along the Celery job identifier responsible for 
        launching this calculation from the etoxwsapi side of things.
        This API jobid is stored in the .dsr and .job files greated
        by eTOXLie and enables API job cancelation by querying for 
        the jobid and touching an empty KILL file in the associated
        eTOXlie job directory.
        """
        
        self.protein='REPLACE_PROTEIN'
        self.version=REPLACE_VERSION
        self.jobid=jobid

        
##################################################################
##    WORKFLOW METHODS
##################################################################

    def buildWorkflow(self, molecules):

        msg="To create or calibrate a new model do not use eTOXlab but the dedicated web interface"

        return (False, msg)


    def predictWorkflow(self, molecules, detail, progress):
        # output of the fuction should be:
        # list of (True, (molPR,molAD,molCI)) for each compound

        # output for each compound is a list of:
        #    (molPR,molAD,molCI)
        # molPR: (True, float with predicted y)
        # molAD: (True, integer-sum of violations)
        # molCI: (True, float- SDEP and violations based value)

        launchDir=os.getcwd()
        logging.info('WORKING DIRECTORY IS %s'%launchDir)
        statEndJob=['FAILED','CANCELLED','DONE']
        completed=False
        modelId={'modelProt':self.protein,'modelProtVer':self.version}
        
        # Check existence file input and check number of molecules for preliminary 'Failed' result setting
        nmols=0
        if os.path.isfile(molecules):
            with open(molecules,'r') as infile:
                sdf=infile.read()
            for line in sdf.splitlines():
                if line=='$$$$':
                    nmols+=1

        listPreds=[ [False,[(False,0),(False,0),(False,0)]] for i in range(nmols)]

        try:
            newSdfFn=molecules
            prediction=True
            fieldExp=''
            # Submit job
            success, results=submitScreen(newSdfFn,modelId,prediction,etoxlie_folder=etoxlie_folder,fieldExp=str(fieldExp),jobid=self.jobid)
            #if success results is the filename id
            if not success:
                raise Exception('Submission failed: %s'%results)

            # Load job as class
            jobFN=results
            jobScreen=jobHandler.jobCpd()
            jobScreen.load(jobFN)

            # Check job status
            while jobScreen.status not in statEndJob:
                updated, status,results=screenSDF(jobScreen.results, jobScreen.prediction, jobScreen.modelProt, jobScreen.modelProtVer,jobScreen.experiment)
                if updated:               
                    jobScreen.status=status
                    jobScreen.results=results
                    jobScreen.update()
                sleep(10)
            # if it is finished, get the results in the proper format
            for i,cpd in enumerate(jobScreen.results):
                doneCpd=False
                molPR=(False,0)
                molAD=(False,0)
                molCI=(False,0)

                if cpd['Status']=='DONE':
                    molPR=(True,cpd['DGcalc'])
                    molAD=(True,cpd['CI'])
                    molCI=(True,cpd['Err'])
                    doneCpd=True

                listPreds[i]=[doneCpd,[molPR,molAD,molCI]]

                # Clean single job
                cpdFileNm=cpd['JobName']
                listScreens=jobHandler.jobsWithCpd(cpdFileNm)
                # if more than one: stop
                if len(listScreens) == 1:
                    # Update screening job
                    cpdJob=jobHandler.jobCpd()
                    cpdJob.load(cpdFileNm)
                    #clean mol folder and index
                    cpdJob.remove()


            jobScreen.remove()
            completed=True
            
        except Exception,e:
            print e
        
        os.chdir(launchDir)

        return (completed,listPreds)
