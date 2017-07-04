# Create template for docking. Perform it wih amber FF (all atoms) and gromacs for minimization

import os, sys
import subprocess32 as sp
import importlib
from tempfile import mkdtemp
from shutil import rmtree

import pybel as pb
import openbabel as ob
from glob import glob

from eTOX_ALLIES.etox import settings
from eTOX_ALLIES.etox.core import modelHandler
from eTOX_ALLIES.etox.topology.amber import refineTau
from eTOX_ALLIES.etox.core.geometry import findOcpd1, normhem, distance as dist

AMBERHOME = settings.get('AMBERHOME')
ACPYPE = settings.get('ACPYPE')
DATADIR = os.path.join(settings.get('etoxlie_root_dir'), 'bin/files/')
GROMACSHOME = settings.get('GROMACSHOME')

ob.obErrorLog.SetOutputLevel(0)


def prepareModel(modelDir,modelData,etoxlie_folder,radiusRes=16):
    '''modelData={name: str, dockSoftware: str, timeSim: float, forceField:str, resSite: None or list
                  pH: float, pHCorr: bool, protConfs: list of dictionaries }
        protConfs: [{'pdb' : str, 'center' : None / [X,Y,Z], 'radius' : float }]
    '''
    
    currDir=os.getcwd()
    success=True
    results=''
    scratchDir=mkdtemp(dir=etoxlie_folder)
    try:
        print modelDir
        for i in modelData:
            print i
            if i != 'protConfs':
                print modelData[i]
       
        ## FROM DOCKING MODULE LOAD REQUIRED FORMAT
        dockFMT=(importlib.import_module("eTOX_ALLIES.etox.docking.%s"%modelData['dockSoftware'])).INPUTFMT 
        ## FROM TOPOLOGY MODULE LOAD appropriate function
        prepareGMX=(importlib.import_module("eTOX_ALLIES.etox.topology.%s"%modelData['forceField'])).prepareGMX

        # preliminary check consistency structures
        seqs=[]
        for conf in modelData['protConfs']:
            seqs.append(getSeq(conf['pdb']))

        if len(seqs)>1:
            for i in range(1, len(seqs)):
                if seqs[i] != seqs[0]:
                    raise Exception('Protein conformation %d has different sequence than conformation 1'%(i+1))

        #process protein conformation
        os.chdir(scratchDir)
        nconf=1
        proteinParams=[]
        protTop=None
        if modelData['resSite'] is None:
            listResSite=[]
            
        for conf in modelData['protConfs']:
            print '##################################################'
            print '########### CONF %d              ##############'%nconf
            print '##################################################'
            preOut="conf_%d"%nconf
            tempPdb="%s.pdb"%preOut
            with open(tempPdb, 'w') as outFile:
                outFile.write(conf['pdb'])
            if nconf==1:
                #get listchanges
                outpdb,listChanges,nhem,ncyp,cypCoord=getListChanges(tempPdb)
                print listChanges, nhem,ncyp,cypCoord
            else:
                outpdb=tempPdb

            # 1. Prepare docking system
            # COMMON: amber to add H, ambpdb for pqr (charges+coords)
            # output is a pybel molecule
            
            print "PREPARE DOCKING"
            refineConfMol,prmtop,inpcrd=refineTau(outpdb, listChanges, nhem, ncyp, "ref_%s"%preOut, cypCoord=cypCoord)
            # pybel molecule is saved in the required format depending on the docking software
            dockTemplate=saveMol(refineConfMol,"DT_%s"%preOut, dockFMT)

            # 1B. Center of docking
            
            print "GET CENTER DOCKING"
            if conf['center'] is None:
                centerDock=CypCenter(refineConfMol)
                centerDock=centerDock.tolist()
            else:
                centerDock=conf['center']

            print "PREPARE TOPOLOGY/ GRO"
            # 2. Prepare topology/ coordinates FF dependent
            temptop,gro, tempattype, tempposre, tempcharge=prepareGMX(refineConfMol,prmtop,inpcrd,'MD_%s'%preOut,protTop)

            if nconf==1:
                protTop=temptop
                attype=tempattype
                posre=tempposre
                chargeProt=tempcharge


            # 3. get residues around center docking    
            print "GET RESLIST"
            if modelData['resSite'] is None:
                tempResList=getSurrResids(refineConfMol,centerDock,radiusRes,listExclude=['WAT','SOL','HOH'])
                listResSite=listResSite+tempResList

            print "GET PROTEINPARAMS"
            #print listResSite
            #print centerDock
            proteinParams.append({
                                  'proteinDock' : dockTemplate,
                                  'proteinCoor' : gro,
                                  'pocket' : centerDock,
                                  'radius' : conf['radius'],
                                  'filter' : centerDock,  #HEME Iron coordinates
                                  })

            nconf+=1
        # when everything is ready: copy files,
        if modelData['resSite'] is None:
            listResSite=','.join(map(str, set(listResSite)))
        else:
            listResSite=modelData['resSite']

        print "PREPARE DICTTOMODEL"
        #prepare dictionary describing model
        dictToModel={'name' : modelData['name'],
                     'forceField' : modelData['forceField'],
                     'pHCorr' : modelData['pHCorr'],
                     'pH' : modelData['pH'],
                     'nprot' : len(proteinParams),
                     'charge' : chargeProt,
                     'timeSim' : modelData['timeSim'],
                     'proteinTop' : protTop,
                     'miscFiles' : [attype],
                     'proteinTopPos' : posre,    
                     'resSite' : listResSite,
                    'dockSoftware' : modelData['dockSoftware'],
                    'proteinParams' : proteinParams,
                    'replicas' : 2
            }

        
        success,results=modelHandler.createModel(dictToModel,modelDir)

    except Exception, e:
        success=False
        results=e

    os.chdir(currDir)
    rmtree(scratchDir)
    return success,results


def getSeq(pdb):
    prot=pb.readstring('pdb',pdb)
    seq=[]
    for res in ob.OBResidueIter(prot.OBMol):
        seq.append([res.GetNum(),res.GetName()])
    
    return seq


def getListChanges(pdbmodel):
    # get modification to apply at pdb to apply amber FF and heme parameters from J Comput Chem. 2012 Jan 15;33(2):119-33
    prefnm,extnm=os.path.splitext(pdbmodel)
    if extnm[1:]!='pdb':
        mutref=pb.readfile(extnm[1:],pdbmodel).next()
        pdbmodel="%s.pdb"%prefnm
        mutref.write('pdb',pdbmodel,overwrite=True)

    outpdb="%s_tau.pdb"%prefnm
    reduceCmd=os.path.join(AMBERHOME,'bin','reduce')
    cmd=[reduceCmd,'-HIS','-FLIPs','-quiet',pdbmodel]

    red=sp.Popen(cmd,stdout=sp.PIPE,stderr=sp.PIPE)
    optimutpdb,err=red.communicate()
    
    optimut=pb.readstring('pdb',optimutpdb)

    listchanges=[]
    listTer=[]
    sCoords={}
    nhemseq=None
    checkCyp=True
    firstRes=True
    isCyp=False
    cypCoord=True
    ncypseq=None
    for res in ob.OBResidueIter(optimut.OBMol):
        resname=res.GetName()
        resNum=res.GetNum()
        if firstRes:
            nres=resNum
            firstRes=False
        else:
            nres+=1

        #FIX HIS TAUTOMERS
        if resname=='HIS':
            atmname=[]
            for atm in ob.OBResidueAtomIter(res):
                atmname.append(res.GetAtomID(atm).strip())
            if ('HE2' in atmname) and ('HD1' in atmname):
                newResName="HIP"
            elif 'HE2' in atmname:
                newResName="HIE" 
            elif 'HD1' in atmname:
                newResName="HID"
            else:
                print "Error in determining the protonation state of HIS %d"%resNum
            listchanges.append((resNum,'HIS',newResName))

        #GET Fe coordinates for coordinating CYS identification
        elif resname=='HEM':
            isCyp=True
            nhem=resNum
            nhemseq=nres
            #nhem=nres
            atmname=[]
            for atm in ob.OBResidueAtomIter(res):
                if atm.GetType()=='FE' or atm.GetType()=='Fe' or res.GetAtomID(atm).strip()=='FE' or res.GetAtomID(atm).strip()=='Fe':
                    Fe_coord=[atm.GetX(),atm.GetY(),atm.GetZ()]
        
        #Collect every cysteine S coordinates
        elif resname=='CYS':
            atmname=[]
            for atm in ob.OBResidueAtomIter(res):
                atmName=res.GetAtomID(atm).strip()
                if atmName=='SG':
                    sCoords[resNum]={}
                    sCoords[resNum]['coor']=[atm.GetX(),atm.GetY(),atm.GetZ()] #resNum
                    sCoords[resNum]['nseq']=nres

        if resname=='CYP':
            ncyp=resNum
            ncypseq=nres
            checkCyp=False

        if resname=='HIH':
            ncypseq=nres
            checkCyp=False
            cypCoord=False 
       
        for atm in ob.OBResidueAtomIter(res):
            atmName=res.GetAtomID(atm).strip()
            if atmName=="OXT":
                listTer.append(resNum)

    if isCyp:
    #Find the coordinating CYS and add to listchanges
        if cypCoord:
           if checkCyp:
                CYP={'nres':0,'dist':100}
                for res in sCoords:
                    distS=dist(sCoords[res]['coor'],Fe_coord)
                    print res, distS
                    if distS<CYP['dist']:
                        CYP['dist']=distS
                        CYP['nres']=res
                        CYP['nseq']=sCoords[res]['nseq']
        
                ncypseq=CYP['nseq']
                listchanges.append((CYP['nres'],'CYS','CYP'))
        

    optimut.removeh()
    optimut.write('pdb',outpdb,overwrite=True)

    return (os.path.join(os.getcwd(),outpdb),listchanges,nhemseq,ncypseq,cypCoord)
     

def saveMol(prot,outPref, outFmt):  
    # apply eventual patches due to amber top to openbabel format differences and output docking format file
    currDir=os.getcwd()

    # adjust H
    for atom in prot:
        if atom.atomicnum==1:
            atom.OBAtom.SetType('H')
    
    for res in ob.OBResidueIter(prot.OBMol):
        resname=res.GetName()
        if resname=='HEM':
            #nhem=res.GetNum()
            for atm in ob.OBResidueAtomIter(res):
                if atm.GetType()=='FE':
                    atm.SetType('Fe')
                if atm.GetType()=='Du':
                    atName=res.GetAtomID(atm).lstrip()
                    atm.SetType(atName[0])

    print "save mol2 for docking"
    templfn=os.path.join(currDir,"%s.%s"%(outPref,outFmt))
    out=pb.Outputfile(outFmt,templfn, overwrite=True)
    out.write(prot)
    out.close()

    return templfn

   
def CypCenter(protMol):
    # Find Fe coordinates
    # Find Virtual O (cpdI) coordinates
    # Find equations of the hem plane
    listheme=['C2B','CHB','C3A','C2A','CHA','C3D','C2D','CHD','C3C','C2C','CHC','C3B']
    Fe_coord=None
    S_coord=None
    hemeplane_coor=dict()
    for res in ob.OBResidueIter(protMol.OBMol):
        resname=res.GetName()
        if resname=='HEM':
            for atm in ob.OBResidueAtomIter(res):
                atmname=res.GetAtomID(atm).strip()
                if atmname in listheme:
                    hemeplane_coor[atmname]=[atm.GetX(), atm.GetY(), atm.GetZ()]
                elif atm.GetType()=='FE' or atm.GetType()=='Fe' or atmname=='FE' or atmname=='FE':
                    Fe_coord=[atm.GetX(),atm.GetY(),atm.GetZ()]

    #get normal of the plane of the heme ring
    hemnorm=normhem(hemeplane_coor)
    catcenter=findOcpd1(hemnorm,Fe_coord,dist=6.0)

    return catcenter


def getSurrResids(pybelMol,centerDock,radiusRes,listExclude=[]):
    listRes=[]
    print 'GETSURRRESIDUE'
    print centerDock, radiusRes
    for res in ob.OBResidueIter(pybelMol.OBMol):
        if not res.GetName() in listExclude:
            inner=False
            for atm in ob.OBResidueAtomIter(res):
                coords=[atm.GetX(),atm.GetY(), atm.GetZ()]
                if dist(coords,centerDock)<radiusRes:
                    inner=True
                    break
            if inner:
                listRes.append(res.GetNum())
    return listRes
    
    