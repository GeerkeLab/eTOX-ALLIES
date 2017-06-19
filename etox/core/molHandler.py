r'''This file contains the core functions needed to do affinity calculations for etoxsys.'''

import logging
import os.path
import sys
import re
import pybel as pb


def file2mol(format,f, addH=False,opt=False):
    '''Translate a file with given format to a pybel molecule
       Does some optimizations if necessary returns pybel molecule'''
    molecule=pb.readfile(format, f).next()  #make molecule object
    if addH:
        molecule.addh()                      #add explicit hydrogens
    if  molecule.dim != 3:   				#unless the molecule is already 3-dimensional
        logging.info( '%d-dimensional molecule, making it 3D.'%molecule.dim)
        molecule.make3D(forcefield='mmff94',steps=100)   	#make it 3-dimensional
    if opt:
        molecule.localopt(forcefield='mmff94', steps=500) #if requested, molecular geometry can be optimized further
    
    return molecule


def rotateMolecule(molecule,xyzangle=[0,0,0,0]):
    '''rotate_molecule(molecule,xyzangle=[0,0,0,0]) -> molecule
       rotate pybel molecule around x,y,z by given angle returns the molecule'''
    x,y,z,angle=xyzangle
    m = pb.ob.matrix3x3()
    m.RotAboutAxisByAngle(pb.ob.vector3(x,y,z), angle) #calculates a rotation matrix around the axes specified, by angle specified
    rotarray = pb.ob.doubleArray(9)
    m.GetArray(rotarray) #convert the matrix to an array
    molecule.OBMol.Rotate(rotarray)   #rotates the supplied molecule object
    
    return molecule


def makeRotations(molecule, rotations=[[1,0,0,90],[1,0,0,-90],[0,1,0,90],[0,1,0,-90],[0,0,1,90],[0,0,1,-90]]):
    '''Takes a pybel molecule and array of rotations to perform on it
       Returns array with rotated pybel molecules'''
    rotated_mols=list()
    for i in range(len(rotations)):
        rotamol = copyMolecule(molecule)             #make copies to store all rotated molecules
        rotateMolecule(rotamol,rotations[i])            #perform rotations on the copies
        rotated_mols.append(rotamol)
        logging.debug("Rotating x,y,z,angle %s"%(rotations[i]))
    
    return rotated_mols


def mols2single(mol_array,output='mol',format='pdb'):
    '''Takes an array of pybel molecules
       Creates a file containing all as separate molecules in the requested format (default as pdb)
       Returns the path of the file created'''
    outname = '%s.%s'%(output,format)
    output_f= pb.Outputfile(format, outname , overwrite=True)
    logging.debug("Molecules will be saved to %s"%outname)
    if os.path.isfile(outname): logging.warning("Overwriting file %s!"%outname) #in case file already exists, make a note in log
    
    for i in mol_array:
        output_f.write(i)
        logging.debug("Added one molecule to %s"%outname)
    output_f.close()
    
    return os.path.abspath(outname)


def copyMolecule(molecule):
  '''Takes a pybel molecule as input.
  Effectively returns a deep copy of the pybel molecule submitted'''
  #using copy.deepcopy didnt return a usable molecule
  newmolecule =pb.Molecule(pb.ob.OBMol(molecule.OBMol)) #converts pybel into OpenBabel and back
  return newmolecule