# Oufan Zhang, Apr 13 2021
# parallels autodock vina docking of a list of smiles strings into the covid spike pocket

# prerequisite: download autodockvina and mgltools (change their paths accordingly),
#               make input and output directories for storing files

import numpy as np
import re
import subprocess
from rdkit.Chem import MolFromSmiles, AddHs, AllChem
from rdkit.Chem.rdmolfiles import MolToPDBBlock, MolToPDBFile

import time
import os
#cwd = os.path.dirname(os.path.abspath(__file__))

# the path to pdbqt conversion script
lig_path = 'mgltools_x86_64Linux2/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'
# smina uses obabel as alternative

# the path to pythonsh (mgltools functionalities written in python2)
py_path = 'mgltools_x86_64Linux2/bin/pythonsh'

#vina = 'autodock_vina_1_1_2_linux_x86/bin/vina'


def smile_pose_generator(smile, nconf, filename):

    '''
    coverts a smile string to pdbqt and runs smina to sample conformations
    ========
    nconf: np.int
       number of conformations to generate
    '''
    #print(smile, filename)
    if not isinstance(smile, str):
        raise TypeError('Input is not a class of string')
        
    m = MolFromSmiles(smile)
    # assert valid smiles
    if m is None:
        raise ValueError(smile, 'is not a valid smile string')
    mh = AddHs(m)
    embed = AllChem.EmbedMolecule(mh, useRandomCoords=False)
    
    #check if rdkit successfully generates structure
    if embed!=0:
        print('RDkit fails to embed molecule', smile, '; file:%s.pdb'%filename)
        return smile, np.nan
        
    # generate mol2 file
    #pdb = MolToPDBFile(mh, 'input/'+filename+'.pdb', flavor=4)
    pdb = MolToPDBBlock(mh, flavor=4)
    open('/tmp/'+filename+'.pdb', 'w').write(pdb)
    
    # convert mol2 to pdbqt by open babel
    try:
        out = subprocess.run(['obabel', '-imol2', '/tmp/'+filename+'.pdb', '-opdbqt', '-O', '/tmp/'+filename+'.pdbqt'])
    except subprocess.CalledProcessError as e:
        print(e.output)
    if not os.path.exists('/tmp/'+filename+'.pdbqt'):
        print("%s does't exist" % (filename+'.pdbqt'))
        return smile, np.nan
    
    filelist = []
    try:
        result = subprocess.run(['sh', 'run_smina_randomize.sh', filename], stdout=subprocess.PIPE)
        result = result.stdout.decode('utf-8')
    except subprocess.CalledProcessError as er:
        print(er.output)
        print(smile, '; file:%s.pdbqt'%filename) 
    
    return smile, filelist

# for testing
#import pandas as pd

if __name__ == '__main__':
    import sys
    import logging



