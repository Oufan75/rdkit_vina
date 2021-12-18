# Oufan Zhang, Dec 17 2021
# autodock vina docking of a list of smiles strings 

# prerequisite: download smina and rdkit

import numpy as np
import re, os
import subprocess
from rdkit.Chem import MolFromSmiles, AddHs, AllChem
from rdkit.Chem.rdmolfiles import MolToPDBBlock

'''
# the path to pdbqt conversion script
lig_path = '/home/oufan/Desktop/autodock/mgltools_x86_64Linux2/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'

# the path to pythonsh (mgltools functionalities written in python2)
py_path = '/home/oufan/Desktop/autodock/mgltools_x86_64Linux2/bin/pythonsh'

#vina = 'autodock_vina_1_1_2_linux_x86/bin/vina'
'''


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
        
    # generate pdb file
    pdb = MolToPDBBlock(mh, flavor=4)
    open(filename+'.pdb', 'w').write(pdb)
    
    # convert pdb to pdbqt by open babel
    try:
        out = subprocess.run(['sh', 'run_obabel.sh', filename])
    except subprocess.CalledProcessError as e:
        print(e.output)
    if not os.path.exists(filename+'.pdbqt'):
        print("%s does't exist" % (filename+'.pdbqt'))
        return smile, np.nan
    
    # generate random conformations
    filelist = []
    for n in range(nconf):
        try:
            result = subprocess.run(['sh', 'run_smina.sh', filename, str(n)], stdout=subprocess.PIPE)
            result = result.stdout.decode('utf-8')
            filelist.append(filename+'_%i.pdbqt'%n)
        except subprocess.CalledProcessError as er:
            print(er.output)
            print(smile, '; file:%s.pdbqt randomize failed'%filename) 
    
    return smile, filelist


# for testing
import pandas as pd

if __name__ == '__main__':

    s, flist = smile_pose_generator('CCCCCC', 2, 'test')
    print(flist)



