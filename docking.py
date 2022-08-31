# Oufan Zhang, Apr 13 2021
# parallels autodock vina docking of a list of smiles strings into the covid spike pocket

# prerequisite: download autodockvina and mgltools (change their paths accordingly),
#               make input and output directories for storing files

import numpy as np
import re
import subprocess
from rdkit.Chem import MolFromSmiles, AddHs, AllChem
from rdkit.Chem.rdmolfiles import MolToPDBBlock, MolToPDBFile

import glob
from multiprocessing import Pool, cpu_count, freeze_support

import time
import os
#cwd = os.path.dirname(os.path.abspath(__file__))

# the path to pdbqt conversion script
lig_path = 'MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'

# the path to pythonsh (mgltools functionalities written in python2)
py_path = 'bin/pythonsh'

#vina = 'autodock_vina_1_1_2_linux_x86/bin/vina'


def parallel(smilelist, cpus=cpu_count()):
    '''
    takes a list of smiles strings as input, 
    and returns zipped tuples of (smiles, binding affinities).
    
    '''
    # record time
    st = time.perf_counter()
    #try:
    #    empty('./input/*')
    #    empty('./output/*')
    #except: 
    #    pass
    freeze_support()
    
    nsmile = len(smilelist)
    args = [(smilelist[n], str(n)) for n in range(nsmile)]
    with Pool(cpus) as pool:   
        out = pool.starmap(docksmile, args)
        
    print('time elapsed: {:0.4f}s'.format(time.perf_counter()-st))
    return out
   
def empty(path):
    for file in glob.glob(path):
        os.remove(file)

def docksmile(smile, filename, gpu=False, return_docked_file=False):

    '''
    coverts a smile string to pdbqt and runs autodock vina,
    returns the binding energy of its top pose
    
    Vina configuration details in config.txt

    gpu: either False (not using GPU), or a string specifying which GPU to use
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
    #pdb = MolToPDBFile(mh, 'input/'+filename+'.pdb', flavor=4)
    pdb = MolToPDBBlock(mh, flavor=4)
    open('/tmp/'+filename+'.pdb', 'w').write(pdb)
    
    # convert pdb to pdbqt
    try:
        out = subprocess.run([py_path, lig_path, '-l', '/tmp/'+filename+'.pdb', '-o','/tmp/'+filename+'.pdbqt'])
    except subprocess.CalledProcessError as e:
        print(e.output)
    if not os.path.exists('/tmp/'+filename+'.pdbqt'):
        print("%s does't exist" % (filename+'.pdbqt'))
        if return_docked_file:
            return smile, np.nan, None
        return smile, np.nan
    
    try:
        if gpu:
            docking_script = "./run_docking_gpu.sh"
            result = subprocess.run(['sh', docking_script, filename, gpu], stdout=subprocess.PIPE)
        else: 
            docking_script = "./run_docking.sh"
            result = subprocess.run(['sh', docking_script, filename], stdout=subprocess.PIPE)
        result = result.stdout.decode('utf-8')
    except subprocess.CalledProcessError as er:
        print(er.output)
        print(smile, '; file:%s.pdbqt'%filename)
        if return_docked_file:
            return smile, np.nan, None
        return smile, np.nan
    
    #print(filename+'.pdbqt','docking success')

    # read energy from output
    energy = np.nan
    strings = re.split('\n', result)
    for line in strings:
        if line[0:4] == '   1':
            energy = float(re.split(' +', line)[2])
    #print(energy )  
    if return_docked_file:
        if os.path.exists('/tmp/'+filename+'_out.pdbqt'):
            with open('/tmp/'+filename+'_out.pdbqt') as f:
                docked_file = f.read()
            return smile, energy, docked_file
        else:
            return smile, energy, None
    return smile, energy

# for testing
#import pandas as pd

if __name__ == '__main__':
    import sys
    import logging
    import argparse

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

    parser = argparse.ArgumentParser()
    parser.add_argument('-smiles', help='the smiles for docking', default="OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N")
    parser.add_argument('-l', '--ligand', help='the ligand .pdbqt file used for docking')
    parser.add_argument('-g', '--gpu', help="GPU index used for docking", default=None)

    args = parser.parse_args()

    logger.info("n_cpus:" + str(cpu_count()))
    t = time.time()

    hashed_smiles = str(hash(args.smiles) % ((sys.maxsize + 1) * 2)) 
    logger.info("accepted smiles: " + args.smiles)

    vina_score = docksmile(args.smiles, hashed_smiles, args.gpu)
    logger.info("time: %.1f" % (time.time() - t))
    print(vina_score[0], vina_score[1])

    with open(f"/tmp/{hashed_smiles}_out.pdbqt") as f:
        print(f.read())




