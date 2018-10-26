#!/usr/bin/env python

#This script reads ft files and generates pairwise-rmsd and cluster files.

import os
import ssdu_utils as sutils
import subprocess
import multiprocessing as mp

#import pdb
#pdb.set_trace()

#rootDirectories = ["dock-std/true", "dock-std/false"]
rootDirectories = ["files-testset/dimer", "files-testset/monomer"]
#The name of ft file to use for each complex.
ftName = "ft.000.00"
rotationFilePath = "/home/athar/Dimer/rot_mol2.prm"

def create_pwrmsd_cluster(compPath):
    
    #Generating pairwise-rmsd file
    pairwiseRMSDName = os.path.join(compPath, "pwrmsd.%s" %ftName)
    ligFileName = os.path.join(compPath, "lig_nmin.pdb")
    recFileName = os.path.join(compPath, "rec_nmin.pdb")
    ftFileName = os.path.join(compPath, ftName)

    subprocess.call(["sblu", "measure", "pwrmsd", ligFileName, ftFileName, rotationFilePath, "-n", "1000", "--only-CA", "--interface-radius", "10", "--rec", recFileName, "--output", pairwiseRMSDName ])
    
    #Generating cluter files
    clusterFileName = os.path.join(compPath, "cluster.%s" %ftName)
    subprocess.call(["sblu", "docking", "cluster", "--output", clusterFileName, pairwiseRMSDName])

#Using multiprocesses to speed up
numProcesses = 4
pool = mp.Pool(numProcesses)
for rootDir in rootDirectories:

    complexList = sutils.comp_list(rootDir)
    
    compCounter = 0
    numComplexes = len(complexList)
    retObjs = numComplexes*[0]
    for comp in complexList:

        compPath = os.path.join(rootDir, comp)
        retObjs[compCounter] = pool.apply_async(create_pwrmsd_cluster, (compPath,))

        compCounter += 1

    for objs in retObjs:
        objs.wait()

pool.close()

