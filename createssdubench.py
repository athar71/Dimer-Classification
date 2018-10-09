#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 11:16:19 2017

@author: athar
"""

"""
#This scrip will generate a benchmark for ssdu on Dimer data,
 it generate files that are corresponding to clusters for our complexes
"""
import os
import ssdu_utils as sutils
import numpy as np
#from sblu.ft import read_ftresults
import subprocess
import  scipy.stats as stats
#import interface_calculation_func
from pathlib import Path



rootDirectories = ["/projectnb/mhcpep/athar/Dimer/files-testset/dimer", "/projectnb/mhcpep/athar/Dimer/files-testset/monomer"]
#rootDirectories = ["files-testset/dimer/dock-std/true", "files-testset/monomer/dock-std/false"]

#The name of ft file to use for each complex.
ftName = "ft.000.00"



num_clusters=5
for file_num,rootDir in enumerate(rootDirectories,0):


    complexList = sutils.comp_list(rootDir)
    
    
    
    
    for counter,comp in enumerate(complexList,0):
        
        compPath = os.path.join(rootDir, comp)
        lig_add=compPath+'/lig.pdb'
        rec_add=compPath+'/rec.pdb'
        
        ft_address=os.path.join(compPath, "%s" %ftName)
        
        #Reading clusters: The output from sutils.read_clusters is a list where
        #each element is a tuple :(center, listOfMembers) and INDICES START FROM ZERO WHEN python_index is set to True
        clusterFileName = os.path.join(compPath, "cluster.%s" %ftName)
        clusters = sutils.read_clusters(clusterFileName, python_index=True)

       
        #ft=np.loadtxt(ft_address)
        with open (ft_address) as ft:
            ft_lines=ft.readlines()
        
        
        # calling the commands for generatinf the file in terminal, nimn file, pbd concate
        os.chdir(compPath)
        subprocess.call(["pdbprep.pl","lig.pdb"])
        subprocess.call(["pdbprep.pl", "rec.pdb"])
        subprocess.call(["pdbnmd.pl", "--smod", "L","--xplor-psf","lig", "?"])
        subprocess.call(["pdbnmd.pl", "--smod", "R","--xplor-psf","rec", "?"])
        subprocess.call(["python","/projectnb/mhcpep/athar/Dimer/shift_ligand.py", "rec_nmin.pdb", "lig_nmin.pdb"])
        subprocess.call(["cp", "new_lig.pdb", "lig_nmin.pdb"])
        subprocess.call(["pdb_psf_concat", "rec_nmin_xplor.psf", "rec_nmin.pdb","lig_nmin_xplor.psf","lig_nmin.pdb", "--prefix=rec_nminlig_nmin"])
        subprocess.call(["obabel", "rec_nminlig_nmin.pdb","-O", "r_l_u_nmin.mol2"])
        
        for cluster_name in range(0,min(num_clusters,len(clusters))):
            filename="ssdu_cluster_%s_%d" %(comp ,cluster_name+1)
            file_name = os.path.join("/projectnb/mhcpep/athar/Dimer/ssdu_files_testset", filename)
            subprocess.call(["mkdir", file_name])
            #subprocess.call(["cp", "new_lig.pdb", os.path.join(compPath,filename,"lig_nmin.pdb")])
        
            subprocess.call(["cp", "rec_nmin.pdb","rec_nmin.psf", "lig_nmin.pdb","lig_nmin.psf","rec_nminlig_nmin.pdb", "rec_nminlig_nmin.psf","r_l_u_nmin.mol2", file_name])
            cluster_members=clusters[cluster_name][1]
            #ft_cluster=ft[cluster_members][:]
            fname= os.path.join(file_name,"ft.%s" %filename)
            f=open(fname,'w')
            for line in cluster_members:
                f.write(ft_lines[line])
            f.close()
            
           
