#!/usr/bin/env python

#This script reads cluster files and saves a feature file

import os
import ssdu_utils as sutils
import numpy as np
from sblu.rmsd import srmsd
from sblu.ft import read_ftresults
import subprocess
import  scipy.stats as stats
#import interface_calculation_func
from pathlib import Path
from math import floor




rootDirectories = ["dock-std/true", "dock-std/false"]
#rootDirectories = ["files-testset/dimer", "files-testset/monomer"]

#The name of ft file to use for each complex.
ftName = "ft.000.00"
rotationFilePath = "/home/athar/Dimer/rot_mol2.prm"

complexList_total=[*sutils.comp_list("dock-std/true"),*sutils.comp_list("dock-std/false")]
#complexList_total=[*sutils.comp_list("files-testset/dimer"),*sutils.comp_list("files-testset/monomer")]
#with open("complex_list_test.txt","w") as f:
with open("complex_list_train.txt","w") as f:
    for complex in complexList_total:
        f.write(complex+"\n")
print (f) 
numComplexes = len(complexList_total)
retObjs = numComplexes*[0]
num_clusters=5
num_features=10*num_clusters+2+(1+8)*num_clusters#coment on features, it can br 2+8 if we add rmsd
features=np.empty((numComplexes,num_features))
Y_label=np.zeros(numComplexes)
not_count=0
tmp_list=[]
for file_num,rootDir in enumerate(rootDirectories,0):


    complexList = sutils.comp_list(rootDir)
    
    
    #if rootDir=='files-testset/dimer':
    if rootDir=='dock-std/true':
        label=1
    else:
        label=-1

    
    for counter,comp in enumerate(complexList,0):
        
        
        compPath = os.path.join(rootDir, comp)
        lig_add=compPath+'/lig_nmin.pdb'
        rec_add=compPath+'/rec_nmin.pdb'

        pairwiseRMSDName = os.path.join(compPath, "pwrmsd.%s" %ftName)
        #np.loadtxt read, root square, reshape(n,n)
        pwd=np.loadtxt(pairwiseRMSDName)
        n=int((len(pwd))**0.5)
        pwr_rmsd_matrix=np.reshape(pwd,(n,n))
        ft_address=os.path.join(compPath, "%s" %ftName)
        
        #Reading clusters: The output from sutils.read_clusters is a list where
        #each element is a tuple :(center, listOfMembers) and INDICES START FROM ZERO WHEN python_index is set to True
        clusterFileName = os.path.join(compPath, "cluster.%s" %ftName)
        clusters = sutils.read_clusters(clusterFileName, python_index=True)

        #reading ft files and extracting the energy column
        #FTRESULT_DTYPE = np.dtype([('roti', 'i4'), ('tv', ('f8', 3)), ('E', 'f8'), ('E1', 'f8'), ('E2', 'f8'), ('E3', 'f8'), ('E4', 'f8'), ('E5', 'f8')])
        FTRESULT_DTYPE = np.dtype([('roti', 'i4'), ('tv', ('f8', 3)), ('E', 'f8'), ('E_components', ('f8',5))])
        #ft=read_ftresults(ft_address,limit=1000)
        ft=np.loadtxt(ft_address,dtype=FTRESULT_DTYPE)
        ft=ft[0:1000][:]
        energy=ft['E']
        energy_components=ft['E_components']
        #Calculate features for the complex
        large_clus_center=clusters[0][0]
        #num_clusters=len(clusters)
        
        cluster_size=np.zeros(num_clusters)
        cluster_center_energy=np.empty(num_clusters)*np.nan
        cluster_center_energy1=np.empty(num_clusters)*np.nan
        cluster_center_energy2=np.empty(num_clusters)*np.nan 
        cluster_center_energy3=np.empty(num_clusters)*np.nan
        cluster_center_energy4=np.empty(num_clusters)*np.nan
        cluster_center_energy5=np.empty(num_clusters)*np.nan                              
        mean_dist_from_center=np.empty(num_clusters)*np.nan
        var_dist_from_center=np.empty(num_clusters)*np.nan
        cluster_interface_area=np.empty(num_clusters)*np.nan
        cluster_PCA_eig = np.empty(num_clusters*5)*np.nan
        cluster_SDU_eig = np.empty(num_clusters*3)*np.nan
        cluster_SSDU_energy = np.empty(num_clusters)*np.nan
        cluster_SSDU_rmsd = np.empty(num_clusters)*np.nan
        #calculating features for each cluster
        
        for cluster_name in range(0,min(num_clusters,len(clusters))):
            
            cluster_size[cluster_name]=len(clusters[cluster_name][1]) 
            cluster_center=clusters[cluster_name][0]
            cluster_center_energy[cluster_name]=energy[cluster_center]
            cluster_center_energy1[cluster_name]=energy_components[cluster_center][0]
            cluster_center_energy2[cluster_name]=energy_components[cluster_center][1]
            cluster_center_energy3[cluster_name]=energy_components[cluster_center][2]
            cluster_center_energy4[cluster_name]=energy_components[cluster_center][3]
            cluster_center_energy5[cluster_name]=energy_components[cluster_center][4]
            #energy should be a matrix of the 5th column of ft.000.00 file that 
            #shows the energy of conformation
            distall_from_center=pwr_rmsd_matrix[cluster_center]
            dist_from_center=list(distall_from_center[i] for i in  clusters[cluster_name][1])
            mean_dist_from_center[cluster_name]=np.mean(dist_from_center)
            var_dist_from_center[cluster_name]=np.var(dist_from_center)
            #finding th einterface area of the cluster centers
            #first we need to make the rotation and the translation of the ligand with respect to the line of the ft for center
            #ftapply_singleline.py ligorig.pdb ft.000.00 ~/prms/rot70k.0.0.6.prm -l 1000 -n 2
#==============================================================================
#             import pdb
#             pdb.set_trace()
#==============================================================================

            subprocess.call(["./ftapply_singleline.py",lig_add,ft_address,rotationFilePath,"-l","1000","-n",str(cluster_center)])
            #cluster_interface_area[cluster_name]=interface_calculation.interface_area(rec_add, 'lig.2.pdb')
            rotate_lig_name=os.path.join("lig_nmin.%s.pdb" %cluster_center)
            
            area=subprocess.check_output(["./interface_calculation.py",rec_add,rotate_lig_name])
            #my_file = Path(rotate_lig_name)
            cluster_interface_area[cluster_name] = float (area)
            
            
            """
            Adding ssdu features, for each complex, consider ssdu features for the clusters with more than 50 member
            average energy of top 20% re-sample conformations,
            PCA and SDU eigenvalues, average rmsd of the top 10%  re-sample conformations with the cluster center 
            """
            folder_ssdu_name = "ssdu_cluster_%s_%d" %(comp ,cluster_name+1)
            eig_ssdu_add = os.path.join("ssdu_dock_std", folder_ssdu_name,"ssdu_run" ,"features_result.txt")
            energy_ssdu_add = os.path.join("ssdu_dock_std", folder_ssdu_name ,"ssdu_run", "conform_energy_out_1_1.txt")
            pdb_ssdu_add = os.path.join("ssdu_dock_std", folder_ssdu_name ,"ssdu_run", "locally_minimized_pdbs")
            
            my_file = Path(eig_ssdu_add)
            if my_file.is_file():
                #reading eigenvalues 
                #5 PCA eigenvalues
                eig = open(eig_ssdu_add)
                lines = eig.readlines()
                PCA_eig_string = lines[1][:]
                PCA_eig = [float(x) for x in PCA_eig_string.split()]
                cluster_PCA_eig[cluster_name*5:5+cluster_name*5] = PCA_eig
                
                #3 SDU eigenvalues
                SDU_eig_string = lines[3][:]
                SDU_eig = [float(x) for x in SDU_eig_string.split()]
                cluster_SDU_eig[cluster_name*3:3+cluster_name*3] = SDU_eig
                
                
                #reading average of energy for top 10% re sampled conformations
                ssdu_energies = np.loadtxt (energy_ssdu_add)
                sort_energy_index = np.argsort(ssdu_energies)
                num_ten_per = floor(0.1*cluster_size[cluster_name])
                ssdu_mean_energy = np.mean (ssdu_energies[sort_energy_index[0:num_ten_per]])
                cluster_SSDU_energy[cluster_name] = ssdu_mean_energy
        #            #sort the energy and pick the top 10% complexes
                    
#            
#            #finding the average of pairwised rmsd of the top 10% with the former cluster center
#            ssdu_rmsd_center = np.empty(num_pdbs_ssdu)
#            for pdb_number_file in range(num_pdbs_ssdu):
#                
#                add_pdb =
#                #at first I want to type it as command, but I founf the function
#                #A = subprocess.Popen(["sblu", "measure", "srmsd", add_pdb, rotate_lig_name], stdout=subprocess.PIPE).communicate()[0])
#                #Out[7]: b'50.2652210043\n'
#                ssdu_rmsd_center[pdb_number_file] = srmsd()
#            
#            #finding the avrage pairwise rmsd of the top 10% less energy complexes
#            ssdu_rmsd_pwr = np.empty(num_ten_per)
#            
#            for i in range(num_ten_per):
#                ssdu_rmsd_pwr_each = np.empty(num_ten_per)
#                conf_num = sort_energy_index[i]
#                pdb_file_name = "model.%d" %conf_num
#                
            
                subprocess.call(["rm",rotate_lig_name])
#===========================================np.loadtxt===================================
#           my_file = Path(rotate_lig_name)
#             if my_file.is_file():
#                 cluster_interface_area[cluster_name]=float (area)
#                 subprocess.call(["rm",rotate_lig_name])
#             else:
#                 tmp_list.append(comp)
#                 not_count=not_count+1
#==============================================================================
        #Global features( for all clusters)
        dist_large_cluster=pwr_rmsd_matrix[large_clus_center]
        #mean distance of cluster centers form the cluster center of the largest cluster
        mean_dist_global=np.mean(dist_large_cluster)
        var_dist_global=np.var(dist_large_cluster)
        #features order: for each complex we have cluster size(number of clusters that we have)+mean and var of 
        #cluster memebers' distance from center+ Global features: mean and var of cluster centers from largest cluster center
        
        
        #index=counter+file_num*len(sutils.comp_list("files-testset/dimer"))
        index=counter+file_num*len(sutils.comp_list("dock-std/true"))
        features[index]=[*cluster_size,*cluster_center_energy,*cluster_center_energy1,*cluster_center_energy2,*cluster_center_energy3,*cluster_center_energy4,*cluster_center_energy5,*mean_dist_from_center,*var_dist_from_center,*cluster_interface_area,*cluster_PCA_eig,*cluster_SDU_eig,*cluster_SSDU_energy,mean_dist_global,var_dist_global]
                    
        Y_label[index]=label

col_mean = np.nanmean(features,axis=0)
#Find indicies that you need to replace
inds = np.where(np.isnan(features))

#Place column means in the indices. Align the arrays using take
features[inds]=np.take(col_mean,inds[1])
#np.savetxt('features_train.txt',features)
#np.savetxt('label_train.txt',Y_label)
#np.savetxt('features_test.txt',features)
#np.savetxt('label_test.txt',Y_label)
np.savetxt('features_train_ssdu.txt',features)
np.savetxt('label_train_ssdu.txt',Y_label)
#f=open('',"w")
#f.writline([])
#write(complex)
#f=close
  
