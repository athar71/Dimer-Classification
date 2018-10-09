#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 12:22:42 2017

@author: athar
"""


#import pdb
#pdb.set_trace()

import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd
pymol.finish_launching()

import sys
cmd.set('dot_solvent', 1)
cmd.set('dot_density', 3)



rec_file = sys.argv[1]
lig_file = sys.argv[2]
lig_pymol=lig_file[:-4]
#==============================================================================
#rec_file = '/home/athar/Dimer/dock-std/true/12as/rec.pdb'
#lig_file = '/home/athar/Dimer/dock-std/true/12as/lig.pdb'
#==============================================================================
#complex_file = sys.argv[3]

cmd.load(rec_file)  # use the name of your pdb file
rec_area=cmd.get_area ('rec_nmin')
cmd.load(lig_file)
#lig_area=cmd.get_area ('lig.2')

lig_area=cmd.get_area (lig_pymol)
cmd.save ('complexfile.pdb')
cmd.delete(all)
cmd.load ('complexfile.pdb')
total_area=cmd.get_area ('complexfile')

area= (abs(rec_area+lig_area-total_area))*0.5# using sasa area 
print area
#return area
