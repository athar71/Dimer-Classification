#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 13:27:41 2017

@author: athar
"""
import ssdu_utils as sutils
import sys

rec_file = sys.argv[1]
lig_file = sys.argv[2]

with open(lig_file) as f:
    ligLines = f.read().splitlines()
    
with open(rec_file) as f:
    recLines = f.read().splitlines()
    
sutils.shiftLigAtomsSerial(recLines, ligLines)
with open('new_lig.pdb', 'w') as out_file:
    out_file.write('\n'.join(ligLines))