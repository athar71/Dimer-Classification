#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 14:34:00 2017

@author: athar
"""

import subprocess
#==============================================================================
# import pdb
# pdb.set_trace()
#==============================================================================

#import interface_calculation
#a=interface_calculation.interface_area('/home/athar/Dimer/dock-std/true/12as/rec.pdb','/home/athar/Dimer/dock-std/true/12as/lig.pdb')
area=subprocess.check_output(["./interface_calculation.py","/home/athar/Dimer/dock-std/true/12as/rec.pdb","/home/athar/Dimer/dock-std/true/12as/lig.pdb"])
area=float (area)