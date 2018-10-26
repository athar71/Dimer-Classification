#!/usr/bin/env python
"""Takes ligand file, applies rotation/translation from specified ft file line,
then saves new ligand"""

import os
import sys
import numpy as np
#from path import path
from itertools import islice
from prody import parsePDB, writePDB, writePDBStream, AtomGroup
from sblu.ft import apply_ftresults_atom_group

FTRESULTS_DTYPE = np.dtype([('roti', 'i4'),
                            ('tv', ('f16', 3)),
                            ('E', 'f8'),
                            ('Ev', ('f8', 5))])

def read_rotations(rotations_path):
    with open(rotations_path, "r") as f:
        rotations = np.loadtxt(f)
    if rotations.shape[-1] == 10:
        rotations = rotations[:, 1:]
    return rotations.reshape(-1, 3, 3)

def read_ftresults(ftresult_path, limit):
    with open(ftresult_path, "r") as f:
        ftresults = np.loadtxt(islice(f, 0, limit),
                               dtype=FTRESULTS_DTYPE)
        ftresults.sort(order='E')
    return ftresults

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument("--ftfile-line", "-n", type=int, default=0, help="0 indexed")
    parser.add_argument("--limit", "-l", type=int, default=None,
                        help="limit for number of ft lines read")
    parser.add_argument("ligand", help="ligand used for docking")
    parser.add_argument("ft_file", help="ftresults file to draw transformations from.")
    parser.add_argument("rot_file", help="Rotation file used during PIPER run to make ft_file.")
    parser.add_argument("--ft_limit", help="Speficy how many ftfile lines to read")

    args = parser.parse_args()

    rotations = read_rotations(args.rot_file)

    #read ft_files
    ftresults = read_ftresults(args.ft_file, args.limit)

    #get center information for original ligands
    lig_orig = parsePDB(args.ligand)
    lig_center = np.mean(lig_orig.getCoords(), axis=0)
    lig_name = args.ligand
    lig_base = lig_name.rsplit('/', 1)[-1][:-4]

    #apply ft lines to the aligned ligands
    coords_apply = apply_ftresults_atom_group(lig_orig, ftresults, rotations, center=lig_center)

    #create new atom group for translated ligand
    lig_new = AtomGroup('%s.%s.pdb'%(lig_base, args.ftfile_line))

    #only choose one coordinate set (0 indexed)
    coords_apply.setACSIndex(args.ftfile_line)
    lig_new.setCoords(coords_apply.getCoords())
    lig_new.setNames(lig_orig.getNames())
    lig_new.setResnames(lig_orig.getResnames())
    lig_new.setResnums(lig_orig.getResnums())

    #write new ligand file
    writePDB("%s.%s.pdb"%(lig_base, args.ftfile_line), lig_new)

if __name__ == "__main__":

    sys.exit(main())
