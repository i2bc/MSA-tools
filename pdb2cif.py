#! /usr/bin/python
"""Short script to convert pdb to cif in light of the colabfold project only (to get a cif compatible with PDBe Molstar plugin)
"""

from Bio.PDB import PDBParser, MMCIFIO
import os
import re
import sys

usage = """usage: pdb2cif.py input.pdb
This script will take a pdb file outputted by AF2 as input and will output the structure in cif format suitable for PDBe Molstar
"""


def get_args(args):

    if len(args) != 2:
        print("\nError: Please provide a single pdb file as input\n")
        print(usage)
        sys.exit()

    pdb_file = args[1]

    if not os.path.exists(pdb_file):
        print("\nError: pdb file doesn't exist\n")
        print(usage)
        sys.exit()

    if os.path.splitext(pdb_file)[1].lower() != '.pdb':
        print("\nError: Please provide a _pdb_ file as input (i.e. *.pdb)\n")
        print(usage)
        sys.exit()
    
    return pdb_file


def convert2cif(pdb_file,cif_file):

    
    p = PDBParser()
    struc = p.get_structure("data_"+os.path.basename(os.path.splitext(cif_file)[0]), pdb_file)
    io = MMCIFIO()
    io.set_structure(struc)
    io.save(cif_file)
    print("pdb2cif done!")


def main():

    pdb_file = get_args(sys.argv)
    overwrite = True
    cif_file = os.path.splitext(pdb_file)[0]+".cif"
    if not os.path.exists(cif_file) or overwrite:
        convert2cif(pdb_file,cif_file)
    else:
        print(cif_file,"already exists")


if __name__ == '__main__':

    main()



