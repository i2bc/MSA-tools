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

    str_scores = """#
loop_
_ma_qa_metric.id
_ma_qa_metric.mode
_ma_qa_metric.name
_ma_qa_metric.software_group_id
_ma_qa_metric.type
1 global pLDDT 1 pLDDT 
2 local  pLDDT 1 pLDDT 
#
_ma_qa_metric_global.metric_id    1
_ma_qa_metric_global.metric_value {:2.2f}
_ma_qa_metric_global.model_id     1
_ma_qa_metric_global.ordinal_id   1
# 
loop_
_ma_qa_metric_local.label_asym_id
_ma_qa_metric_local.label_comp_id
_ma_qa_metric_local.label_seq_id
_ma_qa_metric_local.metric_id
_ma_qa_metric_local.metric_value
_ma_qa_metric_local.model_id
_ma_qa_metric_local.ordinal_id
"""
    
    regex = "ATOM[\s]*[\d]+[\s]*[A-Z] [A-Z0-9]+[\s]*[\.A-Z0-9][\s]*([A-Z]+)[\s]*([\S]+)[\s]*([\S]+)[\s]*([\d]+)[\s]*[\S]+[\s]*[\-\d\.]+[\s]*[\-\d\.]+[\s]*[\-\d\.]+[\s]*[\-\d\.]+[\s]*([\-\d\.]+)[\s]*[\-\d\.]+[\s]*[\S]+[\s]*[\d]+"
    prev_res = ("0","0","0","0","0")
    global_score = 0
    total = 0
    with open(cif_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                res = re.search(regex,line)
                if not res:
                    os.remove(cif_file)
                    print("\nError: something went wrong with the regex...\n")
                    print(pdb_file)
                    print(line)
                    print(regex)
                    return()
                aa, ch, entity, resnum, bfact = res.groups()
                if " ".join((aa, ch, entity, resnum, bfact)) != " ".join(prev_res):
                    str_scores += f"{ch} {aa} {resnum: <3} 2 {bfact: <5} 1 {resnum: <3} \n"
                prev_res = (aa, ch, entity, resnum, bfact)
                global_score += float(bfact)
                total += 1
    str_scores+="#\n"

    with open(cif_file, "a") as f:
        f.write(str_scores.format(global_score*1./total))
    
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



