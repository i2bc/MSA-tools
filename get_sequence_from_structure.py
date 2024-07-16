#!/usr/bin/env python

from pathlib import Path
import pdb2cif
import argparse
from Bio.PDB.MMCIFParser import MMCIFParser

parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True)
parser.add_argument("-o", required=True)

args = parser.parse_args()


class get_sequence_from_structure:
    def __init__(
        self,
        structure_file,
        output,
    ):
        self.structure_file = structure_file
        self.output = output
        self.prefix_structure_file = Path(self.structure_file).stem

    def run(self):
        self.check_structure_file()
        if self.structure_file_type == ".pdb":
            converted_structure_file = self.prefix_structure_file + ".cif"
            pdb2cif.convert2cif(self.structure_file, converted_structure_file)
            self.structure_file = converted_structure_file
        self.write_fasta()

    def check_structure_file(self):
        valid_extensions = [".pdb", ".cif"]
        if Path(self.structure_file).is_file() and Path(self.structure_file).suffix in valid_extensions:
            self.structure_file_type = Path(self.structure_file).suffix
        else:
            print("wrong structure file type, should be either .pdb or .cif")
            exit(1)

    def get_sequence(self):
        parser = MMCIFParser()
        self.structure = parser.get_structure(structure_id="", filename=self.structure_file)

        one_letter = {}
        one_letter["ALA"] = "A"
        one_letter["CYS"] = "C"
        one_letter["CME"] = "C"
        one_letter["CSE"] = "C"
        one_letter["CSD"] = "C"
        one_letter["CSO"] = "C"
        one_letter["CSS"] = "C"
        one_letter["CCS"] = "C"
        one_letter["P1L"] = "C"
        one_letter["CMT"] = "C"
        one_letter["CSZ"] = "C"
        one_letter["CAS"] = "C"
        one_letter["ASP"] = "D"
        one_letter["GLU"] = "E"
        one_letter["PCA"] = "E"
        one_letter["PHE"] = "F"
        one_letter["GLY"] = "G"
        one_letter["GLZ"] = "G"
        one_letter["HIS"] = "H"
        one_letter["DDE"] = "H"
        one_letter["HIC"] = "H"
        one_letter["NEP"] = "H"
        one_letter["ILE"] = "I"
        one_letter["LYS"] = "K"
        one_letter["KCX"] = "K"
        one_letter["MLY"] = "K"
        one_letter["KCX"] = "K"
        one_letter["LLP"] = "K"
        one_letter["LYZ"] = "K"
        one_letter["LEU"] = "L"
        one_letter["MET"] = "M"
        one_letter["MSE"] = "M"
        one_letter["CXM"] = "M"
        one_letter["FME"] = "M"
        one_letter["ASN"] = "N"
        one_letter["MEN"] = "N"
        one_letter["PRO"] = "P"
        one_letter["GLN"] = "Q"
        one_letter["ARG"] = "R"
        one_letter["ARO"] = "R"
        one_letter["SER"] = "S"
        one_letter["SEP"] = "S"
        one_letter["PN2"] = "S"
        one_letter["THR"] = "T"
        one_letter["TPO"] = "T"
        one_letter["VAL"] = "V"
        one_letter["TRP"] = "W"
        one_letter["TRF"] = "W"
        one_letter["TRQ"] = "W"
        one_letter["TYR"] = "Y"
        one_letter["PTR"] = "Y"
        one_letter["PAQ"] = "Y"

        dic_sequence = {}
        dic_sequence["chain"] = []
        for chain in self.structure.get_chains():
            name_chain = chain.get_id()
            dic_sequence["chain"].append(name_chain)
            dic_sequence[name_chain] = {}
            dic_sequence[name_chain]["order"] = []
            sequence = ""
            for residue in chain.get_residues():
                residue_number = residue.get_id()[1]
                name_res = residue.get_resname()
                if name_res in one_letter:
                    amino_acid = one_letter[name_res]
                else:
                    amino_acid = "x"
                sequence += amino_acid
                dic_sequence[name_chain][residue_number] = amino_acid
                dic_sequence[name_chain]["order"].append(residue_number)
            dic_sequence[name_chain]["sequence"] = sequence

        return dic_sequence

    def write_fasta(self):
        dic_sequence = self.get_sequence()
        sequence = ""
        flag_first = True
        for chain in dic_sequence["chain"]:
            if flag_first:
                sequence = dic_sequence[chain]["sequence"]
                flag_first = False
            else:
                sequence += ":" + dic_sequence[chain]["sequence"]
        with open(self.output, "w") as f:
            f.write(f">{self.prefix_structure_file}\n{sequence}\n")


structure_file = args.i
output = args.o
sequence = get_sequence_from_structure(structure_file, output)
sequence.run()
