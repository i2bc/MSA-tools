#!/usr/bin/env python

import os
import sys
from pathlib import Path
import re
import logging
import tempfile
import pdb2cif
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO

RATE4SITE = "rate4site"
MAFFT = "mafft"


class r4s_multi:
    def __init__(
        self,
        msa_input,
        output,
        dir_temp,
        structure_file,
        diverged=30,
        conserved=60,
        threshold_percent_gap=25,
    ):
        self.msa_input = msa_input
        self.output = output
        self.dir_temp = dir_temp
        self.structure_file = structure_file
        self.diverged = diverged
        self.conserved = conserved
        self.threshold_percent_gap = threshold_percent_gap
        self.dic_score = {}
        self.conservation_data = {}
        self.mafft_result = {}
        self.test = {}
        self.setup_logging()
        self.check_structure_file()
        if self.structure_file_type == ".pdb":
            prefix_structure_file = Path(structure_file).stem
            converted_structure_file = prefix_structure_file + ".cif"
            pdb2cif.convert2cif(self.structure_file, converted_structure_file)
            self.structure_file = converted_structure_file

    def run(self):
        self.check_fasta()
        if self.complex:
            self.get_subsequence()
            self.create_fasta()
            for index in self.fasta_sequences:
                fasta = self.fasta_sequences[index]["fasta"]

                self.computeR4S(fasta, self.output + str(index))
                self.get_score_from_r4s(self.output + str(index) + ".grade", index)
                self.compute_bfact(index)
            self.assign_sequence_cif_fasta()
            for group in self.link_sequence_chain:
                self.align_pdb_r4s(group["index"], group["occurrence"], group["chain"])
            self.create_cif_converge_diverge_r4s()
        else:
            self.computeR4S(self.msa_input, self.output)
            self.get_score_from_r4s(self.output + ".grade", 0)
            self.compute_bfact(0)
            self.assign_sequence_cif_fasta()
            for group in self.link_sequence_chain:
                self.align_pdb_r4s(group["index"], group["occurrence"], group["chain"])
            self.create_cif_converge_diverge_r4s()

    def setup_logging(self):
        FORMAT = "%(levelname)s - %(message)s - %(asctime)s"
        DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
        logging.basicConfig(format=FORMAT, datefmt=DATE_FORMAT)
        self.logger = logging.getLogger("rate4site")
        self.logger.setLevel(logging.INFO)
        fh_detailed = logging.FileHandler("rate4site.log")
        fh_detailed.setLevel(logging.INFO)
        fh_detailed.setFormatter(logging.Formatter(fmt=FORMAT, datefmt=DATE_FORMAT))
        self.logger.addHandler(fh_detailed)

    def check_structure_file(self):
        valid_extensions = [".pdb", ".cif"]
        if (
            Path(self.structure_file).is_file()
            and Path(self.structure_file).suffix in valid_extensions
        ):
            self.structure_file_type = Path(self.structure_file).suffix
        else:
            self.logger.error(
                "wrong structure file type, should be either .pdb or .cif"
            )
            sys.exit(1)

    def check_fasta(self):
        """Detect if the fasta has multiple protein in it and store these information (length of these different sequences and their number)"""
        self.complex = False
        with open(self.msa_input, "r") as f:
            firstline = f.readline()
            if firstline.startswith("#"):
                self.complex = True
                firstline_no_comma = firstline[1:]
                pattern_split_in_two_group = r"[\d,]+"
                length_sequences, number_of_sequences = re.findall(
                    pattern=pattern_split_in_two_group, string=firstline_no_comma
                )
                length_sequences = length_sequences.split(",")
                number_of_sequences = number_of_sequences.split(",")
                self.length_number_couple = []
                for index, length in enumerate(length_sequences):
                    self.length_number_couple.append(
                        (length, number_of_sequences[index])
                    )

    def get_subsequence(self):
        """get all the subsequence from a fasta"""
        self.fasta_sequences = {}

        with open(self.msa_input, "r") as f:
            lines = f.readlines()
            current_start_position = 0

            for index, tup in enumerate(self.length_number_couple):
                length = current_start_position + int(tup[0])
                self.fasta_sequences[index] = {}
                self.fasta_sequences[index]["header"] = []
                for line in lines:

                    stripped_line = line.strip()
                    sequence = ""
                    if not line.startswith("#"):
                        if line.startswith(">"):
                            header = stripped_line[1:]
                            header = header.split("\t")
                            current_header = header[index]
                            self.fasta_sequences[index][current_header] = ""

                        else:
                            sequence = stripped_line[current_start_position:length]
                            self.fasta_sequences[index]["header"].append(current_header)
                            self.fasta_sequences[index][current_header] = {}
                            self.fasta_sequences[index][current_header][
                                "sequence"
                            ] = sequence
                            self.fasta_sequences[index][current_header][
                                "occurrence"
                            ] = tup[1]
                current_start_position = length

    def create_fasta(self):
        count = 1
        for index in self.fasta_sequences:
            fasta = f"split_{count}.fasta"
            self.fasta_sequences[index]["fasta"] = fasta
            with open(fasta, "w") as f:
                for header in self.fasta_sequences[index]["header"]:
                    f.write(
                        f"> {header}\n{os.path.join(self.dir_temp, self.fasta_sequences[index][header]['sequence'])}\n"
                    )
            count += 1

    def computeR4S(self, msa_file, result):
        os.system(f"{RATE4SITE}  -s {msa_file} -o {result}.grade -x {result}.ph")

    def get_score_from_r4s(self, input, index):
        sequence = ""
        res_order = []
        score = {"score": {}, "percentage_no_gap": {}}
        with open(input, "r") as f:
            lines = f.readlines()
            for line in lines:
                if not (line.startswith("#") or line == "\n"):
                    parts = line.split()
                    pos = parts[0]
                    seq = parts[1]
                    score_value = float(parts[2])
                    percentage_no_gap = parts[-1].split(
                        "/"
                    )  # get the last character because sometime there are space in qq-interval, messing the positions.
                    a = percentage_no_gap[0]
                    b = percentage_no_gap[1]
                    computed_percentage_no_gap = float(a) / float(b) * 100.0

                    sequence += seq
                    res_order.append(pos)
                    score["score"][pos] = score_value
                    score["percentage_no_gap"][pos] = computed_percentage_no_gap
        self.dic_score[index] = {
            "sequence": sequence,
            "res_order": res_order,
            "score": score["score"],
            "percentage_no_gap": score["percentage_no_gap"],
        }
        if self.complex:
            self.dic_score[index]["occurrence"] = self.fasta_sequences[index][
                self.fasta_sequences[index]["header"][0]
            ]["occurrence"]
        else:
            self.dic_score[index]["occurrence"] = 1

    def compute_bfact(self, index=0):
        conservation_data = {}
        conservation_data["diverged"] = {}
        conservation_data["conserved"] = {}
        conservation_data["percentage_no_gap"] = self.dic_score[index][
            "percentage_no_gap"
        ]
        conservation_data["no_gap"] = {}
        conservation_data["order"] = {}
        conservation_data["bin"] = {}
        scores = self.dic_score[index]["score"].values()
        amplitude = max(scores) - min(scores)
        bfact_min = min(scores)

        for res in self.dic_score[index]["res_order"]:
            MAX_SCORE = 99.0
            MIN_SCORE = 1.0

            bfact = MAX_SCORE - (
                MAX_SCORE
                * (self.dic_score[index]["score"][res] - bfact_min)
                / amplitude
            )

            if bfact > MAX_SCORE:
                conservation_data["bin"][res] = MAX_SCORE
            elif bfact < MIN_SCORE:
                conservation_data["bin"][res] = MIN_SCORE
            else:
                conservation_data["bin"][res] = bfact

            if bfact < self.diverged:
                conservation_data["diverged"][res] = True
            else:
                conservation_data["diverged"][res] = False

            if bfact > self.conserved:
                conservation_data["conserved"][res] = True
            else:
                conservation_data["conserved"][res] = False

            if conservation_data["percentage_no_gap"][res] > self.threshold_percent_gap:
                conservation_data["no_gap"][res] = True
            else:
                conservation_data["no_gap"][res] = False

        self.conservation_data[index] = conservation_data

    def extract_sequence_cif(self):
        parser = MMCIFParser()
        self.structure = parser.get_structure(
            structure_id="", filename=self.structure_file
        )

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
                    amino_acid = "X"
                sequence += amino_acid
                dic_sequence[name_chain][residue_number] = amino_acid
                dic_sequence[name_chain]["order"].append(residue_number)
            dic_sequence[name_chain]["sequence"] = sequence

        return dic_sequence

    def assign_sequence_cif_fasta(self):
        pdb_sequences = self.extract_sequence_cif()
        for index in self.dic_score:
            sequence = self.dic_score[index]["sequence"]

            self.mafft_result[index] = {}
            for chain in pdb_sequences["chain"]:
                input_mafft = tempfile.NamedTemporaryFile(delete=False)
                with open(input_mafft.name, "w") as f:
                    f.write(
                        f"> Fasta\n{sequence}\n> Structure_{self.structure_file_type}\n{pdb_sequences[chain]['sequence']}"
                    )
                output_mafft = tempfile.NamedTemporaryFile(delete=False)

                command = MAFFT + f" --auto {input_mafft.name} > {output_mafft.name}"
                os.system(command)
                self.mafft_result[index][chain] = {}
                self.mafft_result[index][chain]["index_cif"] = pdb_sequences[chain][
                    "order"
                ]
                self.mafft_result[index][chain]["sequences"] = (
                    self.get_sequence_from_alignment(output_mafft.name)
                )
        self.compare_mafft_alignment()

    def compare_mafft_alignment(self):
        pdb_sequences = self.extract_sequence_cif()
        chains = pdb_sequences["chain"]
        self.link_sequence_chain = []
        for index in self.mafft_result:
            for occurrence in range(int(self.dic_score[index]["occurrence"])):
                percentage = {
                    "max": 0,
                    "index": "",
                    "chain": "",
                }
                for chain in chains:
                    # TODO add variable
                    percentage_identity = self.compare_sequences(
                        self.mafft_result[index][chain]["sequences"][
                            self.mafft_result[index][chain]["sequences"]["order"][0]
                        ],
                        self.mafft_result[index][chain]["sequences"][
                            self.mafft_result[index][chain]["sequences"]["order"][1]
                        ],
                    )
                    if percentage_identity > percentage["max"]:
                        percentage["max"] = percentage_identity
                        percentage["index"] = index
                        percentage["chain"] = chain
                        percentage["occurrence"] = occurrence
                self.link_sequence_chain.append(percentage)
                chains.remove(percentage["chain"])

    def compare_sequences(self, seq1, seq2):
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of the same length")

        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)

        percentage_similarity = (matches / len(seq1)) * 100

        return percentage_similarity

    def get_sequence_from_alignment(self, alignment):
        sequences = {}
        sequences["order"] = []
        current_id = None
        seq = ""
        with open(alignment, "r") as f:
            lines = f.readlines()

            for line in lines:
                if line.startswith(">"):
                    if current_id is not None:
                        sequences[current_id] = seq
                        sequences["order"].append(current_id)
                    header = line.strip()
                    current_id = header
                    seq = ""
                else:
                    seq += line.strip()
            sequences[current_id] = seq
            sequences["order"].append(current_id)

        return sequences

    def align_pdb_r4s(self, index, occurrence, chain):
        seq_r4s = self.mafft_result[index][chain]["sequences"][
            self.mafft_result[index][chain]["sequences"]["order"][0]
        ]
        seq_cif = self.mafft_result[index][chain]["sequences"][
            self.mafft_result[index][chain]["sequences"]["order"][1]
        ]
        index_cif = self.mafft_result[index][chain]["index_cif"]

        ite_index = 0  # iterator of index_pdb (true number in the pdb file)
        ite_r4s = 1  # iterator of the r4s_value["bin"] number in the r4s file

        r4s_value_correspondence = {}
        r4s_value_correspondence["diverged"] = {}
        r4s_value_correspondence["conserved"] = {}
        r4s_value_correspondence["percentage_no_gap"] = {}
        r4s_value_correspondence["no_gap"] = {}
        r4s_value_correspondence["order"] = []
        r4s_value_correspondence["bin"] = {}

        for spdb, sr4s in zip(seq_cif, seq_r4s):

            if spdb != "-":

                cur_index_cif = index_cif[ite_index]
                ite_index += 1  # we go to the next index

                if sr4s != "-":
                    r4s_value_correspondence["bin"][cur_index_cif] = (
                        self.conservation_data[index]["bin"][str(ite_r4s)]
                    )
                    r4s_value_correspondence["conserved"][cur_index_cif] = (
                        self.conservation_data[index]["conserved"][str(ite_r4s)]
                    )
                    r4s_value_correspondence["diverged"][cur_index_cif] = (
                        self.conservation_data[index]["diverged"][str(ite_r4s)]
                    )
                    r4s_value_correspondence["percentage_no_gap"][cur_index_cif] = (
                        self.conservation_data[index]["percentage_no_gap"][str(ite_r4s)]
                    )
                    r4s_value_correspondence["no_gap"][cur_index_cif] = (
                        self.conservation_data[index]["no_gap"][str(ite_r4s)]
                    )
                    r4s_value_correspondence["order"].append(cur_index_cif)

                    ite_r4s += 1
                else:
                    r4s_value_correspondence["bin"][cur_index_cif] = 0.0
                    r4s_value_correspondence["conserved"][cur_index_cif] = 0.0
                    r4s_value_correspondence["diverged"][cur_index_cif] = 0.0
                    r4s_value_correspondence["percentage_no_gap"][cur_index_cif] = 0.0
                    r4s_value_correspondence["no_gap"][cur_index_cif] = 0.0
                    r4s_value_correspondence["order"].append(cur_index_cif)

            else:
                if sr4s != "-":
                    ite_r4s += 1
        if index not in self.test:
            self.test[index] = {}
        self.test[index][occurrence] = {}
        self.test[index][occurrence][chain] = r4s_value_correspondence

    def create_cif_converge_diverge_r4s(self):
        with open(self.output + ".cif", "w"):
            for residue in self.structure.get_residues():
                for atom in residue.get_atoms():
                    atom.set_occupancy(0)

            for index in self.test:
                for occurrence in self.test[index]:
                    for chain in self.test[index][occurrence]:
                        struct_chain = self.structure[0][chain]
                        for residue in struct_chain.get_residues():
                            residue_number = residue.get_id()[1]
                            if (
                                residue_number
                                in self.test[index][occurrence][chain]["bin"]
                            ):
                                for atom in residue.get_atoms():
                                    atom.set_occupancy(
                                        self.test[index][occurrence][chain]["bin"][
                                            residue_number
                                        ]
                                    )

        io = MMCIFIO()
        io.set_structure(self.structure)
        io.save(self.output + ".cif")
