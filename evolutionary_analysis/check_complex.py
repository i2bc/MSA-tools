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

import util_cif

RATE4SITE = "rate4site"
MAFFT = "mafft"
HHFILTER = "hhfilter"

# TODO: Delete when done with local testing
# RATE4SITE = "~/programmation/stage/script_msa_tools/tools/bin/rate4site"
# MAFFT = "~/programmation/stage/script_msa_tools/tools/bin/mafft"
# HHFILTER = "~/programmation/stage/script_msa_tools/hhsuite/bin/hhfilter"


class r4s_multi:
    def __init__(
        self,
        msa_input,
        output_directory,
        temporary_directory,
        structure_file,
        array_chains=[],
        diverged=30,
        conserved=60,
        threshold_percent_gap=25,
        maximum_percentage_identity=80,
        maximum_number_sequences=100,
        split_identity=True,
        multiple_msa=False,
        field="occupancy", #occupancy or bfactor
    ):
        self.msa_input = msa_input
        self.structure_file = structure_file
        self.array_chains = array_chains
        self.output_directory = output_directory
        self.temporary_directory = temporary_directory
        self.msa_hhfilter = os.path.join(temporary_directory, "msa_hhfilter.fasta")
        self.msa_filtered = os.path.join(temporary_directory, "msa_filtered.fasta")
        self.diverged = diverged
        self.conserved = conserved
        self.threshold_percent_gap = threshold_percent_gap
        self.maximum_percentage_identity = maximum_percentage_identity
        self.maximum_number_sequences = maximum_number_sequences
        self.split_identity = split_identity
        self.multiple_msa = multiple_msa
        self.field = field
        self.dic_score = {}
        self.conservation_data = {}
        self.mafft_result = {}
        self.dict_r4s_cif_mapping = {}
        self.match_log = ""

        self.setup_logging()
        self.check_structure_file()
        if self.structure_file_type == ".pdb":
            prefix_structure_file = Path(structure_file).stem
            converted_structure_file = prefix_structure_file + ".cif"
            pdb2cif.convert2cif(self.structure_file, converted_structure_file)
            self.structure_file = converted_structure_file

        # Remove h2o and take the first model and first alternate location
        clean_cif = os.path.join(self.output_directory, "cleaned_structure.cif")
        util_cif.clean_cif(self.structure_file, self.structure_file, None, None)
        #self.structure_file = clean_cif

    def setup_logging(self):
        FORMAT = "[%(asctime)s] - %(levelname)s - %(message)s"
        DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
        logging.basicConfig(format=FORMAT, datefmt=DATE_FORMAT)
        self.logger = logging.getLogger("Conservation Pipeline")
        self.logger.setLevel(logging.INFO)

        fh_detailed = logging.FileHandler("rate4site.log")
        fh_detailed.setLevel(logging.INFO)
        fh_detailed.setFormatter(logging.Formatter(fmt=FORMAT, datefmt=DATE_FORMAT))
        self.logger.addHandler(fh_detailed)

        # Create a stream handler for logging to the console
        ch = logging.StreamHandler(stream=sys.stdout)
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter(fmt=FORMAT, datefmt=DATE_FORMAT))
        self.logger.addHandler(ch)


    def run(self):
        if self.multiple_msa:
            self.complex = True
        else:
            self.check_fasta()
        if self.complex:
            fasta_paths = ""
            # Filtering
            if not self.split_identity:
                paired = os.path.join(self.temporary_directory, "paired.fasta")
                self.remove_unpaired(self.msa_input, paired)
                self.run_hhfilter(paired, self.msa_hhfilter)
                self.get_best_match_from_fasta(self.msa_hhfilter, self.msa_filtered)
                self.get_subsequence(self.msa_filtered)
                fasta_paths = self.create_fasta(self.output_directory)
            else:
                if self.multiple_msa:
                    for msa in self.msa_input:
                        print(f"msa: {msa}")
                        self.make_fasta_seq_one_line(msa)
                    self.get_sequence_from_multi_fasta(self.msa_input)
                else:
                    self.make_fasta_seq_one_line(self.msa_input)
                    self.get_subsequence(self.msa_input)
                fasta_paths = self.create_fasta(self.output_directory)

                for index in self.fasta_sequences:
                    self.run_hhfilter(fasta_paths[index], fasta_paths[index])
                    self.get_best_match_from_fasta(fasta_paths[index], fasta_paths[index])

            for index in self.fasta_sequences:
                r4s_output = os.path.join(self.output_directory, "r4s_" + str(index))
                self.computeR4S(fasta_paths[index], r4s_output)
                self.get_score_from_r4s(r4s_output + ".grade", index)
                self.compute_conservation(index)
            self.assign_sequence_cif_fasta()
            for group in self.link_sequence_chain:
                self.align_pdb_r4s(group["index"], group["occurrence"], group["chain"])
            structure_output = os.path.join(self.output_directory, "structure_output.cif")
            self.create_cif_converge_diverge_r4s(structure_output)
        else:
            self.run_hhfilter(self.msa_input, self.msa_hhfilter)
            self.get_best_match_from_fasta(self.msa_hhfilter, self.msa_filtered)
            r4s_output = os.path.join(self.output_directory, "r4s")
            self.computeR4S(self.msa_filtered, r4s_output)
            self.get_score_from_r4s(r4s_output + ".grade", 0)
            self.compute_conservation(0)
            self.assign_sequence_cif_fasta()
            for group in self.link_sequence_chain:
                self.align_pdb_r4s(group["index"], group["occurrence"], group["chain"])
            structure_output = os.path.join(self.output_directory, "structure_output.cif")
            self.create_cif_converge_diverge_r4s(structure_output)

        with open(os.path.join(self.output_directory, "match_log.txt"), "w") as f:
            f.write(self.match_log)

    def run_hhfilter(self, input, output):
        os.system(f"{HHFILTER} -i {input} -o {output} -id {self.maximum_percentage_identity}")
        return output

    def parse_fasta(self, file):
        """Parse the FASTA file and return a list of dictionaries with headers and sequences."""
        sequences = []
        with open(file, "r") as f:
            header = ""
            sequence = ""
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header and sequence:
                        sequences.append({"header": header, "sequence": sequence})
                    header = line
                    sequence = ""
                else:
                    sequence += line
            if header and sequence:
                sequences.append({"header": header, "sequence": sequence})
        return sequences

    def sort_by_best_match(self, sequences):
        """Get the best match compared to the first sequence in the list."""
        if not sequences:
            return None

        main_header = sequences[0]["header"]
        main_sequence = sequences[0]["sequence"]

        for i in range(1, len(sequences)):
            other_sequence = sequences[i]["sequence"]
            match_count = sum(1 for a, b in zip(main_sequence, other_sequence) if a == b)
            sequences[i]["matches"] = match_count
        sorted_by_match = sorted(sequences[1:], reverse=True, key=lambda x: x["matches"])

        sorted_by_match.insert(0, {"header": main_header, "sequence": main_sequence})

        return sorted_by_match

    def write_fasta(self, sequences, output_file):
        with open(output_file, "w") as f:
            for sequence in sequences:
                f.write(sequence["header"] + "\n")
                f.write(sequence["sequence"] + "\n")

    def get_best_match_from_fasta(self, fasta, output):
        """get the closest sequences from the first sequence of a fasta and write them in a file"""
        sequences = self.parse_fasta(fasta)
        sorted_sequences = self.sort_by_best_match(sequences)
        best_sequences = sorted_sequences[: self.maximum_number_sequences]
        self.write_fasta(best_sequences, output)

    def remove_unpaired(self, fasta, output_file):
        """remove unpaired sequences from a fasta which consist of sequences after we see again the header of the main sequence"""
        sequences = self.parse_fasta(fasta)
        main_header = sequences[0]["header"]
        main_header = main_header.split("\t")[0]
        print(main_header)
        for index, sequence in enumerate(sequences[1:]):
            if main_header == sequence["header"]:
                print(f"sequences = {index}")
                self.write_fasta(sequences[:index], output_file)
                break

    def check_structure_file(self):
        """Check if structure_file is a pdb or a cif and if not tell the user about it and stop the pipeline"""
        valid_extensions = [".pdb", ".cif"]
        if Path(self.structure_file).is_file() and Path(self.structure_file).suffix in valid_extensions:
            self.structure_file_type = Path(self.structure_file).suffix
        else:
            self.logger.error("wrong structure file type, should be either .pdb or .cif")
            sys.exit(1)

    def check_fasta(self):
        """Detect if the fasta has multiple protein in it and store these information (length of these different sequences and their number)"""
        self.complex = False
        with open(self.msa_input, "r") as f:
            firstline = f.readline()
            if firstline.startswith("#") and re.search('#[\s]*[0-9,]+',firstline):
                self.logger.info("Fasta with multiple sequences")
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
                    self.length_number_couple.append((length, number_of_sequences[index]))
            else:
                self.logger.info("Fasta with only one sequence")

    def make_fasta_seq_one_line(self, fasta):
        with open(fasta, "r+") as f:
            lines = f.readlines()
            one_line_file = ""
            sequence = ""
            header = ""

            for line in lines:
                if not line.startswith("#"):
                    if line.startswith(">"):
                        if header:
                            one_line_file += header + sequence + "\n"
                        header = line
                        sequence = ""
                    else:
                        sequence += line.strip()
            one_line_file += header + sequence
            f.seek(0)
            f.truncate()
            f.write(one_line_file)

    def get_subsequence(self, fasta):
        """get all the subsequence from a fasta"""
        self.fasta_sequences = {}

        with open(fasta, "r") as f:
            lines = f.readlines()
            current_start_position = 0

            for index, tup in enumerate(self.length_number_couple):
                count_sequences = 0
                length = current_start_position + int(tup[0])
                self.fasta_sequences[index] = {}
                self.fasta_sequences[index]["count_sequences"] = []
                for line in lines:
                    stripped_line = line.strip()
                    sequence = ""
                    if not line.startswith("#"):
                        if line.startswith(">"):
                            count_sequences += 1
                            header = stripped_line
                            self.fasta_sequences[index][count_sequences] = ""

                        else:
                            sequence = stripped_line[current_start_position:length]
                            self.fasta_sequences[index]["count_sequences"].append(count_sequences)
                            self.fasta_sequences[index][count_sequences] = {}
                            self.fasta_sequences[index][count_sequences]["sequence"] = sequence
                            self.fasta_sequences[index][count_sequences]["occurrence"] = tup[1]
                            self.fasta_sequences[index][count_sequences]["header"] = header
                current_start_position = length

    def get_sequence_from_multi_fasta(self, array_fasta):
        """get all the fasta sequences from multiple fasta"""
        self.fasta_sequences = {}
        # dont forget to convert to one liner before
        for index, fasta in enumerate(array_fasta):
            count_sequences = 0
            with open(fasta, "r") as f:
                lines = f.readlines()
                count_sequences = 0
                self.fasta_sequences[index] = {}
                self.fasta_sequences[index]["count_sequences"] = []
                for line in lines:
                    stripped_line = line.strip()
                    sequence = ""
                    if line.startswith(">"):
                        count_sequences += 1
                        header = stripped_line
                        self.fasta_sequences[index][count_sequences] = ""

                    else:
                        sequence = stripped_line
                        # occurrence set to 1 by default
                        occurrence = 1
                        self.fasta_sequences[index]["count_sequences"].append(count_sequences)
                        self.fasta_sequences[index][count_sequences] = {}
                        self.fasta_sequences[index][count_sequences]["sequence"] = sequence
                        self.fasta_sequences[index][count_sequences]["occurrence"] = occurrence
                        self.fasta_sequences[index][count_sequences]["header"] = header

    def create_fasta(self, output_directory):
        """Create a fasta with the subsequence of each protein contained in the original fasta,
        return an array of the output path
        """
        count = 1
        array_path = []
        for index in self.fasta_sequences:
            fasta = f"split_{count}.fasta"
            fasta_output = os.path.join(output_directory, fasta)
            self.fasta_sequences[index]["fasta"] = fasta
            with open(fasta_output, "w") as f:
                for count_sequences in self.fasta_sequences[index]["count_sequences"]:
                    header = self.fasta_sequences[index][count_sequences]["header"]
                    sequence = self.fasta_sequences[index][count_sequences]["sequence"]
                    f.write(f"{header}\n{sequence}\n")
            count += 1
            array_path.append(fasta_output)
        return array_path

    def computeR4S(self, msa_file, result):
        command = f"{RATE4SITE}  -s {msa_file} -o {result}.grade -x {result}.ph"
        self.logger.info("Rate4site command :")
        self.logger.info(command)
        os.system(command)

    def get_score_from_r4s(self, input, index):
        """Compute the score for each, aa, the percentage of no gap"""

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
            first_seq = self.fasta_sequences[index]["count_sequences"][0]
            self.dic_score[index]["occurrence"] = self.fasta_sequences[index][first_seq]["occurrence"]
        else:
            self.dic_score[index]["occurrence"] = 1

    def compute_conservation(self, index=0):
        """Compute the conservation based on the values obtained in rate4site"""
        conservation_data = {}
        conservation_data["diverged"] = {}
        conservation_data["conserved"] = {}
        conservation_data["percentage_no_gap"] = self.dic_score[index]["percentage_no_gap"]
        conservation_data["no_gap"] = {}
        conservation_data["order"] = {}
        conservation_data["bin"] = {}
        scores = self.dic_score[index]["score"].values()
        amplitude = max(scores) - min(scores)
        conservation_mini = min(scores)

        for res in self.dic_score[index]["res_order"]:
            MAX_SCORE = 99.0
            MIN_SCORE = 1.0

            conservation = MAX_SCORE - (MAX_SCORE * (self.dic_score[index]["score"][res] - conservation_mini) / amplitude)

            if conservation > MAX_SCORE:
                conservation_data["bin"][res] = MAX_SCORE
            elif conservation < MIN_SCORE:
                conservation_data["bin"][res] = MIN_SCORE
            else:
                conservation_data["bin"][res] = conservation

            if conservation < self.diverged:
                conservation_data["diverged"][res] = True
            else:
                conservation_data["diverged"][res] = False

            if conservation > self.conserved:
                conservation_data["conserved"][res] = True
            else:
                conservation_data["conserved"][res] = False

            if conservation_data["percentage_no_gap"][res] > self.threshold_percent_gap:
                conservation_data["no_gap"][res] = True
            else:
                conservation_data["no_gap"][res] = False

        self.conservation_data[index] = conservation_data

    def extract_sequence_cif(self, chains=[]):
        """Extract sequence from the cif"""
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
            if chains:
                if name_chain not in chains:
                    continue

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
        pdb_sequences = self.extract_sequence_cif(self.array_chains)
        if self.array_chains: self.logger.info(f"User-specified chains: {self.array_chains}")
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
                self.mafft_result[index][chain]["index_cif"] = pdb_sequences[chain]["order"]
                self.mafft_result[index][chain]["sequences"] = self.get_sequence_from_alignment(output_mafft.name)

        self.compare_mafft_alignment()

    def compare_mafft_alignment(self):
        """Compare and assign the best alignment between r4s and cif (if percentage identity over a certain threshold)"""
        MINIMUM_PERCENTAGE_IDENTITY = 75
        pdb_sequences = self.extract_sequence_cif(self.array_chains)
        chains = pdb_sequences["chain"]
        self.link_sequence_chain = []
        self.match_log = ""

        all_pairs = []
        for index in self.mafft_result:
            for occurrence in range(int(self.dic_score[index]["occurrence"])):
                percentage = {
                    "max": 0,
                    "index": "",
                    "chain": "",
                }
                for chain in chains:
                    header_r4s = self.mafft_result[index][chain]["sequences"]["order"][0]
                    header_cif = self.mafft_result[index][chain]["sequences"]["order"][1]
                    # print(self.mafft_result[index][chain]["sequences"]["order"][0])
                    percentage_identity = self.compare_sequences(
                        self.mafft_result[index][chain]["sequences"][header_r4s],
                        self.mafft_result[index][chain]["sequences"][header_cif],
                    )
                    all_pairs.append([index,occurrence,chain,percentage_identity])
        
        # e.g. all_pairs = [[0, 0, 'A', 76.47058823529412], [1, 0, 'A', 64.85148514851485]]
        seen = {'msa':set(),'chain':set()}
        for index, occurrence, chain, percentage_identity in sorted(all_pairs,key=lambda x: x[-1],reverse=True):
            if percentage_identity < MINIMUM_PERCENTAGE_IDENTITY:
                break
            if (index,occurrence) not in seen['msa'] and chain not in seen['chain']:
                percentage = {'max': percentage_identity, 'index': index, 'chain': chain, 'occurrence': occurrence}
                header_r4s = self.mafft_result[index][chain]["sequences"]["order"][0]
                header_cif = self.mafft_result[index][chain]["sequences"]["order"][1]

                self.logger.info(f"Match found between MSA sequence {index+1} (occurrence {occurrence+1}/{int(self.dic_score[index]['occurrence'])}) and protein chain {chain} in structure with {round(percentage_identity,0)}% identity")
                self.logger.info(f"> MSA {index+1} (occurrence {occurrence+1}/{int(self.dic_score[index]['occurrence'])})")
                self.logger.info(self.mafft_result[index][chain]["sequences"][header_r4s])
                self.logger.info(f"> protein chain {chain} in structure")
                self.logger.info(self.mafft_result[index][chain]["sequences"][header_cif])

                self.match_log += f"Match found between MSA sequence {index+1} (occurrence {occurrence+1}/{int(self.dic_score[index]['occurrence'])}) and protein chain {chain} in structure with {round(percentage_identity,0)}% identity\n"
                self.match_log += f"> MSA {index+1} (occurrence {occurrence+1}/{int(self.dic_score[index]['occurrence'])})\n"
                self.match_log += self.mafft_result[index][chain]["sequences"][header_r4s] + '\n'
                self.match_log += f"> protein chain {chain} in structure\n"
                self.match_log += self.mafft_result[index][chain]["sequences"][header_cif] + '\n'
                self.match_log += "\n"

                self.link_sequence_chain.append(percentage)
                seen['msa'].add((index,occurrence))
                seen['chain'].add(chain)

        for index in self.mafft_result:
            for occurrence in range(int(self.dic_score[index]["occurrence"])):
                if (index,occurrence) not in seen['msa']:
                    self.logger.info(f"No match found for MSA sequence {index+1} (occurrence {occurrence+1}/{int(self.dic_score[index]['occurrence'])})")
                    self.logger.info(f"> MSA {index+1} (occurrence {occurrence+1}/{int(self.dic_score[index]['occurrence'])})")
                    self.logger.info(self.mafft_result[index][chains[0]]["sequences"][self.mafft_result[index][chains[0]]["sequences"]["order"][0]].replace("-",""))

                    self.match_log += f"No match found for MSA sequence {index+1} (occurrence {occurrence+1}/{int(self.dic_score[index]['occurrence'])})\n"
                    self.match_log += f"> MSA {index+1} (occurrence {occurrence+1}/{int(self.dic_score[index]['occurrence'])})\n"
                    self.match_log += self.mafft_result[index][chains[0]]["sequences"][self.mafft_result[index][chains[0]]["sequences"]["order"][0]].replace("-","") + "\n"
                    self.match_log += "\n"

        index = list(self.mafft_result.keys())[0]
        for chain in chains:
            if chain not in seen['chain']:
                self.logger.info(f"No match found for protein chain {chain} in given structure")
                self.logger.info(f"> protein chain {chain}")
                self.logger.info(self.mafft_result[index][chain]["sequences"][self.mafft_result[index][chain]["sequences"]["order"][1]].replace("-",""))

                self.match_log += f"No match found for protein chain {chain} in given structure\n"
                self.match_log += f"> protein chain {chain}\n"
                self.match_log += self.mafft_result[index][chain]["sequences"][self.mafft_result[index][chain]["sequences"]["order"][1]].replace("-","") + "\n"
                self.match_log += "\n"

        # e.g. self.link_sequence_chain = [{'max': 76.47058823529412, 'index': 0, 'chain': 'A', 'occurence': 0}]

    def compare_sequences(self, seq1, seq2):
        """Compare two sequences and return the percentage of identity between both"""
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of the same length")

        matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != "-")

        percentage_similarity = (matches / len(seq1.replace("-", ""))) * 100

        return percentage_similarity

    def get_sequence_from_alignment(self, alignment):
        """Return all the sequences from a fasta"""
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
        """Align the r4s and cif sequence and if match assign the value of the r4s"""
        header_r4s = self.mafft_result[index][chain]["sequences"]["order"][0]
        header_cif = self.mafft_result[index][chain]["sequences"]["order"][1]

        seq_r4s = self.mafft_result[index][chain]["sequences"][header_r4s]
        seq_cif = self.mafft_result[index][chain]["sequences"][header_cif]
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
                    r4s_value_correspondence["bin"][cur_index_cif] = self.conservation_data[index]["bin"][str(ite_r4s)]
                    r4s_value_correspondence["conserved"][cur_index_cif] = self.conservation_data[index]["conserved"][
                        str(ite_r4s)
                    ]
                    r4s_value_correspondence["diverged"][cur_index_cif] = self.conservation_data[index]["diverged"][
                        str(ite_r4s)
                    ]
                    r4s_value_correspondence["percentage_no_gap"][cur_index_cif] = self.conservation_data[index][
                        "percentage_no_gap"
                    ][str(ite_r4s)]
                    r4s_value_correspondence["no_gap"][cur_index_cif] = self.conservation_data[index]["no_gap"][str(ite_r4s)]
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
        if index not in self.dict_r4s_cif_mapping:
            self.dict_r4s_cif_mapping[index] = {}
        self.dict_r4s_cif_mapping[index][occurrence] = {}
        self.dict_r4s_cif_mapping[index][occurrence][chain] = r4s_value_correspondence

    def rna_cleaning(self, cif):
        """Biopython cif parser add '' to rna component, this function trim that"""
        search_text = r"'(\w+)''"
        with open(cif, "r+") as f:
            file = f.read()
            cleaned_file = re.sub(search_text, lambda m: f"{m.group(1)}'", file)
            f.seek(0)
            f.truncate()
            f.write(cleaned_file)

    def add_data(self, cif):
        """Because PyMOL doesn't accept cif files without the data_ information"""
        with open(cif, "r+") as f:
            file = f.read()
            f.seek(0)
            f.truncate()
            f.write("data_output\n" + file)

    def create_cif_converge_diverge_r4s(self, output):
        """Create the a cif with the with setting the value of r4s in occupancy/bfactor column"""
        for residue in self.structure.get_residues():
            for atom in residue.get_atoms():
                if self.field == "occupancy": atom.set_occupancy(0)
                if self.field == "bfactor": atom.set_bfactor(0)

        for index in self.dict_r4s_cif_mapping:
            for occurrence in self.dict_r4s_cif_mapping[index]:
                for chain in self.dict_r4s_cif_mapping[index][occurrence]:
                    struct_chain = self.structure[0][chain]
                    for residue in struct_chain.get_residues():
                        residue_number = residue.get_id()[1]
                        if residue_number in self.dict_r4s_cif_mapping[index][occurrence][chain]["bin"]:
                            for atom in residue.get_atoms():
                                if self.field == "occupancy":
                                    atom.set_occupancy(
                                        self.dict_r4s_cif_mapping[index][occurrence][chain]["bin"][residue_number]
                                    )
                                if self.field == "bfactor":
                                    atom.set_bfactor(
                                        self.dict_r4s_cif_mapping[index][occurrence][chain]["bin"][residue_number]
                                    )

        io = MMCIFIO()
        io.set_structure(self.structure)
        io.save(output)
        temp_dic = io.dic
        # Change the first column who can contain different name chain which ngl have problem interpreting
        if "_atom_site.label_asym_id" in io.dic:
            temp_dic["_atom_site.label_asym_id"] = temp_dic["_atom_site.auth_asym_id"]
            io.set_dict(temp_dic)
            io.save(output)
        self.rna_cleaning(output)
        self.add_data(output)
