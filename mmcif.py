from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO

test = MMCIFParser().get_structure(
    structure_id="cool",
    filename="./Eole3_6VEH-bGFPc_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_59477.cif",
)
test_dic = MMCIF2Dict(
    "./Eole3_6VEH-bGFPc_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_59477.cif"
)


def extract_seq(mmcif_structure):
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
    for chain in mmcif_structure.get_chains():
        name_chain = chain.get_id()
        dic_sequence["chain"] = name_chain
        dic_sequence[chain] = {}
        dic_sequence[chain]["order"] = []
        sequence = ""
        for residue in chain.get_residues():
            # print(dir(residue))
            residue_number = residue.get_id()[1]
            name_res = residue.get_resname()
            if name_res in one_letter:
                amino_acid = one_letter[name_res]
            else:
                amino_acid = "X"
            sequence += amino_acid
            dic_sequence[chain][residue_number] = amino_acid
            dic_sequence[chain]["order"].append(residue_number)
        dic_sequence[chain]["sequence"] = sequence


extract_seq(test)

io = MMCIFIO()
io.set_structure(test)
io.save("result/test_parser.cif")
