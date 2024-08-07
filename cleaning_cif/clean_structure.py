from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO


def parse_cif(cif):
    """return the structures of a cif"""
    parser = MMCIFParser()
    structure = parser.get_structure(structure_id="", filename=cif)

    return structure


def get_model(structure, model=None):
    """Get the specified model from the structure. If no model is specified, return the first model."""
    if not model:
        return structure[0]
    else:
        return structure[model]


def remove_residue(model, residue_name="HOH"):
    for chain in model:
        for residue in chain:
            if residue.get_resname() == residue_name:
                print(residue.get_id())
                chain.detach_child(residue.id)
    return model


def write_cif():
    pass


def clean_cif(cif, model=None, alt_label=None):
    """create a cif without H2O, keep only one alt label and one model"""
    structure = parse_cif(cif)
    model = get_model(structure, model)

    model = remove_residue(model)

    io = MMCIFIO()
    io.set_structure(model)

    io.save("cleaned_structure.cif")


def test(cif_file):
    parser = MMCIFParser()
    structure = parser.get_structure("", cif_file)

    # Get the model (assume single model for simplicity)
    model = structure
    cif_dict = parser._mmcif_dict

    # Iterate over chains and residues to get atoms
    for model in structure:
        # print(chain)
        for residue in model:
            for atom in residue:
                # Print label_alt_id
                label_alt_id = atom.get_full_id()  # Extracting label_alt_id
                # print(label_alt_id)
                pdbx_pdb_model_num = atom.get_parent().get_parent()  # Extracting pdbx_PDB_model_num
                # print(pdbx_pdb_model_num)


clean_cif("./1ssx.cif")
