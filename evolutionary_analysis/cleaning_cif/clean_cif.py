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
    """Remove residues with a specific name from the model."""
    for chain in model:
        residues_to_remove = [residue.id for residue in chain if residue.get_resname() == residue_name]
        for residue_id in residues_to_remove:
            chain.detach_child(residue_id)


def choose_alternative_atom(model, alternative=None):
    """only keep the atoms with the specified alternative"""
    if not alternative:
        for atom in model.get_atoms():
            if atom.is_disordered():
                altlocs = atom.disordered_get_id_list()
                for index in range(1, len(altlocs)):
                    atom.disordered_remove(altlocs[index])
    else:
        for atom in model.get_atoms():
            if atom.is_disordered():
                altlocs = atom.disordered_get_id_list()
                for alt in altlocs:
                    if alt != alternative:
                        atom.disordered_remove(alt)


def write_cif(model, output_file):
    io = MMCIFIO()
    io.set_structure(model)
    io.save(output_file)


def clean_cif(cif, output_cif, model=None, alt_label=None):
    """create a cif without H2O, keep only one alt label and one model"""
    structure = parse_cif(cif)
    model = get_model(structure, model)

    remove_residue(model, "HOH")

    choose_alternative_atom(model, alt_label)

    write_cif(model, output_cif)
