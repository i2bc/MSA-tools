from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB import Select


# https://biopython.org/wiki/Remove_PDB_disordered_atoms
class NotDisordered(Select):

    def __init__(self, alt_label="A"):
        self.alt_label = alt_label

    def accept_atom(self, atom):
        if not atom.is_disordered() or atom.get_altloc() == self.alt_label:
            atom.set_altloc(" ")  # Eliminate alt location ID before output.
            return True
        else:
            return False


def parse_cif(cif):
    """return the structures of a cif"""
    parser = MMCIFParser()
    structure = parser.get_structure(structure_id="", filename=cif)

    return structure


def get_model(structure, model=None):
    """Get the specified model from the structure. If no model is specified, return the first model.
    The model in biopython use an array and not the actual name of the model so you need to specify the model at the position he appear.
    """
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


def choose_alternative_atom(model):
    """find the alternative to keep as it might not always be A, B etc."""
    for atom in model.get_atoms():
        if atom.is_disordered():
            return atom.get_altloc()
    return None

def write_cif(model, output_file, alt_label):
    io = MMCIFIO()
    io.set_structure(model)
    io.save(output_file, select=NotDisordered(alt_label))


def clean_cif(cif, output_cif, model=None, alt_label=None):
    """create a cif without H2O, keep only one alt label and one model"""
    structure = parse_cif(cif)
    model = get_model(structure, model)

    remove_residue(model, "HOH")

    if not alt_label:
        alt_label = choose_alternative_atom(model)

    write_cif(model, output_cif, alt_label)
