# Output form MSA-tools evolutionary analysis

https://bioi2.i2bc.paris-saclay.fr/django/msa-tools

This folder contains:

## Structure & visualisation
- structure_output.cif: input protein structure in cif format with calculated normalised evolutionary rates in the occupancy column
- scripts/launch_chimerax.cxc: ChimeraX script to load structure and colour by evolutionary rate (click right + Open with ChimeraX)
- scripts/launch_pymol.pml: PyMOL script to load structure and colour by evolutionary rate (in Windows: click right + Open with PyMOL, in Linux: open PyMOL then Click on File > Run script and select this file)
- scripts/spectrumany.py: script needed for launch_pymol.pml in PyMOL

## MSA
- structure.a3m: MSA as given by the user or outputted by MMseqs2
- split_?.fasta: filtered and reduced per-chain MSA as inputted to Rate4Site

## Rate4Site scores
- r4s_?.grade: raw Rate4Site output (with scores per residue and confidence measurements)
