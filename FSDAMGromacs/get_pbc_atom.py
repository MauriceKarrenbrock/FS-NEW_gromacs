# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Functions to get the pbc atom
"""

import sys

import MDAnalysis as mda
import numpy as np
import PythonPDBStructures.geometry as _geo
import PythonPDBStructures.pdb.add_chain_id as _chain_id
import PythonPDBStructures.pdb.biopython_utils as _bio


def get_protein_pbc_atom(structure_file,
                         mdanalysis_selection_string='protein'):
    """Get the number of the nearest atom to the geometrical center of the Protein

    it uses `MDAnalysis` and `Biopython`

    Parameters
    ------------
    structure_file : str or Path
        the file (gro, pdb, ...) contining the protein and other stuff
    mdanalysis_selection_string : str, optional, default='protein'
        the sring to use for `MDAnalysis.Universe.select_atoms` method if you want the pbc atom
        of a protein the default is perfecto, otherwise check MDAnalysis documentation online

    Returns
    ----------
    int
        the atom number
    """

    #for now MDAnalysis doesn't deal well with pathlib objects
    # issue 2497 n github https://github.com/MDAnalysis/mdanalysis/issues/2497
    if not isinstance(structure_file, str):

        structure_file = str(structure_file)

    tmp_file = 'TMP_for_pbc_atom.pdb'

    universe = mda.Universe(structure_file)

    atoms = universe.select_atoms(mdanalysis_selection_string)

    atoms.write(tmp_file)

    _chain_id.add_chain_id_pdb(tmp_file, 'A')

    struct = _bio.parse_pdb('idid', tmp_file)

    GEO_CENTER = _geo.get_center_of_mass(struct, geometric=True)

    atoms = struct.get_atoms()

    min_dist = sys.float_info.max

    output_atom = None

    for atom in atoms:

        dist = (atom.coord - GEO_CENTER)**2
        dist = np.sum(dist)
        dist = dist**0.5

        if dist < min_dist:

            min_dist = dist

            output_atom = atom.get_serial_number()

    if output_atom is None:

        raise Exception('Something went very wrong')

    return output_atom
