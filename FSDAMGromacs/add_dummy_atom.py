# -*- coding: utf-8 -*-
# pylint: disable=duplicate-code
# pylint: disable=too-many-lines
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""high level function to add a dummy atom to a gromacs topology and gro file
"""

import mdtraj
import numpy as np
import PythonAuxiliaryFunctions.path as _path

import FSDAMGromacs.itp_files as _itp_files
import FSDAMGromacs.top_files as _top_files


def add_dummy_atom_to_topology(top_file,
                               residue_name='DUM',
                               mass=1.0e20,
                               charge=0.0,
                               sigma=0,
                               epsilon=0):
    """a high level function that adds an atom (and it's residue) to a top file

    will overwrite `top_file`

    Parameters
    -------------
    top_file : str
        the topology file (.top)
    residue_name : str of 3 characters, default=DUM
    mass : float, default=1.0e20
    charge : float, default=0.0
    sigma : float, default=0
        LJ parameter
    epsilon : float, default=0
        LJ parameter
    """

    atom_types, itp_file = _itp_files.create_dummy_atom_itp(mass=mass,
                                                            charge=charge,
                                                            sigma=sigma,
                                                            epsilon=epsilon,
                                                            name=residue_name,
                                                            charge_group=1,
                                                            atom_number=1)

    with open(f'{residue_name}_atomtypes.itp', 'w') as f:

        f.writelines(atom_types)

    with open(f'{residue_name}.itp', 'w') as f:

        f.writelines(itp_file)

    atom_types_path = _path.absolute_filepath(f'{residue_name}_atomtypes.itp')
    itp_file_path = _path.absolute_filepath(f'{residue_name}.itp')

    _top_files.add_include_after_FF(include_line=atom_types_path,
                                    input_top_file=top_file,
                                    output_top_file=top_file)

    _top_files.add_include_after_atomtypes(include_line=itp_file_path,
                                           input_top_file=top_file,
                                           output_top_file=top_file)

    _top_files.add_molecules(name=residue_name,
                             number=1,
                             input_top_file=top_file,
                             output_top_file=top_file)


def add_dummy_atom_to_center_of_mass(structure_file,
                                     select=None,
                                     residue_name='DUM',
                                     atom_name='DU',
                                     xyz=None):
    """Will add a dummy atom (and it's residue) to a the center of mass

    will overwrite the input structure

    Parameters
    -----------
    structure_file : str
        a structure file supported by mdtrj (pdb, gro, ...)
    select : str, default=None
        a mdtraj selection string on which the center of mass will
        be calculated
        if left None the center of mass will be calculated on the
        whole structure
    residue_name : str, default='DUM'
    atom_name : str, default='DU'
    xyz : iterable(float), optional, default=None
        you can choose to give custom coordinates for the dummy atom
        in this case no center of mass will be calculated
    """

    trj = mdtraj.load(structure_file)
    top = trj.top

    if xyz is None:

        com = mdtraj.geometry.distance.compute_center_of_mass(trj,
                                                              select=select)

    else:

        com = np.array([xyz])

    dummy_element = mdtraj.core.element.Element(0, 'dummy_atom', atom_name,
                                                0.0, 0.0)

    dummy_chain = top.add_chain()

    dummy_residue = top.add_residue(residue_name, dummy_chain)

    top.add_atom(atom_name, dummy_element, dummy_residue)

    trj.xyz = np.hstack((trj.xyz, com))

    trj.save(structure_file)
