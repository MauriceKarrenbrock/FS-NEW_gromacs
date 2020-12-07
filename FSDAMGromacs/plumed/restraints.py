# -*- coding: utf-8 -*-
# pylint: disable-msg=too-many-locals
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Functions to introduce restraints with plumed
"""

import PythonAuxiliaryFunctions.files_IO.write_file as write_file


def COM_COM_restraint(atom_groups,
                      restraint_parameters,
                      plumed_file='plumed.dat',
                      distances_file='distances.out',
                      stride=100,
                      geometric=False):
    """create a center of mass (COM) - center of mass restraint with plumed

    creates a plumed input file https://www.plumed.org in order to create
    a restraint between 2 centers of mass

    all distances (input and output) will be in nanometers nm!

    Parameters
    ------------
    atom_groups : dict
        {"group_name" : [<list of atoms>], ...} the keys of the dictionary must
        be the group names (NO SPACES!!!)
        and the value must be a list containing all atom numbers of the group,
        there can be as many groups as
        you like
    restraint_parameters : list
        a nested list: [ ["grop_name_1", "grop_name_2", equilibium_distance_nm,
        harmonic_kappa, linear_slope], ... ]
        the names must be str the other parameters float, remember that you
        can always set a certain parameter to zero 0.
    plumed_file : str
        the name of the plumed input file to be created
        (default 'plumed.dat')
    distances_file : str
        the name of the file that plumed will create and in which all the distances will be
        written out (default 'distances.out')
    stride : int
        `distances_file` will be updated each `stride` MD steps
        (default 100)
    geometric : bool, optional
        if True instead of the center of mass the geometrical center
        will be calulated (default False). Can come in handy with strange
        atoms/residues

    Raises
    ---------
    ValueError
        if `atom_groups` contains only one group
        but there is no check to see if
        `restraint_parameters` doesn't contain a value for
        each possible couple
    """
    def make_atoms_string(atoms):
        """private
        """

        #build atoms string
        atom_string = ''
        for atom in atoms:

            atom_string += f'{atom},'

        #remove last comma
        atom_string = atom_string[:-1]

        return atom_string

    if len(atom_groups) < 2:
        raise ValueError(
            f'need at least 2 atom groups, not {len(atom_groups)}')

    output = []

    #put the units explicitly in order to be always sure what you
    #get in output and what you are expected to put in input
    output.append('UNITS LENGTH=nm TIME=ps \n\n')

    if geometric:
        com = 'CENTER'
    else:
        com = 'COM'

    #WHOLEMOLECULES string
    string = 'WHOLEMOLECULES '
    for i, name in enumerate(atom_groups.keys()):

        #build atoms string
        atom_string = make_atoms_string(atom_groups[name])

        string += f'ENTITY{i}={atom_string} '

    string += '\n'

    output.append(string)

    #define centers of mass (COM or CENTER):
    for name in atom_groups.keys():

        atom_string = make_atoms_string(atom_groups[name])

        string = f'{name}: {com} ATOMS={atom_string}\n'

        output.append(string)

    #define COM-COM distances DISTANCE NOPBC
    names = list(atom_groups.keys())
    distances_list = []
    for i in range(len(names)):

        for j in range(i + 1, len(names)):

            string = f'{names[i]}_{names[j]}_dist: DISTANCE ATOMS={names[i]},{names[j]} NOPBC\n'

            distances_list.append(
                [f'{names[i]}_{names[j]}_dist', names[i], names[j]])

            output.append(string)

    #create restraints RESTRAINT
    for couple in restraint_parameters:

        for distance in distances_list:

            if (distance[1]
                    in (couple[0], couple[1])) and (distance[2]
                                                    in (couple[0], couple[1])):

                string = \
                'RESTRAINT ARG={} AT={:.4f} KAPPA={:.4f} SLOPE={:.4f}\n'.format(
                    distance[0],
                    couple[2],
                    couple[3],
                    couple[4]
                )

                output.append(string)

    #PRINT statement
    string = 'PRINT ARG='

    for distance in distances_list:

        string += distance[0] + ','

    #remove last comma
    string = string[:-1]

    string += f' STRIDE={stride} '

    string += f'FILE={distances_file}\n'

    output.append(string)

    #write file
    write_file.write_file(output, plumed_file)
