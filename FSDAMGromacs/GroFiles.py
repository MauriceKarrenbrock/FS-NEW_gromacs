# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Functions to change or extract gro files
"""


def add_atom_to_gro_file(input_gro_file,
                         output_gro_file,
                         coordinates,
                         velocities=None,
                         atom_name='DU',
                         atom_residue_name='DUM'):
    """adds a given atom to a gro file

    comes in handy to add a dummy atom

    Parameters
    -------------
    input_gro_file : str
        the name (or path) of the input gro file
    output_gro_file : str
        the name (or path) of the input gro file
        can be the same as the input one
    coordinates : iterable of float
        (x,y,z) iterable of any type (list, tuple, ...)
    velocities : iterable of float, optional
        (vx,vy,vz) iterable of any type (list, tuple, ...)
        for default velocities will be left blank
    atom_name : str
        2 characters atom name
    atom_residue_name :
        2 or 3 characters residue name
    """

    with open(input_gro_file, 'r') as f:
        input_lines = f.readlines()

    #adds an atom to the atom count
    input_lines[1] = f'{int(input_lines[1].strip()) + 1}\n'

    for i in range(len(input_lines) - 1, 0, -1):

        if input_lines[i].strip() != '':

            if velocities is not None:

                velocity_string = '{:f8.4}{:f8.4}{:f8.4}'.format(
                    velocities[0], velocities[1], velocities[2])

            else:

                velocity_string = 24 * ' '

            input_lines[
                i -
                1] += '{:>5}{:>5}{:>5}{:>5}{:f8.3}{:f8.3}{:f8.3}{}\n'.format(
                    int(input_lines[i - 1][0:5].strip()) + 1,  #residue number
                    atom_residue_name,  #residue name
                    atom_name,  #atom name
                    int(input_lines[i - 1][15:20].strip()) + 1,  #atom number
                    coordinates[0],
                    coordinates[1],
                    coordinates[2],
                    velocity_string)

            break

    with open(output_gro_file, 'w') as f:
        f.writelines(input_lines)
