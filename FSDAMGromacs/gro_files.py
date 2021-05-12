# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Functions to change or extract gro files
"""

import PythonAuxiliaryFunctions.files_IO.read_file as read_file
import PythonAuxiliaryFunctions.files_IO.write_file as write_file


def add_atom_to_gro_file(input_gro_file,
                         output_gro_file,
                         coordinates,
                         velocities=None,
                         atom_name='DU',
                         atom_residue_name='DUM'):
    """adds a given atom to a gro file

    comes in handy to add a dummy atom
    remember that unlike PDB files GRO files
    are in nanometers nm!!

    Parameters
    -------------
    input_gro_file : str
        the name (or path) of the input gro file
    output_gro_file : str
        the name (or path) of the input gro file
        can be the same as the input one
    coordinates : iterable of float
        (x,y,z) iterable of any type (list, tuple, ...)
        remember that unlike PDB files GRO files
        are in nanometers nm!!
    velocities : iterable of float, optional
        (vx,vy,vz) iterable of any type (list, tuple, ...)
        for default velocities will be left blank
    atom_name : str
        2 characters atom name
    atom_residue_name :
        2 or 3 characters residue name
    """

    input_lines = read_file.read_file(input_gro_file)

    #adds an atom to the atom count
    input_lines[1] = '{:>5}\n'.format(int(input_lines[1].strip()) + 1)

    for i in range(len(input_lines) - 1, 0, -1):

        if input_lines[i].strip() != '':

            if velocities is not None:

                velocity_string = '{:8.4f}{:8.4f}{:8.4f}'.format(
                    velocities[0], velocities[1], velocities[2])

            else:

                velocity_string = 24 * ' '

            input_lines[i - 1] \
                += '{:>5}{:<5}{:>5}{:>5}{:8.3f}{:8.3f}{:8.3f}{}\n'.format(
                    int(input_lines[i - 1][0:5].strip()) + 1,  #residue number
                    atom_residue_name,  #residue name
                    atom_name,  #atom name
                    int(input_lines[i - 1][15:20].strip()) + 1,  #atom number
                    coordinates[0],
                    coordinates[1],
                    coordinates[2],
                    velocity_string)

            break

    write_file.write_file(input_lines, output_gro_file)


def merge_gro_files(input_gro_file_1,
                    input_gro_file_2,
                    output_gro_file,
                    choose_box=1,
                    box_lenghts=None):
    """merge 2 gro files

    can come in handy to mix a particle with a box of water
    for alchemical transformations

    Parameters
    -------------
    input_gro_file_1 : str
        the name (or path) of the input gro file nr 1
    input_gro_file_2 : str
        the name (or path) of the input gro file nr 2
    output_gro_file : str
        the name (or path) of the input gro file
        can be the same as one of the input ones
    choose_box : int, optional
        the last line of a gro file contains the box
        lenghts, if `choose_box` = 1 (default) the output
        gro file will have the box values of `input_gro_file_1`
        if `choose_box` = 1 the box values of `input_gro_file_2`
        if `choose_box` = None the function will read `box_lenghts`
    box_lenghts : iterable, optional
        if `choose_box` = None the box lenghts will be read from
        `box_lenghts`, it must be an iterable with 3 floating point
        values. Gromacs measures in nm
    """

    input_lines_1 = read_file.read_file(input_gro_file_1)

    input_lines_2 = read_file.read_file(input_gro_file_2)

    if choose_box == 1:

        box = input_lines_1[-1].strip().split()

    elif choose_box == 2:

        box = input_lines_2[-1].strip().split()

    elif choose_box is None:

        box = []

        for i in range(len(box_lenghts)):

            box.append('{:.5f}'.format(box_lenghts[i]))

    else:
        raise ValueError(f'{choose_box} is not a valid value for choose_box')

    output_lines = ['Merged gro files\n']

    #adds number of atoms
    output_lines.append(' {:>5}\n'.format(int(input_lines_1[1].strip()) + \
        int(input_lines_2[1].strip())))

    atom = 1
    residue = 1
    pevious_residue = '1'
    for lines in (input_lines_1, input_lines_2):

        for i in range(2, len(lines) - 1):

            if lines[i].strip() != '':

                #check for residue change
                if lines[i][:5].strip() != pevious_residue:

                    pevious_residue = lines[i][:5].strip()

                    residue += 1

                output_lines.append('{:>5}{:<5}{:>5}{:>5}{}\n'.format(
                    residue,  #residue number
                    lines[i][5:10],  #residue name
                    lines[i][10:15],  #atom name
                    atom,  #atom number
                    lines[i][20:].strip('\n')))

                #update atom number
                atom += 1

    del input_lines_1
    del input_lines_2

    output_lines.append(2 * ' ' + '  '.join(box) + '\n')

    write_file.write_file(output_lines, output_gro_file)


def remove_velocities(input_gro_file, output_gro_file, keep_velocities=None):
    """removes the velocities from a gro file

    Some times gromacs has some problems of instability
    if the velocities are given in the gro file,
    in that case it is better to remove them

    Pay attention if you have a heavy dummy atom, that should
    always have zero velocity, for that use `keep_velocities`

    Parameters
    ------------
    input_gro_file : str
        the name (or path) of the input gro file nr 2
    output_gro_file : str
        the name (or path) of the input gro file
        can be the same as one of the input ones
    keep_velocities : list of strings
        the list of residue names (key sensitive) that should
        keep the velocity value (useful for dummy atoms),
        dafalut None no residue keeps its velocity
    """

    input_lines = read_file.read_file(input_gro_file)

    if keep_velocities is None:
        keep_velocities = []

    output_lines = []
    output_lines.append(input_lines[0])
    output_lines.append(input_lines[1])

    for i in range(2, len(input_lines) - 1):

        if input_lines[i][5:10].strip() not in keep_velocities:

            output_lines.append(input_lines[i][:44] + '\n')

        else:

            output_lines.append(input_lines[i])

    output_lines.append(input_lines[-1])

    write_file.write_file(output_lines, output_gro_file)
