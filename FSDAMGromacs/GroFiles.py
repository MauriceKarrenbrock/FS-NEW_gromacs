# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Functions to change or extract gro files
"""

import PythonAuxiliaryFunctions.Run as Run


def extract_gro_from_xtc(xtc_file,
                         tpr_file,
                         output_gro,
                         timestep_ps=100,
                         first_frame_ps=0,
                         last_frame_ps=-1,
                         gromacs_exe='gmx'):
    """Extracts frames from xtc file at given timestep

    Given a xtc and tpr file it extracts trajectories in gro
    format at a certain timestep, starting from a selected frame
    (default the first one) till another selected frame (default the last
    one).
    the output ones will be <output_gro>0.gro <output_gro>1.gro ...

    Parameters
    -----------
    xtc_file : str
        the file name (or path) of the xtc file
    tpr_file : str
        the file name (or path) of the xtc file
    output_gro : str
        output name of the gro files
        the output ones will be numered
        <output_gro>0.gro <output_gro>1.gro ...
    timestep_ps : int or float
        the delta t of the frames to extract in ps
        default 100
    first_frame_ps : int
        first frame to extract in ps
        default 0 (the first one)
    last_frame_ps : int
        last frame to extract in ps
        default -1 (the last one)
    gromacs_exe : str
        the gromacs executable to use
        (absolute path)
        defauld "gmx"

    Notes
    -----------
    This function uses your UNIX shell and specifically `echo` and `|` commands
    """

    if output_gro[-4:] != '.gro':
        output_gro += '.gro'

    string = f'echo System | {gromacs_exe} trjconv -f {xtc_file} -s \
        {tpr_file} -pbc mol -o {output_gro} -sep yes -ur compact -dt \
            {timestep_ps} -b {first_frame_ps}'

    if last_frame_ps != -1:
        string += f'-e {last_frame_ps}'

    error = 'Could not extract configurations, look at stderr stdout printed above'

    Run.subprocess_run(commands=string,
                       shell=True,
                       universal_newlines=True,
                       error_string=error)


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
