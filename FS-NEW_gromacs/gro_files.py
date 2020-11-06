# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Functions to change or extract gro files
"""

import subprocess


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
    """

    if output_gro[-4:] != '.gro':
        output_gro += '.gro'

    string = f'echo System | {gromacs_exe} trjconv -f {xtc_file} -s \
        {tpr_file} -pbc mol -o {output_gro} -sep yes -ur compact -dt \
            {timestep_ps} -b {first_frame_ps}'

    if last_frame_ps != -1:
        string += f'-e {last_frame_ps}'

    r = subprocess.run(string,
                       shell=True,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=False)

    print(r.stdout)
    print(r.stderr)

    if r.returncode != 0:
        raise RuntimeError(
            'Could not extract configurations, look at stderr stdout printed above'
        )
