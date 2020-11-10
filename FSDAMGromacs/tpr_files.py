# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""functions to create tpr files
"""

import PythonAuxiliaryFunctions.run as run


def make_tpr_file(mdp_file, gro_file, top_file, tpr_file, gromacs_path='gmx'):
    """makes a tpr file

    it calls gromacs to make a tpr file locally
    if you want to use gromacs in another computer
    (ex an HPC cluster) this may not be the right function
    because it is better to create the tpr on the machine
    you will do the run

    Parameters
    -------------
    mdp_file : str
        the ABSOLUTE PATH to the file
    gro_file : str
        the ABSOLUTE PATH to the file
    top_file : str
        the ABSOLUTE PATH to the file
    tpr_file : str
        the ABSOLUTE PATH where to save the output file
    gromacs_path : str
        the ABSOLUTE PATH where to the gromacs executable
        (default gmx)
    """

    command = [
        f'{gromacs_path}', 'grompp', '-f', f'{mdp_file}', '-c', f'{gro_file}',
        '-p', f'{top_file}', '-maxwarn', '100', '-o', f'{tpr_file}'
    ]

    error = 'Could not create the TPR file, check the stdout and stderr printed above for more info'

    run.subprocess_run(commands=command,
                       shell=False,
                       universal_newlines=True,
                       error_string=error)
