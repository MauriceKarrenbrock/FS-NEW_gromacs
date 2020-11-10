# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name
# pylint: disable=no-self-use
# pylint: disable=protected-access
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################

import FSDAMGromacs.tpr_files as tpr_files


class Testmake_tpr_file():
    def test_works(self, mocker):

        m_run = mocker.patch('PythonAuxiliaryFunctions.run.subprocess_run')

        mdp = 'mdp.mdp'
        gro = 'gro.gro'
        top = 'top.top'
        tpr = 'tpr.tpr'
        gmx = 'gmx'

        command = [
            f'{gmx}', 'grompp', '-f', f'{mdp}', '-c', f'{gro}', '-p', f'{top}',
            '-maxwarn', '100', '-o', f'{tpr}'
        ]

        error = \
        'Could not create the TPR file, check the stdout and stderr printed above for more info'

        tpr_files.make_tpr_file(mdp, gro, top, tpr, gmx)

        m_run.assert_called_once_with(commands=command,
                                      shell=False,
                                      universal_newlines=True,
                                      error_string=error)
