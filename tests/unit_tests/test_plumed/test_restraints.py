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

import pytest

import FSDAMGromacs.plumed.restraints as restraints


class TestCOM_COM_restraint():
    def test_raises(self):

        with pytest.raises(ValueError):

            restraints.COM_COM_restraint({'aaa': [1, 2]}, [['DUM']])

    @pytest.mark.parametrize('test_type, geometric, com',
                             [('geometric_false', False, 'COM'),
                              ('geometric_true', True, 'CENTER')])
    def test_works(self, mocker, test_type, geometric, com):

        print('Logging test type for visibility: ' + test_type)

        m_write = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.write_file.write_file')

        atom_groups = {
            'group_1': [1, 2, 3, 4],
            'group_2': [5, 6, 7, 8],
            'group_3': [9, 10, 11, 12]
        }

        restraint_parameters = [['group_1', 'group_2', 1., 2., 3.],
                                ['group_1', 'group_3', 4., 5., 6.],
                                ['group_2', 'group_3', 0.7, 0.8, 0.9]]

        plumed_file = 'test_plumed.dat'
        distances_file = 'test_distances.out'
        stride = 10

        expected = [
            'UNITS LENGTH=nm TIME=ps \n\n',
            'WHOLEMOLECULES ENTITY0=1,2,3,4 ENTITY1=5,6,7,8 ENTITY2=9,10,11,12 \n',
            f'group_1: {com} ATOMS=1,2,3,4\n',
            f'group_2: {com} ATOMS=5,6,7,8\n',
            f'group_3: {com} ATOMS=9,10,11,12\n',
            'group_1_group_2_dist: DISTANCE ATOMS=group_1,group_2 NOPBC\n',
            'group_1_group_3_dist: DISTANCE ATOMS=group_1,group_3 NOPBC\n',
            'group_2_group_3_dist: DISTANCE ATOMS=group_2,group_3 NOPBC\n',
            'RESTRAINT ARG=group_1_group_2_dist AT=1.0000 KAPPA=2.0000 SLOPE=3.0000\n',
            'RESTRAINT ARG=group_1_group_3_dist AT=4.0000 KAPPA=5.0000 SLOPE=6.0000\n',
            'RESTRAINT ARG=group_2_group_3_dist AT=0.7000 KAPPA=0.8000 SLOPE=0.9000\n',
            f'PRINT ARG=group_1_group_2_dist,group_1_group_3_dist,group_2_group_3_dist STRIDE={stride} FILE={distances_file}\n'  # pylint: disable=line-too-long
        ]

        restraints.COM_COM_restraint(atom_groups, restraint_parameters,
                                     plumed_file, distances_file, stride,
                                     geometric)

        m_write.assert_called_once_with(expected, plumed_file)
