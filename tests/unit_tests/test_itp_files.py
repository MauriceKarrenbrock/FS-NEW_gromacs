# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name
# pylint: disable=no-self-use
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################

import FSDAMGromacs.itp_files as itp_files


class Testcreate_dummy_atom_itp():
    def test_works(self):

        name = 'DUM'

        atom_number = 22

        mass = 1.E+30

        charge = 0.

        sigma = 0

        epsilon = 0

        charge_group = 1

        atomtypes = [
            '; N.B: SIGMA is in nm and EPSILON in kj mol-1\n',
            '[ atomtypes ]\n',
            ';NAME   AT.NUM  MASS     CHARGE  PTYPE  SIGMA   EPSILON\n',
            f'{name[0:2]}         {atom_number}    {mass}    {charge}   A     {sigma}   {epsilon}\n',  # pylint: disable=line-too-long
            '\n'
        ]

        moleculetype = [
            '; Molecule topology/parameters starts below', '[ moleculetype ]',
            '; Name               nrexcl', f'  {name}                  3', ''
        ]

        atoms = [
            '; Atomic types, pdb names  and groups are defined below',
            '',
            '[ atoms ]',
            ';AT.NUM   TYPE   RESID   RESNAME PDB-NAME  IGRP    CHRGE',
            '    MASS',
            f'  1   {name[0:2]}     1      {name}     {name}      {charge_group}    {charge}    {mass}',  # pylint: disable=line-too-long
            ''
        ]

        itp_file = moleculetype + atoms

        for i in range(len(itp_file)):

            itp_file[i] += '\n'

        output = itp_files.create_dummy_atom_itp(mass, charge, sigma, epsilon,
                                                 name, charge_group,
                                                 atom_number)

        assert output[0] == atomtypes
        assert output[1] == itp_file
