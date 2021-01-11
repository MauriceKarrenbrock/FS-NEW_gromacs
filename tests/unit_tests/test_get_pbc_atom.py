# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name
# pylint: disable=no-self-use
# pylint: disable=protected-access
# pylint: disable=duplicate-code
# pylint: disable=too-many-lines
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################

import numpy as np

import FSDAMGromacs.get_pbc_atom as pbc


class Testget_protein_pbc_atom():
    def test_works_str(self, mocker):

        Geo = np.array([0., 0., 0.])

        class Atom():
            def __init__(self, coord, number):

                self.coord = coord

                self.number = number

            def get_serial_number(self):

                return self.number

        class Molecule():
            def __init__(self, molecule):

                self.molecule = molecule

            def get_atoms(self):

                return self.molecule

        atom_a = Atom(np.array([1., 1., 1.]), 1)

        atom_b = Atom(np.array([2., 2., 2.]), 2)

        atom_c = Atom(np.array([-3., -3., -3.]), 3)

        molecule = Molecule((atom_a, atom_b, atom_c))

        m_mdanalysis = mocker.patch('MDAnalysis.Universe')
        m_bio = mocker.patch(
            'PythonPDBStructures.pdb.biopython_utils.parse_pdb',
            return_value=molecule)
        m_chain = mocker.patch(
            'PythonPDBStructures.pdb.add_chain_id.add_chain_id_pdb')
        m_com = mocker.patch('PythonPDBStructures.geometry.get_center_of_mass',
                             return_value=Geo)

        output = pbc.get_protein_pbc_atom('DUMMY.pdb', 'some_selection')

        assert output == 1

        m_mdanalysis.assert_called_once()
        m_bio.assert_called_once()
        m_chain.assert_called_once()
        m_com.assert_called_once()
