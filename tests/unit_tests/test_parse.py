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

import numpy as np

import FSDAMGromacs.parse as parse


class TestGromacsParseWorkProfile():
    def test_works_parse(self, mocker):

        parsed_stuff = [[0., 1., 2., 3.], [-1., -2., -3., -4.]]

        expected_output = [[0., 1. / 3., 2. * 1. / 3., 3 * 1. / 3.],
                           [-1., -2., -3., -4.]]

        parsed_stuff = np.array(parsed_stuff)

        expected_output = np.array(expected_output)

        mocked_load = mocker.patch('numpy.loadtxt', return_value=parsed_stuff)

        instance = parse.GromacsParseWorkProfile('gromacs')

        output = instance.parse('dhdl.xvg')

        mocked_load.assert_called_once_with('dhdl.xvg',
                                            comments=['#', '@'],
                                            delimiter=' ')

        assert np.testing.assert_allclose(output, expected_output) is None
