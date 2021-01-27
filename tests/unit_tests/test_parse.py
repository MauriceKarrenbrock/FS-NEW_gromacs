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

        parsed_stuff = [[0., -1], [1., -2.], [2., -3.], [3., -4.]]

        expected_output = [[0., 1. / 3., 2. * 1. / 3., 3 * 1. / 3.],
                           [-1., -2., -3., -4.]]

        parsed_stuff = np.array(parsed_stuff)

        expected_output = np.array(expected_output)

        expected_output[1, :] = expected_output[1, :] * 0.23901

        mocked_parse = mocker.patch(
            'FSDAMGromacs.parse.parse_big_gromacs_xvg_files',
            return_value=parsed_stuff)

        instance = parse.GromacsParseWorkProfile('gromacs')

        output = instance.parse('dhdl.xvg',
                                starting_lambda=0.,
                                ending_lambda=1.)

        mocked_parse.assert_called_once_with('dhdl.xvg')

        assert np.testing.assert_allclose(output, expected_output) is None


class TestGromacsParsePullDistances():
    def test_parse(self, mocker):

        parsed_list = [
            '# gmx mdrun is part of G R O M A C S:',
            '#',
            '# GROningen Mixture of Alchemy and Childrens\' Stories',
            '#',
            '@    title "Pull COM"',
            '@    xaxis  label "Time (ps)"',
            '@    yaxis  label "Position (nm)"',
            '@TYPE xy',
            '@ view 0.15, 0.15, 0.75, 0.85',
            '@ legend on',
            '@ legend box on',
            '@ legend loctype view',
            '@ legend 0.78, 0.8',
            '@ legend length 2',
            '@ s0 legend "1"',
            '@ s1 legend "1 g 1 X"',
            '@ s2 legend "1 g 1 Y"',
            '@ s3 legend "1 g 1 Z"',
            '@ s4 legend "1 g 2 X"',
            '@ s5 legend "1 g 2 Y"',
            '@ s6 legend "1 g 2 Z"',
            '@ s7 legend "2"',
            '@ s8 legend "2 g 1 X"',
            '@ s9 legend "2 g 1 Y"',
            '@ s10 legend "2 g 1 Z"',
            '@ s11 legend "2 g 2 X"',
            '@ s12 legend "2 g 2 Y"',
            '@ s13 legend "2 g 2 Z"',
            ' 0.0   1.0   2.0   3.0   4.0   5.0   6.0   7.0   8.0   9.0  10.0  11.0  12.0  13.0  14.',  # pylint: disable=line-too-long
            '15.0  16.0  17.0  18.0  19.0  20.0  21.0  22.0  23.0  24.0  25.0  26.0  27.0  28.0  29.',  # pylint: disable=line-too-long
            '30.0  31.0  32.0  33.0  34.0  35.0  36.0  37.0  38.0  39.0  40.0  41.0  42.0  43.0  44.',  # pylint: disable=line-too-long
            '45.0  46.0  47.0  48.0  49.0  50.0  51.0  52.0  53.0  54.0  55.0  56.0  57.0  58.0  59.'  # pylint: disable=line-too-long
        ]

        parsed_array = [[
            0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14.
        ],
                        [
                            15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                            25., 26., 27., 28., 29.
                        ],
                        [
                            30., 31., 32., 33., 34., 35., 36., 37., 38., 39.,
                            40., 41., 42., 43., 44.
                        ],
                        [
                            45., 46., 47., 48., 49., 50., 51., 52., 53., 54.,
                            55., 56., 57., 58., 59.
                        ]]

        parsed_array = np.array(parsed_array)

        mocked_parse = mocker.patch(
            'FSDAMGromacs.parse.parse_big_gromacs_xvg_files',
            return_value=parsed_array)

        m_list = mocker.patch('FSDAMGromacs.parse.open')

        m_list.return_value.__enter__.return_value = parsed_list

        instance = parse.GromacsParsePullDistances('gromacs')

        output = instance.parse('dummy_pullx.xvg')

        assert isinstance(output, dict)
        assert tuple(output.keys()) == (1, 2)
        assert np.testing.assert_allclose(output[1],
                                          parsed_array[:, 1] * 10.) is None
        assert np.testing.assert_allclose(output[2],
                                          parsed_array[:, 8] * 10.) is None

        mocked_parse.assert_called_once_with('dummy_pullx.xvg')

        m_list.assert_called_once()


class Testparse_big_gromacs_xvg_files():
    def test_works_default(self, mocker):

        parsed_list = [
            '# gmx mdrun is part of G R O M A C S:',
            '#',
            '# GROningen Mixture of Alchemy and Childrens\' Stories',
            '#',
            '@    title "Pull COM"',
            '@    xaxis  label "Time (ps)"',
            '@    yaxis  label "Position (nm)"',
            ' 0.0   1.0   2.0   3.0   4.0   5.0   6.0   7.0   8.0   9.0  10.0  11.0  12.0  13.0  14.',  # pylint: disable=line-too-long
            '15.0  16.0  17.0  18.0  19.0  20.0  21.0  22.0  23.0  24.0  25.0  26.0  27.0  28.0  29.',  # pylint: disable=line-too-long
            '30.0  31.0  32.0  33.0  34.0  35.0  36.0  37.0  38.0  39.0  40.0  41.0  42.0  43.0  44.',  # pylint: disable=line-too-long
            '45.0  46.0  47.0  48.0  49.0  50.0  51.0  52.0  53.0  54.0  55.0  56.0  57.0  58.0  59.'  # pylint: disable=line-too-long
        ]

        m_open = mocker.patch('FSDAMGromacs.parse.open')

        m_open.return_value.__enter__.return_value = parsed_list

        m_pandas_csv = mocker.patch('pandas.read_csv', return_value='csv')

        m_numpy = mocker.patch('numpy.array', return_value='array')

        output = parse.parse_big_gromacs_xvg_files('dummy_file.xvg')

        assert output == 'array'

        m_open.assert_called_once_with('dummy_file.xvg', 'r')

        m_pandas_csv.assert_called_once_with('dummy_file.xvg',
                                             skiprows=7,
                                             skipinitialspace=True,
                                             sep=' ',
                                             header=None)

        m_numpy.assert_called_once_with('csv')
