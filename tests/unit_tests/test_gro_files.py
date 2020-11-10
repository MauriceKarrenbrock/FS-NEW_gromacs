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

import FSDAMGromacs.gro_files as gro_files


class Testgro_files():
    def test_works(self, mocker):

        fake_input_gro = [
            'test-g\n', '    6\n',
            '    1LIG    c01    1   0.125   0.004   0.092  0.0780  0.2389  0.2170\n',
            '    1LIG    h01    2   0.161  -0.096   0.067 -0.0004  0.8299 -2.4035\n',
            '    1LIG   cl01    3   0.227   0.091   0.195  0.0790 -0.1586 -0.2172\n',
            '    1LIG    c02    4   0.018   0.063   0.040 -0.1374  0.4773 -0.1111\n',
            '    1LIG    h02    5  -0.010   0.165   0.067  0.4681  0.2968  1.2689\n',
            '    1LIG   cl02    6  -0.084  -0.012  -0.068 -0.0722 -0.1160  0.2136\n',
            '  50.00000  50.00000  50.00000\n'
        ]

        expected_output_gro = [
            'test-g\n',
            '    7\n',
            '    1LIG    c01    1   0.125   0.004   0.092  0.0780  0.2389  0.2170\n',
            '    1LIG    h01    2   0.161  -0.096   0.067 -0.0004  0.8299 -2.4035\n',
            '    1LIG   cl01    3   0.227   0.091   0.195  0.0790 -0.1586 -0.2172\n',
            '    1LIG    c02    4   0.018   0.063   0.040 -0.1374  0.4773 -0.1111\n',
            '    1LIG    h02    5  -0.010   0.165   0.067  0.4681  0.2968  1.2689\n',
            '    1LIG   cl02    6  -0.084  -0.012  -0.068 -0.0722 -0.1160  0.2136\n    2MUD     MU    7   0.011   0.022   0.033  0.0011  0.0022  0.0033\n',  # pylint: disable=line-too-long
            '  50.00000  50.00000  50.00000\n'
        ]

        mocker.patch('PythonAuxiliaryFunctions.files_IO.read_file.read_file',
                     return_value=fake_input_gro)

        moked_write = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.write_file.write_file')

        gro_files.add_atom_to_gro_file('input_path', 'output_path',
                                       (0.011, 0.022, 0.033),
                                       (0.0011, 0.0022, 0.0033), 'MU', 'MUD')

        moked_write.assert_called_once_with(expected_output_gro, 'output_path')
