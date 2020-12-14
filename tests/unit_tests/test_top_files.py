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

import FSDAMGromacs.top_files as top_files


class Testadd_include_after_FF():
    def test_works(self, mocker):

        top_file = ['\n', ';\n', 'AAAAAA\n', '#include "stuff"\n', 'AAAAAA\n']

        expected = [
            '\n', ';\n', 'AAAAAA\n',
            '#include "stuff"\n\n#include "new_included_stuff"\n', 'AAAAAA\n'
        ]

        m_read = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.read_file.read_file',
            return_value=top_file)

        m_write = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.write_file.write_file')

        top_files.add_include_after_FF('new_included_stuff', 'input', 'output')

        m_read.assert_called_once_with('input')

        m_write.assert_called_once_with(expected, 'output')


class Testadd_include_after_atomtypes():
    def test_works(self, mocker):

        top_file = [
            '\n', ';\n', 'AAAAAA\n', '#include "stuff"\n', 'AAAAAA\n',
            '[ atomtypes ]\n', 'AAAAA\n', '[ moleculetypes ]\n', 'AAAAA\n'
        ]

        expected = [
            '\n', ';\n', 'AAAAAA\n', '#include "stuff"\n', 'AAAAAA\n',
            '[ atomtypes ]\n', 'AAAAA\n',
            '\n#include "new_included_stuff"\n[ moleculetypes ]\n', 'AAAAA\n'
        ]

        m_read = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.read_file.read_file',
            return_value=top_file)

        m_write = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.write_file.write_file')

        top_files.add_include_after_atomtypes('new_included_stuff', 'input',
                                              'output')

        m_read.assert_called_once_with('input')

        m_write.assert_called_once_with(expected, 'output')


class Testadd_molecules():
    def test_works(self, mocker):

        name = 'DUM'
        number = 1

        top_file = ['Stuff\n', 'ssssss\n']

        expected = ['Stuff\n', 'ssssss\n', f'\n{name}     {number}\n']

        m_read = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.read_file.read_file',
            return_value=top_file)

        m_write = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.write_file.write_file')

        top_files.add_molecules(name, number, 'input', 'output')

        m_read.assert_called_once_with('input')

        m_write.assert_called_once_with(expected, 'output')
