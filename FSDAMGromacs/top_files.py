# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""functions to deal with TOP files
"""

import PythonAuxiliaryFunctions.files_IO.read_file as read_file
import PythonAuxiliaryFunctions.files_IO.write_file as write_file


def add_string_after_FF(input_string, input_top_file, output_top_file):
    """adds a string after the FF one

    In case there is no FF include (es a parmed genereted topology) it will
    be added after deafults

    Parameters
    -----------
    input_string : str
    input_top_file : str
    output_top_file : str
        can be the same as the input one
    """

    input_top_lines = read_file.read_file(input_top_file)

    for i in range(len(input_top_lines)):

        if input_top_lines[i].strip() != '':
            if input_top_lines[i].strip()[0] != ';':

                if input_top_lines[i].strip()[0:8] == '#include':

                    input_top_lines[i] += f'\n{input_string}\n'

                    break

                # In case there is no FF include (es a parmed genereted topology)
                if input_top_lines[i].strip()[0] == '[' and \
                    input_top_lines[i].split(';')[0].strip().replace(' ', '') != '[defaults]':

                    input_top_lines[
                        i] = f'\n{input_string}\n' + input_top_lines[
                            i]

                    break

    write_file.write_file(input_top_lines, output_top_file)


def add_include_after_FF(include_line, input_top_file, output_top_file):
    """adds an include statement after the FF one

    In case there is no FF include (es a parmed genereted topology) it will
    be added after deafults

    Parameters
    -----------
    include_line : str
        the thing to include, DON'T write #include
        but only the itp file name (or path)
    input_top_file : str
    output_top_file : str
        can be the same as the input one

    Notes
    ----------
    if you have to include multiple itp files all the [ atomtypes ]
    must be after the force field include and only then you can add the remaining parts of the
    itp files
    """

    add_string_after_FF(f'#include "{include_line}"', input_top_file, output_top_file)


def add_string_after_atomtypes(input_string, input_top_file, output_top_file):
    """adds a string after atomtypes

    It checks for the beginning of a [ ... ] section that is not [ atomtypes ]
    and for the beginning of a #ifdef

    Parameters
    -----------
    input_string : str
    input_top_file : str
    output_top_file : str
        can be the same as the input one
    """
    def is_right_line(line):
        """complex bool expession
        check both for the end of [ atomtypes ] and for the beginning of a #ifdef
        and skips [ defaults ] and [ cmaptypes ]
        """
        _line = line.split(';')[0].strip()

        if _line:
            if _line[0] == '[' and \
                _line.replace(' ', '') not in ('[atomtypes]', '[defaults]', '[cmaptypes]'):

                return True

            if _line[:6] == '#ifdef':
                return True

        return False

    input_top_lines = read_file.read_file(input_top_file)

    for i in range(len(input_top_lines)):

        if input_top_lines[i].strip():

            if is_right_line(input_top_lines[i]):

                input_top_lines[i] = \
                f'\n{input_string}\n{input_top_lines[i].strip()}\n'

                break

    # If end of file was reached try to put it in the end
    # as last resort
    else:
        input_top_lines[-1] += f'\n{input_string}\n'

    write_file.write_file(input_top_lines, output_top_file)


def add_include_after_atomtypes(include_line, input_top_file, output_top_file):
    """adds an include statement after atomtypes

    It checks for the beginning of a [ ... ] section that is not [ atomtypes ]
    and for the beginning of a #ifdef

    Parameters
    -----------
    include_line : str
        the thing to include, DON'T write #include
        but only the itp file name (or path)
    input_top_file : str
    output_top_file : str
        can be the same as the input one

    Notes
    ----------
    if you have to include multiple itp files all the [ atomtypes ]
    must be after the force field include and only then you can add the remaining parts of the
    itp files
    """

    add_string_after_atomtypes(f'#include "{include_line}"', input_top_file, output_top_file)


def add_molecules(name, number, input_top_file, output_top_file):
    """adds a [ molecules ] statement

    Parameters
    -----------
    name : str
        the name of the moelcule
    number : int
        the number of molecules
    input_top_file : str
    output_top_file : str
        can be the same as the input one
    """

    input_top_lines = read_file.read_file(input_top_file)

    input_top_lines.append(f'\n{name}     {number}\n')

    write_file.write_file(input_top_lines, output_top_file)
