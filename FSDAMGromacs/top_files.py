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


def add_include_after_FF(include_line, input_top_file, output_top_file):
    """adds an include statement after the FF one

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

    input_top_lines = read_file.read_file(input_top_file)

    for i in range(len(input_top_lines)):

        if input_top_lines[i].strip() != '':
            if input_top_lines[i].strip()[0] != ';':

                if input_top_lines[i].strip()[0:8] == '#include':

                    input_top_lines[i] += f'\n#include "{include_line}"\n'

                    break

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

    input_top_lines = read_file.read_file(input_top_file)

    for i in range(len(input_top_lines)):

        if input_top_lines[i].strip():

            #complex bool expession
            #check both for the end of [ atomtypes ] and for the beginning of a #ifdef
            is_right_line = (  # pylint: disable=consider-using-ternary
                ((input_top_lines[i].strip()[0] == '[') and
                 (not '[ atomtypes ]' in input_top_lines[i].strip()))
                or (input_top_lines[i].strip()[0:6] == '#ifdef'))

            if is_right_line:

                input_top_lines[i] = \
                f'\n#include "{include_line}"\n{input_top_lines[i].strip()}\n'

                break

    write_file.write_file(input_top_lines, output_top_file)


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
