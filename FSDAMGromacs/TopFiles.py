# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""functions to deal with TOP files
"""


def add_include(include_line, input_top_file, output_top_file):
    """adds an include statement after the FF one

    Parameters
    -----------
    include_line : str
        the thing to include, DON'T write #include
        but only the itp file name (or path)
    input_top_file : str
    output_top_file : str
        can be the same as the input one
    """

    with open(input_top_file, 'r') as f:

        input_top_lines = f.readlines()

    for i in range(len(input_top_lines)):

        if not (input_top_lines[i].strip() != ''
                or input_top_lines[i].strip()[0] != ';'):

            if input_top_lines[i].strip()[0:8] == '#include':

                input_top_lines[i] += f'\n#include "{include_line}"\n'

                break

    with open(output_top_file, 'w') as f:

        f.writelines(input_top_lines)


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

    with open(input_top_file, 'r') as f:

        input_top_lines = f.readlines()

    input_top_lines[-1] += f'{name}     {number}\n'

    with open(output_top_file, 'w') as f:

        f.writelines(input_top_lines)
