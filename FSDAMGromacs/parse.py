# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""functions and classes to parse gromacs output files
"""

import re

import numpy as np
import pandas as pd
import PythonFSDAM.parse.parse_superclasses as superclasses
from simtk import unit


def parse_big_gromacs_xvg_files(file_name, comments=None):
    """parses quickly big gromacs xvg files

    as numpy.loadtxt is too slow to parse houndreds of files
    I did this helper function to use pandas.read_csv

    this function cannot deal with footers but only headers

    Parameters
    ------------
    file_name : str or path
        the file to parse
    comments : interable(str), optional, default=['#', '@']
        lines in the header to jump

    Returns
    -------------
    np.array
        like the one obtained from loadtxt on the same file
    """

    if comments is None:

        comments = ['#', '@']

    with open(file_name, 'r') as f:

        j = 0
        for line in f:

            if line.strip()[0] not in comments:
                break

            j += 1

    output = pd.read_csv(file_name,
                         skiprows=j,
                         skipinitialspace=True,
                         sep=' ',
                         header=None)

    return np.array(output)


class GromacsParseWorkProfile(superclasses.ParseWorkProfileSuperclass):
    """The function to parse gromacs dhdl.xvg file

    this class will usually be called from
   `PythonFSDAM.parse.parse.ParseWorkProfile` factory object
    and inherits from `PythonFSDAM.parse.parse_superclasses.ParseWorkProfileSuperclass`
    so check the documentations for missing informations
    https://github.com/MauriceKarrenbrock/PythonFSDAM

    This class will take for granted that the first column
    is time and the second is dH/dL

    This class will convert the time column in lambda column
    """
    @staticmethod
    def parse(file_name, starting_lambda=1.0, ending_lambda=0.):  # pylint: disable=arguments-differ
        """Parses a dhdl.xvg file

        converts the work from KJ/mol to Kcal/mol (1 Kcal = 0.23901 KJ)

        Parameters
        -----------
        file_name : str
            the file to parse
        starting_lambda : float, optional, default=1.0
            the value of lambda at the beginning of the run
            usually 1.0 for annihilation and 0. for creation
        ending_lambda : float, optional, default=0.0
            the value of lambda at the end of the run
            usually 0. for annihilation and 1.0 for creation


        Returns
        ---------
        a 2-D array containing lambda and dh/dl (converted in Kcal/mol)
        [
            [lambda1, lambda2, ...],
            [dh/dl, dh/dl, ...]
        ]

        Notes
        ---------
        will take for granted that the first column
        is time and the second is dH/dL and that the value of
        lambda is gone from 0 to `abs_lambda_max_val` (or vice-versa) during the run
        with constant speed
        """

        if starting_lambda == ending_lambda:
            raise ValueError(
                'Beginning and end lambda cannot be the same value')

        parsed_file = parse_big_gromacs_xvg_files(file_name)

        #parsed file is time vs dhdl but I want lambda vs dhdl
        delta_lambda = (ending_lambda -
                        starting_lambda) / float(len(parsed_file[:, 0]) - 1)

        tmp = np.arange(len(parsed_file[:, 0]))

        # otherwise parsed_file[:, 0] would go from zero to -1
        if delta_lambda < 0:

            tmp = np.flip(tmp)

            delta_lambda = abs(delta_lambda)

        parsed_file[:, 0] = tmp * delta_lambda

        #in this way I have a line with all lambdas and a line with all
        #dhdl
        parsed_file = parsed_file.transpose()

        #from KJ/mol to Kcal/mol
        tmp = parsed_file[1] * unit.kilojoules_per_mole
        parsed_file[1] = tmp.value_in_unit(unit.kilocalories_per_mole)

        return parsed_file


class GromacsParsePullDistances(superclasses.Parser):
    """Parses the <something>_pullx.xvg file for COM-COM pulls
    """
    @staticmethod
    def parse(file_name):
        """Parses the <something>_pullx.xvg file for COM-COM pulls

        Parameters
        ------------
        file_name : pathlib.Path
            the path to the pullx.xvg file
            it takes for granted that the positions are in nanometers and
            converts them in angstrom (multiplies * 10), if they are in a different unit it's
            your job to deal with the unit change

        Return
        ---------
        dict of numpy.array
            this dictionary has the number of the pull group (int) as key
            and the COM position in ANGSTROMS!! as values (numpy.array)

        Notes
        -------------
        the returned dict will point to a bigger matrix so remember to eliminate it
        when you don't need it anymore to free memory
        """

        #This is an example of the file header:
        # s0 is actulally the second column because the first one is time (xaxis)

        # # gmx mdrun is part of G R O M A C S:
        # #
        # # GROningen Mixture of Alchemy and Childrens' Stories
        # #
        # @    title "Pull COM"
        # @    xaxis  label "Time (ps)"
        # @    yaxis  label "Position (nm)"
        # @TYPE xy
        # @ view 0.15, 0.15, 0.75, 0.85
        # @ legend on
        # @ legend box on
        # @ legend loctype view
        # @ legend 0.78, 0.8
        # @ legend length 2
        # @ s0 legend "1"
        # @ s1 legend "1 g 1 X"
        # @ s2 legend "1 g 1 Y"
        # @ s3 legend "1 g 1 Z"
        # @ s4 legend "1 g 2 X"
        # @ s5 legend "1 g 2 Y"
        # @ s6 legend "1 g 2 Z"
        # @ s7 legend "2"
        # @ s8 legend "2 g 1 X"
        # @ s9 legend "2 g 1 Y"
        # @ s10 legend "2 g 1 Z"
        # @ s11 legend "2 g 2 X"
        # @ s12 legend "2 g 2 Y"
        # @ s13 legend "2 g 2 Z"

        parsed_file = parse_big_gromacs_xvg_files(file_name)

        #from nm to angstrom
        parsed_file = parsed_file * 10.

        #check which columns are useful and should be  returned
        output_dict = {}

        column = 0

        with open(file_name, 'r') as f:

            for line in f:

                line = line.strip()

                if line[0] not in ('#', '@'):

                    break

                if line[0] == '@':

                    if line.split()[1][0] == 's':

                        column += 1

                        #contains only a number
                        if bool(
                                re.match('^[0-9]*$',
                                         line.split()[-1].strip('"'))):

                            output_dict[int(line.split()[-1].strip(
                                '"'))] = parsed_file[:, column]

        return output_dict
