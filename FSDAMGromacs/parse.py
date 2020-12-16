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
import PythonFSDAM.parse.parse_superclasses as superclasses


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
    def parse(file_name, abs_lambda_max_val=1.0):  # pylint: disable=arguments-differ
        """Parses a dhdl.xvg file

        converts the work from KJ/mol to Kcal/mol (1 Kcal = 1/4.148 KJ)

        Parameters
        -----------
        file_name : str
            the file to parse
        abs_lambda_max_val : float, optional
            the maximum absolute value that ambda had during the run
            usually 1.0 (default)

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
        lambda is gone from 0 to 1 (or vice-versa) during the run
        with constant speed
        """

        parsed_file = np.loadtxt(file_name, comments=['#', '@'], delimiter=' ')

        #parsed file is time vs dhdl but I want lambda vs dhdl
        delta_lambda = abs_lambda_max_val / float(len(parsed_file[:, 0]) - 1)

        tmp = np.arange(len(parsed_file[:, 0]))

        parsed_file[:, 0] = tmp * delta_lambda

        #from KJ/mol to Kcal/mol
        parsed_file[:, 1] = parsed_file[:, 1] * (1. / 4.148)

        #in this way I have a line with all lambdas and a line with all
        #dhdl
        return parsed_file.transpose()


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

        parsed_file = np.loadtxt(file_name, comments=['#', '@'], delimiter=' ')

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
