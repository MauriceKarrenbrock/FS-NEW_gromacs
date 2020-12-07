# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""functions and classes to parse gromacs output files
"""

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
    def parse(self, file_name, abs_lambda_max_val=1.0):  # pylint: disable=arguments-differ
        """Parses a dhdl.xvg file

        Parameters
        -----------
        file_name : str
            the file to parse
        abs_lambda_max_val : float, optional
            the maximum absolute value that ambda had during the run
            usually 1.0 (default)

        Returns
        ---------
        a 2-D array containing lambda and dh/dl

        Notes
        ---------
        will take for granted that the first column
        is time and the second is dH/dL and that the value of
        lambda is gone from 0 to 1 (or vice-versa) during the run
        with constant speed
        """

        parsed_file = np.loadtxt(file_name, comments=['#', '@'], delimiter=' ')

        #parsed file is time vs dhdl but I want lambda vs dhdl
        delta_lambda = abs_lambda_max_val / float(parsed_file.shape[1] - 1)

        #only time line
        iterator = np.nditer(parsed_file[0], flags=['c_index'])

        for i in iterator:  # pylint: disable=unused-variable

            parsed_file[0, iterator.index] = delta_lambda * iterator.index

        return parsed_file
