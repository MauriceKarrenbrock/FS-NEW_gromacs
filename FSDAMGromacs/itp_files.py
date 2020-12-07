# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""functions to modify itp files
"""


def create_dummy_atom_itp(mass,
                          charge,
                          sigma,
                          epsilon,
                          name='DUM',
                          charge_group=1,
                          atom_number=1):
    """creates a dummy atom itp file (list of strings)

    Parameters
    ------------
    mass : float
        the mass of the atom
    charge : float
        the charge of the atom
    sigma : float
        the sigma of LJ potential
        unit: nm
    epsilon : float
        the epsilon of LJ potential
        unit: kj mol-1
    name : str, optional
        MAX 3 characters, the residue name of the atom
        default DUM
    charge_group : int, optional
        the charge group for some cutoff schemes
        dafault 1
    atom_number : int, optional
        default 1

    Returns
    -----------
    list of strings
        all the lines already have the newline simbol at the end
    """

    atomtypes = [
        '; N.B: SIGMA is in nm and EPSILON in kj mol-1', '[ atomtypes ]',
        ';NAME   AT.NUM  MASS     CHARGE  PTYPE  SIGMA   EPSILON',
        f'{name[0:2]}         {atom_number}    {mass}    {charge}   A     {sigma}   {epsilon}',
        ''
    ]

    moleculetype = [
        '; Molecule topology/parameters starts below', '[ moleculetype ]',
        '; Name               nrexcl', f'  {name}                  3', ''
    ]

    atoms = [
        '; Atomic types, pdb names  and groups are defined below', '',
        '[ atoms ]',
        ';AT.NUM   TYPE   RESID   RESNAME PDB-NAME  IGRP    CHRGE', '    MASS',
        f'  1   {name[0:2]}     1      {name}     {name}      {charge_group}    {charge}    {mass}',
        ''
    ]

    itp_file = atomtypes + moleculetype + atoms

    for i in range(len(itp_file)):

        itp_file[i] += '\n'

    return itp_file
