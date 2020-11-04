# -*- coding: utf-8 -*-
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""functions to modify itp files
"""

import copy


def create_annihilation_itp_file(input_itp_lines):
    """returns the lines of an ITP file for annihilation

    Given a list of strings (lines) with the newline "\n" at the end
    it returns a new list of lines that can be written on file
    in order to get the ITP file for annihilation.
    There must be only ONE MOLECULE PER ITP FILE!
    In the ITP file every atom will be set for annihilation

    Parameters
    ------------
    input_itp_lines : list
        list of strings (lines of the itp file) with new line
        "\n" in the end

    Returns
    --------
    annihilation_itp : list
        lines of the new ITP file for annihilation with
        newline "\n" in the end

    Notes
    -------
    If you have a system with multiple molecules and you only want to
    annihilate one you have to have different ITP files and ad them to
    the TOP file with #include statements
    """

    annihilation_itp = copy.deepcopy(input_itp_lines)

    nonbond_params = [
        '[ nonbond_params ]\n', '; i  j    func  sigma           epsilon\n'
    ]

    atoms = []

    is_atomtype = False
    for i in range(len(annihilation_itp)):

        #jump empty lines and comments
        if annihilation_itp[i].strip() == '' or \
        annihilation_itp[i].strip()[0] == ';':

            pass

        elif annihilation_itp[i].strip().replace(' ', '') == '[atomtypes]':

            is_atomtype = True

        #create DUM_ atoms (annihilated ones)
        elif is_atomtype:

            line = annihilation_itp[i].strip().split()

            tmp_DUM = [
                'DUM_' + line[0].strip(), line[1].strip(), line[2].strip(),
                '0.0', line[4].strip(), '0.0', '0.0', '\n'
            ]

            tmp_DUM = ' '.join(tmp_DUM)

            #every atom will be daclared with its DUM_ version
            annihilation_itp[i] = annihilation_itp[i].strip() + '\n' + tmp_DUM

            #keeping the normal atoms in order to do the nonbond_params later
            atoms.append(line)

        elif is_atomtype and (annihilation_itp[i].strip()[0] == '['):

            #This will be the last thing to happen

            for k in range(len(atoms)):

                for j in range(k, len(atoms)):

                    #DUM_atom k
                    tmp_nonbond_params = 'DUM_' + atoms[k][0].strip()
                    #DUM_atom j
                    tmp_nonbond_params = tmp_nonbond_params + ' ' + 'DUM_' + atoms[
                        j][0].strip()
                    #func (always 1)
                    tmp_nonbond_params = tmp_nonbond_params + ' ' + '1'

                    #sigma (arithmetic awg of the two sigmas)
                    tmp_nonbond_params = tmp_nonbond_params + ' ' + \
                        f'{( float( atoms[k][5].strip() ) + float( atoms[j][5].strip() ) ) / 2.}'

                    #epsilon (geometric awg of the two epsilons)
                    tmp_nonbond_params = tmp_nonbond_params + ' ' + \
                        f'{( float( atoms[k][6].strip() ) * float( atoms[j][6].strip() ) ) ** 0.5}'

                    #newline
                    tmp_nonbond_params = tmp_nonbond_params + '\n'

                    #append to the cumulative ones
                    nonbond_params.append(tmp_nonbond_params)

            nonbond_params.append('\n\n')

            #add the nonbond_params in annihilation_itp
            annihilation_itp[i:i] = nonbond_params

            break

    is_atoms = False
    for i in range(len(annihilation_itp)):

        #jump empty lines and comments
        if annihilation_itp[i].strip() == '' or \
        annihilation_itp[i].strip()[0] == ';':

            pass

        elif annihilation_itp[i].strip().replace(' ', '') == '[atoms]':

            is_atoms = True

        elif is_atoms and annihilation_itp[i].strip()[0] == '[':

            break

        elif is_atoms:

            tmp_line = annihilation_itp[i].strip().split()

            annihilation_itp[i] = annihilation_itp[i].strip()

            #add DUM_ state b
            annihilation_itp[
                i] = annihilation_itp[i] + ' ' + 'DUM_' + tmp_line[1]

            annihilation_itp[
                i] = annihilation_itp[i] + ' ' + '0.0' + tmp_line[7]

            annihilation_itp[i] = annihilation_itp[i] + '\n'

    return annihilation_itp


def create_creation_itp_file(input_itp_lines):
    """returns the lines of an ITP file for creation

    Given a list of strings (lines) with the newline "\n" at the end
    it returns a new list of lines that can be written on file
    in order to get the ITP file for creation.
    There must be only ONE MOLECULE PER ITP FILE!
    In the ITP file every atom will be set for creation

    Parameters
    ------------
    input_itp_lines : list
        list of strings (lines of the itp file) with new line
        "\n" in the end

    Returns
    --------
    creation_itp : list
        lines of the new ITP file for creation with
        newline "\n" in the end

    Notes
    -------
    If you have a system with multiple molecules and you only want to
    annihilate one you have to have different ITP files and ad them to
    the TOP file with #include statements
    """

    creation_itp = copy.deepcopy(input_itp_lines)

    nonbond_params = [
        '[ nonbond_params ]\n', '; i  j    func  sigma           epsilon\n'
    ]

    atoms = []

    is_atomtype = False

    for i in range(len(creation_itp)):

        #jump empty lines and comments
        if creation_itp[i].strip() == '' or creation_itp[i].strip()[0] == ';':

            pass

        elif creation_itp[i].strip().replace(' ', '') == '[atomtypes]':

            is_atomtype = True

        #creating DUM_ atoms (annihilated)
        elif is_atomtype:

            line = creation_itp[i].strip().split()

            tmp_DUM = [
                'DUM_' + line[0].strip(), line[1].strip(), line[2].strip(),
                '0.0', line[4].strip(), '0.0', '0.0', '\n'
            ]

            tmp_DUM = ' '.join(tmp_DUM)

            #every atom will be daclared with its DUM_ version
            creation_itp[i] = creation_itp[i].strip() + '\n' + tmp_DUM

            #keeping the normal atoms in order to do the nonbond_params later
            atoms.append(line)

        elif is_atomtype and (creation_itp[i].strip()[0] == '['):

            #This will be the last thing to happen

            for k in range(len(atoms)):

                for j in range(k, len(atoms)):

                    #DUM_atom k
                    tmp_nonbond_params = 'DUM_' + atoms[k][0].strip()
                    #DUM_atom j
                    tmp_nonbond_params = tmp_nonbond_params + ' ' + 'DUM_' + atoms[
                        j][0].strip()
                    #func (always 1)
                    tmp_nonbond_params = tmp_nonbond_params + ' ' + '1'

                    #sigma (arithmetic awg of the two sigmas)
                    tmp_nonbond_params = tmp_nonbond_params + ' ' + \
                        f'{( float( atoms[k][5].strip() ) + float( atoms[j][5].strip() ) ) / 2.}'

                    #epsilon (geometric awg of the two epsilons)
                    tmp_nonbond_params = tmp_nonbond_params + ' ' + \
                        f'{( float( atoms[k][6].strip() ) * float( atoms[j][6].strip() ) ) ** 0.5}'

                    #newline
                    tmp_nonbond_params = tmp_nonbond_params + '\n'

                    #append to the cumulative ones
                    nonbond_params.append(tmp_nonbond_params)

            nonbond_params.append('\n\n')

            #add the nonbond_params in creation_itp
            creation_itp[i:i] = nonbond_params

            break

    is_atoms = False
    for i in range(len(creation_itp)):

        #jump empty lines and comments
        if creation_itp[i].strip() == '' or creation_itp[i].strip()[0] == ';':

            pass

        elif creation_itp[i].strip().replace(' ', '') == '[atoms]':

            is_atoms = True

        elif is_atoms and creation_itp[i].strip()[0] == '[':
            #end of [atoms]
            break

        elif is_atoms:

            tmp_line = creation_itp[i].strip().split()

            #add DUM_ state a
            creation_itp[i] = tmp_line[0]

            creation_itp[i] = creation_itp[i] + ' ' + 'DUM_' + tmp_line[1]

            creation_itp[i] = creation_itp[i] + ' ' + tmp_line[2]

            creation_itp[i] = creation_itp[i] + ' ' + tmp_line[3]

            creation_itp[i] = creation_itp[i] + ' ' + tmp_line[4]

            creation_itp[i] = creation_itp[i] + ' ' + tmp_line[5]

            creation_itp[i] = creation_itp[i] + ' ' + '0.0'

            creation_itp[i] = creation_itp[i] + ' ' + tmp_line[7]

            creation_itp[i] = creation_itp[i] + ' ' + tmp_line[2]

            creation_itp[i] = creation_itp[i] + ' ' + tmp_line[6]

            creation_itp[i] = creation_itp[i] + ' ' + tmp_line[7]

            creation_itp[i] = creation_itp[i] + '\n'

    return creation_itp
