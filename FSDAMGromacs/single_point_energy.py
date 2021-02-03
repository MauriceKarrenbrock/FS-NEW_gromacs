# -*- coding: utf-8 -*-
"""calculate single point energies

it is useful when you want to concatenate different MD programs
"""

import PythonAuxiliaryFunctions.run as _run

import FSDAMGromacs.mdp_files as _mdp
import FSDAMGromacs.parse as _parse
import FSDAMGromacs.tpr_files as _tpr


def get_single_point_energies(trajectory,
                              topology,
                              gromacs_path='gmx',
                              omp_threads=None):
    """get the single point energies of a structure or trajectory with gromacs

    Parameters
    ------------
    trajectory : str or path
        any kind of file that is recognised by gmx mdrun -rerun
        like pdb gro trr xtc
    topology : str or path
        a gromacs topology file
    gromacs_path : str or path, optional, default=gmx
        the path to the gromacs executable
    omp_threads : int, optional
        the number of OpenMP threads that gromacs should use
        default=gromacs choses automatically

    Returns
    ------------
    single_point_energies : numpy.array
        an array with the value of the single point energy for each structure in the
        trajectory in the same order

    Raises
    -----------------
    ValueError
        if omp_threads <= 0 and omp_threads is not None
    """

    if omp_threads is None:

        omp_string = ''

    elif omp_threads > 0:

        omp_string = f' -ntomp {omp_threads} '

    else:

        raise ValueError(
            f'omp_threads must be bigger than zero or None not: {omp_threads}')

    mdp_file_name = 'single_point_energy.mdp'

    #make mdp_file
    input_dict = {'continuation': 'yes', 'nsteps': 0, 'constraints': 'none'}

    mdp_maker = _mdp.BasicMdpFileMixIn()

    mdp_lines = mdp_maker.create_basic_template(input_dict)

    with open(mdp_file_name, 'w') as f:

        for line in mdp_lines:

            f.write(line + '\n')

    tpr_file_name = 'single_point_energy.tpr'

    _tpr.make_tpr_file(mdp_file=mdp_file_name,
                       gro_file=trajectory,
                       top_file=topology,
                       tpr_file=tpr_file_name,
                       gromacs_path=gromacs_path,
                       shell=True)

    gromacs_command_string = (
        f'{gromacs_path} mdrun -rerun {trajectory} -s {tpr_file_name} -deffnm single_point_energy'
        + omp_string)

    _run.subprocess_run(
        gromacs_command_string,
        shell=True,
        error_string=
        'error during single point calculation through gmx mdrun -rerun')

    with open('tmp_file.txt', 'w') as f:

        f.write('Potential\n0\n')

    energy_file = 'single_point_energy.xvg'

    gromacs_command_string = (
        f'< tmp_file.txt {gromacs_path} energy '
        f'-f single_point_energy.edr -s {tpr_file_name} -o {energy_file}')

    _run.subprocess_run(gromacs_command_string,
                        shell=True,
                        error_string='error during gmx energy')

    single_point_energies = _parse.parse_big_gromacs_xvg_files(energy_file)

    single_point_energies = single_point_energies.transpose()

    single_point_energies = single_point_energies[1]

    return single_point_energies
