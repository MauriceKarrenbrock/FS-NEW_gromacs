# -*- coding: utf-8 -*-
# pylint: disable=missing-docstring
# pylint: disable=redefined-outer-name
# pylint: disable=no-self-use
# pylint: disable=protected-access
# pylint: disable=duplicate-code
# pylint: disable=too-many-lines
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################

import pytest

import FSDAMGromacs.mdp_files as mdp_files


class Testmake_free_energy_lines():
    def test_works(self):

        expected = [
            '; Free energy control stuff', 'free-energy              = yes',
            'init-lambda              = 1', 'delta-lambda             = -0.1',
            'couple-moltype           = ALC', 'couple-lambda0           =vdw',
            'couple-lambda1           =none', 'couple-intramol          =no',
            'sc-alpha                 = 0.0', 'sc-coul                  = no',
            'sc-sigma                 = 0.25', 'sc-power                 = 1',
            'nstdhdl                  = 100', 'separate-dhdl-file       = yes',
            ''
        ]

        output = mdp_files.make_free_energy_lines(condition_lambda0='vdw',
                                                  condition_lambda1='none',
                                                  alchemical_molecule='ALC',
                                                  lambda_step=-0.1,
                                                  starting_lambda=1,
                                                  couple_intramol='no',
                                                  sc_alpha=0.0,
                                                  sc_coul='no',
                                                  sc_sigma=0.25,
                                                  sc_power=1,
                                                  nstdhdl=100,
                                                  separate_dhdl_file='yes',
                                                  free_energy='yes')

        assert output == expected


class Testcreate_COMCOM_pulling_strings():
    def test_works_pbc_none(self):

        COM_pull_goups = ['DU1', 'Protein', 'DU2']

        harmonic_kappa = [['DU1', 'Protein', 120], ['DU1', 'DU2', 121],
                          ['Protein', 'DU2', 0]]

        pbc = None

        output = mdp_files.create_COMCOM_pulling_strings(
            COM_pull_goups=COM_pull_goups,
            harmonic_kappa=harmonic_kappa,
            pbc_atoms=pbc)

        expected_output = [
            ';COM PULLING', 'pull                     = yes',
            'pull-print-com           = yes', 'pull-print-components    = no',
            f'pull-ncoords            = {len(COM_pull_goups)-1}',
            'pull-nstxout            = 10',
            f'pull-ngroups            = {len(COM_pull_goups)}',
            f'pull-group1-name        = {COM_pull_goups[0]}\n' +
            f'pull-group2-name        = {COM_pull_goups[1]}\n' +
            f'pull-group3-name        = {COM_pull_goups[2]}\n',
            'pull-pbc-ref-prev-step-com  = yes',
            'pull-coord1-geometry    = distance',
            'pull-coord1-type        = umbrella',
            'pull-coord1-dim         = Y Y Y', 'pull-coord1-groups      = 1 2',
            'pull-coord1-start       = yes', 'pull-coord1-init       = 0.0',
            'pull-coord1-rate       = 0', 'pull-coord1-k          = 120',
            'pull-coord2-geometry    = distance',
            'pull-coord2-type        = umbrella',
            'pull-coord2-dim         = Y Y Y', 'pull-coord2-groups      = 1 3',
            'pull-coord2-start       = yes', 'pull-coord2-init       = 0.0',
            'pull-coord2-rate       = 0', 'pull-coord2-k          = 121'
        ]

        assert output == expected_output

    def test_works_given_pbc(self):

        COM_pull_goups = ['DU1', 'Protein', 'DU2']

        harmonic_kappa = [['DU1', 'Protein', 120], ['DU1', 'DU2', 121],
                          ['Protein', 'DU2', 0]]

        pbc = (124, 0, 0)

        output = mdp_files.create_COMCOM_pulling_strings(
            COM_pull_goups=COM_pull_goups,
            harmonic_kappa=harmonic_kappa,
            pbc_atoms=pbc)

        expected_output = [
            ';COM PULLING',
            'pull                     = yes',
            'pull-print-com           = yes',
            'pull-print-components    = no',
            f'pull-ncoords            = {len(COM_pull_goups)-1}',
            'pull-nstxout            = 10',
            f'pull-ngroups            = {len(COM_pull_goups)}',
            f'pull-group1-name        = {COM_pull_goups[0]}\n' +
            f'pull-group2-name        = {COM_pull_goups[1]}\n' +
            f'pull-group3-name        = {COM_pull_goups[2]}\n',
            'pull-pbc-ref-prev-step-com  = yes',
            'pull-coord1-geometry    = distance',
            'pull-coord1-type        = umbrella',
            'pull-coord1-dim         = Y Y Y',
            'pull-coord1-groups      = 1 2',
            'pull-coord1-start       = yes',
            'pull-coord1-init       = 0.0',
            'pull-coord1-rate       = 0',
            'pull-coord1-k          = 120',
            'pull-group1-pbcatom     = 124',
            'pull-coord2-geometry    = distance',
            'pull-coord2-type        = umbrella',
            'pull-coord2-dim         = Y Y Y',
            'pull-coord2-groups      = 1 3',
            'pull-coord2-start       = yes',
            'pull-coord2-init       = 0.0',
            'pull-coord2-rate       = 0',
            'pull-coord2-k          = 121',
            'pull-group2-pbcatom     = 0',
            'pull-group3-pbcatom     = 0',
        ]

        assert output == expected_output


class TestMdpFile():
    def test_init(self):

        mdp_file = 'mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = ['a', 'b', 'c']

        instance = mdp_files.MdpFile(mdp_file, alchemical_molecule,
                                     timestep_ps, number_of_steps, temperature,
                                     lambda_steps, COM_pull_goups)

        output = [
            instance.mdp_file, instance.alchemical_molecule,
            instance.timestep_ps, instance.number_of_steps,
            instance.temperature, instance.lambda_steps,
            instance.COM_pull_goups, instance._template
        ]

        expected = [
            mdp_file + '.mdp', alchemical_molecule, timestep_ps,
            number_of_steps, temperature, lambda_steps, COM_pull_goups, []
        ]

        assert output == expected

    def test__create_free_energy_strings(self, mocker):

        mdp_file = 'mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = None

        instance = mdp_files.MdpFile(mdp_file, alchemical_molecule,
                                     timestep_ps, number_of_steps, temperature,
                                     lambda_steps, COM_pull_goups)

        m_free = mocker.patch('FSDAMGromacs.mdp_files.make_free_energy_lines',
                              return_value=['A', 'B'])

        output = instance._create_free_energy_strings(condition_lambda0='vdw',
                                                      condition_lambda1='none',
                                                      starting_lambda=1,
                                                      couple_intramol='no',
                                                      sc_alpha=0.0,
                                                      sc_coul='no',
                                                      sc_sigma=0.25,
                                                      sc_power=1,
                                                      nstdhdl=100,
                                                      separate_dhdl_file='yes',
                                                      free_energy='yes')

        assert output == ['A', 'B']

        m_free.assert_called_once_with(condition_lambda0='vdw',
                                       condition_lambda1='none',
                                       alchemical_molecule='alc',
                                       lambda_step=None,
                                       starting_lambda=1,
                                       couple_intramol='no',
                                       sc_alpha=0.0,
                                       sc_coul='no',
                                       sc_sigma=0.25,
                                       sc_power=1,
                                       nstdhdl=100,
                                       separate_dhdl_file='yes',
                                       free_energy='yes')

    def test__create_COMCOM_pulling_strings_None(self):

        mdp_file = 'mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = None

        instance = mdp_files.MdpFile(mdp_file, alchemical_molecule,
                                     timestep_ps, number_of_steps, temperature,
                                     lambda_steps, COM_pull_goups)

        assert instance._create_COMCOM_pulling_strings() == ['']

    def test__create_COMCOM_pulling_strings_empty_list(self):

        mdp_file = 'mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = []

        instance = mdp_files.MdpFile(mdp_file, alchemical_molecule,
                                     timestep_ps, number_of_steps, temperature,
                                     lambda_steps, COM_pull_goups)

        assert instance._create_COMCOM_pulling_strings() == ['']

    def test__create_COMCOM_pulling_strings_works(self, mocker):

        mdp_file = 'mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = ['DU1', 'Protein', 'DU2']
        harmonic_kappa = [['DU1', 'Protein', 120], ['DU1', 'DU2', 121],
                          ['Protein', 'DU2', 0]]

        m_COM = mocker.patch(
            'FSDAMGromacs.mdp_files.create_COMCOM_pulling_strings',
            return_value=['A', 'B'])

        instance = mdp_files.MdpFile(mdp_file, alchemical_molecule,
                                     timestep_ps, number_of_steps, temperature,
                                     lambda_steps, COM_pull_goups,
                                     harmonic_kappa)

        assert instance._create_COMCOM_pulling_strings() == ['A', 'B']

        m_COM.assert_called_once_with(COM_pull_goups=COM_pull_goups,
                                      harmonic_kappa=harmonic_kappa,
                                      pbc_atoms=None)

    @pytest.mark.parametrize('test_type, harmonic_kappa',
                             [('harmonic_kappa None', None),
                              ('harmonic_kappa []', [])])
    def test__create_COMCOM_pulling_strings_raises_valueerror(
            self, test_type, harmonic_kappa):

        print('Logging test type for visibility: ' + test_type)

        mdp_file = 'mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = ['DU1', 'Protein', 'DU2']

        instance = mdp_files.MdpFile(mdp_file, alchemical_molecule,
                                     timestep_ps, number_of_steps, temperature,
                                     lambda_steps, COM_pull_goups,
                                     harmonic_kappa)

        with pytest.raises(ValueError):

            instance._create_COMCOM_pulling_strings()

    def test_get_template(self):

        mdp_file = 'mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = ['DU1', 'Protein', 'DU2']

        instance = mdp_files.MdpFile(mdp_file, alchemical_molecule,
                                     timestep_ps, number_of_steps, temperature,
                                     lambda_steps, COM_pull_goups)

        instance._template = ['AAAAA', 'BBBBB\n', 'ccccc']

        expected = ['AAAAA\n', 'BBBBB\n', 'ccccc\n']

        assert instance._get_template() == expected

    def test_write_template(self, mocker):

        mocked_write = mocker.patch(
            'PythonAuxiliaryFunctions.files_IO.write_file.write_file')

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = ['DU1', 'Protein', 'DU2']

        instance = mdp_files.MdpFile(mdp_file, alchemical_molecule,
                                     timestep_ps, number_of_steps, temperature,
                                     lambda_steps, COM_pull_goups)

        instance._write_template(['template'])

        mocked_write.assert_called_once_with(['template'], mdp_file)

    def test_execute(self, mocker):

        m_hook = mocker.patch.object(mdp_files.MdpFile, '_hook')
        m_get = mocker.patch.object(mdp_files.MdpFile, '_get_template')
        m_write = mocker.patch.object(mdp_files.MdpFile, '_write_template')

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = ['DU1', 'Protein', 'DU2']

        instance = mdp_files.MdpFile(mdp_file, alchemical_molecule,
                                     timestep_ps, number_of_steps, temperature,
                                     lambda_steps, COM_pull_goups)

        instance.execute()

        m_hook.assert_called_once()
        m_get.assert_called_once()
        m_write.assert_called_once()


class TestAnnihilateVdwMdpBoundState():
    def test__create_template(self, mocker):

        mocked_COM = mocker.patch.object(mdp_files.AnnihilateVdwMdpBoundState,
                                         '_create_COMCOM_pulling_strings',
                                         return_value=['COM'])

        mocked_free = mocker.patch.object(mdp_files.AnnihilateVdwMdpBoundState,
                                          '_create_free_energy_strings',
                                          return_value=['FREE'])

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = None

        instance = mdp_files.AnnihilateVdwMdpBoundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._create_template()

        mocked_COM.assert_called_once()

        mocked_free.assert_called_once()

        expected_mdp = [
            '; VARIOUS PREPROCESSING OPTIONS',
            '; Preprocessor information: use cpp syntax.',
            '; e.g.: -I/home/joe/doe -I/home/mary/roe',
            'include                  =',
            '; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)',
            'define                   =', '', '; RUN CONTROL PARAMETERS',
            'integrator               = md', '; Start time and timestep in ps',
            'tinit                    = 0',
            f'dt                       = {timestep_ps}',
            f'nsteps                   = {number_of_steps}',
            '; For exact run continuation or redoing part of a run',
            'init-step                = 0',
            '; Part index is updated automatically on checkpointing (keeps files separate)',
            'simulation-part          = 1',
            '; mode for center of mass motion removal',
            'comm-mode                = Linear',
            '; number of steps for center of mass motion removal',
            'nstcomm                  = 100',
            '; group(s) for center of mass motion removal',
            'comm-grps                =', '',
            '; TEST PARTICLE INSERTION OPTIONS',
            'rtpi                     = 0.05', '', '; OUTPUT CONTROL OPTIONS',
            '; Output frequency for coords (x), velocities (v) and forces (f)',
            'nstxout                  = 10000',
            'nstvout                  = 10000',
            'nstfout                  = 10000',
            '; Output frequency for energies to log file and energy file',
            'nstlog                   = 1000',
            'nstcalcenergy            = 100',
            'nstenergy                = 1000',
            '; Output frequency and precision for .xtc file',
            'nstxtcout                = 2000',
            'xtc-precision            = 1000',
            '; This selects the subset of atoms for the .xtc file. You can',
            '; select multiple groups. By default all atoms will be written.',
            'xtc-grps                 =', '; Selection of energy groups',
            'energygrps               = System', '',
            '; NEIGHBORSEARCHING PARAMETERS',
            '; cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)',
            '; nblist update frequency', 'cutoff-scheme            = Verlet',
            'nstlist                  = 20',
            'verlet-buffer-tolerance  = 0.0001',
            '; ns algorithm (simple or grid)',
            'ns_type                  = grid',
            '; Periodic boundary conditions: xyz, no, xy',
            'pbc                      = xyz', 'periodic-molecules       = no',
            '; Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,',
            '; a value of -1 means: use rlist', '; nblist cut-off',
            'rlist                    = 1.0',
            '; long-range cut-off for switched potentials',
            'rlistlong                = -1', '',
            '; OPTIONS FOR ELECTROSTATICS AND VDW',
            '; Method for doing electrostatics',
            'coulombtype              = PME', 'rcoulomb-switch          = 0',
            'rcoulomb                 = 1.0',
            '; Relative dielectric constant for the medium and the reaction field',
            'epsilon-r                = 1', 'epsilon-rf               = 0',
            '; Method for doing Van der Waals',
            'vdw-type                 = Cut-off', '; cut-off lengths',
            'rvdw-switch              = 0', 'rvdw                     = 1.0',
            '; Apply long range dispersion corrections for Energy and Pressure',
            'DispCorr                 = EnerPres',
            '; Extension of the potential lookup tables beyond the cut-off',
            'table-extension          = 1',
            '; Separate tables between energy group pairs',
            'energygrp-table          =',
            '; Spacing for the PME/PPPM FFT grid',
            'fourierspacing           = 0.1',
            '; FFT grid size, when a value is 0 fourierspacing will be used',
            'fourier-nx               = 0', 'fourier-ny               = 0',
            'fourier-nz               = 0', '; EWALD/PME/PPPM parameters',
            'pme-order                = 4', 'ewald-rtol               = 1e-05',
            'ewald-geometry           = 3d', 'epsilon-surface          =',
            'optimize-fft             = no', '',
            '; IMPLICIT SOLVENT ALGORITHM', 'implicit-solvent         = No',
            '', '; OPTIONS FOR WEAK COUPLING ALGORITHMS',
            '; Temperature coupling', 'tcoupl                   = v-rescale',
            'nsttcouple               = -1', 'nh-chain-length          = 1',
            '; Groups to couple separately',
            'tc-grps                  = System',
            '; Time constant (ps) and reference temperature (K)',
            'tau-t                    = 0.2',
            f'ref-t                    = {temperature}', '; pressure coupling',
            'pcoupl                   = Parrinello-Rahman',
            'pcoupltype               = Isotropic',
            'nstpcouple               = -1',
            '; Time constant (ps), compressibility (1/bar) and reference P (bar)',
            'tau-p                    = 1.0',
            'compressibility          = 4.6e-5',
            'ref-p                    = 1',
            '; Scaling of reference coordinates, No, All or COM',
            'refcoord-scaling         = COM', '',
            '; GENERATE VELOCITIES FOR STARTUP RUN',
            'gen-vel                  = no', 'gen-temp                 = 500',
            'gen-seed                 = 173529', '', '; OPTIONS FOR BONDS',
            'constraints              = all-bonds',
            '; Type of constraint algorithm',
            'constraint-algorithm     = Lincs',
            '; Do not constrain the start configuration',
            'continuation             = no',
            '; Use successive overrelaxation to reduce the number of shake iterations',
            'Shake-SOR                = no', '; Relative tolerance of shake',
            'shake-tol                = 0.00001',
            '; Highest order in the expansion of the constraint coupling matrix',
            'lincs-order              = 5',
            '; Number of iterations in the final step of LINCS. 1 is fine for',
            '; normal simulations, but use 2 to conserve energy in NVE runs.',
            '; For energy minimization with constraints it should be 4 to 8.',
            'lincs-iter               = 2',
            '; Lincs will write a warning to the stderr if in one step a bond',
            '; rotates over more degrees than',
            'lincs-warnangle          = 30',
            '; Convert harmonic bonds to morse potentials',
            'morse                    = no', '', 'FREE', 'COM'
        ]

        assert instance._template == expected_mdp

    def test_hook_lambda_None(self, mocker):

        mocked_template = \
            mocker.patch.object(mdp_files.AnnihilateVdwMdpBoundState, '_create_template')

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = None

        instance = mdp_files.AnnihilateVdwMdpBoundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._hook()

        mocked_template.assert_called_once()

        assert instance.lambda_steps == (-1. / number_of_steps)

    def test_hook_lambda_input(self, mocker):

        mocked_template = \
            mocker.patch.object(mdp_files.AnnihilateVdwMdpBoundState, '_create_template')

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = 33
        COM_pull_goups = None

        instance = mdp_files.AnnihilateVdwMdpBoundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._hook()

        mocked_template.assert_called_once()

        assert instance.lambda_steps == lambda_steps


class TestAnnihilateQMdpBoundState():
    def test__create_template(self, mocker):

        mocked_COM = mocker.patch.object(mdp_files.AnnihilateQMdpBoundState,
                                         '_create_COMCOM_pulling_strings',
                                         return_value=['COM'])

        mocked_free = mocker.patch.object(mdp_files.AnnihilateQMdpBoundState,
                                          '_create_free_energy_strings',
                                          return_value=['FREE'])

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = None

        instance = mdp_files.AnnihilateQMdpBoundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._create_template()

        mocked_COM.assert_called_once()
        mocked_free.assert_called_once()

        expected_mdp = [
            '; VARIOUS PREPROCESSING OPTIONS',
            '; Preprocessor information: use cpp syntax.',
            '; e.g.: -I/home/joe/doe -I/home/mary/roe',
            'include                  =',
            '; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)',
            'define                   =', '', '; RUN CONTROL PARAMETERS',
            'integrator               = md', '; Start time and timestep in ps',
            'tinit                    = 0',
            f'dt                       = {timestep_ps}',
            f'nsteps                   = {number_of_steps}',
            '; For exact run continuation or redoing part of a run',
            'init-step                = 0',
            '; Part index is updated automatically on checkpointing (keeps files separate)',
            'simulation-part          = 1',
            '; mode for center of mass motion removal',
            'comm-mode                = Linear',
            '; number of steps for center of mass motion removal',
            'nstcomm                  = 100',
            '; group(s) for center of mass motion removal',
            'comm-grps                =', '',
            '; TEST PARTICLE INSERTION OPTIONS',
            'rtpi                     = 0.05', '', '; OUTPUT CONTROL OPTIONS',
            '; Output frequency for coords (x), velocities (v) and forces (f)',
            'nstxout                  = 10000',
            'nstvout                  = 10000',
            'nstfout                  = 10000',
            '; Output frequency for energies to log file and energy file',
            'nstlog                   = 1000',
            'nstcalcenergy            = 100',
            'nstenergy                = 1000',
            '; Output frequency and precision for .xtc file',
            'nstxtcout                = 2000',
            'xtc-precision            = 1000',
            '; This selects the subset of atoms for the .xtc file. You can',
            '; select multiple groups. By default all atoms will be written.',
            'xtc-grps                 =', '; Selection of energy groups',
            'energygrps               = System', '',
            '; NEIGHBORSEARCHING PARAMETERS',
            '; cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)',
            '; nblist update frequency', 'cutoff-scheme            = Verlet',
            'nstlist                  = 20',
            'verlet-buffer-tolerance  = 0.0001',
            '; ns algorithm (simple or grid)',
            'ns_type                  = grid',
            '; Periodic boundary conditions: xyz, no, xy',
            'pbc                      = xyz', 'periodic-molecules       = no',
            '; Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,',
            '; a value of -1 means: use rlist', '; nblist cut-off',
            'rlist                    = 1.0',
            '; long-range cut-off for switched potentials',
            'rlistlong                = -1', '',
            '; OPTIONS FOR ELECTROSTATICS AND VDW',
            '; Method for doing electrostatics',
            'coulombtype              = PME', 'rcoulomb-switch          = 0',
            'rcoulomb                 = 1.0',
            '; Relative dielectric constant for the medium and the reaction field',
            'epsilon-r                = 1', 'epsilon-rf               = 0',
            '; Method for doing Van der Waals',
            'vdw-type                 = Cut-off', '; cut-off lengths',
            'rvdw-switch              = 0', 'rvdw                     = 1.0',
            '; Apply long range dispersion corrections for Energy and Pressure',
            'DispCorr                 = EnerPres',
            '; Extension of the potential lookup tables beyond the cut-off',
            'table-extension          = 1',
            '; Separate tables between energy group pairs',
            'energygrp-table          =',
            '; Spacing for the PME/PPPM FFT grid',
            'fourierspacing           = 0.1',
            '; FFT grid size, when a value is 0 fourierspacing will be used',
            'fourier-nx               = 0', 'fourier-ny               = 0',
            'fourier-nz               = 0', '; EWALD/PME/PPPM parameters',
            'pme-order                = 4', 'ewald-rtol               = 1e-05',
            'ewald-geometry           = 3d', 'epsilon-surface          =',
            'optimize-fft             = no', '',
            '; IMPLICIT SOLVENT ALGORITHM', 'implicit-solvent         = No',
            '', '; OPTIONS FOR WEAK COUPLING ALGORITHMS',
            '; Temperature coupling', 'tcoupl                   = v-rescale',
            'nsttcouple               = -1', 'nh-chain-length          = 1',
            '; Groups to couple separately',
            'tc-grps                  = System',
            '; Time constant (ps) and reference temperature (K)',
            'tau-t                    = 0.2',
            f'ref-t                    = {temperature}', '; pressure coupling',
            'pcoupl                   = Parrinello-Rahman',
            'pcoupltype               = Isotropic',
            'nstpcouple               = -1',
            '; Time constant (ps), compressibility (1/bar) and reference P (bar)',
            'tau-p                    = 1.0',
            'compressibility          = 4.6e-5',
            'ref-p                    = 1',
            '; Scaling of reference coordinates, No, All or COM',
            'refcoord-scaling         = COM', '',
            '; GENERATE VELOCITIES FOR STARTUP RUN',
            'gen-vel                  = no', 'gen-temp                 = 500',
            'gen-seed                 = 173529', '', '; OPTIONS FOR BONDS',
            'constraints              = all-bonds',
            '; Type of constraint algorithm',
            'constraint-algorithm     = Lincs',
            '; Do not constrain the start configuration',
            'continuation             = no',
            '; Use successive overrelaxation to reduce the number of shake iterations',
            'Shake-SOR                = no', '; Relative tolerance of shake',
            'shake-tol                = 0.00001',
            '; Highest order in the expansion of the constraint coupling matrix',
            'lincs-order              = 5',
            '; Number of iterations in the final step of LINCS. 1 is fine for',
            '; normal simulations, but use 2 to conserve energy in NVE runs.',
            '; For energy minimization with constraints it should be 4 to 8.',
            'lincs-iter               = 2',
            '; Lincs will write a warning to the stderr if in one step a bond',
            '; rotates over more degrees than',
            'lincs-warnangle          = 30',
            '; Convert harmonic bonds to morse potentials',
            'morse                    = no', '', 'FREE', 'COM'
        ]

        assert instance._template == expected_mdp

    def test_hook_lambda_None(self, mocker):

        mocked_template = \
            mocker.patch.object(mdp_files.AnnihilateQMdpBoundState, '_create_template')

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = None

        instance = mdp_files.AnnihilateQMdpBoundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._hook()

        mocked_template.assert_called_once()

        assert instance.lambda_steps == (-1. / number_of_steps)

    def test_hook_lambda_input(self, mocker):

        mocked_template = \
            mocker.patch.object(mdp_files.AnnihilateQMdpBoundState, '_create_template')

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = 33
        COM_pull_goups = None

        instance = mdp_files.AnnihilateQMdpBoundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._hook()

        mocked_template.assert_called_once()

        assert instance.lambda_steps == lambda_steps


class TestCreateVdwMdpUnboundState():
    def test__create_template(self, mocker):

        mocked_COM = mocker.patch.object(mdp_files.CreateVdwMdpUnboundState,
                                         '_create_COMCOM_pulling_strings',
                                         return_value=['COM'])

        mocked_free = mocker.patch.object(mdp_files.CreateVdwMdpUnboundState,
                                          '_create_free_energy_strings',
                                          return_value=['FREE'])

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = None

        instance = mdp_files.CreateVdwMdpUnboundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._create_template()

        mocked_COM.assert_called_once()
        mocked_free.assert_called_once()

        expected_mdp = [
            '; VARIOUS PREPROCESSING OPTIONS',
            '; Preprocessor information: use cpp syntax.',
            '; e.g.: -I/home/joe/doe -I/home/mary/roe',
            'include                  =',
            '; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)',
            'define                   =', '', '; RUN CONTROL PARAMETERS',
            'integrator               = md', '; Start time and timestep in ps',
            'tinit                    = 0',
            f'dt                       = {timestep_ps}',
            f'nsteps                   = {number_of_steps}',
            '; For exact run continuation or redoing part of a run',
            'init-step                = 0',
            '; Part index is updated automatically on checkpointing (keeps files separate)',
            'simulation-part          = 1',
            '; mode for center of mass motion removal',
            'comm-mode                = Linear',
            '; number of steps for center of mass motion removal',
            'nstcomm                  = 100',
            '; group(s) for center of mass motion removal',
            'comm-grps                =', '',
            '; TEST PARTICLE INSERTION OPTIONS',
            'rtpi                     = 0.05', '', '; OUTPUT CONTROL OPTIONS',
            '; Output frequency for coords (x), velocities (v) and forces (f)',
            'nstxout                  = 10000',
            'nstvout                  = 10000',
            'nstfout                  = 10000',
            '; Output frequency for energies to log file and energy file',
            'nstlog                   = 500', 'nstcalcenergy            = 100',
            'nstenergy                = 1000',
            '; Output frequency and precision for .xtc file',
            'nstxtcout                = 2000',
            'xtc-precision            = 1000',
            '; This selects the subset of atoms for the .xtc file. You can',
            '; select multiple groups. By default all atoms will be written.',
            'xtc-grps                 =', '; Selection of energy groups',
            'energygrps               = System', '',
            '; NEIGHBORSEARCHING PARAMETERS',
            '; cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)',
            '; nblist update frequency', 'cutoff-scheme            = Verlet',
            'nstlist                  = 20',
            'verlet-buffer-tolerance  = 0.0001',
            '; ns algorithm (simple or grid)',
            'ns_type                  = grid',
            '; Periodic boundary conditions: xyz, no, xy',
            'pbc                      = xyz', 'periodic-molecules       = no',
            '; Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,',
            '; a value of -1 means: use rlist', '; nblist cut-off',
            'rlist                    = 1.0',
            '; long-range cut-off for switched potentials',
            'rlistlong                = -1', '',
            '; OPTIONS FOR ELECTROSTATICS AND VDW',
            '; Method for doing electrostatics',
            'coulombtype              = PME', 'rcoulomb-switch          = 0',
            'rcoulomb                 = 1.0',
            '; Relative dielectric constant for the medium and the reaction field',
            'epsilon-r                = 1', 'epsilon-rf               = 0',
            '; Method for doing Van der Waals',
            'vdw-type                 = Cut-off', '; cut-off lengths',
            'rvdw-switch              = 0', 'rvdw                     = 1.0',
            '; Apply long range dispersion corrections for Energy and Pressure',
            'DispCorr                 = EnerPres',
            '; Extension of the potential lookup tables beyond the cut-off',
            'table-extension          = 1',
            '; Separate tables between energy group pairs',
            'energygrp-table          =',
            '; Spacing for the PME/PPPM FFT grid',
            'fourierspacing           = 0.1',
            '; FFT grid size, when a value is 0 fourierspacing will be used',
            'fourier-nx               = 0', 'fourier-ny               = 0',
            'fourier-nz               = 0', '; EWALD/PME/PPPM parameters',
            'pme-order                = 4', 'ewald-rtol               = 1e-05',
            'ewald-geometry           = 3d', 'epsilon-surface          =',
            'optimize-fft             = no', '',
            '; IMPLICIT SOLVENT ALGORITHM', 'implicit-solvent         = No',
            '', '; OPTIONS FOR WEAK COUPLING ALGORITHMS',
            '; Temperature coupling', 'tcoupl                   = v-rescale',
            'nsttcouple               = -1', 'nh-chain-length          = 1',
            '; Groups to couple separately',
            'tc-grps                  = System',
            '; Time constant (ps) and reference temperature (K)',
            'tau-t                    = 0.2',
            f'ref-t                    = {temperature}', '; pressure coupling',
            'pcoupl                   = Parrinello-Rahman',
            'pcoupltype               = Isotropic',
            'nstpcouple               = -1',
            '; Time constant (ps), compressibility (1/bar) and reference P (bar)',
            'tau-p                    = 1.0',
            'compressibility          = 4.6e-5',
            'ref-p                    = 1',
            '; Scaling of reference coordinates, No, All or COM',
            'refcoord-scaling         = COM', '',
            '; GENERATE VELOCITIES FOR STARTUP RUN',
            'gen-vel                  = no', 'gen-temp                 = 500',
            'gen-seed                 = 173529', '', '; OPTIONS FOR BONDS',
            'constraints              = all-bonds',
            '; Type of constraint algorithm',
            'constraint-algorithm     = Lincs',
            '; Do not constrain the start configuration',
            'continuation             = no',
            '; Use successive overrelaxation to reduce the number of shake iterations',
            'Shake-SOR                = no', '; Relative tolerance of shake',
            'shake-tol                = 0.00001',
            '; Highest order in the expansion of the constraint coupling matrix',
            'lincs-order              = 5',
            '; Number of iterations in the final step of LINCS. 1 is fine for',
            '; normal simulations, but use 2 to conserve energy in NVE runs.',
            '; For energy minimization with constraints it should be 4 to 8.',
            'lincs-iter               = 2',
            '; Lincs will write a warning to the stderr if in one step a bond',
            '; rotates over more degrees than',
            'lincs-warnangle          = 30',
            '; Convert harmonic bonds to morse potentials',
            'morse                    = no', '', 'FREE', 'COM'
        ]

        assert instance._template == expected_mdp

    def test_hook_lambda_None(self, mocker):

        mocked_template = \
            mocker.patch.object(mdp_files.CreateVdwMdpUnboundState, '_create_template')

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = None

        instance = mdp_files.CreateVdwMdpUnboundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._hook()

        mocked_template.assert_called_once()

        assert instance.lambda_steps == (1. / number_of_steps)

    def test_hook_lambda_input(self, mocker):

        mocked_template = \
            mocker.patch.object(mdp_files.CreateVdwMdpUnboundState, '_create_template')

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = 33
        COM_pull_goups = None

        instance = mdp_files.CreateVdwMdpUnboundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._hook()

        mocked_template.assert_called_once()

        assert instance.lambda_steps == lambda_steps


class TestCreateQMdpUnboundState():
    def test__create_template(self, mocker):

        mocked_COM = mocker.patch.object(mdp_files.CreateQMdpUnboundState,
                                         '_create_COMCOM_pulling_strings',
                                         return_value=['COM'])

        mocked_free = mocker.patch.object(mdp_files.CreateQMdpUnboundState,
                                          '_create_free_energy_strings',
                                          return_value=['FREE'])

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = None

        instance = mdp_files.CreateQMdpUnboundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._create_template()

        mocked_COM.assert_called_once()
        mocked_free.assert_called_once()

        expected_mdp = [
            '; VARIOUS PREPROCESSING OPTIONS',
            '; Preprocessor information: use cpp syntax.',
            '; e.g.: -I/home/joe/doe -I/home/mary/roe',
            'include                  =',
            '; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)',
            'define                   =', '', '; RUN CONTROL PARAMETERS',
            'integrator               = md', '; Start time and timestep in ps',
            'tinit                    = 0',
            f'dt                       = {timestep_ps}',
            f'nsteps                   = {number_of_steps}',
            '; For exact run continuation or redoing part of a run',
            'init-step                = 0',
            '; Part index is updated automatically on checkpointing (keeps files separate)',
            'simulation-part          = 1',
            '; mode for center of mass motion removal',
            'comm-mode                = Linear',
            '; number of steps for center of mass motion removal',
            'nstcomm                  = 100',
            '; group(s) for center of mass motion removal',
            'comm-grps                =', '',
            '; TEST PARTICLE INSERTION OPTIONS',
            'rtpi                     = 0.05', '', '; OUTPUT CONTROL OPTIONS',
            '; Output frequency for coords (x), velocities (v) and forces (f)',
            'nstxout                  = 10000',
            'nstvout                  = 10000',
            'nstfout                  = 10000',
            '; Output frequency for energies to log file and energy file',
            'nstlog                   = 500', 'nstcalcenergy            = 100',
            'nstenergy                = 1000',
            '; Output frequency and precision for .xtc file',
            'nstxtcout                = 2000',
            'xtc-precision            = 1000',
            '; This selects the subset of atoms for the .xtc file. You can',
            '; select multiple groups. By default all atoms will be written.',
            'xtc-grps                 =', '; Selection of energy groups',
            'energygrps               = System', '',
            '; NEIGHBORSEARCHING PARAMETERS',
            '; cut-off scheme (group: using charge groups, Verlet: particle based cut-offs)',
            '; nblist update frequency', 'cutoff-scheme            = Verlet',
            'nstlist                  = 20',
            'verlet-buffer-tolerance  = 0.0001',
            '; ns algorithm (simple or grid)',
            'ns_type                  = grid',
            '; Periodic boundary conditions: xyz, no, xy',
            'pbc                      = xyz', 'periodic-molecules       = no',
            '; Allowed energy drift due to the Verlet buffer in kJ/mol/ps per atom,',
            '; a value of -1 means: use rlist', '; nblist cut-off',
            'rlist                    = 1.0',
            '; long-range cut-off for switched potentials',
            'rlistlong                = -1', '',
            '; OPTIONS FOR ELECTROSTATICS AND VDW',
            '; Method for doing electrostatics',
            'coulombtype              = PME', 'rcoulomb-switch          = 0',
            'rcoulomb                 = 1.0',
            '; Relative dielectric constant for the medium and the reaction field',
            'epsilon-r                = 1', 'epsilon-rf               = 0',
            '; Method for doing Van der Waals',
            'vdw-type                 = Cut-off', '; cut-off lengths',
            'rvdw-switch              = 0', 'rvdw                     = 1.0',
            '; Apply long range dispersion corrections for Energy and Pressure',
            'DispCorr                 = EnerPres',
            '; Extension of the potential lookup tables beyond the cut-off',
            'table-extension          = 1',
            '; Separate tables between energy group pairs',
            'energygrp-table          =',
            '; Spacing for the PME/PPPM FFT grid',
            'fourierspacing           = 0.1',
            '; FFT grid size, when a value is 0 fourierspacing will be used',
            'fourier-nx               = 0', 'fourier-ny               = 0',
            'fourier-nz               = 0', '; EWALD/PME/PPPM parameters',
            'pme-order                = 4', 'ewald-rtol               = 1e-05',
            'ewald-geometry           = 3d', 'epsilon-surface          =',
            'optimize-fft             = no', '',
            '; IMPLICIT SOLVENT ALGORITHM', 'implicit-solvent         = No',
            '', '; OPTIONS FOR WEAK COUPLING ALGORITHMS',
            '; Temperature coupling', 'tcoupl                   = v-rescale',
            'nsttcouple               = -1', 'nh-chain-length          = 1',
            '; Groups to couple separately',
            'tc-grps                  = System',
            '; Time constant (ps) and reference temperature (K)',
            'tau-t                    = 0.2',
            f'ref-t                    = {temperature}', '; pressure coupling',
            'pcoupl                   = Parrinello-Rahman',
            'pcoupltype               = Isotropic',
            'nstpcouple               = -1',
            '; Time constant (ps), compressibility (1/bar) and reference P (bar)',
            'tau-p                    = 1.0',
            'compressibility          = 4.6e-5',
            'ref-p                    = 1',
            '; Scaling of reference coordinates, No, All or COM',
            'refcoord-scaling         = COM', '',
            '; GENERATE VELOCITIES FOR STARTUP RUN',
            'gen-vel                  = no', 'gen-temp                 = 500',
            'gen-seed                 = 173529', '', '; OPTIONS FOR BONDS',
            'constraints              = all-bonds',
            '; Type of constraint algorithm',
            'constraint-algorithm     = Lincs',
            '; Do not constrain the start configuration',
            'continuation             = no',
            '; Use successive overrelaxation to reduce the number of shake iterations',
            'Shake-SOR                = no', '; Relative tolerance of shake',
            'shake-tol                = 0.00001',
            '; Highest order in the expansion of the constraint coupling matrix',
            'lincs-order              = 5',
            '; Number of iterations in the final step of LINCS. 1 is fine for',
            '; normal simulations, but use 2 to conserve energy in NVE runs.',
            '; For energy minimization with constraints it should be 4 to 8.',
            'lincs-iter               = 2',
            '; Lincs will write a warning to the stderr if in one step a bond',
            '; rotates over more degrees than',
            'lincs-warnangle          = 30',
            '; Convert harmonic bonds to morse potentials',
            'morse                    = no', '', 'FREE', 'COM'
        ]

        assert instance._template == expected_mdp

    def test_hook_lambda_None(self, mocker):

        mocked_template = \
            mocker.patch.object(mdp_files.CreateQMdpUnboundState, '_create_template')

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = None
        COM_pull_goups = None

        instance = mdp_files.CreateQMdpUnboundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._hook()

        mocked_template.assert_called_once()

        assert instance.lambda_steps == (1. / number_of_steps)

    def test_hook_lambda_input(self, mocker):

        mocked_template = \
            mocker.patch.object(mdp_files.CreateQMdpUnboundState, '_create_template')

        mdp_file = 'mdp.mdp'
        alchemical_molecule = 'alc'
        timestep_ps = 0.002
        number_of_steps = 1000
        temperature = 297.20
        lambda_steps = 33
        COM_pull_goups = None

        instance = mdp_files.CreateQMdpUnboundState(
            mdp_file, alchemical_molecule, timestep_ps, number_of_steps,
            temperature, lambda_steps, COM_pull_goups)

        instance._hook()

        mocked_template.assert_called_once()

        assert instance.lambda_steps == lambda_steps
