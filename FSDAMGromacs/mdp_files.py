# -*- coding: utf-8 -*-
# pylint: disable=duplicate-code
# pylint: disable=too-many-lines
#############################################################
# Copyright (c) 2020-2020 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Classes to create the mdp files for NE alchemical transformation

Creation and annihilation of VdW (Van der Waals contribution)
and Q (charges)
"""

import PythonAuxiliaryFunctions.files_IO.write_file as write_file


def make_free_energy_lines(condition_lambda0,
                           condition_lambda1,
                           alchemical_molecule,
                           lambda_step,
                           starting_lambda=0,
                           couple_intramol='no',
                           sc_alpha=0.0,
                           sc_coul='no',
                           sc_sigma=0.25,
                           sc_power=1,
                           nstdhdl=100,
                           separate_dhdl_file='yes',
                           free_energy='yes'):
    """Get the free energy lines of a mdp_file

    Parameters
    ------------
    condition_lambda0 : str
        the condition of the sistem for lambda=0
        the options are 'vdw-q', 'vdw', 'none'
    condition_lambda1 : str
        the condition of the sistem for lambda=1
        the options are 'vdw-q', 'vdw', 'none'
    alchemical_molecule : str
        the residue name of the molecule that will be
        annihilated/created
    lambda_step : float
        how much lambda will increase/decrease each timestep
        (can be both positive and negative)
    starting_lambda : int, default=0
        the starting value of lambda (usually 1 or 0)
    couple_intramol : str, optional, default='no'
        check the mdp documentation on the gromacs website
        don't use this options if you don't know what you are doing
    sc_alpha : float, optional, default=0.0
        check the mdp documentation on the gromacs website
        don't use this options if you don't know what you are doing
    sc_coul : str, optional, default='no'
        check the mdp documentation on the gromacs website
        don't use this options if you don't know what you are doing
    sc_sigma : float, optional, default=0.25
        check the mdp documentation on the gromacs website
        don't use this options if you don't know what you are doing
    sc_power : int, optional, default=1
        check the mdp documentation on the gromacs website
        don't use this options if you don't know what you are doing
    nstdhdl : int, optional, default=100
        check the mdp documentation on the gromacs website
        don't use this options if you don't know what you are doing
    separate_dhdl_file : str, optional, default='yes'
        check the mdp documentation on the gromacs website
        don't use this options if you don't know what you are doing
    free_energy : str, optional, default='yes'
        check the mdp documentation on the gromacs website
        don't use this options if you don't know what you are doing

    Returns
    ------------
    list of str
        the lines of the mdp file (newline missing)

    Notes
    ----------
    A good convention but it's not compulsory is that during creation
    lambda goes from 0 -> 1 and during annihilation 1 -> 0
    """

    lines = [
        '; Free energy control stuff',
        f'free-energy              = {free_energy}',
        f'init-lambda              = {starting_lambda}',
        f'delta-lambda             = {lambda_step}',
        f'couple-moltype           = {alchemical_molecule}',
        f'couple-lambda0           ={condition_lambda0}',
        f'couple-lambda1           ={condition_lambda1}',
        f'couple-intramol          ={couple_intramol}',
        f'sc-alpha                 = {sc_alpha}',
        f'sc-coul                  = {sc_coul}',
        f'sc-sigma                 = {sc_sigma}',
        f'sc-power                 = {sc_power}',
        f'nstdhdl                  = {nstdhdl}',
        f'separate-dhdl-file       = {separate_dhdl_file}', ''
    ]

    return lines


def create_COMCOM_pulling_strings(COM_pull_goups,
                                  harmonic_kappa,
                                  pbc_atoms=None):
    """The strings for COM-COM pulling, in mdp file

    Parameters
    ------------
    COM_pull_goups : list of strings
        the list of gromacs groups (System,Protein,<ligand_resname>,...)
        that will have an harmonic COM-COM (center of mass) constrain
        Usually for unbound state no pulling is needed (default)
        for bound state it is ["Protein", "<ligand residue name>", "<dummy heavy atom>"]
    harmonic_kappa : list
        [ ["group_1", "group_2", harmonic_kappa_value], ... ] (str, str, float)
        it is a nested list containing the couple-couple harmonic kappa value
        for the umbrella COM-COM pulling, a good number may be 120
        if you don't want to groups to pull each other set kappa to 0
    pbc_atoms : iterable of int, optional
        an iterable containing the number of each nearest atom to the
        geometric center of the `COM_pull_goups` molecule
        must be as long as `COM_pull_goups`, for small molecules
        you can set it to zero (gromacs will guess it) for big ones (proteins)
        it must be given, if you keep it None it will gess it on any molecule

    Returns
    ----------
    list of strings
        missing newlines
    """

    pull_ngroups = len(COM_pull_goups)

    pull_ncoords = 0
    for i in range(pull_ngroups - 1, 0, -1):

        pull_ncoords += i

    pull_groups_name_number = {}
    pull_group_name = ''
    for i, group in enumerate(COM_pull_goups):
        pull_group_name += f'pull-group{i + 1}-name        = {group}\n'

        pull_groups_name_number[group] = i + 1

    pull_coord = []
    for i, couple in enumerate(harmonic_kappa):

        #don't write a pull for things with zero harmonic constant
        if couple[2] not in (0, 0., '0', '0.'):

            pull_coord += [
                f'pull-coord{i + 1}-geometry    = distance',
                f'pull-coord{i + 1}-type        = umbrella',
                f'pull-coord{i + 1}-dim         = Y Y Y',

                f'pull-coord{i + 1}-groups      = ' + \
                f'{pull_groups_name_number[couple[0]]} {pull_groups_name_number[couple[1]]}',

                f'pull-coord{i + 1}-start       = yes',
                f'pull-coord{i + 1}-init       = 0.0',
                f'pull-coord{i + 1}-rate       = 0',
                f'pull-coord{i + 1}-k          = {couple[2]}'
            ]

        else:

            pull_ncoords -= 1

        if pbc_atoms is not None:

            pull_coord += [f'pull-group{i + 1}-pbcatom     = {pbc_atoms[i]}']

    COM_pulling_strings = [
        ';COM PULLING', 'pull                     = yes',
        'pull-print-com           = yes', 'pull-print-components    = no',
        f'pull-ncoords            = {pull_ncoords}',
        'pull-nstxout            = 10',
        f'pull-ngroups            = {pull_ngroups}', f'{pull_group_name}',
        'pull-pbc-ref-prev-step-com  = yes'
    ]

    COM_pulling_strings += pull_coord

    return COM_pulling_strings


class MdpFile(object):
    """Super class for MDP file creation 4 alchemical transformations

    Interface for MDP creation classes for alchemical
    transformations

    Attributes
    -----------
    _template : list
        Private, list of strings of the template
        that will be written on the mdp file
        ("\n" new line caracters may be missing)

    Parameters
    ----------
    mdp_file : str
        the name of the output mdp file (or the path)
    alchemical_molecule : str
        the residue name (resname) or gromacs group
        (System, Protein, Other, ...) that identifies
        the alchemical atoms/molecules
    timestep_ps : float
        the timestep of the MD run (ps)
        default 0.001
    number_of_steps : int
        number of MD steps
        default 1000000
    temperature : float
        temperature in Kelvin (K)
    lambda_steps : float, optional
        DON'T USE IF NOT SURE OF WHAT YOU ARE DOING
        default will be good 99% of the times.
        how much for each step should the achemical
        parameter (lambda) increase/decrease,
        default (if None) each subclass will have its default
        usually `lambda_steps`= +/- 1/`number_of_steps`
        in order to have a complete creation/annihilation
        of the alchemical molecule (it is subclass dependent)
    COM_pull_goups : list of strings
        the list of gromacs groups (System,Protein,<ligand_resname>,...)
        that will have an harmonic COM-COM (center of mass) constrain
        Usually for unbound state no pulling is needed (default)
        for bound state it is ["Protein", "<ligand residue name>", "<dummy heavy atom>"]
    harmonic_kappa : list
        [ ["group_1", "group_2", harmonic_kappa_value], ... ] (str, str, float)
        it is a nested list containing the couple-couple harmonic kappa value
        for the umbrella COM-COM pulling, a good numbe may be 120
        if you don't want to groups to pull each other set kappa to 0
    pbc_atoms : iterable of int, optional
        an iterable containing the number of each nearest atom to the
        geometric center of the `COM_pull_goups` molecule
        must be as long as `COM_pull_goups`, for small molecules
        you can set it to zero (gromacs will guess it) for big ones (proteins)
        it must be given, if you keep it None it will gess it on any molecule


    Methods
    ---------
    execute()
        the only pubblic method, writes the mdp file
        with the right template

    Notes
    -----------
    `COM_pull_goups` in the bound state may need a dummy heavy atom
    (no LJ nor Q but infinite mass) because
    COM COM pulling in gromacs crashes because of poor PBC implemetation if the protein
    crosses the box

    TIPS for input values:
    This are some good values if you are creating/annihilating a small
    organic molecule

    * vdw creation:
        * `timestep_ps` = 0.001
        * `number_of_steps` = 160000
    * q creation:
        * `timestep_ps` = 0.001
        * `number_of_steps` = 160000
    * vdw annihilation:
        * `timestep_ps` = 0.0015
        * `number_of_steps` = 1000000
    * q annihilation:
        * `timestep_ps` = 0.0015
        * `number_of_steps` = 500000
    """
    def __init__(self,
                 mdp_file,
                 alchemical_molecule,
                 timestep_ps=0.001,
                 number_of_steps=1000000,
                 temperature=298.15,
                 lambda_steps=None,
                 COM_pull_goups=None,
                 harmonic_kappa=None,
                 pbc_atoms=None):

        if mdp_file[-4:] != '.mdp':
            mdp_file += '.mdp'

        self.mdp_file = mdp_file

        self.alchemical_molecule = alchemical_molecule

        self.timestep_ps = timestep_ps

        self.number_of_steps = number_of_steps

        self.temperature = temperature

        self.lambda_steps = lambda_steps

        self.COM_pull_goups = COM_pull_goups

        self.harmonic_kappa = harmonic_kappa

        self.pbc_atoms = pbc_atoms

        self._template = []

    def _hook(self):
        """a hook for subclasses
        """

    def _create_free_energy_strings(self,
                                    condition_lambda0,
                                    condition_lambda1,
                                    starting_lambda,
                                    couple_intramol='no',
                                    sc_alpha=0.0,
                                    sc_coul='no',
                                    sc_sigma=0.25,
                                    sc_power=1,
                                    nstdhdl=100,
                                    separate_dhdl_file='yes',
                                    free_energy='yes'):
        """Makes the free energy strings

        wrapper of `make_free_energy_lines` function
        """

        return make_free_energy_lines(
            condition_lambda0=condition_lambda0,
            condition_lambda1=condition_lambda1,
            alchemical_molecule=self.alchemical_molecule,
            lambda_step=self.lambda_steps,
            starting_lambda=starting_lambda,
            couple_intramol=couple_intramol,
            sc_alpha=sc_alpha,
            sc_coul=sc_coul,
            sc_sigma=sc_sigma,
            sc_power=sc_power,
            nstdhdl=nstdhdl,
            separate_dhdl_file=separate_dhdl_file,
            free_energy=free_energy)

    def _create_COMCOM_pulling_strings(self):
        """Private creates the strings for COM-COM pulling

        wrapper of `create_COMCOM_pulling_strings` function
        """

        if self.COM_pull_goups is None:
            return ['']
        elif len(self.COM_pull_goups) == 0:
            return ['']

        # pylint: disable=no-else-raise
        if self.harmonic_kappa is None:
            raise ValueError(
                'If you give some COM pull groups you shall give some harmonic constants'
            )
        elif len(self.harmonic_kappa) == 0:
            raise ValueError(
                'If you give some COM pull groups you shall give some harmonic constants'
            )

        return create_COMCOM_pulling_strings(
            COM_pull_goups=self.COM_pull_goups,
            harmonic_kappa=self.harmonic_kappa,
            pbc_atoms=self.pbc_atoms)

    def _get_template(self):
        """PRIVATE the template to write on MDP file

        Returns
        ----------
        tempalte : list
            list of strings with final newline
            already added ("\n")
        """

        for i in range(len(self._template)):

            if self._template[i] == '':
                self._template[i] += '\n'

            elif self._template[i][-1] != '\n':

                self._template[i] += '\n'

        return self._template

    def _write_template(self, template):
        """PRIVATE writes the template on a file
        """

        write_file.write_file(template, self.mdp_file)

    def execute(self):
        """writes the mdp

        it is the only pubblic method of the class
        and writes you the mdp file
        """

        #a hook for the subclasses
        self._hook()

        template = self._get_template()

        self._write_template(template=template)


class AnnihilateVdwMdpBoundState(MdpFile):
    """Creates mdp for VdW annihilation

    Remember to first annihilate the charges Q with
    `AnnihilateQMdpBoundState` class!
    See `MdpFile` class (the superclass) documentation
    """
    def _create_template(self):
        """PRIVATE creates the template

        it modifies self._template
        """

        self._template = [
            '; VARIOUS PREPROCESSING OPTIONS',
            '; Preprocessor information: use cpp syntax.',
            '; e.g.: -I/home/joe/doe -I/home/mary/roe',
            'include                  =',
            '; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)',
            'define                   =', '', '; RUN CONTROL PARAMETERS',
            'integrator               = md', '; Start time and timestep in ps',
            'tinit                    = 0',
            f'dt                       = {self.timestep_ps}',
            f'nsteps                   = {self.number_of_steps}',
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
            f'ref-t                    = {self.temperature}',
            '; pressure coupling',
            'pcoupl                   = Parrinello-Rahman',
            'pcoupltype               = Isotropic',
            'nstpcouple               = -1',
            '; Time constant (ps), compressibility (1/bar) and reference P (bar)',
            'tau-p                    = 0.5',
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
            'morse                    = no', ''
        ]

        self._template += self._create_free_energy_strings(
            condition_lambda0='none',
            condition_lambda1='vdw',
            starting_lambda=1,
            couple_intramol='no',
            sc_alpha=0.0,
            sc_coul='no',
            sc_sigma=0.25,
            sc_power=1,
            nstdhdl=100,
            separate_dhdl_file='yes',
            free_energy='yes')

        self._template += self._create_COMCOM_pulling_strings()

    def _hook(self):
        """defines self.lambda_steps if left None & creates self._template
        """

        if self.lambda_steps is None:
            self.lambda_steps = -1. / self.number_of_steps

        self._create_template()


class AnnihilateQMdpBoundState(MdpFile):
    """Creates mdp for Q annihilation

    annihilates charges
    See `MdpFile` class (the superclass) documentation
    """
    def _create_template(self):
        """PRIVATE creates the template

        it modifies self._template
        """

        self._template = [
            '; VARIOUS PREPROCESSING OPTIONS',
            '; Preprocessor information: use cpp syntax.',
            '; e.g.: -I/home/joe/doe -I/home/mary/roe',
            'include                  =',
            '; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)',
            'define                   =', '', '; RUN CONTROL PARAMETERS',
            'integrator               = md', '; Start time and timestep in ps',
            'tinit                    = 0',
            f'dt                       = {self.timestep_ps}',
            f'nsteps                   = {self.number_of_steps}',
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
            f'ref-t                    = {self.temperature}',
            '; pressure coupling',
            'pcoupl                   = Parrinello-Rahman',
            'pcoupltype               = Isotropic',
            'nstpcouple               = -1',
            '; Time constant (ps), compressibility (1/bar) and reference P (bar)',
            'tau-p                    = 0.5',
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
            'morse                    = no', ''
        ]

        self._template += self._create_free_energy_strings(
            condition_lambda0='vdw',
            condition_lambda1='vdw-q',
            starting_lambda=1,
            couple_intramol='no',
            sc_alpha=0.0,
            sc_coul='no',
            sc_sigma=0.25,
            sc_power=1,
            nstdhdl=100,
            separate_dhdl_file='yes',
            free_energy='yes')

        self._template += self._create_COMCOM_pulling_strings()

    def _hook(self):
        """defines self.lambda_steps if left None & creates self._template
        """

        if self.lambda_steps is None:
            self.lambda_steps = -1. / self.number_of_steps

        self._create_template()


class CreateVdwMdpUnboundState(MdpFile):
    """Creates mdp for VdW creation

    See `MdpFile` class (the superclass) documentation
    """
    def _create_template(self):
        """PRIVATE creates the template

        it modifies self._template
        """

        self._template = [
            '; VARIOUS PREPROCESSING OPTIONS',
            '; Preprocessor information: use cpp syntax.',
            '; e.g.: -I/home/joe/doe -I/home/mary/roe',
            'include                  =',
            '; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)',
            'define                   =', '', '; RUN CONTROL PARAMETERS',
            'integrator               = md', '; Start time and timestep in ps',
            'tinit                    = 0',
            f'dt                       = {self.timestep_ps}',
            f'nsteps                   = {self.number_of_steps}',
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
            f'ref-t                    = {self.temperature}',
            '; pressure coupling',
            'pcoupl                   = Parrinello-Rahman',
            'pcoupltype               = Isotropic',
            'nstpcouple               = -1',
            '; Time constant (ps), compressibility (1/bar) and reference P (bar)',
            'tau-p                    = 0.5',
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
            'morse                    = no', ''
        ]

        self._template += self._create_free_energy_strings(
            condition_lambda0='none',
            condition_lambda1='vdw',
            starting_lambda=0,
            couple_intramol='no',
            sc_alpha=0.5,
            sc_coul='no',
            sc_sigma=0.25,
            sc_power=1,
            nstdhdl=100,
            separate_dhdl_file='yes',
            free_energy='yes')

        self._template += self._create_COMCOM_pulling_strings()

    def _hook(self):
        """defines self.lambda_steps if left None & creates self._template
        """

        if self.lambda_steps is None:
            self.lambda_steps = 1. / self.number_of_steps

        self._create_template()


class CreateQMdpUnboundState(MdpFile):
    """Creates mdp for Q creation
    creates charges
    REMEMBER TO FIRST CREATE VdW with `CreateVdwMdpUnboundState`
    See `MdpFile` class (the superclass) documentation
    """
    def _create_template(self):
        """PRIVATE creates the template

        it modifies self._template
        """

        self._template = [
            '; VARIOUS PREPROCESSING OPTIONS',
            '; Preprocessor information: use cpp syntax.',
            '; e.g.: -I/home/joe/doe -I/home/mary/roe',
            'include                  =',
            '; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)',
            'define                   =', '', '; RUN CONTROL PARAMETERS',
            'integrator               = md', '; Start time and timestep in ps',
            'tinit                    = 0',
            f'dt                       = {self.timestep_ps}',
            f'nsteps                   = {self.number_of_steps}',
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
            f'ref-t                    = {self.temperature}',
            '; pressure coupling',
            'pcoupl                   = Parrinello-Rahman',
            'pcoupltype               = Isotropic',
            'nstpcouple               = -1',
            '; Time constant (ps), compressibility (1/bar) and reference P (bar)',
            'tau-p                    = 0.5',
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
            'morse                    = no', ''
        ]

        self._template += self._create_free_energy_strings(
            condition_lambda0='vdw',
            condition_lambda1='vdw-q',
            starting_lambda=0,
            couple_intramol='no',
            sc_alpha=0.5,
            sc_coul='no',
            sc_sigma=0.25,
            sc_power=1,
            nstdhdl=100,
            separate_dhdl_file='yes',
            free_energy='yes')

        self._template += self._create_COMCOM_pulling_strings()

    def _hook(self):
        """defines self.lambda_steps if left None & creates self._template
        """

        if self.lambda_steps is None:
            self.lambda_steps = 1. / self.number_of_steps

        self._create_template()
