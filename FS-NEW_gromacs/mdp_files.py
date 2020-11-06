# -*- coding: utf-8 -*-
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
    """
    def __init__(self,
                 mdp_file,
                 alchemical_molecule,
                 timestep_ps=0.001,
                 number_of_steps=1000000,
                 temperature=298.15,
                 lambda_steps=None,
                 COM_pull_goups=None):

        if mdp_file[-4:] != '.mdp':
            mdp_file += '.mdp'

        self.mdp_file = mdp_file

        self.alchemical_molecule = alchemical_molecule

        self.timestep_ps = timestep_ps

        self.number_of_steps = number_of_steps

        self.temperature = temperature

        self.lambda_steps = lambda_steps

        self.COM_pull_goups = COM_pull_goups

        self._template = []

    def _hook(self):
        """a hook for subclasses
        """

    def _create_COMCOM_pulling_strings(self):
        """Private creates the strings for COM-COM pulling

        Returns
        ----------
        string
            will return an empty string if `self.COM_pull_goups` is empty or None
        """

        if self.COM_pull_goups is None:
            return ''
        elif len(self.COM_pull_goups) == 0:
            return ''

        pull_ngroups = len(self.COM_pull_goups)

        pull_group_name = ''
        for i, group in enumerate(self.COM_pull_goups):
            pull_group_name += f'pull-group{i}-name        = {group}\n'

        COM_pulling_strings = [
            ';COM PULLING', 'pull                     = yes',
            'pull-print-com           = yes', 'pull-ncoords            = 2',
            'pull-nstxout            = 10',
            f'pull-ngroups            = {pull_ngroups}', f'{pull_group_name}',
            'pull-pbc-ref-prev-step-com  = yes',
            'pull-group1-pbcatom     = 2373',
            'pull-coord1-geometry    = distance',
            'pull-coord1-type        = umbrella',
            'pull-coord1-dim         = Y Y Y', 'pull-coord1-groups      = 1 2',
            'pull-coord1-start       = yes', 'pull-coord1-init       = 0.0',
            'pull-coord1-rate       = 0', 'pull-coord1-k          = 120',
            'pull-coord2-geometry    = distance',
            'pull-coord2-type        = umbrella',
            'pull-coord2-dim         = Y Y Y', 'pull-coord2-groups      = 1 3',
            'pull-coord2-start       = yes', 'pull-coord2-init       = 0.0',
            'pull-coord2-rate       = 0', 'pull-coord2-k          = 120'
        ]

        return '\n'.join(COM_pulling_strings)

    def _get_template(self):
        """PRIVATE the template to write on MDP file

        Returns
        ----------
        tempalte : list
            list of strings with final newline
            already added ("\n")
        """

        for i in range(len(self._template)):

            if self._template[i][-1] != '\n':

                self._template[i] += '\n'

        return self._template

    def _write_template(self, template):
        """PRIVATE writes the template on a file
        """

        with open(self.mdp_file, 'w') as output:
            for line in template:
                output.write(line)

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
            'morse                    = no', '', '; Free energy control stuff',
            'free-energy              = yes', 'init-lambda              = 1',
            f'delta-lambda             = {self.lambda_steps}',
            f'couple-moltype           = {self.alchemical_molecule}',
            'couple-lambda0           =none', 'couple-lambda1           =vdw',
            'couple-intramol          =no', 'sc-alpha                 = 0.0',
            'sc-coul                  = no', 'sc-sigma                 = 0.25',
            'sc-power                 = 1', ''
        ]

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
            'morse                    = no', '', '; Free energy control stuff',
            'free-energy              = yes', 'init-lambda              = 1',
            f'delta-lambda             = {self.lambda_steps}',
            f'couple-moltype           = {self.alchemical_molecule}',
            'couple-lambda0           =vdw', 'couple-lambda1           =vdw-q',
            'couple-intramol          =no', 'sc-alpha                 = 0.0',
            'sc-coul                  = no', 'sc-sigma                 = 0.25',
            'sc-power                 = 1', ''
        ]

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
            'morse                    = no', '', '; Free energy control stuff',
            'free-energy              = yes', 'init-lambda              = 0',
            f'delta-lambda             = {self.lambda_steps}',
            f'couple-moltype           = {self.alchemical_molecule}',
            'couple-lambda0           =none', 'couple-lambda1           =vdw',
            'couple-intramol          =no', 'sc-alpha                 = 0.5',
            'sc-coul                  = no', 'sc-sigma                 = 0.25',
            'sc-power                 = 1', ''
        ]

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
            'morse                    = no', '', '; Free energy control stuff',
            'free-energy              = yes', 'init-lambda              = 0',
            f'delta-lambda             = {self.lambda_steps}',
            f'couple-moltype           = {self.alchemical_molecule}',
            'couple-lambda0           =vdw', 'couple-lambda1           =vdw-q',
            'couple-intramol          =no', 'sc-alpha                 = 0.5',
            'sc-coul                  = no', 'sc-sigma                 = 0.25',
            'sc-power                 = 1', ''
        ]

    def _hook(self):
        """defines self.lambda_steps if left None & creates self._template
        """

        if self.lambda_steps is None:
            self.lambda_steps = 1. / self.number_of_steps

        self._create_template()
