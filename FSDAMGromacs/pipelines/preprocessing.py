# -*- coding: utf-8 -*-
# pylint: disable=too-many-locals
#############################################################
# Copyright (c) 2020-2021 Maurice Karrenbrock               #
#                                                           #
# This software is open-source and is distributed under the #
# BSD 3-Clause "New" or "Revised" License                   #
#############################################################
"""Pipeline classes for preprocessing
"""

import pathlib

import PythonAuxiliaryFunctions.path as _path
import PythonFSDAM.pipelines.superclasses as superclasses

import FSDAMGromacs.get_pbc_atom as get_pbc_atom
import FSDAMGromacs.mdp_files as mdp_files


class PreprocessGromacsFSDAM(superclasses.PreProcessingPipeline):
    """Preprocesses files fot non-equilibrium alchemical transformations

    It is the first step to make an absolute binding free energy, you must instantiate
    and run one instance for the bound and one for the unbound state to use Jarzynki

    it inherits from `PythonFSDAM.pipelines.superclasses.PreProcessingPipeline`

    Parameters
    -----------
    COM_pull_groups : list of str
        check the `mdp_files.MdpFile` for more info
    harmonic_kappa : nested list
        check the `mdp_files.MdpFile` for more info
    pbc_atoms : iterable
        check the `mdp_files.MdpFile` for more info
    creation : bool, default=True
        if True it will make a creation of the alchemical molecule
        if false an annihilation
    vdw_timestep_ps : float, optional
        timestep in ps for vdw alchemical transformation
        the default will be ok most of the times
    q_timestep_ps : float, optional
        timestep in ps for q alchemical transformation
        the default will be ok most of the times
    vdw_number_of_steps : int, optional
        number of md steps for vdw alchemical transformation
        the default will be ok most of the times
    q_number_of_steps : int, optional
        number of md steps for q alchemical transformation
        the default will be ok most of the times
    constrains : str, default all-bonds
        possible values none h-bonds all-bonds

    Attributes
    ------------
    DEFAULTS : dict
        the defaults for the alchemical tranformations
    """

    DEFAULTS = {
        'vdw': {
            'annihilation': {
                'timestep_ps': 0.0015,
                'number_of_steps': 800000
            },
            'creation': {
                'timestep_ps': 0.001,
                'number_of_steps': 250000
            }
        },
        'q': {
            'annihilation': {
                'timestep_ps': 0.0015,
                'number_of_steps': 500000
            },
            'creation': {
                'timestep_ps': 0.001,
                'number_of_steps': 250000
            }
        },
    }

    MDP_CLASSES = {
        'q': {
            'annihilation': mdp_files.AnnihilateQMdpBoundState,
            'creation': mdp_files.CreateQMdpUnboundState,
        },
        'vdw': {
            'annihilation': mdp_files.AnnihilateVdwMdpBoundState,
            'creation': mdp_files.CreateVdwMdpUnboundState,
        }
    }

    def __init__(self,
                 topology_files,
                 md_program_path,
                 alchemical_residue,
                 structure_files=None,
                 COM_pull_groups=None,
                 harmonic_kappa=None,
                 temperature=298.15,
                 pbc_atoms=None,
                 constraints=None,
                 creation=True,
                 vdw_timestep_ps=None,
                 q_timestep_ps=None,
                 vdw_number_of_steps=None,
                 q_number_of_steps=None,
                 constrains=None):

        super().__init__(topology_files=topology_files,
                         md_program_path=md_program_path,
                         alchemical_residue=alchemical_residue,
                         structure_files=structure_files,
                         temperature=temperature,
                         constrains=constrains)

        self.COM_pull_groups = COM_pull_groups

        self.harmonic_kappa = harmonic_kappa

        self.pbc_atoms = pbc_atoms

        self.constraints = constraints

        self.creation = creation

        self.vdw_timestep_ps = vdw_timestep_ps

        self.q_timestep_ps = q_timestep_ps

        self.vdw_number_of_steps = vdw_number_of_steps

        self.q_number_of_steps = q_number_of_steps

        self.md_program = 'gromacs'

    def _update_pbc_atoms(self):
        """Works only for standard and good behaving
        protein ligand systems
        """

        if self.pbc_atoms is None and self.COM_pull_groups is not None:

            if self.COM_pull_groups.count('Protein') != 1:

                raise Exception(
                    'Could not figure out the pbc atoms, please give them manually at instantiation'
                )

            tmp_pbc_atoms = []

            for item in self.COM_pull_groups:

                if item == 'Protein':

                    tmp_pbc_atoms.append(
                        get_pbc_atom.get_protein_pbc_atom(
                            self.structure_files[0]))

                else:

                    tmp_pbc_atoms.append(0)

            self.pbc_atoms = tuple(tmp_pbc_atoms)

    def _make_mdp_files(self, mdp_type, mdp_file):
        """Make a given mdp file

        Parameters
        ------------
        mdp_type : list of strings
            first string is 'vdw' or 'q'
            the second is 'annihilation'
            or 'creation'
        mdp_file : str or path
            the name of the file that will be created
        """

        self._update_pbc_atoms()

        ##############################################
        #set some defaults if no input has been given
        if mdp_type[0] == 'vdw':

            if self.vdw_timestep_ps is None:

                ts = self.DEFAULTS[mdp_type[0]][mdp_type[1]]['timestep_ps']

            else:

                ts = self.vdw_timestep_ps

            if self.vdw_number_of_steps is None:

                n = self.DEFAULTS[mdp_type[0]][mdp_type[1]]['number_of_steps']

            else:

                n = self.vdw_timestep_ps

        elif mdp_type[0] == 'q':

            if self.q_timestep_ps is None:

                ts = self.DEFAULTS[mdp_type[0]][mdp_type[1]]['timestep_ps']

            else:

                ts = self.q_timestep_ps

            if self.q_number_of_steps is None:

                n = self.DEFAULTS[mdp_type[0]][mdp_type[1]]['number_of_steps']

            else:

                n = self.q_timestep_ps
        ##############################################

        creator = self.MDP_CLASSES[mdp_type[0]][mdp_type[1]](
            mdp_file=mdp_file,
            timestep_ps=ts,
            number_of_steps=n,
            alchemical_molecule=self.alchemical_residue,
            temperature=self.temperature,
            COM_pull_goups=self.COM_pull_groups,
            harmonic_kappa=self.harmonic_kappa,
            pbc_atoms=self.pbc_atoms,
            constraints=self.constraints)

        creator.execute()

    def _make_tpr_lines(self, mdp_type, top_file, vdw_mdp, q_mdp):
        """PRIVATE

        makes the lines to create the TPR files

        RETURNS
        ----------
        make_vdw_tpr_commands, make_q_tpr_commands, vdw_tpr_files, q_tpr_files : list(str), list(str), list(str), list(str)
        """  # pylint: disable=line-too-long

        if self.creation:

            q_files = []

            for i in range(len(self.structure_files)):

                q_files.append(f'vdw_{mdp_type}_{i}.gro')

            vdw_files = [
                pathlib.Path(i).resolve().name for i in self.structure_files
            ]

        else:

            q_files = [
                pathlib.Path(i).resolve().name for i in self.structure_files
            ]

            vdw_files = []

            for i in range(len(self.structure_files)):

                vdw_files.append(f'q_{mdp_type}_{i}.gro')

        #create the create tpr commands in the right order and way

        #make vdw tpr files

        make_vdw_tpr_commands = []

        vdw_tpr_files = []

        for i, item in enumerate(vdw_files):

            command_string = \
                f'{self.md_program_path} grompp -f {vdw_mdp} ' + \
                    f'-c {item} -p {top_file} -o vdw_{mdp_type}_{i}.tpr -maxwarn 100'

            make_vdw_tpr_commands.append(command_string)

            vdw_tpr_files.append(f'vdw_{mdp_type}_{i}.tpr')

        #make q tpr files

        make_q_tpr_commands = []

        q_tpr_files = []

        for i, item in enumerate(q_files):

            command_string = \
                f'{self.md_program_path} grompp -f {q_mdp} ' + \
                    f'-c {item} -p {top_file} -o q_{mdp_type}_{i}.tpr -maxwarn 100'

            make_q_tpr_commands.append(command_string)

            q_tpr_files.append(f'q_{mdp_type}_{i}.tpr')

        return make_vdw_tpr_commands, make_q_tpr_commands, vdw_tpr_files, q_tpr_files

    @staticmethod
    def _get_suffix(path):
        """helper function to get file suffixes

        both for str and for pathlib.Path
        """

        if isinstance(path, str):

            return path[-4:]

        return path.suffix

    def _get_relative_path(self):
        """makes all paths relative to the CWD
        """

        for path_list in (self.topology_files, self.structure_files):

            for i, path in enumerate(path_list):

                if not isinstance(path, pathlib.Path):

                    path = pathlib.Path(path)

                path = _path.absolute_filepath(path)

                path_list[i] = path.relative_to(pathlib.Path().cwd())

    def execute(self):
        # pylint: disable=line-too-long
        """This method preprocesses the bound or unbound NE Alchemical transformations

        As gromacs forces you to divide the creation/annihilation in different runs this pipeline
        doesn't only preprocess but also run the MD simulation

        Returns
        ---------
        dict
            {
                'vdw_mdp' : path_to_vdw_mdp,

                'q_mdp' : path_to_q_mdp,

                'make_vdw_tpr' : list of commands to use to make the vdw tpr,

                'make_q_tpr' : list of commands to make the q tpr,

                'run_vdw' : list of commands to run the vdw tpr,

                'run_q' : list of commands to run the q tpr,

                'vdw_dhdl' : list of names of the vdw dhdl files that will be created (if you use the run commands unaltered),

                'q_dhdl' : list of names of the q dhdl files that will be created (if you use the run commands unaltered),
            }

        Notes
        -------
        If you are doing an annihilation first create and run the q tpr and later create and run the vdw one
        If you are doing a creation do it viceversa
        """

        #create the output dictionary
        output_dictionary = {}

        try:
            self._get_relative_path()
        except ValueError:
            pass

        if self.creation:

            mdp_type = 'creation'

        else:

            mdp_type = 'annihilation'

        #make vdw mdp files
        self._make_mdp_files(mdp_type=['vdw', mdp_type],
                             mdp_file=f'vdw_{mdp_type}.mdp')

        output_dictionary['vdw_mdp'] = f'vdw_{mdp_type}.mdp'

        #make q mdp files
        self._make_mdp_files(mdp_type=['q', mdp_type],
                             mdp_file=f'q_{mdp_type}.mdp')

        output_dictionary['q_mdp'] = f'q_{mdp_type}.mdp'

        #find the .top file (there are also .itp)
        for i in self.topology_files:

            if self._get_suffix(i) == '.top':

                top_file = i

                break

        else:
            raise ValueError('You gave no .top file as input')

        #create the create tpr commands in the right order and way
        (output_dictionary['make_vdw_tpr'], output_dictionary['make_q_tpr'],
         vdw_tpr_files, q_tpr_files) = (self._make_tpr_lines(
             mdp_type=mdp_type,
             top_file=top_file,
             vdw_mdp=output_dictionary['vdw_mdp'],
             q_mdp=output_dictionary['q_mdp']))

        #create the right run strings

        #create run_vdw

        run_vdw = []
        vdw_dhdl = []
        for i, item in enumerate(vdw_tpr_files):

            run_vdw.append(
                f'{self.md_program_path} mdrun -deffnm vdw_{mdp_type}_{i} -s {item}'
            )

            vdw_dhdl.append(f'vdw_{mdp_type}_{i}.xvg')

        output_dictionary['run_vdw'] = run_vdw

        output_dictionary['vdw_dhdl'] = vdw_dhdl

        #create run_q

        run_q = []
        q_dhdl = []
        for i, item in enumerate(q_tpr_files):

            run_q.append(
                f'{self.md_program_path} mdrun -deffnm q_{mdp_type}_{i} -s {item}'
            )

            q_dhdl.append(f'q_{mdp_type}_{i}.xvg')

        output_dictionary['run_q'] = run_q

        output_dictionary['q_dhdl'] = q_dhdl

        return output_dictionary
