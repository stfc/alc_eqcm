!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that generates:
! - atomistic models to match a target stoichiometry, which can be either
!   the stoichiometry of a pristine sample (defined in &block_species) or the
!   stoichiometries solutions from the charge and mass balance equations using
!   EQCM data. The user must provide a reference geometry as input structure
!   (file INPUT_GEOM/INPUT_STRUCTURE)
! - generates backup informtion of the atomistic models  
! - input files for the simulation of the atomistic models
! - scripts for job submission to HPC 
!
! In addition, this module can use the generated atomistic models to rebuild 
! input files for simulations and/or HPC scripts if the user decides any 
! parameter must be adjusted without the need to re-generate atomistic models    
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author       - i.scivetti  January-March 2021
! Contribution - i.scivetti  June-August   2021 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module atomistic_models 
  Use constants,          Only : avogadro,&
                                 Bohr_to_A,&
                                 chemsymbol, & 
                                 cm_to_Ang, &
                                 max_components,&
                                 NPTE,&
                                 twopi
                                 
  Use code_castep,        Only : print_castep_settings
  Use code_cp2k,          Only : print_cp2k_settings
  Use code_onetep,        Only : print_onetep_settings
  Use code_vasp,          Only : print_vasp_settings
                                 
  Use electrode,          Only : electrode_type
  Use eqcm,               Only : eqcm_type
  Use fileset,            Only : file_type,              &
                                 FILE_INTERCALATION_OX,  &
                                 FILE_INTERCALATION_RED, &
                                 FILE_INPUT_STRUCTURE,   &
                                 FILE_MODEL_SUMMARY,     &
                                 FILE_OUTPUT_STRUCTURE,  &
                                 FILE_SELECTED_SOLUTIONS,&
                                 FILE_SET_EQCM,          &
                                 FILE_HPC_SETTINGS, &
                                 FILE_KPOINTS, &
                                 FILE_RECORD_MODELS,&
                                 FILE_SET_SIMULATION, &
                                 FOLDER_ANALYSIS, &
                                 FOLDER_RESTART,&
                                 FOLDER_DFT, &
                                 FOLDER_INPUT_GEOM, &
                                 FOLDER_MODELS, &
                                 refresh_out_eqcm
  Use hpc,                Only : hpc_type, &
                                 summary_hpc_settings 
  Use input_types,        Only : in_integer, &
                                 in_integer_array, & 
                                 in_logic,   &
                                 in_string,  &
                                 in_param,   &
                                 in_scalar
  Use numprec,            Only : li, &
                                 wi, &
                                 wp
  Use process_data,       Only : capital_to_lower_case 
  Use redox,              Only : redox_type
  Use simulation_output,  Only : summary_simulation_settings,&
                                 warning_simulation_settings
  Use simulation_setup,   Only : simul_type
  Use stoichiometry,      Only : stoich_type
  Use unit_output,        Only : error_stop,&
                               info 

  Implicit None
  Private

  ! Maximum number of species
  Integer(Kind=wi), Parameter, Public  :: max_species=100
  ! Maximum number of atoms per species
  Integer(Kind=wi), Parameter  :: max_at_species= 500
  ! Maximum number of units per species 
  Integer(Kind=wi), Parameter  :: max_num_species_units=  100000
  ! Maximum number of models generated
  Integer(Kind=wi), Parameter  :: max_models=1000
  ! Maximum limit for the cell size to account for convergence adjustments
  Real(Kind=wp), Parameter, Public :: large_cell_limit = 50.0_wp
  ! Tolerance for length
  Real(Kind=wp), Parameter, Public :: length_tol = 1.0E-8
  ! Manimum amount of nearest neighbours for a single atom
  Integer(Kind=wi), Parameter, Public :: max_nn = 12 

  ! Limit for the maximum bond distance
  Real(Kind=wp), Parameter, Public   :: max_intra=2.30_wp  
  ! Limit for the minimum bond distance
  Real(Kind=wp), Parameter, Public   :: min_intra=0.70_wp  
  ! Limit for the minimum distance between atoms of different species 
  Real(Kind=wp), Parameter, Public   :: min_inter=1.70_wp  

  ! multiple of the total number of attempts to intercalate species 
  Integer(Kind=wi), Parameter  :: times_number_attempts=100


  ! Type to describe the atoms
  Type :: atom_type
     Real(Kind=wp)    :: r(3)
     Character(Len=8) :: tag
     Character(Len=2) :: element
     Integer(Kind=wi) :: atomic_number
     Logical          :: dynamics(3)
     Logical          :: in_species
     Logical          :: vanish
  End Type 
 
  ! Types related to format for VASP files 
  Type :: list_type
    Character(Len=8)  :: tag(max_components)
    Character(Len=2)  :: element(max_components) 
    Integer(Kind=wi)  :: N0(max_components)   
    Integer(Kind=wi)  :: num_elements
    Integer(Kind=wi)  :: net_elements
    Character(Len=32) :: coord_type 
  End Type 
 
  ! Type for the  &block_input_composition
  Type :: component_in_block
    Character(Len=8)  :: tag(max_components)
    Character(Len=2)  :: element(max_components)
    Integer(Kind=wi)  :: atomic_number(max_components)
    Integer(Kind=wi)  :: N0(max_components)
    Integer(Kind=wi)  :: numtot
  End Type 

  ! Type for the components of a given species 
  Type :: component_model
    Character(Len=8)  :: tag(max_components)
    Character(Len=2)  :: element(max_components)
    Integer(Kind=wi)  :: atomic_number(max_components)
    Integer(Kind=wi)  :: N(max_components)
    Integer(Kind=wi)  :: N0(max_components)
  End Type  
 
  ! Type for the definition of geometry of the species   
  Type :: define_species
     Character(Len=8) :: tag(max_at_species)
     Character(Len=2) :: element(max_at_species)
     Integer(Kind=wi) :: atomic_number(max_components)
     Real(Kind=wp)    :: r0(max_at_species,3)
     Logical          :: in_species(max_at_species)
  End Type
 
  ! Type for the species units
  Type :: species_units
     Integer(Kind=wi) :: list(max_at_species)
     Logical            :: vanish
  End Type 

  ! Type to describe species within the models
  Type :: species_model
    Character(Len=32)     :: tag
    Character(Len=32)     :: topology
    Integer(Kind=wi)      :: num
    Integer(Kind=wi)      :: num_show
    Integer(Kind=wi)      :: num_extra  
    Integer(Kind=wi)      :: D_num
    Integer(Kind=wi)      :: atoms_per_species 
    Integer(Kind=wi)      :: num_components
    Logical               :: change_content
    Real(Kind=wp)         :: s
    Type(component_model) :: component
    Type(define_species)  :: definition
    Type(species_units)   :: units(max_num_species_units)
  End Type 

  ! Type for the sample models 
  Type :: sample_type
    ! Path to subfolders
    Character(Len=256) :: path 
    ! Cell vectors 
    Real(Kind=wp)      :: cell(3,3)
    ! Inverse cell vectors
    Real(Kind=wp)      :: invcell(3,3)
    ! Length of cell vectors
    Real(Kind=wp)      :: cell_length(3)
    ! Scale factor vasp
    Real(Kind=wp)      :: scale_factor_vasp=1.0_wp
    ! Atoms in the model
    Type(atom_type), Allocatable :: atom(:)
    ! Maximum number of atom
    Integer(Kind=wi) :: max_atoms
    ! Total number of atoms  
    Integer(Kind=wi) :: num_atoms
    ! Total number of atoms (needed if inserting species) 
    Integer(Kind=wi) :: num_atoms_extra
    ! list elements
    Type(list_type)  :: list
    ! Minimum number of species units among all the species
    Integer(Kind=wi) :: min_species
    ! Minimum number of atomic  components among all the fixed species
    Integer(Kind=wi) :: min_components
    ! Stoichiometry of the atomic component with min_component
    Real(Kind=wp)      :: min_stoich
    ! Information for the species  
    Type(species_model), Allocatable :: species(:)
    ! Slab Area
    Real(Kind=wp)             :: slab_area
  End Type  

  ! Type for the modelling related variables 
  Type, Public :: model_type
    Private
    ! Multiple for the amount of output_atom
    Type(in_integer), Public :: multiple_output_atoms 
    ! Multiple for the amount of input_atoms
    Type(in_integer), Public :: multiple_input_atoms 
    ! Constrained dynamics
    Logical                  :: selective_dyn=.False.
    ! Change format
    Logical                  :: change_format
    ! Old format 
    Character(Len=32)        :: old_format 
    ! Flag to rotate the species
    Type(in_logic),  Public  :: rotate_species 
    ! Format of input model 
    Type(in_string), Public  :: input_model_format
    ! Format of output model 
    Type(in_string), Public  :: output_model_format
    ! Vector normal to the slab surface (only for electrodeposition) 
    Type(in_string), Public  :: normal_vector
    ! Flag to insert species
    Logical :: insert_species
    ! Flag to remove species
    Logical :: remove_species
    ! Total retition of the input model 
    Integer(Kind=wi) :: input_times
    ! Arrays for sampling the spatial region within the simulation cell 
    Integer(Kind=wi), Dimension(3)  :: scan
    Integer(Kind=wi), Dimension(3)  :: npoints 
    ! Distance cutoff (minium value for inter-species separation)
    Type(in_param), Public   :: distance_cutoff
    ! Repetition of the input model to build the output sample
    Type(in_integer_array), Public  :: repeat_input_model 
    ! Number of different species types
    Integer(Kind=wi) :: types_species 
    ! Type for those components defined in block_input_composition
    Type(component_in_block), Public  ::  block_component
    ! Block with vectors for the input structure model 
    Type(in_string), Public  :: block_input_cell 
    ! Block for tags and amount of atoms that compose the input model
    Type(in_string), Public  :: block_input_composition 
    ! Input model 
    Type(sample_type), Public  :: input 
    ! Generated atomic sample
    Type(sample_type), Public  :: sample
    ! Accepted stoichiometry error for models
    Type(in_scalar), Public   :: stoichiometry_error 
    ! Target number of models to sample the stoichiometric space when the solutions are multiple
    Type(in_integer), Public  :: targeted_num_models
    ! Number of generated models to sample the stoichiometric space when the solutions are multiple
    Integer(Kind=wi) :: num_models_sampling
    ! Arrays with the stoichiometry solutions to build models
    Real(Kind=wp)             ::  stoich_target(max_species, max_models)
    ! Reference stoichiometry
    Real(Kind=wp)             :: ref_stoich
    ! Decide if the structure is shifted or not
    Type(in_logic), Public    :: shift_structure
    ! Level where to start deposition from
    Type(in_param), Public    :: deposition_level
    ! Discretization of the spatial coordinates
    Type(in_param), Public   :: delta_space
    ! scale_cell
    Type(in_scalar), Public  :: scale_cell
  Contains
     Private
     Procedure         :: atomic_arrays_input   => allocate_input_atomic_arrays
     Procedure         :: species_arrays        => allocate_species_arrays
     Procedure         :: atomic_arrays_model   => allocate_atomic_arrays_model
     Procedure, Public :: init_input_variables  => allocate_input_variables
     Final             :: cleanup
  End Type model_type

  Public :: check_atomic_settings, build_atomistic_models, generate_hpc_simulation_files 

Contains

  Subroutine allocate_input_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential variables to build atomistic models
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type), Intent(InOut)  :: T

    Integer(Kind=wi)     :: fail(1)
    Character(Len=256)   :: message

    Allocate(T%repeat_input_model%value(3), Stat=fail(1))

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for input variables needed for&
                                & atomistic modelling (subroutine allocate_input_variables)'
      Call error_stop(message)
    End If

  End Subroutine allocate_input_variables

  Subroutine allocate_species_arrays(T, stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate species array for characterization of 
    ! the participating species 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type), Intent(InOut)  :: T
    Type(stoich_type), Intent(In   )  :: stoich_data

    Integer(Kind=wi)    :: fail(2), i, j, num_species
    Character(Len=256)  :: message

    fail=0
    num_species= stoich_data%num_species%value

    Allocate(T%sample%species(num_species),   Stat=fail(1))
    Allocate(T%input%species(num_species),     Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for species arrays&
                               & (subroutine allocate_species_arrays)'
      Call error_stop(message)
    End If

    ! Copy variables from the stoich_data%species
    Do i = 1, num_species
      T%input%species(i)%tag = stoich_data%species(i)%tag
      T%input%species(i)%topology = stoich_data%species(i)%topology
      T%input%species(i)%atoms_per_species = stoich_data%species(i)%atoms_per_species
      T%input%species(i)%num_components = stoich_data%species(i)%num_components
      T%input%species(i)%component%tag = stoich_data%species(i)%component%tag
      T%input%species(i)%component%element = stoich_data%species(i)%component%element
      T%input%species(i)%component%atomic_number = stoich_data%species(i)%component%atomic_number
      T%input%species(i)%component%N0=stoich_data%species(i)%component%N0
      If (Trim(stoich_data%species(i)%vartype0)=='fixed') Then
        T%input%species(i)%change_content=.False.
      Else
        T%input%species(i)%change_content=.True.
      End If 
    End Do
    
    ! Copy from input to sample
    Do i = 1, num_species
      T%sample%species(i)%tag=T%input%species(i)%tag
      T%sample%species(i)%topology=T%input%species(i)%topology
      T%sample%species(i)%atoms_per_species=T%input%species(i)%atoms_per_species
      T%sample%species(i)%component=T%input%species(i)%component
      T%sample%species(i)%num_components=T%input%species(i)%num_components
      T%sample%species(i)%definition=T%input%species(i)%definition
      T%sample%species(i)%change_content=T%input%species(i)%change_content
    End Do

    Do i = 1, num_species
      Do j = 1, max_num_species_units
        T%sample%species(i)%units(j)%vanish=.False.
        T%input%species(i)%units(j)%vanish=.False.
      End Do
    End Do

  End Subroutine allocate_species_arrays

  Subroutine allocate_atomic_arrays_model(T,num_atoms)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate atomic arrays for the sample model
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),   Intent(InOut)  :: T
    Integer(Kind=wi),    Intent(In   )  :: num_atoms

    Integer(Kind=wi)     :: fail(1)
    Character(Len=256)   :: message

    fail=0
    T%sample%max_atoms=T%multiple_output_atoms%value*T%input_times*num_atoms  

    Allocate(T%sample%atom(T%sample%max_atoms),        Stat=fail(1))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for the atomic arrays&
                               & of output structure model (subroutine allocate_atomic_arrays_model)'
      Call error_stop(message)
    End If

    T%sample%atom(:)%in_species=.False.
    T%sample%atom(:)%vanish=.False.

  End Subroutine allocate_atomic_arrays_model

  Subroutine allocate_input_atomic_arrays(T,num_atoms)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate atomic arrays for the input model
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),  Intent(InOut)  :: T
    Integer(Kind=wi),   Intent(In   )  :: num_atoms

    Integer(Kind=wi)    :: fail(1)
    Character(Len=256)  :: message

    T%input%max_atoms=T%multiple_input_atoms%value*num_atoms  
    fail=0

    Allocate(T%input%atom(T%input%max_atoms),   Stat=fail(1))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for the atomic arrays of input model&
                                & (subroutine allocate_input_atomic_arrays). It seems the input model is&
                                & too large. Reduce the value of "multiple_input_atoms" (set to 10000 by default)'
      Call error_stop(message)
    End If

    T%input%atom(:)%in_species=.False.
    T%input%atom(:)%vanish=.False.

  End Subroutine allocate_input_atomic_arrays

  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(model_type) :: T

    If (Allocated(T%input%atom)) Then
      Deallocate(T%input%atom)
    End If 

    If (Allocated(T%sample%atom)) Then
      Deallocate(T%sample%atom)
    End If 

    If (Allocated(T%input%species)) Then
      Deallocate(T%input%species)
    End If 

    If (Allocated(T%sample%species)) Then
      Deallocate(T%sample%species)
    End If 

    If (Allocated(T%repeat_input_model%value)) Then
      Deallocate(T%repeat_input_model%value)
    End If 

  End Subroutine cleanup 

  Subroutine build_atomistic_models(files, electrode_data, eqcm_data, redox_data,&
                                  & stoich_data, model_data, simulation_data, hpc_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to build atomistic models and input files for simulation
    !
    ! author    - i.scivetti Jan 2021
    ! contribution - i.scivetti  June-August   2021 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(electrode_type),Intent(In   ) :: electrode_data
    Type(eqcm_type),     Intent(In   ) :: eqcm_data
    Type(redox_type),    Intent(In   ) :: redox_data
    Type(stoich_type),   Intent(In   ) :: stoich_data
    Type(model_type),    Intent(InOut) :: model_data
    Type(simul_type),    Intent(InOut) :: simulation_data
    Type(hpc_type),      Intent(In   ) :: hpc_data

    Logical :: activate_random, alloc_sample_arrays, activate_loop
    Character(Len= 16) :: name_leg, numchar
    Character(Len=256) :: messages(5)
    Character(Len=256) :: message
    Character(Len=256) :: type_sim, process
    Real(Kind=wp), Dimension(3)  :: v1, v2
    Real(Kind=wp), Allocatable  :: target_stoich(:)
    Integer(Kind=wi)   :: i, ic, jc, kc, record_unit
    Integer(Kind=wi)   :: cycles, legs
    Integer(Kind=wi)   :: fail(1)
    Logical            :: fortho
    !Time related variables
    Real(Kind=wp)      :: t_ini, t_final
  
    activate_random=.True.
    activate_loop=.True.
    alloc_sample_arrays=.True.
    type_sim=eqcm_data%analysis%type
    process= eqcm_data%process%type

    model_data%input_times=1
    Do i = 1,3 
      model_data%input_times=model_data%input_times*model_data%repeat_input_model%value(i)
    End Do

    Call info(' ', 1)
    If (Trim(type_sim)=='model_pristine_sample' .Or. Trim(type_sim)=='model_disordered_system') Then
      Write (messages(1),'(1x,a)') '====================='
      Write (messages(2),'(1x,a)') 'Build atomistic model'
      Write (messages(3),'(1x,a)') '====================='
      Call info(messages, 3)
      cycles=1
      legs=1
    Else If (Trim(type_sim)=='model_cycled_sample') Then
      Write (messages(1),'(1x,a)') '======================'
      Write (messages(2),'(1x,a)') 'Build atomistic models'
      Write (messages(3),'(1x,a)') '======================'
      Call info(messages, 3)
      ! Define limits
      cycles=redox_data%limit_cycles
      legs=redox_data%legs
    End If

    ! Refresh out_eqcm
    Call refresh_out_eqcm(files)
     
    Call print_atomistic_settings(model_data)
    Call read_input_model(files, stoich_data, model_data)
    If (model_data%shift_structure%stat) Then
      Call shift_structure(model_data, process)
    End If
    If (simulation_data%solvation%info%stat) Then
      Call check_orthorhombic_cell(model_data%input%cell, fortho) 
      If (.Not.fortho) Then
        Write (message,'(1x,1a)') '***ERROR: computation with implicit solvent is only possible for&
                                 & orthorhombic cells.'
        Call error_stop(message)
      End If        
    End If
    Call check_cell_consistency(model_data)
    Call model_data%species_arrays(stoich_data)  
    Call compute_stoichiometry_input(stoich_data, model_data)

    model_data%types_species=stoich_data%num_species%value
    Allocate(target_stoich(model_data%types_species), Stat=fail(1))

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for "target_stoich" array in subroutine&
                                & "build_atomistic_models". Please check'
      Call error_stop(message)
    End If

    Call identify_species_input(stoich_data,model_data)

    ! Scale atomic coordinates only if scale_cell =/ 1.0 
    If (Abs(model_data%scale_cell%value-1.0_wp)>epsilon(model_data%scale_cell%value)) Then 
      Call scale_model(model_data)
    End If

    If (Trim(type_sim)=='model_pristine_sample' .Or. Trim(type_sim)=='model_disordered_system') Then
      Call stoichiometry_input_model(files, stoich_data, model_data)
    End If

    ! Define simulation cell
    Call define_cell(model_data, simulation_data)

    Call info(' ', 1)
    If (Trim(type_sim)=='model_pristine_sample' .Or. Trim(type_sim)=='model_disordered_system') Then
      Call info(' Atomistic details for the generated model of the pristine sample', 1)
    Else
      Call info(' Atomistic details for the generated models of the cycled samples', 1)
    End If
      Call info(' ----------------------------------------------------------------', 1)
    ! Print a message with the size of the generated model
    If (model_data%repeat_input_model%fread) Then
      If (model_data%repeat_input_model%value(1)==1 .And. &
          model_data%repeat_input_model%value(2)==1 .And. &
          model_data%repeat_input_model%value(3)==1) Then
          Write (message,'(1x,a,3(i3,a))') 'The size for the generated model(s) is equal to the input model.'
      Else
          Write (message,'(1x,a,3(i3,a))') 'The size for the generated model(s) is (', &
                                      & model_data%repeat_input_model%value(1), ',',&
                                      & model_data%repeat_input_model%value(2), ',',&
                                      & model_data%repeat_input_model%value(3),     &
                                      & ') times the input model along its cell vectors.'
      End If 
    Else
      Write (message,'(1x,a)') 'By default, the generated model(s) have the same as the input model.&
                               & The user can increase the size with directive&
                               & "repeat_input_model".'
    End If 
    Call info(message,1)
 
    Write (message,'(1x,a)') 'Model cell vectors:'
    Call info(message,1)
    Write (message,'(1x,a,3f20.12,a)') 'c1 = (', (model_data%sample%cell(1,i), i=1,3), ' )'
    Call info(message,1)
    Write (message,'(1x,a,3f20.12,a)') 'c2 = (', (model_data%sample%cell(2,i), i=1,3), ' )'
    Call info(message,1)
    Write (message,'(1x,a,3f20.12,a)') 'c3 = (', (model_data%sample%cell(3,i), i=1,3), ' )'
    Call info(message,1)
    Call info(' ',1)

    ! Compute the surface area in case of electrodeposition 
    If (Trim(process)=='electrodeposition') Then
      If (Trim(model_data%normal_vector%type)=='c3') Then
        v1(:)=model_data%input%cell(1,:)
        v2(:)=model_data%input%cell(2,:)
      ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
        v1(:)=model_data%input%cell(1,:)
        v2(:)=model_data%input%cell(3,:)
      ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
        v1(:)=model_data%input%cell(2,:)
        v2(:)=model_data%input%cell(3,:)
      End If
      Call compute_area_slab(v1,v2,model_data%input%slab_area)
    End If

    ! Opening file for record of relevant modelling settings
    Open(Newunit=files(FILE_RECORD_MODELS)%unit_no, File=files(FILE_RECORD_MODELS)%filename,Status='Replace')
    record_unit=files(FILE_RECORD_MODELS)%unit_no
    Write (record_unit, '(2(3x,a))') Trim(type_sim), Trim(model_data%output_model_format%type) 
 
    ic=1
    Do While (ic <= cycles .And. activate_loop)
      jc=1
      Do While (jc <= legs .And. activate_loop)

        If (Trim(type_sim)=='model_pristine_sample' .Or. Trim(type_sim)=='model_disordered_system') Then
          model_data%num_models_sampling=1
        Else If (Trim(type_sim)=='model_cycled_sample') Then        
          name_leg=redox_data%label_leg(jc,ic)
          Call stoichiometry_solutions(ic, jc , name_leg, process, files, electrode_data, stoich_data, model_data)
          If (Trim(process)=='electrodeposition') Then
            model_data%num_models_sampling=1
          End If
        End If

        If (model_data%num_models_sampling/=1) Then
          activate_loop=.False. 
          Write (messages(1), '(1x,a,i3,2a)') 'Cycle', ic, '-', Trim(name_leg)
          Write (messages(2), '(1x,a)')&
                   & '====================================================================================================='
          Write (messages(3), '(1x,a,i3,a)')  '***IMPORTANT: multiple stoichiometric solutions computed. A total of ',&
                                             & model_data%num_models_sampling, ' models for this redox part'
          Write (messages(4), '(1x,a)')       'will be generated and the execution will stop afterwards.&
                                             & The set of selected stoichiometries for'
          Write (messages(5), '(1x,3a)')      'these models are printed to file ',&
                                             & Trim(FOLDER_ANALYSIS)//'/'//Trim(files(FILE_SELECTED_SOLUTIONS)%filename)
          Call info(messages, 5)             
          If (model_data%targeted_num_models%value /= model_data%num_models_sampling) Then
            Write (messages(1), '(1x,a)')     '***WARNING: Note that the number of generated models differs from the value&
                                             & specified in directive'
             Write (messages(2), '(1x,a)')    '"targeted_num_models". This is due to the implemented algorithm (see manual)'
            Call info(messages,2)
          End If 
          Call info(' =====================================================================================================', 1)
          Call info(' ', 1)       
        End If
 
        Do kc=1, model_data%num_models_sampling 

          ! Record initial time
          Call cpu_time(t_ini) 

          If (Trim(process)=='intercalation' .Or. (ic==1 .And. jc==1)) Then   
            !Assign the initial value for the new number of atoms to be equal to the number of atoms in then structure 
             model_data%input%num_atoms_extra=model_data%input%num_atoms
          Else
             Call identify_species_input(stoich_data,model_data)
          End If 
          !Assign the target stoichiometry 
          If (Trim(type_sim)=='model_pristine_sample' .Or. Trim(type_sim)=='model_disordered_system') Then
             target_stoich(:) = stoich_data%species(:)%s0_pristine
          Else If (Trim(type_sim)=='model_cycled_sample') Then
             Do i=1, model_data%types_species 
               target_stoich(i) = model_data%stoich_target(i,kc) 
             End Do
             If (model_data%num_models_sampling==1) Then
               Write (messages(1), '(a,i3,2a)') ' Cycle', ic, '-', Trim(name_leg)
             Else
               Write (numchar, '(i3)') kc
               Write (messages(1), '(a,i3,5a)') ' Cycle', ic, '-', Trim(name_leg), ' (model ', Trim(Adjustl(numchar)), ')'
             End If
             Call info(messages,1)
           End If

           Call match_target_stoichiometry(model_data, target_stoich, 'input')

          ! Set initialization for random number: This is done only once when activate_random=.True.
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
          If ((model_data%remove_species .Or. model_data%insert_species) .And. activate_random) Then
            Call init_random_seed()         
            activate_random=.False.
          End If
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

          ! Remove species in the input model to match target stoichiometry
          If (model_data%remove_species) Then
            Call remove_species(model_data%input, model_data%types_species, process, 'input')
          End If
  
          ! Insert species to match target stoichiometry. This is done after having removed species to make some space
          If (model_data%insert_species) Then
            Call insert_species(process, stoich_data, model_data)
          End If

          ! Calculate the amount of atoms of the modelled sample
          model_data%sample%num_atoms=model_data%input_times*model_data%input%num_atoms_extra

          ! Allocate atomic arrays for the output model, only once
          If (alloc_sample_arrays) Then
            Call model_data%atomic_arrays_model(model_data%sample%num_atoms)  
            alloc_sample_arrays=.False.
          End If

          Call define_repeated_model(model_data)
          Call match_target_stoichiometry(model_data, target_stoich, 'output')

          ! Remove species in the input model to match target stoichiometry
          If (model_data%remove_species) Then
            Call remove_species(model_data%sample, model_data%types_species, process, 'sample')
          End If

          ! Print modelled stoichiometry
          Call stoichiometry_sample_model(files, target_stoich, model_data)      

          ! Create list
          Call create_list_net_elements(model_data%sample, model_data%types_species)

          ! Print output file
          Call print_output_files(type_sim, files, model_data, simulation_data, hpc_data, ic, kc, name_leg)

          ! Record relevant settings to file RECORD_MODELS
          Call record_models(record_unit, model_data)  

          ! Reset values 
          If (Trim(process)=='electrodeposition') Then 
            model_data%input%num_atoms=model_data%input%num_atoms_extra
            model_data%input%species(:)%num=model_data%input%species(:)%num_extra
          End If
          Call reset_to_input(model_data)

          ! Record initial time
          Call cpu_time(t_final) 

          ! Refresh out_eqcm
          Call refresh_out_eqcm(files)
        End Do
        jc=jc+1
      End Do
      ic=ic+1
    End Do

    ! If the user has requested the generation of simulation files, the following will provide summary of the 
    ! settings as well as aspects to be taken into consideration
    ! Print summary of HPC settings requested by the user
    If (simulation_data%generate) Then
      Call info(' ', 1)
      Call info(' =======================================================', 1)
      Call info(' Summary of the generated input settings for simulations', 1)
      Call info(' =======================================================', 1)
      Call summary_simulation_settings(simulation_data)
    Else
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'INFO: No input files for atomistic level simulations have been generated.'
      Write (messages(2),'(1x,a)') 'To generate input files, the user must define&
                                 & "&block_simulation_settings" in the SET_EQCM file.'
      Call info(messages,2)
    End If
    
    If (hpc_data%generate) Then
      Write (messages(1),'(1x,a)')  'In addition, the user has requested to build HPC script files.'
      Write (messages(2),'(1x,3a)') 'Each sub-folder contains the generated "', Trim(hpc_data%script_name),&
                                   &'" file for job submission.' 
      Call info(messages,2)
      Call summary_hpc_settings(hpc_data)
    Else
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'INFO: No scripts for HPC job submission has been generated.'
      Write (messages(2),'(1x,a)') 'To generate HPC scripts, the user must define "&block_hpc_settings"&
                                 & in the SET_EQCM file.'
      Call info(messages,2)
    End If

    If ((.Not.simulation_data%generate) .Or. (.Not. hpc_data%generate)) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)')  'INFO: Input files for atomistic level simulations and/or HPC scripts can be&
                                  & generated/updated without' 
      Write (messages(2),'(1x,a)')  'the need to re-build atomistic models using the option "hpc_simulation_files"&
                                  & for directive "analysis".'
      Call info(messages,2)
    End If

    If (simulation_data%generate) Then
      ! Print warnings
      Call warning_simulation_settings(simulation_data)
    End If
    
    ! Close RECORD_MODELS
    Close(record_unit)

    ! Deallocate 
    Deallocate(target_stoich)

  End Subroutine build_atomistic_models

  Subroutine scale_model(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to scale the cell of the atomistic model 
    !
    ! author    - i.scivetti June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),  Intent(InOut) :: model_data

    Integer(Kind=wi) :: i
    Real(Kind=wp)    :: v_direct(3), v_cart(3)

    Do i=1, model_data%input%num_atoms
      v_direct= MatMul(model_data%input%atom(i)%r, model_data%input%invcell)
      model_data%input%atom(i)%r=v_direct
    End Do

    model_data%input%cell=model_data%scale_cell%value * model_data%input%cell 
    Call about_cell(model_data%input%cell,model_data%input%invcell,model_data%input%cell_length)

    Do i=1, model_data%input%num_atoms
      v_cart= MatMul(model_data%input%atom(i)%r, model_data%input%cell)
      model_data%input%atom(i)%r=v_cart
    End Do

  End Subroutine scale_model 


  Subroutine define_cell(model_data, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to define the cell of the atomistic model 
    !
    ! author    - i.scivetti June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),  Intent(InOut) :: simulation_data
    Type(model_type),  Intent(InOut) :: model_data

    Integer(Kind=wi) :: i

    ! Set size for the output sample 
    Do i = 1, 3 
       model_data%sample%cell(i,:)=model_data%input%cell(i,:)*model_data%repeat_input_model%value(i) 
    End Do

    !Compute cell properties and inverse
    Call about_cell(model_data%sample%cell,model_data%sample%invcell,model_data%sample%cell_length)    

    ! In case the user wants to generate input files for simulations
    If (simulation_data%generate) Then
     !Call large_cell(model_data, simulation_data) 
      ! copy cell info to simulation_data%cell 
      simulation_data%cell=model_data%sample%cell
      simulation_data%cell_length=model_data%sample%cell_length
      ! If the size of the cell is large...
      simulation_data%large_cell=.False.
      Do i=1,3
        If (model_data%sample%cell_length(i)>large_cell_limit) Then
          simulation_data%large_cell=.True.
        End If
      End Do
    End If

  End Subroutine define_cell

  Subroutine record_models(record_unit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to keep relevant settings from the generation 
    ! of atomistic models 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(InOut) :: record_unit
    Type(model_type),  Intent(InOut) :: model_data

    Integer(Kind=wi) :: i

    Write (record_unit,*) Trim(model_data%sample%path)
    Write (record_unit,*) model_data%sample%list%net_elements
    Write (record_unit,'(*(6x,a4))') (Trim(model_data%sample%list%tag(i)), i=1, model_data%sample%list%net_elements)
    Write (record_unit,'(*(6x,a4))') (model_data%sample%list%element(i), i=1, model_data%sample%list%net_elements)
    Write (record_unit,'(*(2x,i8))')    (model_data%sample%list%N0(i),  i=1, model_data%sample%list%net_elements)

  End Subroutine record_models
   
  Subroutine check_atomic_settings(files, eqcm_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives for building atomistic models
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(eqcm_type),     Intent(In   ) :: eqcm_data
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256)  :: messages(8)
    Character(Len=64 )  :: error_set_eqcm
    Integer(Kind=wi)    :: i

    error_set_eqcm = '***ERROR in file '//Trim(files(FILE_SET_EQCM)%filename)//' -'

    ! Check "input_model_format"
    If (model_data%input_model_format%fread) Then
      If (model_data%input_model_format%fail) Then
        Write (messages(1),'(2a)')  Trim(error_set_eqcm), ' Wrong specification for directive "input_model_format"'
        Call info(messages,1)
        Call error_stop(' ') 
      Else
        If (Trim(model_data%input_model_format%type) /= 'vasp'    .And. &
          Trim(model_data%input_model_format%type) /= 'onetep'   .And. &
          Trim(model_data%input_model_format%type) /= 'castep'   .And. &
          Trim(model_data%input_model_format%type) /= 'cif'   .And. &
          Trim(model_data%input_model_format%type) /= 'xyz') Then
          Write (messages(1),'(2a)') Trim(error_set_eqcm), ' Specification for directive "input_model_format" should either be:'
          Write (messages(2),'(1x,a)') '- vasp'
          Write (messages(3),'(1x,a)') '- onetep'
          Write (messages(4),'(1x,a)') '- castep'
          Write (messages(5),'(1x,a)') '- cif'
          Write (messages(6),'(1x,a)') '- xyz'
          Call info(messages, 6)
          Call error_stop(' ') 
        End If
      End If
    Else
      Write (messages(1),'(2a)')  Trim(error_set_eqcm), 'Atomistic analysis requires specification of directive&
                                                      & "input_model_format".'
      Call info(messages,1)
      Call error_stop(' ') 
    End If

    ! Check "output_model_format"
    If (model_data%output_model_format%fread) Then
      If (model_data%output_model_format%fail) Then
        Write (messages(1),'(2a)')  Trim(error_set_eqcm), ' Wrong specification for directive "output_model_format"'
        Call info(messages,1)
        Call error_stop(' ') 
      Else
        If (Trim(model_data%output_model_format%type) /= 'vasp' .And. &
          Trim(model_data%output_model_format%type) /= 'onetep'  .And. &
          Trim(model_data%output_model_format%type) /= 'castep'   .And. &
          Trim(model_data%output_model_format%type) /= 'cp2k'   .And. &
          Trim(model_data%output_model_format%type) /= 'siesta'   .And. &
          Trim(model_data%output_model_format%type) /= 'cif'   .And. &
          Trim(model_data%output_model_format%type) /= 'xyz') Then
          Write (messages(1),'(2a)') Trim(error_set_eqcm), ' Specification for directive "output_model_format" should either be:'
          Write (messages(2),'(1x,a)') '- vasp'
          Write (messages(3),'(1x,a)') '- onetep'
          Write (messages(4),'(1x,a)') '- castep'
          Write (messages(5),'(1x,a)') '- xyz'
          Write (messages(6),'(1x,a)') '- cif'
          Write (messages(7),'(1x,a)') '- cp2k'
          Write (messages(8),'(1x,a)') '- siesta`'
          Call info(messages,8)
          Call error_stop(' ') 
        End If
      End If
    Else
      Write (messages(1),'(2a)')  Trim(error_set_eqcm), ' Requested analysis needs the specification of&
                                                      & directive "output_model_format"'
      Call info(messages,1)
      Call error_stop(' ') 
    End If

    ! Need to read &block_input_cell only if input_model_format is "xyz"
    If (Trim(model_data%input_model_format%type) == 'xyz' ) Then
      If (.Not. model_data%block_input_cell%fread) Then
        Write (messages(1),'(4a)')  Trim(error_set_eqcm), ' Format ', Trim(model_data%input_model_format%type),&
                                                      & ' for the input model&
                                                      & requires the specification of the cell vectors via &block_input_cell'
        Call info(messages,1)
        Call error_stop(' ') 
      End If
    Else
      If (model_data%block_input_cell%fread) Then 
        Write (messages(1),'(4a)') ' Specification for cell vectors in &block_input_cell will be neglected. Cell vectors must&
                                 & be specified in file ', Trim(files(FILE_SET_EQCM)%filename), ' according to the format',&
                                 & Trim(model_data%input_model_format%type)  
        Call info(messages,1)
      End If
    End If

    ! Check if directive "perpendicular_cell_vector" has been defined correctly
    If (Trim(eqcm_data%process%type) == 'electrodeposition') Then
      If (model_data%normal_vector%fread) Then
        If (model_data%normal_vector%fail) Then
          Write (messages(1),'(2a)')  Trim(error_set_eqcm), ' Wrong specification for directive "normal_vector"'
          Call info(messages,1)
          Call error_stop(' ')
        Else
          If (Trim(model_data%normal_vector%type) /= 'c1' .And. &
             Trim(model_data%normal_vector%type) /= 'c2' .And. &
             Trim(model_data%normal_vector%type) /= 'c3') Then
            Write (messages(1),'(2a)') Trim(error_set_eqcm), ' Specification for directive "normal_vector" should either be &
                                     &"c1", "c2" or "c3", which refer to the three cell vectors'
            Call info(messages,1)
            Call error_stop(' ')
          End If
        End If
      Else
        Write (messages(1),'(2a)')  Trim(error_set_eqcm), 'Building models for electrodeposition requires specification of&
                                  & directive "normal_vector"'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      If (model_data%normal_vector%fread) Then
         Write (messages(1),'(2a)') Trim(error_set_eqcm), ' Specification of directive "normal_vector" is only meaningful&
                                  & to electrodeposition processes. Otherwise, it must be removed.'
         Call info(messages,1)
         Call error_stop(' ')
      End If
    End If

    If (.Not. model_data%block_input_composition%fread) Then
      Write (messages(1),'(2a)')  Trim(error_set_eqcm),&
                              & 'Building atomistic models requires the specification of &block_input_composition'
      Call info(messages,1)
      Call error_stop(' ')
    End If

    If (Trim(model_data%input_model_format%type) == 'cif') Then
      Call execute_command_line('[ ! -e cif-input-geom.py ]', Exitstat=i)
      If (i==0) Then
        Call info(' ', 1)
        Write (messages(1),'(4(1x,a))') Trim(error_set_eqcm), 'File cif-input-geom.py not found.&
                                 & Please copy the file from the "scripts" folder located in the ALC_EQCM root directory&
                                 & and execute again'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    End If

    If (Trim(model_data%output_model_format%type) == 'cif') Then
      Call execute_command_line('[ ! -e cif-output-geom.py ]', Exitstat=i)
      If (i==0) Then
        Call info(' ', 1)
        Write (messages(1),'(4(1x,a))') Trim(error_set_eqcm), 'File cif-output-geom.py not found. &
                                 & Please copy the file from the "scripts" folder located in the ALC_EQCM root directory&
                                 & and execute again'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    End If

  End Subroutine check_atomic_settings


  Subroutine print_atomistic_settings(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print relevant settings used to building atomistic models
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256)  :: message

    ! Print stoichiometry error
    If (model_data%stoichiometry_error%fread) Then
      Write (message, '(1x,a,f8.3)') 'Error tolerance for the stoichiometry of the generated models:', &
                                   & model_data%stoichiometry_error%value  
    Else
      Write (message, '(1x,a,f8.3)') 'Error tolerance for the stoichiometry of the generated models (default):',&
                                   & model_data%stoichiometry_error%value  
    End If
    Call info(message,1) 

    ! Print delta_space 
    If (model_data%delta_space%fread) Then
      Write (message,'(1x,a,f5.3,1x,a)') 'Spatial discretization of the simulation cell: ',&
                                           & model_data%delta_space%value, '[Angstrom]'
    Else
      Write (message,'(1x,a,f5.3,1x,a)') 'Spatial discretization of the simulation cell (default): ', &
                                          & model_data%delta_space%value, '[Angstrom]'
    End If
    Call info(message,1) 

    ! Print distance_cutoff
    If (model_data%distance_cutoff%fread) Then
      Write (message,'(1x,a,f4.2,1x,a)') 'Minimum separation distance between species: ',&
                                           & model_data%distance_cutoff%value, '[Angstrom]'
    Else
      Write (message,'(1x,a,f4.2,1x,a)') 'Minimum separation distance between species (default): ', &
                                          & model_data%distance_cutoff%value, '[Angstrom]'
    End If
    Call info(message,1) 

    ! Print format for input/output
    Write (message, '(2(1x,a))') 'Format for input model:    ', Trim(model_data%input_model_format%type)
    Call info(message,1) 
     
    Write (message, '(2(1x,a))') 'Format for output model(s):', Trim(model_data%output_model_format%type)
    Call info(message,1) 

    ! Print in case rotations are prevented 
    If (.Not. model_data%rotate_species%stat) Then
      Write (message, '(1x,a)') 'In case there are molecular species to being incorparated to the model, they will&
                               & not be randomly rotated but keep their orientation, as defined in the corresponding xyz files'
      Call info(message,1) 
    End If

  End Subroutine print_atomistic_settings 

  Subroutine stoichiometry_input_model(files, stoich_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the comparison between the target stoichiometry of the 
    ! pristine sample and the stoichiometry of the input model (table format) 
    !
    ! author    - i.scivetti Feb 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(stoich_type),   Intent(In   ) :: stoich_data
    Type(model_type),    Intent(In   ) :: model_data

    Integer(Kind=wi)   :: i
    Character(Len=256) :: message
    Character(Len= 64) :: fmt0, fmt1, fmt2, fmt3
    Real(Kind=wp)      :: error
    
    fmt0='(15x,a,13x,a))'
    fmt1='(6x,a,2x,4(a,6x))'
    fmt2='(1x,a12,4x,i7,x,f12.3,x,f12.3,x,f12.3)'
    fmt3='(1x,a12,4x,i7,x,f12.3,x,f12.3,x,f12.3,a)'

    Call info(' ', 1)
    Write (message,'(1x,2a)') 'Atomistic details for the input structure provided in file '//Trim(FOLDER_INPUT_GEOM)//&
                            &'/'//Trim(files(FILE_INPUT_STRUCTURE)%filename) 
    Call info(message, 1)
  
    Call info(' ----------------------------------------------------------------', 1)
    Write (message, fmt0)            '|  Total  |',      'Stoichiometry' 
    Call info(message,1)
    Write (message, fmt1) 'Species', '|  number |','Input', 'Target', 'Difference'
    Call info(message,1)
    Call info(' ----------------------------------------------------------------', 1)
    Do i=1, stoich_data%num_species%value
      error=model_data%input%species(i)%s-stoich_data%species(i)%s0_pristine
      If (Abs(error) < model_data%stoichiometry_error%value) Then
        Write (message, fmt2)  Trim(stoich_data%species(i)%tag), model_data%input%species(i)%num, model_data%input%species(i)%s, &
                            stoich_data%species(i)%s0_pristine, error
      Else
        Write (message, fmt3)  Trim(stoich_data%species(i)%tag), model_data%input%species(i)%num, model_data%input%species(i)%s, &
                            stoich_data%species(i)%s0_pristine, error, ' **** EXCEEDS ERROR TOLERANCE'
      End If
      Call info(message,1)
    End Do
    Call info(' ----------------------------------------------------------------', 1)

  End Subroutine stoichiometry_input_model

  Subroutine stoichiometry_sample_model(files, target_stoich, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the comparison between the target stoichiometry 
    ! and the stoichiometry of the generated sample model 
    !
    ! author       - i.scivetti Feb 2021
    ! contibution  - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),  Intent(InOut) :: files(:)
    Real(Kind=wp),    Intent(In   ) :: target_stoich(:)
    Type(model_type), Intent(InOut) :: model_data

    Integer(Kind=wi)   :: i, num, summary
    Character(Len=256) :: message, line
    Character(Len= 64) :: fmt0, fmt1, fmt2, fmt3
    Real(Kind=wp)      :: error
    Logical            :: accept

    accept=.True.
    fmt0='(15x,a,13x,a))'
    fmt1='(6x,a,2x,4(a,6x))'
    fmt2='(1x,a12,4x,i7,x,f12.3,x,f12.3,x,f12.3)'
    fmt3='(1x,a12,4x,i7,x,f12.3,x,f12.3,x,f12.3,a)'
   
    ! Open MODEL_SUMMARY  
    Open(Newunit=files(FILE_MODEL_SUMMARY)%unit_no, File=files(FILE_MODEL_SUMMARY)%filename, Status='Replace')
    summary=files(FILE_MODEL_SUMMARY)%unit_no
 
    Write (line,'(a)') ' ----------------------------------------------------------------'

    Call info(line,1);     Write (summary,'(a)') Trim(line)
    Write (message, fmt0)            '|  Total  |',      'Stoichiometry' 
    Call info(message,1);  Write (summary,'(a)') Trim(message)
    Write (message, fmt1) 'Species', '|  number |','Model', 'Target', 'Difference'
    Call info(message,1);  Write (summary,'(a)') Trim(message)
    Call info(line,1);     Write (summary,'(a)') Trim(line) 

    Do i=1, model_data%types_species 
      error=model_data%sample%species(i)%s-target_stoich(i)
      num=model_data%sample%species(i)%D_num+ model_data%sample%species(i)%num_show
      If (Abs(error) < model_data%stoichiometry_error%value) Then
        Write (message, fmt2)  Trim(model_data%sample%species(i)%tag), num,&
                            model_data%sample%species(i)%s, target_stoich(i), error
      Else
        accept=.False.
        Write (message, fmt3)  Trim(model_data%sample%species(i)%tag), num,&
                            model_data%sample%species(i)%s, target_stoich(i), error, ' **** EXCEEDS ERROR TOLERANCE'
      End If
      Call info(message,1); Write (summary,'(a)') Trim(message)
    End Do
    Call info(line,1); Write (summary,'(a)') Trim(line)
    If (.Not.accept) Then
      If (model_data%stoichiometry_error%fread) Then 
        Write (message, '(a,f5.2,a)') ' *** WARNING: Differences between target and modelled stoichiometries are larger than ', &
                                  & model_data%stoichiometry_error%value, ' (set in the SET_EQCM file).&
                                  & Model is printed, regardless.'
      Else 
        Write (message, '(a,f5.2,a)') ' *** WARNING: Differences between target and modelled stoichiometries are larger than ', &
                                  & model_data%stoichiometry_error%value, ' (default). Model is printed, regardless.'
      End If
      Call info(message, 1); Write (summary,'(a)') Trim(message)
      Write (message, '(a)')    ' To improve the matching between modelled and target stoichiometries either increase:'
      Call info(message, 1); Write (summary,'(a)') Trim(message)
      Write (message, '(a)')    ' 1) the size of the generated model via directive "repeat_input_model", or'
      Call info(message, 1); Write (summary,'(a)') Trim(message)
      Write (message, '(a)')    ' 2) the value of the error tolerance via directive "error_stoichiometry".'
      Call info(message, 1); Write (summary,'(a)') Trim(message)
    Else
      Write (message, '(a)')    ' ***SUCCESS: target and modelled stoichiometries agree within the error tolerance.'
      Call info(message, 1); Write (summary,'(a)') Trim(message)
    End If

    Close(summary)

  End Subroutine stoichiometry_sample_model


  Subroutine compute_area_slab(a, b, area)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the surface area from the 
    ! two vectors that define the 2D periodicity 
    !
    ! author    - i.scivetti Feb 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   ) :: a(3)
    Real(Kind=wp), Intent(In   ) :: b(3)
    Real(Kind=wp), Intent(InOut) :: area 

    Real(Kind=wp) :: cross(3)

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

    ! Calculate area
    area= norm2(cross)

  End Subroutine compute_area_slab


  Subroutine stoichiometry_solutions(ic, jc ,leg, process, files, electrode_data, stoich_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to select stoichiometry solutions among a set of multiple
    ! solutions from the ANALYSIS_EQCM/INTERCALATION_OX or ANALYSIS_EQCM/INTERCALATION_RED files
    !
    ! author    - i.scivetti Feb 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),    Intent(In   ) :: ic
    Integer(Kind=wi),    Intent(In   ) :: jc
    Character(Len=*),    Intent(In   ) :: leg
    Character(Len=*),    Intent(In   ) :: process
    Type(file_type),     Intent(InOut) :: files(:)
    Type(electrode_type),Intent(In   ) :: electrode_data
    Type(stoich_type),   Intent(In   ) :: stoich_data
    Type(model_type),    Intent(InOut) :: model_data


    Character(Len=32 ) :: input_file, set_error
    Character(Len=256) :: word
    Character(Len=256) :: messages(3)
    Character(Len=256) :: message

    Logical          :: safe
    Integer(Kind=wi) :: i, j, k, kf, ix, indx
    Integer(Kind=wi) :: iunit, numsol, io
    Integer(Kind=wi) :: ndim, nlim, types_species
    Integer(Kind=wi) :: fail(7)

    Logical          :: loop
    Real(Kind=wp)    :: s_test, mol, area
    Real(Kind=wp)    :: ds, ds_match
    Integer(Kind=wi) :: nsamples_fine,  nsamples_coarse
    Integer(Kind=wi) :: nmodels

    ! Allocataable arrays
    Real(Kind=wp),    Allocatable  :: solutions(:,:)
    Real(Kind=wp),    Allocatable  :: minmax(:,:)
    Real(Kind=wp),    Allocatable  :: s(:)
    Logical,          Allocatable  :: touched(:)
    Logical,          Allocatable  :: s_comparison(:)
    Integer(Kind=wi), Allocatable  :: list(:)
    Integer(Kind=wi), Allocatable  :: maxlist(:)

    numsol=stoich_data%numsol(jc,ic)
    types_species=stoich_data%num_species%value
    set_error = '***ERROR -'
    Write (messages(1),'(3a,i3,a)') Trim(set_error), Trim(leg), ' of cycle ', ic, ' contains multiple stoichiometry solutions.'

    !!!!!!!!!!!!!!!
    ! INTERCALATION 
    !!!!!!!!!!!!!!!
    If (Trim(process)=='intercalation') Then
      If (numsol==1) Then
        model_data%num_models_sampling=1
        Do i=1, types_species  
          model_data%stoich_target(i,model_data%num_models_sampling)=stoich_data%solution_coeff(i,jc,ic)  
        End Do     
      Else
        ! Prepare the subset of stoichiometry solutions to map the stoichiometry space

        ! Define variables
        ndim=stoich_data%num_indep
        nlim=2**ndim
        ! Define the number of models
        If (model_data%targeted_num_models%fread) Then
          nmodels=model_data%targeted_num_models%value
        Else
          Write (messages(2),'(a,i3,a)') 'The user must define the required number of models to sample the stoichimetric&
                                    & domain via directive "targeted_num_models" in the SET_EQCM file'
          Call info(messages,2)
          Call error_stop(' ')
        End If
       
        ! Allocate variables
        Allocate(solutions(types_species,numsol), Stat=fail(1)) 
        Allocate(touched(numsol),    Stat=fail(2)) 
        Allocate(minmax(2,ndim),     Stat=fail(3))
        Allocate(list(ndim),         Stat=fail(4))
        Allocate(maxlist(ndim),      Stat=fail(5))
        Allocate(s(ndim),            Stat=fail(6))
        Allocate(s_comparison(ndim), Stat=fail(7))

        If (Any(fail > 0)) Then
          Write (message,'(1x,1a)') '***ERROR: Allocation problems for allocatable arrays in subroutine&
                                    & "stoichiometry_solutions". Please check'
          Call error_stop(message)
        End If

        touched=.False.

        !Define is oxidation or reduction
        If (Trim(leg)=='oxidation') Then
           kf=FILE_INTERCALATION_OX
        Else
           kf=FILE_INTERCALATION_RED 
        End If

        input_file = Trim(files(kf)%filename)
        
        ! Open the file
        Inquire(File=Trim(FOLDER_ANALYSIS)//'/'//Trim(files(kf)%filename), Exist=safe)

        If (.not.safe) Then
          Call info(' ', 1)
          Write (messages(1),'(4(1x,a))') Trim(set_error), 'File', Trim(input_file), ' is not accesible'
          Call error_stop(messages(1))
        Else
          Open(Newunit=files(kf)%unit_no, File=Trim(FOLDER_ANALYSIS)//'/'//Trim(files(kf)%filename),Status='old')
          iunit=files(kf)%unit_no
        End If

        ! Read and pass the header
        Read (iunit, Fmt=*, iostat=io) word
        Read (iunit, Fmt=*, iostat=io) word
        Read (iunit, Fmt=*, iostat=io) word
    
        ! Read solutions and store them in a grid
        Do i=1, numsol
         solutions(:,i)=stoich_data%species(:)%s0_pristine
         Read (iunit, Fmt=*, iostat=io) (solutions(stoich_data%index_indep(j),i), j=1, stoich_data%num_indep), &
                                            (solutions(stoich_data%index_dep(j),i),   j=1, 2) 
        End Do

        Close(iunit)

        ! Define minimum and maximum for each independent stoichioetry species
        minmax(1,:)= Huge(1.0_wp)
        minmax(2,:)=-Huge(1.0_wp)
        Do i=1, numsol
          Do j=1, stoich_data%num_indep
            If (solutions(stoich_data%index_indep(j),i) < minmax(1,j)) Then
              minmax(1,j)=solutions(stoich_data%index_indep(j),i)
            End If
            If (solutions(stoich_data%index_indep(j),i) > minmax(2,j)) Then
              minmax(2,j)=solutions(stoich_data%index_indep(j),i)
            End If
          End Do
        End Do

        ! Compute the number of sampling points within the domain determined by minmax
        nsamples_fine=1
        Do j=1, ndim
          nsamples_fine=nsamples_fine*(Int((minmax(2,j)-minmax(1,j))/stoich_data%discretization%value)+1)
        End Do

        ! Define the coarse discretization
        If (nmodels<=numsol) Then
          nsamples_coarse= Nint(1.0_wp*nmodels*nsamples_fine/numsol)-1
          ds = ((nsamples_fine * stoich_data%discretization%value**ndim)/nsamples_coarse)**(1.0_wp/ndim)
          i=1 
          loop=.True.
          Do While (loop)
            ds_match=i*stoich_data%discretization%value
            If ((ds_match-ds)>0) Then
              ds=ds_match
              loop=.False.  
            End If
            i=i+1
          End Do 
        Else
          Write (messages(2),'(3a)') 'The requested number of models is larger than the computed number of solutions&
                                    & in file ', Trim(input_file), '. Please reduce the value of directive "targeted_num_models"'
          Call info(messages,2)
          Call error_stop(' ')
        End If

        ! Define the coarser number of points to be sampled in the space of stoichiometric solutions
        maxlist=1
        Do i=1, ndim
          loop=.True.
          j=1
          Do While (loop)
            s_test=j*ds+minmax(1,i)
            If (s_test<minmax(2,i)) Then
              maxlist(i)=maxlist(i)+1
              j=j+1
            Else
              loop=.False.
            End If
          End Do
        End Do

        ! Find the coarser set of solutions
        list=0
        ix=1
        model_data%num_models_sampling=0
        Do While (ix /= 0)
          list(ix)=list(ix)+1 
          indx= stoich_data%index_indep(ix)
          If (list(ix) <= maxlist(ix) .And. list(ix) >= 1 ) Then
          s(ix)=minmax(1,ix)+(list(ix)-1)*ds
            If (ix == stoich_data%num_indep) Then
              i=1
              loop=.True.
              Do While (i<= numsol .And. loop)
                If (.Not. touched(i)) Then
                  s_comparison=.False.
                  Do j=1, ndim 
                    If (Abs(solutions(stoich_data%index_indep(j),i)-s(j))<stoich_data%discretization%value/2.0_Wp) Then
                      s_comparison(j)=.True.
                    End If
                  End Do
                  If (All(s_comparison)) Then
                    touched(i)=.True.
                    model_data%num_models_sampling= model_data%num_models_sampling+1
                    loop=.False. 
                    Do j=1, types_species                
                      model_data%stoich_target(j,model_data%num_models_sampling)=solutions(j,i)
                    End Do 
                  End If
                End If
                i=i+1
              End Do
            Else
              ix=ix+1
            End If
          Else
            list(ix) = 0
            ix = ix - 1
          End If
        End Do

       ! Add solutions from the corners if any
        maxlist=2
        list=0
        ix=1

        Do While (ix /= 0)
          list(ix)=list(ix)+1 
          indx= stoich_data%index_indep(ix)
          If (list(ix) <= maxlist(ix) .And. list(ix) >= 1 ) Then
          s(ix)=minmax(list(ix),ix)
          
            If (ix == stoich_data%num_indep) Then
              i=1
              loop=.True.
              Do While (i<= numsol .And. loop)
                If (.Not. touched(i)) Then
                  s_comparison=.False.
                  Do j=1, ndim 
                    If (Abs(solutions(stoich_data%index_indep(j),i)-s(j))<stoich_data%discretization%value/2.0_wp) Then
                      s_comparison(j)=.True.
                    End If
                  End Do
                  If (All(s_comparison)) Then
                    touched(i)=.True.
                    model_data%num_models_sampling= model_data%num_models_sampling+1
                    loop=.False. 
                    Do j=1, types_species                
                      model_data%stoich_target(j,model_data%num_models_sampling)=solutions(j,i)
                    End Do 
                  End If
                End If
                i=i+1
              End Do
            Else
              ix=ix+1
            End If
          Else
            list(ix) = 0
            ix = ix - 1
          End If
        End Do
        ! Print selected solutons to file SELECTED_SOLUTIONS
        kf=FILE_SELECTED_SOLUTIONS

        input_file = Trim(files(kf)%filename)

        Open(Newunit=files(kf)%unit_no, File=files(kf)%filename,Status='Replace')
        iunit=files(kf)%unit_no

        ! Write the header
        If (stoich_data%block_constraints%fread) Then
          Write (iunit,'(a1,*(a12,2x))') '#',&
                                   (Trim(stoich_data%species(stoich_data%index_indep(k))%tag), k=1,stoich_data%num_indep),&
                                   (Trim(stoich_data%species(stoich_data%index_dep(k))%tag),   k=1,stoich_data%num_dep),&
                                   (Trim(stoich_data%species(stoich_data%index_const(k))%tag), k=1,stoich_data%net_const)
        Else
          Write (iunit,'(a1,*(a12,2x))') '#',&
                                   (Trim(stoich_data%species(stoich_data%index_indep(k))%tag), k=1,stoich_data%num_indep),&
                                   (Trim(stoich_data%species(stoich_data%index_dep(k))%tag),   k=1,stoich_data%num_dep)
        End If
        Write (iunit,'(a,i3)')     '# Cycle number: ', ic
        Write (iunit,'(a,f12.6)')  '# Electrons: ', stoich_data%electrons 
        Do i=1, model_data%num_models_sampling
          If (stoich_data%block_constraints%fread) Then
            Write (iunit,'(1x,*(f12.6,2x))') (model_data%stoich_target(stoich_data%index_indep(j),i),  j=1,stoich_data%num_indep),&
                                 &          (model_data%stoich_target(stoich_data%index_dep(j),i),    j=1,stoich_data%num_dep),  &
                                 &          (model_data%stoich_target(stoich_data%index_const(j),i),  j=1,stoich_data%net_const)
          Else
            Write (iunit,'(1x,*(f12.6,2x))') (model_data%stoich_target(stoich_data%index_indep(j),i),  j=1,stoich_data%num_indep),&
                                          & (model_data%stoich_target(stoich_data%index_dep(j),i),    j=1,stoich_data%num_dep)
   
           End If
        End Do 
        Close(iunit)
        ! Move file
         Call execute_command_line('mv '//Trim(Trim(files(FILE_SELECTED_SOLUTIONS)%filename))//' '//Trim(FOLDER_ANALYSIS)) 

        ! Deallocate arrays
        Deallocate(s_comparison)
        Deallocate(s)
        Deallocate(maxlist)
        Deallocate(minmax)
        Deallocate(touched)
        Deallocate(solutions)
      End If

    !!!!!!!!!!!!!!!!!!!
    ! ELECTRODEPOSITION
    !!!!!!!!!!!!!!!!!!!
    ElseIf (Trim(process)=='electrodeposition') Then
      model_data%num_models_sampling=1
      ! Calcuate the amount of moles for the modelled area
      If (ic==1 .And. jc==1) Then 
        Do i=1, types_species
          model_data%stoich_target(i,1)=stoich_data%species(i)%s0_pristine
        End Do
      End If 
      Do i=1, stoich_data%num_indep
        j=stoich_data%index_indep(i)
        mol=stoich_data%solution_moles(jc,ic)
        area=model_data%input%slab_area
        ds=mol*area/electrode_data%area_geom%value*avogadro/cm_to_Ang**2
        ds=ds/model_data%input%min_species/model_data%input%min_components*model_data%input%min_stoich
        model_data%stoich_target(j,1)=model_data%stoich_target(j,1)+ds
      End Do
    End If

  End Subroutine stoichiometry_solutions

  Subroutine reset_to_input(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Reset to values and lists found after the execution of subroutine
    ! identify_species_input 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),  Intent(InOut) :: T
 
    Integer(Kind=wi) :: i, j

    Do i = 1, T%types_species
      Do j = 1, T%input%species(i)%num_extra
        T%input%species(i)%units(j)%vanish=.False.
      End Do
    End Do

    Do i=1, T%input%num_atoms_extra
      T%input%atom(i)%vanish=.False.
      T%input%atom(i)%in_species=.False.
    End Do

  End Subroutine reset_to_input


  Subroutine print_output_files(type_sim, files, model_data, simulation_data, hpc_data, ic, kc, leg)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print output files with:
    !  - the generated sample models
    !  - all related files with settings for the simulation
    !
    ! author         - i.scivetti Jan 2021
    ! contribution   - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=32), Intent(In   ) :: type_sim
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(simul_type),  Intent(InOut) :: simulation_data
    Type(hpc_type),    Intent(In   ) :: hpc_data
    Integer(Kind=wi),  Intent(In   ) :: ic
    Integer(Kind=wi),  Intent(In   ) :: kc
    Character(Len=*),  Intent(In   ) :: leg

    Character(Len=256) :: exec_mkdir, exec_cp, exec_mv, exec_cat
    Character(Len=256) :: fileformat
    Character(Len=256) :: filename, path, modchar
    Character(Len=3)   :: ic_char

    fileformat=model_data%output_model_format%type
  
    filename='SAMPLE.'//Trim(fileformat)

    If (Trim(fileformat)=='vasp') Then
      Call print_vasp_output(files, model_data%sample, model_data%types_species, model_data%selective_dyn, simulation_data)
    Else If (Trim(fileformat)=='cp2k') Then
      Call print_cp2k_output(files, model_data%sample, model_data%types_species, simulation_data)
    Else If (Trim(fileformat)=='castep') Then
      Call print_castep_output(files, model_data%sample, model_data%types_species, simulation_data)
    Else If (Trim(fileformat)=='onetep') Then
      Call print_onetep_output(files, model_data%sample, model_data%types_species, simulation_data)
    Else If (Trim(fileformat)=='siesta') Then
      Call print_siesta_output(files, model_data%sample, model_data%types_species)
    Else If (Trim(fileformat)=='xyz') Then 
      Call print_xyz_output(files, model_data%sample, model_data%types_species)
    Else If (Trim(fileformat)=='cif') Then 
      Call print_cif_output(files, model_data%sample, model_data%types_species, model_data%selective_dyn)
    End If 

    Call execute_command_line('[ ! -d '//Trim(FOLDER_MODELS)//' ] && '//'mkdir '//Trim(FOLDER_MODELS))

    If (Trim(type_sim)=='model_pristine_sample') Then
      Write (path,'(a)') Trim(FOLDER_MODELS)//'/pristine' 
      exec_mkdir='mkdir '//Trim(path)
      Call execute_command_line('[ ! -d '//Trim(path)//' ] && '//Trim(exec_mkdir))
    Else If (Trim(type_sim)=='model_disordered_system') Then
      Write (path,'(a)') Trim(FOLDER_MODELS)//'/disordered' 
      exec_mkdir='mkdir '//Trim(path)
      Call execute_command_line('[ ! -d '//Trim(path)//' ] && '//Trim(exec_mkdir))
    Else If (Trim(type_sim) == 'model_cycled_sample') Then
      Write (ic_char,'(i3)') ic
      Write (path,'(a)') Trim(FOLDER_MODELS)//'/'//Trim(Adjustl(ic_char))//'cycle-'//Trim(leg)
      exec_mkdir='mkdir '//Trim(path)
      Call execute_command_line('[ ! -d '//Trim(path)//' ] && '//Trim(exec_mkdir))
      If (model_data%num_models_sampling/=1) Then
        Write (modchar,'(i4)') kc
        Write (path,'(a)') Trim(FOLDER_MODELS)//'/'//Trim(Adjustl(ic_char))//&
                             & 'cycle-'//Trim(leg)//'/model'//Adjustl(Trim(modchar))
        exec_mkdir='mkdir '//Trim(path)
        Call execute_command_line('[ ! -d '//Trim(path)//' ] && '//Trim(exec_mkdir))
      End If
    End If

    Call info(' Atomistic model is printed to file '//Trim(path)//'/'//Trim(filename), 1)

    ! Moving simulation setting files to the corresponding directory
    If (simulation_data%generate) Then
      Call info(' Files for the simulation of the generated model are printed to folder '//Trim(path), 1)
      If (Trim(fileformat)=='vasp') Then
        ! Generate POSCAR
        exec_cp='cp '//Trim(files(FILE_OUTPUT_STRUCTURE)%filename)//' '//Trim(path)//'/POSCAR'
        Call execute_command_line(exec_cp)
      Else If (model_data%output_model_format%type=='castep') Then
        exec_cat='cat '//Trim(files(FILE_SET_SIMULATION)%filename)//' '//&
              &Trim(files(FILE_OUTPUT_STRUCTURE)%filename)//' > model.cell'
        Call execute_command_line(exec_cat)
        exec_mv ='mv  '//Trim(files(FILE_SET_SIMULATION)%filename)//' CI-castep.cell' 
        Call execute_command_line(exec_mv)
      Else If (model_data%output_model_format%type=='onetep') Then
        exec_cat='cat '//Trim(files(FILE_SET_SIMULATION)%filename)//' '//&
              &Trim(files(FILE_OUTPUT_STRUCTURE)%filename)//' > model.dat'
        Call execute_command_line(exec_cat)
        exec_mv ='mv  '//Trim(files(FILE_SET_SIMULATION)%filename)//' CI-onetep.dat'  
        Call execute_command_line(exec_mv)
      End If

      Call print_simulation_files(path, files, model_data, simulation_data, hpc_data)
    End If

    ! MOving structure
    exec_mv='mv '//Trim(files(FILE_OUTPUT_STRUCTURE)%filename)//' '//Trim(path)//'/'//Trim(filename)
    Call execute_command_line(exec_mv)
      
    ! Move MODEL_SUMMARY files
    exec_mv='mv '//Trim(files(FILE_MODEL_SUMMARY)%filename)//' '//Trim(path)//'/'//Trim(files(FILE_MODEL_SUMMARY)%filename)
    Call execute_command_line(exec_mv)

    model_data%sample%path=path
    If (Trim(type_sim) == 'model_cycled_sample')Call info(' ',1)

  End Subroutine print_output_files 

  Subroutine print_simulation_files(path, files, model_data, simulation_data, hpc_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print output files with settings for the simulation of the models
    ! This abstraction is convenient to generated files with simulation settings
    ! without the need of re-generating atomistic models (the most expensive bit)
    !
    ! author         - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=256), Intent(In   ) :: path
    Type(file_type),    Intent(InOut) :: files(:)
    Type(model_type),   Intent(InOut) :: model_data
    Type(simul_type),   Intent(In   ) :: simulation_data
    Type(hpc_type),     Intent(In   ) :: hpc_data

    Character(Len=256) :: exec_cp, exec_mv, fileformat
    Character(Len=256) :: path_pp

    fileformat=model_data%output_model_format%type

    !!!! VASP
    If (Trim(fileformat)=='vasp') Then
      ! Generate INCAR
      exec_mv='mv '//Trim(files(FILE_SET_SIMULATION)%filename)//' '//Trim(path)//'/INCAR'   
      Call execute_command_line(exec_mv)
      ! Generate KPOINTS 
      exec_mv='mv '//Trim(files(FILE_KPOINTS)%filename)//' '//Trim(path)//'/KPOINTS'   
      Call execute_command_line(exec_mv)
      ! Generate POTCAR
      If (simulation_data%dft%pp_info%stat)Then
        exec_mv='mv POTCAR'//' '//Trim(path)//'/POTCAR'   
        Call execute_command_line(exec_mv)
      End If
    !!!! CP2K
    Else If (Trim(fileformat)=='cp2k') Then
      ! Generated input.cp2k 
      exec_mv='mv '//Trim(files(FILE_SET_SIMULATION)%filename)//' '//Trim(path)//'/input.cp2k'
      Call execute_command_line(exec_mv)
      ! copy the potential
      If (simulation_data%dft%pp_info%stat)Then
        path_pp=Trim(FOLDER_DFT)//'/PPs/'//Trim(simulation_data%dft%pseudo_pot(1)%file_name)  
        exec_cp='cp '//Trim(path_pp)//' '//Trim(path)
        Call execute_command_line(exec_cp)
      End If
      ! copy basis set 
      If (simulation_data%dft%basis_info%stat)Then
        path_pp=Trim(FOLDER_DFT)//'/BASIS_SET'
        exec_cp='cp '//Trim(path_pp)//' '//Trim(path)
        Call execute_command_line(exec_cp)
      End If
    !!!! CASTEP
    Else If (Trim(fileformat)=='castep') Then
      exec_mv='mv '//'model.param'//' '//Trim(path) 
      Call execute_command_line(exec_mv)
      exec_mv='mv '//'model.cell'//' '//Trim(path) 
      Call execute_command_line(exec_mv)
      exec_mv='mv '//'CI-castep.cell'//' '//Trim(path) 
      Call execute_command_line(exec_mv)
      If (simulation_data%dft%pp_info%stat)Then
        exec_cp='cp DFT/PPs/* '//Trim(path)
        Call execute_command_line(exec_cp)
      End If
    !!!! ONETEP
    Else If (Trim(fileformat)=='onetep') Then
      exec_mv='mv '//'model.dat'//' '//Trim(path) 
      Call execute_command_line(exec_mv)
      exec_mv='mv '//'CI-onetep.dat'//' '//Trim(path) 
      Call execute_command_line(exec_mv)
      If (simulation_data%dft%pp_info%stat) Then
        exec_cp='cp DFT/PPs/* '//Trim(path)
        Call execute_command_line(exec_cp)
      End If
    End If

    ! If vdW correction uses kernel, copy the file
    If (simulation_data%dft%need_vdw_kernel) Then
      exec_cp='cp '// Trim(FOLDER_DFT)//'/'//Trim(simulation_data%dft%vdw_kernel_file)//' '//Trim(path)//'/'   
      Call execute_command_line(exec_cp)
    End If

    ! If hpc settings are defined, copy to the corresponding HPC file
    If (hpc_data%generate) Then   
      exec_cp='cp '//Trim(files(FILE_HPC_SETTINGS)%filename)//' '//Trim(path)//'/'//Trim(hpc_data%script_name)
      Call execute_command_line(exec_cp)
    End If

  End Subroutine print_simulation_files 

  Subroutine create_list_net_elements(T, types_species)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Create list for the atoms 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     

    Integer(Kind=wi) :: icount, isp, j, i, inat

    ! Define the amount of net elements
    icount=0
    T%list%net_elements=0
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        inat=0
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. (Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) ) Then
            inat=inat+1
          End If
        End Do
        If (inat/=0) Then
          icount=icount+1
          T%list%N0(icount)=inat
          T%list%tag(icount)=Trim(T%species(isp)%component%tag(j))
          T%list%element(icount)=Trim(T%species(isp)%component%element(j))   
          T%list%net_elements=T%list%net_elements+1
        End If
      End Do
    End Do

  End Subroutine create_list_net_elements

  Subroutine print_cif_output(files, T, types_species, dynamics)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print file in cif format 
    !
    ! author    - i.scivetti Jan 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     
    Logical,           Intent(In   ) :: dynamics

    Integer(Kind=wi) :: isp, i, j
    Integer(Kind=wi) :: iunit
    Character(Len=256)   :: message, exec_cif, exec_rm

    exec_rm='rm POSCAR'
    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File='POSCAR', Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    Write (iunit,'(a)') 'Model for sample' 
    Write (iunit,'(f19.16)') T%scale_factor_vasp
    Do i= 1, 3
      Write (iunit,'(3f20.12)') (T%cell(i,j), j=1,3)
    End Do
    Write (iunit,'(*(6x,a2))') (Trim(T%list%element(j)), j=1, T%list%net_elements)   
    Write (iunit,'(*(i8))') (T%list%N0(j), j=1, T%list%net_elements)
    If (dynamics) Then
      Write (iunit, '(a)')  'Selective dynamics'
    End If   
    Write (iunit, '(a)')    'Cartesian'
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            If (dynamics) Then
              Write (iunit, '(3f20.12,3l3)') T%atom(i)%r, T%atom(i)%dynamics
            Else
              Write (iunit, '(3f20.12)') T%atom(i)%r
            End If
          End If
        End Do
      End Do
    End Do

    ! Close file
    Close(iunit)

    exec_cif= 'python cif-output-geom.py'
    Call execute_command_line(exec_cif, Exitstat=i)
    If (i/=0) Then
      Call info(' ', 1)
      Write (message,'(1x,a)') '***ERROR: Are Python and ASE installed in your machine?'
      Call execute_command_line(exec_rm)
      Call error_stop(message)      
    End If   

    Call execute_command_line(exec_rm)

  End Subroutine print_cif_output

  Subroutine print_vasp_output(files, T, types_species, dynamics, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print files in vasp format 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     
    Logical,           Intent(In   ) :: dynamics
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: isp, i, j
    Integer(Kind=wi) :: iunit


    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    Write (iunit,'(a)') 'Model for sample' 
    Write (iunit,'(f19.16)') T%scale_factor_vasp
    Do i= 1, 3
      Write (iunit,'(3f20.12)') (T%cell(i,j), j=1,3)
    End Do
    Write (iunit,'(*(6x,a2))') (Trim(T%list%element(j)), j=1, T%list%net_elements)   
    Write (iunit,'(*(i8))') (T%list%N0(j), j=1, T%list%net_elements)
    If (dynamics) Then
      Write (iunit, '(a)')  'Selective dynamics'
    End If   
    Write (iunit, '(a)')    'Cartesian'
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            If (dynamics) Then
              Write (iunit, '(3f20.12,3l3)') T%atom(i)%r, T%atom(i)%dynamics
            Else
              Write (iunit, '(3f20.12)') T%atom(i)%r
            End If
          End If
        End Do
      End Do
    End Do

    ! Close file
    Close(iunit)

    If (simulation_data%generate) Then
      Call print_vasp_settings(files, T%list%net_elements, T%list%element, T%list%tag, T%list%N0, simulation_data) 
    End If

  End Subroutine print_vasp_output

  Subroutine print_cp2k_output(files, T, types_species, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print files in CP2K format 
    !
    ! author        - i.scivetti May  2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     
    Type(simul_type),  Intent(In   ) :: simulation_data

    Integer(Kind=wi) :: isp, i, j, ntot
    Integer(Kind=wi) :: iunit

    ntot=0
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            ntot=ntot+1
          End If
        End Do
      End Do
    End Do

    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    ! print atomic position
    Write (iunit,'(i10)') ntot
    Write (iunit,*) '  ' 

    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            Write (iunit,'(a,2x,3f20.12)') Trim(T%atom(i)%tag), T%atom(i)%r
          End If
        End Do
      End Do
    End Do

    Close(iunit)

    If (simulation_data%generate) Then
      Call print_cp2k_settings(files, T%list%net_elements, T%list%element, T%list%tag, T%list%N0, simulation_data)
    End If

  End Subroutine print_cp2k_output 

  Subroutine print_castep_output(files, T, types_species, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print files in CASTEP format 
    !
    ! author        - i.scivetti June 2021 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: isp, i, j, k
    Integer(Kind=wi) :: iunit
    Logical :: loop 
    Character(Len=8) :: atom_tag  

    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    Write (iunit,'(a)') '%BLOCK LATTICE_CART'
    Do i= 1, 3
      Write (iunit,'(3f20.12)') (T%cell(i,j), j=1,3)
    End Do  
    Write (iunit,'(a)') '%ENDBLOCK LATTICE_CART' 
    Write (iunit, *) '  '

    ! print atomic position
    Write (iunit,'(a)') '%BLOCK POSITIONS_ABS'

    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            If (simulation_data%dft%mag_info%fread) Then
                k=1
                loop=.True.
                Do While (k <= simulation_data%total_tags .And. loop)
                  If (Trim(T%atom(i)%tag)==Trim(simulation_data%dft%magnetization(k)%tag)) Then
                    If (simulation_data%dft%pp_info%fread) Then
                      atom_tag=Trim(T%atom(i)%tag)
                    Else
                      atom_tag=Trim(T%atom(i)%element)
                    End If
                    Write (iunit,'(a,2x,3f20.12,4x,a,f5.2)') Trim(atom_tag), T%atom(i)%r, 'spin=', &
                                                        & simulation_data%dft%magnetization(k)%value 
                    loop=.False.
                  End If
                  k=k+1
                End Do
            Else
              If (simulation_data%dft%pp_info%fread) Then
                atom_tag=Trim(T%atom(i)%tag)
              Else
                atom_tag=Trim(T%atom(i)%element)
              End If
              Write (iunit,'(a,2x,3f20.12)') Trim(atom_tag), T%atom(i)%r
            End If
          End If
        End Do
      End Do
    End Do
    Write (iunit,'(a)') '%ENDBLOCK POSITIONS_ABS'


    Close(iunit)

    If (simulation_data%generate) Then
      Call print_castep_settings(files, T%list%net_elements, T%list%element, T%list%tag, T%list%N0, simulation_data)
    End If

  End Subroutine print_castep_output 

  Subroutine print_onetep_output(files, T, types_species, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print files in CASTEP format
    !
    ! author        - i.scivetti Jan  2021 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: isp, i, j
    Integer(Kind=wi) :: iunit

    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    ! Print simulation cell
    Write (iunit,'(a)') '%block lattice_cart'
    Write (iunit,'(a)') 'ang'
    Do i= 1, 3
      Write (iunit,'(3f12.6)') (T%cell(i,j), j=1,3)
    End Do  
    Write (iunit,'(a)') '%endblock lattice_cart' 
    Write (iunit, *) '  '

    ! print atomic position
    Write (iunit,'(a)') '%block positions_abs'
    Write (iunit,'(a)') 'ang'
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
              Write (iunit,'(a,2x,3f12.6)') Trim(T%atom(i)%tag), T%atom(i)%r
          End If
        End Do
      End Do
    End Do
    Write (iunit,'(a)') '%endblock positions_abs' 

    Close(iunit)

    If (simulation_data%generate) Then
      Call print_onetep_settings(files, T%list%net_elements, T%list%tag, T%list%N0, simulation_data)
    End If

  End Subroutine print_onetep_output 

  Subroutine print_siesta_output(files, T, types_species)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print files in SIESTA format 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     

    Integer(Kind=wi) :: isp, i, j, ntot, nat 
    Integer(Kind=wi) :: iunit
    Logical :: lpass

    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    ! Print lattice vectors
    Write (iunit, *) '  '
    Write (iunit,'(a)') 'LatticeConstant 1.00 Ang'
    Write (iunit, *) '  '
    Write (iunit,'(a)') '%block LatticeVectors'
    Do i= 1, 3
      Write (iunit,'(3f12.6)') (T%cell(i,j), j=1,3)
    End Do  
    Write (iunit,'(a)') '%endblock LatticeVectors'
    Write (iunit, *) '  '
    Write (iunit, *) '  '

    Write (iunit,'(a)') '%block ChemicalSpeciesLabel'
    ntot=0
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        nat=0
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            nat=nat+1 
          End If
        End Do
        If (nat/=0) Then
          ntot=ntot+1
          Write (iunit,'(i3,4x,i4,4x,a2)')  ntot, T%species(isp)%component%atomic_number(j), T%species(isp)%component%element(j) 
        End If  
      End Do
    End Do
    Write (iunit,'(a)') '%block ChemicalSpeciesLabel'

    Write (iunit, *) '  '
    Write (iunit, *) '  '
    Write (iunit,'(a)') 'AtomicCoordinatesFormat Ang'
    Write (iunit,'(a)') 'AtomCoorFormatOut       Ang'
    Write (iunit, *) '  '
    Write (iunit,'(a)') '%block positions_abs'
    
    ntot=0
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        lpass=.True.
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
              If (lpass) Then
                lpass=.False.
                ntot=ntot+1
              EndIf
              Write (iunit, '(3f16.8,3x,i3,5x,a2)') T%atom(i)%r, ntot, T%atom(i)%element
            End If
        End Do
      End Do
    End Do

    Write (iunit,'(a)') '%endblock positions_abs'
    Write (iunit, *) '  '
    Close(iunit)

  End Subroutine print_siesta_output


  Subroutine print_xyz_output(files, T, types_species)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prints model in xyz format. 
    !
    ! author        - i.scivetti July  2021 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(sample_type), Intent(InOut) :: T
    Integer(Kind=wi),  Intent(In   ) :: types_species     

    Integer(Kind=wi) :: isp, i, j, ntot
    Integer(Kind=wi) :: iunit

    ntot=0
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
            ntot=ntot+1
          End If
        End Do
      End Do
    End Do

    ! Open OUTPUT_STRUCTURE file
    Open(Newunit=files(FILE_OUTPUT_STRUCTURE)%unit_no, File=files(FILE_OUTPUT_STRUCTURE)%filename,Status='Replace')
    iunit=files(FILE_OUTPUT_STRUCTURE)%unit_no

    ! print atomic position
    Write (iunit,'(i10)') ntot
    Write (iunit,'(a, 9f10.4, a)') 'Lattice = "' , ((T%cell(i,j), j=1,3), i=1,3), '"' 
    Do isp = 1, types_species
      Do j=1, T%species(isp)%num_components
        Do i=1, T%num_atoms
          If ((.Not. T%atom(i)%vanish) .And. &
            & Trim(T%atom(i)%tag)==Trim(T%species(isp)%component%tag(j))) Then
              Write (iunit,'(a,2x,3f12.6)') Trim(T%atom(i)%element), T%atom(i)%r
          End If
        End Do
      End Do
    End Do

    Close(iunit)

  End Subroutine print_xyz_output 

  Subroutine match_target_stoichiometry(T, target_stoich, option)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Find the difference in the number of species between the target stoichiometry
    ! and the input/sample model
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),  Intent(InOut) ::  T
    Real(Kind=wp),     Intent(In   ) ::  target_stoich(T%types_species)
    Character(Len=*),  Intent(In   ) ::  option
  
    Integer(Kind=wi) :: i, num_target, min_comp

    T%remove_species=.False.
    T%insert_species=.False.

    ! Find the number of:
    ! 1) species that the pristine model would have
    ! 2) species to add/remove from the input model to match the units of the target model 
    Do i = 1, T%types_species
      If (Trim(option)=='input') Then
        !Initialis
        T%input%species(i)%D_num=0
        min_comp=1
        If (T%input%species(i)%change_content) Then
          min_comp=T%input%min_components
        End If
        num_target=Floor(target_stoich(i)*T%input%min_species*min_comp/T%ref_stoich)
        If (T%input_times*Abs(T%input%species(i)%s - target_stoich(i)) >  T%stoichiometry_error%value)Then
         ! If (Abs(target_stoich(i)*T%min_species*min_comp-num_target) > epsilon(1.0_wp)) Then 
            num_target=num_target+1
         ! End If
        End If
        T%input%species(i)%D_num=num_target-T%input%species(i)%num
        If (T%input%species(i)%D_num<0) Then
           T%remove_species=.True.
        End If
        If (T%input%species(i)%D_num>0) Then
           T%insert_species=.True.
        End If
      ElseIf (Trim(option)=='output') Then
        !Initialis
        T%sample%species(i)%D_num=0
        min_comp=1
        If (T%sample%species(i)%change_content) Then
          min_comp=T%sample%min_components
        End If
        num_target=Nint(target_stoich(i)*T%sample%min_species*min_comp/T%ref_stoich)
        T%sample%species(i)%D_num=num_target-T%sample%species(i)%num_show
        T%sample%species(i)%s=Real(T%sample%species(i)%D_num+T%sample%species(i)%num_show)/T%sample%min_species*T%ref_stoich
        T%sample%species(i)%s=T%sample%species(i)%s/min_comp
        If (T%sample%species(i)%D_num<0) Then
           T%remove_species=.True.
        End If
        If (T%sample%species(i)%D_num>0) Then
           T%insert_species=.True.
        End If
      End If


    End Do   

  End Subroutine match_target_stoichiometry
   

  Subroutine define_repeated_model(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Uses the modified input model to build the output sample
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: i, j, k
    Integer(Kind=wi) :: ic(3), ic1, ic2, ic3
    Integer(Kind=wi) :: na, nat_extra, nsp_extra
    Integer(Kind=wi) :: ind_at, ind_sp, isp 
    Integer(Kind=wi) :: ris(3)
    Real(Kind=wp)    :: shift(3)
    Character(Len=256) :: messages(3)

    ris=model_data%repeat_input_model%value
    nat_extra=model_data%input%num_atoms_extra

    model_data%sample%min_species=model_data%input%min_species*model_data%input_times
    model_data%sample%min_components=model_data%input%min_components

    ! Copy species
    Do i = 1, model_data%types_species    
      model_data%sample%species(i)%num=model_data%input%species(i)%num_extra*model_data%input_times
    End Do

    ! Copy the first num_atoms_extra atoms
    If (model_data%input%num_atoms_extra>model_data%sample%max_atoms) Then
       Write (messages(1),'(1x,1a)') '***ERROR: the initial allocated dimensions for the maximum number of atoms&
                                    & to generate models has been exceeded. '      
       Write (messages(2),'(1x,1a)') '   The user must define a larger value for the "multiple_output_atoms"&
                                    & directive, which is set to 2 by default.'
       Write (messages(3),'(1x,1a)') '   Try first by increasing this directive to 5 first. Increase it more if needed'
       Call info(messages, 3)
       Call error_stop(' ')
    End If            
    Do j =1 , model_data%input%num_atoms_extra
      model_data%sample%atom(j)=model_data%input%atom(j)      
    End Do
  
    Do ic3 = 1, ris(3)
      ic(3)=ic3
      Do ic2 = 1, ris(2)
        ic(2)=ic2
        Do ic1 = 1, ris(1)
          ic(1)=ic1
          ind_at=nat_extra*(ic1-1)+nat_extra*ris(1)*(ic2-1)+nat_extra*ris(1)*ris(2)*(ic3-1)
          ! Set atomic properties
          shift=0.0_wp
          Do j=1,3
            Do k=1,3
              shift(j)=shift(j)+model_data%input%cell(k,j)*(ic(k)-1)
            End Do
          End Do
          Do i=1, nat_extra
            na=i+ind_at
            model_data%sample%atom(na)%r=model_data%input%atom(i)%r+shift
            model_data%sample%atom(na)%element=model_data%input%atom(i)%element
            model_data%sample%atom(na)%atomic_number=model_data%input%atom(i)%atomic_number
            model_data%sample%atom(na)%tag=model_data%input%atom(i)%tag
            model_data%sample%atom(na)%in_species=model_data%input%atom(i)%in_species
            model_data%sample%atom(na)%vanish=model_data%input%atom(i)%vanish
            model_data%sample%atom(na)%dynamics=model_data%input%atom(i)%dynamics
          End Do
          !Set species
          Do isp = 1, model_data%types_species
            nsp_extra=model_data%input%species(isp)%num_extra
            ind_sp=nsp_extra*(ic1-1)+nsp_extra*ris(1)*(ic2-1)+nsp_extra*ris(1)*ris(2)*(ic3-1)
            If (model_data%input%species(isp)%change_content) Then
              Do j=1, model_data%input%species(isp)%num_extra
                If ((j+ind_sp) > max_num_species_units) Then
                  Call error_stop('***ERROR: requested size for the model structure exceeeds the stipulated limits.&
                                 & Please reduce the size using directive "repeat_input_model"') 
                Else
                   model_data%sample%species(isp)%units(j+ind_sp)%vanish = &
                        model_data%input%species(isp)%units(j)%vanish
                End If
                Do k= 1, model_data%input%species(isp)%atoms_per_species
                  model_data%sample%species(isp)%units(j+ind_sp)%list(k) = &
                                                    model_data%input%species(isp)%units(j)%list(k)+ind_at 
                End Do
              End Do
            End If
          End Do 
  
        End Do  
      End Do  
    End Do 
 
    Do isp = 1, model_data%types_species
      model_data%sample%species(isp)%num_show=0
      Do j = 1, model_data%sample%species(isp)%num
        If (.Not. model_data%sample%species(isp)%units(j)%vanish) Then
          model_data%sample%species(isp)%num_show=model_data%sample%species(isp)%num_show+1
        End If
      End Do
    End Do

  End Subroutine define_repeated_model


  Subroutine insert_species(process, stoich_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to insert species into the model  
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=32),   Intent(In   ) :: process
    Type(stoich_type),   Intent(In   ) :: stoich_data
    Type(model_type),    Intent(InOut) :: model_data

    Logical :: loop
    Integer(Kind=wi) :: i, j, index
    Real(Kind=wp)    :: min_diff, diff, length

    Do i=1, model_data%types_species  
      model_data%input%species(i)%num_extra=model_data%input%species(i)%num
    End Do

    ! We first try to insert molecules, then atoms
    Do i=1, model_data%types_species 
      If (model_data%input%species(i)%D_num>0) Then
        If (model_data%input%species(i)%topology=='molecule') Then
          Call read_molecular_units(model_data%input%species(i), stoich_data%species(i)%bond_cutoff)  
        ElseIf (model_data%input%species(i)%topology=='atom') Then
          model_data%input%species(i)%definition%r0(1,:)=0.0_wp 
          model_data%input%species(i)%definition%tag(1)=model_data%input%species(i)%component%tag(1)
          model_data%input%species(i)%definition%element(1)=model_data%input%species(i)%component%element(1)
          model_data%input%species(i)%definition%atomic_number(1)=model_data%input%species(i)%component%atomic_number(1)
        End If
      End If
    End Do

    ! Find the total number of grid points based on the delta_space
    Do i=1,3
      model_data%npoints(i)=Int(model_data%input%cell_length(i)/model_data%delta_space%value)
      model_data%scan(i)=1 
    End Do

    If (Trim(process)=='electrodeposition') Then
      If (Trim(model_data%normal_vector%type)=='c3') Then
        index=3
        length=model_data%input%cell_length(3)
      ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
        index=2
        length=model_data%input%cell_length(2)
      ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
        index=1
        length=model_data%input%cell_length(1)
      End If
      
      min_diff=Huge(1.0_wp)
      Do i=1, model_data%npoints(index)
        diff=Abs((i-1)*length/model_data%npoints(index)-Abs(model_data%deposition_level%value))
        If (diff < min_diff) Then
          min_diff=diff
          model_data%scan(index)=i
        End If
      End Do

    End If

    j=1
    loop=.True.
    Do While (loop)
      loop=.False.
      Do i=1, model_data%types_species 
        If (model_data%input%species(i)%change_content) Then
           If (j<= model_data%input%species(i)%D_num) Then  
              loop=.True.
              If (Trim(process)=='intercalation') Then 
                Call fitting_species_intercalation(model_data%input, model_data%npoints,&
                                  & i, j, model_data%distance_cutoff, model_data%rotate_species%stat)
              ElseIf (Trim(process)=='electrodeposition') Then
                Call fitting_species_electrodeposition(model_data%input, model_data%scan, model_data%npoints, &
                                  & i, j, model_data%distance_cutoff, model_data%rotate_species%stat)
              End If
           End If
        End If
      End Do
      j=j+1
    End Do

  End Subroutine insert_species 

  Subroutine fitting_species_electrodeposition(T, scan, npoints, isp, junit, distance_cutoff, rotate)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to fit species over the input electrode surface
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(sample_type),  Intent(InOut) :: T
    Integer(Kind=wi),   Intent(InOut) :: scan(3)
    Integer(Kind=wi),   Intent(InOut) :: npoints(3)
    Integer(Kind=wi),   Intent(In   ) :: isp, junit
    Type(in_param),     Intent(In   ) :: distance_cutoff
    Logical,            Intent(In   ) :: rotate 

    Character(Len=256)   :: messages(5)
    Integer(Kind=wi), Dimension(3) :: ic
    Integer(Kind=wi) :: i, j, l
    Real(Kind=wp) :: centre(3), a(3), b(3)
    Real(Kind=wp) :: r0_rot(max_at_species,3)
    Real(Kind=wp) :: dist
    Logical       :: nofit, too_short

    nofit=.True.
    ic=scan

    Do While (ic(3) < npoints(3) .And. nofit)
      Do While (ic(2) < npoints(2) .And. nofit)
        Do While (ic(1) < npoints(1) .And. nofit)
          Do l=1,3
            centre(l)=(ic(l)-1.0_wp)/(npoints(l))
          End Do
          centre=MatMul(centre,T%cell)
          If (T%species(isp)%topology=='molecule') Then
            If (rotate) Then  
              Call rotate_species(T%species(isp)%definition%r0, T%species(isp)%atoms_per_species, r0_rot)
            End If
          Else
           r0_rot(1,:)=T%species(isp)%definition%r0(1,:)
          End If

          Do l=1,3
           Do i=1, T%species(isp)%atoms_per_species
             r0_rot(i,l)=r0_rot(i,l)+centre(l)
           End Do
          End Do

          i=1
          too_short=.False.
          Do While (i <= T%species(isp)%atoms_per_species .And. .Not.(too_short))
            j=1
            Do While (j <= T%num_atoms_extra .And. .Not.(too_short))
              If (.Not.(T%atom(j)%vanish)) Then
              a(:)=r0_rot(i,:)
              b(:)=T%atom(j)%r(:)
              Call compute_distance_PBC(a, b, T%cell, T%invcell, dist)
              If (dist<distance_cutoff%value) Then
                too_short=.True.
              End If
              End If
              j=j+1
            End Do
            i=i+1
          End Do

          If (.Not. too_short) Then
            scan=ic
            ! Assign arrays
            nofit=.False.
            Do i=1,T%species(isp)%atoms_per_species
              T%atom(T%num_atoms_extra+i)%r(:)=r0_rot(i,:)
              T%atom(T%num_atoms_extra+i)%element=T%species(isp)%definition%element(i)
              T%atom(T%num_atoms_extra+i)%atomic_number=T%species(isp)%definition%atomic_number(i)
              T%atom(T%num_atoms_extra+i)%tag=T%species(isp)%definition%tag(i)
              T%atom(T%num_atoms_extra+i)%vanish=.False.
              T%atom(T%num_atoms_extra+i)%in_species=.True.
              T%atom(T%num_atoms_extra+i)%dynamics=.True.        
            End Do
            Do i=1, T%species(isp)%atoms_per_species
              T%species(isp)%units(T%species(isp)%num_extra+1)%list(i)=T%num_atoms_extra+i
            End Do
            T%num_atoms_extra=T%num_atoms_extra+T%species(isp)%atoms_per_species
            T%species(isp)%num_extra=T%species(isp)%num_extra+1
          End If
          ic(1)=ic(1)+1
        End Do
        ic(1)=1
        ic(2)=ic(2)+1
      End Do
        ic(2)=1
        ic(3)=ic(3)+1
    End Do

    If (nofit) Then
      Write (messages(1),'(2(a,i3),3a)') '***ERROR: fail to deposit unit ', junit, ' (out of ',&
                               T%species(isp)%D_num, ') for species "', Trim(T%species(isp)%tag), '"'
      Write (messages(2),'(a)') 'This is probably due to:'
      Write (messages(3),'(a)') '1) an incorrect geometry for the input model. Either the surface area is not large enough to&
                              & accommodate the required species or there is not enough vacuum region.'
      If (distance_cutoff%fread) Then
        Write (messages(4),'(a,f6.2,a)') '2) a rather large value for the cutoff distance between species, set to ',&
                                       & distance_cutoff%value, ' Angstrom.&
                                       & Please reduce the value of directive "distance_cutoff".'
      Else
        Write (messages(4),'(a,f4.2,a)') '2) the default cutoff distance between species, set to ',&
                                       & distance_cutoff%value, ' is not enough to fit the species. If the is the case, &
                                       & please increase the vaccum region.'
      End If
      If (.Not. rotate) Then
        Write (messages(5),'(a)') '3) for large molecular species, there could be an inconsistency between the orientation/size&
                                & of the species to deposit and the surface of the modelled substrate. Please check.'
      Else
        Write (messages(5),'(a)') '3) for large molecular species, it is convenient to set "rotate_species" to .False. and adapt&
                                & the orientation of the molecule to the substrate.'
      End If
      Call info(messages,5)
      Call execute_command_line('rm  RECORD_MODELS')
      Call error_stop(' ')
    End If

  End Subroutine fitting_species_electrodeposition

  Subroutine fitting_species_intercalation(T, npoints, isp, junit, distance_cutoff, rotate)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to fit species within the input struture 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(sample_type),  Intent(InOut) :: T
    Integer(Kind=wi),   Intent(InOut) :: npoints(3)
    Integer(Kind=wi),   Intent(In   ) :: isp, junit
    Type(in_param),     Intent(In   ) :: distance_cutoff
    Logical,            Intent(In   ) :: rotate 
    
    Integer(Kind=wi), Dimension(3) :: ic
    Integer(Kind=wi) :: i, j, l
    Real(Kind=wp) :: centre(3), centre2(3), a(3), b(3)
    Real(Kind=wp) :: r0_rot(max_at_species,3)
    Real(Kind=wp) :: dist, rn
    Integer(Kind=li) :: nlimit, ncount 
    Logical       :: nofit, too_short
    Character(Len=256)   :: messages(6)

    nofit=.True.
    nlimit=1
    Do i=1,3 
      nlimit=nlimit*npoints(i)
    End Do
    nlimit=nlimit*times_number_attempts
    ncount=0

    Do While (nofit .And. ncount<=nlimit)
      ncount=ncount+1
      Do l=1, 3
        Call random_number(rn)
        ic(l)=floor(npoints(l)*rn+1)
      End Do
      Do l=1,3
        centre2(l)=(ic(l)-1.0_wp)/(npoints(l))
      End Do
      centre=MatMul(centre2, T%cell)
      If (T%species(isp)%topology=='molecule') Then
        If (rotate) Then
          Call rotate_species(T%species(isp)%definition%r0, T%species(isp)%atoms_per_species, r0_rot)
        Else
          r0_rot=T%species(isp)%definition%r0
        End If
      Else
       r0_rot(1,:)=T%species(isp)%definition%r0(1,:)
      End If
   
      Do i=1, T%species(isp)%atoms_per_species
        Do l=1,3
          r0_rot(i,l)=r0_rot(i,l)+centre(l)
        End Do
      End Do
  
      i=1
      too_short=.False. 
      Do While (i <= T%species(isp)%atoms_per_species .And. .Not.(too_short))
        j=1
        Do While (j <= T%num_atoms_extra .And. .Not.(too_short))
          If (.Not.(T%atom(j)%vanish)) Then
          a(:)=r0_rot(i,:)
          b(:)=T%atom(j)%r(:)
          Call compute_distance_PBC(a, b, T%cell, T%invcell, dist)
          If (dist<distance_cutoff%value) Then
            too_short=.True.
          End If
          End If
          j=j+1
        End Do
        i=i+1
      End Do
      
      If (.Not. too_short) Then
        ! Assign arrays  
        nofit=.False.
        Do i=1,T%species(isp)%atoms_per_species
          If (T%num_atoms_extra+i > T%max_atoms) Then 
            Write(messages(1),'(a)')    '***ERROR: it seems the user wants to build a model with a large number of atoms.' 
            Write(messages(2),'(a,i5)') '   The maximum number of atoms to build a model using the input structure is set to ', &
                                     & T%max_atoms
            Write(messages(3),'(a)')    '   However, the required number of atoms as specified in the SET_EQCM settings is larger&
                                     & than this limit.'
            Write(messages(4),'(a)')    '   The user must define a larger value for "multiple_input_atoms" (10000 by default)' 
            Call info(messages, 4) 
            Call error_stop(' ') 
          End If 
          T%atom(T%num_atoms_extra+i)%r(:)=r0_rot(i,:)
          T%atom(T%num_atoms_extra+i)%element=T%species(isp)%definition%element(i)
          T%atom(T%num_atoms_extra+i)%atomic_number=T%species(isp)%definition%atomic_number(i)
          T%atom(T%num_atoms_extra+i)%tag=T%species(isp)%definition%tag(i)
          T%atom(T%num_atoms_extra+i)%vanish=.False.
          T%atom(T%num_atoms_extra+i)%in_species=.True.
          T%atom(T%num_atoms_extra+i)%dynamics=.True.
        End Do    
        Do i=1, T%species(isp)%atoms_per_species
          T%species(isp)%units(T%species(isp)%num_extra+1)%list(i)=T%num_atoms_extra+i
        End Do
        T%num_atoms_extra=T%num_atoms_extra+T%species(isp)%atoms_per_species
        T%species(isp)%num_extra=T%species(isp)%num_extra+1
      End If
      ic(1)=ic(1)+1    
    End Do

    If (ncount>nlimit) Then
      Write (messages(1),'(2(a,i6),3a)') '***ERROR: fail to insert unit ', junit, ' (out of ',&
                               T%species(isp)%D_num, ') for species "', Trim(T%species(isp)%tag), '"'
      Write (messages(2),'(a)') 'The maximum limit of attemps to insert the species has been reached. Possible reasons:' 
      Write (messages(3),'(a)') '1) an incorrect geometry for the input structure.' 
      Write (messages(4),'(a)') '2) the input structure is not large enough to accommodate the required species. The user might&
                               & consider enlarging the input structure by settig directive "scale_cell" to larger than 1.' 
      If (distance_cutoff%fread) Then
        Write (messages(5),'(a,f5.2,a)') '3) a rather large value for the cutoff distance between species, set to ',&
                                       & distance_cutoff%value, ' Angstrom. The user could try reducing the value of directive&
                                       & "distance_cutoff".'
      Else 
        Write (messages(5),'(a,f4.2,a)') '3) the default minimum cutoff distance between species, set to ',&
                                       & distance_cutoff%value, ' Angstrom is not enough. If this is the case, set &
                                       & "scale_cell" to larger than 1 and try again.'
      End If
      If (.Not. rotate) Then
        Write (messages(6),'(a)') '4) for large molecular species, there could be an inconsistency between the orientation/size&
                                 & of the species to intercalate and the structure of the input model. Please check.'
      Else
        Write (messages(6),'(a)') '4) for large molecular species, it is convenient to set "rotate_species" to .False. and adapt&
                                & the orientation of the molecule to the input structure.'
      End If
      Call info(messages,6)
      Call execute_command_line('rm  RECORD_MODELS')
      Call error_stop(' ')
    End If

  End Subroutine fitting_species_intercalation

  Subroutine rotate_species(r0, nat, r0_rot)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to rotate molecular species 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp),    Intent(In   ) :: r0(max_at_species,3)
    Integer(Kind=wi), Intent(In   ) :: nat
    Real(Kind=wp),    Intent(  Out) :: r0_rot(max_at_species,3)

    Integer(Kind=wi) :: i
    Real(Kind=wp) :: rn, alpha, beta, gamma
    Real(Kind=wp) :: rot_matrix(3,3), a(3), b(3)

    !Define the Euler angles
    Call random_number(rn)
    alpha=twopi*rn
    Call random_number(rn)
    beta =twopi*rn
    Call random_number(rn)
    gamma=twopi*rn
    ! Build rotation matrix
    rot_matrix(1,1)=cos(alpha)*cos(beta)
    rot_matrix(2,1)=sin(alpha)*cos(beta)
    rot_matrix(3,1)=-sin(beta)
    rot_matrix(1,2)=cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma)
    rot_matrix(2,2)=sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma)
    rot_matrix(3,2)=cos(beta)*sin(gamma)
    rot_matrix(1,3)=cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma)
    rot_matrix(2,3)=sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma)
    rot_matrix(3,3)=cos(beta)*cos(gamma)

    Do i= 1, nat      
      a(:)=r0(i,:)
      b=MatMul(a, rot_matrix)
      r0_rot(i,:)=b(:)
    End Do

  End Subroutine rotate_species 


  Subroutine read_molecular_units(T, bond_cutoff)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read species units from file 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(species_model),     Intent(InOut) :: T 
    Real(Kind=wp),           Intent(In   ) :: bond_cutoff 
 
    Logical  :: safe, error, loop
    Integer(Kind=wi)  :: i, j, m
    Integer(Kind=wi)  :: iunit, io, num_atoms
    Character(Len=256)   :: set_error, title
    Character(Len=256)   :: messages(4), error_format(10)
    Character(Len= 64)   :: filename 
   
    Integer(Kind=wi)  :: ncount(T%num_components) 

    Integer(Kind=wi) :: num_nn, num_tot
    Integer(Kind=wi), Dimension(max_nn) :: list_nn, list_nn_next
    Integer(Kind=wi) :: ka, k, nacc, knni, il, knn, nlist
    Real(Kind=wp), Dimension(3) :: a, b, dr, geom_centre
    Real(Kind=wp) :: dist
    Logical       :: loop_sp, is_nn 

    error=.False.
    filename=Trim(FOLDER_INPUT_GEOM)//'/'//Trim(T%tag)//'.xyz' 
    set_error='***ERROR in file '//Trim(filename)//':'

    error_format(1) ='In ALC_EQCM, the model to define the units for each species must be:'
    error_format(2) = ' '
    error_format(3) = 'Number of atoms per species (Nsp)'
    error_format(4) = 'Comment with the description of the system (compulsory)'
    error_format(5) = 'Element_1         X_1      Y_1      Z_1       Tag_1'
    error_format(6) = 'Element_2         X_2      Y_2      Z_2       Tag_2'
    error_format(7) = '...........       .....    .....    .....     ...... '
    error_format(8) = 'Element_Nsp       X_Nsp    Y_Nsp    Z_Nsp     Tag_Nsp'    
    error_format(9) = ' '
    error_format(10) = 'where the Tag correspond to the tag used for the component of each species,&
                     & as defined in &Block_Species_Components'

    ! Open the SET_EQCM file with EQCM settings
    Inquire(File=filename, Exist=safe)

    ! Check if file exists
    If (.not.safe) Then
      Write (messages(1),'(4(1x,a))') '***ERROR - File', Trim(filename), 'not found for the&
                                 & specification of species unit', Trim(T%tag)
      Call info(messages,1)
      Call error_stop(' ')
    Else
      Open(Newunit=iunit, File=filename,Status='old')
    End If


    Read (iunit, Fmt=*, iostat=io) num_atoms
    If (io/=0) Then
      Write (messages(1),'(4a)') Trim(set_error), ' First line must be the number of atoms (integer) that constitute species "',&
                                 Trim(T%tag), '". Please check'
      Call info(messages,1)
      error=.True.
    Else
      If (num_atoms/=T%atoms_per_species) Then
        Write (messages(1),'(4a,i3)') Trim(set_error), ' First line must be the number of atoms that constitute species "',&
                               & Trim(T%tag), '", and equal to ', T%atoms_per_species
        Call info(messages,1)
        error=.True.
      End If
    End If

    If (error) Then
      Call info(error_format,10)
      Call error_stop(' ')
    End If

    ! Read title
    Read (iunit, Fmt=*, iostat=io) title


    ! Read specification for the species unit 
    i=0
    Do While (i < T%atoms_per_species)
        i=i+1
        Read (iunit, Fmt=*, iostat=io) T%definition%element(i), (T%definition%r0(i,m), m=1,3), T%definition%tag(i)
        If (io/=0) Then
          If (i== T%atoms_per_species) Then
            Write (messages(1),'(2a,i5)') Trim(set_error), ' Missing specification for atom ', i
          Else
            Write (messages(1),'(2a,i5)') Trim(set_error), ' Wrong specification for atom', i
          End If
          Call info(messages, 1)
          Call info(error_format,10)
          Call error_stop(' ')
        End If
    End Do

    Close(iunit)

    ! Count the number of atoms por component
    ncount=0
    Do i=1, T%atoms_per_species
      loop=.True.
      j=1
      Do While (j <= T%num_components .And. loop)
        If (Trim(T%definition%tag(i))==Trim(T%component%tag(j))) Then
          T%definition%atomic_number(i)=T%component%atomic_number(j)
          ncount(j)=ncount(j)+1
          loop=.False.
          If (Trim(T%definition%element(i))/=Trim(T%component%element(j))) Then
            Write (messages(1),'(4a,i3,3a,i3,3a)') Trim(set_error), ' Setting "', Trim(T%definition%element(i)), &
                                        & '" for the chemical element of atom ', i, ' does not agree with&
                                        & setting "', Trim(T%component%element(j)), '" specified for component ', &
                                        & j, ' of species "', Trim(T%tag), '" in  &Block_Species_Components'
            Call info(messages,1)
            Call error_stop(' ')  
          End If
        End If
        j=j+1
      End Do
      If (loop) Then
        Write (messages(1),'(4a,i3,3a)') Trim(set_error), ' Setting "', Trim(T%definition%tag(i)), &
                                    & '" for the tag of atom ', i, ' has not been assigned any component tag&
                                    & of species "', Trim(T%tag), '" in  &Block_Species_Components'
        Call info(messages,1)
        Call error_stop(' ')  
      End If
    End Do
    
    Do j=1, T%num_components
        If (ncount(j)/=T%component%N0(j)) Then
           Write (messages(1),'(4a,i3,a,i3,a)') Trim(set_error), ' The number of components with tag "', Trim(T%component%tag(j)),&
                                         & '" is equal to ', ncount(j), ', which is different from the amout of ', &
                                         & T%component%N0(j), ' specified in &Block_Species_Components. Please check the&
                                         & definition for the molecular geometry.'    
          Call info(messages,1)
          Call error_stop(' ')  
        End If
    End Do

    ! Check it is a sensible geometry
    T%definition%in_species=.False.
    num_nn=1
    list_nn(1)=1
    ka=1
    loop=.True. 
    num_tot=0
    list_nn_next=0 
    nlist=1


    Do While (loop)
      nacc=0
      Do knn=1, num_nn
        ka=list_nn(knn)
        Do k=1, T%atoms_per_species
          loop_sp=.True.
          If (k/=ka .And. (.Not. T%definition%in_species(k)) .And. (.Not. T%definition%in_species(ka))) Then
             ! Compute distances
             a=T%definition%r0(k,:)
             b=T%definition%r0(ka,:)
             dr=a-b
             dist=Norm2(dr)
             ! Check if distances are within a given cutoff
             If (dist<bond_cutoff) Then 
               If (dist<min_intra) Then
                 Write (messages(1),'(2a)') Trim(set_error), ' wrong input geometry'
                 Write (messages(2),'(2(a,i4))') 'Intermolecular distance between atoms ', k, ' and ', ka
                 Write (messages(3),'(a,f4.2,a)') 'is shorter than the input minimum distance criteria of ', min_intra,&
                                                 & ' Angstrom for bonding. Please review the input geometry.'  
                 Call info(messages,3)
                 Call error_stop(' ')
               Else
                 loop_sp=.False.
                 is_nn=.True.
                 il=1
                 Do While (il <= max_nn .And. is_nn) 
                   If (k==list_nn_next(il)) Then
                     is_nn=.False.
                   End If
                   il=il+1
                 End Do 
                 If (is_nn) Then    
                   nacc=nacc+1
                   list_nn_next(nacc)=k
                 End If
               End If  
             End If
          End If
        End Do
        If (.Not. T%definition%in_species(ka)) Then
          T%definition%in_species(ka)=.True.
          num_tot=num_tot+1
          If (num_tot+nacc == T%atoms_per_species) Then
            Do knni= 1, nacc
              T%definition%in_species(list_nn_next(knni))=.True.
            End Do
            loop=.False.
            num_tot=T%atoms_per_species
          End If
        End If
      End Do
      num_nn=nacc  
      list_nn=list_nn_next
      list_nn_next=0 
      If (nacc==0) loop=.False.
    End Do

    If (T%atoms_per_species/=num_tot) Then
      Write (messages(1),'(2a)') Trim(set_error), ' Problems with the input geometry'
      Write (messages(2),'(a)')  'Possible reasons:' 
      Write (messages(3),'(a)')          ' 1) wrong input geometry, where species are dissociated (e.g. H3O instead of H2O).&
                                        & If this is the case, the user must provide an input model consistent with the&
                                        & definition of species.' 
      Write (messages(4),'(a,f4.2,3a)')  ' 2) the value for the bond cutoff (', bond_cutoff, &
                                        & ' Angstrom) as specified for species ',  Trim(T%tag) ,' in  &block_species_components&
                                        & is not sufficient. The user must adjust this value' 
      Call info(messages, 4)
      Call error_stop(' ')
    End If

    ! Compute the geometric centre of the species and refer atomic coordinates to such point in space 
    geom_centre=0.0_wp
    Do i =1, T%atoms_per_species
      Do k=1,3
        geom_centre(k)=geom_centre(k)+T%definition%r0(i,k)
      End Do
    End Do

    geom_centre=geom_centre/T%atoms_per_species

    Do i =1, T%atoms_per_species
      Do k=1,3
        T%definition%r0(i,k)=T%definition%r0(i,k)-geom_centre(k)
      End Do
    End Do  
    
  End Subroutine read_molecular_units


  Subroutine remove_species(T, types_species, process, stage)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to remove atoms from the model  
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(sample_type), Intent(InOut) :: T 
    Integer(Kind=wi),  Intent(in   ) :: types_species
    Character(Len=*),  Intent(In   ) :: process
    Character(Len=*),  Intent(In   ) :: stage 

    Integer(Kind=wi) :: i, j, k, m, l, na

    Integer(Kind=wi) :: list_rm(max_num_species_units)
    Real(Kind=wp) :: rn
    Logical       :: hit(max_at_species)

    list_rm=0
    Do i = 1, types_species
      If (T%species(i)%change_content) Then
        If (T%species(i)%D_num<0) Then
          If (Trim(process)=='intercalation' .Or. &
            (Trim(process)=='electrodeposition' .And. Trim(stage)=='sample')) Then
            hit=.False.
            k=1 
            Do While (k <= (-T%species(i)%D_num)) 
              Call random_number(rn)
              na=floor(T%species(i)%num*rn)+1
              If (.Not. T%species(i)%units(na)%vanish) Then
                list_rm(k)=na 
                T%species(i)%units(na)%vanish=.True.
                k=k+1
              End If
            End Do
          ElseIf (Trim(process)=='electrodeposition' .And. Trim(stage)=='input') Then
            k=1
            Do j=T%species(i)%num, T%species(i)%num+T%species(i)%D_num+1, -1
              list_rm(k)=j
              T%species(i)%units(j)%vanish=.True.
              k=k+1
            End Do
          End If
          ! remove 
          Do k =1, -T%species(i)%D_num
             j=list_rm(k)
             Do m= 1, T%species(i)%atoms_per_species
               l=T%species(i)%units(j)%list(m)
               T%atom(l)%vanish=.True.
             End Do  
          End Do

        End IF
      End If
    End Do

  End Subroutine remove_species

  Subroutine read_input_model(files, stoich_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read input model from file INPUT_GEOM/INPUT_STRUCTURE
    ! Format Options: 
    ! - VASP
    ! - xyz 
    ! - ONETEP
    ! - CASTEP  
    ! - CIF
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(stoich_type), Intent(In   ) :: stoich_data
    Type(model_type),  Intent(InOut) :: model_data

    Logical            :: safe
    Integer(Kind=wi)   :: i, ifolder
    Character(Len=256) :: message, exec_cif
    Character(Len=256) :: messages(3)
    Character(Len=32 ) :: input_file, set_error, path_to_file

    input_file=Trim(files(FILE_INPUT_STRUCTURE)%filename)

    ! Check if folder INPUT_GEOM exists
    Call execute_command_line('[ -d '//Trim(FOLDER_INPUT_GEOM)//' ]', exitstat=ifolder)
    If (ifolder/=0) Then
      Call info(' ', 1)
      Write (messages(1), '(1x,3a)') '***ERROR: folder ', Trim(FOLDER_INPUT_GEOM), ' cannot be found.'
      Write (messages(2), '(1x,3a)') 'This folder must contain file ', Trim(input_file),' and the xyz&
                                   & files only for the participant MOLECULAR species (if any).'
      Write (messages(3), '(1x,a)') 'The requested analysis cannot be conducted. Please create the folder&
                                   & and add the required information.'
      Call info(messages, 3)
      Call error_stop(' ')
    End If
 
    path_to_file = Trim(FOLDER_INPUT_GEOM)//'/'//Trim(files(FILE_INPUT_STRUCTURE)%filename)
    set_error = '***ERROR -'
    model_data%change_format=.False.

    Inquire(File=path_to_file, Exist=safe)

    If (.not.safe) Then
      Call info(' ', 1)
      Write (message,'(4(1x,a))') Trim(set_error), 'File', Trim(path_to_file), ' not found'
      Call error_stop(message)
    End If

    If (Trim(model_data%input_model_format%type) == 'cif') Then
      Inquire(File=Trim(FOLDER_INPUT_GEOM)//'/'//Trim(files(FILE_INPUT_STRUCTURE)%filename)//'.cif', Exist=safe)
      If (safe) Then
        exec_cif='mv '//Trim(path_to_file)//'.cif  '//Trim(path_to_file)
        Call execute_command_line(exec_cif, Exitstat=i) 
      End If
      exec_cif='cp '//Trim(path_to_file)//' input.cif'
      Call execute_command_line(exec_cif, Exitstat=i) 
      exec_cif='mv '//Trim(path_to_file)//' '//Trim(path_to_file)//'.cif'
      Call execute_command_line(exec_cif, Exitstat=i) 
      exec_cif= 'python cif-input-geom.py'
      Call execute_command_line(exec_cif, Exitstat=i)
      If (i/=0) Then
        Call info(' ', 1)
        Write (message,'(1x,4a)') Trim(set_error), 'Either ', Trim(path_to_file), ' file does not comply with the "cif"&
                                   & format or there is an error in the file. Otherwise, is Python and ASE installed in&
                                   & your machine?'
        exec_cif='rm input.cif '//Trim(input_file)
        Call execute_command_line(exec_cif)
        exec_cif='mv '//Trim(path_to_file)//'.cif '//Trim(path_to_file)
        Call execute_command_line(exec_cif, Exitstat=i) 
        Call error_stop(message)      
      End If     
      model_data%input_model_format%type = 'vasp'
      model_data%change_format=.True.
      model_data%old_format='cif'
      exec_cif='mv '//Trim(input_file)//' '//Trim(path_to_file)
      Call execute_command_line(exec_cif, Exitstat=i)
      exec_cif='rm input.cif '
      Call execute_command_line(exec_cif)
    End If

    ! Open the INPUT_STRUCTURE file
    Open(Newunit=files(FILE_INPUT_STRUCTURE)%unit_no, File=Trim(path_to_file),Status='old')

    ! Select the calling of the subroutine according to the specification of input_model_format directive
    If (Trim(model_data%input_model_format%type) == 'vasp') Then
        Call read_input_vasp_format(files, stoich_data, model_data)
    Else If (Trim(model_data%input_model_format%type) == 'xyz') Then
            Call read_input_xyz_format(files, model_data)
    Else If (Trim(model_data%input_model_format%type) == 'castep' .Or. & 
            Trim(model_data%input_model_format%type) == 'onetep') Then
            Call read_input_geom_format(files, model_data)  
    End If

    Call about_cell(model_data%input%cell,model_data%input%invcell,model_data%input%cell_length)

    ! Close file
    Close(files(FILE_INPUT_STRUCTURE)%unit_no) 

    If (model_data%change_format) Then
      If (Trim(model_data%old_format)=='cif') Then
        exec_cif='mv '//Trim(path_to_file)//'.cif '//Trim(path_to_file)
        Call execute_command_line(exec_cif)
      End If
    End If

 
  End Subroutine read_input_model

  Subroutine shift_structure(model_data, process)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to shift the set of atomic coordinates to the origin 
    !
    ! author    - i.scivetti March 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data
    Character(Len=*),    Intent(In   ) :: process

    Integer(Kind=wi) :: i, indx
    Real(Kind=wp)    :: min_dist, max_dist
    Real(Kind=wp)    :: dist
    Real(Kind=wp)    :: r(3), shift(3), normal(3), level(3)

    Character(Len=256) :: message

    min_dist=Huge(1.0_wp)
    normal=0.0_wp

    If (Trim(process)=='intercalation') Then
      ! Find the closest atom to the (0,0,0) 
      Do i = 1, model_data%input%num_atoms
        r(:)=model_data%input%atom(i)%r(:)
        dist=norm2(r)
        If (dist < min_dist) Then
          min_dist = dist
          shift=r
        End If 
      End Do
  
      ! Shift all the atoms 
      Do i = 1, model_data%input%num_atoms
        model_data%input%atom(i)%r(:)=model_data%input%atom(i)%r(:)-shift(:)    
      End Do
 
    Else If (Trim(process)=='electrodeposition') Then
      If (Trim(model_data%normal_vector%type)=='c3') Then
        indx=3
      ElseIf (Trim(model_data%normal_vector%type)=='c2') Then
        indx=2
      ElseIf (Trim(model_data%normal_vector%type)=='c1') Then
        indx=1
      End If

      r(:)=model_data%input%cell(indx,:)
      Call normal_along_vector(r,normal)

      ! Find the atom with the lowest component along the normal vector to the surface 
      Do i = 1, model_data%input%num_atoms
        r(:)=model_data%input%atom(i)%r(:)
        dist=Dot_product(r,normal)
        If (dist < min_dist) Then
          min_dist = dist
          shift=r
        End If 
      End Do

      ! Shift atoms
      Do i = 1, model_data%input%num_atoms
        model_data%input%atom(i)%r(:)=model_data%input%atom(i)%r(:)-shift(:)
      End Do

      If (model_data%deposition_level%fread) Then
        level=0.0_wp
        model_data%deposition_level%value=model_data%deposition_level%value-shift(indx)
        level(indx)=model_data%deposition_level%value
        r(:)=model_data%input%cell(indx,:)
        If ( Abs(model_data%deposition_level%value) >= norm2(r) .Or. Dot_product(level,normal) < 0.0_wp ) Then
           Write (message, '(1x,a)') '***ERROR: innapropriate value of "deposition_level" for the input system.&
                                    & Please check model and simulation box.'
           Call error_stop(message)      
        End If 
        Write (message, '(1x,a)') 'IMPORTAT: The user has defined the point along the dimension perpendicular to the surface from&
                                & which species will be added. Please check this setting is reasonable for the model.'
        Call info(message, 1) 
      Else
        ! Find the atom with the largest component along the normal vector to the surface 
        max_dist=-Huge(1.0_wp)
        Do i = 1, model_data%input%num_atoms
          r(:)=model_data%input%atom(i)%r(:)
          dist=Dot_product(r,normal)
          If (dist > max_dist) Then
            max_dist = dist
          End If 
        End Do
        model_data%deposition_level%value=max_dist 
      End If

    End If
    
  End Subroutine shift_structure

  Subroutine normal_along_vector(a,normal)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the normal corresponding to a lattice vector 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   )   :: a(3)
    Real(Kind=wp), Intent(  Out)   :: normal(3)
 
    Real(Kind=wp) :: norm
    Integer(Kind=wi) :: i

    norm=0.0_wp

    Do i=1, 3
      norm=norm+a(i)**2
    End Do

    norm=sqrt(norm)
    normal=a/norm

  End Subroutine normal_along_vector 

  Subroutine check_cell_consistency(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the consistency between the simulation cell and the
    ! atomic coordinates  
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data

    Integer(Kind=wi) :: i, j
    Real(Kind=wp)    :: dist, dist0
    Logical          :: changed_geo, error 

    Real(Kind=wp)       :: a(3), b(3), r(3)
    Character(Len=256)  :: message, messages(2)
    Real(Kind=wp)       :: r0(model_data%input%num_atoms,3)
    
    Do i=1, model_data%input%num_atoms
      r0(i,:)=model_data%input%atom(i)%r(:)
    End Do

    error=.False.

    Write (messages(1),'(a)') '***ERROR: Inconsistency found between the atomic coordinates&
                            & and the simulation cell vectors. Please verify input coordinates and cell'

    ! Move all inside the ce
    Do i = 1, model_data%input%num_atoms
      r(:)=model_data%input%atom(i)%r(:)
      Call check_PBC(r, model_data%input%cell, model_data%input%invcell, 1.0_wp, changed_geo)
      If (changed_geo) Then
        model_data%input%atom(i)%r(:)=r(:)
      End If
    End Do

    ! Try to move atoms again....if any atom is now moved, there is an inconsistency
    ! between the coordinates and the cell 
    Do i =  1, model_data%input%num_atoms
      r(:)=model_data%input%atom(i)%r(:)
      Call check_PBC(r, model_data%input%cell, model_data%input%invcell, 1.0_wp, changed_geo)
      If (changed_geo) Then
        ! Once again. If any atom is now moved, there is an inconsistency
        Call check_PBC(r, model_data%input%cell, model_data%input%invcell, 1.0_wp, changed_geo)
      End If  
      If (changed_geo) Then
        Write (message,'(a,i4)') '***PROBLEMS with atom',  i
        Call info(message, 1)
        error=.True.
      End If
    End Do

    If (error) Then
       If (.Not. model_data%shift_structure%stat) Then                 
         Write (messages(2),'(a)')    '*** The user should try by setting the "shift_structure"&
                                     & directive to .True.'
         Call info(messages,2)
       Else
         Call info(messages,1)
       End If
       Call error_stop(' ')
    End If

    Do i = 1, model_data%input%num_atoms-1
      Do j = i+1, model_data%input%num_atoms
        a(:)=model_data%input%atom(i)%r(:)
        b(:)=model_data%input%atom(j)%r(:)
        Call compute_distance_PBC(a, b, model_data%input%cell, model_data%input%invcell, dist)
        a(:)=r0(i,:)
        b(:)=r0(j,:)
        Call compute_distance_PBC(a, b, model_data%input%cell, model_data%input%invcell, dist0)
        If (Abs(dist-dist0)>length_tol) Then
          Write (message,'(a,2(i7,a))') '***PROBLEMS: Distance between atom ', i, ' and ', j, &
                                   &' does not comply with the crystal symmetry imposed by the cell vectors.'
          Call info(message,1)
          Call info(messages,1)
          Call error_stop(' ')
        End If
      End Do
    End Do
 
  End Subroutine check_cell_consistency


  Subroutine read_input_geom_format(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read file INPUT_GEOM/INPUT_STRUCTURE according to either 
    ! CASTEP or ONETEP format 
    !
    ! author    - i.scivetti March 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: io, num_geo, num, iunit
    Integer(Kind=wi) :: internal
    Integer(Kind=wi) :: i, j, k, m
    Character(Len=256)    :: messages(2)
    
    Character(Len=256)    :: set_error, input_file, path_to_file
    Character(Len=256)    :: exec_numgeo, word

    input_file = Trim(files(FILE_INPUT_STRUCTURE)%filename)
    path_to_file=Trim(FOLDER_INPUT_GEOM)//'/'//Trim(input_file)
    iunit=files(FILE_INPUT_STRUCTURE)%unit_no

    set_error = '***ERROR in file '//Trim(path_to_file)//&
               &' (inconsistency with CASTEP/ONETEP format)'

    ! Check the correctness of the file

    exec_numgeo='grep "<-- E" '//Trim(path_to_file)//' > number_geom.dat'
    Call execute_command_line(exec_numgeo) 
    Call execute_command_line('wc -l number_geom.dat > nlines.dat')
    Open(Newunit=internal, File='nlines.dat' ,Status='old')
    Read (internal, Fmt=*, iostat=io) num_geo
    If (num_geo == 0) Then
      Call error_stop(set_error)
    Else If (num_geo > 1) Then
      Write (messages(1),'(2a,i3,a)') Trim(set_error), ': number of geometrical configurations indentified&
                                   & is equal to', num_geo, '. Please provide only one configuration.'
      Call error_stop(messages(1)) 
    End If
    Close(internal)
    Call execute_command_line('rm nlines.dat number_geom.dat')

    exec_numgeo='grep "<-- R" '//Trim(path_to_file)//' > number_atoms.dat'
    Call execute_command_line(exec_numgeo) 
    Call execute_command_line('wc -l number_atoms.dat > nat.dat')
    Open(Newunit=internal, File='nat.dat' ,Status='old')
    Read (internal, Fmt=*, iostat=io) model_data%input%num_atoms
    If (model_data%input%num_atoms == 0) Then
      Call error_stop(set_error)
    End If
    Close(internal)
    Call execute_command_line('rm nat.dat number_atoms.dat')

    If (model_data%input%num_atoms /= model_data%block_component%numtot) Then
     Write (messages(1),'(3a)') 'ERROR***: Inconsistency between the number of atoms in file ',&
                             & Trim(path_to_file),&
                             &' and the amount of atoms defined in &Block_input_composition. Please check'   
     Call error_stop(messages(1))
    End If

    ! Read the file
    Read (iunit, Fmt=*, iostat=io) word  
    Read (iunit, Fmt=*, iostat=io) word  
    Do i= 1, 3  
      Read (iunit, Fmt=*, iostat=io) (model_data%input%cell(i,j), j=1,3) 
      If (io/=0) Then
        Write (messages(1),'(2a,i1,a)') Trim(set_error), ': Problems with the specification of cell vector ', i,&
                                     & '. See those lines that end with "<-- h"'  
        Call error_stop(messages(1))
      End If 
      model_data%input%cell(i,:)=Bohr_to_A*model_data%input%cell(i,:) 
    End Do

    ! Allocate atomic arrays for the input model   
    Call model_data%atomic_arrays_input(model_data%input%num_atoms)  

    ! Read atomic coordinates
    j=1; k=0; i=0
    Do While (i < model_data%input%num_atoms)
      If (model_data%block_component%N0(j)/=0) Then
        i=i+1
        Read (iunit, Fmt=*, iostat=io) model_data%input%atom(i)%element, num, (model_data%input%atom(i)%r(m), m=1,3)
        If (io/=0) Then
          Write (messages(1),'(2a,i5)') Trim(set_error), ' Wrong specification for atom', i
          Call info(messages, 1)
          Call error_stop(' ')
        Else
          If (Trim(model_data%input%atom(i)%element)/=Trim(model_data%block_component%element(j))) Then
            Write (messages(1),'(a,i5,5a)') '***ERROR: Chemical element of atom ', i, &
                                          & ' is set to "', Trim(model_data%input%atom(i)%element), '", but& 
                                          & according to the definition of &block_input_composition it must be "',&
                                          & Trim(model_data%block_component%element(j)), '".'
            Write (messages(2),'(3a)')      'Check the consistency between the data in ', &
                                          &  Trim(path_to_file), ' and &block_input_composition'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End If
        ! Transform to Angstom
        model_data%input%atom(i)%r=Bohr_to_A*model_data%input%atom(i)%r 
        k=k+1
        model_data%input%atom(i)%tag=model_data%block_component%tag(j)
        If (k==model_data%block_component%N0(j)) Then
          k=0
          j=j+1
        End If
      Else
        k=0
        j=j+1
      End If
    End Do

  End Subroutine read_input_geom_format 


  Subroutine read_input_xyz_format(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read file INPUT_STRUCTURE according to xyz format 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: io, iunit
    Integer(Kind=wi) :: i, j, k, m
    Character(Len=256)    :: title
    Character(Len=256)    :: messages(5)
    
    Character(Len=256)    :: set_error, error_format(10), input_file, path_to_file

    input_file = Trim(files(FILE_INPUT_STRUCTURE)%filename) 
    path_to_file = Trim(FOLDER_INPUT_GEOM)//'/'//Trim(input_file)    
    iunit=files(FILE_INPUT_STRUCTURE)%unit_no

    set_error = '***ERROR in file '//Trim(path_to_file)//' (inconsistency with xyz format):'
    error_format(1) ='In ALC_EQCM, the structure of the '//Trim(input_file)//' file in xyz format must be:'
    error_format(2) = ' '
    error_format(3) = 'Number of atoms (Nat)'
    error_format(4) = 'Comment with the description of the system (compulsory)'
    error_format(5) = 'Element_1         X_1      Y_1      Z_1'
    error_format(6) = 'Element_2         X_2      Y_2      Z_2'
    error_format(7) = '...........       .....    .....    .....'
    error_format(8) = 'Element_Nat       X_Nat    Y_Nat    Z_Nat'
    error_format(9) = ' ' 
    error_format(10) = 'Please check consistency between the structure of the file and&
                       & directive "output_model_format".'
    ! Start reading the file
    !!!!!!!!!!!!!!!!!!!!!!!!

    ! Read number of atoms 
    Read (iunit, Fmt=*, iostat=io) model_data%input%num_atoms
    If (io/=0) Then
      Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the number of atoms, which must be&
                               & specified in the first line.'
      Call info(messages,1)
      Call info(error_format,10)
      Call error_stop(' ')
    End If    

    If (model_data%input%num_atoms /= model_data%block_component%numtot) Then
     Write (messages(1),'(3a)') 'ERROR***: Inconsistency between the number of atoms in file ', Trim(path_to_file),&
                             &' and the amount of atoms defined in &Block_input_composition. Please check.'   
     Call error_stop(messages(1))
    End If

    ! Read title
    Read (iunit, Fmt=*, iostat=io) title

    ! Allocate atomic arrays for the input model   
    Call model_data%atomic_arrays_input(model_data%input%num_atoms)  

    ! Read atomic coordinates
    j=1; k=0; i=0
    Do While (i < model_data%input%num_atoms)
      If (model_data%block_component%N0(j)/=0) Then
        i=i+1
        Read (iunit, Fmt=*, iostat=io) model_data%input%atom(i)%element, (model_data%input%atom(i)%r(m), m=1,3)
        If (io/=0) Then
          If (i== model_data%input%num_atoms) Then
            Write (messages(1),'(2a)') Trim(set_error), ' Missing input coordinates somewhere in the list.'
          Else
            Write (messages(1),'(2a,i5)') Trim(set_error), ' Wrong specification for atom', i
          End If
          Call info(messages, 1)
          Call info(error_format,8)
          Call error_stop(' ')
        Else
          If (Trim(model_data%input%atom(i)%element)/=Trim(model_data%block_component%element(j))) Then
            Write (messages(1),'(a,i5,5a)') '***ERROR: Chemical element of atom ', i, &
                                          & ' is set to "', Trim(model_data%input%atom(i)%element), '", but& 
                                          & according to the definition of &block_input_composition it must be "',&
                                          & Trim(model_data%block_component%element(j)), '".'
            Write (messages(2),'(3a)')  'Check the consistency between the data in ', Trim(path_to_file),&
                                          & ' and &block_input_composition.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End If
        k=k+1
        model_data%input%atom(i)%tag=model_data%block_component%tag(j)
        If (k==model_data%block_component%N0(j)) Then
          k=0
          j=j+1
        End If
      Else
        k=0
        j=j+1
      End If
    End Do

  End Subroutine read_input_xyz_format


  Subroutine read_input_vasp_format(files, stoich_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read file INPUT_STRUCTURE according to VASP/POSCAR format 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(stoich_type),   Intent(In   ) :: stoich_data
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: io, iunit
    Integer(Kind=wi) :: i, j, k, m
    Character(Len=256)  :: title, word, exec
    Character(Len=256)  :: messages(5)
    Logical             :: loop, loop2, error_vasp, error_elements
    Character(Len=  2)  :: buffer
    
    Character(Len=256)  :: input_file, path_to_file, set_error, webpage

    Real(Kind=wp)       :: v_cart(3)

    Character(Len= 2), Allocatable :: element_file(:)
    Integer(Kind=wi),  Allocatable :: amount_file(:)
    Integer(Kind=wi)  :: fail(2)

    error_vasp = .False.
    input_file = Trim(files(FILE_INPUT_STRUCTURE)%filename)
    path_to_file= Trim(FOLDER_INPUT_GEOM)//'/'//Trim(input_file)   
    iunit=files(FILE_INPUT_STRUCTURE)%unit_no

    If (model_data%change_format) Then
      If (Trim(model_data%old_format)=='cif') Then
        set_error = '***ERROR in file '//Trim(path_to_file)//' (problems with CIF format):'
      End If
    Else  
      set_error = '***ERROR in file '//Trim(path_to_file)//' (inconsistency with POSCAR format):'
      webpage='For correct format and details see the following link: https://www.vasp.at/wiki/index.php/POSCAR'
    End If
    

    ! Determine vasp elements and amounts of atoms from the data in &block_input_composition and &block_species_components
    model_data%input%list%num_elements=1 
    j=1
    model_data%input%list%N0(j)=model_data%block_component%N0(j)
    model_data%input%list%element(j)=model_data%block_component%element(j)   
    buffer=model_data%block_component%element(j)
    Do i = 2,  stoich_data%total_tags
      If (Trim(model_data%block_component%element(i)) /= Trim(buffer)) Then
        buffer=model_data%block_component%element(i)
        If (model_data%block_component%N0(i)/=0) Then
          model_data%input%list%num_elements=model_data%input%list%num_elements+1
          j=j+1
          model_data%input%list%N0(j)=model_data%block_component%N0(i)
          model_data%input%list%element(j)=model_data%block_component%element(i)
        End If
      Else
        model_data%input%list%N0(j)=model_data%input%list%N0(j)+model_data%block_component%N0(i)
      End If
    End Do 

    model_data%input%num_atoms=0
    Do i=1, model_data%input%list%num_elements
       model_data%input%num_atoms=model_data%input%num_atoms + model_data%input%list%N0(i)
    End Do

    ! Allocate atomic arrays for the input model   
    Call model_data%atomic_arrays_input(model_data%input%num_atoms)  

    ! Allocate vasp related variables for reading the file (deallocated below)
    Allocate(element_file(model_data%input%list%num_elements), Stat=fail(1))
    Allocate(amount_file(model_data%input%list%num_elements),  Stat=fail(2))

    If (Any(fail > 0)) Then
      Write (messages(1),'(1x,1a)') '***ERROR: Allocation problems for arrays in subroutine "read_input_vasp_format"'
      Call error_stop(messages(1))
    End If

    ! Start reading the file
    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Read elements or title
    Read (iunit, Fmt=*, iostat=io) (element_file(i) , i=1, model_data%input%list%num_elements)
    loop2=.True.
    i=1
    Do While (i<= model_data%input%list%num_elements .And. loop2)
      loop=.True.
      k=1
      Do While (k <= NPTE .And. loop)
        If (Trim(chemsymbol(k))==Trim(element_file(i))) Then
          loop=.False.
        End If
        k=k+1
      End Do
      If (loop) Then
        loop2=.False.
      End If
      i=i+1
    End Do

    If (loop2) Then
      error_elements=.False.
    Else
      error_elements=.True.
      Rewind iunit 
      Read (iunit, Fmt=*, iostat=io) title
    End If   
 
    ! Read scale factor
    Read (iunit, Fmt=*, iostat=io) model_data%input%scale_factor_vasp
    If (io/=0) Then
      Write (messages(1),'(2a)') Trim(set_error), ' Scale factor for cell vectors is invalid.'
      Write (messages(2),'(a)') webpage
      Call info(messages,2)
      Call error_stop(' ')
    End If 

    Do i= 1, 3
      Read (iunit, Fmt=*, iostat=io) (model_data%input%cell(i,j), j=1,3)
      If (io/=0) Then
        Write (messages(1),'(2a,i1,a)') Trim(set_error), ' Definition for cell vector ', i, ' is incorrect.'
        Write (messages(2),'(a)') webpage
        Call info(messages,2)
        Call error_stop(' ')
      End If 
    End Do

    If (error_elements) Then
      Read (iunit, Fmt=*, iostat=io) (element_file(i) , i=1, model_data%input%list%num_elements)
      If (io/=0) Then
        Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the list of atomic species.'
        Write (messages(2),'(a)') webpage
        Call info(messages,2)
        Call error_stop(' ') 
      End If
    End If

    Do i=1, model_data%input%list%num_elements
      loop=.True.
      k=1
      Do While (k <= NPTE .And. loop)
        If (Trim(chemsymbol(k))==Trim(element_file(i))) Then
          loop=.False.
        End If
        k=k+1
      End Do
      If (loop) Then
        Write (messages(1),'(4a,i2,a)') Trim(set_error), ' Chemical element "' , Trim(element_file(i)), &
                                       & '" defined for component ',  i, ' of the list does not correspond to an element&
                                       & of the Periodic Table. Please use a valid element.'
        Write (messages(2),'(a)')   webpage
        Write (messages(3),'(3a)') 'IMPORTANT: The user should also check for inconsistencies between the list/number of atoms&
                                & in file ', Trim(input_file),' and the settings in &Block_input_composition.&
                                & Have you missed info?'
        Call info(messages,3)
        Call error_stop(' ')
      End If

      If (Trim(element_file(i)) /= Trim(model_data%input%list%element(i))) Then
        error_vasp=.True.
      End If 
    End Do


    Read (iunit, Fmt=*, iostat=io) (amount_file(i) , i=1, model_data%input%list%num_elements)
    If (io/=0) Then
      Read (iunit, Fmt=*, iostat=io) (amount_file(i) , i=1, model_data%input%list%num_elements)
      If (io/=0) Then
        Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the list with the number of atoms&
                                               & for each atomic species'
        Write (messages(2),'(a)') webpage
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End If

    If (model_data%change_format) Then
      If (Trim(model_data%old_format)=='cif') Then
        Write (messages(1),'(3a)') '***ERROR: Inconsistent between the list of atomic species in ', Trim(path_to_file),&
                             & ' (CIF format) and the data specified in &block_input_composition.'
      End If
    Else
      Write (messages(1),'(3a)') '***ERROR: Inconsistent between the list of atomic species in ',  Trim(path_to_file),&
                             & ' (VASP format) and the data specified in &block_input_composition.'
    End If
    Write (messages(2),'(a)') 'The list of elements and amount of atoms per element resulting from the definition&
                             & in &block_input_composition should be:'
          
    Write (messages(3),'(*(4x,a2))') (model_data%input%list%element(j), j=1, model_data%input%list%num_elements)
    Write (messages(4),'(*(1x,i5))') (model_data%input%list%N0(j), j=1, model_data%input%list%num_elements)
    Write (messages(5),'(3a)') 'Which DO NOT MATCH the order/amount of atoms in file ', Trim(path_to_file), & 
                            & '. Please review the settings of &block_input_composition (or the geometry of the input structure)'

    Do i=1, model_data%input%list%num_elements
      If (amount_file(i) /= model_data%input%list%N0(i)) Then
        error_vasp=.True.
      End If
    End Do

    If (error_vasp) Then
      Call info(messages, 5)
      If (model_data%change_format) Then
        If (Trim(model_data%old_format)=='cif') Then
          exec='mv input.cif '//Trim(files(FILE_INPUT_STRUCTURE)%filename)
        End If
        Call execute_command_line(exec)
      End If
      Call error_stop(' ')
    End If

    Read (iunit, Fmt=*, iostat=io) word
    Call capital_to_lower_case(word)
    If (word(1:4)=='sele') Then
      model_data%selective_dyn=.True. 
    Else       
      Backspace iunit
    End If

    Read (iunit, Fmt=*, iostat=io) model_data%input%list%coord_type
    If (io/=0) Then
      Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the type of coordinates. Valid options are either& 
                                             & "Cartesian" or "Direct"'
      Write (messages(2),'(a)') webpage
      Call info(messages,2)
      Call error_stop(' ')
    End If
    Call capital_to_lower_case(model_data%input%list%coord_type)   
  
    If (Trim(model_data%input%list%coord_type) /= 'cartesian' .And. &
      Trim(model_data%input%list%coord_type) /= 'direct') Then
      Write (messages(1),'(2a)') Trim(set_error), ' Wrong option for the type of coordinates. Valid options are either& 
                                             & "Cartesian" or "Direct"'
      Write (messages(2),'(a)') webpage
      Call info(messages,2)
      Call error_stop(' ')
    End If

    ! Read atomic coordinates
    j=1; k=0; i=0
    Do While (i < model_data%input%num_atoms)
      If (model_data%block_component%N0(j)/=0) Then
        i=i+1
        If (model_data%selective_dyn) Then
          Read (iunit, Fmt=*, iostat=io) (model_data%input%atom(i)%r(m), m=1,3), &
                                              (model_data%input%atom(i)%dynamics(m), m=1,3)
        Else 
          Read (iunit, Fmt=*, iostat=io) (model_data%input%atom(i)%r(m), m=1,3)
        End If
        If (io/=0) Then
          If (i== model_data%input%num_atoms) Then
            Write (messages(1),'(2a)') Trim(set_error), ' Missing input coordinates somewhere in the list'
          Else
            Write (messages(1),'(2a,i5)') Trim(set_error), ' Wrong specification for the input coordinates of atom ', i
          End If
          Call info(messages, 1)
          Call error_stop(' ')
        End If
        k=k+1
        model_data%input%atom(i)%tag=model_data%block_component%tag(j)
        model_data%input%atom(i)%element=model_data%block_component%element(j)
        If (k==model_data%block_component%N0(j)) Then
          k=0
          j=j+1 
        End If
      Else
        k=0
        j=j+1
      End If
    End Do

    ! Deallocate arrays
    Deallocate(amount_file)
    Deallocate(element_file)

    ! Transform using the scaling factor
    If (Trim(model_data%input%list%coord_type) == 'cartesian') Then
      Do i = 1 , model_data%input%num_atoms
        model_data%input%atom(i)%r=model_data%input%scale_factor_vasp*model_data%input%atom(i)%r
      End Do
    Else If (Trim(model_data%input%list%coord_type) == 'direct') Then
      Do i = 1, model_data%input%num_atoms 
        v_cart=MatMul(model_data%input%atom(i)%r, model_data%input%cell)
        model_data%input%atom(i)%r=model_data%input%scale_factor_vasp*v_cart
      End Do
    End If
    
    model_data%input%cell= model_data%input%scale_factor_vasp * model_data%input%cell

  End Subroutine read_input_vasp_format

  Subroutine compute_stoichiometry_input(stoich_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to 
    ! 1) identify the number of the species (defined in &block_species)
    ! within the input model 
    ! 2) computes the resulting stoichiometry  
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),   Intent(In   ) :: stoich_data
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: i, j, k
    Character(Len=256) :: messages(2)
    Logical            :: loop
    Real(Kind=wp)      :: nsp, nsp_round
    Integer(Kind=wi)   :: nsp_int, imin

    ! Determine the number of components for each species of the input model
    Do i=1, stoich_data%total_tags
      loop=.True.
      j=1
      Do While (j <= stoich_data%num_species%value .And. loop)
        k=1
        Do While (k <= stoich_data%species(j)%num_components .And. loop)
          If (Trim(stoich_data%species(j)%component%tag(k))==model_data%block_component%tag(i)) Then
            model_data%input%species(j)%component%N(k)=model_data%block_component%N0(i)
            loop=.False.
          End If
          k=k+1
        End Do
        j=j+1
      End Do
    End Do

    ! For each species:
    ! 1) check that the number of atoms provided makes an interget number of units
    ! 2) find the number of units. 
    Do j=1, stoich_data%num_species%value
      model_data%input%species(j)%num=0
      Do k= 1, stoich_data%species(j)%num_components
        If ((Trim(model_data%input%species(j)%topology) == 'crystal'  .Or. &
             Trim(model_data%input%species(j)%topology) == 'solute') .And. &
            model_data%input%species(j)%component%N(k) == 0) Then
          Write (messages(1),'(5a)') '***ERROR: Component "',  Trim(model_data%input%species(j)%component%tag(k)),&
                                   & '" of "fixed" species "', Trim(model_data%input%species(j)%tag), &
                                   & '" cannot have zero atoms in the input structure. Please review settings.'          
          Call info(messages,1)
          Call error_stop(' ')    
        End If
        model_data%input%species(j)%num=model_data%input%species(j)%num+model_data%input%species(j)%component%N(k)
      End Do
      nsp=Real(model_data%input%species(j)%num,Kind=wp)/stoich_data%species(j)%atoms_per_species
      nsp_int=Int(nsp)
      nsp_round=Real(nsp_int,Kind=wp)
      If (Abs(nsp-nsp_round) > epsilon(nsp)) Then
        Write (messages(1),'(3a,f12.5)') '***ERROR: Amount of atoms provided for species "',&
                                 & Trim(model_data%input%species(j)%tag),&
                                 & '" does not make an integer number of units for the input model.&
                                 & Computed number of species is ', nsp    
        Write (messages(2),'(a)')  'Please check consistent between the definition of &block_species_component and&
                                  & &block_input_composition, as well as the atoms provided in input structure'  
        Call info(messages,2)
        Call error_stop(' ')  
      End If
        model_data%input%species(j)%num=model_data%input%species(j)%num/stoich_data%species(j)%atoms_per_species
    End Do

    ! For all the species of topology "crystal" or "solute" find that particular species with minimum number of units 
    model_data%input%min_species=Huge(1)
    Do i = 1, stoich_data%num_species%value
      If (Trim(model_data%input%species(i)%topology)=='crystal' .Or.&
         Trim(model_data%input%species(i)%topology) == 'solute') Then
        If (model_data%input%species(i)%num==0) Then
          Write (messages(1),'(3a)') '***ERROR: the input number for "fixed" species "', Trim(model_data%input%species(i)%tag),&   
                                 & '" is zero. This is a wrong setting for a "fixed" species. Please change.'
          Call info(messages,1)
          Call error_stop(' ')    
        End If
        If (model_data%input%species(i)%num < model_data%input%min_species) Then
          imin=i
          model_data%input%min_species=model_data%input%species(i)%num
        End If
      End If
    End Do     

    model_data%input%min_stoich=stoich_data%species(imin)%s0

    ! For all the species of topology "crystal" or "solute" find the minimum number of atomic components
    model_data%input%min_components=Huge(1)
    Do i = 1, stoich_data%num_species%value
      If (Trim(model_data%input%species(i)%topology)=='crystal' .Or.&
          Trim(model_data%input%species(i)%topology) == 'solute') Then
        Do k=1, model_data%input%species(i)%num_components
          If (stoich_data%species(i)%component%N0(k) < model_data%input%min_components) Then
            model_data%input%min_components=stoich_data%species(i)%component%N0(k)
          End If
        End Do
      End If
    End Do     


    ! Find the stoichiometry of the involved species within the input model 
    Do i = 1, stoich_data%num_species%value
       model_data%input%species(i)%s=Real(model_data%input%species(i)%num,Kind=wp)/model_data%input%min_species
       model_data%input%species(i)%s=stoich_data%species(imin)%s0_pristine*model_data%input%species(i)%s
       If (Trim(model_data%input%species(i)%topology)/='crystal' .Or. &
           Trim(model_data%input%species(i)%topology) == 'solute') Then
         model_data%input%species(i)%s=model_data%input%species(i)%s/model_data%input%min_components
       End If
    End Do     

    model_data%ref_stoich=stoich_data%species(imin)%s0_pristine

    ! Check that the stochiometry of the fixed species is consistent with the sotichiometry defined in &block_species
    Do i = 1, stoich_data%num_species%value
      If (Trim(model_data%input%species(i)%topology)=='crystal' .Or. &
          Trim(model_data%input%species(i)%topology) == 'solute') Then
        If (Abs(model_data%input%species(i)%s-stoich_data%species(i)%s0_pristine)>model_data%stoichiometry_error%value) Then
          Write (messages(1),'(a)') '***ERROR: the stoichiometry for the "fixed" species of the input model does not match&
                                 & the stoichiometry defined in &block_species. Please review the input data.&
                                 & If needed, the user could change the species from "fixed" to "dependent/independent".'
          Call info(messages,1)
          Call error_stop(' ')    
        End If
      End If
    End Do

  End Subroutine compute_stoichiometry_input

  Subroutine identify_species_input(stoich_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to group the species by creating a list with the involved atoms
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),   Intent(In   ) :: stoich_data
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: i, j, k, m
    Logical               :: loop

    Do i = 1, model_data%types_species
      If (model_data%input%species(i)%change_content) Then
        Do j =1, model_data%input%species(i)%num
          loop=.True.
          k=1
          Do While (k <= model_data%input%num_atoms .And. loop)
            m=1
            Do While (m <= model_data%input%species(i)%num_components .And. loop)
              If (Trim(model_data%input%atom(k)%tag)==Trim(model_data%input%species(i)%component%tag(m)) .And. &
                   (.Not. model_data%input%atom(k)%in_species)) Then
                loop=.False.
                model_data%input%species(i)%units(j)%list(1)=k
                If (model_data%input%species(i)%atoms_per_species > 1) Then
                  Call find_neighbours(stoich_data%species(i)%bond_cutoff, model_data, k, j, i)
                Else
                  model_data%input%atom(k)%in_species=.True. 
                End If
              End If
              m=m+1
            End Do
            k=k+1
          End Do
        End Do
      End If
    End Do

  End Subroutine identify_species_input 


  Subroutine find_neighbours(bond_cutoff, model_data, kin, jin, iin)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to find the all the neighbours-bonded atoms that form a species 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp),     Intent(In   ) :: bond_cutoff
    Type(model_type),  Intent(InOut) :: model_data
    Integer(Kind=wi),  Intent(In   ) :: kin, jin, iin

    Integer(Kind=wi) :: m, k, ka, il
    Integer(Kind=wi) :: knn, nacc, num_nn, num_tot, nlist, knni
    Logical             :: loop, loop_sp, is_nn, is_list
    Real(Kind=wp)       :: a(3), b(3)
    Real(Kind=wp)       :: dist
    Character(Len=256)  :: messages(6)

    Integer(Kind=wi), Dimension(max_nn) :: list_nn, list_nn_next

    num_nn=1
    list_nn(1)=kin

    loop=.True. 
    num_tot=0
    list_nn_next=0 
    nlist=1

    Do While (loop)
      nacc=0
      Do knn=1, num_nn
        ka=list_nn(knn)
        Do k=1, model_data%input%num_atoms
          loop_sp=.True.
          m=1
          Do While (m <= model_data%input%species(iin)%num_components .And. loop_sp)
            If (Trim(model_data%input%atom(k)%tag)==Trim(model_data%input%species(iin)%component%tag(m)) .And. k/=ka .And.&
               (.Not. model_data%input%atom(k)%in_species) .And. (.Not. model_data%input%atom(ka)%in_species) ) Then
               ! Compute distances
               a=model_data%input%atom(k)%r(:)
               b=model_data%input%atom(ka)%r(:)
               Call compute_distance_PBC(a, b, model_data%input%cell, model_data%input%invcell, dist)
               ! Check if distances are within a given cutoff
               If (dist<bond_cutoff) Then 
                 If (dist<min_intra) Then
                   Write (messages(1),'(a)') '***ERROR in the input structure: '
                   Write (messages(2),'(3(a,i4),3a)') 'Intermolecular distance between atoms ', k, ' and ', ka, &
                                                    ' in unit ', jin, ' of species "',&
                                                    & Trim(model_data%input%species(iin)%tag),'"'
                   Write (messages(3),'(a,f4.2,a)') 'is shorter than the input minimum distance criteria of ', &
                                                   &  min_intra,&
                                                   & ' Angstrom for bonding. Please review the input geometry.'
                   Call info(messages,3)
                   Call error_stop(' ')
                 Else
                   loop_sp=.False.
                   is_nn=.True.
                   il=1
                   Do While (il <= max_nn .And. is_nn) 
                     If (k==list_nn_next(il)) Then
                       is_nn=.False.
                     End If
                     il=il+1
                   End Do 
                   If (is_nn) Then    
                     nacc=nacc+1
                     If (nacc==max_nn+1) Then
                        Write (messages(1),'(a,i4,3a)') '***ERROR: Trouble to identify unit ', jin, ' of species "',&
                                                 & Trim(model_data%input%species(iin)%tag),&
                                                 & '" from the input coordinates of the input structure.'    
                        Write (messages(2),'(a)') ' It is likely the species represents a fixed crystal structure with bonds&
                                                 & periodically repeated. If this is the case, the user must define the&
                                                 & species as "fixed" in &Block_species'
                        Call info(messages, 2)
                        Call error_stop(' ')
             
                     End If
                     list_nn_next(nacc)=k
                     il=1
                     is_list=.True.
                     Do While (il <= model_data%input%species(iin)%atoms_per_species .And. is_list) 
                       If (k==model_data%input%species(iin)%units(jin)%list(il)) Then
                         is_list=.False.
                       End If 
                       il=il+1
                     End Do
                     If (is_list) Then
                       nlist=nlist+1
                       model_data%input%species(iin)%units(jin)%list(nlist)=k
                     End If
                   End If
                 End If  
               End If
            End If
            m=m+1
          End Do
        End Do
        If (.Not. model_data%input%atom(ka)%in_species) Then
          model_data%input%atom(ka)%in_species=.True.
          num_tot=num_tot+1
          If (num_tot+nacc == model_data%input%species(iin)%atoms_per_species) Then
            Do knni= 1, nacc
              model_data%input%atom(list_nn_next(knni))%in_species=.True.
            End Do
            loop=.False.
            num_tot=model_data%input%species(iin)%atoms_per_species
          End If
        End If
      End Do
      num_nn=nacc  
      list_nn=list_nn_next
      list_nn_next=0 
      If (nacc==0) loop=.False.
    End Do

    If (model_data%input%species(iin)%atoms_per_species/=num_tot) Then
      Write (messages(1),'(a,i4,3a)') '***ERROR: Trouble to identify unit ', jin, ' of species "',&
                               & Trim(model_data%input%species(iin)%tag),&
                               &'" from the input coordinates of the input structure.'
      Write (messages(2),'(a)')  'Possible reasons:' 
      Write (messages(3),'(a)')  ' 1) wrong input geometry, where one(or more) species are dissociated (e.g. H3O instead of H2O).&
                           & If this is the case, the user must provide an input model consistent with the definition of species.'

      Write (messages(4),'(a,f4.2,a)')  ' 2) the value for the bond cutoff (', bond_cutoff, &
                                    & ' Angstrom) as specified in &block_species_components&
                                    & is not adequate. The user must adjust this value'

      Write (messages(5),'(a)') ' 3) the species represents a fixed crystal structure with bonds periodically repeated.&
                               & If this is the case, the user must define the species as "fixed" in &Block_species'
      Write (messages(6),'(a)') ' 4) the cell vectors are not strictly consistent with the geometry of the input structure.'
      Call info(messages, 6)
      Call error_stop(' ')
    End If

  End Subroutine find_neighbours

  Subroutine compute_distance_PBC(a, b, cell, invcell, dist)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the distance between atoms with PCB 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   ) :: a(3), b(3)
    Real(Kind=wp), Intent(In   ) :: cell(3,3)
    Real(Kind=wp), Intent(In   ) :: invcell(3,3)
    Real(Kind=wp), Intent(  Out) :: dist

    Real(Kind=wp) :: Dr_cart(3)
    Logical :: modified
    Integer :: ir

    ! Vector difference
    Do ir=1,3
      Dr_cart(ir)=a(ir)-b(ir) 
    End Do

    ! Find the vector difference for the nearest neighbours (NN)
    Call check_PBC(Dr_cart, cell, invcell, 0.5_wp, modified)
    ! Calculate norm
    dist=norm2(Dr_cart)

  End Subroutine compute_distance_PBC

  Subroutine check_PBC(v_cart, basis, inv_basis, ratio, changed_geo)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the components of given vector "v_cart" 
    ! (given in cartesian coordinates of the Euclidean space) in terms of a general
    ! 3D vector basis, specified by the 3x3 matrix "basis".
    ! The inverse of the basis matrix is named "inv_basis".
    ! The factor "ratio" is used to evaluate how large these components are, 
    ! and modify the "vect" accordingly to accound for PBC.
    !
    ! The input value of "ratio" depends on the quantity to evalue. Thus,
    ! - ratio=0.5 is used for nearest-neighbour distances 
    ! - ratio=1.0 is used to evaluate if atomic positions lie within the volume 
    !   defined by the basis.
    !
    ! Logical variable "changed_geo" is used to check is the vector has been modified 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Intent(InOut) :: v_cart(3)
    Real(Kind=wp), Intent(In   ) :: basis(3,3)
    Real(Kind=wp), Intent(In   ) :: inv_basis(3,3)
    Real(Kind=wp), Intent(In   ) :: ratio 
    Logical,       Intent(  Out) :: changed_geo                     

    Real(Kind=wp) :: v_direct(3), limit1, limit2
    Integer(Kind=wi) :: ir, i
    Logical          :: flag

    changed_geo=.False.
 
    If (Abs(ratio-0.5_wp) < epsilon(ratio)) Then
      limit1=ratio
      limit2=-ratio
    Else If (Abs(ratio-1.0_wp) < epsilon(ratio)) Then
      limit1=ratio+length_tol
      limit2=-length_tol
    End If

    ! Express vector difference in terms of the cell vectors
    v_direct= MatMul(v_cart, inv_basis)

    ! PCB effect
    i=1
    flag=.True.
    Do While (i < 4 .And. flag)
      Do ir = 1, 3
        If (v_direct(ir) > limit1) Then
           v_cart(:)= v_cart(:) - basis(ir,:)
           changed_geo=.True.
        Else If (v_direct(ir) < limit2) Then
           v_cart(:)= v_cart(:) + basis(ir,:)
           changed_geo=.True.
        End If
        If (changed_geo) Then
          v_direct= MatMul(v_cart, inv_basis)
        End If
      End Do
      If (.Not. changed_geo) Then
        flag=.False.      
      End If        
      i=i+1
    End Do

  End Subroutine check_PBC 

  Subroutine check_orthorhombic_cell(A, flag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check if the simulation cell 'A' is orthorhombic or not. 
    ! vectors. If A is not orthorhombic, flag will be set to False
    !
    ! author    - i.scivetti Aug 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   )  :: A(3,3)
    Logical,       Intent(  Out)  :: flag
   
    Integer(Kind=wi) :: i, j

    i=1
    flag=.True.
    Do i = 1, 2
     Do j = i+1, 3
       If (Abs(Dot_product(A(i,:), A(j,:)))>epsilon(1.0_wp)) Then
         flag=.False.      
       End If        
     End Do
    End Do 

  End Subroutine check_orthorhombic_cell

  Subroutine about_cell(A,invA,length)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to invert a general 3x3 matrix and compute the length of the cell
    ! vectors.
    ! Matrix A is the input matrix.
    ! Matrix invA is the output matrix (inverse of A) 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Intent(In   )  :: A(3,3)
    Real(Kind=wp), Intent(  Out)  :: invA(3,3)
    Real(Kind=wp), Intent(  Out)  :: length(3)

    Real(Kind=wp) :: Det
    Real(Kind=wp) :: Cofactor(3,3)
   
    Integer(Kind=wi) :: i, j

    length = 0.0_wp     

    Do i = 1, 3 
      Do j= 1, 3
        length(i) = length(i)+ A(i,j)**2   
      End Do
      length(i)=sqrt(length(i))
    End Do 

    Det =   A(1,1)*A(2,2)*A(3,3)  &
          - A(1,1)*A(2,3)*A(3,2)  &
          - A(1,2)*A(2,1)*A(3,3)  &
          + A(1,2)*A(2,3)*A(3,1)  &
          + A(1,3)*A(2,1)*A(3,2)  &
          - A(1,3)*A(2,2)*A(3,1)

    If (Abs(Det) <= epsilon(det)) Then
      Call error_stop('***ERROR: The determinant of the simulation cell is zero. Please check the definition of the cell vectors')
    End If

    Cofactor(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    Cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    Cofactor(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    Cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    Cofactor(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    Cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    Cofactor(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    Cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    Cofactor(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

    invA = Transpose(Cofactor) / Det

  End Subroutine about_cell

  Subroutine init_random_seed()
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to generate the seed for random number generation 
    !
    ! author    - i.scivetti Feb 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi) :: i, n, clock
    Integer(Kind=wi), DIMENSION(:), Allocatable :: seed
                                      
    Call random_seed(size = n)
    
    Allocate(seed(n))
    Call system_clock(Count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    Call random_seed(put = seed)
    Deallocate(seed)

  End Subroutine init_random_seed

  Subroutine generate_hpc_simulation_files(code_format, files, stoich_data, model_data, simulation_data, hpc_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to generate simulation/HPC files without having to generate 
    ! atomisitic models 
    !
    ! author    - i.scivetti June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*),  Intent(In   ) :: code_format 
    Type(file_type),   Intent(InOut) :: files(:) 
    Type(stoich_type), Intent(InOut) :: stoich_data
    Type(model_type),  Intent(InOut) :: model_data
    Type(simul_type),  Intent(InOut) :: simulation_data
    Type(hpc_type),    Intent(In   ) :: hpc_data     
  
    Character(Len=256) :: messages(5)
    Character(Len=256) :: type_sim, file_record, model_name
    Character(Len=256) :: saved_format, set_error
    Integer(Kind=wi)   :: i, j, iunit, io, ifolder, ifile
    Logical             :: safe, loop, fhpc, fsim, loop_pp
    Character(Len=256)  :: exec_cat 
    Character(Len=256)  :: pseudo_list(max_components) 

    file_record=Trim(files(FILE_RECORD_MODELS)%filename)
    set_error='***ERROR: file '//Trim(FOLDER_RESTART)//'/'//Trim(file_record)
    fhpc=.False.
    fsim=.False.
    
    Call info(' ', 1)
    Call info(' Starting the calculation', 1)
    Call info(' ========================', 1)
    Inquire(File=Trim(FOLDER_RESTART)//'/'//Trim(file_record), Exist=safe)
    If (.Not. safe) Then
      Write (messages(1), '(1x,2a)') Trim(set_error), ' does not exist.'
      Write (messages(2), '(1x,a)')  'This file is generated when the atomistic models are built.'
      Write (messages(3), '(1x,a)')  'Have you generated the models previously?&
                                    & Have you deleted/moved/renamed this file by mistake?'
      Write (messages(4), '(1x,a)')  'In any case, the requested analysis is not possible. If any simulation or hpc&
                                    & setting needs adjustment, we recommend to generate the models again.'

      Call info(messages, 4)
      Call error_stop(' ')
    End If

    Open(Newunit=files(FILE_RECORD_MODELS)%unit_no, File=Trim(FOLDER_RESTART)//'/'//Trim(file_record),Status='Old')
    iunit=files(FILE_RECORD_MODELS)%unit_no

    Read (iunit, Fmt=*, iostat=io) type_sim, saved_format
    If (is_iostat_end(io)) Then
      Write (messages(1), '(1x,2a)') Trim(set_error), ' is empty! The requested analysis is not possible'
      Call info(messages, 1)
      Call error_stop(' ')
    End If

    Read (iunit, Fmt='(a)', iostat=io) model_data%sample%path
    If (is_iostat_end(io)) Then
      Write (messages(1), '(1x,2a)') Trim(set_error),' has been corrupted or modified!&
                                & The requested analysis is not possible'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    Backspace iunit

    If (Trim(saved_format)/=Trim(code_format))then
      Write (messages(1), '(1x,3a)') '***ERROR: atomistic models were generated according to code "', Trim(saved_format), '".'
      Write (messages(2), '(1x,3a)') 'However, the user has now selected option "', Trim(code_format), '" for directive&
                                    & "output_model_format". Please adjust.' 
      Write (messages(3), '(1x,3a)')  'If the user still wants to generate files for option "', Trim(code_format),&
                                    '", atomistic models must be re-built.'
      Call info(messages, 3) 
      Call error_stop(' ')
    End If

    ! simulation cell
    Call read_input_model(files, stoich_data, model_data)
    Call define_cell(model_data, simulation_data)

    loop=.True.
    Do While (loop) 
      ! Read path
      Read (iunit, Fmt='(a)', iostat=io) model_data%sample%path
      ! Read number of net elements
      Read (iunit, Fmt=*, iostat=io) model_data%sample%list%net_elements
      Call check_record_file(io, file_record)
      ! Read species tag 
      Read (iunit, Fmt=*, iostat=io) (model_data%sample%list%tag(i), i=1, model_data%sample%list%net_elements)
      Call check_record_file(io, file_record)
      ! Read species elements 
      Read (iunit, Fmt=*, iostat=io) (model_data%sample%list%element(i), i=1, model_data%sample%list%net_elements)
      Call check_record_file(io, file_record)
      ! Read number of atoms per species
      Read (iunit, Fmt=*, iostat=io) (model_data%sample%list%N0(i), i=1, model_data%sample%list%net_elements)
      Call check_record_file(io, file_record)

      Call execute_command_line('[ -d '//Trim(model_data%sample%path)//' ]', exitstat=ifolder)

      If (ifolder/=0) Then 
        Write (messages(1), '(1x,3a)') '***WARNING: folder ', Trim(model_data%sample%path), ' cannot be found.'
        Write (messages(2), '(1x,3a)') 'The requested analysis cannot be conducted for this particular model.'
        Call info(messages, 2)
      Else
        Write (messages(1), '(1x,2a)') '=== Folder ', Trim(Adjustl(model_data%sample%path))
        Call info(messages, 1)
        model_name='SAMPLE.'//Trim(code_format)
        Call execute_command_line('[ -f '//Trim(Adjustl(model_data%sample%path))//'/'//Trim(model_name)//&
                                 &' ]', exitstat=ifile)
        If (ifile/=0) Then
          Write (messages(1), '(1x,3a)') '***PROBLEMS: Model file ', Trim(model_name), ' cannot be found.'
          Write (messages(2), '(1x,3a)') '             The requested analysis cannot be conducted for this particular model.'
          Write (messages(3), '(1x,3a)') '             The user must re-build the atomistic model for this composition.'
          Call info(messages, 3)
        Else 
          If (simulation_data%generate) Then    
             If (simulation_data%dft%need_vdw_kernel) Then
                Call ammend_files(model_data%sample%path, Trim(FOLDER_DFT)//'/'//simulation_data%dft%vdw_kernel_file,&
                                  simulation_data%dft%vdw_kernel_file, '&dft_settings', fsim)
             End If     
 
            If (Trim(code_format)=='vasp') Then
              Call print_vasp_settings(files, model_data%sample%list%net_elements, model_data%sample%list%element,&
                                     & model_data%sample%list%tag, model_data%sample%list%N0, simulation_data)
              Call ammend_files(model_data%sample%path, files(FILE_SET_SIMULATION)%filename, 'INCAR', &
                                & '&block_simulation_settings', fsim) 
              Call execute_command_line('rm '//Trim(files(FILE_SET_SIMULATION)%filename))
              Call ammend_files(model_data%sample%path, files(FILE_KPOINTS)%filename, 'KPOINTS', &
                                & '&dft_settings', fsim)
              Call execute_command_line('rm '//Trim(files(FILE_KPOINTS)%filename))
              If (simulation_data%dft%pp_info%stat) Then           
                Call ammend_files(model_data%sample%path, 'POTCAR', 'POTCAR', '&dft_settings', fsim)
                Call execute_command_line('rm POTCAR')
              End If

            Else If (Trim(code_format)=='cp2k') Then
              Call print_cp2k_settings(files, model_data%sample%list%net_elements, model_data%sample%list%element, & 
                                     & model_data%sample%list%tag, model_data%sample%list%N0, simulation_data)
              Call ammend_files(model_data%sample%path, files(FILE_SET_SIMULATION)%filename, 'input.cp2k', &
                                & '&block_simulation_settings', fsim)
              If (simulation_data%dft%basis_info%stat)Then          
                Call ammend_files(model_data%sample%path, Trim(FOLDER_DFT)//'/'//'BASIS_SET', 'BASIS_SET', &
                                & '&dft_settings', fsim)
              End If    
              ! Delete temporary files
              Call execute_command_line('rm '//Trim(files(FILE_SET_SIMULATION)%filename))
              ! Checking PPs
              If (simulation_data%dft%pp_info%stat) Then
                Do i=1,  model_data%sample%list%net_elements 
                  j=1
                  loop_pp=.True.
                  Do While (j <= simulation_data%total_tags .And. loop_pp)
                    If (Trim(model_data%sample%list%element(i))==Trim(simulation_data%dft%pseudo_pot(j)%element)) Then
                      pseudo_list(i)= Trim(simulation_data%dft%pseudo_pot(j)%file_name)
                      loop_pp=.False.
                      Call execute_command_line('[ -f '//Trim(Adjustl(model_data%sample%path))//'/'//Trim(pseudo_list(i))//&
                                 &' ]', exitstat=ifile) 
                      If (ifile /= 0) Then
                        Write (messages(1), '(1x,a)') 'Generating PP file '//Trim(pseudo_list(i))//&
                                                     &' according to new specifications'
                        Call info(messages, 1)
                        Call execute_command_line('cp DFT/PPs/'//Trim(pseudo_list(i))//' '//Trim(Adjustl(model_data%sample%path)))
                        fsim=.True.
                      End If
                    End If
                    j=j+1
                  End Do
                End Do              
              End If ! End checking PPs               

            Else If (Trim(code_format)=='castep') Then
              Call print_castep_settings(files, model_data%sample%list%net_elements, model_data%sample%list%element,&
                                     & model_data%sample%list%tag, model_data%sample%list%N0, simulation_data)
              exec_cat='cat '//Trim(files(FILE_SET_SIMULATION)%filename)//' '//Trim(model_data%sample%path)//'/SAMPLE.castep '&
                      &//'> model.cell'
              Call execute_command_line(exec_cat)
              Call ammend_files(model_data%sample%path, 'model.cell', 'model.cell', '&block_simulation_settings', fsim) 
              Call ammend_files(model_data%sample%path, 'model.param', 'model.param', '&block_simulation_settings', fsim)
              ! Delete temporary files
              Call execute_command_line('rm '//Trim(files(FILE_SET_SIMULATION)%filename)//' model.cell model.param' ) 
              
              ! Checking PPs
              If (simulation_data%dft%pp_info%stat) Then
                Do i=1,  model_data%sample%list%net_elements 
                  j=1
                  loop_pp=.True.
                  Do While (j <= simulation_data%total_tags .And. loop_pp)
                    If (Trim(model_data%sample%list%element(i))==Trim(simulation_data%dft%pseudo_pot(j)%element)) Then
                      pseudo_list(i)= Trim(simulation_data%dft%pseudo_pot(j)%file_name)
                      loop_pp=.False.
                      Call execute_command_line('[ -f '//Trim(Adjustl(model_data%sample%path))//'/'//Trim(pseudo_list(i))//&
                                 &' ]', exitstat=ifile)
                      If (ifile /= 0) Then
                        Write (messages(1), '(1x,a)') 'Generating PP file '//Trim(pseudo_list(i))//&
                                                    & ' according to new specifications'
                        Call info(messages, 1)
                        Call execute_command_line('cp DFT/PPs/'//Trim(pseudo_list(i))//' '//Trim(Adjustl(model_data%sample%path)))
                        fsim=.True.
                      End If
                    End If
                    j=j+1
                  End Do
                End Do              
              End If ! End checking PPs               

            Else If (Trim(code_format)=='onetep') Then
              Call print_onetep_settings(files, model_data%sample%list%net_elements, model_data%sample%list%tag,&
                                         model_data%sample%list%N0, simulation_data)
              exec_cat='cat '//Trim(files(FILE_SET_SIMULATION)%filename)//' '//Trim(model_data%sample%path)//'/SAMPLE.onetep '&
                      &//'> model.dat'
              Call execute_command_line(exec_cat)
              Call ammend_files(model_data%sample%path, 'model.dat', 'model.dat', '&block_simulation_settings', fsim)
              ! Delete temporary files
              Call execute_command_line('rm '//Trim(files(FILE_SET_SIMULATION)%filename)//' model.dat' )
              ! Checking PPs
              If (simulation_data%dft%pp_info%stat) Then
                Do i=1,  model_data%sample%list%net_elements 
                  j=1
                  loop_pp=.True.
                  Do While (j <= simulation_data%total_tags .And. loop_pp)
                    If (Trim(model_data%sample%list%element(i))==Trim(simulation_data%dft%pseudo_pot(j)%element)) Then
                      pseudo_list(i)= Trim(simulation_data%dft%pseudo_pot(j)%file_name)
                      loop_pp=.False.
                      Call execute_command_line('[ -f '//Trim(Adjustl(model_data%sample%path))//'/'//Trim(pseudo_list(i))//&
                                 &' ]', exitstat=ifile) 
                      If (ifile /= 0) Then
                        Write (messages(1), '(1x,a)') 'Generating PP file '//Trim(pseudo_list(i))//&
                                                    & ' according to new specifications'
                        Call info(messages, 1)
                        Call execute_command_line('cp DFT/PPs/'//Trim(pseudo_list(i))//' '//Trim(Adjustl(model_data%sample%path)))
                        fsim=.True.
                      End If
                    End If
                    j=j+1
                  End Do
                End Do              
              End If ! End checking PPs               
            End If
          End If 
 
          If (hpc_data%generate) Then
            Call ammend_files(model_data%sample%path, files(FILE_HPC_SETTINGS)%filename, hpc_data%script_name, &
                             & '&block_hpc_settings', fhpc) 
          End If

        End If
      End If

      ! Check for the end of file
      Read (iunit, Fmt='(a)', iostat=io) model_data%sample%path
      If (is_iostat_end(io)) Then
        loop=.False.
      Else
        Backspace iunit
      End If
    End Do


    If (simulation_data%generate) Then
      Call info(' ',1)
      If (.Not. fsim) Then
        Write (messages(1), '(1x,a)') 'INFO: No simulation file has been generated/updated.&
                                     & This means there was no change with respect to the previous simulation settings. '
        Call info(messages, 1)
      Else
        Call summary_simulation_settings(simulation_data)
      End If
    End If
    If (hpc_data%generate) Then
      If (.Not. fhpc) Then
        Call info(' ',1)
        Write (messages(1), '(1x,a)') 'INFO: No HPC script file has been generated/updated.&
                                     & This means there was no change with respect to the previous HPC settings. '
        Call info(messages, 1)
      Else
        Call summary_hpc_settings(hpc_data)
      End If
    End If

    If ((.Not. fhpc) .And. hpc_data%generate .And. (.Not. fsim) .And. simulation_data%generate) Then 
      Call info(' ',1)
      Write (messages(1), '(1x,a)') '***WARNING: No simulation nor HCP script file has generated/changed.&
                                  & Is the user sure about having added/modified the blocks for&
                                  & simulations/HPC settings?'
      Call info(messages, 1)
    End If

    If (simulation_data%generate) Then
      ! Print warnings
      Call warning_simulation_settings(simulation_data)
    End If

  End Subroutine generate_hpc_simulation_files

  Subroutine check_record_file(io, filename)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check RECORD_MODELS file
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: io
    Character(Len=*), Intent(In   ) :: filename

    Character(Len=256)    :: message

    If (io /= 0 .Or. is_iostat_end(io)) Then
      Write (message, '(1x,3a)') '***ERROR: file ', Trim(filename), ' has been corrupted or modified!&
                              & The requested analysis is not possible.'
      Call error_stop(message)
    End If
 
  End Subroutine check_record_file 


  Subroutine ammend_files(folder, file_new, file_ref, block, flag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check between files and decide what to do
    !
    ! author    - i.scivetti June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: folder
    Character(Len=*), Intent(In   ) :: file_new
    Character(Len=*), Intent(In   ) :: file_ref
    Character(Len=*), Intent(In   ) :: block
    Logical,          Intent(InOut) :: flag 

    Character(Len=256)  :: message
    Character(Len=256)  :: path, exec_cp
    Integer(Kind=wi)    :: ifile

    path=Trim(folder)//'/'//Trim(file_ref)
    exec_cp='cp '//Trim(file_new)//' '//Trim(path)
    Call execute_command_line('[ -f '//Trim(path)//' ]', exitstat=ifile)
    If (ifile/=0) Then
      Write (message, '(1x,a)') 'Generate file '//Trim(Adjustl(file_ref))//' according to '//Trim(block)
      Call info(message, 1)
      Call execute_command_line(exec_cp)
      flag=.True.
    Else
      Call execute_command_line('cmp -s '//Trim(file_new)//Trim(path), exitstat=ifile) 
      If (ifile/=0) Then
        Write (message, '(1x,a)') 'Updating file '//Trim(Adjustl(file_ref))//' according to changes in '//Trim(block)
        Call info(message, 1)
        Call execute_command_line(exec_cp)
        flag=.True.
      Else
        Write (message, '(1x,a)') 'No need to change file '//Trim(Adjustl(file_ref))
        Call info(message, 1)
      End If             
    End If


  End Subroutine ammend_files

End Module atomistic_models
