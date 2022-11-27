!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that defines simulation type and procedure
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author    - i.scivetti Oct  2021
! Refact    - i.scivetti Sept 2022
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module simulation_setup

  Use constants,        Only : max_components
  Use input_types,      Only : in_integer, &
                               in_integer_array,&
                               in_logic,   &
                               in_string,  &
                               in_param,   &
                               in_param_array, &
                               in_scalar
  Use numprec,          Only : wi, &
                               wp
  Use stoichiometry,    Only : stoich_type
  Use process_data,     Only : capital_to_lower_case, &
                               remove_symbols  
  Use unit_output,      Only : error_stop,&
                               info 

  Implicit None
  Private

  ! Maximum directives for simulations
  Integer(Kind=wi), Parameter, Public  :: max_directives=100 
  ! Maximum number of Boltzmann ions
  Integer(Kind=wi), Parameter, Public  :: max_number_boltzmann_ions=20

  ! Components inherited from model_data
  Type :: component_in_block
    Character(Len=8) :: tag
    Character(Len=2) :: element
    Integer(Kind=wi) :: atomic_number
    Real(Kind=wp)    :: valence_electrons
  End Type component_in_block

  ! NGWF type 
  Type :: type_ngwf
    Character(Len=8) :: tag
    Character(Len=2) :: element
    Integer(Kind=wi) :: ni
    Real(Kind=wp)    :: radius
  End Type type_ngwf

  ! Type for pseudopotentials 
  Type :: type_pseudo
    Character(Len=256) :: file_name
    Character(Len=256) :: potential 
    Character(Len=8)   :: tag
    Character(Len=2)   :: element
  End Type type_pseudo

  ! Type for extra directives 
  Type, Public :: type_extra 
    Character(Len=256) :: array(max_directives)
    Character(Len=256) :: key(max_directives)
    Character(Len=256) :: set(max_directives)
    Integer(Kind=wi)   :: N0
  End Type type_extra

  ! type for reference data
  Type, Public :: type_ref_data
    Character(Len=256) :: key
    Character(Len=16)  :: keytype
    Character(Len=256) :: msn
    Character(Len=256) :: set_default
    Character(Len=16)  :: units
    Integer(Kind=wi)   :: N0
  End Type type_ref_data

  ! Type for basis_set
  Type :: type_basis
    Character(Len=8)   :: tag
    Character(Len=2)   :: element
    Character(Len=256) :: basis
    Character(Len=256) :: type
  End Type type_basis

  ! Type for magnetization
  Type :: type_mag
    Character(Len=8)  :: tag
    Character(Len=2)  :: element
    Real(Kind=wp)     :: value
  End Type type_mag

  ! Type for Hubbard corrections 
  Type :: type_hubbard
    Character(Len=8)  :: tag
    Character(Len=2)  :: element
    Integer(Kind=wi)  :: l_orbital  
    Real(Kind=wp)     :: U 
    Real(Kind=wp)     :: J 
  End Type type_hubbard

  ! Type to read quantities from a list of species from two-rows blocks
  Type, Public :: species_list
    Character(Len=8)  :: tag
    Character(Len=2)  :: element
    Real(Kind=wp)     :: value
  End Type species_list

  ! Type for GC-DFT
  Type :: type_gcdft 
    Type(in_logic)  ::  activate     
    Type(in_param)  ::  reference_potential    
    Type(in_param)  ::  electrode_potential
    Type(in_scalar) ::  electron_threshold    
  End Type
 
  ! Type for boltzmann_ions 
  Type :: boltzmann_ions
    Character(Len=8)  :: tag
   !Character(Len=2)  :: element
    Real(Kind=wp)     :: charge
    Real(Kind=wp)     :: conc
    Real(Kind=wp)     :: necs_shift
  End Type boltzmann_ions 

  ! Type for DFT settings
  Type :: dft_type
    ! Flag to ensure DFT block is not defined more than once 
    Logical  :: generate=.False.
    ! Type XC functional
    Type(in_string)  :: xc_level     
    ! Type XC version 
    Type(in_string)  :: xc_version     
    ! Staring XC base approach 
    Character(Len=256) :: xc_base
    ! XC reference
    Character(Len=256) :: xc_ref
    ! vdW
    Type(in_string)    :: vdw   
    ! vdW reference
    Character(Len=256) :: vdw_ref
    ! vdW kernel
    Logical            :: need_vdw_kernel 
    Character(Len=256) :: vdw_kernel_file
    ! Flag to set spin polarised simulation 
    Type(in_logic)  :: spin_polarised
    ! Energy cutoff
    Type(in_param)  :: encut
    ! Precision
    Type(in_string) :: precision   
    ! Smearing
    Type(in_string) :: smear
    ! Width for smearing
    Type(in_param)  :: width_smear
    ! Mixing 
    Type(in_string) :: mixing
    ! SFC steps
    Type(in_integer) :: scf_steps
    ! Energy convergence
    Type(in_param)   :: delta_e
    ! k-point sampling
    Integer(Kind=wi)       :: total_kpoints 
    Type(in_integer_array) :: kpoints
    ! basis set
    Type(in_logic)   :: basis_info 
    Type(type_basis), Allocatable :: basis_set(:)
    ! pseudo-potentials
    Type(in_logic)   :: pp_info 
    Type(type_pseudo), Allocatable :: pseudo_pot(:)
    ! Maximum l_orbital
    Type(in_integer) :: max_l_orbital
    ! Total magnetization 
    Type(in_param)   :: total_magnetization 
    ! magnetization
    Type(in_logic)   :: mag_info 
    Type(type_mag),  Allocatable :: magnetization(:)
    ! hubbard
    Type(in_logic)   :: hubbard_info
    Logical          :: hubbard_all_U_zero
    Type(type_hubbard),  Allocatable :: hubbard(:)
    ! Orbital Transformation (OT), only valid for CP2K
    Type(in_logic)   :: ot
    ! Ensemble DFT (EDFT), only valid for CASTEP and ONETEP
    Type(in_logic)   :: edft 
    ! Bands paralellization, only for VASP
    Type(in_integer) :: npar
    ! kpoints paralellization, only for VASP
    Type(in_integer) :: kpar
    ! Bands
    Type(in_integer) :: bands
    ! NGWF, compulsory only for ONETEP
    Type(in_logic)   :: ngwf_info
    Type(type_ngwf),  Allocatable  :: ngwf(:)
    ! PAW for onetep
    Logical :: onetep_paw
    
    ! GC-DFT
    Type(type_gcdft) :: gc

  End Type

  ! Type for motion settings
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  Type :: motion_type
    ! Flag to ensure motion block is not defined more than once 
    Logical  :: generate=.False.
    ! Relaxation method 
    Type(in_string) :: relax_method
    ! force convergence
    Type(in_param_array) :: delta_f
    ! Time step
    Type(in_param) :: timestep
    ! Number of ionic step, either for relaxation or MD
    Type(in_integer) :: ion_steps
    ! Change simulation cell volume
    Type(in_logic)  :: change_cell_volume
    ! Change simulation cell shape 
    Type(in_logic)  :: change_cell_shape
    ! Ensemble 
    Type(in_string) :: ensemble
    ! Temperature 
    Type(in_param)  :: temperature 
    ! Pressure 
    Type(in_param)  :: pressure 
    ! Thermostat 
    Type(in_string) :: thermostat
    ! Thermostat relaxation time 
    Type(in_param)  :: relax_time_thermostat 
    ! Barostat 
    Type(in_string) :: barostat
    ! Barostat relaxation time 
    Type(in_param)  :: relax_time_barostat 
    ! Masses 
    Type(in_logic)   :: mass_info
    Type(species_list), Allocatable :: mass(:) 
  End Type

  ! Type for solvation 
  Type :: solvation_onetep
    Type(in_logic)    :: info
    Type(in_logic)    :: in_vacuum_first 
    Type(in_string)   :: cavity_model
    Type(in_string)   :: dielectric_function
    Type(in_scalar)   :: density_threshold
    Type(in_scalar)   :: density_min_threshold
    Type(in_scalar)   :: density_max_threshold
    Type(in_scalar)   :: beta 
    Type(in_scalar)   :: permittivity_bulk
    Character(Len=256) :: bib_epsilon
    Type(in_logic)    :: soft_radii_info  
    Type(in_scalar)   :: soft_sphere_scale    
    Type(in_scalar)   :: soft_sphere_delta    
    Type(species_list), Allocatable :: soft_radii(:)
    !Apolar terms
    Type(in_string)   :: apolar_terms
    Type(in_string)   :: sasa_definition
    Type(in_scalar)   :: apolar_scaling 
    Type(in_param)    :: solvent_pressure      
    Type(in_param_array) :: surf_tension 
    Type(in_param)    :: smear_ion_width 
  End Type solvation_onetep

  ! Type for Poisson-Boltzmann 
  Type :: electrolyte
    Type(in_logic)    :: info
    Type(in_string)   :: solver
    Type(in_string)   :: neutral_scheme 
    Type(in_string)   :: steric_potential 
    Type(in_param)    :: boltzmann_temp
    Type(in_param)    :: steric_isodensity 
    Type(in_param)    :: steric_smearing
    Type(in_param)    :: capping 
    ! Boltzman ions
    Logical           :: set_necs_shift
    Type(in_logic)    :: boltzmann_ions_info
    Integer           :: number_boltzmann_ions
    Type(in_logic)    :: solvent_radii_info
    Type(species_list), Allocatable :: solvent_radii(:)
    Type(boltzmann_ions) :: boltzmann_ions(max_number_boltzmann_ions)
  End Type electrolyte

  ! Type for multi-grid 
  Type :: multigrid
    Type(in_logic)    :: info
  End Type multigrid

  ! Type for the modelling related variables 
  Type, Public :: simul_type
    Private
    ! General
    !!!!!!!!!
    ! Flag to generate simulation files
    Logical, Public  :: generate=.False.
    ! Details for the components
    Type(component_in_block), Public :: component(max_components)
    !number of total tags
    Integer(Kind=wi),   Public  :: total_tags
    ! Code 
    Character(Len=256), Public  :: code_format
    ! Code version 
    Character(Len=256), Public  :: code_version
    ! physical process 
    Character(Len=256), Public  :: process
    ! vector normal to the surface 
    Character(Len=256), Public  :: normal_vector
    ! Simulation cell
    Real(Kind=wp),      Public  :: cell(3,3)
    ! Length of cell vectors
    Real(Kind=wp),      Public  :: cell_length(3)
    ! Large cell
    Logical,            Public  :: large_cell

    ! Specific variables 
    !!!!!!!!!!!!!!!!!!!!
    ! Type of the simulation to be performed
    Type(in_string), Public :: simulation     
    ! Level of theory 
    Type(in_string), Public :: theory_level
    ! Total charge
    Type(in_scalar), Public :: net_charge
    ! DFT directives
    Type(dft_type),    Public :: dft 
    ! Ions related variables
    Type(motion_type), Public :: motion
    ! Extra directives
    Type(in_logic),    Public :: extra_info
    Type(type_extra),  Public :: extra_directives
    ! Set directives
    Type(type_extra),  Public :: set_directives

    ! Solvation
    Type(solvation_onetep), Public :: solvation
    ! Poisson-Boltzmann 
    Type(electrolyte), Public :: electrolyte 

  Contains
    Private
    Procedure, Public  :: init_input_dft_variables    =>  allocate_input_dft_variables
    Procedure, Public  :: init_input_motion_variables =>  allocate_input_motion_variables
    Procedure, Public  :: init_input_solvation_variables =>  allocate_input_solvation_variables
    Procedure, Public  :: init_input_electrolyte_variables =>  allocate_input_electrolyte_variables
    Final              :: cleanup

  End Type simul_type

Contains

  Subroutine allocate_input_dft_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential DFT variables to build input files for simulations 
    !
    ! author    - i.scivetti April-June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(simul_type), Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(6)
    Character(Len=256)  :: message

    Allocate(T%dft%pseudo_pot(T%total_tags),     Stat=fail(1))
    Allocate(T%dft%magnetization(T%total_tags),  Stat=fail(2))
    Allocate(T%dft%hubbard(T%total_tags),        Stat=fail(3))
    Allocate(T%dft%basis_set(T%total_tags),      Stat=fail(4))
    Allocate(T%dft%kpoints%value(3),             Stat=fail(5))
    Allocate(T%dft%ngwf(T%total_tags),           Stat=fail(6))

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of "DFT" variables to build input files&
                               & for simulations'
      Call error_stop(message)
    End If

    !Set to False just in case
    T%dft%basis_info%stat=.False. 
    T%dft%pp_info%stat=.False. 
    T%dft%mag_info%stat=.False.
    T%dft%hubbard_info%stat=.False.
    T%dft%ngwf_info%stat=.False.

  End Subroutine allocate_input_dft_variables

  Subroutine allocate_input_motion_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential motion variables to build input files for simulations 
    !
    ! author    - i.scivetti April-June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(simul_type), Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(3)
    Character(Len=256)  :: message

    Allocate(T%motion%delta_f%value(1),               Stat=fail(1))
    Allocate(T%motion%delta_f%units(2),               Stat=fail(2))
    Allocate(T%motion%mass(T%total_tags),             Stat=fail(3))  

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of "motion" variables to build input files&
                               & for simulations'
      Call error_stop(message)
    End If

    !Set to False just in case
    T%motion%mass_info%stat=.False.

  End Subroutine allocate_input_motion_variables

  Subroutine allocate_input_solvation_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential solvation input variable 
    !
    ! author    - i.scivetti August 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(simul_type), Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(3)
    Character(Len=256)  :: message

    Allocate(T%solvation%surf_tension%value(1),    Stat=fail(1))
    Allocate(T%solvation%surf_tension%units(2),    Stat=fail(2))
    Allocate(T%solvation%soft_radii(T%total_tags), Stat=fail(3))  

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of "solvation" variables to build input files&
                               & for simulations'
      Call error_stop(message)
    End If

    !Set to False just in case
    T%solvation%soft_radii_info%stat=.False.

  End Subroutine allocate_input_solvation_variables

  Subroutine allocate_input_electrolyte_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential elctrolyte input variable 
    !
    ! author    - i.scivetti September 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(simul_type), Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(1)
    Character(Len=256)  :: message

    Allocate(T%electrolyte%solvent_radii(T%total_tags), Stat=fail(1))  

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of "electrolyte" variables to build input files&
                               & for simulations'
      Call error_stop(message)
    End If

    !Set to False just in case
    T%electrolyte%solvent_radii_info%stat=.False.

  End Subroutine allocate_input_electrolyte_variables

  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti April 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type) :: T

    If (Allocated(T%motion%delta_f%units)) Then
      Deallocate(T%motion%delta_f%units)
    End If

    If (Allocated(T%motion%delta_f%value)) Then
      Deallocate(T%motion%delta_f%value)
    End If

    If (Allocated(T%motion%mass)) Then
      Deallocate(T%motion%mass)
    End If


    If (Allocated(T%dft%kpoints%value)) Then
      Deallocate(T%dft%kpoints%value)
    End If

    If (Allocated(T%dft%pseudo_pot)) Then
      Deallocate(T%dft%pseudo_pot)
    End If

    If (Allocated(T%dft%hubbard)) Then
      Deallocate(T%dft%hubbard)
    End If

    If (Allocated(T%dft%basis_set)) Then
      Deallocate(T%dft%basis_set)
    End If

    If (Allocated(T%dft%ngwf)) Then
      Deallocate(T%dft%ngwf)
    End If

    If (Allocated(T%solvation%soft_radii)) Then
      Deallocate(T%solvation%soft_radii)
    End If

  End Subroutine cleanup

End Module simulation_setup
