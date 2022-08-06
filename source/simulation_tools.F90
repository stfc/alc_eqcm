!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that defines simulation type and required subroutines to
! process data and build input files for simulation. This module supports
! simulation_setup.  
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author    - i.scivetti Oct 2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module simulation_tools

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
  Integer(Kind=wi), Parameter, Public  :: max_directives=200 
  ! Error in the initialization of magnetization                            
  Real(Kind=wp), Parameter, Public :: error_mag = 0.00001_wp

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
    Character(Len=256) :: basis
    Character(Len=256) :: potential 
    Character(Len=8)   :: tag
    Character(Len=2)   :: element
  End Type type_pseudo

  ! Type for extra directives 
  Type, Public :: type_extra 
    Character(Len=256) :: array(max_directives)
    Integer(Kind=wi)   :: N0
  End Type type_extra

  ! Type for basis_set
  Type :: type_basis
    Character(Len=8)   :: tag
    Character(Len=2)   :: element
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

  ! Type for masses
  Type :: type_mass
    Character(Len=8)  :: tag
    Character(Len=2)  :: element
    Real(Kind=wp)     :: value
  End Type type_mass

 
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
    Type(type_mass), Allocatable :: mass(:) 
  End Type

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


  Contains
    Private
    Procedure, Public  :: init_input_dft_variables    =>  allocate_input_dft_variables
    Procedure, Public  :: init_input_motion_variables =>  allocate_input_motion_variables
    Final              :: cleanup

  End Type simul_type

  Public :: check_extra_directives, record_directive, scan_extra_directive, check_initial_magnetization
  Public :: print_warnings
  
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

  End Subroutine allocate_input_motion_variables


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

  End Subroutine cleanup

  
  Subroutine check_extra_directives(sentence, symbol, code, ic)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check the structure of sub-block &extra_directives 
    !
    ! author    - i.scivetti Jul 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=256), Intent(In   ) :: sentence
    Character(Len=*),   Intent(In   ) :: symbol 
    Character(Len=*),   Intent(In   ) :: code
    Integer(Kind=wi),   Intent(In   ) :: ic

    Character(Len=256)  :: messages(4)
    Logical             :: error

    error=.False.

    Write (messages(1),'(1x,a,i3)') '***ERROR: in block &extra_directives: user-defined directive ', ic
    Write (messages(2),'(1x,a)')     Trim(Adjustl(sentence))
    Write (messages(4),'(1x,a)')    'Please check. Sentences starting with "#" are assumed as comments.'

    If (Index(Trim(Adjustl(sentence)), '%') /= 1) Then
      If (Index(Trim(Adjustl(sentence)), '#') > 1) Then
        If (Index(Trim(Adjustl(sentence)), Trim(symbol)) == 0 ) Then
          Write (messages(3),'(1x,5a)')    'does not contain the symbol "', Trim(symbol), '", which is needed to specify&
                                       & the directive according to the ', Trim(code), ' format.'
          error=.True.
        Else If (Index(Trim(Adjustl(sentence)), Trim(symbol)) == 1 ) Then
          Write (messages(3),'(1x,3a)')    'contains the symbol "', Trim(symbol), '" at the beginning of the declaration.'
          error=.True.
        End If
  
        If (Index(sentence, "=") > Index(sentence, '#')) Then
           Write (messages(3),'(1x,3a)')    'contains the symbol "', Trim(symbol),&
                                          & '" but it is wrongly used to define the directive.'
           error=.True.
        End If
  
      Else If (Index(Trim(Adjustl(sentence)), '#') == 0) Then
  
        If (Index(Trim(Adjustl(sentence)), Trim(symbol)) == 0 ) Then
          Write (messages(3),'(1x,5a)')    'does not contain the symbol "', Trim(symbol), '", which is needed to specify&
                                       & the directive according to ', Trim(code), ' format.'
          error=.True.
        Else If (Index(Trim(Adjustl(sentence)), Trim(symbol)) == 1 ) Then
          Write (messages(3),'(1x,a)')    'contains the symbol "', Trim(symbol), '" at the beginning of the declaration.'
          error=.True.
        End If
  
      End If
  
      If (error) Then
        Call info(messages, 4)
        Call error_stop(' ')
      End If
    Else
      If (Trim(code) /= 'VASP') Then
        Write (messages(3),'(2(1x,a))')    'Definition of blocks are not allowed in format', Trim(code)
      Else
        Write (messages(3),'(2(1x,a))')    'Sentences starting with "%" are not allowed in format', Trim(code)
      End If
      Call info(messages, 3)
      Call error_stop(' ')
    End If

  End Subroutine check_extra_directives  

  Subroutine record_directive(iunit, message, tag, name_dir, ic) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print and keep record of the directives for the
    ! atomisitc simulations 
    ! 
    ! author    - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit 
    Character(Len=*), Intent(In   ) :: message  
    Character(Len=*), Intent(In   ) :: tag      
    Character(Len=*), Intent(  Out) :: name_dir
    Integer(Kind=wi), Intent(InOut) :: ic

    Write(iunit, '(a)') Trim(message)
    name_dir= Trim(tag)
    ic=ic+1

  End Subroutine record_directive

  Subroutine scan_extra_directive(sentence, set_directives, symbol, found, tag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to scan the directives defined in block &extra_directives against
    ! the directives set from the definitions of &block_simulation_settings 
    ! 
    ! author    - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: sentence 
    Type(type_extra), Intent(In   ) :: set_directives
    Character(Len=*), Intent(In   ) :: symbol
    Logical,          Intent(InOut) :: found
    Character(Len=*), Intent(  Out) :: tag
   
    Character(Len=256) :: word, word2, dir
    Integer(Kind=wi)   :: j

    word=Trim(Adjustl(sentence))
    Call remove_symbols(word,Trim(symbol))        
    Read(word,Fmt=*) tag
    Call capital_to_lower_case(word)
    Read(word,Fmt=*) dir
    j=1
    Do While (j <= set_directives%N0 .And. (.Not. found))
      word2=set_directives%array(j)
      Call capital_to_lower_case(word2)
      If (Trim(dir) == Trim(word2)) Then 
        found=.True.
      End If
      j=j+1
    End Do

  End Subroutine scan_extra_directive
  
  Subroutine check_initial_magnetization(net_elements, list_tag, N0, mag_ini, target_mag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the initial magnetic moments of the resulting models. 
    ! This subroutine is only invoked if there are differences between the assigned
    ! total magnetization and the initial magnetization of the generated model
    ! 
    ! author    - i.scivetti May-July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: net_elements
    Character(Len=8),  Intent(In   ) :: list_tag(max_components) 
    Integer(Kind=wi),  Intent(In   ) :: N0(net_elements)
    Real(Kind=wp),     Intent(In   ) :: mag_ini(max_components)
    Real(Kind=wp),     Intent(In   ) :: target_mag

    Real(Kind=wp)      :: tot_mag
    Character(Len=256) :: message
    Character(Len=256) :: messages(3)
    Integer(Kind=wi) :: i
  
    tot_mag=0.0_wp
    Do i=1, net_elements
      tot_mag=tot_mag+mag_ini(i)*N0(i)
    End Do
 
    If (Abs(tot_mag-target_mag) > error_mag) Then
      Call info(' ', 1)
      Call info(' Summary of the total amount and initial magnetic moment of species', 1)
      Call info(' -------------------------------------------------', 1)
      Write (message, '(1x,a,2(6x,a))') 'Tag', 'Amount', 'Initial magnetic moment/spin'
      Call info(message, 1)
      Call info(' -------------------------------------------------', 1)
      Do i=1, net_elements
        Write (message, '(1x,a,2x,i5,f10.2)') list_tag(i), N0(i), mag_ini(i) 
        Call info(message, 1)    
      End Do
      Call info(' -------------------------------------------------', 1)
      Write (messages(1),'(1x,a,f10.2)') 'Total initial magnetic moment/spin: ', tot_mag
      Write (messages(2),'(1x,a,f10.2,a)') 'Targeted magnetic moment/spin:      ',&
                                     & target_mag,&
                                     & ' (from the value of "total_magnetization")' 
      Write (messages(3),'(1x,a)') '***ERROR: "total_magnetization" must be&
                                   & equal to the total initial magnetization (or viceversa). Please change the&
                                   & value of "total_magnetization" (or values in sub-block &magnetization).&
                                   & If problems persist, remove "total_magnetization".'
      Call info(messages, 3)
      Call error_stop(' ')
    End If

  End Subroutine check_initial_magnetization

  Subroutine print_warnings(header, print_header, message, dim)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Auxiliary subroutine to print warning headers 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=256), Intent(In   ) :: header
    Logical,            Intent(InOut) :: print_header
    Character(Len=256), Intent(In   ) :: message(*) 
    Integer(Kind=wi),   Intent(In   ) :: dim     

    If (print_header) Then
      Call info(header,1)
      print_header=.False.
    End If
    Call info(message,dim)

  End Subroutine print_warnings
 
  
End Module simulation_tools
