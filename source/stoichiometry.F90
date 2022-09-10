!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Module that finds the stoichiometric coefficients from the charge and mass balance equations
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author        - i.scivetti Nov   2020
! Contribution  - i.scivetti March 2021 - Fixes to allow the generation of atomistic models
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module stoichiometry
  
  Use constants,        Only : Farad, &
                               chemsymbol, &
                               g_to_ng, &
                               NPTE,&
                               max_components
  Use electrode,        Only : electrode_type
  Use eqcm,             Only : eqcm_type

  Use fileset,          Only : file_type, &
                               FILE_CHARACT,&
                               FILE_ELECTRO, &
                               FILE_INTERCALATION_OX,&
                               FILE_INTERCALATION_RED,& 
                               FILE_SET_EQCM,&
                               FOLDER_ANALYSIS, &
                               refresh_out_eqcm

  Use input_types,      Only : in_integer, &
                               in_logic,  & 
                               in_string, &
                               in_scalar
  Use numprec,          Only : wi, &
                               wp
  Use process_data,     Only : check_for_symbols, &
                               get_word_length, &
                               remove_symbols
  Use redox,            Only : redox_type
  Use unit_output,      Only : error_stop,&
                               info

  Implicit None
  Private
   
  ! Default efficiency
  Real(Kind=wp), Public              ::  efficiency_default = 1.0_wp
  ! Default discretization of the stoichiometric space
  Real(Kind=wp), Public              ::  ds_default = 0.01_wp
  ! Paramter for the maximum length of the species name
  Integer(Kind=wi), Parameter :: max_length_name_species = 12

  ! Type for the compoentents
  Type :: components_type
    Character(Len=8)  :: tag(max_components) 
    Character(Len=2)  :: element(max_components)
    Integer(Kind=wi)  :: atomic_number(max_components)
    Integer(Kind=wi)  :: N0(max_components)
    Logical           :: fread= .False.
  End Type  

  Type :: species_stoich
    !Description
    Real(Kind=wp)     :: bond_cutoff
    Real(Kind=wp)     :: s0_pristine
    Real(Kind=wp)     :: s0
    Real(Kind=wp)     :: s=0.0_wp
    Real(Kind=wp)     :: mass
    Real(Kind=wp)     :: ox
    Character(Len=32) :: tag    
    Character(Len=32) :: vartype0
    Character(Len=32) :: vartype
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Integer(Kind=wi)  :: num_components
    Integer(Kind=wi)  :: atoms_per_species 
    Character(Len=32) :: topology
    Type(components_type) :: component
    !Constraints-related
    Real(Kind=wp)  :: minvalue 
    Real(Kind=wp)  :: maxvalue 
    Real(Kind=wp)  :: range(2)
    Real(Kind=wp)  :: limit(2)
    Logical        :: fmin  
    Logical        :: fmax  
    Logical        :: frange  
  End Type


  Type :: const_unit
    Character(Len=32)  :: type 
    Character(Len=32)  :: leg
    Character(Len=32)  :: tag_species 
    Character(Len=32)  :: tag0_species 
    Character(Len=32)  :: species(2) 
    Real(Kind=wp)      :: value(2) 
    Logical            :: fread= .False.
    Logical            :: fail = .False.
  End Type


  ! Variables for Sotichiometric analysis
  Type, Public :: stoich_type
  Private
    ! Logical for one or more solutions
    Logical                            :: single_sol
    ! Logical to check if solution is found
    Logical                            :: sol_exist =.True.
    ! molar mass
    Real(Kind=wp)                      :: molar_mass=0.0_wp
    ! number of mole of the host material for intercalation
    Real(Kind=wp), Allocatable         ::  mole(:)
    ! area ratio for oxidation/reduction
    Real(Kind=wp), Allocatable         ::  elec_depo_ratio_ox(:)
    Real(Kind=wp), Allocatable         ::  elec_depo_ratio_red(:)
    ! number of mole of material to being electrodeposited
    Real(Kind=wp), Allocatable         ::  mole_ox(:)
    Real(Kind=wp), Allocatable         ::  mole_red(:)
    ! Number of collected electrons
    Real(Kind=wp), Public              ::  electrons
    ! Discretization
    Type(in_scalar), Public            ::  discretization
    ! fragment of the CV cycle (oxidation or reduction)
    Character(Len=64)                  ::  cv_leg
    Character(Len=64)                  ::  first_leg
    ! Stoichiometry cycle
    Integer(Kind=wi)                   ::  cv_cycle
    ! Stoichiometric flags for CV cycling
    Logical                            ::  continue_cv=.True.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Variables for species 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Indexes for dependent variables   
    Integer(Kind=wi),              Public       :: index_dep(2)   
    ! Indexes, list, counters and boundaries for independent variables   
    Integer(Kind=wi), Allocatable, Public       :: index_indep(:)   
    Integer(Kind=wi), Allocatable               :: list(:)
    Integer(Kind=wi), Allocatable               :: maxlist(:)
    Integer(Kind=wi), Allocatable               :: list_increment(:)
    ! matrices for charge and mass balance
    Real(Kind=wp), Allocatable, Public  ::  matrix_balance0(:,:) 
    Real(Kind=wp), Allocatable, Public  ::  matrix_balance(:,:) 
    ! 2x1 matrix with integrated redox data
    Real(Kind=wp)                       ::  matrix_redox(2) 
    ! 2x2 matrix for calculation of dependent variables 
    Real(Kind=wp)                       ::  matrix_eq(2,2) 
    ! 2x2 matrix for calculation of dependent variables 
    Real(Kind=wp)                       ::  inv_matrix_eq(2,2) 
    ! Total number of stoichiometric solutions
    Integer(Kind=wi), Allocatable, Public   :: numsol(:,:)
    ! Number of species involved in the reaction (fixed and variable)
    Type(in_integer) ,   Public    :: num_species
    ! Number of variable stoichiometry coefficients
    Integer(Kind=wi),             Public    :: num_variables
    ! Number of variable independent stoichiometry coefficients
    Integer(Kind=wi),             Public    :: num_indep
    ! Number of variable independent stoichiometry coefficients
    Integer(Kind=wi),             Public    :: num_dep
    ! Stoichiometric characterization of the species
    Type(species_stoich), Allocatable, Public :: species(:)
    ! Number of different tags used to name the involved species 
    Integer(Kind=wi),             Public    :: total_tags
    ! Stoichiometric solutions as a function of cycle
    Real(Kind=wp), Allocatable, Public  :: solution_coeff(:,:,:)
    ! Number of moles as a function of cycle
    Real(Kind=wp), Allocatable, Public  :: solution_moles(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Variables for constraints
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Directive for beginning of block_species
    Type(in_string), Public        :: block_species
    ! Directive for beginning of block_constraints
    Type(in_string), Public        :: block_constraints
    ! Directive for beginning of block_species_components
    Type(in_string), Public        :: block_species_components
    ! Number of species involved in the reaction (fixed and variable)
    Type(in_integer) ,   Public    :: num_constraints
    ! Stoichiometric characterization of the species
    Type(const_unit), Allocatable, Public :: constraints(:)
    ! Matrix to introduce contraints
    Real(Kind=wp), Allocatable, Public  ::  matrix_constraints(:,:) 
    ! Indexes for constrainted variables   
    Integer(Kind=wi), Allocatable, Public        :: index_const(:)
    ! Number of net constrainted stoichiometry coefficients
    Integer(Kind=wi), Public                     :: net_const
    ! Number of constraints per fragment
    Integer(Kind=wi)                        :: nconst_leg
  Contains
    Private
    Procedure, Public :: init_species_arrays     => allocate_species_arrays
    Procedure, Public :: init_constraints_arrays => allocate_constraints_arrays
    Procedure         :: species_extra           => allocate_species_extra_arrays
    Final             :: cleanup
  End Type stoich_type

  Public :: stoichometry_analysis, check_stoich_settings 
  Public :: check_constraints_settings
  Public :: check_components_species

Contains

  Subroutine allocate_constraints_arrays(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays for stoichiometry analysis 
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stoich_type),  Intent(InOut)  :: T

    Integer(Kind=wi) :: fail(3)
    Integer(Kind=wi) :: number_constraints
    Integer(Kind=wi) :: number_species
    Character(Len=256)   :: message

    number_constraints= T%num_constraints%value
    number_species    = T%num_species%value

    fail =0

    Allocate(T%constraints(number_constraints),                   Stat=fail(1))
    Allocate(T%matrix_constraints(number_species,number_species), Stat=fail(2))
    Allocate(T%index_const(number_constraints),                   Stat=fail(3))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for constraints (stoichimetric) arrays'
      Call error_stop(message)
    End If

    T%constraints(:)%value(1)=0.0_wp

  End Subroutine allocate_constraints_arrays


  Subroutine allocate_species_arrays(T, what)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays for stoichiometry analysis 
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stoich_type), Intent(InOut)  :: T
    Character(Len=*),   Intent(In   )  :: what

    Integer(Kind=wi) :: fail(6)
    Integer(Kind=wi) :: number_species
    Character(Len=256)   :: message

    number_species=T%num_species%value
    fail =0

    If (what=='block') Then
      Allocate(T%species(number_species)          , Stat=fail(1))
      If (Any(fail > 0)) Then
        Write (message,'(1x,1a)') '***ERROR: Allocation problems for stoichiometric species arrays'
        Call error_stop(message)
      End If
    ElseIf (what=='matrices') Then
      Allocate(T%index_indep(T%num_indep),          Stat=fail(1))
      Allocate(T%list(T%num_indep),                 Stat=fail(2))
      Allocate(T%maxlist(T%num_indep),              Stat=fail(3))
      Allocate(T%matrix_balance0(2,number_species), Stat=fail(4))
      Allocate(T%matrix_balance(2,number_species),  Stat=fail(5))
      Allocate(T%list_increment(T%num_indep),       Stat=fail(6))
      If (Any(fail > 0)) Then
        Write (message,'(1x,1a)') '***ERROR: Allocation problems for stoichiometric working matrices'
        Call error_stop(message)
      End If
    End If

  End Subroutine allocate_species_arrays

  Subroutine allocate_species_extra_arrays(T, process, cycles)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays for stoichiometric mole
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(stoich_type), Intent(InOut) :: T
    Type(in_string),    Intent(In   ) :: process  
    Integer(Kind=wi),   Intent(In   ) :: cycles  
 
    Integer(Kind=wi)   :: fail(6)
    Integer(Kind=wi)   :: number_species
    Character(Len=256) :: message
 
    fail=0
    number_species=T%num_species%value
 
    If (Trim(process%type) == 'intercalation') Then
      Allocate(T%mole(1)          , Stat=fail(1))
      Allocate(T%numsol(2,cycles) , Stat=fail(2))
      Allocate(T%solution_coeff(number_species,2,cycles),  Stat=fail(3))
      If (Any(fail > 0)) Then
        Write (message,'(1x,1a)') '***ERROR: Allocation problems for extra stoichometric arrays (intercalation)'
        Call error_stop(message)
      End If
      T%numsol=1
    Else If (Trim(process%type) == 'electrodeposition') Then
      Allocate(T%mole_ox(cycles),             Stat=fail(1))
      Allocate(T%mole_red(cycles),            Stat=fail(2))
      Allocate(T%elec_depo_ratio_ox(cycles),  Stat=fail(3))
      Allocate(T%elec_depo_ratio_red(cycles), Stat=fail(4))
      Allocate(T%solution_moles(2,cycles),    Stat=fail(5))
      Allocate(T%numsol(2,cycles),            Stat=fail(6))
      If (Any(fail > 0)) Then
        Write (message,'(1x,1a)') '***ERROR: Allocation problems for extra stoichometric arrays (electrodeposition)'
        Call error_stop(message)
      End If
    End If
 
  End Subroutine allocate_species_extra_arrays        


  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type) :: T

    If (Allocated(T%species)) Then
      Deallocate(T%species)
    End If

    If (Allocated(T%matrix_balance)) Then
      Deallocate(T%matrix_balance)
    End If

    If (Allocated(T%matrix_balance0)) Then
      Deallocate(T%matrix_balance0)
    End If

    If (Allocated(T%index_indep)) Then
      Deallocate(T%index_indep) 
    End If

    If (Allocated(T%list)) Then
      Deallocate(T%list) 
    End If

    If (Allocated(T%list_increment)) Then
      Deallocate(T%list_increment) 
    End If

    If (Allocated(T%maxlist)) Then
      Deallocate(T%maxlist) 
    End If

    If (Allocated(T%mole)) Then
      Deallocate(T%mole) 
    End If

    If (Allocated(T%numsol)) Then
      Deallocate(T%numsol) 
    End If
    
    If (Allocated(T%solution_coeff)) Then
      Deallocate(T%solution_coeff)
    End If

    If (Allocated(T%solution_moles)) Then
      Deallocate(T%solution_moles)
    End If 

    If (Allocated(T%mole_ox)) Then
      Deallocate(T%mole_ox) 
    End If

    If (Allocated(T%mole_red)) Then
      Deallocate(T%mole_red) 
    End If

    If (Allocated(T%elec_depo_ratio_ox)) Then
      Deallocate(T%elec_depo_ratio_ox) 
    End If

    If (Allocated(T%elec_depo_ratio_red)) Then
      Deallocate(T%elec_depo_ratio_red) 
    End If

    If (Allocated(T%constraints)) Then
      Deallocate(T%constraints) 
    End If

    If (Allocated(T%matrix_constraints)) Then
      Deallocate(T%matrix_constraints) 
    End If

    If (Allocated(T%index_const)) Then
      Deallocate(T%index_const) 
    End If

  End Subroutine cleanup


  Subroutine stoichometry_analysis(files, eqcm_data, electrode_data, redox_data, stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Perform the stoichiometric analysis 
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut) :: files(:)
    Type(eqcm_type),      Intent(InOut) :: eqcm_data
    Type(electrode_type), Intent(InOut) :: electrode_data
    Type(redox_type),     Intent(InOut) :: redox_data
    Type(stoich_type),    Intent(InOut) :: stoich_data

    Real(Kind=wp)       :: DM, DQ
    Character(Len=256)  :: messages(4), message_end(2)
    Character(Len=256)  :: iunit_name(2), error_msg
    Integer(Kind=wi)    :: iunit_stoich(2)
    Integer(Kind=wi)    :: i, j, k
    Logical             :: fprint, ferror
    Logical             :: fsol, fsol_dep
    Integer(Kind=wi)    :: legs, redox_cycles 

    fsol=.True.
    fsol_dep=.True. 
    fprint=.True.
    ferror=.False.

    Call info(' ',1)
    Call stoich_data%species_extra(eqcm_data%process, redox_data%limit_cycles)
        
    Write (messages(1),'(1x,a)') 'Stoichiometry analysis'
    Write (messages(2),'(1x,a)') '======================'
    Call info(messages, 2)

    If (Trim(eqcm_data%process%type) == 'intercalation') Then
      ! Determine the number of moles for the host material to participate in the reaction
      Call set_mole(eqcm_data, electrode_data, stoich_data)
      If (redox_data%label_leg(1,1)=='reduction') Then
        Open(Newunit=files(FILE_INTERCALATION_RED)%unit_no, File=files(FILE_INTERCALATION_RED)%filename, Status='Replace')
        iunit_stoich(1)=files(FILE_INTERCALATION_RED)%unit_no
        iunit_name(1)=Trim(files(FILE_INTERCALATION_RED)%filename) 
      ElseIf (redox_data%label_leg(1,1)=='oxidation') Then
        Open(Newunit=files(FILE_INTERCALATION_OX)%unit_no, File=files(FILE_INTERCALATION_OX)%filename, Status='Replace')
        iunit_stoich(1)=files(FILE_INTERCALATION_OX)%unit_no
        iunit_name(1)=Trim(files(FILE_INTERCALATION_OX)%filename) 
      End If
    End If

    Call print_stoichiometric_settings(eqcm_data, stoich_data)

    If (stoich_data%num_variables >= 2 .And. Trim(eqcm_data%process%type) == 'intercalation') Then
      Call set_matrix_balance(stoich_data, 'initial')
    End If 

    ! Define limits
    redox_cycles=redox_data%limit_cycles
    legs=2 

    ! Find if the first leg is oxidation or reduction    
    stoich_data%first_leg = Trim(redox_data%label_leg(1,1))
     
    i=1
    ! Here is where the loop starts
    Do While (i <= redox_cycles .And. stoich_data%continue_cv .And. stoich_data%sol_exist)
      stoich_data%cv_cycle=i
      j=1
      Do While (j <= legs .And. stoich_data%continue_cv .And. stoich_data%sol_exist)
        If (redox_data%label_leg(j,i)=='oxidation') Then
          DM=redox_data%DM_ox(i)
          DQ=redox_data%DQ_ox(i)
        Else If (redox_data%label_leg(j,i)=='reduction') Then
          DM=redox_data%DM_red(i)
          DQ=redox_data%DQ_red(i)
        End If
        stoich_data%cv_leg = Trim(redox_data%label_leg(j,i))
        ! If intercalation....
        If (Trim(eqcm_data%process%type) == 'intercalation') Then
          If (fprint) Then
            If (j == 2) Then
              If (redox_data%label_leg(2,1)=='reduction') Then
                Open(Newunit=files(FILE_INTERCALATION_RED)%unit_no, File=files(FILE_INTERCALATION_RED)%filename, Status='Replace')
                iunit_stoich(2)=files(FILE_INTERCALATION_RED)%unit_no
                iunit_name(2)=Trim(files(FILE_INTERCALATION_RED)%filename) 
              ElseIf (redox_data%label_leg(2,1)=='oxidation') Then
                Open(Newunit=files(FILE_INTERCALATION_OX)%unit_no, File=files(FILE_INTERCALATION_OX)%filename, Status='Replace')
                iunit_stoich(2)=files(FILE_INTERCALATION_OX)%unit_no
                iunit_name(2)=Trim(files(FILE_INTERCALATION_OX)%filename) 
              End If
              Write (message_end(2),'(1x,4a)') 'Print stoichiometric solutions to files ',&
                                               & Trim(FOLDER_ANALYSIS)//'/'//Trim(iunit_name(1)), ' and ',&
                                               & Trim(FOLDER_ANALYSIS)//'/'//Trim(iunit_name(2))
              fprint=.False.
            Else
              Write (message_end(1),'(1x,2a)') 'Print stoichiometric solutions to file ',&
                                              & Trim(FOLDER_ANALYSIS)//'/'//Trim(iunit_name(1))
            End If
          End If
          !Set matrix DM-DQ 
          !!!!!!!!!!!!!!!!! 
          stoich_data%matrix_redox(1)=DM/stoich_data%mole(1)
          stoich_data%matrix_redox(2)=DQ/Farad/stoich_data%mole(1)
          !!!!!!!!!!!!!!!!! 
          If (stoich_data%block_constraints%fread) Then
            If (i==1 .And. j==1) Then 
              Call info('CONSTRAINTS (from &block_constraints_species)',1)  
            End If
          End If 
          !Find which variables are dependent, independent and fixed, including contraints
          Call classify_variables(stoich_data)
          ! Build the charge and mass balance matrix
          If (stoich_data%num_variables >= 2) Then
            Call set_matrix_balance(stoich_data, 'cv')
          End If 
          ! Find the domains for the involved stoichiometric variables
          Call set_stoich_domains(stoich_data)
          !Print headings to files INTERCALATION_X
          If (i==1) Then
            Call print_stoich_file(stoich_data, iunit_stoich(j), 'header')
          End If
          If (.Not. stoich_data%single_sol) Call print_stoich_file(stoich_data, iunit_stoich(j), 'cycle_only')
          ! Solve the stoichiometry according to the particular case
          If (stoich_data%num_indep >= 1)then   
            Call sample_stoich_phase_space(stoich_data, 'find_solutions', iunit_stoich(j), i, j)
          ElseIf (stoich_data%num_indep == 0)then
            If (stoich_data%num_dep == 2) Then
              Call solve_stoichiometry(stoich_data, fsol_dep)
              If (fsol_dep) Then
                If (stoich_data%block_constraints%fread) Then
                  Call setcheck_constrained_solutions(stoich_data, fsol)
                Else
                  fsol=.True.
                End If
              Else
                fsol=.False.
              End If
            ElseIf (stoich_data%num_dep == 1) Then
              Call solve_stoichiometry(stoich_data, fsol)
            End If
            If (fsol) Then
              Call obtain_electrons(stoich_data) 
              Call print_stoich_file(stoich_data, iunit_stoich(j), 'coefficients')
            Else
              stoich_data%sol_exist=.False.
            End If
            Do k=1, stoich_data%num_species%value
               stoich_data%species(k)%s0= stoich_data%species(k)%s + stoich_data%species(k)%s0
            End Do
          End If
          ! Store stoichiometry solution only if there is a single solution to the charge/mass balance equation
          If (stoich_data%single_sol) Then
            stoich_data%solution_coeff(:,j,i)=stoich_data%species(:)%s0
          End If
        ! If electrodeposition....
        ElseIf (Trim(eqcm_data%process%type) == 'electrodeposition') Then
          !Find which variables are dependent, independent and fixed, including contraints
          Call classify_variables(stoich_data)
          Call electro_deposition(i, j, DM, DQ, stoich_data, eqcm_data%mass_frequency%fread, eqcm_data%mass%fread)
        End If
        j=j+1
      End Do
      i=i+1
    End Do

    If (Trim(eqcm_data%process%type) == 'electrodeposition') Then
      Call print_deposition(stoich_data, files, redox_data%limit_cycles, eqcm_data%mass_frequency%fread, eqcm_data%mass%fread) 
    Else If (Trim(eqcm_data%process%type) == 'intercalation') Then
      If (.Not. stoich_data%sol_exist) Then
         ferror=.True.
         Write (messages(1),'(3a,i3)') '***ERROR:  No solution is found during the ', Trim(stoich_data%cv_leg),&
                                      &' part of cycle ', stoich_data%cv_cycle
         Write (messages(2),'(1x,2a)')  'Please review i) stoichiometric settings in &block_species,& 
                                      & ii) the efficiency set for the reaction and&
                                      & iii) computed quantities in the table above and in file ', &
                                      & Trim(FOLDER_ANALYSIS)//'/'//Trim(files(FILE_CHARACT)%filename)

         If (stoich_data%nconst_leg == 0) Then
           Write (messages(3),'(a)')     ' '
         Else
           Write (messages(3),'(1x,3a)')  '***Constraints seetings for ', Trim(stoich_data%cv_leg), ' should&
                                        & also be reviewed***' 
         End If
        
         Write (messages(4),'(1x,a)') 'IMPORTANT: the reaction might not involve only the&
                                      & species defined in &block_species.'

        If (stoich_data%cv_cycle==1) Then
          If ((j-1) == 1) Then ! First cycle has failed
            Write (messages(1),'(a)') '***ERROR: No stoichiometric solution found for the initial fragment of cycle 1 !!!'
            Call info(messages,4)
            Write (error_msg, '(1x,3a)')  'File '//Trim(FOLDER_ANALYSIS)//'/'//Trim(iunit_name(1))//' does NOT CONTAIN values'
          Else
            Call info(messages,4)
            Call info(message_end(1),1)
            Write (error_msg, '(1x,3a)')  'File '//Trim(FOLDER_ANALYSIS)//'/'//Trim(iunit_name(2))//' does NOT CONTAIN values'
          End If 
        Else 
          Call info(messages,4)
          Write (error_msg, '(a)') Trim(message_end(2))
        End If
      End If 

      If (.Not. ferror) Call info(' ', 1)

      ! Close STOICHIOMETRY  files
      Close(iunit_stoich(1))
      ! Move file
      Call execute_command_line('mv '//Trim(iunit_name(1))//' '//Trim(FOLDER_ANALYSIS)) 
      If (.Not. fprint) Then
        If (.Not. ferror) Call info(message_end(2),1)
        Close(iunit_stoich(2))
        ! Move file
        Call execute_command_line('mv '//Trim(iunit_name(2))//' '//Trim(FOLDER_ANALYSIS)) 
      Else
        If (.Not. ferror) Call info(message_end(1),1)
      End If
 
      If (ferror) Call error_stop(error_msg) 

    End If 
 
    ! refresh OUT_EQCM
    Call refresh_out_eqcm(files)
 
  End Subroutine stoichometry_analysis


  Subroutine set_stoich_domains(stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set the domain of involved stoichiometric species 
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),  Intent(InOut) :: stoich_data

    Integer(Kind=wi) :: i
    
    !Define lower limit 
    Do i=1, stoich_data%num_species%value
      If (stoich_data%species(i)%vartype /= 'fixed') Then
        stoich_data%species(i)%limit(1)=-stoich_data%species(i)%s0
        stoich_data%species(i)%limit(2)=-Huge(1.0_wp)
      End If
    End Do

    !Define upper limit for independent variables
    If (stoich_data%num_indep >= 1)then
      Call sample_stoich_phase_space(stoich_data, 'upper_limit')
    End If 

  End Subroutine set_stoich_domains


  Subroutine sample_stoich_phase_space(stoich_data, operation, iunit, icycle, jleg)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Sample the stoichiometric space to determine the range of solutions
    ! compatible with EQCM data
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),          Intent(InOut) :: stoich_data
    Character(Len=*),           Intent(In   ) :: operation
    Integer(Kind=wi), Optional, Intent(In   ) :: iunit
    Integer(Kind=wi), Optional, Intent(In   ) :: icycle, jleg

    Integer(Kind=wi) :: i, k
    Integer(Kind=wi) :: ix, indx

    Real(Kind=wp) :: maxs0, Delta_stoich
    Integer(Kind=wi) :: isol_found, isol_domain
    Integer(Kind=wi) :: itot
    Character(Len=256)   :: message
    Logical :: fsol_dep, fsol, felect

    fsol=.True. 
    fsol_dep=.True.
    felect=.True.

    
    ! Define variables before looping
    ! Variables for upper limits of independent variables
    If (operation=='upper_limit') Then
      itot=1
      isol_domain=0
      Delta_stoich=Ds_default
      maxs0=-Huge(0.0_wp)

      Do i=1,stoich_data%num_species%value
        If (stoich_data%species(i)%vartype /= 'fixed') Then
          If (stoich_data%species(i)%s0 > maxs0) Then
            maxs0=stoich_data%species(i)%s0
          End If  
        End If
      End Do

      maxs0=4.0

      Do i=1, stoich_data%num_indep
        indx=stoich_data%index_indep(i)
        stoich_data%maxlist(i)=nint((maxs0-stoich_data%species(indx)%limit(1))/Ds_default+1.0_wp)
        If (stoich_data%species(indx)%fmax) Then
          stoich_data%list(i)=stoich_data%maxlist(i)+1
          stoich_data%list_increment(i)=-1
        Else
          stoich_data%list(i)=0
          stoich_data%list_increment(i)=1 
        End If
        itot=itot*stoich_data%maxlist(i)
      End Do


    ! Set mesh details for computation of the sotichiometric space
    ElseIf (operation=='find_solutions') Then

      Delta_stoich=stoich_data%discretization%value
      itot=1
      isol_found=0
      Do i=1, stoich_data%num_indep
        indx=stoich_data%index_indep(i)
        stoich_data%maxlist(i)=nint((stoich_data%species(indx)%limit(2)-stoich_data%species(indx)%limit(1))/Delta_stoich+1.0_wp)
        If (stoich_data%species(indx)%fmax) Then
          stoich_data%species(indx)%limit(1)=stoich_data%species(indx)%limit(2)
          stoich_data%maxlist(i)=1
          stoich_data%list(i)=2
          stoich_data%list_increment(i)=-1
        ElseIf (stoich_data%species(indx)%fmin) Then
          stoich_data%species(indx)%limit(1)=stoich_data%species(indx)%limit(2)
          stoich_data%maxlist(i)=1
          stoich_data%list(i)=0
          stoich_data%list_increment(i)=1  
        Else
          stoich_data%list(i)=0
          stoich_data%list_increment(i)=1  
        End If
        itot=itot*stoich_data%maxlist(i)
      End Do

    End If

    ix=1
    Do While (ix /= 0)
      stoich_data%list(ix)=stoich_data%list(ix)+ stoich_data%list_increment(ix)
      indx= stoich_data%index_indep(ix)
      stoich_data%species(indx)%s= stoich_data%species(indx)%limit(1) + (stoich_data%list(ix)-1)*Delta_stoich
      If (stoich_data%list(ix) <= stoich_data%maxlist(ix) .And. &
         stoich_data%list(ix) >= 1 ) Then
        If (ix == stoich_data%num_indep) Then
          Call solve_stoichiometry(stoich_data, fsol_dep)
          If (fsol_dep) Then
            Call check_validity_solution(stoich_data, fsol)
            If (fsol) Then
              If (operation=='upper_limit') Then
                isol_domain=isol_domain+1
                Do i=1, stoich_data%num_indep
                  indx=stoich_data%index_indep(i)
                  If (stoich_data%species(indx)%limit(2)<stoich_data%species(indx)%s) Then
                    stoich_data%species(indx)%limit(2)=stoich_data%species(indx)%s
                  End If
                End Do
              ElseIf (operation=='find_solutions') Then
                Call obtain_electrons(stoich_data)
                isol_found=isol_found+1
                If (felect) Then
                  Call print_stoich_file(stoich_data, iunit, 'electrons')
                   felect=.False.
                End If
                Call print_stoich_file(stoich_data, iunit, 'coefficients')
              End If
            End If
          End If
        Else 
         ix=ix+1
        End If

      Else 
      ! end of loop - reset counter
        If (stoich_data%species(indx)%fmax) Then
          stoich_data%list(ix) = stoich_data%maxlist(ix)+1 
        Else
          stoich_data%list(ix) = 0 
        End If
      ! Move to previous index
        ix = ix - 1
      End If

    End Do


    If (operation=='find_solutions') Then
      If (isol_found==1) Then
        Do k=1, stoich_data%num_species%value 
          stoich_data%species(k)%s0= stoich_data%species(k)%s + stoich_data%species(k)%s0
        End Do
      Else If (isol_found > 1) Then
        Call info(' ', 1)
        stoich_data%continue_cv=.False.
        Write (message, '(1x,2a,a,i2,a)') 'Multiple solutions generated during ', Trim(stoich_data%cv_leg), &
                                    &' of cycle ', stoich_data%cv_cycle, '. The stoichiometric analysis&
                                    & stops here.'
         Call info(message,1)
        If (stoich_data%single_sol) Then
          Write (message, '(a,i2,3a)')  '***PROBLEMS: A total of ', stoich_data%nconst_leg, ' constraints have&
                                   & been set for ', Trim(stoich_data%cv_leg), &
                                   &' which should have led to A SINGLE stoichiometric solution.'
          Call info(message, 1)
          Write (message, '(1x,3a)') ' Please check the definition of constraints for ', Trim(stoich_data%cv_leg),&
                               & ' (see Table with types of variables above).'
          Call info(message, 1)
        End If
      Else If (isol_found==0) Then
        Write (message, '(2a,a,i2)') '***WARNING:  No stoichiometric solution has been found during the ',&
                                  & Trim(stoich_data%cv_leg), ' of cycle ', stoich_data%cv_cycle   
        Call info(message, 1)
        stoich_data%sol_exist=.False.
      End If
      stoich_data%numsol(jleg, icycle)=isol_found
    End If

    If (operation=='upper_limit' .And. isol_domain==0) Then
       Write (message, '(2a,a,i2)') '***ERROR: No stoichiometric solution has been found when computing the DOMAINS&
                                  & of variables for ', Trim(stoich_data%cv_leg), ' of cycle ', stoich_data%cv_cycle   
       Call info(message, 1)
       If (stoich_data%nconst_leg/=0) Then
         Write (message, '(1x,2a)') ' Please check the definition of constraints for ', Trim(stoich_data%cv_leg)
         Call info(message, 1)
       End If
       Call error_stop(' ')
    End If 

  End Subroutine sample_stoich_phase_space

  Subroutine check_validity_solution(stoich_data, fsol)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if the computed solution is within a realistic range and complies 
    ! with the contraints (if set)
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),    Intent(InOut) :: stoich_data
    Logical,              Intent(InOut) :: fsol
 
    Integer(Kind=wi)       :: i

    fsol=.True.

    If (stoich_data%nconst_leg/=0) Then
      ! Set if constrained variables are within a physical range  
      Call setcheck_constrained_solutions(stoich_data, fsol)

      ! Check is variable has reached minimum
      i=1
      Do While (i <= stoich_data%num_species%value .And. fsol)
        If (stoich_data%species(i)%fmin) Then
          If (stoich_data%species(i)%s <= stoich_data%species(i)%minvalue) Then
            stoich_data%species(i)%minvalue=stoich_data%species(i)%s
            fsol=.True.
          Else 
            fsol=.False. 
          End If
        End If  
        i=i+1
      End Do

      ! Check is variable has reached maximum
      i=1
      Do While (i <= stoich_data%num_species%value .And. fsol)
        If (stoich_data%species(i)%fmax) Then
          If (stoich_data%species(i)%s >= stoich_data%species(i)%maxvalue) Then
            stoich_data%species(i)%maxvalue=stoich_data%species(i)%s
            fsol=.True.
          Else 
            fsol=.False. 
          End If
        End If  
        i=i+1
      End Do

      ! Check if variable is within the range
      i=1
      Do While (i <= stoich_data%num_species%value .And. fsol)
        If (stoich_data%species(i)%frange) Then
          If ((stoich_data%species(i)%s+stoich_data%species(i)%s0)>= stoich_data%species(i)%range(1) .And. &
             (stoich_data%species(i)%s+stoich_data%species(i)%s0)<= stoich_data%species(i)%range(2)  ) Then
            fsol=.True.
          Else 
            fsol=.False. 
          End If
        End If  
        i=i+1
      End Do

    End If


   End Subroutine check_validity_solution


  Subroutine setcheck_constrained_solutions(stoich_data, fsol)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute solutions under constraints and check if they are consistent 
    ! with the domain 
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),    Intent(InOut) :: stoich_data
    Logical,              Intent(InOut) :: fsol            

    Integer(Kind=wi) :: i, j
    Real(Kind=wp) :: svec(stoich_data%num_species%value)

    ! Set solutions
    svec(:)=stoich_data%species(:)%s
    stoich_data%species(:)%s=0.0_wp
    Do i=1, stoich_data%num_species%value
      Do j=1, stoich_data%num_species%value
        stoich_data%species(i)%s=stoich_data%species(i)%s + &
        stoich_data%matrix_constraints(i,j)*svec(j)
      EndDo
    End Do

    ! Check solutions
    Do i=1, stoich_data%num_species%value
      If (stoich_data%species(i)%vartype == 'constrained') Then
        If (stoich_data%species(i)%s<stoich_data%species(i)%limit(1)) Then
          fsol=.False.
        End If
      End If
    End Do   

  End Subroutine setcheck_constrained_solutions 


  Subroutine print_stoich_file(stoich_data, iunit, fragment)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print stoichiometric solutions
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type), Intent(In   ) :: stoich_data
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Character(Len=*),  Intent(In   ) :: fragment

    Integer(Kind=wi) :: k
    
    If (stoich_data%block_constraints%fread) Then
      If (fragment=='header') Then
        If (stoich_data%single_sol) Then   
          Write (iunit,'(a,*(a12,2x))') '# Cycle',&
                                          (Trim(stoich_data%species(stoich_data%index_indep(k))%tag), k=1,stoich_data%num_indep),&
                                          (Trim(stoich_data%species(stoich_data%index_dep(k))%tag),   k=1,stoich_data%num_dep),&
                                          (Trim(stoich_data%species(stoich_data%index_const(k))%tag), k=1,stoich_data%net_const), &
                                         'Electrons'
        Else
          Write (iunit,'(a1,*(a12,2x))') '#',&
                                          (Trim(stoich_data%species(stoich_data%index_indep(k))%tag), k=1,stoich_data%num_indep),&
                                          (Trim(stoich_data%species(stoich_data%index_dep(k))%tag),   k=1,stoich_data%num_dep),&
                                          (Trim(stoich_data%species(stoich_data%index_const(k))%tag), k=1,stoich_data%net_const)
        End If
      Else If (fragment=='cycle_only') Then
        Write (iunit,'(a,i3)')   '# Cycle number: ', stoich_data%cv_cycle

      Else If (fragment=='electrons' .And. (.Not. stoich_data%single_sol)) Then
        Write (iunit,'(a,f12.6)')   '# Electrons: ', stoich_data%electrons 

      ElseIf (fragment=='coefficients') Then
        If (stoich_data%single_sol) Then
           Write (iunit,'(4x,i3,*(f12.6,2x))') stoich_data%cv_cycle, &
                                         (stoich_data%species(stoich_data%index_indep(k))%s+&
                                           stoich_data%species(stoich_data%index_indep(k))%s0,        k=1,stoich_data%num_indep),& 
                                          (stoich_data%species(stoich_data%index_dep(k))%s  + &
                                           stoich_data%species(stoich_data%index_dep(k))%s0,          k=1,stoich_data%num_dep),&
                                           (stoich_data%species(stoich_data%index_const(k))%s+  &
                                           stoich_data%species(stoich_data%index_const(k))%s0,        k=1,stoich_data%net_const),&
                                           stoich_data%electrons
        Else
           Write (iunit,'(1x,*(f12.6,2x))') (stoich_data%species(stoich_data%index_indep(k))%s+&
                                           stoich_data%species(stoich_data%index_indep(k))%s0,        k=1,stoich_data%num_indep),& 
                                          (stoich_data%species(stoich_data%index_dep(k))%s  +&
                                           stoich_data%species(stoich_data%index_dep(k))%s0,          k=1,stoich_data%num_dep),&
                                          (stoich_data%species(stoich_data%index_const(k))%s  +&
                                          stoich_data%species(stoich_data%index_const(k))%s0,         k=1,stoich_data%net_const)
        End If 
      End If  
    Else
      If (fragment=='header') Then
        If (stoich_data%single_sol) Then   
          Write (iunit,'(a,(*(a12,2x)))') '# Cycle',&
                                          (Trim(stoich_data%species(stoich_data%index_indep(k))%tag), k=1,stoich_data%num_indep),&
                                          (Trim(stoich_data%species(stoich_data%index_dep(k))%tag),   k=1,stoich_data%num_dep),  &
                                          'Electrons'
        Else
          Write (iunit,'(a1,(*(a12,2x)))') '#',&
                                          (Trim(stoich_data%species(stoich_data%index_indep(k))%tag), k=1,stoich_data%num_indep),&
                                          (Trim(stoich_data%species(stoich_data%index_dep(k))%tag),   k=1,stoich_data%num_dep)
        End If
      Else If (fragment=='cycle_only') Then
        Write (iunit,'(a,i3)')           '# Cycle number: ', stoich_data%cv_cycle

      Else If (fragment=='electrons' .And. (.Not. stoich_data%single_sol)) Then
        Write (iunit,'(a,f12.6)')   '# Electrons: ', stoich_data%electrons 


      ElseIf (fragment=='coefficients') Then
        If (stoich_data%single_sol) Then
           Write (iunit,'(4x,i3,*(f12.6,2x))') stoich_data%cv_cycle,  (stoich_data%species(stoich_data%index_indep(k))%s+&
                                           stoich_data%species(stoich_data%index_indep(k))%s0,        k=1,stoich_data%num_indep),& 
                                          (stoich_data%species(stoich_data%index_dep(k))%s  + &
                                           stoich_data%species(stoich_data%index_dep(k))%s0,          k=1,stoich_data%num_dep),  &
                                           stoich_data%electrons
        Else
           Write (iunit,'(1x,*(f12.6,2x))') (stoich_data%species(stoich_data%index_indep(k))%s+&
                                           stoich_data%species(stoich_data%index_indep(k))%s0,        k=1,stoich_data%num_indep),& 
                                          (stoich_data%species(stoich_data%index_dep(k))%s  +&
                                          stoich_data%species(stoich_data%index_dep(k))%s0,           k=1,stoich_data%num_dep)
        End IF   
      End If  
    End If
   
  End Subroutine print_stoich_file

  Subroutine solve_stoichiometry(stoich_data, fsol_dep)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Determine the stoichiometric coefficients from the solution of mass and
    ! charge balance equations. This subroutine is called only if the number 
    ! of involved variables is >= 2
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),            Intent(InOut) :: stoich_data
    Logical,                      Intent(  Out) :: fsol_dep

    Integer(Kind=wi)   :: i, j, k
    Real(Kind=wp)      :: sol(2), aux(2)
    Character(Len=256) :: message

    If (stoich_data%num_variables == 1) Then
      Write (message,'(a)') '***ERROR: stoichiometry analysis consistent with EQCM data is NOT possible using&
                           & a single (only one) variable chemical species for intercalation processes. Please review settings.'
      Call error_stop(message)
    Else
      aux=stoich_data%matrix_redox
      Do j=1,stoich_data%num_species%value
        If (stoich_data%species(j)%vartype == 'independent') Then
          aux(:)=aux(:)-stoich_data%matrix_balance(:,j)*stoich_data%species(j)%s
        End If
      End Do
      
      sol = 0.0_wp
      Do k=1, 2
       Do j=1,2
         sol(k)=sol(k)+stoich_data%inv_matrix_eq(k,j)*aux(j)
       End Do
      End Do
    End If 
    i=1
    fsol_dep=.True.
    Do While (i <= stoich_data%num_dep .And. fsol_dep)
      fsol_dep= (stoich_data%species(stoich_data%index_dep(i))%limit(1) <= sol(i))
      i=i+1
    End Do
  
    If (fsol_dep) Then
      Do i=1, stoich_data%num_dep
        stoich_data%species(stoich_data%index_dep(i))%s=sol(i)
      End Do  
    End If
  
  End Subroutine solve_stoichiometry

  Subroutine obtain_electrons(stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute number of accumulated electrons for reaction
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),    Intent(InOut) :: stoich_data

    Integer(Kind=wi) :: j

    stoich_data%electrons= 0.0_wp
    Do j=1,stoich_data%num_species%value
      If (stoich_data%species(j)%vartype /= 'fixed') Then
        stoich_data%electrons=stoich_data%electrons-stoich_data%species(j)%ox*stoich_data%species(j)%s
      End If
    End Do

  End Subroutine obtain_electrons 


  Subroutine classify_variables(stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Re-classifies stoichiometric variables when constraints are applied
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),    Intent(InOut) :: stoich_data

    Integer(Kind=wi) :: i, j, k, m
    Integer(Kind=wi) :: lc
    Integer(Kind=wi) :: ic1, ic2
    Integer(Kind=wi) :: kdep, kindep, kconst
    Character(Len=256) :: message
    Real(Kind=wp)      :: fix_value

    lc=0
    ic1=0 ; ic2=0
    stoich_data%species(:)%fmax=.False.
    stoich_data%species(:)%fmin=.False.
    stoich_data%species(:)%frange=.False.
    stoich_data%single_sol=.False.
    stoich_data%nconst_leg=0

    ! Set the type of variables back to the definition of &block_species
    stoich_data%species(:)%vartype=stoich_data%species(:)%vartype0 

    ! If contraints are set, define the new vartype for the chemical species
    If (stoich_data%block_constraints%fread) Then
      !Define de contraints matrix
      stoich_data%matrix_constraints=0.0_wp
      Do i = 1, stoich_data%num_species%value
        stoich_data%matrix_constraints(i,i)=1.0_wp
      End Do   
      !Set array to indicate which atoms have changed label
      ! Loop over contraints
      Do i=1, stoich_data%num_constraints%value
        ! Only consider contraints for oxidation (reduction) only if the portion of the CV is oxidation (reduction)
        If (Trim(stoich_data%constraints(i)%leg)==Trim(stoich_data%cv_leg)) Then
          lc=lc+1
          ! If one fixes the stoichiometry...
          If (Trim(stoich_data%constraints(i)%type)=='target_value') Then
            Do j=1, stoich_data%num_species%value
              If (Trim(stoich_data%species(j)%tag)==Trim(stoich_data%constraints(i)%species(1))) Then
                stoich_data%species(j)%vartype='constrained'
                fix_value=stoich_data%constraints(i)%value(1)-stoich_data%species(j)%s0
                stoich_data%matrix_redox(1)=stoich_data%matrix_redox(1)-stoich_data%matrix_balance0(1,j)*fix_value
                stoich_data%matrix_redox(2)=stoich_data%matrix_redox(2)-stoich_data%matrix_balance0(2,j)*fix_value
                stoich_data%species(j)%s = fix_value
              End If
            End Do

          ! Settings to consider solutions within a given range 
          ElseIf (Trim(stoich_data%constraints(i)%type)=='target_range') Then
            Do j=1, stoich_data%num_species%value
              If (Trim(stoich_data%species(j)%tag)==Trim(stoich_data%constraints(i)%species(1))) Then
                stoich_data%species(j)%vartype='independent'
                stoich_data%species(j)%frange=.True.  
                stoich_data%species(j)%range(1)=stoich_data%constraints(i)%value(1)
                stoich_data%species(j)%range(2)=stoich_data%constraints(i)%value(2)
                If (Abs(stoich_data%species(j)%range(1)-stoich_data%species(j)%range(2))< epsilon(1.0_wp) .Or. &
                  (stoich_data%species(j)%range(1) > stoich_data%species(j)%range(2)) ) Then  
                  Write (message, '(1x,3a)') '***ERROR in &blocK_constraints_species: Problems in the definition&
                                           & "target_range" for species ',&
                                           & Trim(stoich_data%species(j)%tag), '. The first value (lower bound)&
                                           & must be lower the second bound (upper bound).'
                  Call info(message, 1)
                  Call error_stop(' ')
                End If          
              End If
            End Do

          ! Settings to target min and maximum values for stoichiometric coefficients 
          ElseIf (Trim(stoich_data%constraints(i)%type)=='target_min' .Or. &
                 Trim(stoich_data%constraints(i)%type)=='target_max') Then
            Do j=1, stoich_data%num_species%value
              If (Trim(stoich_data%species(j)%tag)==Trim(stoich_data%constraints(i)%species(1))) Then
                stoich_data%species(j)%vartype='independent'
                If (Trim(stoich_data%constraints(i)%type)=='target_min') Then
                  stoich_data%species(j)%fmin=.True.  
                  stoich_data%species(j)%minvalue=Huge(1.0_wp)
                End If
                If (Trim(stoich_data%constraints(i)%type)=='target_max') Then
                  stoich_data%species(j)%fmax=.True.  
                  stoich_data%species(j)%maxvalue=-Huge(1.0_wp)
                End If
              End If
            End Do


          ! Settings for linear constraints between variables 
          ElseIf (Trim(stoich_data%constraints(i)%type)=='keep_ratio' .Or. &
                 Trim(stoich_data%constraints(i)%type)=='ratio_fixed') Then
            Do j=1, stoich_data%num_species%value
              If (Trim(stoich_data%species(j)%tag)==Trim(stoich_data%constraints(i)%species(2))) Then
                ic2=j
                stoich_data%species(j)%vartype='dependent'
              End If
              If (Trim(stoich_data%species(j)%tag)==Trim(stoich_data%constraints(i)%species(1))) Then
                ic1=j
                stoich_data%species(j)%vartype='constrained'
              End If
            End Do
            stoich_data%matrix_constraints(ic1,ic1)=0.0_wp
            If (Trim(stoich_data%constraints(i)%type)=='ratio_fixed') Then         
              stoich_data%matrix_constraints(ic1,ic2)=stoich_data%constraints(i)%value(1)
            ElseIf (Trim(stoich_data%constraints(i)%type)=='keep_ratio') Then
              If (Trim(stoich_data%constraints(i)%leg) /= Trim(stoich_data%first_leg)) Then
                stoich_data%matrix_constraints(ic1,ic2)=stoich_data%species(ic1)%s/stoich_data%species(ic2)%s
              Else
                Write (message,'(5a)')'***ERROR: definition of "keep_ratio" type of constraint during ',&
                                  & Trim(stoich_data%constraints(i)%leg), ' is wrong, as ', Trim(stoich_data%constraints(i)%leg),& 
                                  &' takes place first in the CV cycling. Use of "keep_ratio" MUST be'
                Call info(message,1)
                Write (message,'(a)') 'applied to the portion of the CV cycle that follows the first fragment.&
                                  & For example, if oxidation occurs first, the definition of "keep_ratio" '
                Call info(message,1)
                Write (message,'(a)') 'is only valid during reduction, and vice versa. Please review the definition&
                                   & of constraints in &block_constraint_species'
                Call info(message,1)
                Call error_stop(' ')
              End If
            End If
          End If

        End If
      End Do

      !Count number of dependent variables
      kdep=0 
      Do j= 1, stoich_data%num_species%value
        If (Trim(stoich_data%species(j)%vartype)=='dependent') Then
          kdep=kdep+1
        End If
      End Do 


      ! Check over the species that have changed type
      j=1; k=0
      Do While (j <= stoich_data%num_species%value) 
          If (Trim(stoich_data%species(j)%vartype)/='fixed' .And. &
           Trim(stoich_data%species(j)%vartype)/='constrained' .And. kdep/=2) Then
           If (kdep<2) Then
             If (Trim(stoich_data%species(j)%vartype)=='independent') Then
               stoich_data%species(j)%vartype='dependent'
               k=k+1
               kdep=kdep+1
             End If
           ElseIf (kdep>2) Then
             If (Trim(stoich_data%species(j)%vartype)=='dependent') Then
               stoich_data%species(j)%vartype='independent'
               k=k+1
               kdep=kdep-1
             End If
           End If
         End If
       j=j+1
      End Do

    ! Count how many dependent, independent and total stoichimetric variables are left
    ! Constrained variables are also numbered 
      kdep=0 ; kindep=0; kconst=0
      Do j= 1, stoich_data%num_species%value
        If (Trim(stoich_data%species(j)%vartype)=='independent') Then
          kindep=kindep+1
        End If
        If (Trim(stoich_data%species(j)%vartype)=='dependent') Then
          kdep=kdep+1
        End If
        If (Trim(stoich_data%species(j)%vartype)=='constrained') Then
          kconst=kconst+1
        End If
      End Do
      
      ! Define new quantities
      stoich_data%num_dep=kdep
      stoich_data%num_indep=kindep
      stoich_data%net_const=kconst

      ! Number of constraints for a oxidation/reduction
      stoich_data%nconst_leg=lc

      ! Print Contraints to file, both for oxidation and reduction, depending on what happens first
      If (stoich_data%cv_cycle==1) Then
        If (stoich_data%nconst_leg == 0) Then
          Write (message,'(1x,3a)') '-', Trim(stoich_data%cv_leg), ': No constraints set'
        Else 
          Write (message,'(1x,3a,i2)') '-', Trim(stoich_data%cv_leg), ': Number of constraints is equal to ', lc
        End If
        Call info(message,1)
        If (stoich_data%nconst_leg == 0) Then
          Call info('  The type of variable for each stoichiometric species remains as for &block_species above', 1)  
        Else
          Write(message,'(a)') '  The type of variable for each species is set as follows (it might remain same as&
                              & the definition of &block_species above):'
          Call info(message, 1) 
          Call info(' ---------------------------------', 1)
          Write (message, '(6x,2(a,5x))') 'Species', 'Type of variable'
          Call info(message,1)
          Call info(' ---------------------------------', 1)
          Do i=1, stoich_data%num_species%value
            Write (message, '(1x,a12,5x,a)')  Trim(stoich_data%species(i)%tag), Trim(stoich_data%species(i)%vartype)
            Call info(message,1)
          End Do
          Call info(' ---------------------------------', 1)
        End If
      End If

    End If
    
    ! Define new indexes
    k=0; m=0; i=0
    Do j= 1, stoich_data%num_species%value
      If (stoich_data%species(j)%vartype=='dependent') Then
        k=k+1 
        stoich_data%index_dep(k)=j
      End If
      If (stoich_data%species(j)%vartype=='independent') Then
        m=m+1
        stoich_data%index_indep(m)=j
      End If
      If (stoich_data%species(j)%vartype=='constrained') Then
        i=i+1
        stoich_data%index_const(i)=j
      End If
    End Do

    If ((stoich_data%num_variables-lc)<=2) Then
      stoich_data%single_sol=.True.
      If (stoich_data%block_constraints%fread) Then
        i=1
        Do While ((i<= stoich_data%num_constraints%value) .And. stoich_data%single_sol)
          If (Trim(stoich_data%constraints(i)%leg)==Trim(stoich_data%cv_leg)) Then
            If (Trim(stoich_data%constraints(i)%type)=='target_range') Then
              stoich_data%single_sol=.False.
            End If
          End If 
          i=i+1
        End Do
      End If
    End If

    ! Set initial values of stoichiometric solutions to zero
    Do j=1,stoich_data%num_species%value
      If (stoich_data%species(j)%vartype /= 'fixed' .And. &
         stoich_data%species(j)%vartype /= 'constrained'  ) Then
         stoich_data%species(j)%s=0.0_wp
      End If
    End Do   


  End Subroutine classify_variables  

  Subroutine set_matrix_balance(stoich_data, stage)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set matrix for charge and mass equations  
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),    Intent(InOut) :: stoich_data
    Character(Len=*),     Intent(In   ) :: stage

    Real(Kind=wp)       :: mol, det
    Integer(Kind=wi)    :: i, j, k
    Character(Len=256)  :: message
    
    mol=stoich_data%mole(1)

    If (stage=='initial') Then
      !Matrix balance 
      Do j= 1, stoich_data%num_species%value
        If (stoich_data%species(j)%vartype0=='fixed') Then
          stoich_data%matrix_balance0(:,j) = 0.0_wp
        Else 
          stoich_data%matrix_balance0(1,j) = stoich_data%species(j)%mass*g_to_ng
          stoich_data%matrix_balance0(2,j) = -stoich_data%species(j)%ox
        End If 
      End Do 

    Else If (stage=='cv') Then
      If (stoich_data%block_constraints%fread) Then 
        Do i=1, 2
          Do j=1,stoich_data%num_species%value
            stoich_data%matrix_balance(i,j)=0.0_wp
            Do k=1, stoich_data%num_species%value
              stoich_data%matrix_balance(i,j)=stoich_data%matrix_balance(i,j) + &
              & stoich_data%matrix_balance0(i,k)*stoich_data%matrix_constraints(k,j)
            End Do
          End Do
        End Do
      Else
        stoich_data%matrix_balance=stoich_data%matrix_balance0
      End If

      k=0
      !Matrix of dependent variables
      Do j= 1, stoich_data%num_species%value
        If (stoich_data%species(j)%vartype=='dependent') Then
          k=k+1
          stoich_data%matrix_eq(:,k) = stoich_data%matrix_balance(:,j)
        End If
      End Do
      ! Inverse matrix
      det=stoich_data%matrix_eq(1,1)*stoich_data%matrix_eq(2,2)-stoich_data%matrix_eq(2,1)*stoich_data%matrix_eq(1,2)
      If (Abs(det)<epsilon(det)) Then
        Write (message,'(a)') '***PROBLEMS: Determinant of dependent matrix is zero! Check details of block_species!'
        Call info(message, 1)
        Write (message,'(a)') '   At least one of the dependent variables must have a non-zero oxidation number!'
        Call info(message, 1)
        Call error_stop(' ')
      End If
      stoich_data%inv_matrix_eq(1,1)= stoich_data%matrix_eq(2,2)
      stoich_data%inv_matrix_eq(2,2)= stoich_data%matrix_eq(1,1)
      stoich_data%inv_matrix_eq(2,1)=-stoich_data%matrix_eq(2,1)
      stoich_data%inv_matrix_eq(1,2)=-stoich_data%matrix_eq(1,2)
      stoich_data%inv_matrix_eq(:,:)=stoich_data%inv_matrix_eq(:,:)/det
    End If

  End Subroutine set_matrix_balance  

  Subroutine set_mole(eqcm_data, electrode_data, stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate mole of material from the formula unit and electrode mass,
    ! only if electrode_mole is not set in the SET_EQCM file
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(electrode_type), Intent(In   ) :: electrode_data
    Type(eqcm_type),      Intent(InOut) :: eqcm_data 
    Type(stoich_type),    Intent(InOut) :: stoich_data

    Integer(Kind=wi) :: i
    Character(Len=256)  :: messages(2)

    If ((.Not. electrode_data%mass%fread .And. .Not. electrode_data%mole%fread)) Then
      Write (messages(1),'(a,i2,a)') '***PROBLEMS: Stoichometric analysis for this case has ', stoich_data%num_variables, &
                                   &' variable species.'
      Write (messages(2),'(a)') 'Solution requires the specification of either "electrode_mass" or "electrode_moles" directives &
                               &in file SET_EQCM.'
      Call info(messages, 2)
      Call error_stop(' ')
    End If

    If (electrode_data%mole%fread) Then
      stoich_data%mole(1)=electrode_data%mole%value
      stoich_data%mole(1)=stoich_data%mole(1)*eqcm_data%efficiency%value
    Else
    Do  i=1, stoich_data%num_species%value
      stoich_data%molar_mass=stoich_data%molar_mass+stoich_data%species(i)%mass*stoich_data%species(i)%s0
    End Do
      Write (messages(1),'(1x,a,f12.3,a)') 'Molar mass of pristine electrode:    ', stoich_data%molar_mass, ' [g/mol]'
      Call info(messages, 1)
      stoich_data%mole(1)=electrode_data%mass%value/stoich_data%molar_mass/g_to_ng
      Write (messages(1),'(1x,a,26x,E14.7)')   'Calculated moles:', stoich_data%mole(1)
      Call info(messages, 1)

      If (eqcm_data%efficiency%fread) Then
        Write (messages(1),'(1x,a,10x,f6.3)') 'Efficiency set for the reaction: ', eqcm_data%efficiency%value 
      Else
        Write (messages(1),'(1x,a,f6.3)') 'Efficiency set for the reaction (default): ', eqcm_data%efficiency%value 
      End If
      Call info(messages,1)

      stoich_data%mole(1)=stoich_data%mole(1)*eqcm_data%efficiency%value
      Write (messages(1),'(1x,a,4x,E14.7)')   'Moles that participate in the reaction:', stoich_data%mole(1)
      Call info(messages, 1)
    End If
    
  End Subroutine set_mole

  Subroutine print_stoichiometric_settings(eqcm_data, stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print relevant input information for each species as specified in 
    ! &block_species 
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),    Intent(InOut) :: eqcm_data
    Type(stoich_type),  Intent(InOut) :: stoich_data

    Character(Len=256)  :: message

     Call info(' ',1)
    Write (message,'(1x,a,i2)') 'Species that characterise the reaction:    ', &
                               &stoich_data%num_species%value
    Call info(message,1)
    Write (message,'(1x,a,i2)') 'Species that participate in the reaction:  ', &
                               &stoich_data%num_variables
    Call info(message,1)

    If (Trim(eqcm_data%process%type)=='intercalation' .And. stoich_data%num_indep > 0) Then
      If (stoich_data%discretization%fread) Then
        Write (message,'(1x,a,f6.3)') 'Discretization of the stoichiometric space: ', stoich_data%discretization%value 
      Else
        Write (message,'(1x,a,f6.3)') 'Discretization of the stoichiometric space (default):', stoich_data%discretization%value
      End if
      Call info(message,1)
    End If

    Call print_pristine_stoichiometry(stoich_data)

  End Subroutine print_stoichiometric_settings 

  Subroutine print_pristine_stoichiometry(stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print info of &block_species 
    ! 
    ! author    - i.scivetti Feb 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type),  Intent(InOut) :: stoich_data

    Integer(Kind=wi) :: i
    Character(Len=256)  :: message
    Character(Len= 64)  :: fmt1, fmt2

    fmt1='(6x,5(a,6x))'
    fmt2='(1x,a12,5x,f12.3,18x,f5.2,22x,f6.3,6x,a)' 

    Call info(' ', 1)
    Call info(' Summary of the information defined in &block_species:', 1)
    Call info(' ------------------------------------------------------------------------------------------------------', 1)     
    Write (message, fmt1) 'Species', 'Mass [g/mol]', 'Oxidation number', 'Pristine stoichiometry', 'Type of variable'   
    Call info(message,1)
    Call info(' ------------------------------------------------------------------------------------------------------', 1)     
    Do i=1, stoich_data%num_species%value
      Write (message, fmt2)  Trim(stoich_data%species(i)%tag), stoich_data%species(i)%mass,  stoich_data%species(i)%ox, &
                           stoich_data%species(i)%s0_pristine, Trim(stoich_data%species(i)%vartype0)
      Call info(message,1)
    End Do
    Call info(' ------------------------------------------------------------------------------------------------------', 1)     

  End Subroutine print_pristine_stoichiometry

  Subroutine print_deposition(stoich_data, files, icycle, record_mass_freq, record_mass)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print ELECTRODEPOSITION file
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(stoich_type), Intent(InOut) :: stoich_data
    Type(file_type),   Intent(InOut) :: files(:)
    Integer(Kind=wi),  Intent(In   ) :: icycle 
    Logical,           Intent(In   ) :: record_mass_freq
    Logical,           Intent(In   ) :: record_mass    

    Integer(Kind=wi) :: iunit, i
    Character(Len=256)  :: message

    Open(Newunit=files(FILE_ELECTRO)%unit_no, File=files(FILE_ELECTRO)%filename, Status='Replace')
    iunit=files(FILE_ELECTRO)%unit_no
    Write (message,'(1x,2a)') 'Print results to file ', Trim(FOLDER_ANALYSIS)//'/'//Trim(files(FILE_ELECTRO)%filename)
    Call info(' ', 1)
    If (record_mass) Then
        Write (iunit,'(a,(*(a18,2x)))') '#Cycle', '[Dm(Q)/Dm_eff]_ox', '[Dm(Q)/Dm_eff]_red', 'mole_ox [nmol]', 'mole_red [nmol]'
    Else
      If (record_mass_freq) Then       
        Write (iunit,'(a,(*(a18,2x)))') '#Cycle', 'A_eff/A_geom(ox)', 'A_eff/A_geom(red)', 'mole_ox [nmol]', 'mole_red [nmol]'
      End If  
    End If

    Call info(message,1)
    Do i=1, icycle
      Write (iunit,'(1x,i2,(*(f18.3,2x)))')  i, stoich_data%elec_depo_ratio_ox(i), stoich_data%elec_depo_ratio_red(i),&
                                            &  stoich_data%mole_ox(i)/1.0E-9, stoich_data%mole_red(i)/1.0E-9 
    End Do

    Close(iunit)

    ! Move file
    Call execute_command_line('mv '//Trim(files(FILE_ELECTRO)%filename)//' '//Trim(FOLDER_ANALYSIS))

  End Subroutine print_deposition

 
  Subroutine electro_deposition(ic, jc, DM, DQ, stoich_data, record_mass_freq, record_mass)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute electrodeposition quantities
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: ic        
    Integer(Kind=wi),  Intent(In   ) :: jc       
    Real(Kind=wp),     Intent(In   ) :: DM       
    Real(Kind=wp),     Intent(In   ) :: DQ       
    Type(stoich_type), Intent(InOut) :: stoich_data
    Logical,           Intent(In   ) :: record_mass_freq
    Logical,           Intent(In   ) :: record_mass

    Real(Kind=wp)      :: ratio, mol, mtot, oxtot, dMQ
    Character(Len=256) :: message
    Integer(Kind=wi)   ::  i 

    mtot=0.0_wp
    oxtot=0.0_wp

    Do i=1,stoich_data%num_indep
      mtot=mtot+stoich_data%species(stoich_data%index_indep(i))%mass
      oxtot=oxtot+stoich_data%species(stoich_data%index_indep(i))%ox
    End Do

    If (Abs(oxtot)<epsilon(1.0_wp)) Then
      Write (message,'(a)') '***ERROR: the sum of the oxidation states for all the independent variables is zero!!!&
                          & Electrodepositon analysis is not possible. Please check'
      Call error_stop(message)
    End If

    If (record_mass) Then
      dMQ=-DQ*(mtot*g_to_ng)/(oxtot*Farad)
      ratio=DMQ/DM
      If (ratio < 0.0_wp) Then
        Write (message,'(3a,i2,a)') '***ERROR: Computed DMeff/DMQ ratio for ', Trim(stoich_data%cv_leg), &
                                  &' of cycle ', stoich_data%cv_cycle, ' is negative. Please review the setting&
                                  & for the selected sign of the cathodic/anodic current'
        Call error_stop(message)
      End If
      mol=-DQ/(oxtot*Farad)    
    Else
      If (record_mass_freq) Then
        ratio=-mtot*g_to_ng/Farad 
        ratio= ratio/oxtot/(DM/DQ) 
        If (ratio < 0.0_wp) Then
          Write (message,'(3a,i2,a)') '***ERROR: Computed A_eff/A_geom ratio for ', Trim(stoich_data%cv_leg), &
                                    &' of cycle ', stoich_data%cv_cycle, ' is negative. Please review the setting&
                                    & for the selected sign of the cathodic/anodic current'
          Call error_stop(message)
        End If
        mol=DM*ratio/(mtot*g_to_ng)
      End If
    End If        

    If (Trim(stoich_data%cv_leg)=='oxidation') Then
      stoich_data%elec_depo_ratio_ox(stoich_data%cv_cycle)=ratio
      stoich_data%mole_ox(stoich_data%cv_cycle)=mol
      stoich_data%solution_moles(jc,ic)=mol
    ElseIf (Trim(stoich_data%cv_leg)=='reduction') Then
      stoich_data%elec_depo_ratio_red(stoich_data%cv_cycle)=ratio
      stoich_data%mole_red(stoich_data%cv_cycle)=mol
      stoich_data%solution_moles(jc,ic)=mol
    End If

    stoich_data%numsol(jc,ic)=1

  End Subroutine electro_deposition


  Subroutine check_components_species(files, stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives for the atoms defined for 
    ! each component in &block_species_components
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(stoich_type),   Intent(InOut) :: stoich_data

    Integer(Kind=wi)  ::  i, j, k, m
    Character(Len=256)  :: messages(4)
    Character(Len=64 )  :: error_block_species_components
    Logical :: loop

    error_block_species_components = '***ERROR in &block_species_components of file '//Trim(files(FILE_SET_EQCM)%filename)

    Write (messages(1),'(a)') error_block_species_components

    ! Check correctness in the definition for the atomic components inside the same species
    Do i=1, stoich_data%num_species%value
      Do j= 1, (stoich_data%species(i)%num_components-1)
       Do k=j+1, stoich_data%species(i)%num_components
          If (Trim(stoich_data%species(i)%component%tag(j)) == &
            Trim(stoich_data%species(i)%component%tag(k)) ) Then
            Write (messages(2),'(2(a,i2),3a)') 'Tags for components ', j, ' and ', k, ' of species "', &
                                             & Trim(stoich_data%species(i)%tag), '" are the identical.'
            Write (messages(3),'(a)')          'Please use different local tags for components of the same species'
            Call info(messages, 3)
            Call error_stop(' ')
          End If 
          If (Trim(stoich_data%species(i)%component%element(j)) == &
            Trim(stoich_data%species(i)%component%element(k)) ) Then
            Write (messages(2),'(2(a,i2),3a)') 'Chemical element for components ', j, ' and ', k, ' of species "', &
                                             & Trim(stoich_data%species(i)%tag), '" are the identical.'
            Write (messages(3),'(a)')          'Please group atoms of the same type in only one component'
            Call info(messages, 3)
            Call error_stop(' ')
          End If 
        End Do
      End Do

      If (Trim(stoich_data%species(i)%topology) /= 'crystal' .And. & 
        Trim(stoich_data%species(i)%topology) /= 'solute' .And. &  
        Trim(stoich_data%species(i)%topology) /= 'molecule' .And. &  
        Trim(stoich_data%species(i)%topology) /= 'atom' ) Then
        Write (messages(2),'(4a)') 'Incorrect definition for the topology of species "', Trim(stoich_data%species(i)%tag), &
                                 & '": ', Trim(stoich_data%species(i)%topology)
                      
        Write (messages(3),'(a)')  'Valid options are: crystal, solute, molecule or atom'
        Call info(messages, 3)
        Call error_stop(' ')
      End If

      If (Trim(stoich_data%species(i)%topology) == 'crystal' .Or. Trim(stoich_data%species(i)%topology) == 'solute') Then
        If (Trim(stoich_data%species(i)%vartype0)/='fixed') Then
          Write (messages(2),'(5a)') 'Inconsistent definition for species "', Trim(stoich_data%species(i)%tag), &
                                 & '".  If the topology is chosen to be "', Trim(stoich_data%species(i)%topology),&
                                 & '", it must correspond to a "fixed" variable during the reaction'
                      
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End If

 
      !Check if global tags are consistent with the elements of the periodic table  
      Do j= 1, stoich_data%species(i)%num_components
        loop=.True.
        k=1
        Do While (k <= NPTE .And. loop)
          If (Trim(chemsymbol(k))==Trim(stoich_data%species(i)%component%element(j))) Then
           stoich_data%species(i)%component%atomic_number(j)=k
           loop=.False.
          End If
          k=k+1
        End Do
        If (loop) Then 
          Write (messages(2),'(3a,i2,3a)') 'Chemical element "' , Trim(stoich_data%species(i)%component%element(j)), &
                                         & '" defined for component ',  j, ' of species "', Trim(stoich_data%species(i)%tag), &
                                         & '" does not correspond to any element of the Periodic Table. Please use a valid chemical&
                                         & element'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End Do
    End Do

    ! Check correctness in the definition for the atomic components inside the same species
    Do i=1, stoich_data%num_species%value-1
      Do j = i+1, stoich_data%num_species%value 
        Do k=1, stoich_data%species(i)%num_components 
          Do m=1, stoich_data%species(j)%num_components
            If (Trim(stoich_data%species(i)%component%tag(k)) == &
              Trim(stoich_data%species(j)%component%tag(m)) ) Then
              Write (messages(2),'(3a,i2,a,i2,5a)')  'Tag ', Trim(stoich_data%species(i)%component%tag(k)), &
                                          & ' is the same for components ', k, ' and ' , m , ' of species ',            &
                                          & Trim(stoich_data%species(i)%tag), ' and ', Trim(stoich_data%species(j)%tag),&
                                          & ', respectively.' 
              Write (messages(3),'(a)')   'Please ammend! Tags for the components of species must be unique and not used&
                                          & in multiple species.'
              Call info(messages, 3)
              Call error_stop(' ')      
            End If
          End Do
        End Do        
      End Do
    End Do
 
  End Subroutine check_components_species


  Subroutine check_constraints_settings(files, stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives for the contraints on the
    ! solution of the stoichiometric problem
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(stoich_type),   Intent(InOut) :: stoich_data

    Integer(Kind=wi)  ::  i, j
    Integer(Kind=wi)  ::  ic_ox, ic_red
    Logical  ::  error, cond1, cond2
    Character(Len=256)  :: messages(4)
    Character(Len=64 )  :: error_set_eqcm
    Character(Len=64 )  :: error_block_constraints

    error_set_eqcm = '***ERROR in file '//Trim(files(FILE_SET_EQCM)%filename)//' -'
    error_block_constraints = '***ERROR in &block_constraints_species of file`'//Trim(files(FILE_SET_EQCM)%filename)

    error=.False.

    ! Check correctness of constraints for stoichiometric analysis
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Write (messages(1),'(a)') error_block_constraints

      Do i=1, stoich_data%num_constraints%value
        If (stoich_data%constraints(i)%fail) Then
          Write (messages(2),'(a,i2,a)') 'Problems with the specificaton for constraint ', i, &
                                      & ' (see manual for the correct format/syntax)'
          error=.True.
          Call info(messages,2)
        Else
          If (Trim(stoich_data%constraints(i)%leg) /= 'oxidation'  .And. &
            Trim(stoich_data%constraints(i)%leg) /= 'reduction') Then
            Write (messages(2),'(a,i2,3a)') 'The portion of the CV selected for constraint ', i,&
                                         & ' has been wrongly labelled as "', Trim(stoich_data%constraints(i)%leg),&
                                         & '".  It MUST be either "oxidation" or "reduction".'
            error=.True.
            Call info(messages,2)

          End If
        End If

       If (stoich_data%constraints(i)%value(1) < 0.0_wp .And. Trim(stoich_data%constraints(i)%type)=='target_value') Then
            Write (messages(2),'(3a)') 'Assigned value for any "', Trim(stoich_data%constraints(i)%type), '" type of constraint&
                                          & MUST ALWAYS BE >= 0.0: it is unphysical to derive cycled structures with negative&
                                          & stoichometric coefficients!'
            error=.True.
            Call info(messages,2)
       End If

       Call obtain_species_from_contraints(files,stoich_data, i)
     End Do


     Do i=1, stoich_data%num_constraints%value
       Do j=i+1, stoich_data%num_constraints%value
         !Check to prevent target_value constraint for a given species (and CV leg) to be duplicated
         If ((Trim(stoich_data%constraints(i)%type)=='target_value' .Or. Trim(stoich_data%constraints(i)%type)=='target_min' .Or. &
             Trim(stoich_data%constraints(i)%type)=='target_max' .Or. Trim(stoich_data%constraints(i)%type)=='target_range').And.&
            (Trim(stoich_data%constraints(j)%type)=='target_value' .Or. Trim(stoich_data%constraints(j)%type)=='target_min' .Or. &
             Trim(stoich_data%constraints(j)%type)=='target_max' .Or. Trim(stoich_data%constraints(j)%type)=='target_range')) Then
           If (Trim(stoich_data%constraints(i)%leg)==Trim(stoich_data%constraints(j)%leg)) Then
             If (Trim(stoich_data%constraints(i)%species(1))==Trim(stoich_data%constraints(j)%species(1))) Then
               Write (messages(2),'(5(a,4x))') 'The following constraints are redundant:'
               Write (messages(3),'(5(a,4x))') 'Constraint:', Trim(stoich_data%constraints(i)%type), &
                                               Trim(stoich_data%constraints(i)%leg),   &
                                               Trim(stoich_data%constraints(i)%tag_species)
               Write (messages(4),'(4(a,4x))') 'Constraint:', Trim(stoich_data%constraints(j)%type), &
                                               Trim(stoich_data%constraints(j)%leg),   &
                                               Trim(stoich_data%constraints(j)%tag_species)
               Call info(messages,4)
               error=.True.
               Call info('Please remove duplication', 1)
             End If
           End If
         End If

         !Check to prevent duplication of ratio_fixed constraint for a given pair of species (and CV leg)
         If ((Trim(stoich_data%constraints(i)%type)=='ratio_fixed' .Or. Trim(stoich_data%constraints(i)%type)=='keep_ratio') .And. &
            (Trim(stoich_data%constraints(j)%type)=='ratio_fixed' .Or. Trim(stoich_data%constraints(j)%type)=='keep_ratio') ) Then
           If (Trim(stoich_data%constraints(i)%leg)==Trim(stoich_data%constraints(j)%leg)) Then
             cond1=(Trim(stoich_data%constraints(i)%species(1))==Trim(stoich_data%constraints(j)%species(1))) .Or. &
                   (Trim(stoich_data%constraints(i)%species(1))==Trim(stoich_data%constraints(j)%species(2)))
             cond2=(Trim(stoich_data%constraints(i)%species(2))==Trim(stoich_data%constraints(j)%species(1))) .Or. &
                   (Trim(stoich_data%constraints(i)%species(2))==Trim(stoich_data%constraints(j)%species(2)))
             If (cond1 .And. cond2) Then
               Write (messages(2),'(4(a,4x))') 'Constraint:', Trim(stoich_data%constraints(i)%type), &
                                               Trim(stoich_data%constraints(i)%leg),   &
                                               Trim(stoich_data%constraints(i)%tag0_species)
               Call info(messages,2)
               error=.True.
               Call info('appears to be redundantly defined. Please review the settings for the defined contraints', 1)
             End If
           End If
         End If
       End Do
     End Do

     ic_ox=0
     ic_red=0

     Do i=1, stoich_data%num_constraints%value
       If (Trim(stoich_data%constraints(i)%leg) == 'oxidation') Then
         ic_ox=ic_ox+1
       Else If (Trim(stoich_data%constraints(i)%leg) == 'reduction') Then
         ic_red=ic_red+1
       End If
     End Do

     If ((stoich_data%num_variables-ic_ox)< 2) Then
        Write (messages(2),'((a,i2))') 'Number of defined contraints for "oxidation" is ', ic_ox
        Call info(messages,2)
        Write (messages(2),'((a,i2))') 'The number of (in)dependent stoichiometric variables defined in &block_species is ',&
                                      & stoich_data%num_variables
        Call info(messages(2),1)
        Write (messages(2),'(a)')      'The resulting system of equations is overdetermined. Please reduce the number of&
                                     & contraints for "oxidation".'
        Call info(messages(2),1)
        Write (messages(2),'((a,i2))') 'Maximum number of allowed constraints for "oxidation" is ', stoich_data%num_variables-2
        Call info(messages(2),1)
        error=.True.
     End If

     If ((stoich_data%num_variables-ic_red)< 2) Then
        Write (messages(2),'((a,i2))') 'Number of defined contraints for "reduction" is ', ic_red
        Call info(messages,2)
        Write (messages(2),'((a,i2))') 'The number of (in)dependent stoichiometric variables defined in &block_species is ',&
                                      & stoich_data%num_variables
        Call info(messages(2),1)
        Write (messages(2),'(a)')      'The resulting system of equations is overdetermined. Please reduce the number of&
                                     & contraints for "reduction".'
        Call info(messages(2),1)
        Write (messages(2),'((a,i2))') 'Maximum number of allowed constraints for "reduction" is ', stoich_data%num_variables-2
        Call info(messages(2),1)
        error=.True.
     End If

     If (error) Then
       Call error_stop(' ')
     End If

  End Subroutine check_constraints_settings

  Subroutine obtain_species_from_contraints(files,stoich_data, ic)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to:
    !  - indentify the involved constrained species
    !  - verify constraints are defined correctly
    !  - corroborate consistency with block_species
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(stoich_type),  Intent(InOut) :: stoich_data
    Integer(Kind=wi),   Intent(In   ) :: ic

    Integer(Kind=wi) :: i, j, k
    Integer(Kind=wi) :: nc, io, len
    Character(Len=256)  :: messages(2), string
    Character(Len=64 )  :: error_block_constraints_species

    error_block_constraints_species = '***ERROR in &block_constraints_species of file '&
                                      &//Trim(files(FILE_SET_EQCM)%filename)
    Write (messages(1),'(a)') error_block_constraints_species

    nc=0
    io=0

    If (Trim(stoich_data%constraints(ic)%type)=='target_value' .Or. &
       Trim(stoich_data%constraints(ic)%type)=='target_min'   .Or. &
       Trim(stoich_data%constraints(ic)%type)=='target_max'   .Or. &
       Trim(stoich_data%constraints(ic)%type)=='target_range'  ) Then
      stoich_data%constraints(ic)%species(1)=stoich_data%constraints(ic)%tag_species
      Do i = 1, stoich_data%num_species%value
        If (Trim(stoich_data%constraints(ic)%species(1))==Trim(stoich_data%species(i)%tag)) Then
          If (Trim(stoich_data%species(i)%vartype0)=='fixed') Then
            Write (messages(2),'(3a,i2,a)') 'Species ', Trim(stoich_data%constraints(ic)%species(1)),&
                             & ' (set for constraint ', ic,&
                             & ' of type "target_value") is defined as "fixed" &block_species.&
                             & Constraints MUST ONLY apply to (in)dependent variables'
            Call info(messages,2)
            Call error_stop(' ')
          Else
            nc=nc+1
          End If
        End If
      End Do
      If (nc == 0) Then
        Write (messages(2),'(3a,i2,a)') 'Species ', Trim(stoich_data%constraints(ic)%species(1)),&
                                      & ' (set for constraint ', ic,&
                                      & ' of type "target_value") has not been defined in &block_species'
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End If

    If (Trim(stoich_data%constraints(ic)%type)=='ratio_fixed' .Or. &
      Trim(stoich_data%constraints(ic)%type)=='keep_ratio') Then
      Call remove_symbols(stoich_data%constraints(ic)%tag_species, '/|\')
      string=Trim(stoich_data%constraints(ic)%tag_species)
      Read (string, Fmt=*, iostat=io) stoich_data%constraints(ic)%species(1), stoich_data%constraints(ic)%species(2)
      If (io/=0) Then
        Write (messages(2),'(a,i2,a)') 'Wrong definition for the involved species of constraint ', ic,&
                                     & ' (type "target_value"). Format MUST be d(species_X)|d(species_Y)'
        Call info(messages,2)
        Call error_stop(' ')
      End If
      Do j=1,2
        If (stoich_data%constraints(ic)%species(j)(1:1) /= 'd' .And. &
          stoich_data%constraints(ic)%species(j)(1:1) /= 'D') Then
          Write (messages(2),'(a,i2,a)') 'Missing "d" or "D" symbol (delta) in the definition of the species ratio&
                                       & for constraints', ic, ' (type "target_value"). Format MUST be&
                                       & d(species_X)|d(species_Y)'
          Call info(messages,2)
          Call error_stop(' ')
        Else
          len=Len_Trim(stoich_data%constraints(ic)%species(j))
          stoich_data%constraints(ic)%species(j)=stoich_data%constraints(ic)%species(j)(2:len)
        End If
        Call remove_symbols(stoich_data%constraints(ic)%species(j),'{[()]}')
        nc=0
        Do k = 1, stoich_data%num_species%value
          If (Trim(stoich_data%constraints(ic)%species(j))==Trim(stoich_data%species(k)%tag)) Then
            If (Trim(stoich_data%species(k)%vartype0)=='fixed') Then
              Write (messages(2),'(3a,i2,a)') 'Species ', Trim(stoich_data%constraints(ic)%species(j)),&
                               & ' (set for constraint ', ic,&
                               & ' of type "', Trim(stoich_data%constraints(ic)%type), &
                               & '") is defined as "fixed" &block_species.&
                               & Constraints MUST ONLY apply to (in)dependent variables'
              Call info(messages,2)
              Call error_stop(' ')
            End If
          Else
            nc=nc+1
          End If
        End Do
        If (nc == stoich_data%num_species%value) Then
          Write (messages(2),'(3a,i2,3a)') 'Species ', Trim(stoich_data%constraints(ic)%species(j)),&
                                & ' (set for constraint ', ic, &
                                & ' of type "', Trim(stoich_data%constraints(ic)%type), &
                                & '") has not been defined in &block_species'
          Call info(messages,2)
          Call error_stop(' ')
        End If
      End Do
      If (Trim(stoich_data%constraints(ic)%species(1))==Trim(stoich_data%constraints(ic)%species(2))) Then
        Write (messages(2),'(a,i2,5a)')  'Constraint ', ic, ' of type "', Trim(stoich_data%constraints(ic)%type), &
                                      & '" involves the same chemical species ', Trim(stoich_data%constraints(ic)%species(1)),&
                                      & '. Please correct the specification for this constraint.'
        Call info(messages,2)
        Call error_stop(' ')
      End If

    End If

  End Subroutine obtain_species_from_contraints

  Subroutine check_stoich_settings(files, eqcm_data, stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives for stoichimetric analysis
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(eqcm_type),      Intent(InOut) :: eqcm_data 
    Type(stoich_type),    Intent(InOut) :: stoich_data

    Integer(Kind=wi)  ::  i, j, length
    Integer(Kind=wi)  ::  ic_var, ic_dep
    Logical  ::  error
    Character(Len=256)  :: messages(3)
    Character(Len=64 )  :: error_block_species
    Character(Len=64 )  :: error_set_eqcm

    error_set_eqcm = '***ERROR in file '//Trim(files(FILE_SET_EQCM)%filename)//' -'
    error_block_species = '***ERROR in &block_species of file '//Trim(files(FILE_SET_EQCM)%filename)

    error=.False.

    ! check directive process
    If (eqcm_data%process%fread) Then
      If (eqcm_data%process%fail) Then
        Write (messages(1),'(2a)')  Trim(error_set_eqcm), ' Wrong (or missing) specification for directive "process"'
        Call error_stop(messages(1))
      Else
        If (Trim(eqcm_data%analysis%type) /= 'model_disordered_system') Then
          If (Trim(eqcm_data%process%type) /= 'electrodeposition' .And. &
            Trim(eqcm_data%process%type) /= 'intercalation') Then
            Write (messages(1),'(2a)') Trim(error_set_eqcm), ' Wrong specification for directive "process". It should&
                                     & be either "electrodeposition" or "intercalation"'
            Call error_stop(messages(1))
          End If
          Call info(' ', 1) 
          Write (messages(1),'(1x,a)')  '===================================================='
          Write (messages(2),'(1x,2a)') '=== Type of physical process: ', Trim(eqcm_data%process%type)
          Write (messages(3),'(1x,a)')  '===================================================='
          Call info(messages, 3)
        Else
          Write (messages(1),'(4a)')  Trim(error_set_eqcm), ' Directive "process" must not be set for "',&
                                    & Trim(eqcm_data%analysis%type), '" type of analysis'
          Call error_stop(messages(1))
        End If
      End If
    Else
      If (Trim(eqcm_data%analysis%type) /= 'model_disordered_system') Then
        Write (messages(1),'(a,1x,2a)')  Trim(error_set_eqcm), Trim(eqcm_data%analysis%type), &
                                ' analysis the requires specification of directive "process"'
        Call error_stop(messages(1))
      Else
        eqcm_data%process%type='intercalation'
      End If
    End If

    ! Check paramentes in blokc_species
    If (.Not. stoich_data%block_species%fread) Then
      Write (messages(1),'(4a)') Trim(error_set_eqcm), ' &block_species, needed for "', Trim(eqcm_data%analysis%type),&
                               &'" analysis, not found'
      error=.True.
      Call info(messages,1)
      Call error_stop(' ')
    Else
      Write (messages(1),'(a)') error_block_species

      Do i=1, stoich_data%num_species%value
        If (stoich_data%species(i)%fail) Then
            Write (messages(2),'(2a)') 'Problems with the specification of chemical species ', Trim(stoich_data%species(i)%tag)
          error=.True.
          Call info(messages,2)
        Else
          Call check_for_symbols(stoich_data%species(i)%tag, '{[()]}', stoich_data%species(i)%fail)
          If (stoich_data%species(i)%fail) Then
            Write (messages(2),'(3a)') 'Species name ', Trim(stoich_data%species(i)%tag), ' contains brackets that MUST be&
                                     & avoided. Please remove all bracket or rename the chemical species.'

            error=.True.
            Call info(messages,2)
          End If

          Call get_word_length(stoich_data%species(i)%tag,length)
          If (length > max_length_name_species) Then
            Write (messages(2),'(3a,i2,a)') 'Species name ', Trim(stoich_data%species(i)%tag), ' exceeds the maximum number of ', &
                                   max_length_name_species, ' set by default. Please rename the species using a shorter string.'

            error=.True.
            Call info(messages,2)
          End If

          If (stoich_data%species(i)%mass <= epsilon(stoich_data%species(i)%mass)) Then
            Write (messages(2),'(3a)') 'Mass for species ', Trim(stoich_data%species(i)%tag), ' must be larger than zero.'
            error=.True.
            Call info(messages,2)
          End If

          If (stoich_data%species(i)%s0 < 0.0_wp) Then
            Write (messages(2),'(3a)') 'Stoichiometry for species ', Trim(stoich_data%species(i)%tag), ' cannot be negative.'
            error=.True.
            Call info(messages,2)
          End If

          If (Abs(stoich_data%species(i)%s0) < epsilon(0.0_wp) .And. Trim(stoich_data%species(i)%vartype0) == 'fixed') Then
            Write (messages(2),'(3a)') 'Species name ', Trim(stoich_data%species(i)%tag), ' is set as "fixed" with zero&
                                     & pristine content. This is a wrong setting for a "fixed" species.'
            error=.True.
            Call info(messages,2)

          End If 


          If (Trim(stoich_data%species(i)%vartype0) /= 'fixed' .And. Trim(stoich_data%species(i)%vartype0) /= 'dependent' &
             .And. Trim(stoich_data%species(i)%vartype0) /= 'independent' ) Then
            Write (messages(2),'(3a)') 'The specification type for species ', Trim(stoich_data%species(i)%tag), ' must be &
                                     &either "fixed", "dependent" or "independent".'
            error=.True.
            Call info(messages,2)
          End If

        End If
      End Do


      If (error) Then
        Call error_stop(' ')
      End If

      ! Analyse if the information provided is sufficient for the requested stoichiometric analysis
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ic_var=0
      ic_dep=0

      Do i=1, stoich_data%num_species%value
       If (Trim(stoich_data%species(i)%vartype0) == 'dependent' .Or. &
          Trim(stoich_data%species(i)%vartype0) == 'independent') Then
         ic_var=ic_var+1
         If (Trim(stoich_data%species(i)%vartype0) == 'dependent') Then
           ic_dep=ic_dep+1
         End If
       End If
      End Do
   
      stoich_data%num_variables=ic_var
      stoich_data%num_indep=ic_var-ic_dep
      stoich_data%num_dep=ic_dep

      If (ic_var==0) Then
        Write (messages(2),'(3a)') 'All stoichiometric species are defined to be "fixed". It is NOT possible to perform the "',&
                                 Trim(eqcm_data%analysis%type), '" analysis'  
        Call info(messages,2)
        Call error_stop(' ')
      End If

      If (eqcm_data%analysis%type == 'model_pristine_sample'  .Or. &
         eqcm_data%analysis%type == 'model_disordered_system' .Or. &
         eqcm_data%analysis%type == 'model_cycled_sample') Then
         If (stoich_data%num_species%value==stoich_data%num_variables) Then
           Write (messages(2),'(a)') 'Modelling analysis necessarily requires the definition of "fixed" species in &block_species'
           Call info(messages,2)
           Call error_stop(' ')
         End If   
         If ((stoich_data%num_species%value-stoich_data%num_variables)==1) Then
           Do i=1, stoich_data%num_species%value
             If (Trim(stoich_data%species(i)%vartype0) == 'fixed') Then
               If (Abs(stoich_data%species(i)%s0_pristine-1.0_wp)>0.0_wp) Then
                 Write (messages(2),'(a)') 'Only one "fixed" species is defined in &block_species. The stoichiometry of this&
                                         & species must be set to 1.0 to build the atomistic model'
                 Call info(messages,2)
                 Call error_stop(' ')
               End If 
             End If 
           End Do
         End If
      End If        


      If (Trim(eqcm_data%process%type) == 'electrodeposition') Then
        If (eqcm_data%analysis%type == 'model_cycled_sample'  .Or. &
            eqcm_data%analysis%type == 'stoichiometry') Then
          If (ic_dep>0) Then
            Write (messages(2),'(a)') 'For electro-deposition processes, the stoichiometric solution is ONLY possible&
                                  & if the  participating species (to be deposited and/or removed) are all&
                                  & defined as "independent" in &Block_species.'
            Write (messages(3),'(a)') 'All independent species must participate&
                                  & with the same number of moles. Species for the substrate must all be defined as "fixed"'
            Call info(messages,3)
            Call error_stop(' ')
          End If

          ! First of all, all independent variables must have the same stoichiometry
          Do i=1, stoich_data%num_species%value-1
            Do j=i+1, stoich_data%num_species%value
              If (stoich_data%species(i)%vartype0=='independent' .And. &
                 stoich_data%species(j)%vartype0=='independent') Then

                If (Abs(stoich_data%species(i)%s0_pristine - stoich_data%species(j)%s0_pristine)>epsilon(1.0_wp)) Then
                  Write (messages(2),'(5a)') 'Subspecies ',  Trim(stoich_data%species(i)%tag), ' and ', &
                                      & Trim(stoich_data%species(j)%tag), ' have different values for the pristine stoichiometry.'
                  Write (messages(3),'(a)') 'Requested analysis for the eletrodeposition of more than one independent subspecies&
                                           & is ONLY possible if they participate with the same amount of moles'
                  Call info(messages,3)
                  Call error_stop(' ')
                End If
              End If
            End Do
          End Do
        End If

      Else If (Trim(eqcm_data%process%type) == 'intercalation') Then
        If (stoich_data%num_species%value == 1 .And. ic_var >0) Then
          Write (messages(2),'(a,i2,a)') 'Intercalation processes require the definition of (fixed) species corresponding&
                                    & to the host material. Please review the settings'
          Call info(messages,2)
          Call error_stop(' ')
        End If
 
        If (eqcm_data%analysis%type == 'model_cycled_sample' .Or. &
            eqcm_data%analysis%type == 'stoichiometry') Then
          If (ic_dep > 2) Then
            Write (messages(2),'(a)') 'The number of dependent stoichiometric variables MUST ALWAYS BE <= 2'
            Call info(messages,2)
            Call error_stop(' ')
          End If

          If (ic_var >= 2 .And. ic_dep < 2) Then
            Write (messages(2),'(a,i2,a)') 'The number of dependent/independent species in &block_species is equal to ', ic_var,&
                                          & '. Nevertheless, it is necessary to set two of these species as "dependent".'
            Call info(messages,2)
            Call error_stop(' ')
          End If
        End If
      End If

      ! allocate relevant matrices only if the stoichiometric problem will be solved
      If (eqcm_data%analysis%type == 'stoichiometry' .Or. &
         eqcm_data%analysis%type == 'model_cycled_sample') Then 
         Call stoich_data%init_species_arrays('matrices')
      End If

    End If

  End Subroutine check_stoich_settings

End Module stoichiometry
