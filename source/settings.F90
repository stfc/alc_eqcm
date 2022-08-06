!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Module that:
! - reads SET_EQCM and defines the settings for EQCM analysis 
! - checks correctness of directives
! - checks feasibility of the setup according to the experimental values (read from DATA_EQCM)
! - checks feasibility of stoichiometric calculation 
! - checks input directives to build atomistic models 
! - checks input directives to build scripts for HPC execution 
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author        - i.scivetti May     2020
! Contribution  - i.scivetti Dec     2020 (Stoichiometric settings)
! Contribution  - i.scivetti Mar     2021 (Atomistic models) + rearrangement
! Contribution  - i.scivetti Apr-Jul 2021 Reading directives for building files for DFT simulation
! Contribution  - i.scivetti May     2021 Reading directives for building scripts for HPC execution 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module settings 

  Use atomistic_models,  Only : check_atomic_settings, &
                                max_species,&
                                max_intra,  &
                                min_intra,  &
                                min_inter,  &
                                model_type
  Use constants,         Only : Bohr_to_A,  &
                                chemsymbol, & 
                                NPTE  
  Use electrode,         Only : electrode_type
  Use eqcm,              Only : eqcm_type
  Use fileset,           Only : file_type, &
                                FILE_SET_EQCM, &  
                                FILE_DATA_EQCM,&
                                FILE_OUT_EQCM, & 
                                FILE_RECORD_MODELS,& 
                                FOLDER_ANALYSIS, &
                                refresh_out_eqcm
  Use fft,               Only : fft_type
  Use filtering,         Only : filter_type
  Use hpc_setup,         Only : hpc_type, &
                                check_hpc_settings
  Use numprec,           Only : wi, &
                                wp
  Use process_data,      Only : capital_to_lower_case, &
                                check_for_rubbish, &
                                get_word_length
  Use sauerbrey,         Only : sauer_type
  Use simulation_setup,  Only : check_simulation_settings     
  Use simulation_tools,  Only : simul_type
  Use system,            Only : system_type
  Use stoichiometry,     Only : stoich_type, &
                                check_components_species,&
                                check_constraints_settings,&
                                check_stoich_settings, &
                                ds_default,&
                                efficiency_default 
  Use unit_output,       Only : error_stop,&
                                info

  Implicit None
  Private

  Public :: check_feasibility, read_settings, check_all_settings

Contains

  Subroutine duplication_error(directive)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Aborts execution when duplication for
   ! a directive is found
   !
   ! author - i. scivetti  June 2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: directive

    Character(Len=256)  :: message

    Write (message,'(4a)') '***ERROR - Directive "', Trim(directive), '" is duplicated!'
    Call error_stop(message)

  End Subroutine duplication_error  

  Subroutine set_read_status(word, io, fread, fail, string)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to:
    !  - prevent duplication
    !  - define input directive is read by setting fread=.True. 
    !  - test if there was a problem with reading a directive, indicated by io/=0. This sets fail=.True.
    !
    ! author    - i.scivetti June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: word
    Integer(Kind=wi), Intent(In   ) :: io
    Logical,          Intent(  Out) :: fread 
    Logical,          Intent(InOut) :: fail
    Character(Len=*), Optional, Intent(InOut) :: string

    If (fread)then
      Call duplication_error(word)
    Else
      fread=.True.
      If (io /= 0) Then
        fail=.True.
      End If
    End If

    If (present(string)) then
      Call capital_to_lower_case(string)
    End If

  End Subroutine set_read_status 


  Subroutine read_settings(files, eqcm_data, fft_data, filter_data, system_data,&
                         & electrode_data, sauer_data, stoich_data, model_data, &
                         & simulation_data, hpc_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read EQCM settings from SET_EQCM file.
    ! Lines starting with # are ignored and assumed as comments. 
    ! If a directive is identified during the reading of the file, subroutine "set_read_fail" 
    ! assigns fread=.True. On the contrary, the subroutine assigns fail=.True. (fail=.False.) 
    ! if the format/syntax for the directive is correct (incorrect)
    ! If the directive is repeated the execution is aborted via subroutine duplication 
    ! 
    ! author        - i.scivetti May 2020
    ! contribution  - i.scivetti Nov 2020      
    ! contribution  - i.scivetti Jan 2021
    ! contribution  - i.scivetti Apr 2021 (setting simulation files) 
    ! contribution  - i.scivetti May 2021 (setting script files for submission of jobs) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),       Intent(InOut) :: files(:)
    Type(eqcm_type),       Intent(InOut) :: eqcm_data 
    Type(fft_type),        Intent(InOut) :: fft_data
    Type(filter_type),     Intent(InOut) :: filter_data
    Type(system_type),     Intent(InOut) :: system_data
    Type(electrode_type),  Intent(InOut) :: electrode_data
    Type(sauer_type),      Intent(InOut) :: sauer_data
    Type(stoich_type),     Intent(InOut) :: stoich_data 
    Type(model_type),      Intent(InOut) :: model_data 
    Type(simul_type),      Intent(InOut) :: simulation_data 
    Type(hpc_type),        Intent(InOut) :: hpc_data 
 
    Logical            :: safe
    Character(Len=256) :: word
    Integer(Kind=wi)   :: i, j
    Integer(Kind=wi)   :: length, io, iunit
  
    Character(Len=256)  :: message

    Character(Len=32 )  :: eqcm_file
    Character(Len=32 )  :: set_error

    eqcm_file = Trim(files(FILE_SET_EQCM)%filename)
    set_error = '***ERROR -'

    ! Initialise relevant arrays
    Call eqcm_data%init_input_variables()
    Call sauer_data%init_input_variables()
    Call model_data%init_input_variables()
 
    ! Open the SET_EQCM file with EQCM settings
    Inquire(File=files(FILE_SET_EQCM)%filename, Exist=safe)
    
    If (.not.safe) Then
      Call info(' ', 1)
      Write (message,'(4(1x,a))') Trim(set_error), 'File', Trim(eqcm_file), '(settings for EQCM analysis) not found'
      Call error_stop(message)
    Else
      Open(Newunit=files(FILE_SET_EQCM)%unit_no, File=Trim(eqcm_file), Status='old')
      iunit=files(FILE_SET_EQCM)%unit_no 
    End If

     Read (iunit, Fmt=*, iostat=io) word
     ! If nothing is found, complain and abort
     If (is_iostat_end(io)) Then
       Write (message,'(3(1x,a))') Trim(set_error), Trim(eqcm_file), 'file seems to be empty?. Please check'
       Call error_stop(message)
     End If
     ! Check header has "#" as the first character 
     If (word(1:1)/='#') Then
       Write (message,'(4(1x,a))') Trim(set_error), 'Heading comment in file', Trim(eqcm_file), & 
                                  'is required and MUST be preceded with the symbol "#"'
       Call error_stop(message)
     End If

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Exit
      end If
      Call check_for_rubbish(iunit, Trim(eqcm_file)) 
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
        ! Do nothing if line is a comment of we have an empty line
        Read (iunit, Fmt=*, iostat=io) word

      !!!!!!!!!!!!!!!!!!!!!!! 
      ! EQCM directives
      !!!!!!!!!!!!!!!!!!!!!!! 
      Else If (word(1:length) == 'analysis') Then
        Read (iunit, Fmt=*, iostat=io) word, eqcm_data%analysis%type
        Call set_read_status(word, io, eqcm_data%analysis%fread, eqcm_data%analysis%fail, eqcm_data%analysis%type)

      Else If (word(1:length) == 'software') Then
        Read (iunit, Fmt=*, iostat=io) word, eqcm_data%software%type
        Call set_read_status(word, io, eqcm_data%software%fread, eqcm_data%software%fail, eqcm_data%software%type)
  
      Else If (word(1:length) == 'scan_rate') Then
        Read (iunit, Fmt=*, iostat=io) word, eqcm_data%scan_rate%value(1), (eqcm_data%scan_rate%units(j), j=1,2)
        Call set_read_status(word, io, eqcm_data%scan_rate%fread, eqcm_data%scan_rate%fail)
        Call capital_to_lower_case(eqcm_data%scan_rate%units(1))
        Call capital_to_lower_case(eqcm_data%scan_rate%units(2))

      Else If (word(1:length) == 'cycles') Then
        Read (iunit, Fmt=*, iostat=io) word, (eqcm_data%range_cycles%value(i), i=1,2)
        Call set_read_status(word, io, eqcm_data%range_cycles%fread, eqcm_data%range_cycles%fail)

      Else If (word(1:length) == 'voltage_range') Then
        Read (iunit, Fmt=*, iostat=io) word,&
                                    & (eqcm_data%voltage_range%value(i), i=1,2), eqcm_data%voltage_range%units(1)
        Call set_read_status(word, io, eqcm_data%voltage_range%fread, eqcm_data%voltage_range%fail)
        Call capital_to_lower_case(eqcm_data%voltage_range%units(1))

      Else If (word(1:length) == 'current_offset') Then
        Read (iunit, Fmt=*, iostat=io) word, eqcm_data%current_offset%stat
        Call set_read_status(word, io, eqcm_data%current_offset%fread, eqcm_data%current_offset%fail)

      Else If (word(1:length) == 'process') Then 
        Read (iunit, Fmt=*, iostat=io) word, eqcm_data%process%type
        Call set_read_status(word, io, eqcm_data%process%fread, eqcm_data%process%fail, eqcm_data%process%type)

      Else If (word(1:length) == 'efficiency') Then
        Read (iunit, Fmt=*, iostat=io) word, eqcm_data%efficiency%value
        Call set_read_status(word, io, eqcm_data%efficiency%fread, eqcm_data%efficiency%fail)

      !!!!!!!!!!!!!!!!!!!!!!! 
      ! Filter directives
      !!!!!!!!!!!!!!!!!!!!!!! 
      Else If (word(1:length) == 'filter_cutoff') Then
        Read (iunit, Fmt=*, iostat=io) word, filter_data%cutoff%value, filter_data%cutoff%units
        Call set_read_status(word, io, filter_data%cutoff%fread, filter_data%cutoff%fail)
        Call capital_to_lower_case(filter_data%cutoff%units)

      Else If (word(1:length) == 'endpoints_mass_frequency') Then
        Read (iunit, Fmt=*, iostat=io) word, fft_data%end_mass_frequency%value
        Call set_read_status(word, io, fft_data%end_mass_frequency%fread, fft_data%end_mass_frequency%fail)

      Else If (word(1:length) == 'endpoints_current') Then
        Read (iunit, Fmt=*, iostat=io) word, fft_data%end_current%value
        Call set_read_status(word, io, fft_data%end_current%fread, fft_data%end_current%fail)

      !!!!!!!!!!!!!!!!!!!!!!! 
      ! Directives for conversion from Mass-Frequency to mass
      !!!!!!!!!!!!!!!!!!!!!!! 
      Else If (word(1:length) == 'sauerbrey') Then
        Read (iunit, Fmt=*, iostat=io) word, sauer_data%factor%value(1), (sauer_data%factor%units(j), j=1,3)
        Call set_read_status(word, io, sauer_data%factor%fread, sauer_data%factor%fail)
        Do j=1,3
          Call capital_to_lower_case(sauer_data%factor%units(j))
        End Do

      Else If (word(1:length) == 'v_to_hz') Then
        Read (iunit, Fmt=*, iostat=io) word, eqcm_data%V_to_Hz%value
        Call set_read_status(word, io, eqcm_data%V_to_Hz%fread, eqcm_data%V_to_Hz%fail)

      !!!!!!!!!!!!!!!!!!!!!!! 
      ! System specifications
      !!!!!!!!!!!!!!!!!!!!!!! 
      Else If (word(1:length) == 'electrode_area') Then
        Read (iunit, Fmt=*, iostat=io) word, electrode_data%area_geom%value, electrode_data%area_geom%units 
        Call set_read_status(word, io, electrode_data%area_geom%fread, electrode_data%area_geom%fail)
        Call capital_to_lower_case(electrode_data%area_geom%units)

      Else If (word(1:length) == 'area_scale') Then
        Read (iunit, Fmt=*, iostat=io) word, electrode_data%area_scale%value
        Call set_read_status(word, io, electrode_data%area_scale%fread, electrode_data%area_scale%fail)

      Else If (word(1:length) == 'electrode_mass') Then
        Read (iunit, Fmt=*, iostat=io) word, electrode_data%mass%value, electrode_data%mass%units
        Call set_read_status(word, io, electrode_data%mass%fread, electrode_data%mass%fail)
        Call capital_to_lower_case(electrode_data%mass%units)

      Else If (word(1:length) == 'electrode_mole') Then
        Read (iunit, Fmt=*, iostat=io) word, electrode_data%mole%value
        Call set_read_status(word, io, electrode_data%mole%fread, electrode_data%mole%fail)

      Else If (word(1:length) == 'quartz_freq') Then
        Read (iunit, Fmt=*, iostat=io) word, system_data%quartz_freq%value, system_data%quartz_freq%units
        Call set_read_status(word, io, system_data%quartz_freq%fread, system_data%quartz_freq%fail)
        Call capital_to_lower_case(system_data%quartz_freq%units)

      !!!!!!!!!!!!!!!!!!!!!!!
      ! Stoichiometry specifications
      !!!!!!!!!!!!!!!!!!!!!!!
      Else If (word(1:length) == '&block_species') Then
        Read (iunit, Fmt=*, iostat=io) stoich_data%block_species%type
        Call set_read_status(word, io, stoich_data%block_species%fread, stoich_data%block_species%fail)
        ! Read information inside the block
        Call read_block_species(iunit, stoich_data)  

      Else If (word(1:length) == '&block_constraints_species') Then
        If (.Not.stoich_data%block_species%fread) Then
          Write (message,'(3(1x,a))') Trim(set_error), '&block_constraints_species must be specified after', &
                                     '&block_species. Have you specified &block_species?'
          Call error_stop(message) 
        End If
        Read (iunit, Fmt=*, iostat=io) stoich_data%block_constraints%type
        Call set_read_status(word, io, stoich_data%block_constraints%fread, stoich_data%block_constraints%fail)
        ! Read information inside the block
        Call read_block_constraints_species(iunit, stoich_data)  

      Else If (word(1:length) == 'delta_stoich') Then
        Read (iunit, Fmt=*, iostat=io) word, stoich_data%discretization%value
        Call set_read_status(word, io, stoich_data%discretization%fread, stoich_data%discretization%fail)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Specifications for building atomistic models
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Else If (word(1:length) == '&block_species_components') Then
        If (.Not. eqcm_data%analysis%fread) Then
          Write (message,'(3(1x,a))') Trim(set_error), 'Directive "Analysis" must be specified before', &
                                     '&block_species_components'
          Call error_stop(message) 
        End If
        If (.Not.stoich_data%block_species%fread) Then
          Write (message,'(3(1x,a))') Trim(set_error), '&block_species_components must be specified after', &
                                     '&block_species. Have you specified &block_species?'
          Call error_stop(message) 
        End If
        Read (iunit, Fmt=*, iostat=io) stoich_data%block_species_components%type
        Call set_read_status(word, io, stoich_data%block_species_components%fread, stoich_data%block_species_components%fail)
        ! Read information inside the block
        Call read_block_species_components(iunit, stoich_data, eqcm_data%analysis%type)  

      Else If (word(1:length) == '&block_input_composition') Then
        If (.Not.stoich_data%block_species_components%fread) Then
          Write (message,'(3(1x,a))') Trim(set_error), '&block_input_composition must be specified after', &
                                     '&block_species_components. Have you specified &block_species_components at all?'
          Call error_stop(message)
        End If
        Read (iunit, Fmt=*, iostat=io) model_data%block_input_composition%type
        Call set_read_status(word, io, model_data%block_input_composition%fread, model_data%block_input_composition%fail)
        ! Read information inside the block
        Call read_block_input_composition(iunit, stoich_data, model_data)

      Else If (word(1:length) == '&block_input_cell') Then
        Read (iunit, Fmt=*, iostat=io) model_data%block_input_cell%type
        Call set_read_status(word, io, model_data%block_input_cell%fread, model_data%block_input_cell%fail)
        ! Read information inside the block
        Call read_block_input_cell(iunit, model_data)  

      Else If (word(1:length) == 'input_model_format') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%input_model_format%type
        Call set_read_status(word, io, model_data%input_model_format%fread, model_data%input_model_format%fail, &
                           & model_data%input_model_format%type)

      Else If (word(1:length) == 'output_model_format') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%output_model_format%type
        Call set_read_status(word, io, model_data%output_model_format%fread, model_data%output_model_format%fail, &
                           & model_data%output_model_format%type)

      Else If (word(1:length) == 'delta_space') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%delta_space%value, model_data%delta_space%units 
        Call set_read_status(word, io, model_data%delta_space%fread, model_data%delta_space%fail)
        Call capital_to_lower_case(model_data%delta_space%units)

      Else If (word(1:length) == 'scale_cell') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%scale_cell%value 
        Call set_read_status(word, io, model_data%scale_cell%fread, model_data%scale_cell%fail)

      Else If (word(1:length) == 'distance_cutoff') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%distance_cutoff%value, model_data%distance_cutoff%units 
        Call set_read_status(word, io, model_data%distance_cutoff%fread, model_data%distance_cutoff%fail)
        Call capital_to_lower_case(model_data%distance_cutoff%units)

      Else If (word(1:length) == 'deposition_level') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%deposition_level%value, model_data%deposition_level%units 
        Call set_read_status(word, io, model_data%deposition_level%fread, model_data%deposition_level%fail)
        Call capital_to_lower_case(model_data%deposition_level%units)

      Else If (word(1:length) == 'repeat_input_model') Then 
        Read (iunit, Fmt=*, iostat=io) word, (model_data%repeat_input_model%value(j), j=1,3)
        Call set_read_status(word, io, model_data%repeat_input_model%fread, model_data%repeat_input_model%fail)

      Else If (word(1:length) == 'rotate_species') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%rotate_species%stat
        Call set_read_status(word, io, model_data%rotate_species%fread, model_data%rotate_species%fail)

      Else If (word(1:length) == 'stoichiometry_error') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%stoichiometry_error%value
        Call set_read_status(word, io, model_data%stoichiometry_error%fread, model_data%stoichiometry_error%fail)

      Else If (word(1:length) == 'targeted_num_models') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%targeted_num_models%value
        Call set_read_status(word, io, model_data%targeted_num_models%fread, model_data%targeted_num_models%fail)

      Else If (word(1:length) == 'multiple_input_atoms') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%multiple_input_atoms%value
        Call set_read_status(word, io, model_data%multiple_input_atoms%fread, model_data%multiple_input_atoms%fail)

      Else If (word(1:length) == 'normal_vector') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%normal_vector%type
        Call set_read_status(word, io, model_data%normal_vector%fread, model_data%normal_vector%fail,&
                           & model_data%normal_vector%type)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Specifications of settings to build simulation files 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Else If (word(1:length) == '&block_simulation_settings') Then
        If (.Not.model_data%block_input_composition%fread) Then
          Write (message,'(2(1x,a))') Trim(set_error), '"&block_simulation_settings" must be specified after&
                                    & "&block_input_composition"'
          Call error_stop(message)
        End If
 
        Read (iunit, Fmt=*, iostat=io) word 
        If (simulation_data%generate) Then
          Call duplication_error(word)
        End If
        simulation_data%generate=.True.
        simulation_data%total_tags=stoich_data%total_tags
      
        ! Asssign variables before reading simulation settings
        Do i=1, simulation_data%total_tags
          simulation_data%component(i)%tag=model_data%block_component%tag(i)
          simulation_data%component(i)%element=model_data%block_component%element(i)
          simulation_data%component(i)%atomic_number=model_data%block_component%atomic_number(i)
        End Do 
        ! Now, it is ready to read information inside &block_simulation
        Call read_simulation_settings(iunit, simulation_data)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Specifications of HPC settings 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Else If (word(1:length) == '&block_hpc_settings') Then
        Read (iunit, Fmt=*, iostat=io) word 
        If (hpc_data%generate) Then
          Call duplication_error(word)
        End If
        hpc_data%generate=.True.
        Call hpc_data%init_input_variables()

        ! Now, it is ready to read information inside &block_hpc_settings
        Call read_hpc_settings(iunit, hpc_data)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Directive not recognised. Inform and kill 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      Else
        If (word(1:1)=='&') Then
          Write (message,'(1x,6a)') Trim(set_error), ' unknown directive found in ', Trim(eqcm_file), &
                                  &' file: ', Trim(word),'. Do you use "&" to define a block? If so,&
                                  & make sure the block is located in the right place, as it might be&
                                  & a sub-block actually. Please also check for the right syntax.'
        Else
          Write (message,'(5(1x,a))') Trim(set_error), 'unknown directive found in', Trim(eqcm_file), &
                                  &'file: ', Trim(word)
        End If 
        Call error_stop(message)
      End If

    End Do

    ! Close file
    Close(files(FILE_SET_EQCM)%unit_no)


  End Subroutine read_settings

  Subroutine check_all_settings(files, eqcm_data, fft_data, filter_data, system_data, electrode_data, &
                               &sauer_data, stoich_data, model_data, simulation_data, hpc_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the correctness of all directive read in file SET_EQCM
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut) :: files(:)
    Type(eqcm_type),      Intent(InOut) :: eqcm_data 
    Type(fft_type),       Intent(InOut) :: fft_data
    Type(filter_type),    Intent(InOut) :: filter_data
    Type(system_type),    Intent(InOut) :: system_data
    Type(electrode_type), Intent(InOut) :: electrode_data
    Type(sauer_type),     Intent(InOut) :: sauer_data
    Type(stoich_type),    Intent(InOut) :: stoich_data 
    Type(model_type),     Intent(InOut) :: model_data 
    Type(simul_type),     Intent(InOut) :: simulation_data
    Type(hpc_type),       Intent(InOut) :: hpc_data
 
    Integer(Kind=wi)   ::  j
 
    Character(Len=256) :: message
    Character(Len=256) :: messages(12)
    Character(Len=64 ) :: error_set_eqcm 
    Logical            :: fprint, error

    error_set_eqcm = '***ERROR in file '//Trim(files(FILE_SET_EQCM)%filename)//' -'   
  
    ! EQCM directives
    !!!!!!!!!!!!!!!!!!!!!!! 
    ! Check the type of analysis selected
    If (eqcm_data%analysis%type /= 'spectra'                .And.  &
       eqcm_data%analysis%type /= 'print_eqcm_raw'          .And.  &
       eqcm_data%analysis%type /= 'print_eqcm_filter'       .And.  &
       eqcm_data%analysis%type /= 'mass_calibration'        .And.  &
       eqcm_data%analysis%type /= 'massogram'               .And.  &
       eqcm_data%analysis%type /= 'characterization'        .And.  &
       eqcm_data%analysis%type /= 'stoichiometry'           .And.  &
       eqcm_data%analysis%type /= 'model_pristine_sample'   .And.  &
       eqcm_data%analysis%type /= 'model_cycled_sample'     .And.  &
       eqcm_data%analysis%type /= 'model_disordered_system' .And.  &
       eqcm_data%analysis%type /= 'hpc_simulation_files' ) Then 
      Write (messages(1),'(1x,a)')  'Possibe settings for option "analysis" are:'
      Write (messages(2),'(1x,a)')  '- Spectra'
      Write (messages(3),'(1x,a)')  '- Print_EQCM_raw'
      Write (messages(4),'(1x,a)')  '- Print_EQCM_filter'
      Write (messages(5),'(1x,a)')  '- Mass_Calibration'
      Write (messages(6),'(1x,a)')  '- Massogram'
      Write (messages(7),'(1x,a)')  '- Characterization'
      Write (messages(8),'(1x,a)')  '- Stoichiometry'  
      Write (messages(9),'(1x,a)')  '- Model_pristine_sample'
      Write (messages(10),'(1x,a)') '- Model_cycled_sample' 
      Write (messages(11),'(1x,a)') '- Model_disordered_system' 
      Write (messages(12),'(1x,a)') '- hpc_simulation_files'
      Call info(messages,12)
      Call info(' ',1)
      Write (message,'(2(1x,a))')  Trim(error_set_eqcm), 'No (or wrong) specification for directive "analysis"' 
      Call error_stop(message) 
    Else
      Call info(' ',1)
      Write (message,'(1x,2a)') 'Requested analysis: ', Trim(eqcm_data%analysis%type)
      Call info(message,1)
              Call info(' ----------------------------------------------------------------------------------------', 1)
      If (eqcm_data%analysis%type == 'spectra' ) Then
        Write (messages(1),'(1x,a)')  'This option computes the FFT for current and/or mass density. Magnitudes are reported'
        Write (messages(2),'(1x,a)')  'for the different components of the reciprocal space in units of 1/V (if the elapsed'
        Write (messages(3),'(1x,a)')  'is not recorded) or Hz (is time is recorded). The user is referred to the manual for'
        Write (messages(4),'(1x,a)')  'details. This feature is convenient to identify those regions of the reciprocal space'
        Write (messages(5),'(1x,a)')  'that might be responsible for the presence of strong noise in the EQCM data.'
        Call info(messages,5)
      Else If (eqcm_data%analysis%type == 'print_eqcm_raw' ) Then
        Write (messages(1),'(1x,a)')  'This option prints the raw current and/or mass density reported in file DATA_EQCM'
        Write (messages(2),'(1x,a)')  'for the requested range of CV cycles, as specified by directive "cycles". If "cycles"'
        Write (messages(3),'(1x,a)')  'is not set, the analysis is carried out for all the CV cycles found in the set of data.'
        Write (messages(4),'(1x,a)')  'This capability is convenient to visualise experimental data and identify any problem'
        Write (messages(5),'(1x,a)')  'such as the presence of strong noise.'
        Call info(messages,5)
      Else If (eqcm_data%analysis%type == 'print_eqcm_filter' ) Then
        Write (messages(1),'(1x,a)')  'This option filters and prints current and/or mass density reported in file DATA_EQCM'
        Write (messages(2),'(1x,a)')  'for the requested range of CV cycles, as specified by directive "cycles". If "cycles"'
        Write (messages(3),'(1x,a)')  'is not set, the analysis is carried out for all the CV cycles found in the set of data.'
        Write (messages(4),'(1x,a)')  'Filtering is conducted by applying a low-pass Gaussian filter unsing a cutoff value'
        Write (messages(5),'(1x,a)')  'provided by the user via directive "filter_cutoff". See manual for more information.' 
        Call info(messages,5)
      Else If (eqcm_data%analysis%type == 'mass_calibration' ) Then
        Write (messages(1),'(1x,a)')  'This option prints the mass-frequency and (computed) accumulated charge for each point'
        Write (messages(2),'(1x,a)')  'of the CV data. If directive "cycles" is not set, the analysis is carried out for all'
        Write (messages(3),'(1x,a)')  'CV cycles found in the set of the data. Plots of mass-frequency as a function of charge'
        Write (messages(4),'(1x,a)')  'are often used to fit calibration coefficients for EQCM devices (see manual).' 
        Call info(messages,4)
      Else If (eqcm_data%analysis%type == 'massogram' ) Then
        Write (messages(1),'(1x,a)')  'This option prints the rate of mass density as a function of the applied voltage'
        Write (messages(2),'(1x,a)')  'for the requested CV cycle range, as specified by directive "cycles". If "cycles"'
        Write (messages(3),'(1x,a)')  'is not set, the analysis is carried out for all CV cycles found in the set of data.'
        Write (messages(4),'(1x,a)')  'Massograms provide a sensitive tool to identify the presence of side reactions' 
        Call info(messages,4)
      Else If (eqcm_data%analysis%type == 'characterization' ) Then
        Write (messages(1),'(1x,a)')  'This option computes relevant quantities associated with the redox process from the'
        Write (messages(2),'(1x,a)')  'data reported in file DATA_EQCM. The range of cycles for analysis is given by directive'
        Write (messages(3),'(1x,a)')  '"cycles". If "cycles" is not set, all cycles found in the set of data are considered.'
        Call info(messages,3)
      Else If (eqcm_data%analysis%type == 'stoichiometry' ) Then
        Write (messages(1),'(1x,a)')  'This option computes the stoichiometric related solutions from the accumulated charge'
        Write (messages(2),'(1x,a)')  'and mass, as well as the variables defined in &block_species. The outcome depends on'
        Write (messages(3),'(1x,a)')  'the type of process under consideration: "intercalation" or "electrodeposition".'
        Call info(messages,3)
      Else If (eqcm_data%analysis%type == 'model_cycled_sample' ) Then
        Write (messages(1),'(1x,a)')  'This option generates atomistic models compatible with EQCM data and the definition'
        Write (messages(2),'(1x,a)')  'of the participating species in &block_species.'
        Call info(messages,2)
      Else If (eqcm_data%analysis%type == 'model_pristine_sample' ) Then
        Write (messages(1),'(1x,a)')  'This option generates a single atomistic model, whose stoichiometry will target the'
        Write (messages(2),'(1x,a)')  'stoichiometric values defined in &block_species. The term pristine refers to the'
        Write (messages(3),'(1x,a)')  'system before EQCM cycling. Thus, the generated model does not depend on the EQCM data.'
        Write (messages(4),'(1x,a)')  'If defined, directives related to EQCM specification and data analysis will be ignored.'
        Call info(messages,4)
      Else If (eqcm_data%analysis%type == 'model_disordered_system' ) Then
        Write (messages(1),'(1x,a)')  'This option generates a single atomistic model of a disordered system, whose stoichiometry'
        Write (messages(2),'(1x,a)')  'will be set to approximate the values defined in &block_species. The generated'
        Write (messages(3),'(1x,a)')  'model does not depend on the EQCM data. If defined, directives related to EQCM'
        Write (messages(4),'(1x,a)')  'specification and data analysis will be ignored. See manual for more information'
        Call info(messages,4)
      Else If (eqcm_data%analysis%type == 'hpc_simulation_files' ) Then
        Write (messages(1),'(1x,a)')  'This option generates (or regenerates):'
        Write (messages(2),'(1x,a)')  ' - input files for simulation (if "&block_simulation_settings" is present)'
        Write (messages(3),'(1x,a)')  ' - script files for HPC submission (if "&block_hpc_settings" is present)'
        Write (messages(4),'(1x,2a)') 'without the need to generate the atomistic models again. This analysis&
                                    & strictly needs the file ', Trim(files(FILE_RECORD_MODELS)%filename)
        Write (messages(5),'(1x,a)')  'which is created after the generation of atomistic models and contains&
                                    & all the relevant settings.'
        Write (messages(6),'(1x,a)')  'The output with information of the atomistic models has been saved to&
                                    & the OUT_EQCM_BACKUP file.'
        Write (messages(7),'(1x,3a)') 'Thus, the new ', Trim(files(FILE_OUT_EQCM)%filename), ' file must be&
                                    & interpreted as a complement of the OUT_EQCM_BACKUP file.' 
        Call info(messages,7)
      End If
              Call info(' ----------------------------------------------------------------------------------------', 1)
      Call info(' ', 1)
    End If

    fprint=.True.
    If (eqcm_data%analysis%type== 'model_pristine_sample'    .Or. &
        eqcm_data%analysis%type== 'model_disordered_system'  .Or. &
        eqcm_data%analysis%type== 'hpc_simulation_files') Then
      fprint=.False.
    End If
    
    ! Create folder ANALYSIS
    If (eqcm_data%analysis%type /= 'hpc_simulation_files'    .And. &
        eqcm_data%analysis%type /= 'model_disordered_system' .And. &
        eqcm_data%analysis%type /= 'model_pristine_sample') Then
       Call execute_command_line('[ ! -d '//Trim(FOLDER_ANALYSIS)//' ] && '//'mkdir '//Trim(FOLDER_ANALYSIS)) 
    End If

    If (fprint) Then
      Write (messages(1),'(1x,a)') 'Relevant EQCM settings'
      Write (messages(2),'(1x,a)') '======================'
      Call info(messages,2)
    End If


    ! Check the type of software used for collecting EQCM data
    If (eqcm_data%software%fread) Then
      If (eqcm_data%software%type/='qcm200'  .And. &
         eqcm_data%software%type/='metrohm'  .And. &
         eqcm_data%software%type/='ch-inst') Then
         eqcm_data%software%warn=.True.
      End If
    Else
      eqcm_data%software%type='Not specified'
    End If

    Write (message,'(1x,2a)') 'EQCM software:      ', Trim(eqcm_data%software%type)
    If (fprint) Then
      Call info(message,1)
    End If

    ! Voltage scan rate
    If (eqcm_data%scan_rate%fread) Then 
      If (eqcm_data%scan_rate%fail) Then 
        If (fprint)Then
          Call info(' ',1)
        End If
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong specification for directive "scan_rate"'
        Call error_stop(message)
      Else     
        If (eqcm_data%scan_rate%value(1) <= epsilon(eqcm_data%scan_rate%value(1))) Then
          If (fprint) Then
            Call info(' ',1)
          End If
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), '"scan_rate" should be larger than zero.'
          Call error_stop(message)
        End If
        error=.False.
        If (eqcm_data%scan_rate%units(1) == 'v') Then
          eqcm_data%scan_rate%convert = 1.0_wp
        Else If (eqcm_data%scan_rate%units(1) == 'mv') Then
          eqcm_data%scan_rate%convert = 0.001_wp 
        Else
          error=.True.
        End If
        If (eqcm_data%scan_rate%units(2) /= 's-1') Then
          error=.True.
        End If

        If(error)Then
          If (fprint) Then
            Call info(' ',1)
          End If
          Write (messages(1),'(5(1x,a))') Trim(error_set_eqcm), 'Invalid units for "scan_rate" set to: ',&
                                         & (Trim(eqcm_data%scan_rate%units(j)), j=1,2)
          Write (messages(2),'(1x,a)')   'Valid units must be "mV s-1" or "V s-1"'
          Call info(messages, 2)
          Call error_stop(' ')
        End If

        eqcm_data%scan_rate%value(1) = eqcm_data%scan_rate%convert * eqcm_data%scan_rate%value(1)
        Write (message,'(1x,a,f8.4,a)') 'Scan rate:        ', eqcm_data%scan_rate%value(1), ' [V/s]' 
        If (fprint) Then
          Call info(message,1)
        End If
      End If
    End If

    ! Cycling specification
    If (eqcm_data%range_cycles%fread) Then 
      If (eqcm_data%range_cycles%fail) Then 
        If (fprint) Then
          Call info(' ',1)
        End If
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong specification for the range of CV cycles to be &
                                                         &considered. Check specification for "cycles" directive'
        Call error_stop(message)
      Else
        If (eqcm_data%range_cycles%value(1) > eqcm_data%range_cycles%value(2) ) Then
          If (fprint) Then
            Call info(' ',1)
          End If
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Flag "cycles": lowest bound for the cycle range &
                                                           &MUST be set before the upmost bound'
          Call error_stop(message)
        Else If (eqcm_data%range_cycles%value(1)<=0 .Or. eqcm_data%range_cycles%value(2)<=0) Then
          If (fprint) Then
            Call info(' ',1)
          End If
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Flag "cycles": Cycle numbers must be larger than zero'
          Call error_stop(message)
        End If 
      End If
      If (eqcm_data%range_cycles%value(1)==eqcm_data%range_cycles%value(2)) Then
        Write (message,'(1x,a,i3,a)') 'CV cycle range:     Only cycle ', eqcm_data%range_cycles%value(1), &
                                      ' has been selected for the requested analysis'
        If (fprint) Then
          Call info(message,1)
        End If
      Else 
        Write (message,'(1x,2(a,i3))') 'CV cycle range:    ', eqcm_data%range_cycles%value(1), ' to ',&
                                      & eqcm_data%range_cycles%value(2)
        If (fprint) Then
          Call info(message,1)
        End If

      End If
    Else
      eqcm_data%range_cycles%warn=.True.
    End If    

    If (eqcm_data%voltage_range%fread) Then 
      If (eqcm_data%voltage_range%fail) Then 
        If (fprint) Then
          Call info(' ',1)
        End If
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong specification for the voltage_range directive'
        Call error_stop(message)
      Else
        If (eqcm_data%voltage_range%value(1) > eqcm_data%voltage_range%value(2) ) Then
          If (fprint) Then
            Call info(' ',1)
          End If
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Flag "voltage_range": lowest potential MUST be set before&
                                                          & the highest potential'
          Call error_stop(message)
        Else If (Abs(eqcm_data%voltage_range%value(1)-eqcm_data%voltage_range%value(2)) &
              & < epsilon(eqcm_data%voltage_range%value(2)) ) Then
          If (fprint) Then
            Call info(' ',1)
          End If 
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Values for voltage_range potential are the same'
          Call error_stop(message)
        End If 
        ! convert units
        If (eqcm_data%voltage_range%units(1)=='v') Then
          eqcm_data%voltage_range%convert= 1.0_wp 
        ElseIf (eqcm_data%voltage_range%units(1)=='mv') Then
          eqcm_data%voltage_range%convert= 0.001_wp 
        Else 
         Write (messages(1),'(3(1x,a))')  Trim(error_set_eqcm), 'Invalid units for voltage range:', &
                                    & Trim(eqcm_data%voltage_range%units(1))
         Write (messages(2),'(1(1x,a))') 'Valid options: V or mV'
         Call info(messages, 2)
         Call error_stop(' ')
        End If
        
        eqcm_data%voltage_range%value=eqcm_data%voltage_range%convert*eqcm_data%voltage_range%value

        Write (message,'(1x,2(a,f6.3),a)') 'Voltage range:   ', eqcm_data%voltage_range%value(1), '  to  ', &
                                          & eqcm_data%voltage_range%value(2), ' [V]'
        If (fprint) Call info(message,1)
        If (eqcm_data%analysis%type == 'spectra') Then
           Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Specification of "voltage_range" is incompatible&
                                       & with "spectra" analysis. Please remove it and rerun.'
           Call error_stop(message)
        End If 
      End If
    End If    

    If (eqcm_data%V_to_Hz%fread) Then
      If (eqcm_data%V_to_Hz%fail) Then 
        If (fprint) Then
          Call info(' ',1)
        End If
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Missing (or wrong) settings for "V_to_Hz" directive'
        Call error_stop(message)
      Else
        If (Abs(eqcm_data%V_to_Hz%value) <= epsilon(eqcm_data%V_to_Hz%value)) Then
         If (fprint) Then 
           Call info(' ',1)
         End If
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Value for "V_to_Hz" must different from zero!'
          Call error_stop(message)
        End If
      End If
    End If 

    If (eqcm_data%current_offset%fread) Then
      If (eqcm_data%current_offset%fail) Then
        If (fprint) Then
          Call info(' ',1)
        End If
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Missing (or wrong) specification for directive&
                                  & "current_offset" (choose either .True. or .False.)'
        Call error_stop(message)
      End If
    Else
      eqcm_data%current_offset%warn=.True.
      eqcm_data%current_offset%stat=.True. 
    End If

           
    ! System directives 
    !!!!!!!!!!!!!!!!!!!
    If (electrode_data%area_geom%fread) Then
      If (electrode_data%area_geom%fail) Then 
        If (fprint) Then
          Call info(' ',1)
        End If
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Missing (or wrong) settings for "electrode_area" directive'
        Call error_stop(message)
      Else
        If (electrode_data%area_geom%value <= epsilon(electrode_data%area_geom%value)) Then
          If (fprint) Then
           Call info(' ',1)
          End If
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Value for "electrode_area" must be positive!'
          Call error_stop(message)
        End If
        
        ! Convert units
        If (electrode_data%area_geom%units  == 'inch2') Then
          electrode_data%area_geom%convert = 6.4516_wp
        Else If (electrode_data%area_geom%units  == 'cm2') Then
          electrode_data%area_geom%convert = 1.0_wp
        Else
          If (fprint) Then
            Call info(' ',1)
          End If
          Write (message,'(2(1x,a))')  Trim(error_set_eqcm), &
                                      &'Specified units for electrode area is not valid. Options: cm2 or inch2'
          Call error_stop(message)
        End If
        
        electrode_data%area_geom%value = electrode_data%area_geom%convert*electrode_data%area_geom%value
        Write (message,'(1x,a,f5.3,1x,a)') 'Electrode area:     ', electrode_data%area_geom%value, '[cm2]'
        If (fprint) Then
          Call info(message,1)
        End If
      End If
    End If 

    If (electrode_data%mass%fread) Then
      If (electrode_data%mass%fail) Then 
        If (fprint) Call info(' ', 1)
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Missing (or wrong) settings for "electrode_mass" directive'
        Call error_stop(message)
      Else
        If (electrode_data%mass%value <= epsilon(electrode_data%mass%value)) Then
         If (fprint) Call info(' ',1)
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Value for "electrode_mass" must be positive!'
          Call error_stop(message)
        End If
        
        ! Convert units
        If (electrode_data%mass%units  == 'g') Then
          electrode_data%mass%convert = 1.0E+9
        Else If (electrode_data%mass%units  == 'mg') Then
          electrode_data%mass%convert = 1.0E+6
        Else If (electrode_data%mass%units  == 'ug') Then
          electrode_data%mass%convert = 1.0E+3
        Else If (electrode_data%mass%units  == 'ng') Then
          electrode_data%mass%convert = 1.0_wp
        Else
          If (fprint) Call info(' ',1)
          Write (message,'(2(1x,a))')  Trim(error_set_eqcm), &
                                      &'Specified units for electrode mass is not valid. Options: g, mg, ug or ng'
          Call error_stop(message)
        End If
        
        electrode_data%mass%value = electrode_data%mass%convert*electrode_data%mass%value
        Write (message,'(1x,a,f12.3,1x,a)') 'Electrode mass:   ', electrode_data%mass%value, ' [ng]'
        If (fprint) Call info(message,1)
      End If
    End If 


    If (electrode_data%area_scale%fread) Then
      If (electrode_data%area_scale%fail) Then 
        If (fprint) Call info(' ',1)
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Missing (or wrong) settings for "area_scale" directive'
        Call error_stop(message)
      Else
        If (electrode_data%area_scale%value <= epsilon(electrode_data%area_scale%value)) Then
         If (fprint) Call info(' ',1)
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Value for "area_scale" must be larger than zero!'
          Call error_stop(message)
        End If
        
        electrode_data%area_geom%value = electrode_data%area_scale%value * electrode_data%area_geom%value
        Write (message,'(1x,a,f7.3,1x,a)') 'Effective area:   ', electrode_data%area_geom%value, '[cm2]'
        If (fprint) Call info(message,1)
      End If
    Else
      electrode_data%area_scale%value=1.0_wp
    End If 
            
    If (electrode_data%mole%fread) Then
      If (electrode_data%mole%fail) Then 
        If (fprint) Call info(' ',1)
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Missing (or wrong) settings for "electrode_mole" directive'
        Call error_stop(message)
      Else
        If (electrode_data%mole%value <= epsilon(electrode_data%mole%value)) Then
         If (fprint) Call info(' ',1)
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Value for "electrode_mole" must be positive!'
          Call error_stop(message)
        End If
        Write (message,'(1x,a,E15.7)') 'Number of moles:  ', electrode_data%mole%value
        If (fprint) Call info(message,1)
      End If
    End If 

    If (system_data%quartz_freq%fread) Then
      If (system_data%quartz_freq%fail) Then 
        If (fprint) Call info(' ',1)
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Missing (or wrong) settings for "quartz_frequency" directive'
        Call error_stop(message)
      Else
       If (system_data%quartz_freq%value <= epsilon(system_data%quartz_freq%value)) Then
         If (fprint) Call info(' ',1)
         Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Value for "quartz_frequency" must be larger than zero!'
         Call error_stop(message)
       End If
       !Convert units 
       If (system_data%quartz_freq%units  == 'mhz') Then
         system_data%quartz_freq%convert = 1000000.0_wp
       Else If (system_data%quartz_freq%units  == 'khz') Then
         system_data%quartz_freq%convert = 1000.0_wp
       Else If (system_data%quartz_freq%units  == 'hz') Then
         system_data%quartz_freq%convert = 1.0_wp
       Else
         If (fprint) Call info(' ',1)
         Write (message,'(2(1x,a))')  Trim(error_set_eqcm), &
                                     &'Specified units for quartz frequency is not valid. Options: MHz, kHz or Hz'
         Call error_stop(message)
       End If
       system_data%quartz_freq%value=system_data%quartz_freq%convert * system_data%quartz_freq%value
       Write (message,'(1x,a,f7.3,1x,a)') 'Quartz frequency: ', system_data%quartz_freq%value/1000000.0_wp, '[MHz]'
       If (fprint) Call info(message,1)
      End If
    End If 


    ! Frequency to mass conversion directives
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    If (sauer_data%factor%fread) Then
      If (sauer_data%factor%fail) Then 
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong (or missing) settings for sauerbrey directive.'
        Call error_stop(message)
      Else 
        If (Abs(sauer_data%factor%value(1)) < epsilon(sauer_data%factor%value(1))) Then
          If (fprint) Call info(' ',1)
            Write (message,'(2(1x,a))') Trim(error_set_eqcm), &
                                    &'Input value for Sauerbrey factor MUST be different from zero'
            Call error_stop(message)
        End If

        If (sauer_data%factor%units(1) /= 'hz'   .Or. &
            sauer_data%factor%units(2) /= 'ng-1' .Or. &
            sauer_data%factor%units(3) /= 'cm2'   ) Then
           Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Units for sauberbrey directive MUST be [Hz ng-1 cm2]'
           Call error_stop(message)
        End If
        Write (message,'(1x,a,f10.5,1x,a)') 'Sauerbrey factor: ', sauer_data%factor%value(1), '[Hz ng-1 cm2]'
        If (fprint) Call info(message,1)
      End If
    End If 

    ! Filtering directives 
    !!!!!!!!!!!!!!!!!!!!!!! 
    If (filter_data%cutoff%fread) Then
      If (filter_data%cutoff%fail) Then   
        If (fprint) Call info(' ',1)
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong specification for directive "filter_cutoff" '
        Call error_stop(message)
      Else 
        If (filter_data%cutoff%value <= epsilon(filter_data%cutoff%value)) Then
          If (fprint) Call info(' ',1)
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Value of "filter_cutoff" should be larger than zero'
          Call error_stop(message)
        End If
        If (filter_data%cutoff%units /= 'v-1' .And. filter_data%cutoff%units /= 'hz') Then
          If (fprint) Call info(' ',1)
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong units for "filter_cutoff" '
          Call error_stop(message)
        End If
        If (fprint) Call info(' ',1)
        If (filter_data%cutoff%units == 'v-1')then
          Write (message,'(1x,a,f12.3,1x,a)') 'Low-pass filter cutoff:', filter_data%cutoff%value, &
                                           & '[1/V]'
        Else
          Write (message,'(1x,a,f12.3,1x,a)') 'Low-pass filter cutoff:', filter_data%cutoff%value, &
                                           & '[Hz]'
        End If
        If (fprint) Call info(message,1)
      End If
    Else 
      If (eqcm_data%analysis%type/='spectra' .And. & 
          eqcm_data%analysis%type/='print_eqcm_raw') Then 
        Write (message,'(1x,a)') 'EQCM analysis with raw data (no filtering)'
        If (fprint) Then
          Call info(' ', 1)
          Call info(message,1)
        End If
      End If 
    End If 

    ! Required strings in case of error for endpoints (avoids duplication)
    Write (messages(2),'(1x,a)') 'of the set of values. The average of these points is used to set the value of'
    Write (messages(3),'(1x,a)') 'the padding points (outside the range of meassured values) prior to the FFT.'
    Write (messages(4),'(1x,a)') 'This number MUST be always greater than 0. Default value is 1. '

    If (fft_data%end_current%fread) Then
      Write (messages(1),'(1x,a)') '***Problems: "endpoints_current" is the number of current points taken from each extreme'
      If (fft_data%end_current%fail) Then
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong (or missing) specification for directive "endpoints_current"'
        If (fprint) Call info(' ',1)
        If (fprint) Call info(messages,4)
        If (fprint) Call info(' ',1)
        Call error_stop(message)
      Else
        If (fft_data%end_current%value < 1) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), '"endpoints_current" should be larger than 0'
          If (fprint) Call info(' ',1)
          If (fprint) Call info(messages,4)
          If (fprint) Call info(' ',1)
          Call error_stop(message)
        End If
      End If
    Else
      fft_data%end_current%value=1  
    End If
 
    If (fft_data%end_mass_frequency%fread) Then
      Write (messages(1),'(1x,a)') '***Problems: "endpoints_mass_frequency" is the number of mass-frequency '&
                                  'points taken from each extreme'
      If (fft_data%end_mass_frequency%fail) Then
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong (or missing) specification for directive &
                                   &"endpoints_mass_frequency"'
        If (fprint) Call info(' ',1)
        If (fprint) Call info(messages,4)
        If (fprint) Call info(' ',1)
        Call error_stop(message)
      Else
        If (fft_data%end_mass_frequency%value < 1) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), '"endpoints_mass_frequency" should be larger than 0'
          If (fprint) Call info(' ',1)
          If (fprint) Call info(messages,4)
          If (fprint) Call info(' ',1)
          Call error_stop(message)
        End If
      End If 
    Else 
      fft_data%end_mass_frequency%value=1
    End If

  ! Stoichiometric directives
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
    If (eqcm_data%analysis%type == 'stoichiometry' .Or. & 
       eqcm_data%analysis%type == 'model_pristine_sample' .Or. &
       eqcm_data%analysis%type == 'model_disordered_system' .Or. &
       eqcm_data%analysis%type == 'model_cycled_sample') Then
       Call check_stoich_settings(files, eqcm_data, stoich_data)
    End If
    
    If (eqcm_data%analysis%type == 'stoichiometry' .Or. & 
       eqcm_data%analysis%type == 'model_cycled_sample') Then
        If (eqcm_data%process%type == 'intercalation') Then
          Write (messages(1),'(1x,a)')  'ALC_EQCM computes the set of stoichiometric solutions compatible with EQCM data.'
          Write (messages(2),'(1x,a)')  'If the solution is unequivocally defined, the algorithm provides stoichiometric'
          Write (messages(3),'(1x,a)')  'resolution as a function of the EQCM cycling.'
          Call info(messages,3)
        Else If (eqcm_data%process%type == 'electrodeposition') Then
          Write (messages(1),'(1x,a)')  'The implemented algorithm only allows a single species to be deposited/removed'
          Write (messages(2),'(1x,a)')  'onto/from the electrode surface. Results give the ratio between the effective'
          Write (messages(3),'(1x,a)')  'geometrical areas, as well the amount of moles of the species that participate'
          Write (messages(4),'(1x,a)')  'in the reaction. Species can be composed of multiple subspecies.'
          Write (messages(5),'(1x,a)')  'Each subspecies must contribute with the same number of moles.'
          Call info(messages,5)
        End If
    End If

    If (eqcm_data%analysis%type == 'stoichiometry' .Or. & 
       eqcm_data%analysis%type == 'model_cycled_sample') Then

      If (eqcm_data%efficiency%fread) Then
        If (eqcm_data%efficiency%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Missing (or wrong) settings for "efficiency" directive'
          Call error_stop(message)
        Else
          If (eqcm_data%efficiency%value <= epsilon(eqcm_data%efficiency%value) .Or.  &
             eqcm_data%efficiency%value  > 1.0_wp ) Then
             Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Value for "efficiency" must be positive and <= 1.0'
             Call error_stop(message)
          End If
        End If
      Else
        eqcm_data%efficiency%value = 1.0_wp
      End If

      If (stoich_data%discretization%fread) Then
        If (stoich_data%discretization%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Missing (or wrong) settings for "delta_stoich" directive'
          Call error_stop(message)
        Else
          If (stoich_data%discretization%value <= epsilon(stoich_data%discretization%value)) Then
             Write (messages,'(2(1x,a))') Trim(error_set_eqcm), 'Value for "sotich_discretization" must be positive'
             Call error_stop(message)
          End If
        End If
      Else
        stoich_data%discretization%value = ds_default
      End If

       If (electrode_data%mole%fread .And. electrode_data%mass%fread) Then
         If (fprint) Call info(' ',1)
         Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Both "electrode_mole" and "electrode_mass" have been specified. &
                                   &This information is redundant and could lead to mistakes.'
         If (fprint) Call info(message,1)
         Write (message,'(1x,a)') 'Review the settings and select only one of these directives'
         If (fprint) Call info(message,1)
         Call error_stop(' ')
       End If

       If (eqcm_data%range_cycles%fread .And. eqcm_data%range_cycles%value(1) /= 1) Then
         If (fprint) Call info(' ',1)
         Write (message,'(1x,4a)') Trim(error_set_eqcm), ' The low limit of directive "cycle" for ', &
                                    Trim(eqcm_data%analysis%type), ' analysis MUST BE set to 1.'
         If (fprint) Call info(message,1)
         Write (message,'(2(1x,a))') 'Only the initial stoichiometric composition can be unequivocally determined.',& 
                                    'Please either change the specification for "cycle" or remove the directive&
                                    & (to consider all cycles)'
         If (fprint) Call info(message,1)
         Call error_stop(' ')
        
       End If

       If (stoich_data%block_constraints%fread) Then
         If (Trim(eqcm_data%process%type)=='intercalation') Then
           Call check_constraints_settings(files, stoich_data) 
         If (stoich_data%num_variables < 3) Then
           stoich_data%block_constraints%warn=.True.
         End If
         Else
           Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Definition of &block_constraints_species&
                                              & is incompatible with electrodeposition analysis. See manual'
           Call error_stop(message)
         End If
       End If

    End If

  ! Settings for atomistic modelling
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    If (Trim(eqcm_data%analysis%type) == 'model_pristine_sample' .Or. &
       Trim(eqcm_data%analysis%type) == 'model_cycled_sample'   .Or. &
       Trim(eqcm_data%analysis%type) == 'model_disordered_system' .Or. &
       Trim(eqcm_data%analysis%type) == 'hpc_simulation_files' ) Then

      ! rotate_species
      If (model_data%rotate_species%fread) Then
        If (model_data%rotate_species%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong specification for directive "rotate_species". It must be &
                                    & either .True. or .False. (.True. by default)'
          Call error_stop(message)
        End If
      Else
        model_data%rotate_species%stat=.True.
      End If

      ! repeat_input_model
      If (model_data%repeat_input_model%fread) Then
        If (model_data%repeat_input_model%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm),&
                                  & 'Wrong (or missing) settings for "repeat_input_model" directive (see manual).'
          Call error_stop(message)
        Else
          Do j = 1, 3
            If (model_data%repeat_input_model%value(j) < 1) Then
              Write (message,'(1x,2a, i2,a)') Trim(error_set_eqcm), &
                                      &' Input value', j, ' for directive "repeat_input_model" must be positive'
              Call error_stop(message)
            End If
          End Do 
        End If
      Else  
        model_data%repeat_input_model%value=1
      End If

      ! targeted_num_models
      If (model_data%targeted_num_models%fread) Then
        If (model_data%targeted_num_models%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm),&
                                  & 'Wrong (or missing) settings for "targeted_num_models" directive.'
          Call error_stop(message)
        Else
          If (model_data%targeted_num_models%value <= 1) Then
            Write (message,'(1x,2a)') Trim(error_set_eqcm), ' Directive "targeted_num_models" must be larger than 1.'
            Call error_stop(message)
          End If
        End If
      End If

      ! multiple_input_atoms 
      If (model_data%multiple_input_atoms%fread) Then
        If (model_data%multiple_input_atoms%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm),&
                                  & 'Wrong (or missing) settings for "multiple_input_atoms" directive.&
                                  & Is it not too large? It should not exceed 1 millon.'
          Call error_stop(message)
        Else
          If (model_data%multiple_input_atoms%value <= 0) Then
            Write (message,'(1x,2a)') Trim(error_set_eqcm), ' Directive "multiple_input_atoms" must be larger than 0.'
            Call error_stop(message)
          End If
        End If
      Else
        model_data%multiple_input_atoms%value=10000
      End If

      ! block species components
      If (stoich_data%block_species_components%fread) Then
        Call check_components_species(files, stoich_data) 
      Else 
        Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Missing &block_species_components, which is required to specify the&
                                                        & atomistic details for the components of each chemical species.' 
        Call error_stop(message)
      End If
      Call check_atomic_settings(files, eqcm_data, model_data)

      ! delta_space
      If (model_data%delta_space%fread) Then
        If (model_data%delta_space%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong settings for "delta_space" directive.'
          Call error_stop(message)
        Else
          If (model_data%delta_space%value < epsilon(model_data%delta_space%value)) Then
            If (fprint) Call info(' ',1)
            Write (message,'(2(1x,a))') Trim(error_set_eqcm), &
                                      &'Input value for delta_space MUST be larger than zero'
            Call error_stop(message)
          End If
          If (Trim(model_data%delta_space%units) == 'angstrom' ) Then
            model_data%delta_space%convert=1.0_wp
          Else If (Trim(model_data%delta_space%units) == 'bohr' ) Then
            model_data%delta_space%convert= Bohr_to_A
          Else
            Write (message,'(2a)')  Trim(error_set_eqcm), 'Invalid units for directive "delta_space". Options: Angstrom or Bohr'
            Call info(message, 1)
            Call error_stop(' ') 
          End If
          model_data%delta_space%units='Angstrom' 
          model_data%delta_space%value=model_data%delta_space%value*&
                                                     model_data%delta_space%convert
        End If
      Else
        model_data%delta_space%units='Angstrom'
        model_data%delta_space%value= 0.1_wp 
      End If

      ! distance cutoff
      If (model_data%distance_cutoff%fread) Then
        If (model_data%distance_cutoff%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong settings for "distance_cutoff" directive.'
          Call error_stop(message)
        Else
          If (model_data%distance_cutoff%value < epsilon(model_data%distance_cutoff%value)) Then
            If (fprint) Call info(' ',1)
            Write (message,'(2(1x,a))') Trim(error_set_eqcm), &
                                      &'Input value for "distance_cutoff" MUST be larger than zero'
            Call error_stop(message)
          End If
          If (Trim(model_data%distance_cutoff%units) == 'angstrom' ) Then
            model_data%distance_cutoff%convert=1.0_wp
          Else If (Trim(model_data%distance_cutoff%units) == 'bohr' ) Then
            model_data%distance_cutoff%convert= Bohr_to_A
          Else
            Write (message,'(2a)')  Trim(error_set_eqcm), 'Invalid units for directive "distance_cutoff".&
                                  & Options: Angstrom or Bohr'
            Call info(message, 1)
            Call error_stop(' ') 
          End If
          model_data%distance_cutoff%units='Angstrom' 
          model_data%distance_cutoff%value=model_data%distance_cutoff%value*&
                                                     model_data%distance_cutoff%convert
          If (model_data%distance_cutoff%value < min_inter) Then
            Write (message,'(1x,a,2(a,f4.2),a)') Trim(error_set_eqcm), &
                                            & ' The minimum limit for the inter-species distance of ',&
                                            & model_data%distance_cutoff%value,&
                                            &' Angstrom (set via "distance_cutoff" directive) is unphysical.&
                                            & The value MUST NOT be lower than ', min_inter, ' Angstrom. Please change'
            Call info(message,1)
            Call error_stop(' ')
          End If
        End If
      Else
        model_data%distance_cutoff%units='Angstrom'
        model_data%distance_cutoff%value= min_inter 
      End If

      ! Deposition level 
      If (model_data%deposition_level%fread) Then
        If (eqcm_data%process%type /= 'electrodeposition') Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Definition of "deposition_level" directive is only&
                                     & compatible with option "electrodeposition" for directive "process".'
          Call error_stop(message)
        End If
        If (model_data%deposition_level%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm), 'Wrong settings for "deposition_level" directive.'
          Call error_stop(message)
        Else
          If (Trim(model_data%deposition_level%units) == 'angstrom' ) Then
            model_data%deposition_level%convert=1.0_wp
          Else If (Trim(model_data%deposition_level%units) == 'bohr' ) Then
            model_data%deposition_level%convert= Bohr_to_A
          Else
            Write (message,'(2a)')  Trim(error_set_eqcm), 'Invalid units for directive "deposition_level".&
                                  & Options: Angstrom or Bohr'
            Call info(message, 1)
            Call error_stop(' ') 
          End If
          model_data%deposition_level%units='Angstrom' 
          model_data%deposition_level%value=model_data%deposition_level%value*model_data%deposition_level%convert
        End If
      End If

      ! stoichiometry error
      If (model_data%stoichiometry_error%fread) Then
        If (model_data%stoichiometry_error%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm),&
                                  & 'Wrong (or missing) settings for "stoichiometry_error" directive (see manual).'
          Call error_stop(message)
        Else
          If (model_data%stoichiometry_error%value<0.0_wp) Then
            Write (message,'(2(1x,a))') Trim(error_set_eqcm),&
                                  & 'Value for directive "stoichiometry_error" must be positive'
            Call error_stop(message)
          Else If (model_data%stoichiometry_error%value>0.1_wp) Then
            Write (message,'(2(1x,a))') Trim(error_set_eqcm),&
                                  & 'Value for directive "stoichiometry_error" MUST NOT be larger than 0.1, as this&
                                  & limit already represent a rather large error margin'
            Call error_stop(message)
          End If
        End If
      Else
        model_data%stoichiometry_error%value= 0.01_wp  
      End If

      ! scale cell 
      If (model_data%scale_cell%fread) Then
        If (model_data%scale_cell%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set_eqcm),&
                                  & 'Wrong (or missing) settings for "scale_cell" directive (see manual).'
          Call error_stop(message)
        Else
          If (model_data%scale_cell%value<0.0_wp) Then
            Write (message,'(2(1x,a))') Trim(error_set_eqcm),&
                                  & 'Value for directive "scale_cell" must be positive'
            Call error_stop(message)
          End If
        End If
      Else
        model_data%scale_cell%value= 1.0_wp  
      End If

    End If

  ! Settings for building simulation files 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    If (Trim(eqcm_data%analysis%type) /= 'hpc_simulation_files') Then
      If (Trim(eqcm_data%analysis%type) /= 'model_pristine_sample'   .And. &
         Trim(eqcm_data%analysis%type) /= 'model_disordered_system' .And. &
         Trim(eqcm_data%analysis%type) /= 'model_cycled_sample') Then
        If (simulation_data%generate .Or. hpc_data%generate) Then
           Write (messages(1),'(2(1x,a))') Trim(error_set_eqcm), &
                                      & 'The user has requested the generation of files for simulations (dft_settings)&
                                      & and/or the generation of HPC files (block_HPC_settings).'
           Write (messages(2),'(1x,a)') 'Any of these features, however, requires that directive "Analysis" is set to either:' 
           Write (messages(3),'(1x,a)')  ' - model_pristine_sample'
           Write (messages(4),'(1x,a)')  ' - model_cycled_sample'
           Write (messages(5),'(1x,a)')  ' - model_disordered_system'
           Write (messages(6),'(1x,a)')  ' - hpc_simulation_files'
           Write (messages(7),'(1x,a)')  'Please correct the settings.'
           Call info(messages, 7)
           Call error_stop(' ')
        End If
      Else If (Trim(eqcm_data%analysis%type) == 'model_pristine_sample' .Or. &
              Trim(eqcm_data%analysis%type) /= 'model_disordered_system' .Or. &
              Trim(eqcm_data%analysis%type) == 'model_cycled_sample') Then  
        If (simulation_data%generate) Then
          simulation_data%code_format=model_data%output_model_format%type
          simulation_data%process=eqcm_data%process%type
          simulation_data%normal_vector=model_data%normal_vector%type
          Call check_simulation_settings(files, stoich_data, simulation_data)
        End If
        If (hpc_data%generate) Then
          Call check_hpc_settings(files, model_data%output_model_format%type, hpc_data)
        End If
      End If

    Else
      If ((.Not. simulation_data%generate) .And. (.Not. hpc_data%generate)) Then  
        Write (messages(1),'(2(1x,a))') Trim(error_set_eqcm), &
                                   & 'The user has requested the generation of files for simulations (dft_settings)&
                                   & and/or the generation of HPC files (block_HPC_settings).'
        Write (messages(2),'(1x,a)') 'However, none of these blocks is defined. Requested analysis is not possible.&
                                   & Please correct.'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
      If (simulation_data%generate) Then
        simulation_data%code_format=model_data%output_model_format%type
        simulation_data%process=eqcm_data%process%type
        simulation_data%normal_vector=model_data%normal_vector%type
        Call check_simulation_settings(files, stoich_data, simulation_data)
      End If
      If (hpc_data%generate) Then
        Call check_hpc_settings(files, model_data%output_model_format%type, hpc_data)
      End If
    End If

    ! Write warnings
    !!!!!!!!!!!!!!!!
    If (eqcm_data%current_offset%warn) Then
       If (eqcm_data%analysis%type /= 'print_eqcm_raw'    .And. &
           eqcm_data%analysis%type /= 'print_eqcm_filter' .And. &
           eqcm_data%analysis%type /= 'spectra'  )Then
         If (fprint) Call info(' ',1)
         Write (messages(1),'(1x,a)') '***IMPORTANT: By default, the offset for EQCM current is activated. If the user wants to'
         Write (messages(2),'(1x,a)') '   deactivate the offset, directive "current_offset" must be set to .False.'
         If (fprint) Call info(messages,2)
       End If
    End If

    If (eqcm_data%software%warn) Then
      If (fprint) Call info(' ',1)
      Write (message,'(1x,a)') '***WARNING: Check specification for directive "software".' 
      If (fprint) Call info(message,1)
    end If

    If (eqcm_data%range_cycles%warn) Then
      If (fprint) Call info(' ',1)
      Write (messages(1),'(1x,a)') '***IMPORTANT: Directive "cycles" (cycle range) has not been specified.'
      Write (messages(2),'(1x,a)') '   Thus, all CV cycles will be considered in the analysis.'
      If (fprint) Call info(messages, 2)
    End If

    If (stoich_data%block_constraints%warn) Then
      Write (messages(1),'(2(1x,a))') '***WARNING: The number of stoichiometric variables is <= 2. Information for constraints,&
                                    & as specified by &block_constraints_species, is irrelevant and it will be ignored.'
      If (fprint) Then
        Call info(messages,1)
      End If 
    End If

    ! refresh out_eqcm
    Call refresh_out_eqcm(files)

  End Subroutine check_all_settings


  Subroutine check_feasibility(files, eqcm_data, filter_data, system_data, electrode_data, sauer_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the feasibility of the requested analysis.
    ! Criteria is based on the settings read in SET_EQCM as well as the variables
    ! and the number of cycles identified from DATA_EQCM file, with the exception of 
    ! model_pristine_sample and model_disordered_system type of analysis, which do 
    ! not require EQCM data.
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut) :: files(:)
    Type(eqcm_type),      Intent(InOut) :: eqcm_data
    Type(filter_type),    Intent(InOut) :: filter_data
    Type(system_type),    Intent(InOut) :: system_data
    Type(electrode_type), Intent(InOut) :: electrode_data
    Type(sauer_type),     Intent(InOut) :: sauer_data

    Character(Len=256)                :: message
    Character(Len=256)                :: messages(6)
    Character(Len= 32)                :: error_set
    Character(Len= 32)                :: set_eqcm
    Character(Len= 32)                :: data_eqcm

    error_set = '***ERROR -'
    set_eqcm   = Trim(files(FILE_SET_EQCM)%filename) 
    data_eqcm  = Trim(files(FILE_DATA_EQCM)%filename)

    If (.Not. eqcm_data%voltage%fread) Then
      Write (message,'(4(1x,a))') Trim(error_set), 'No values for applied potential/voltage are found in file',&
                                &Trim(data_eqcm), '. EQCM analysis is not possible'
      Call error_stop(message)
    End If

    If (eqcm_data%range_cycles%fread) Then 
      If (eqcm_data%ncycles < eqcm_data%range_cycles%value(2)) Then
        Call info(' ',1)
        Write (message,'(3(1x,a),2(a,i3),2(1x,a))') Trim(error_set), 'Number of cycles identified in', Trim(data_eqcm),&
                                                 &' is ', eqcm_data%ncycles, ', which is lower than the upper bound of ', &
                                                 & eqcm_data%range_cycles%value(2), ' specified by directive "cycles" in file', &
                                                 & Trim(set_eqcm)  
        Call error_stop(message)
      End If
    Else
      ! if directive cycle is missing, use all the CV cycles found in DATA_EQCM
      eqcm_data%range_cycles%value(1)=1
      eqcm_data%range_cycles%value(2)=eqcm_data%ncycles 
    End If

    If (eqcm_data%analysis%type == 'mass_calibration' .Or. eqcm_data%analysis%type == 'massogram') Then
      If ((.Not. eqcm_data%mass_frequency%fread) .Or. (.Not. eqcm_data%current%fread)) Then 
         Write (message,'(8(1x,a))') Trim(error_set), 'No values for mass-frequency and/or current found in', &
                                 & Trim(data_eqcm),'file. Requested', Trim(eqcm_data%analysis%type),'analysis in file',  & 
                                 & Trim(set_eqcm),'is not possible.'
         Call error_stop(message)
      End If
    End If

    If (eqcm_data%analysis%type == 'characterization') Then
      If ((.Not. eqcm_data%current%fread)) Then
         Write (message,'(7(1x,a))') Trim(error_set), 'No values for current found in', Trim(data_eqcm),&
                                 &'file. Requested', Trim(eqcm_data%analysis%type),'analysis is not possible.'
         Call error_stop(message)
      End If    
    End If

    If (.Not. eqcm_data%mass_frequency%fread) Then
      If (eqcm_data%analysis%type == 'stoichiometry'    .Or.  &
         eqcm_data%analysis%type == 'model_cycled_sample' ) Then
         Write (message,'(8(1x,a))') Trim(error_set), 'No values for mass-frequency found in', Trim(data_eqcm),&
                                 &'file. Requested', Trim(eqcm_data%analysis%type),'analysis in file', Trim(set_eqcm),&
                                 &'is not possible.'
        Call error_stop(message)
      End If
    End If

    If ((.Not. eqcm_data%current%fread)) Then
      If (eqcm_data%analysis%type == 'stoichiometry'    .Or.  &
         eqcm_data%analysis%type == 'model_cycled_sample' ) Then
         Write (message,'(8(1x,a))') Trim(error_set), 'No values for current found in', Trim(data_eqcm),&
                                 &'file. Requested', Trim(eqcm_data%analysis%type),'analysis in file', Trim(set_eqcm),&
                                 &'is not possible.'
        Call error_stop(message)
      End If
    End If

    If (eqcm_data%mass_frequency%fread) Then
      If (sauer_data%factor%fread)then  
        Write (messages(1),'(1x,a)')    '****PROBLEMS'
        Write (messages(4),'(3(1x,a))') 'using the Sauerbrey factor defined in', Trim(set_eqcm), 'file.'

        If (system_data%quartz_freq%value <= epsilon(system_data%quartz_freq%value)) Then
          Write (messages(2),'(3(1x,a))') 'The', Trim(eqcm_data%analysis%type), ' option chosen for EQCM analysis requires the'
          Write (messages(3),'(3(1x,a))') 'convertion of mass-frequency values found in file', Trim(data_eqcm),&
                                        &'to mass-density mass'
          Write (messages(5),'(1x,1a)') 'To corroborate the validation of the Sauerbrey approximation,&
                                       & the change in mass-frequency'
          Write (messages(6),'(1x,1a)') 'must be compared with the quartz frequency of the EQCM device.'
          Call info(messages,6)
          Call info(' ',1)
          Write (message,'(3(1x,a))') Trim(error_set), 'No specification for "quartz_freq" in file', Trim(set_eqcm)
          Call error_stop(message)
        End If
  
        If (eqcm_data%analysis%type == 'characterization' .Or.  &
           eqcm_data%analysis%type == 'stoichiometry'    .Or.  &
           eqcm_data%analysis%type == 'model_cycled_sample' ) Then
          If (electrode_data%area_geom%value <= epsilon(electrode_data%area_geom%value)) Then
            Write (messages(2),'(3(1x,a))') 'The', Trim(eqcm_data%analysis%type), &
                                          &'option chosen for EQCM analysis requires conversion'
            Write (messages(3),'(3(1x,a))') 'of mass-frequency values found in file', Trim(data_eqcm), 'to total mass'
            Write (messages(5),'(1x,1a)')   'However, this requires information of the electrode area, under the implicit'
            Write (messages(6),'(1x,1a)')   'assumption the reaction occurs evenly throghout the sample.'
            Call info(messages,6)
            Call info(' ',1)
            Write (message,'(3(1x,a))') Trim(error_set), 'No specification for "electrode_area" in file', Trim(set_eqcm)
            Call error_stop(message)
          End If
        End If

      End If
    End If
  
    ! Stop the execution if filter is the user wants to print raw data but filter_cutoff is turned on
    If (eqcm_data%analysis%type == 'print_eqcm_raw') Then 
      If (filter_data%cutoff%fread) Then 
        Write (message,'(4(1x,a))') Trim(error_set), 'Analysis', Trim(eqcm_data%analysis%type), 'is not possible if&
                                 & "Filter_cutoff" is set. Please remove directive "Filter_cutoff" and run again'
        Call error_stop(message)
      End If 
    End If


    ! Check conditions for spectra and print_eqcm_raw (or print_eqcm_filter)
    If (eqcm_data%analysis%type == 'spectra'  .Or. eqcm_data%analysis%type == 'print_eqcm_raw' .Or. &
                                             eqcm_data%analysis%type == 'print_eqcm_filter') Then
      If ((.Not. eqcm_data%current%fread) .And. (.Not. eqcm_data%mass_frequency%fread) ) Then
         Write (message,'(6(1x,a))') Trim(error_set), 'No values for current nor mass-frequency found in', Trim(data_eqcm),&
                                    &'file. Requested', Trim(eqcm_data%analysis%type), 'analysis is not possible.'
         Call error_stop(message)
      End If
    End If

    ! Check conditions for printing and filtering
    If (eqcm_data%analysis%type == 'print_eqcm_filter'  .And. (.Not. filter_data%cutoff%fread) ) Then
      Call info(' ',1)
      Write (message,'(3(1x,a))') Trim(error_set), 'No specification for "filter_cutoff" in file', Trim(set_eqcm)
      Call error_stop(message)
    End If
   
    If (eqcm_data%analysis%type == 'massogram') Then
      If (.Not.  eqcm_data%time%fread) Then
        If (.Not.  eqcm_data%scan_rate%fread) Then
        Write (messages(1),'(3(1x,a))') 'File', Trim(data_eqcm), 'contains no elapsed time.'
        Write (messages(2),'(1x,a)')    'To compute the massogram (dM/dt) from the applied voltage, it is required to'
        Write (messages(3),'(3(1x,a))') 'specify the scan_rate in file', Trim(set_eqcm), ' used in the experiment.'
        Call info(messages,3)
        Call info(' ',1)
        Write (message,'(3(1x,a))') Trim(error_set), 'Wrong (or No) value of scan_rate in the file', Trim(set_eqcm)
        Call error_stop(message)
        End If
      End If
    End If

    If (eqcm_data%analysis%type == 'characterization' .Or.  eqcm_data%analysis%type == 'mass_calibration') Then
      If (.Not. eqcm_data%time%fread .And. eqcm_data%current%fread) Then
        If (.Not.  eqcm_data%scan_rate%fread) Then
        Write (messages(1),'(3(1x,a))') 'File', Trim(data_eqcm), 'contains the meassured current but no elapsed time.' 
        Write (messages(2),'(1x,a)')    'Thus, to integrate the current and compute the total charge involved'
        Write (messages(3),'(2(1x,a))') 'in the process, the sc_rate directive MUST be specified in file', Trim(set_eqcm)
        Call info(messages,3)
        Call info(' ',1)
        Write (message,'(3(1x,a))') Trim(error_set), 'Wrong (or No) value of scan_rate in the file', Trim(set_eqcm)
        Call error_stop(message)
        End If
      End If
    End If

   !!!!!!!!!!!!!!!!!!!!!!!
   ! Redefine reciprocal space units
   !!!!!!!!!!!!!!!!!!!!!!!

   If ( (.Not. eqcm_data%time%fread) .And. filter_data%cutoff%fread ) Then
     If (filter_data%cutoff%units == 'hz') Then
       Write (messages(1),'(1x,a)')    '***WARNING'
       Write (messages(2),'(3(1x,a))') 'Cutoff frequency for filtering is set in units of Hz in', Trim(set_eqcm), 'file &
                                     &(see above).'  
       Write (messages(3),'(3(1x,a))') 'However, there is no elapsed time recorded in', Trim(data_eqcm), 'file. The value'
       Write (messages(4),'(1x,a)') 'for the cutoff in the 1/V domain is obtained by filter_cutoff/scan_rate, assuming the'
       Write (messages(5),'(1x,a)') 'scan_rate was kept constant during the CV (standard protocol). If the scan rate was changed'
       Write (messages(6),'(1x,a)') 'during the experiment, the cutoff MUST be given in units of 1/V (need of spectra analysis).'
       Call info(messages,6)
       If (eqcm_data%scan_rate%value(1) < epsilon(eqcm_data%scan_rate%value(1))) Then
         Call info(' ',1)
         Write (message,'(3(1x,a))') Trim(error_set), 'Wrong (or No) value of scan_rate in the file', Trim(set_eqcm)
         Call error_stop(message)
       End If 
       filter_data%cutoff%value= filter_data%cutoff%value/eqcm_data%scan_rate%value(1)
       Write (message,'(1x,a,f12.5,1x,a)') 'New cutoff: ', filter_data%cutoff%value, '[1/V]'  
       Call info(message,1)
     End If
   End If

    ! refresh out_eqcm
    Call refresh_out_eqcm(files)

  End Subroutine check_feasibility

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutines to read parameters and settings for atomistic models
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine check_end(io, string)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check if there is missing data and the end of the file
    ! has been reached
    !
    ! author    - i.scivetti Dec 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer,          Intent(In   ) :: io
    Character(Len=*), Intent(In   ) :: string

    Character(Len=256) :: messages(2)

    If (is_iostat_end(io))Then
      Call info(' ', 1)
      Write (messages(1),'(1x,2a)') '*** ERROR in ', Trim(string)
      Write (messages(2),'(1x,2a)') 'End of file is detected. It seems there is missing data or the block is not&
                                  & closed properly. Please check'
      Call info(messages, 2)
      Call error_stop(' ')
    End If

  End Subroutine check_end

  Subroutine read_block_input_cell(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the simulation cell vectors of the input structure
    ! defined in &block_input_cell
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Type(model_type),   Intent(InOut) :: model_data
   
    Integer(Kind=wi)   :: io, i, j
    Character(Len=64 ) :: error_block_input_cell
    Character(Len=256) :: messages(2), word
    Logical            :: endblock

    error_block_input_cell = '***ERROR in &block_input_cell of SET_EQCM file'
    Write (messages(1),'(a)') error_block_input_cell

    i=1
   
    Do While (i <= 3)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&block_input_cell')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&block_input_cell') 
          Read (iunit, Fmt=*, iostat=io) (model_data%input%cell(i,j), j=1,3)
          If (io/=0) Then
            Write (messages(2),'(a,i2)') 'Problems with the definition of cell vector', i
            Call info(messages, 2)
            Call error_stop(' ')
          End If
          i=i+1
        Else
          Write (messages(2),'(1x,a)') 'End of block found! Not all the cell vectors for the&  
                                     & input structure have been defined. Please check.'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End If
    End Do
 
    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&block_input_cell')
      Call capital_to_lower_case(word)
      If (word /= '&end_block_input_cell') Then
        If (word(1:1) /= '#') Then
          Write (messages(2),'(a)') 'Block for cell vectors must be closed with&
                                  & sentence &end_block_input_cell.'
          Call info(messages,2)
          Call error_stop(' ')
        End If
      Else
          endblock=.True.
      End If
    End Do

  End Subroutine read_block_input_cell

  Subroutine read_block_input_composition(iunit, stoich_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the amount of atoms of each atomic components of the 
    ! participating species from &block_input_composition
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(stoich_type), Intent(InOut) :: stoich_data
    Type(model_type),  Intent(InOut) :: model_data

    Integer(Kind=wi)   :: io, i, j, k
    Character(Len=256) :: messages(8), word
    Character(Len=64 ) :: error_block_input_composition
    Logical  :: endblock, loop
    Logical  :: header(2), error_duplication

    error_block_input_composition = '***ERROR in &block_input_composition of SET_EQCM file'
    Write (messages(1),'(a)') error_block_input_composition

    header=.False.
    error_duplication=.False.

    i=1

    Do While (i <= 2)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&block_input_composition')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call capital_to_lower_case(word) 
          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Call check_for_rubbish(iunit, '&block_input_composition') 
              Read (iunit, Fmt=*, iostat=io) word, (model_data%block_component%tag(j), j = 1, stoich_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, stoich_data%total_tags-1
                 Do k=j+1, stoich_data%total_tags
                   If (Trim(model_data%block_component%tag(j))==Trim(model_data%block_component%tag(k))) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(model_data%block_component%tag(j)), 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the atomic components of the species must be declared,&
                                                    & each tag only once'
                     Call info(messages, 3)
                     Call error_stop(' ')
                   End If
                 End Do
              End Do  
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='amounts') Then
            If (.Not. header(2)) Then 
              i=i+1
              Call check_for_rubbish(iunit, '&block_input_composition') 
              Read (iunit, Fmt=*, iostat=io) word, (model_data%block_component%N0(j), j = 1, stoich_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the amount of atoms for each atomic tag'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),'". Valid options are "tags" and "amounts".&
                                         & Please refer to the manual'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')    'The correct structure for the block must be:'
          Write (messages(3),'(1x,a)')    '&block_input_composition'
          Write (messages(4),'(1x,a)')    '  tags      tg1    tg2    tg3   .... tgNsp'
          Write (messages(5),'(1x,a)')    '  amounts   N_tg1  N_tg2  N_tg3 .... N_tgNsp'
          Write (messages(6),'(1x,a)')    '&end_block_input_composition'
          Write (messages(7),'(1x,a,i3)') 'where, from &block_species_components, in this case Nsp = ', stoich_data%total_tags
          Write (messages(8),'(1x,a)')    'See manual for details'
          Call info(messages, 8)
          Call error_stop(' ')
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within the block'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&block_input_composition') 
      Call capital_to_lower_case(word)
      If (word /= '&end_block_input_composition') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Descriptors "tags" and "amounts" have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_block_input_composition.' 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

    ! Check consistency of the tags with those defined in block_species_components
    ! Assing the chemical elements to model_data%input%element
    Do i=1, stoich_data%total_tags
      loop=.True.
      j=1
      Do While (j <= stoich_data%num_species%value .And. loop)
        k=1
        Do While (k <= stoich_data%species(j)%num_components .And. loop)
          If (Trim(stoich_data%species(j)%component%tag(k))==model_data%block_component%tag(i)) Then
            model_data%block_component%element(i)=Trim(stoich_data%species(j)%component%element(k))
            loop=.False.
          End If   
          k=k+1
        End Do
        j=j+1        
      End Do 
      If (loop) Then
        Write (messages(2),'(3a)') 'Tag "', Trim(model_data%block_component%tag(i)), '" has not been defined in&
                                 & &block_species_components. Please check'
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End Do

 
    ! Check if the number of atoms are correct
    Do i=1, stoich_data%total_tags
      If (model_data%block_component%N0(i)< 0) Then
        Write (messages(2),'(3a)') 'Tag "', Trim(model_data%block_component%tag(i)), '" CANNOT have a negative number of&
                                 & atoms within the input structure! Please correct'
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End Do

   ! Calculate the number of total atoms set in the block
    model_data%block_component%numtot=0
    Do i=1, stoich_data%total_tags
      model_data%block_component%numtot=model_data%block_component%numtot+model_data%block_component%N0(i)
    End Do

   ! Assing atomic numbers
   Do i=1, stoich_data%total_tags
     loop=.True.
     k=1
     Do While (k <= NPTE .And. loop)
       If (Trim(chemsymbol(k))==Trim(model_data%block_component%element(i))) Then
         loop=.False.
         model_data%block_component%atomic_number(i)=k
       End If
       k=k+1
     End Do
   End Do

  End Subroutine read_block_input_composition


  Subroutine read_block_species_components(iunit, stoich_data, analysis)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read atomic composition from &block_species_components for 
    ! each chemical species defined in &block_species
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(stoich_type), Intent(InOut) :: stoich_data
    Character(Len=*),  Intent(In   ) :: analysis

    Integer(Kind=wi)   :: io, i, j, k, l, num_comp, length
    Character(Len=256) :: messages(8), word
    Character(Len=64 ) :: error_block_species_components
    Character(Len=64)  ::  separator
    Logical  :: endblock, species_read, loop, fread

    error_block_species_components = '***ERROR in &block_species_components of SET_EQCM file'
    Write (messages(1),'(a)') error_block_species_components
    Write (messages(3),'(1x,a)') ' '
    Write (messages(5),'(1x,a)')  'tag_1  element_1  N0_1  "separator"&
                                & tag_2  element_2  N0_2 .... "separator"&  
                                & tag_Nc  element_Nc  N0_Nc' 
    Write (messages(6),'(1x,a)') ' '
    Write (messages(7),'(1x,a)') 'where N0_i is the amount of atoms for each of the Nc component of the species "i". "separator"&
                                & can be any string to separate the specifications of each components for a better visualization'
    Write (messages(8),'(1x,a)') 'IMPORTANT: All tags must have maximum of 4 characters. Have you missed any of the NC components?&
                                & No comment is allowed in between the two senteces.'

    i=1
    Do While (i <= stoich_data%num_species%value)

      Read (iunit, Fmt=*, iostat=io) word
      If (io/=0) Then
        Exit 
      End If

      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          species_read=.False.
          j=1
          loop=.True.
          Do While (j <= stoich_data%num_species%value .And. loop)
            If (Trim(word)==Trim(stoich_data%species(j)%tag)) Then 
              species_read=.True.
              If (.Not. stoich_data%species(j)%component%fread) Then
                stoich_data%species(j)%component%fread = .True.
                loop=.False.
                Call check_for_rubbish(iunit, '&block_species_components')
                Read (iunit, Fmt=*, iostat=io) word, stoich_data%species(j)%num_components,&
                                                   stoich_data%species(j)%topology

                If (stoich_data%species(j)%num_components==1 .And. Trim(stoich_data%species(j)%topology)=='molecule')Then
                  Write (messages(2),'(1x,3a)') 'Wrong format for the atomic specification of species "', Trim(word), '"' 
                  Write (messages(3),'(1x,a)')  'A molecular species must have a number of atoms larger than 1'
                  Write (messages(4),'(1x,a)')  'If this species has only 1 atom, change its topology to "atom"'
                  Call info(messages, 4)
                  Call error_stop(' ')
                End If
 
 
                Write (messages(2),'(1x,3a)') 'Wrong format for the atomic specification of species "', Trim(word), &
                                             '". Format for each species must be:' 
                Write (messages(4),'(1x,a)')  'Species_tag  number_components(Nc)   topology   &
                                              &bond_cutoff (bond_cutoff only needed if topology is "molecule")' 
                If (io /= 0) Then
                  Call info(messages, 8)
                  Call error_stop(' ')
                End If 

                Call capital_to_lower_case(stoich_data%species(j)%topology) 

                If (stoich_data%species(j)%topology=='molecule') Then
                  Backspace iunit
                  Read (iunit, Fmt=*, iostat=io) word, stoich_data%species(j)%num_components,&
                                                   stoich_data%species(j)%topology, stoich_data%species(j)%bond_cutoff
                  If (io /= 0) Then
                    Write (messages(2),'(1x,3a)') 'Wrong (or missing) setting for bond_cutoff of species "', Trim(word), &
                                             '". Format for each species must be:' 
                        
                    Call info(messages, 8)
                    Call error_stop(' ')
                  End If 
                  ! Check if the bond_cutoff value is sensible
                  If (stoich_data%species(j)%bond_cutoff > max_intra) Then
                    Write (messages(2),'(1x,a,f4.2,3a,f4.2,a)') 'The maximum limit for the bonding distance of ',&
                                                          & stoich_data%species(j)%bond_cutoff,&
                                                          &' Angstrom set for species "', Trim(word), '" MUST NOT be larger than ',&
                                                          & max_intra, ' Angstrom. Please change'
                    Call info(messages,2)
                    Call error_stop(' ')
                  End If
                
                  If (stoich_data%species(j)%bond_cutoff < min_intra) Then
                    Write (messages(2),'(1x,a,f4.2,3a,f4.2,a)') 'The minimum limit for the bonding distance of ',&
                                                          & stoich_data%species(j)%bond_cutoff,&
                                                          &' Angstrom set for species "', Trim(word), '" MUST NOT be lower than ',&
                                                          & min_intra, ' Angstrom. Please change'
                    Call info(messages,2)
                    Call error_stop(' ')
                  End If

                End If

                If (Trim(stoich_data%species(j)%vartype0)=='fixed') Then 
                  If (Trim(analysis)=='model_cycled_sample' .Or. Trim(analysis)=='model_pristine_sample') Then 
                    If (Trim(stoich_data%species(j)%topology) /= 'crystal') Then
                      Write (messages(2),'(1x,3a)') 'The requested analysis is "', Trim(analysis), '"'
                      Write (messages(3),'(1x,3a)') 'Species ', Trim(word), ' has been defined as "fixed" in &Block_species. The&
                                              & topology for this species in &block_species_components must be set to "crystal"' 
                      Call info(messages, 3)
                      Call error_stop(' ')
                    End If 
                  Else If (Trim(analysis)=='model_disordered_system') Then
                    If (Trim(stoich_data%species(j)%topology) /= 'solute') Then
                      Write (messages(2),'(1x,3a)') 'The requested analysis is "', Trim(analysis), '"'
                      Write (messages(3),'(1x,3a)') 'Species ', Trim(word), ' has been defined as "fixed" in &Block_species. The&
                                              & topology for this species in &block_species_components must be set to "solute"' 
                      Call info(messages, 3)
                      Call error_stop(' ')
                    End If 
                  End If 
                End If 
     
                If ((Trim(stoich_data%species(j)%vartype0)=='dependent' .Or. Trim(stoich_data%species(j)%vartype0)=='independent')&
                  .And. Trim(stoich_data%species(j)%topology)=='crystal') Then
                  Write (messages(2),'(1x,5a)') 'Species ', Trim(word), ' has been defined as "',&
                                              & Trim(stoich_data%species(j)%vartype0),&
                                              & '" in &Block_species. The topology for this species in &block_species_components&
                                              & must be either "atom" or "molecule"' 
                  Call info(messages, 2)
                  Call error_stop(' ')
                End If 
 
                If (stoich_data%species(j)%num_components <= 0) Then
                  Write (messages(2),'(1x,3a)') 'The number of atomic components for species ', Trim(word), &
                                               ' must be larger than zero!'
                  Call info(messages, 2)
                  Call error_stop(' ')
                End If 

                num_comp=stoich_data%species(j)%num_components
                fread= .True.
                Do While (fread)
                  Read (iunit, Fmt=*, iostat=io) word
                  If (word(1:1)/='#') Then
                    fread=.False.
                    Call check_for_rubbish(iunit, '&block_species_components')
                  End If
                End Do

                Read (iunit, Fmt=*, iostat=io) (stoich_data%species(j)%component%tag(k), &
                                                   stoich_data%species(j)%component%element(k), &
                                                   stoich_data%species(j)%component%N0(k), separator, &   
                                                   k = 1, num_comp-1), &
                                                   stoich_data%species(j)%component%tag(num_comp), &
                                                   stoich_data%species(j)%component%element(num_comp), &
                                                   stoich_data%species(j)%component%N0(num_comp)
 
                If (io /= 0) Then
                  Call info(messages, 8)
                  Call error_stop(' ')
                End If 

                If (Trim(stoich_data%species(j)%vartype0)/='fixed') Then
                  If (stoich_data%species(j)%num_components==1 .And. stoich_data%species(j)%component%N0(1)==1 .And. &
                     stoich_data%species(j)%topology/='atom') Then
                    Write (messages(2),'(1x,3a)') 'It appears that the topology of species ', Trim(word), ' should be set&
                                                 & to "atom" in &block_species_components. Please check'
                    Call info(messages, 2)
                    Call error_stop(' ')
                  End If   

                  If (stoich_data%species(j)%component%N0(1)/=1 .And. stoich_data%species(j)%topology/='molecule') Then
                    Write (messages(2),'(1x,3a)') 'It appears that the topology of species ', Trim(word), ' should be set&
                                                 & to "molecule" in &block_species_components. Please check'
                    Call info(messages, 2)
                    Call error_stop(' ')
                  End If   
                End If

              Else
                Write (messages(2),'(3(1x,a))') 'Atomic description for species', &
                                               Trim(word), 'has been defined more than once'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
            End If
            j=j+1
          End Do

          If (.Not. species_read) Then
            Write (messages(2),'(3a)') 'Species "', Adjustl(Trim(word)), '" has NOT been defined in&
                             & &Block_species. Please check for typo (capitalization DOES matter&
                             & for the definition of the species)'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
          If (.Not. stoich_data%species(j-1)%component%fread) Then
            If (word(1:1)=="&") Then
              Write (messages(2),'(a)') 'Missing atomic description. Not all species defined&
                                       & in &Block_species have been included' 
              Call info(messages, 2)
              Call error_stop(' ')
            Else 
              Write (messages(2),'(3(1x,a))') 'Atomic description for species', &
                                             Trim(word), 'has been defined more than once'
              Call info(messages, 2)
              Call error_stop(' ')
            End IF
          End If
          i=i+1
        Else
          Write (messages(2),'(1x,a)') 'End of block found! Missing atomistic description for&
                                     & the species defined in &Block_species'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End If
    End Do 

    endblock=.False.
    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call capital_to_lower_case(word)
      If (word /= '&end_block_species_components') Then
        If (word(1:1) /= '#') Then
          If ((i-1)/=stoich_data%num_constraints%value) Then
            Write (messages(2),'(a)') 'Atomistic details for the involved species must be closed with&
                                    & sentence &end_block_species_components. Check 1) syntax 2) if the&
                                    & number of specifications is larger than the value assigned to "number_species"'
          Else
            Write (messages(2),'(a)') 'The number of specifications for the composition of the involved&
                                     & species is larger than the number of species (see&
                                     & value for "number_species" in &Block_species)'
          End If
          Call info(messages,2)
          Call error_stop(' ')
        End If
      Else
          endblock=.True.
      End If
    End Do

    !Check that tags for compoenents do not exceed 4 characters
    Do i=1, stoich_data%num_species%value
      Do j=1, stoich_data%species(i)%num_components
        Call get_word_length(stoich_data%species(i)%component%tag(j),length)
        If (length > 4) Then
          Write (messages(2),'(5a)') 'Tag "', Trim(stoich_data%species(i)%component%tag(j)), '" for atomic component&
                                   & of species "', Trim(stoich_data%species(i)%tag), '" must have a maximum of 4 characters.&
                                   & Please use a different tag'
          Call info(messages,2)
          Call error_stop(' ')  
        End If
      End Do
    End Do
  
    ! Within the same species
    Do i=1, stoich_data%num_species%value
      Do j=1, stoich_data%species(i)%num_components
        Do k=1, stoich_data%species(i)%num_components
          If (j/=k) Then
            If (Trim(stoich_data%species(i)%component%tag(j))==Trim(stoich_data%species(i)%component%tag(k))) Then
              Write (messages(2),'(1x,5a)') 'Component tag "', Trim(stoich_data%species(i)%component%tag(k)), &
                             &'" for species "', Trim(stoich_data%species(i)%tag), '" has been used more than once.&  
                             & Tags for components must be unequivocaly defined.'
              Call info(messages, 2)
              Call error_stop(' ')  
            End If 
          End If
        End Do
      End Do 
    End Do
    ! Check if tags for components are unequivocaly defined
    ! In other species
    Do i=1, stoich_data%num_species%value
      Do l=1, stoich_data%num_species%value
        If (i/=l) Then
          Do j=1, stoich_data%species(i)%num_components
            Do k=1, stoich_data%species(l)%num_components
                If (Trim(stoich_data%species(i)%component%tag(j))==Trim(stoich_data%species(l)%component%tag(k))) Then
                  Write (messages(2),'(1x,7a)') 'Component tag "', Trim(stoich_data%species(l)%component%tag(k)), &
                                 &'" has been used for species "', Trim(stoich_data%species(i)%tag), '" and species "',&  
                                 & Trim(stoich_data%species(l)%tag), '". Tags for components must be unequivocaly defined.'
                  Call info(messages, 2)
                  Call error_stop(' ')  
                End If 
            End Do
          End Do 
        End If
      End Do
    End Do

    ! Calculate the total number of different tags
    stoich_data%total_tags=0
    Do i=1, stoich_data%num_species%value
      stoich_data%total_tags=stoich_data%total_tags+stoich_data%species(i)%num_components
    End Do

    ! Calculate total number of atoms per species
    Do i=1, stoich_data%num_species%value
      stoich_data%species(i)%atoms_per_species=0
      Do j=1,stoich_data%species(i)%num_components
        stoich_data%species(i)%atoms_per_species = stoich_data%species(i)%atoms_per_species + &
                                                   stoich_data%species(i)%component%N0(j)
      End Do
    End Do


  End Subroutine read_block_species_components

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutines to read parameters and settings for stoichiometric analysis 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine read_individual_constraint(iunit, stoich_data, ic)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read each individual constraint 
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(stoich_type), Intent(InOut) :: stoich_data
    Integer(Kind=wi),  Intent(In   ) :: ic

    Integer(Kind=wi) :: io
    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_block_constraints_species

    error_block_constraints_species = '***ERROR in &block_constraints_species of SET_EQCM file'
    Write (messages(1),'(a)') error_block_constraints_species

    Call capital_to_lower_case(stoich_data%constraints(ic)%type)

    If (Trim(stoich_data%constraints(ic)%type)=='target_value' .Or. &
       Trim(stoich_data%constraints(ic)%type)=='ratio_fixed') Then
      Read (iunit, Fmt=*, iostat=io) stoich_data%constraints(ic)%type, stoich_data%constraints(ic)%leg,   &
                                        stoich_data%constraints(ic)%tag0_species, &
                                        stoich_data%constraints(ic)%value(1)
   
    Else If (Trim(stoich_data%constraints(ic)%type)=='target_min' .Or. &
            Trim(stoich_data%constraints(ic)%type)=='target_max') Then
      Read (iunit, Fmt=*, iostat=io) stoich_data%constraints(ic)%type, stoich_data%constraints(ic)%leg,   &
                                        stoich_data%constraints(ic)%tag0_species
 
    ElseIf (Trim(stoich_data%constraints(ic)%type)=='target_range') Then
      Read (iunit, Fmt=*, iostat=io) stoich_data%constraints(ic)%type, stoich_data%constraints(ic)%leg,   &
                                        stoich_data%constraints(ic)%tag0_species, &
                                        stoich_data%constraints(ic)%value(1), stoich_data%constraints(ic)%value(2)  
  
    ElseIf (Trim(stoich_data%constraints(ic)%type)=='keep_ratio') Then
      Read (iunit, Fmt=*, iostat=io) stoich_data%constraints(ic)%type, stoich_data%constraints(ic)%leg,   &
                                        stoich_data%constraints(ic)%tag0_species
 
    ElseIf (stoich_data%constraints(ic)%type(1:1)=='#') Then
      Write (messages(2),'(a)') 'Comments are not allowed within blocks'
      Call info(messages,2)
      Call error_stop(' ')
    Else If (Trim(stoich_data%constraints(ic)%type) == '&end_block_constraints_species') Then
      Write (messages(2),'(2(a,i2),a)') 'Missing specification for constraints: only ',  ic-1,&
                       &' constraint set out of ', stoich_data%num_constraints%value,&
                        ' (see number_constraints).'      
      Call info(messages,2) 
      Call error_stop(' ')
 
    Else 
      Write (messages(2),'(3a,i2,a)') 'Setting "', Trim(stoich_data%constraints(ic)%type), '" for constraint ', ic, &
                                     &' is not valid. Possible options are: target_value, target_min, target_max,&
                                     & ratio_fixed and keep_ratio (see manual).'
      Call info(messages,2)
      Call error_stop(' ') 
    End If

    Call capital_to_lower_case(stoich_data%constraints(ic)%type)
    Call capital_to_lower_case(stoich_data%constraints(ic)%leg)

    stoich_data%constraints(ic)%tag_species=stoich_data%constraints(ic)%tag0_species
    Call set_read_status('constraint', io, stoich_data%constraints(ic)%fread, stoich_data%constraints(ic)%fail)


  End Subroutine read_individual_constraint

  Subroutine read_block_constraints_species(iunit, stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read settings inside block_constraints, in case the user
    ! wants to reduce the phase space of solutions
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Type(stoich_type),  Intent(InOut) :: stoich_data 

    Integer(Kind=wi)   ::  io, i
    Character(Len=256) :: word, messages(2)
    Character(Len=64 ) :: error_block_constraints_species
    Logical  :: error, endblock, fread

    error= .False.

    error_block_constraints_species = '***ERROR in &block_constraints_species of SET_EQCM file'
    Write (messages(1),'(a)') error_block_constraints_species

    fread= .True.
    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      If (word(1:1)/='#') Then 
        fread=.False.
        Call check_for_rubbish(iunit, '&block_constraints_species')
      End If
    End Do

    Read (iunit, Fmt=*, iostat=io) stoich_data%num_constraints%tag, stoich_data%num_constraints%value
    If (io/=0) Then
      Write (messages(2),'(a)') 'Problems to read "number_constraints". It must be the first defined&
                                & directive. See manual'
      Call info(messages,2)
      Call error_stop(' ')
    End If  
    Call set_read_status(word, io, stoich_data%num_constraints%fread, stoich_data%num_constraints%fail)

    If (Trim(stoich_data%num_constraints%tag) /= 'number_constraints') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(stoich_data%num_constraints%tag), &
                         & '" has been found, but directive "number_constraints" is expected.'
      error=.True.
    End If

    If (stoich_data%num_constraints%fail) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "number_constraints"'
      error=.True.
    Else
      If (stoich_data%num_constraints%value<1) Then
        Write (messages(2),'(a)') 'Number of constraints for species MUST be >= 1'
        error=.True.
      End If
    End If

    If (error) Then
      Call info(messages,2)
      Call error_stop(' ')
    End If
 
    ! Initialise arrays     
    Call stoich_data%init_constraints_arrays()

    i=1
    Do While (i <= stoich_data%num_constraints%value)
      Read (iunit, Fmt=*, iostat=io) word
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&block_constraints_species')
        Read (iunit, Fmt=*, iostat=io) stoich_data%constraints(i)%type
        Backspace iunit
        Call read_individual_constraint(iunit,stoich_data,i)
        i=i+1
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call capital_to_lower_case(word)
      If (word /= '&end_block_constraints_species') Then
        If (word(1:1) /= '#') Then
          If ((i-1)/=stoich_data%num_constraints%value) Then 
            Write (messages(2),'(a)') 'Number of defined constraints is larger than&
                                & the value given by directive "number_constraints"'
          Else
            Write (messages(2),'(a)') 'Constraints block (for species) must be closed with&
                                    & sentence &end_block_constraints_species. Please check if the&
                                    & number of contraints defined is larger than the value assigned&
                                    & to "number_constraints".'
          End If   
          Call info(messages,2) 
          Call error_stop(' ')
        End If  
      Else
          endblock=.True.
      End If
    End Do

  End Subroutine read_block_constraints_species


  Subroutine read_block_species(iunit, stoich_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read settings inside block_species, needed for 
    ! stoichimetric analysis 
    !
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(stoich_type), Intent(InOut) :: stoich_data 

    Integer(Kind=wi)  ::  io, i

    Character(Len=256)  :: word, messages(2)
    Character(Len=64 )  :: error_block_species
    Logical  :: error, endblock, fread 

    error= .False.
    error_block_species = '***ERROR in &block_species of SET_EQCM file'
    Write (messages(1),'(a)') error_block_species 

    fread= .True.
    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&block_species')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&block_species')
      End If
    End Do

    ! Read number of species
    Read (iunit, Fmt=*, iostat=io) stoich_data%num_species%tag, stoich_data%num_species%value
    Call set_read_status(word, io, stoich_data%num_species%fread, stoich_data%num_species%fail)

    If (Trim(stoich_data%num_species%tag) /= 'number_species') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(stoich_data%num_species%tag), &
                         & '" has been found, but directive "number_species" is expected.'
      error=.True.
    End If 

    If (stoich_data%num_species%fail) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "number_species"'
      error=.True.
    Else  
      If (stoich_data%num_species%value<1) Then
        Write (messages(2),'(a)') 'Number of stoichiometric species MUST be >= 1'
        error=.True.
      ElseIf (stoich_data%num_species%value>max_species) Then
        Write (messages(2),'(a,i3,a)') 'Are you sure you want to consider a system involving more than ', max_species,&
                                    & ' different types of species?. We have set this value as a sensible limit...'
        error=.True.
      End If
    End If

    If (error) Then
      Call info(messages,2) 
      Call error_stop(' ')
    End If


    ! Initialise arrays
    Call stoich_data%init_species_arrays('block')   

    i=1
    Do While (i <= stoich_data%num_species%value)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&block_species')
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&block_species')
  
        Read (iunit, Fmt=*, iostat=io) stoich_data%species(i)%tag, stoich_data%species(i)%mass,   &
                                            stoich_data%species(i)%ox, stoich_data%species(i)%s0, &
                                            stoich_data%species(i)%vartype0

        If (io/=0) Then
          If (Trim(stoich_data%species(i)%tag) == '&end_block_species') Then
            Write (messages(2),'(2(a,i2),a)') 'Missing specification for chemical species. Only ',  i-1,&
                                    &' species set out of ', stoich_data%num_species%value, ' (see number_species)'      
            Call info(messages,2) 
            Call error_stop(' ')
          End If 
        End If

        Call capital_to_lower_case(stoich_data%species(i)%vartype0)

        !Assign s0 to the pristine sample
        stoich_data%species(i)%s0_pristine=stoich_data%species(i)%s0
        Call set_read_status(word, io, stoich_data%species(i)%fread, stoich_data%species(i)%fail)
        i=i+1
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&block_species')
      Call capital_to_lower_case(word)
      If (word /= '&end_block_species') Then
        If (word(1:1) /= '#') Then
          If ((i-1)/=stoich_data%num_species%value) Then 
            Write (messages(2),'(a)') 'Number of chemical species specified is larger than&
                                & the value given by directive "number_species"'
          Else
            Write (messages(2),'(a)') 'Block must be closed with sentence &end_block_species. Please check.&
                                     & Is the number of defined species the same as set for directive "number_species"?'
          End If   
          Call info(messages,2) 
          Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If  
    End Do

  End Subroutine read_block_species

  Subroutine read_hpc_settings(iunit, hpc_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read hpc directives, which will be used to generate the 
    ! script files to submit the simulations 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Type(hpc_type),   Intent(InOut) :: hpc_data

    Character(Len=256) :: word, word2
    Integer(Kind=wi)   :: length, io, i
    Integer(Kind=wi)   :: wl, wl2
  
    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &block_hpc_settings -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&End_Block_hpc_settings" to close the block. Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_block_hpc_settings') Exit
      Backspace iunit


      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'number_mpi_tasks') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%number_mpi_tasks%value
        Call set_read_status(word, io, hpc_data%number_mpi_tasks%fread, hpc_data%number_mpi_tasks%fail)

      Else If (word(1:length) == 'number_nodes') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%number_nodes%value
        Call set_read_status(word, io, hpc_data%number_nodes%fread, hpc_data%number_nodes%fail)

      Else If (word(1:length) == 'cpus_per_node') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%cpus_per_node%value
        Call set_read_status(word, io, hpc_data%cpus_per_node%fread, hpc_data%cpus_per_node%fail)

      Else If (word(1:length) == 'threads_per_process') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%threads_per_process%value
        Call set_read_status(word, io, hpc_data%threads_per_process%fread, hpc_data%threads_per_process%fail)

      Else If (word(1:length) == 'machine_name') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%machine_name%type
        Call set_read_status(word, io, hpc_data%machine_name%fread, hpc_data%machine_name%fail,&
                           & hpc_data%machine_name%type)

      Else If (word(1:length) == 'platform') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%platform%type
        Call set_read_status(word, io, hpc_data%platform%fread, hpc_data%platform%fail, hpc_data%platform%type)

      Else If (word(1:length) == 'parallelism_type') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%parallelism_type%type
        Call set_read_status(word, io, hpc_data%parallelism_type%fread, hpc_data%parallelism_type%fail, &
                           & hpc_data%parallelism_type%type)

      Else If (word(1:length) == 'job_name') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%job_name%type
        Call set_read_status(word, io, hpc_data%job_name%fread, hpc_data%job_name%fail)

      Else If (word(1:length) == 'project_name') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%project_name%type
        Call set_read_status(word, io, hpc_data%project_name%fread, hpc_data%project_name%fail)

      Else If (word(1:length) == 'queue') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%queue%type
        Call set_read_status(word, io, hpc_data%queue%fread, hpc_data%queue%fail, hpc_data%queue%type)

      Else If (word(1:length) == 'memory_per_cpu') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%memory_per_cpu%value
        Call set_read_status(word, io, hpc_data%memory_per_cpu%fread, hpc_data%memory_per_cpu%fail)

      Else If (word(1:length) == 'time_limit') Then 
        Read (iunit, Fmt=*, iostat=io) word, (hpc_data%time_limit%value(i), i=1,3)
        Call set_read_status(word, io, hpc_data%time_limit%fread, hpc_data%time_limit%fail)

      Else If (word(1:length) == 'mkl') Then
        Read (iunit, Fmt=*, iostat=io) word, hpc_data%mkl%stat
        Call set_read_status(word, io, hpc_data%mkl%fread, hpc_data%mkl%fail)

      Else If (word(1:length) == '&modules') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, hpc_data%modules_info%fread, hpc_data%modules_info%fail)
        Call read_hpc_modules(iunit, hpc_data)

      Else If (word(1:length) == 'executable') Then
        Read (iunit, Fmt=*, iostat=io) word2
        Backspace iunit
        Read (iunit, Fmt='(a)', iostat=io) word
        Call set_read_status(word, io, hpc_data%executable%fread, hpc_data%executable%fail)

        word=Adjustl(Trim(word))
        wl2= Len(Trim(word2))
        Do i=1, wl2
          word(i:i)=' '
        End Do
        word=Adjustl(Trim(word))
        wl=Index(word, ' ')
        wl2=Len(word)
        Do i= wl, wl2
          word(i:i)=' '
        End Do
        hpc_data%executable%type=Trim(word)

      Else If (word(1:length) == 'exec_options') Then
        Read (iunit, Fmt=*, iostat=io) word2
        Backspace iunit
        Read (iunit, Fmt='(a)', iostat=io) word
        Call set_read_status(word, io, hpc_data%exec_options%fread, hpc_data%exec_options%fail)

        word=Adjustl(Trim(word))
        wl2= Len(Trim(word2))
        Do i=1, wl2
          word(i:i)=' '
        End Do
        word=Adjustl(Trim(word))
        hpc_data%exec_options%type=Trim(word)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word), '" is not recognised as a valid hpc settings.',&
                                & ' See manual. Have you properly closed the block with "&End_Block_hpc_settings"?'
        Call error_stop(message)
      End If

     End Do

  End Subroutine read_hpc_settings

  Subroutine read_hpc_modules(iunit, hpc_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the names of the modules needed to setup the script
    ! file for HPC job submission
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Type(hpc_type),   Intent(InOut) :: hpc_data

    Character(Len=256) :: word
    Integer(Kind=wi)   :: i, j
    Integer(Kind=wi)   :: io, length, total_length

    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &modules (inside &block_hpc_settings):'

    i=1
    Do 
      Read (iunit, Fmt=*, iostat=io) word
      If (is_iostat_end(io)) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'End of file? Have you closed the block with &end_modules?'
        Call error_stop(message)
      End If

      Call capital_to_lower_case(word) 

      If (Trim(word) == '&end_modules') Then
        Exit
      Else
        If (word(1:1)/='#') Then
          If (word(1:1)=='&') Then
            Write (message,'(2a)') Trim(set_error), ' It appears that &modules has not been closed properly. Please use&
                                 & &end_modules.'
            Call error_stop(message) 
          Else 
            Backspace iunit 
            Read (iunit, Fmt='(a)', iostat=io) word 
            hpc_data%modules%element(i)=Adjustl(Trim(word))
            i=i+1
          End If  
        End If
      End If
    End Do
    hpc_data%modules%num=i-1

    total_length=Len(hpc_data%modules%element(1))
    ! Clean modules if there were comments inserted
    Do i=1, hpc_data%modules%num
      Call get_word_length(hpc_data%modules%element(i),length)
      Do j= length+1, total_length
        hpc_data%modules%element(i)(j:j)=' '
      End Do
    End Do

  End Subroutine read_hpc_modules

  Subroutine read_simulation_settings(iunit, simulation_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read simulation directives, which will be used to generate input 
    ! files, required for atomistic level simulation
    !
    ! author    - i.scivetti June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io
  
    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &Block_simulation_settings -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&End_Block_simulation_settings" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_block_simulation_settings') Exit
      Call check_for_rubbish(iunit, '&block_simulation_settings')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'simulation_type') Then 
        Read (iunit, Fmt=*, iostat=io) word, simulation_data%simulation%type
        Call set_read_status(word, io, simulation_data%simulation%fread, simulation_data%simulation%fail, &
                           & simulation_data%simulation%type)

      Else If (word(1:length) == 'theory_level') Then 
        Read (iunit, Fmt=*, iostat=io) word, simulation_data%theory_level%type
        Call set_read_status(word, io, simulation_data%theory_level%fread, simulation_data%theory_level%fail, &
                           & simulation_data%theory_level%type)

      Else If (word(1:length) == 'net_charge') Then
       Read (iunit, Fmt=*, iostat=io) word, simulation_data%net_charge%value
       Call set_read_status(word, io, simulation_data%net_charge%fread, simulation_data%net_charge%fail)

      Else If (word(1:length) == '&dft_settings') Then
        Read (iunit, Fmt=*, iostat=io) word 
        If (simulation_data%dft%generate) Then
          Call duplication_error(word)
        End If
        simulation_data%dft%generate=.True.
        Call simulation_data%init_input_dft_variables()
        ! Now, it is ready to read information inside &dft_settings
        Call read_dft_settings(iunit, simulation_data)

      Else If (word(1:length) == '&motion_settings') Then
        Read (iunit, Fmt=*, iostat=io) word 
        If (simulation_data%motion%generate) Then
          Call duplication_error(word)
        End If
        simulation_data%motion%generate=.True.
        ! Now, it is ready to read information inside &dft_settings
        Call simulation_data%init_input_motion_variables()
        Call read_motion_settings(iunit, simulation_data)

      ! Extra directives
      Else If (word(1:length) == '&extra_directives') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, simulation_data%extra_info%fread, simulation_data%extra_info%fail)
        simulation_data%extra_info%stat = .True.
        ! Read extra input 
        Call read_extra_directives(iunit, simulation_data)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid simulation settings.',&
                                & ' See manual. Have you properly closed the block with "&End_block_simulation_settings"?'
        Call error_stop(message)
      End If

    End Do
  End Subroutine read_simulation_settings

  Subroutine read_extra_directives(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read extra directives defined by the user
    !
    ! author    - i.scivetti Jul 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Character(Len=256)  :: word
    Integer(Kind=wi)    :: i, io

    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &modules (inside &extra_directives):'

    i=1
    Do 
      Read (iunit, Fmt='(a)', iostat=io) word
      If (is_iostat_end(io)) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'End of file? Have you closed the block with "&end_extra_directives"?'
        Call error_stop(message)
      End If

      word=Trim(Adjustl(word))
      Call capital_to_lower_case(word) 

      If (Trim(word) == '&end_extra_directives') Then
        Exit
      Else
        If (word(1:1)=='&') Then
          If (word(1:7)/= '&block ' .And. word(1:10)/='&endblock ') Then
            Write (message,'(2a)') Trim(set_error), ' It appears that &extra_directives has not been closed properly. Please use&
                               & "&end_extra_directives".'
            Call error_stop(message) 
          Else
            Backspace iunit 
            Read (iunit, Fmt='(a)', iostat=io) word 
            simulation_data%extra_directives%array(i)=Adjustl(Trim(word))
            i=i+1
          End If
        Else 
          Backspace iunit 
          Read (iunit, Fmt='(a)', iostat=io) word 
          simulation_data%extra_directives%array(i)=Adjustl(Trim(word))
          i=i+1
        End If  
      End If
    End Do
    simulation_data%extra_directives%N0=i-1

  End Subroutine read_extra_directives


  Subroutine read_motion_settings(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read directives from &motion_settings
    !
    ! author    - i.scivetti June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: T

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io
 
    Integer             :: j
    Character(Len=256)  :: message
    Character(Len=265)  :: set_error

    set_error = '***ERROR in &motion_settings -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_motion_settings" to close the block. Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_motion_settings') Then
        Exit
      End If
      Call check_for_rubbish(iunit, '&motion_settings')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'relax_method') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%motion%relax_method%type
        Call set_read_status(word, io, T%motion%relax_method%fread, T%motion%relax_method%fail, T%motion%relax_method%type)

      Else If (word(1:length) == 'ensemble') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%motion%ensemble%type
        Call set_read_status(word, io, T%motion%ensemble%fread, T%motion%ensemble%fail, T%motion%ensemble%type)

      Else If (word(1:length) == 'change_cell_volume') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%motion%change_cell_volume%stat
        Call set_read_status(word, io, T%motion%change_cell_volume%fread, T%motion%change_cell_volume%fail)

      Else If (word(1:length) == 'change_cell_shape') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%motion%change_cell_shape%stat
        Call set_read_status(word, io, T%motion%change_cell_shape%fread, T%motion%change_cell_shape%fail)

      Else If (word(1:length) == 'force_tolerance') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%delta_f%value(1), (T%motion%delta_f%units(j), j=1,2)
        Call set_read_status(word, io, T%motion%delta_f%fread, T%motion%delta_f%fail)
        Call capital_to_lower_case(T%motion%delta_f%units(1))
        Call capital_to_lower_case(T%motion%delta_f%units(2))

      Else If (word(1:length) == 'timestep') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%timestep%value, T%motion%timestep%units
        Call set_read_status(word, io, T%motion%timestep%fread, T%motion%timestep%fail)
        Call capital_to_lower_case(T%motion%timestep%units)

      Else If (word(1:length) == 'ion_steps') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%ion_steps%value
        Call set_read_status(word, io, T%motion%ion_steps%fread, T%motion%ion_steps%fail)

      Else If (word(1:length) == 'pressure') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%pressure%value, T%motion%pressure%units
        Call set_read_status(word, io, T%motion%pressure%fread, T%motion%pressure%fail)
        Call capital_to_lower_case(T%motion%pressure%units)

      Else If (word(1:length) == 'temperature') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%temperature%value, T%motion%temperature%units
        Call set_read_status(word, io, T%motion%temperature%fread, T%motion%temperature%fail)
        Call capital_to_lower_case(T%motion%temperature%units)

      Else If (word(1:length) == 'thermostat') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%thermostat%type
        Call set_read_status(word, io, T%motion%thermostat%fread, T%motion%thermostat%fail, T%motion%thermostat%type)

      Else If (word(1:length) == 'relax_time_thermostat') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%relax_time_thermostat%value, T%motion%relax_time_thermostat%units
        Call set_read_status(word, io, T%motion%relax_time_thermostat%fread, T%motion%relax_time_thermostat%fail)
        Call capital_to_lower_case(T%motion%relax_time_thermostat%units)

      Else If (word(1:length) == 'barostat') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%barostat%type
        Call set_read_status(word, io, T%motion%barostat%fread, T%motion%barostat%fail, T%motion%barostat%type)

      Else If (word(1:length) == 'relax_time_barostat') Then
        Read (iunit, Fmt=*, iostat=io) word, T%motion%relax_time_barostat%value, T%motion%relax_time_barostat%units
        Call set_read_status(word, io, T%motion%relax_time_barostat%fread, T%motion%relax_time_barostat%fail)
        Call capital_to_lower_case(T%motion%relax_time_barostat%units)

      ! Masses
      Else If (word(1:length) == '&masses') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%motion%mass_info%fread, T%motion%mass_info%fail)
        T%motion%mass_info%stat = .True.
        !Read masses from block &mass
        Call read_masses(iunit, T)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word), '" is not recognised as a valid motion settings.',&
                                & ' See manual. Have you properly closed the block with "&end_motion_settings"?'
        Call error_stop(message)
      End If

     End Do

  End Subroutine read_motion_settings

  Subroutine read_masses(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the masses for each atomic tags (this is strictly needed for 
    ! MD simulations with isotopes)
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Type(simul_type),   Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: io, i, j, k
    Character(Len=256)  :: messages(8), word
    Character(Len=64 )  :: set_error
    Logical  :: endblock
    Logical  :: header(2), error_duplication

    set_error = '***ERROR in &masses (inside &motion_settings):'
    Write (messages(1),'(a)') set_error

    header=.False.
    error_duplication=.False.

    i=1

    Do While (i <= 2)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&masses (inside &motion_settings)')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&masses (inside &motion_settings)')
          Call capital_to_lower_case(word) 
          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%motion%mass(j)%tag, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, simulation_data%total_tags-1
                 Do k=j+1, simulation_data%total_tags
                   If (Trim(simulation_data%motion%mass(j)%tag)==Trim(simulation_data%motion%mass(k)%tag)) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(simulation_data%motion%mass(j)%tag),&
                                                  & 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the components of the species must be declared,&
                                                   & each tag only once'
                     Call info(messages, 3)
                     Call error_stop(' ')
                   End If
                 End Do
              End Do  
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='values') Then
            If (.Not. header(2)) Then 
              i=i+1
              Read (iunit,Fmt=*,iostat=io) word, (simulation_data%motion%mass(j)%value, j=1,simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read masses values for the defined tags. Please check if tags&
                                            & are consistent with those defined in "&Block_Species_Components".'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),'". Valid options are "tags" and "values"&
                                         & Please refer to the manual.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')      'The correct structure must be:'
          Write (messages(3),'(1x,a)')      '&masses'
          Write (messages(4),'(1x,a)')      '  tags      tg1        tg2      tg3    ....   tgNsp'
          Write (messages(5),'(1x,a)')      '  values  mass_tg1  mass_tg2  mass_tg3 .... mass_tgNsp'
          Write (messages(6),'(1x,a)')      '&end_masses'
          Write (messages(7),'(1x,a,i3,a)') 'where in this case Nsp = ', simulation_data%total_tags,&
                                          & ', which corresponds to the number of tags defined in "&Block_Species_Components".'
          Write (messages(8),'(1x,a)')      'See manual for details'
          Call info(messages, 8)
          Call error_stop(' ')
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within &masses'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&masses (inside &motion_settings)')
      Call capital_to_lower_case(word)
      If (word /= '&end_masses') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Descriptors "tags" and "values" have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_masses.' 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

  End Subroutine read_masses

  Subroutine read_dft_settings(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read DFT directives from &DFT_settings
    !
    ! author         - i.scivetti May  2021
    ! contribution   - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: T 

    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io
    Integer(Kind=wi)   :: i
  
    Character(Len=256) :: message
    Character(Len=265) :: set_error

    set_error = '***ERROR in &dft_settings -'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&dft_settings')
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&End_dft_settings" to close the block. Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_dft_settings') Then
        Exit
      End If
      Call check_for_rubbish(iunit, '&dft_settings')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'xc_level') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%xc_level%type
        Call set_read_status(word, io, T%dft%xc_level%fread, T%dft%xc_level%fail, T%dft%xc_level%type)

      Else If (word(1:length) == 'xc_version') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%xc_version%type
        Call set_read_status(word, io, T%dft%xc_version%fread, T%dft%xc_version%fail, T%dft%xc_version%type)

      Else If (word(1:length) == 'vdw') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%vdw%type
        Call set_read_status(word, io, T%dft%vdw%fread, T%dft%vdw%fail, T%dft%vdw%type)

      Else If (word(1:length) == 'spin_polarised') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%spin_polarised%stat
        Call set_read_status(word, io, T%dft%spin_polarised%fread, T%dft%spin_polarised%fail)

      Else If (word(1:length) == 'energy_cutoff') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%encut%value, T%dft%encut%units
        Call set_read_status(word, io, T%dft%encut%fread, T%dft%encut%fail)
        Call capital_to_lower_case(T%dft%encut%units)
      
      Else If (word(1:length) == 'precision') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%precision%type
        Call set_read_status(word, io, T%dft%precision%fread, T%dft%precision%fail, T%dft%precision%type)

      Else If (word(1:length) == 'smearing') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%smear%type
        Call set_read_status(word, io, T%dft%smear%fread, T%dft%smear%fail, T%dft%smear%type)

      Else If (word(1:length) == 'width_smear') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%width_smear%value, T%dft%width_smear%units
        Call set_read_status(word, io, T%dft%width_smear%fread, T%dft%width_smear%fail)
        Call capital_to_lower_case(T%dft%width_smear%units)
      
      Else If (word(1:length) == 'scf_energy_tolerance') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%delta_e%value, T%dft%delta_e%units
        Call set_read_status(word, io, T%dft%delta_e%fread, T%dft%delta_e%fail)
        Call capital_to_lower_case(T%dft%delta_e%units)

      Else If (word(1:length) == 'scf_steps') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%scf_steps%value
        Call set_read_status(word, io, T%dft%scf_steps%fread, T%dft%scf_steps%fail)

      Else If (word(1:length) == 'kpoints') Then 
        Read (iunit, Fmt=*, iostat=io) word, T%dft%kpoints%tag, (T%dft%kpoints%value(i), i=1,3)
        Call set_read_status(word, io, T%dft%kpoints%fread, T%dft%kpoints%fail, T%dft%kpoints%tag)

      Else If (word(1:length) == 'max_l_orbital') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%max_l_orbital%value
        Call set_read_status(word, io, T%dft%max_l_orbital%fread, T%dft%max_l_orbital%fail)

      Else If (word(1:length) == 'total_magnetization') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%total_magnetization%value
        Call set_read_status(word, io, T%dft%total_magnetization%fread, T%dft%total_magnetization%fail)

      Else If (word(1:length) == 'ot') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%ot%stat
        Call set_read_status(word, io, T%dft%ot%fread, T%dft%ot%fail)

      Else If (word(1:length) == 'edft') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%edft%stat
        Call set_read_status(word, io, T%dft%edft%fread, T%dft%edft%fail)

      Else If (word(1:length) == 'npar') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%npar%value
        Call set_read_status(word, io, T%dft%npar%fread, T%dft%npar%fail)

      Else If (word(1:length) == 'kpar') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%kpar%value
        Call set_read_status(word, io, T%dft%kpar%fread, T%dft%kpar%fail)

      Else If (word(1:length) == 'bands') Then
        Read (iunit, Fmt=*, iostat=io) word, T%dft%bands%value
        Call set_read_status(word, io, T%dft%bands%fread, T%dft%bands%fail)

      ! Pseudopotentials
      Else If (word(1:length) == '&pseudo_potentials') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%pp_info%fread, T%dft%pp_info%fail)
        T%dft%pp_info%stat=.True.
        ! Read names of the pseudo potentials 
        Call read_pseudo_poptentials(iunit, T)

      ! Basis set
      Else If (word(1:length) == '&basis_set') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%basis_info%fread, T%dft%basis_info%fail)
        T%dft%basis_info%stat=.True.
        ! Read type of basis 
        Call read_basis_set(iunit, T)

      ! Magnetization    
      Else If (word(1:length) == '&magnetization') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%mag_info%fread, T%dft%mag_info%fail)
        T%dft%mag_info%stat = .True.
        ! Read names of the pseudo potentials 
        Call read_dft_magnetization(iunit, T)

      ! Hubbard    
      Else If (word(1:length) == '&hubbard') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%hubbard_info%fread, T%dft%hubbard_info%fail)
        T%dft%hubbard_info%stat = .True.
        ! Read Hubbard corrections
        Call read_dft_hubbard(iunit, T)

      ! NGWF
      Else If (word(1:length) == '&ngwf') Then 
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, T%dft%ngwf_info%fread, T%dft%ngwf_info%fail)
        T%dft%ngwf_info%stat = .True.
        ! Read NGWF setings 
        Call read_dft_ngwf(iunit, T)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word), '" is not recognised as a valid DFT settings.',&
                                & ' See manual. Have you properly closed the block with "&end_dft_settings"?'
        Call error_stop(message)
      End If

     End Do

  End Subroutine read_dft_settings

  Subroutine read_basis_set(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the type of basis assigned to each atomic site
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Character(Len=256) :: word
    Integer(Kind=wi)   :: io
    Integer(Kind=wi)   :: i
    Logical            :: endblock

    Character(Len=256) :: message
    Character(Len=265) :: set_error

    set_error = '***ERROR in &basis_set (inside &dft_settings):'

    i=1
    Do While (i<= simulation_data%total_tags)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&basis_set (inside &dft_settings)')
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&basis_set (inside &dft_settings)')
        Read (iunit, Fmt=*, iostat=io) simulation_data%dft%basis_set(i)%tag, simulation_data%dft%basis_set(i)%type
        Call capital_to_lower_case(simulation_data%dft%basis_set(i)%type)
        i=i+1  
      End If
    End Do

    If (i-1 < simulation_data%total_tags) Then
      Write (message,'(2a)') Trim(set_error), ' The number of declared atomic tags is less than those defined in&
                          & "&block_input_composition". Please define all tags with the corresponding&
                          & basis set type'
      Call error_stop(message)
    End If

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&basis_set (inside &dft_settings)')
      Call capital_to_lower_case(word)
      If (word /= '&end_basis_set') Then
        If (word(1:1) /= '#') Then
          Write (message,'(2a)') Trim(set_error), ' It seems the user has provided wrong or additional information.&
                               & This must be closed with sentence "&end_basis_set"'
          Call error_stop(message)
        End If
      Else
          endblock=.True.
      End If
    End Do

  End Subroutine read_basis_set

  Subroutine read_pseudo_poptentials(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the names of the pseudo potentials files, which
    ! are required to build files for atomistic level simulation
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Character(Len=256) :: word
    Integer(Kind=wi)   :: io
    Integer(Kind=wi)   :: i
    Logical            :: endblock

    Character(Len=256) :: message
    Character(Len=265) :: set_error

    set_error = '***ERROR in &pseudo_potentials (inside &dft_settings):'

    i=1
    Do While (i<= simulation_data%total_tags)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&pseudo_potentials (inside &dft_settings)')
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&pseudo_potentials (inside &dft_settings)')
        Read (iunit, Fmt=*, iostat=io) simulation_data%dft%pseudo_pot(i)%tag, simulation_data%dft%pseudo_pot(i)%file_name
        i=i+1
      End If
    End Do

    If (i-1 < simulation_data%total_tags) Then
      Write (message,'(2a)') Trim(set_error), ' The number of declared atomic tags is less than those defined in&
                          & "&block_input_composition". Please define all tags with the corresponding file name for the&
                          & pseudopotential'
      Call error_stop(message)
    End If

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&pseudo_potentials (inside &dft_settings)')
      Call capital_to_lower_case(word)
      If (word /= '&end_pseudo_potentials') Then
        If (word(1:1) /= '#') Then
          Write (message,'(2a)') Trim(set_error), ' It seems the user has provided wrong or additional information.&
                               & This block must be closed with sentence "&end_pseudo_potentials"'
          Call error_stop(message)
        End If
      Else
          endblock=.True.
      End If
    End Do

  End Subroutine read_pseudo_poptentials

  Subroutine read_dft_magnetization(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the initial magnetization for the atomic tags 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Type(simul_type),   Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: io, i, j, k
    Character(Len=256)  :: messages(8), word
    Character(Len=64 )  :: set_error
    Logical  :: endblock
    Logical  :: header(2), error_duplication

    set_error = '***ERROR in &magnetization (inside &dft_settings):'
    Write (messages(1),'(a)') set_error

    header=.False.
    error_duplication=.False.

    i=1

    Do While (i <= 2)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&magnetization (inside &dft_settings)')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&magnetization (inside &dft_settings)')
          Call capital_to_lower_case(word) 

          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%magnetization(j)%tag, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, simulation_data%total_tags-1
                 Do k=j+1, simulation_data%total_tags
                   If (Trim(simulation_data%dft%magnetization(j)%tag)==Trim(simulation_data%dft%magnetization(k)%tag)) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(simulation_data%dft%magnetization(j)%tag),&
                                                  & 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the components of the species must be declared,&
                                                   & each tag only once'
                     Call info(messages, 3)
                     Call error_stop(' ')
                   End If
                 End Do
              End Do  
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='values') Then
            If (.Not. header(2)) Then 
              i=i+1
              Read (iunit,Fmt=*,iostat=io) word, (simulation_data%dft%magnetization(j)%value, j=1,simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the magnetization values for the defined tags. Please check if tags&
                                            & are consistent with those defined in "&Block_Species_Components".'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),'". Valid options are "tags" and "values"&
                                         & Please refer to the manual.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')      'The correct structure must be:'
          Write (messages(3),'(1x,a)')      '&magnetization'
          Write (messages(4),'(1x,a)')      '  tags      tg1     tg2     tg3    .... tgNsp'
          Write (messages(5),'(1x,a)')      '  values    mu_tg1  mu_tg2  mu_tg3 .... mu_tgNsp'
          Write (messages(6),'(1x,a)')      '&end_magnetization'
          Write (messages(7),'(1x,a,i3,a)') 'where in this case Nsp = ', simulation_data%total_tags,&
                                          & ', which corresponds to the number of tags defined in "&Block_Species_Components".'
          Write (messages(8),'(1x,a)')      'See manual for details'
          Call info(messages, 8)
          Call error_stop(' ')
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within &magnetization'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call capital_to_lower_case(word)
      If (word /= '&end_magnetization') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Descriptors "tags" and "values" have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_magnetization.' 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

  End Subroutine read_dft_magnetization

  Subroutine read_dft_hubbard(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the hubbard specification for each atomic tags 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: io, i, j, k
    Character(Len=256)  :: messages(10), word
    Character(Len=64 )  :: set_error
    Logical  :: endblock
    Logical  :: header(4), error_duplication

    set_error = '***ERROR in &hubbard (inside &dft_settings):'
    Write (messages(1),'(a)') set_error

    header=.False.
    error_duplication=.False.

    i=1

    Do While (i <= 4)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&hubbard (inside &dft_settings)')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&hubbard (inside &dft_settings)')
          Call capital_to_lower_case(word) 

          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%hubbard(j)%tag, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, simulation_data%total_tags-1
                 Do k=j+1, simulation_data%total_tags
                   If (Trim(simulation_data%dft%hubbard(j)%tag)==Trim(simulation_data%dft%hubbard(k)%tag)) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(simulation_data%dft%hubbard(j)%tag),&
                                                    & 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the components of the species must be declared,&
                                                    & each tag only once'
                     Call info(messages, 3)
                     Call error_stop(' ')
                   End If
                 End Do
              End Do  
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='l_orbital') Then
            If (.Not. header(2)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%hubbard(j)%l_orbital, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the l_orbital (which orbital to apply the correction)&
                                           & for each tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='u') Then
            If (.Not. header(3)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%hubbard(j)%U, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the "U" values for tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(3)=.True.        
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='j') Then
            If (.Not. header(4)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%hubbard(j)%J, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the "J" values for tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(4)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),'". Valid options are: "tags", "l_orbital", "U" and "J".&
                                         & Please refer to the manual.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')    'The correct structure must be:'
          Write (messages(3),'(1x,a)')    '&hubbard'
          Write (messages(4),'(1x,a)')    '  tags       tg1     tg2     tg3    .... tgNsp'
          Write (messages(5),'(1x,a)')    '  l_orbital  l_tg1   l_tg2   l_tg3  .... l_tgNsp'
          Write (messages(6),'(1x,a)')    '  U          U_tg1   U_tg2   U_tg3  .... U_tgNsp'
          Write (messages(7),'(1x,a)')    '  J          J_tg1   J_tg2   J_tg3  .... J_tgNsp'
          Write (messages(8),'(1x,a)')    '&end_hubbard'
          Write (messages(9),'(1x,a,i3,a)') 'where in this case Nsp = ', simulation_data%total_tags, &
                                          & ', which corresponds to the number of tags defined in "&Block_Species_Components".'
          Write (messages(10),'(1x,a)')    'See manual for details'
          Call info(messages, 10)
          Call error_stop(' ')
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within &hubbard'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&hubbard (inside &dft_settings)')
      Call capital_to_lower_case(word)
      If (word /= '&end_hubbard') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Descriptors "tags", "l_orbital", "U" and "J" have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_hubbard' 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

  End Subroutine read_dft_hubbard


  Subroutine read_dft_ngwf(iunit, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the NGWF specification for each atomic tags (only for ONETEP)
    !
    ! author    - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(simul_type),  Intent(InOut) :: simulation_data

    Integer(Kind=wi) :: io, i, j, k
    Character(Len=256)  :: messages(9), word
    Character(Len=64 )  :: set_error
    Logical  :: endblock
    Logical  :: header(4), error_duplication

    set_error = '***ERROR in &ngwf (inside &dft_settings):'
    Write (messages(1),'(a)') set_error

    header=.False.
    error_duplication=.False.

    i=1

    Do While (i <= 3)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&ngwf')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&ngwf (inside &dft_settings)')
          Call capital_to_lower_case(word) 
          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%ngwf(j)%tag, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              Do j=1, simulation_data%total_tags-1
                 Do k=j+1, simulation_data%total_tags
                   If (Trim(simulation_data%dft%ngwf(j)%tag)==Trim(simulation_data%dft%ngwf(k)%tag)) Then
                     Write (messages(2),'(3(1x,a))') 'Tag', Trim(simulation_data%dft%ngwf(j)%tag),&
                                                    & 'is repeated in the list!'
                     Write (messages(3),'((1x,a))') 'All tags for the components of the species must be declared,&
                                                    & each tag only once'
                     Call info(messages, 3)
                     Call error_stop(' ')
                   End If
                 End Do
              End Do  
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='number') Then
            If (.Not. header(2)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%ngwf(j)%ni, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the number of ngwf functions (row labelled as number)&
                                           & for each tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(2)=.True.        
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='radius') Then
            If (.Not. header(3)) Then 
              i=i+1
              Read (iunit, Fmt=*, iostat=io) word, (simulation_data%dft%ngwf(j)%radius, j = 1, simulation_data%total_tags)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the valus of "radius" for the tags. Please check'
                Call info(messages, 2)
                Call error_stop(' ')
              End If  
              header(3)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word),&
                                        & '". Valid options are: "tags", "number" and "radius". Please refer to the manual.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')    'The correct structure must be:'
          Write (messages(3),'(1x,a)')    '&ngwf'
          Write (messages(4),'(1x,a)')    '  tags         tg1     tg2     tg3    .... tgNsp'
          Write (messages(5),'(1x,a)')    '  number       n_tg1   n_tg2   n_tg3  .... n_tgNsp'
          Write (messages(6),'(1x,a)')    '  radius       r_tg1   r_tg2   r_tg3  .... r_tgNsp'
          Write (messages(7),'(1x,a)')    '&end_ngwf'
          Write (messages(8),'(1x,a,i3,a)') 'where in this case Nsp = ', simulation_data%total_tags, &
                                        & ', which corresponds to the number of tags defined in "&Block_Species_Components".'
          Write (messages(9),'(1x,a)')    'See manual for details.'
          Call info(messages, 9)
          Call error_stop(' ')
        End If
      End If
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within &ngwf'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word

      Call capital_to_lower_case(word)
      If (word /= '&end_ngwf') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'Descriptors "tags", "number" and "radius" have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be closed with sentence &end_ngwf' 
            Call info(messages,2)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

  End Subroutine read_dft_ngwf

End module settings

