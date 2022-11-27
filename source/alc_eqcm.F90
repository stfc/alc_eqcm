!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! Welcome to ALC_EQCM
!
! This program performs post-processing of EQCM data to the 
! purpose of electrochemical characterization and building
! of atomistic models compatible with experiments. This code
! also generates input files for atomistic simulations and 
! scripts for computation in High Performance Computing 
! facilities.
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)  
!               
!
! Author:     Ivan Scivetti (i.scivetti)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program alc_eqcm 

  Use atomistic_models, Only: model_type, &
                              build_atomistic_models, &
                              generate_hpc_simulation_files
  
  Use electrode,         Only: electrode_type
  
  Use eqcm,              Only: eqcm_cycling_summary, &
                               eqcm_check_voltage_range, &
                               eqcm_filter,     & 
                               eqcm_mass_calibration,    &
                               eqcm_massogram,  &
                               eqcm_spectra,    &
                               eqcm_type,       &
                               get_eqcm_charge, & 
                               print_eqcm_data, &
                               read_eqcm_data 
  Use fileset,           Only: file_type, &
                               FILE_HPC_SETTINGS, &
                               FILE_OUT_EQCM, &
                               FILE_RECORD_MODELS,&
                               FOLDER_RESTART,& 
                               NUM_FILES, &
                               print_header_out_eqcm, &
                               set_system_files, &
                               wrapping_up
  Use fft,               Only: fft_type, & 
                               fftnmax 
  Use filtering,         Only: filter_type
  Use hpc,               Only: hpc_type, &
                               build_hpc_script 
  Use numprec,           Only: wi,& 
                               wp
  Use redox,             Only: redox_type, &
                               redox_characterization     
  Use sauerbrey,         Only: sauer_type, &
                               sauerbrey_correlation, &
                               sauerbrey_transformation
  Use simulation_setup,  Only: simul_type 
  Use settings,          Only: read_settings, &
                               check_all_settings,&
                               check_feasibility
  Use stoichiometry,     Only: stoich_type, &
                               stoichometry_analysis
  Use system,            Only: system_type
  Use unit_output,       Only: info

Implicit None

! all simulation variables
  Type(eqcm_type)      :: eqcm_data
  Type(file_type)      :: files(NUM_FILES)
  Type(filter_type)    :: filter_data 
  Type(fft_type)       :: fft_data
  Type(redox_type)     :: redox_data
  Type(sauer_type)     :: sauer_data 
  Type(system_type)    :: system_data
  Type(electrode_type) :: electrode_data
  Type(stoich_type)    :: stoich_data
  Type(model_type)     :: model_data
  Type(simul_type)     :: simulation_data
  Type(hpc_type)       :: hpc_data

  !Time related variables
  Character(Len=256) :: message
  Integer(kind=wi)   :: start,finish,rate

  Call system_clock(count_rate=rate)
  ! Record initial time
  Call system_clock(start)
  ! Initialise settings for input/output files
  Call set_system_files(files)
  ! Print header of OUT_EQCM
  Call print_header_out_eqcm(files) 
  ! Read settings from SET_EQCM
  Call read_settings(files, eqcm_data, fft_data, filter_data, system_data, electrode_data, & 
                   & sauer_data, stoich_data, model_data, simulation_data, hpc_data)
! Check the specification of settings in SET_EQCM
  Call check_all_settings(files, eqcm_data, fft_data, filter_data, system_data, electrode_data, &
                   & sauer_data, stoich_data, model_data, simulation_data, hpc_data)

  ! HPC related directives
  If (hpc_data%generate) Then
    Call build_hpc_script(files, model_data%output_model_format%type, hpc_data)
  End If

  If (Trim(eqcm_data%analysis%type) /= 'hpc_simulation_files') Then
    If (eqcm_data%analysis%type     /= 'model_pristine_sample' .And. &
        eqcm_data%analysis%type     /= 'model_disordered_system') Then 
      !Read settings from DATA_EQCM
      Call read_eqcm_data(files, eqcm_data, fftnmax)
      !Check feasibility of the requested simulation from EQCM settings and EQCM data 
      Call check_feasibility(files, eqcm_data, filter_data, system_data, electrode_data, sauer_data)

      If (Trim(eqcm_data%analysis%type) /= 'spectra' .And. &
          Trim(eqcm_data%analysis%type) /= 'print_eqcm_raw'    .And. &
          Trim(eqcm_data%analysis%type) /= 'print_eqcm_filter') Then
         !If current is reported, compute the accummulated charge
         If (eqcm_data%current%fread) Then
           Call get_eqcm_charge(eqcm_data) 
         End If 
      End If 
    
      !Check voltage range and report summary of EQCM cycles
      If (eqcm_data%voltage_range%fread) Then
        Call eqcm_check_voltage_range(eqcm_data)
      End If
        
      ! Summarise cycling
      Call eqcm_cycling_summary(files, eqcm_data)
       
      !Print mass-frequency and accumulated charge for mass calibration
      If (eqcm_data%analysis%type=='mass_calibration') Then
        Call eqcm_mass_calibration(eqcm_data, files) 
      Else
      !Perform mass-frequency to mass conversion   
        If (eqcm_data%mass_frequency%fread .And. (.Not. eqcm_data%mass%fread)) Then
          Call sauerbrey_correlation(files, eqcm_data, system_data, sauer_data) 
          Call sauerbrey_transformation(eqcm_data, sauer_data)
        End If
      !Obtain spectrum in the recirocal space if "spectra" type of analysis is chosen 
        If (eqcm_data%analysis%type=='spectra') Then 
          Call eqcm_spectra(eqcm_data, fft_data, files)
        Else
          !Smooth EQCM data if filter is active
          If (filter_data%cutoff%fread)then
            Call eqcm_filter(eqcm_data, fft_data, filter_data)
          End If
      
          !Print EQCM data either raw or filtered, according 
          !to the definition in the  SET_EQCM file
          If (eqcm_data%analysis%type=='print_eqcm_raw' .Or. eqcm_data%analysis%type=='print_eqcm_filter') Then
            Call print_eqcm_data(eqcm_data, fft_data, files)
          Else If (eqcm_data%analysis%type=='massogram') Then
            Call eqcm_massogram(eqcm_data, files)
          Else
            !Analyse collected data 
            Call redox_characterization(files, redox_data, eqcm_data, electrode_data)
            If (eqcm_data%analysis%type =='stoichiometry' .Or. eqcm_data%analysis%type=='model_cycled_sample') Then
              !Call stoichiometry     
              Call stoichometry_analysis(files, eqcm_data, electrode_data, redox_data, stoich_data)
            End If
          End If
        End If
      End If
    End If
    ! Build atomistic_models. Option "model_cycled_sample" needs all the previous analysis
    ! If &block_simulation_settings is defined, ALC_EQCM generates input files for simulations 
    ! If &block_hpc_settings is defined, ALC_EQCM generates script files for job submission
    If (eqcm_data%analysis%type == 'model_pristine_sample'   .Or. &
        eqcm_data%analysis%type == 'model_disordered_system' .Or. &
        eqcm_data%analysis%type=='model_cycled_sample') Then
      Call build_atomistic_models(files, electrode_data, eqcm_data, redox_data,&
                                & stoich_data, model_data, simulation_data, hpc_data) 
    End If
  Else
    ! Building simulations files and hpc settings
    ! No EQCM data analysis nor modelling.
     Call generate_hpc_simulation_files(model_data%output_model_format%type, files, stoich_data,&
                                      & model_data, simulation_data, hpc_data)
  End If

  ! Record final time
  Call system_clock(finish)

  Call info(' ', 1)
  Call info(' ==========================================', 1)
  Write (message, '(1x,a,f9.3,a)') 'Total execution time = ',  Real(finish-start,Kind=wp)/rate,  ' seconds.' 
  Call info(message, 1)
  Call info(' ==========================================', 1)

  ! Print appendix to OUT_EQCM file
  Call wrapping_up(files)

  ! Create restart for subsequent runs, only if atomic models have been generated
  If (Trim(eqcm_data%analysis%type) == 'model_cycled_sample' .Or. &
     Trim(eqcm_data%analysis%type) == 'model_disordered_system' .Or. &
     Trim(eqcm_data%analysis%type) == 'model_pristine_sample') Then
     Call execute_command_line('[ ! -d '//Trim(FOLDER_RESTART)//' ] && '//'mkdir '//Trim(FOLDER_RESTART))
     Call execute_command_line('cp '//Trim(files(FILE_OUT_EQCM)%filename)//' '//Trim(FOLDER_RESTART)//'/OUT_EQCM_BACKUP') 
     Call execute_command_line('mv '//Trim(files(FILE_RECORD_MODELS)%filename)//' '//Trim(FOLDER_RESTART)) 
  End If
  
  ! Remove temporary HPC_SETTINGS file
  If (hpc_data%generate) Then
    Call execute_command_line('rm '//Trim(files(FILE_HPC_SETTINGS)%filename))
  End If

End Program alc_eqcm
