!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to check, define and print DFT directives for simulations
! with CASTEP. This module also warns the user about aspects to take
! into consideration when performing simulations
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author    - i.scivetti Oct 2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module dft_castep
  
  Use constants,        Only : date_RELEASE, &
                               chemsymbol, & 
                               max_components, &
                               K_to_eV, &
                               NPTE

  Use fileset,          Only : file_type, &
                               FILE_SET_EQCM, &
                               FILE_SET_SIMULATION,&
                               FOLDER_DFT 
  Use numprec,          Only : wi, &
                               wp
  Use process_data,     Only : capital_to_lower_case
  Use references,       Only : bib_blyp, bib_g06, bib_jchs, bib_mbd, bib_obs, bib_pbe, bib_pbesol, bib_pw91, bib_pz, &
                               bib_rp, bib_ts, bib_wc, bib_blyp, bib_revpbe, web_castep
  Use simulation_tools, Only : simul_type,&
                               check_extra_directives, &
                               check_initial_magnetization, &
                               print_warnings, &
                               record_directive, &
                               scan_extra_directive 
  Use stoichiometry,    Only : stoich_type
  Use unit_output,      Only : error_stop,&
                               info 

  Implicit None
  Private

  Public :: define_castep_settings, print_castep_settings, warnings_castep
  
Contains

  Subroutine define_castep_settings(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define (and check) settings for CASTEP directives (set defaults values)
    !
    ! author    - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data   

    Character(Len=256)  :: message, messages(15)
    Character(Len=256)  :: pp_path, pp_file, pp_name
    Character(Len=256)  :: error_dft, error_motion
    Character(Len=256)  :: ref_pp_extension, pp_extension 
    Integer(Kind=wi)    :: i, k, internal
    Logical             :: error, loop

    error_dft    = '***ERROR in &dft_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
    error_motion = '***ERROR in &motion_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
    pp_path   = Trim(FOLDER_DFT)//'/PPs/'    

    ! latest version of the code
    simulation_data%code_version= '18.1' 

    ! Pseudopotentials
    !!!!!!!!!!!!!!!!!!
    ! Check that atomic tags do not have more than 3 characters (only if &pseudo_potentials block is defined)
    If (simulation_data%dft%pp_info%fread)Then 
      Do i=1, simulation_data%total_tags
        If(Len_Trim(simulation_data%component(i)%tag)>3)Then
          Write (messages(1),'(1x,a)') '***ERROR: In CASTEP, tags for atomic species must not contain&
                                      & more than 3 characters. Please rename atomic tags' 
        End If
      End Do

      ! Check if all PP files correspond to the same format
      Do i =1, simulation_data%total_tags
        ! check the extension of the file (recpot vs. abinit vs. unknown)
        Write (pp_name, '(a)')  Trim(simulation_data%dft%pseudo_pot(i)%file_name)
        If (Index(pp_name,'.usp') /= 0 .Or. Index(pp_name,'.recpot') /= 0) Then

          If (Index(pp_name,'.recpot') /= 0) Then
            pp_extension='.recpot'
            If (i==1) Then
              ref_pp_extension='.recpot'
            End If         

          Else If (Index(pp_name,'.usp') /= 0) Then
            pp_extension='.usp'
            If (i==1) Then
              ref_pp_extension='.usp'
            End If

          End If

          If (Trim(pp_extension) /= Trim(ref_pp_extension)) Then
             Write (message, '(1x,3a)') '***ERROR: pseudo potential files defined in &pseudo_potentials contain&
                                       & different extensions. For ONETEP, all files must be either ".recpot"&
                                       & or ".usp" type.'
             Call error_stop(message)
          End If

        Else
          Write (message,'(1x,5a)') '***ERROR: The extension of file ',  Trim(pp_name), ' for the pseudo potential of species "',&
                                  & Trim(simulation_data%dft%pseudo_pot(i)%tag), '" is not recognised by CASTEP.&
                                  & Valid extensions are ".recpot" and ".usp".'
          Call error_stop(message)
        End If 

      End Do 

      ! Check consistency between pseudpotentials and XC directive
      Do i=1, simulation_data%total_tags
        pp_file=Trim(pp_path)//Trim(simulation_data%dft%pseudo_pot(i)%file_name)
        ! Open PP file
        Open(Newunit=internal, File=Trim(pp_file), Status='old')
        If (Trim(ref_pp_extension) =='.recpot') Then
          Call check_recpot_castep(internal, simulation_data, i)  
        Else If (Trim(ref_pp_extension) == '.usp') Then
          Call check_usp_castep(internal, simulation_data, i) 
        End If
        ! Close PP file
        Close(internal)
      End Do
    End If

    ! Electronic structure
    !!!!!!!!!!!!!!!!!!!!!!

    ! Check XC_version
    If (Trim(simulation_data%dft%xc_version%type) /= 'pz'    .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'pw91'   .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'pbe'    .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'revpbe' .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'pbesol' .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'blyp'   .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'wc' ) Then
      Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                  &'Invalid specification for directive "XC_version" for CASTEP.&
                                  & Implement options for CASTEP are:'
      Write (messages(2),'(1x,a)')   '==== LDA-level =================='
      Write (messages(3),'(1x,2a)')  '- PZ      Perdew-Zunger          ', Trim(bib_pz) 
      Write (messages(4),'(1x,a)')   '==== GGA-level =================='
      Write (messages(5),'(1x,2a)')  '- PW91    Perdew-Wang 91         ', Trim(bib_pw91)
      Write (messages(6),'(1x,2a)')  '- PBE     Perdew-Burke-Ernzerhof ', Trim(bib_pbe)
      Write (messages(7),'(1x,2a)')  '- revPBE  revPBE                 ', Trim(bib_revpbe)
      Write (messages(8),'(1x,2a)')  '- PBEsol  PBE for solids         ', Trim(bib_pbesol)
      Write (messages(9),'(1x,2a)')  '- WC      Wu-Cohen               ', Trim(bib_wc)
      Write (messages(10),'(1x,2a)') '- BLYP    Becke-Lee-Young-Parr   ', Trim(bib_blyp)
      Write (messages(11),'(1x,a)')  '================================='
      Call info(messages, 11)
      Call error_stop(' ')
    End If

    ! Check if basis set was defined, complain and abort
    If (simulation_data%dft%basis_info%fread) Then
      Write (message,'(1x,2a)') Trim(error_dft), &
                        &' Definition of basis sets is not required for CASTEP.&
                        & Please remove &block_basis_set and rerun.'
      Call error_stop(message)
    End If

    ! vdW settings
    !!!!!!!!!!!!!!!!!!!!!!
    simulation_data%dft%need_vdw_kernel=.False.
    If (simulation_data%dft%vdw%fread) Then
      If (Trim(simulation_data%dft%vdw%type) /= 'g06'  .And.&
         Trim(simulation_data%dft%vdw%type)  /= 'obs'  .And.&
         Trim(simulation_data%dft%vdw%type)  /= 'jchs' .And.&
         Trim(simulation_data%dft%vdw%type)  /= 'ts'   .And.&
         Trim(simulation_data%dft%vdw%type)  /= 'mbd'   ) Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                  &'Invalid specification of directive "vdW" for CASTEP. Valid options are:'
        Write (messages(2),'(1x,2a)')  '- G06    Grimme 2006 ', Trim(bib_g06)
        Write (messages(3),'(1x,2a)')  '- OBS    Ortmann-Bechstedt-Schmidt ', Trim(bib_obs)
        Write (messages(4),'(1x,2a)')  '- JCHS   Jurecka-Cerny-Hobza-Salahub ', Trim(bib_jchs)
        Write (messages(5),'(1x,2a)')  '- TS     Tkatchenko-Scheffler method ', Trim(bib_ts)
        Write (messages(6),'(1x,2a)')  '- MBD    Many-body dispersion energy method ', Trim(bib_mbd)
        Call info(messages, 6)
        Call error_stop(' ')
      End If

      ! vdW only for GGA type
      If (Trim(simulation_data%dft%xc_level%type) /= 'gga') Then
        Write (message,'(1x,4a)') Trim(error_dft), &
                                &' Dispersion correction type "', Trim(simulation_data%dft%vdw%type), '" requires of GGA&
                                & option for directive XC_level. Please change'
        Call error_stop(message)
      End If

      If (Trim(simulation_data%dft%vdw%type) /= 'obs') Then
         If (Trim(simulation_data%dft%xc_version%type) /= 'pbe') Then
           Write (message,'(1x,4a)') Trim(error_dft), &
                                  &' Dispersion correction type "', Trim(simulation_data%dft%vdw%type), '" requires of&
                                  & option PBE for XC_version.'
           Call error_stop(message)
         End If
      Else
        If (Trim(simulation_data%dft%xc_version%type) /= 'pw91') Then 
           Write (message,'(1x,4a)') Trim(error_dft), &
                                  &' Dispersion correction type "', Trim(simulation_data%dft%vdw%type), '" requires of&
                                  & option PW91 for XC_version.'

        Call error_stop(message)
        End If
      End If
    End If

    ! Orbital transformation
    If (simulation_data%dft%ot%stat) Then
      Write (message,'(2(1x,a))') Trim(error_dft), 'Requested Orbital Transformation via directive "OT"&
                                 & is not possible for CASTEP simulations. Please remove it'
      Call error_stop(message)
    End If

   ! Energy cutoff 
    If (Trim(simulation_data%dft%encut%units)/='ev') Then
       Write (message,'(2(1x,a))') Trim(error_dft), &
                                   &'Units for directive "energy_cutoff" for CASTEP simulations must be in eV'
       Call error_stop(message)
    End If
      
   ! precision
    If (simulation_data%dft%precision%fread) Then
      Write (message,'(2(1x,a))') Trim(error_dft), 'For CASTEP, "precision" directive is not needed, as it is&
                                 & determined by the energy cutoff. Please remove "precision" and rerun.'
      Call error_stop(message)
    End If

   ! Smearing
   If (simulation_data%dft%smear%fread) Then
     If (Trim(simulation_data%dft%smear%type) /= 'gaussian' .And.&
        Trim(simulation_data%dft%smear%type) /= 'fermi'     .And.& 
        Trim(simulation_data%dft%smear%type) /= 'fix_occupancy') Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                     &'Invalid specification of directive "smearing" for CASTEP. Options are:'
        Write (messages(2),'(1x,a)') '- Gaussian       (Gaussian distribution)' 
        Write (messages(3),'(1x,a)') '- Fermi          (Fermi-Dirac distribution)'
        Write (messages(4),'(1x,a)') '- Fix-Occupancy  (Fix the occupancies of the bands)' 
        Write (messages(5),'(1x,a)') ' '
        Write (messages(6),'(1x,a)') 'IMPORTANT: options GaussianSplines, HermitePolynomials and ColdSmearing&
                                     & are not implemented in ALC_EQCM.' 
        Call info(messages, 6)
        Call error_stop(' ')
     End If
   Else
     Write (message,'(2(1x,a))') Trim(error_dft), 'The user must specify directive "smearing" for CASTEP simulations'
     Call error_stop(message)
   End If 
 
   ! Mixing 
   If (simulation_data%dft%mixing%fread) Then
     If (simulation_data%dft%edft%stat) Then 
       Write (message,'(2(1x,a))') Trim(error_dft), 'For CASTEP, the definition of "mixing_scheme" is incompatible with&
                                                & eDFT. Please review settings.'
       Call error_stop(message)
     End If        
     If (Trim(simulation_data%dft%mixing%type)   /= 'kerker'      .And.&
        Trim(simulation_data%dft%mixing%type)    /= 'linear'      .And.& 
        Trim(simulation_data%dft%mixing%type)    /= 'broyden-2nd' .And.& 
        Trim(simulation_data%dft%mixing%type)    /= 'pulay') Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                     &'Invalid specification of directive "mixing_scheme" for CASTEP. Options are:'
        Write (messages(2),'(1x,a)') '- Kerker' 
        Write (messages(3),'(1x,a)') '- Linear'
        Write (messages(4),'(1x,a)') '- Broyden-2nd (input option is "Broyden")' 
        Write (messages(5),'(1x,a)') '- Pulay' 
        Call info(messages, 5)
        Call error_stop(' ')
     End If
   Else
     simulation_data%dft%mixing%type='broyden-2nd'      
   End If 

   ! Fix occupancy
   If (Trim(simulation_data%dft%smear%type) == 'fix_occupancy') Then
     If (simulation_data%dft%edft%stat) Then 
       Write (message,'(2(1x,a))') Trim(error_dft), 'For CASTEP, the choice of "fix_occupancy" for the "smearing" directive&
                                  & is incompatible with the request of EDFT. Please review settings'
       Call error_stop(message)
     End If        
           
     If (simulation_data%dft%bands%fread) Then
       Write (message,'(2(1x,a))') Trim(error_dft), 'For CASTEP, the definition of "bands" is incompatible with&
                                 & the choice of "fix_occupancy" for the "smearing" directive. Please review settings'
       Call error_stop(message) 
     End If
     If (simulation_data%dft%width_smear%fread) Then
       Write (message,'(2(1x,a))') Trim(error_dft), 'For CASTEP, the definition of "width_smear" is incompatible with&
                                 & the choice of "fix_occupancy" for the "smearing" directive. Please review settings'
       Call error_stop(message) 
     End If
     If (simulation_data%dft%mixing%fread) Then
       Write (message,'(2(1x,a))') Trim(error_dft), 'Definition of "mixing_scheme" is incompatible with&
                                 & the choice of "fix_occupancy" for the "smearing" directive. Please review settings'
       Call error_stop(message) 
     End If
   End If
        
    ! Width smearing
    If (Trim(simulation_data%dft%smear%type) /= 'fix_occupancy') Then
      If (.Not. simulation_data%dft%bands%fread) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Requested settings for computation of the electronic problem&
                                  & needs the specification of extra number of bands via "bands".'
        Call error_stop(message) 
      End If
      If (.Not. simulation_data%dft%width_smear%fread) Then
        simulation_data%dft%width_smear%value=0.20_wp
        simulation_data%dft%width_smear%units='eV'
      Else
        If (Trim(simulation_data%dft%width_smear%units) /= 'ev') Then
           Write (message,'(2a)')  Trim(error_dft), ' Units of directive "width_smear" for CASTEP must be in eV'
          Call error_stop(message)
        End If
      End If
    End If
   
    ! SCF energy tolerance 
    If (simulation_data%dft%delta_e%fread) Then
      If (Trim(simulation_data%dft%delta_e%units) /= 'ev' ) Then
         Write (message,'(2a)')  Trim(error_dft), ' Units for directive "SCF_energy_tolerance" in CASTEP must be eV'
         Call info(message, 1)
         Call error_stop(' ')
      End If
    End If

    ! kpoints
    If (Trim(simulation_data%dft%kpoints%tag) == 'automatic' ) Then
       Write (messages(1),'(2a)')  Trim(error_dft), 'The only option available in CASTEP for the kpoint mesh is&
                                  & "MPack". Please change and rerun.'
      Call info(messages, 1)
      Call error_stop(' ')
    End If

    ! Ions related settings
    !!!!!!!!!!!!!!!!!!!!!!!
    !Relaxation method
    If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      If (Trim(simulation_data%motion%relax_method%type) /= 'bfgs'    .And. &
        Trim(simulation_data%motion%relax_method%type)  /= 'lbfgs'  ) Then
        Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                &'Invalid specification of directive "relax_method" for CASTEP. Implemented options are:'
        Write (messages(2),'(1x,a)') '- BFGS  (Broyden-Fletcher-Goldfarb-Shanno)'
        Write (messages(3),'(1x,a)') '- LBFGS (Linear Broyden-Fletcher-Goldfarb-Shanno)'
        Write (messages(4),'(1x,a)') ' '
        Write (messages(5),'(1x,a)') 'WARNING: CASTEP options Delocalized, DampedMD and TPSD are not implemented in ALC_EQCM.'
        Call info(messages, 5)
        Call error_stop(' ')
      End If
    
      If (simulation_data%motion%ion_steps%fread) Then
        If (simulation_data%motion%ion_steps%value == 1) Then
          simulation_data%simulation%type='singlepoint'   
          Call info(' ***WARNING: since the number of ionic steps was set to 1, the simulation was changed to SinglePoint', 1)
        Else If (simulation_data%motion%ion_steps%value == 2) Then
          Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                &' In CASTEP, "ion_steps" for geometry relaxation must be larger than 2. Please change'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If 
     
    End If

    ! Force tolerance
    error=.False.
    If (simulation_data%motion%delta_f%fread) Then
      If (Trim(simulation_data%motion%delta_f%units(1)) /= 'ev' ) Then
        error=.True.
      End If
      If (Trim(simulation_data%motion%delta_f%units(2)) /= 'angstrom-1' ) Then
        error=.True.
      End If
    Else
      simulation_data%motion%delta_f%units(1)='eV' 
      simulation_data%motion%delta_f%units(2)='Angstrom-1'
      simulation_data%motion%delta_f%value(1)= 0.01_wp 
    End If
 
    If (error) Then
      Write (messages(1),'(2a)')  Trim(error_motion), 'Invalid units of directive "force_tolerance" for CASTEP.&
                                & Units must be "eV Angstrom-1"'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
 
    ! Timestep
    If (.Not. simulation_data%motion%timestep%fread) Then
      simulation_data%motion%timestep%units='fs'
      simulation_data%motion%timestep%value= 1.0_wp 
    End If

    ! Ensemble 
    If (simulation_data%motion%ensemble%fread) Then
        If (Trim(simulation_data%motion%ensemble%type) /= 'nve'  .And.&
           Trim(simulation_data%motion%ensemble%type) /= 'nvt'  .And.&
           Trim(simulation_data%motion%ensemble%type) /= 'npt'  .And.&
           Trim(simulation_data%motion%ensemble%type) /= 'nph'  ) Then
           Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                    &'Invalid specification of "ensemble" for CASTEP. Options are:'
           Write (messages(2),'(1x,a)') '- NVE (Microcanonical ensemble)'
           Write (messages(3),'(1x,a)') '- NVT (Canonical ensemble)'
           Write (messages(4),'(1x,a)') '- NpT (Isothermal-Isobaric ensemble)'
           Write (messages(5),'(1x,a)') '- NpH (Isoenthalpic-Isobaric ensemble)'
           Call info(messages, 5)
           Call error_stop(' ')
        End If
    End If
   
    ! Thermostat
    If (Trim(simulation_data%simulation%type) == 'md') Then    
      If (simulation_data%motion%thermostat%fread) Then
        If (Trim(simulation_data%motion%thermostat%type) /= 'langevin'  .And. &
           Trim(simulation_data%motion%thermostat%type) /= 'nose-hoover'  ) Then
          Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                  &'Specification for "thermostat" is not supported by CASTEP. Options are:'
          Write (messages(2),'(1x,a)') '- Langevin'
          Write (messages(3),'(1x,a)') '- Nose-Hoover'
          Call info(messages, 3)
          Call error_stop(' ')
        End If
      End If
    End If

    ! Barostat
    If (Trim(simulation_data%simulation%type) == 'md') Then
      If (Trim(simulation_data%motion%ensemble%type) == 'npt'  .Or.&
          Trim(simulation_data%motion%ensemble%type) == 'nph'  ) Then
        If (Trim(simulation_data%motion%barostat%type)/='andersen-hoover' .And. &
           Trim(simulation_data%motion%barostat%type)/='parrinello-rahman')Then
           Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                            &'Specification for "barostat" is not supported by CASTEP. Options are:'
           Write (messages(2),'(1x,a)') '- Andersen-Hoover'
           Write (messages(3),'(1x,a)') '- Parrinello-Rahman'
           Call info(messages, 3)
           Call error_stop(' ')
        End If
      End If
    End If

    ! Relaxation time for the thermostat
    If (Trim(simulation_data%motion%ensemble%type) == 'nvt' .Or. &
       Trim(simulation_data%motion%ensemble%type) == 'npt' ) Then
      If (.Not. simulation_data%motion%relax_time_thermostat%fread) Then
          Write (messages(1),'(1x,4a)') Trim(error_motion), ' For CASTEP, thermostat "', &
                                       Trim(simulation_data%motion%thermostat%type), '" requires the specification&
                                       & of "relax_time_thermostat", which is missing.'
          Call info(messages, 1)
          Call error_stop(' ')
      End If
    End If

    ! Relaxation time for the barostat
    If (Trim(simulation_data%motion%ensemble%type) == 'npt' .Or. &
       Trim(simulation_data%motion%ensemble%type) == 'nph') Then 
      If (.Not. simulation_data%motion%relax_time_barostat%fread) Then
        Write (messages(1),'(1x,4a)') Trim(error_motion), ' For "', &
                                     Trim(simulation_data%motion%ensemble%type), '" simulations in CASTEP, the&
                                     & user must specify "relax_time_barostat", which is missing.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If
      
    If (simulation_data%extra_info%stat) Then
      ! Check if user defined directives contain only symbol "="
      Do i = 1, simulation_data%extra_directives%N0
        Call check_extra_directives(simulation_data%extra_directives%array(i), &
                                    simulation_data%extra_directives%key(i),   &
                                    simulation_data%extra_directives%set(i), ':', 'CASTEP')
      End Do
    End If

    If (simulation_data%motion%mass_info%fread) Then
      ! Assing atomic numbers
      Do i=1, simulation_data%total_tags
        loop=.True.
        k=1
        Do While (k <= NPTE .And. loop)
          If (Trim(chemsymbol(k))==Trim(simulation_data%component(i)%tag)) Then
            loop=.False.
          End If
          k=k+1
        End Do
        If (loop) Then
          If (Index(Trim(simulation_data%component(i)%tag),':') == 0) Then
             Write (messages(1),'(1x,3a)') '***ERROR: Definition of atomic tag "', &
                                   & Trim(simulation_data%component(i)%tag), &
                                   & '" must include the symbol ":" for proper execution in CASTEP'
             Write (messages(2),'(1x,a)')  'To the left of ":" the user must set the chemical symbol&
                                         & and to the right any specification.'
             Write (messages(3),'(1x,a)') 'For example, deuterium can be defined as "H:D". Bear in mind&
                                         & one cannot exceeds 4 characters'
             Call info(messages, 3)  
             Call error_stop(' ')            
          End If  
        End If
      End Do
 
    End If

  End Subroutine define_castep_settings

  Subroutine print_castep_settings(files, net_elements, list_element, list_tag, list_number, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define settings for CASTEP directives
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Integer(Kind=wi),   Intent(In   ) :: net_elements 
    Character(Len=2),   Intent(In   ) :: list_element(max_components) 
    Character(Len=8),   Intent(In   ) :: list_tag(max_components) 
    Integer(Kind=wi),   Intent(In   ) :: list_number(max_components)
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256), Dimension(max_components)  ::  l_orbital
    Integer(Kind=wi)   ::  iunit, i, j, ic
    Character(Len=256) :: pseudo_list(max_components)
    Character(Len=256) :: exec_mv
    Character(Len=256) :: pp_path
    Character(Len=256) :: message, messages(2)
    Real(Kind=wp)      :: mag_ini(max_components)
    Real(Kind=wp), Dimension(max_components) :: uc, jc  

    Logical :: found
    Logical :: loop
 
    ic=1
    pp_path   = Trim(FOLDER_DFT)//'/PPs/'    

    ! Open FILE_SET_SIMULATION file
    Open(Newunit=files(FILE_SET_SIMULATION)%unit_no, File=files(FILE_SET_SIMULATION)%filename,Status='Replace')
    iunit=files(FILE_SET_SIMULATION)%unit_no

    Write (iunit,'(a)')  '##############################'
    Write (iunit,'(a)')  '# File generated with ALC_EQCM'
    Write (iunit,'(a)')  '##############################'
    Write (iunit,'(a)') ' '
 
    Write (iunit,'(a)') '##### Type of calculation'
    If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      Write (message,'(a)') 'TASK :  GeometryOptimization' 
    Else If (Trim(simulation_data%simulation%type) == 'md') Then
      Write (message,'(a)') 'TASK :  MolecularDynamics' 
    Else If (Trim(simulation_data%simulation%type) == 'singlepoint') Then
      Write (message,'(a)') 'TASK :  SinglePoint' 
    End If
    Call record_directive(iunit, message, 'TASK', simulation_data%set_directives%array(ic), ic)

    Write (iunit,'(a)')    '    '
    Write (iunit,'(a)') '##### Electronic parameters'
    Write (iunit,'(a)') '#=========================='
    ! Define the XC part
     If (Trim(simulation_data%dft%xc_version%type) == 'pz') Then
       Write (iunit,'(2a)') '# Perdew-Zunger (PZ) LDA functional ', Trim(bib_pz) 
       Write (message,'(a)')  'XC_FUNCTIONAL : LDA'
     Else If (Trim(simulation_data%dft%xc_version%type) == 'pw91') Then
       Write (iunit,'(2a)') '# Perdew-Wang 91 (PW91) XC functional ', Trim(bib_pw91)
       Write (message,'(a)')  'XC_FUNCTIONAL : PW91'
     Else If (Trim(simulation_data%dft%xc_version%type) == 'pbe') Then
       Write (iunit,'(2a)') '# Perdew-Burke-Ernzerhof (PBE) XC functional ', Trim(bib_pbe)
       Write (message,'(a)')  'XC_FUNCTIONAL : PBE'
     Else If (Trim(simulation_data%dft%xc_version%type) == 'revpbe') Then 
       Write (iunit,'(2a)') '# Hammer-Hansen-Norskov (RPEB) XC functional ', Trim(bib_rp)
       Write (message,'(a)')  'XC_FUNCTIONAL : RPBE'
     Else If (Trim(simulation_data%dft%xc_version%type) == 'pbesol') Then 
       Write (iunit,'(2a)') '# PBE for solids (PBEsol) XC functional ', Trim(bib_pbesol)
       Write (message,'(a)')  'XC_FUNCTIONAL : PBESOL'
     Else If (Trim(simulation_data%dft%xc_version%type) == 'wc') Then 
       Write (iunit,'(6x,2a)') '# Wu-Cohen (WC) exchange (no correlation) ', Trim(bib_wc) 
       Write (message,'(a)')  'XC_FUNCTIONAL : WC'
     Else If (Trim(simulation_data%dft%xc_version%type) == 'blyp') Then
       Write (iunit,'(2a)') '# Becke-Lee-Young-Parr (BLYP) XC functional ', Trim(bib_blyp)
       Write (message,'(a)')  'XC_FUNCTIONAL : BLYP'
     End If 
     Call record_directive(iunit, message, 'XC_FUNCTIONAL', simulation_data%set_directives%array(ic), ic)
 
    If (simulation_data%dft%vdw%fread) Then
      Write (iunit,'(a)') '#==== vdW corrections'
      Write (message,'(a)') 'SEDC_APPLY : True'
      Call record_directive(iunit, message, 'SEDC_APPLY', simulation_data%set_directives%array(ic), ic)
      If (Trim(simulation_data%dft%vdw%type)     == 'g06'   ) Then
        Write (iunit,'(2a)') '# Grimme 2006 ', Trim(bib_g06)
        Write (message,'(a)')  'SEDC_SCHEME : g06'
      Else If (Trim(simulation_data%dft%vdw%type)     == 'ts'  ) Then
        Write (iunit,'(2a)') '# Tkatchenko-Scheffler method (TS) ', Trim(bib_ts)
        Write (message,'(a)')  'SEDC_SCHEME : TS'
      Else If (Trim(simulation_data%dft%vdw%type)     == 'obs' ) Then
        Write (iunit,'(2a)') '# Ortmann-Bechstedt-Schmidt (OBS) ', Trim(bib_obs)
        Write (message,'(a)')  'SEDC_SCHEME : OBS'
      Else If (Trim(simulation_data%dft%vdw%type)     == 'jchs'   ) Then
        Write (iunit,'(2a)') '# Jurecka-Cerny-Hobza-Salahub (JCHS) ', Trim(bib_jchs)
        Write (message,'(a)')  'SEDC_SCHEME : JCHS'
      Else If (Trim(simulation_data%dft%vdw%type)     == 'mbd'   ) Then
        Write (iunit,'(2a)') '# Many-body dispersion energy method MBD ', Trim(bib_mbd) 
        Write (message,'(a)')  'SEDC_SCHEME : MBD*'
      End If
      Call record_directive(iunit, message, 'SEDC_SCHEME', simulation_data%set_directives%array(ic), ic)
    End If
    
    If (Trim(simulation_data%dft%smear%type) == 'fix_occupancy') Then
      Write (message,'(a)')   'FIX_OCCUPANCY   :  True '
      Call record_directive(iunit, message, 'FIX_OCCUPANCY', simulation_data%set_directives%array(ic), ic)
    Else
      If (simulation_data%dft%edft%stat) Then
        Write (iunit,'(a)')             '#=== Ensemble DFT'      
        Write (message,'(a)')           'METALS_METHOD  :   EDFT'
        Call record_directive(iunit, message, 'METALS_METHOD', simulation_data%set_directives%array(ic), ic)
      Else
        Write (iunit,'(a)')           '#=== Density mixing'      
        Write (message,'(a)')         'METALS_METHOD  :   DM'
        Call record_directive(iunit, message, 'METALS_METHOD', simulation_data%set_directives%array(ic), ic)
        If (Trim(simulation_data%dft%mixing%type)    == 'broyden-2nd') Then
          Write (message,'(a)')  'MIXING_SCHEME  :  Broyden   # This is actually the Broyden-2nd method of Johnson.'
        Else
          Write (message,'(2a)') 'MIXING_SCHEME  :  ', Trim(simulation_data%dft%mixing%type)
        End If
        Call record_directive(iunit, message, 'MIXING_SCHEME', simulation_data%set_directives%array(ic), ic)
      End If
      If (Trim(simulation_data%dft%smear%type) == 'gaussian') Then
        Write (message,'(a)') 'SMEARING_SCHEME :    Gaussian   # Gaussian smearing'
        Call record_directive(iunit, message, 'SMEARING_SCHEME', simulation_data%set_directives%array(ic), ic)
      Else If (Trim(simulation_data%dft%smear%type) == 'fermi') Then
        Write (message,'(a)') 'SMEARING_SCHEME :    FermiDirac   # Fermi smearing'
        Call record_directive(iunit, message, 'SMEARING_SCHEME', simulation_data%set_directives%array(ic), ic)
      End If
      Write (message,'(a,f6.2,1x,a)') 'SMEARING_WIDTH : '  , simulation_data%dft%width_smear%value,&
                                     Trim(simulation_data%dft%width_smear%units)
      Call record_directive(iunit, message, 'SMEARING_WIDTH', simulation_data%set_directives%array(ic), ic)
      If (simulation_data%dft%bands%fread) Then
        Write (message,'(a,i4,a)') 'NEXTRA_BANDS : ', simulation_data%dft%bands%value,&
                                & ' # add extra bands to improve convergence'
        Call record_directive(iunit, message, 'NEXTRA_BANDS', simulation_data%set_directives%array(ic), ic)
      End If
    End If

    Write (message,'(a,i4,a)') 'MAX_SCF_CYCLES : ', simulation_data%dft%scf_steps%value, ' # maximum number of SC steps'
    Call record_directive(iunit, message, 'MAX_SCF_CYCLES', simulation_data%set_directives%array(ic), ic)
    Write (message,'(a,f8.2,1x,2a)') 'CUT_OFF_ENERGY : '  , simulation_data%dft%encut%value,&
                                     & Trim(simulation_data%dft%encut%units), '  # Energy cutoff'
    Call record_directive(iunit, message, 'CUT_OFF_ENERGY', simulation_data%set_directives%array(ic), ic)
    Write (message,'(a,e10.3,1x,2a)') 'ELEC_ENERGY_TOL : '  , simulation_data%dft%delta_e%value, &
                                   & Trim(simulation_data%dft%delta_e%units), ' # Energy tolerance'
    Call record_directive(iunit, message, 'ELEC_ENERGY_TOL', simulation_data%set_directives%array(ic), ic)

    If (simulation_data%dft%spin_polarised%stat) Then
      Write (message,'(a)')      'SPIN_POLARIZED : True  # Spin-polarised' 
      Call record_directive(iunit, message, 'SPIN_POLARIZED', simulation_data%set_directives%array(ic), ic)
    End If

   If (simulation_data%dft%mag_info%fread) Then
     If (simulation_data%dft%total_magnetization%fread) Then
       Write (iunit,'(a)') ' '
       Write (iunit,'(a)') '#==== Total magnetization'
       Write (message,'(a,f6.2)') 'SPIN : ', simulation_data%dft%total_magnetization%value
       Call record_directive(iunit, message, 'SPIN', simulation_data%set_directives%array(ic), ic)
       Do i=1, net_elements
         j=1
         loop=.True.
         Do While (j <= simulation_data%total_tags .And. loop)
           If (Trim(list_tag(i))==Trim(simulation_data%dft%magnetization(j)%tag)) Then
             mag_ini(i)=simulation_data%dft%magnetization(j)%value
             loop=.False.
           End If
           j=j+1
         End Do
       End Do
       Call check_initial_magnetization(net_elements, list_tag, list_number, mag_ini,&
                                      & simulation_data%dft%total_magnetization%value)

       If (Trim(simulation_data%dft%smear%type) /= 'fix_occupancy') Then
         Write (message,'(a)') 'SPIN_FIX      : -10' 
         Call record_directive(iunit, message, 'SPIN_FIX', simulation_data%set_directives%array(ic), ic)
         Write (message,'(a)') 'GEOM_SPIN_FIX : -10'
         Call record_directive(iunit, message, 'GEOM_SPIN_FIX', simulation_data%set_directives%array(ic), ic)
       End If

     End If
   End If

   If (simulation_data%net_charge%fread) Then
     Write (iunit,'(a)') ' '
     Write (iunit,'(a)') '#==== Net charge'
     Write (message,'(a,f20.4)') 'CHARGE : ', simulation_data%net_charge%value
     Call record_directive(iunit, message, 'CHARGE', simulation_data%set_directives%array(ic), ic)
   End If
   
   Write (iunit,'(a)') ' '
   If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
     Write (iunit,'(a)') '#### Geometry relaxation'
     Write (iunit,'(a)') '#======================='
     ! atoms      
     Write (message,'(a,i4,a)') 'GEOM_MAX_ITER : ', simulation_data%motion%ion_steps%value, ' # Number of ionic steps'
     Call record_directive(iunit, message, 'GEOM_MAX_ITER', simulation_data%set_directives%array(ic), ic)
     If (Trim(simulation_data%motion%relax_method%type) == 'bfgs') Then
       Write (message, '(a)') 'GEOM_METHOD :  BFGS   # geometry relaxation with Broyden-Fletcher-Goldfarb-Shanno method'
     Else If (Trim(simulation_data%motion%relax_method%type) == 'lbfgs') Then
       Write (message, '(a)') 'GEOM_METHOD :  LBFGS  # geometry relaxation with the linear Broyden-Fletcher-Goldfarb-Shanno method'
     End If
     Call record_directive(iunit, message, 'GEOM_METHOD', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,f6.2,a)') 'GEOM_FORCE_TOL : ', simulation_data%motion%delta_f%value(1), ' ev/ang'   
     Call record_directive(iunit, message, 'GEOM_FORCE_TOL', simulation_data%set_directives%array(ic), ic)

   Else If (Trim(simulation_data%simulation%type) == 'md') Then
     Write (iunit,'(a)') '#### Molecular dynamics'
     Write (iunit,'(a)') '#======================'
     Write (message,'(a,f6.2,a)') 'MD_DELTA_T     : ', simulation_data%motion%timestep%value,    ' fs'
     Call record_directive(iunit, message, 'MD_DELTA_T', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,f6.2,a)') 'MD_TEMPERATURE : ', simulation_data%motion%temperature%value, ' K'
     Call record_directive(iunit, message, 'MD_TEMPERATURE', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,i4,a)')   'MD_NUM_ITER    : ', simulation_data%motion%ion_steps%value, ' # Number of MD steps'
     Call record_directive(iunit, message, 'MD_NUM_ITER', simulation_data%set_directives%array(ic), ic)
     
     If (Trim(simulation_data%motion%ensemble%type) == 'nve') Then
       Write (iunit, '(a)') '#==== NVE ensemble'
       Write (message, '(a)') 'MD_ENSEMBLE : NVE'
       Call record_directive(iunit, message, 'MD_ENSEMBLE', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%motion%ensemble%type) == 'nvt' .Or. &
            & Trim(simulation_data%motion%ensemble%type) == 'npt') Then
       If (Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
         Write (iunit, '(a)') '#==== NVT ensemble'
         Write (message, '(a)') 'MD_ENSEMBLE : NVT'
       Else If (Trim(simulation_data%motion%ensemble%type) == 'npt') Then
         Write (iunit, '(a)') '#==== NPT ensemble'
         Write (message, '(a)') 'MD_ENSEMBLE : NPT'
       End If 
       Call record_directive(iunit, message, 'MD_ENSEMBLE', simulation_data%set_directives%array(ic), ic)

       Write (iunit, '(a)') '#Thermostat'
       If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover') Then
         Write (message, '(a)') 'MD_THERMOSTAT  : Nose-Hoover'
         Call record_directive(iunit, message, 'MD_THERMOSTAT', simulation_data%set_directives%array(ic), ic)
       Else If (Trim(simulation_data%motion%thermostat%type) == 'langevin') Then
         Write (message, '(a)') 'MD_THERMOSTAT  : Langevin'
         Call record_directive(iunit, message, 'MD_THERMOSTAT', simulation_data%set_directives%array(ic), ic)
       End If
       Write (message, '(a,f6.2,a)') 'MD_ION_T  : ', simulation_data%motion%relax_time_thermostat%value, ' fs'
       Call record_directive(iunit, message, 'MD_ION_T', simulation_data%set_directives%array(ic), ic)

     Else If (Trim(simulation_data%motion%ensemble%type) == 'nph') Then
        Write (iunit, '(a)') '#==== NPH ensemble'
        Write (message, '(a)') 'MD_ENSEMBLE : NPH'
        Call record_directive(iunit, message, 'MD_ENSEMBLE', simulation_data%set_directives%array(ic), ic)
     End If 

     If (Trim(simulation_data%motion%ensemble%type) == 'npt' .Or. &
           &  Trim(simulation_data%motion%ensemble%type) == 'nph') Then
        Write (iunit, '(a)') '#Barostat'
        If (Trim(simulation_data%motion%barostat%type) == 'andersen-hoover') Then
          Write (message, '(a)') 'MD_BAROSTAT : Andersen-Hoover'
        Else If (Trim(simulation_data%motion%barostat%type) == 'parrinello-rahman') Then
          Write (message, '(a)') 'MD_BAROSTAT : Parrinello-Rahman'
        End If
        Call record_directive(iunit, message, 'MD_BAROSTAT', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a,f10.2,a)') 'MD_CELL_T : ', simulation_data%motion%relax_time_barostat%value, ' fs'
        Call record_directive(iunit, message, 'MD_CELL_T', simulation_data%set_directives%array(ic), ic)
     End If
   End If

   simulation_data%set_directives%N0=ic-1

   ! Extra-directives
   If (simulation_data%extra_info%stat) Then
     Write (iunit,'(a)') ' '
     Write (iunit,'(a)') '##### Extra directives'
     Write (iunit,'(a)') '#====================='
     found=.False.
     Do i=1, simulation_data%extra_directives%N0
       Write (iunit,'(a)') Trim(Adjustl(simulation_data%extra_directives%array(i)))
       If (Index(Trim(Adjustl(simulation_data%extra_directives%array(i))), '#') /= 1 ) Then
         Call scan_extra_directive(simulation_data%extra_directives%key(i), simulation_data%set_directives, found)         
         If (found)Then
           Close(iunit)
           Call info(' ', 1)
           Write (messages(1), '(1x,3a)') '***ERROR in sub-bock &extra_directives: directive "', &
                                         & Trim(simulation_data%extra_directives%key(i)),&
                                         & '" has already been defined.'
           Write (messages(2), '(1x,3a)') 'Please check temporary file ', Trim(files(FILE_SET_SIMULATION)%filename), &
                                     & '. The user must review the settings of &extra_directives.&
                                     & Directives must not be duplicated.'
           Call info(messages, 2)
           Call error_stop(' ')
         Else
           simulation_data%set_directives%N0=simulation_data%set_directives%N0+1
           simulation_data%set_directives%array(simulation_data%set_directives%N0)=Trim(simulation_data%extra_directives%key(i))
         End If
       End If
     End Do
   End If

   Close(iunit)
   ! create model.param file 
   exec_mv= 'mv '//Trim(files(FILE_SET_SIMULATION)%filename)//' model.param'
   Call execute_command_line(exec_mv)


   ! Open FILE_SET_SIMULATION file (again)
   Open(Newunit=files(FILE_SET_SIMULATION)%unit_no, File=files(FILE_SET_SIMULATION)%filename,Status='Replace')
   iunit=files(FILE_SET_SIMULATION)%unit_no

    Write (iunit,'(a)')  '##############################'
    Write (iunit,'(a)')  '# File generated with ALC_EQCM'
    Write (iunit,'(a)')  '##############################'
    Write (iunit,'(a)') ' '
 
   ! Cell
   If (simulation_data%motion%change_cell_volume%stat .And. simulation_data%motion%change_cell_shape%stat) Then
     Write (iunit, '(a)') 'FIX_ALL_CELL  False  # Allow simulation cell (volume and shape) to change '
   Else If ((.Not. simulation_data%motion%change_cell_volume%stat) .And.&
            (.Not. simulation_data%motion%change_cell_shape%stat)) Then
     Write (iunit, '(a)') 'FIX_ALL_CELL  True   # Fix simulation cell'
   End If

   ! Define constraints for the relaxation of the simulation cell 
   If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
     If ( ((.Not. simulation_data%motion%change_cell_volume%stat) .And. simulation_data%motion%change_cell_shape%stat) .Or. &
           (simulation_data%motion%change_cell_volume%stat .And. (.Not. simulation_data%motion%change_cell_shape%stat)) ) Then
       Write (iunit, '(a)') '#==== Constraints for cell relaxation'
       Write (iunit, '(a)') '%BLOCK CELL_CONSTRAINTS'
       If ( (.Not. simulation_data%motion%change_cell_volume%stat) .And. simulation_data%motion%change_cell_shape%stat) Then 
         Write (iunit, '(a)') '0    0    0'
         Write (iunit, '(a)') '1    2    3'
       ElseIf (simulation_data%motion%change_cell_volume%stat .And. (.Not. simulation_data%motion%change_cell_shape%stat)) Then
         Write (iunit, '(a)') '1    2    3'
         Write (iunit, '(a)') '0    0    0'
       End If 
       Write (iunit, '(a)') '%ENDBLOCK CELL_CONSTRAINTS'
       Write (iunit, '(a)') ' '
     End If
   End If

   ! External pressure 
   If (simulation_data%motion%pressure%fread) Then
     Write (iunit,'(a)') '#==== External pressure'
     Write (iunit, '(a)') '%BLOCK EXTERNAL_PRESSURE'
       Write (iunit, '(a)') 'Mbar'
       Write (iunit,'(    3e20.12)')   simulation_data%motion%pressure%value/1000.0_wp, 0.0_wp, 0.0_wp
       Write (iunit,'(20x,2e20.12)') simulation_data%motion%pressure%value/1000.0_wp, 0.0_wp
       Write (iunit,'(40x,e20.12)')  simulation_data%motion%pressure%value/1000.0_wp
     Write (iunit, '(a)') '%ENDBLOCK EXTERNAL_PRESSURE'
     Write (iunit, '(a)') ' '
   End If  


   ! Hubbard correction
   If (simulation_data%dft%hubbard_info%fread) Then
     Write (iunit,'(a)') '#==== Hubbard corrections'
     Write (iunit,'(a)') '%BLOCK HUBBARD_U'
     Do i=1, net_elements
       j=1
       loop=.True.
       Do While (j <= simulation_data%total_tags .And. loop)
         If (Trim(list_element(i))==Trim(simulation_data%dft%hubbard(j)%element)) Then
           If (Abs(simulation_data%dft%hubbard(j)%U)<epsilon(1.0_wp) .And.&
               Abs(simulation_data%dft%hubbard(j)%J)<epsilon(1.0_wp)) Then
           Else
             If (simulation_data%dft%hubbard(j)%l_orbital == 0) Then
               l_orbital(i)='s'
             Else If (simulation_data%dft%hubbard(j)%l_orbital == 1) Then
               l_orbital(i)='p'
             Else If (simulation_data%dft%hubbard(j)%l_orbital == 2) Then
               l_orbital(i)='d'
             Else If (simulation_data%dft%hubbard(j)%l_orbital == 3) Then
               l_orbital(i)='f'
             End If 
             uc(i)=simulation_data%dft%hubbard(j)%U 
             jc(i)=simulation_data%dft%hubbard(j)%j
             If (simulation_data%dft%pp_info%fread) Then
               Write (iunit,'(2(2x,a),f5.2)') Trim(simulation_data%dft%hubbard(j)%tag), Trim(l_orbital(i))//':  ',&
                                      &  uc(i)-jc(i)   
             Else
               Write (iunit,'(2(2x,a),f5.2)') Trim(simulation_data%dft%hubbard(j)%element), Trim(l_orbital(i))//':  ',&
                                      &  uc(i)-jc(i)   
             End If
           End If
           loop=.False.
         End If
         j=j+1
       End Do      
     End Do
     Write (iunit,'(a)') '%ENDBLOCK HUBBARD_U'
     Write (iunit, '(a)') ' '
   End If

   ! k-point sampling
   Write (iunit,'(a)') '#==== Monkhorst-Pack mesh for k-point sampling'
   Write (iunit,'(a,3i3)') 'KPOINTS_MP_GRID  ',  (simulation_data%dft%kpoints%value(i), i= 1, 3)
   Write (iunit, '(a)') ' '

   ! Pseudo potentials
   Write (iunit,'(a)') '#==== Pseudopotentials'
   Write (iunit,'(a)') '%BLOCK SPECIES_POT'
   Do i=1, net_elements
     If (simulation_data%dft%pp_info%stat) Then
       j=1
       loop=.True.
       Do While (j <= simulation_data%total_tags .And. loop)
         If (Trim(list_element(i))==Trim(simulation_data%dft%pseudo_pot(j)%element)) Then
           pseudo_list(i)= Trim(simulation_data%dft%pseudo_pot(j)%file_name)
           loop=.False.
           If (simulation_data%dft%pp_info%fread) Then  
             Write (iunit,'(1x,a4,4x,a)') Trim(list_tag(i)), Trim(pseudo_list(i))
           Else  
             Write (iunit,'(1x,a2,4x,a)') Trim(list_element(i)), Trim(pseudo_list(i))
           End If
         End If
         j=j+1
       End Do
     Else
       If (simulation_data%dft%pp_info%fread) Then
         Write (iunit,'(1x,a4,1x,a)') Trim(list_tag(i)), 'NCP'
       Else
         Write (iunit,'(1x,a2,1x,a)') Trim(list_element(i)), 'NCP'
       End If 
     End If 
   End Do
   Write (iunit,'(a)') '%ENDBLOCK SPECIES_POT'
   Write (iunit, '(a)') ' '

   If (simulation_data%motion%mass_info%stat) Then
     Write (iunit,'(a)') '#==== Atomic masses'
     Write (iunit,'(a)') '%BLOCK SPECIES_MASS'
     Do i=1, net_elements
       j=1
       loop=.True.
       Do While (j <= simulation_data%total_tags .And. loop)
         If (Trim(list_tag(i))==Trim(simulation_data%motion%mass(j)%tag)) Then
           loop=.False.
           Write (iunit,'(1x,a3,4x,f8.3)') Trim(list_tag(i)), simulation_data%motion%mass(j)%value 
         End If
         j=j+1
       End Do
     End Do
     Write (iunit,'(a)') '%ENDBLOCK SPECIES_MASS'
     Write (iunit, '(a)') ' '
   End If

   Write (iunit,'(a)') '#==== Simulation cell and atomic coordinates'

   Close(iunit)

  End Subroutine print_castep_settings

  Subroutine warnings_castep(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about CASTEP settings 
    !
    ! author    - i.scivetti June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256) :: messages(10), header
    Character(Len=256) :: in_extra
    Logical            :: warning, error, print_header
    Integer(Kind=wi)   :: i
    
    in_extra='using the &extra_directives block'
    warning=.False.
    print_header=.True.

    If (simulation_data%extra_info%fread) Then
      Call info(' ', 1)
      Call info(' ***WARNING****************', 1)
      Write (messages(1), '(1x,a)') ' - ALC_EQCM only checked that the information added in "&extra_directives" has not been'
      Write (messages(2), '(1x,a)') '   already defined from the settings provided in "&block_simulation_settings"'   
      Write (messages(3), '(1x,a)') ' - specification of directives in "&extra_directives" for functionalities that are not'
      Write (messages(4), '(1x,a)') '   implemented in ALC_EQCM have been printed, but its correctness is full responsibility&
                                    & of the user'
      Write (messages(5), '(1x,a)') ' - CASTEP directives added via &extra_directives must only correspond to those of the'
      Write (messages(6), '(1x,a)') '   ".param" file (see CASTEP manual). This is also responsability of the user.'
      Call info(messages, 6)
      Call info(' **************************', 1)
    End If

    Call info(' ', 1)
    Write (messages(1), '(1x,3a)') 'The efficiency in the parallelization can be optimised by changing directives&
                                & OPT_STRATEGY and DATA_DISTRIBUTION ', Trim(in_extra), '.'
    Call info(messages, 1)
    Write (messages(1), '(1x,a)') 'In case of problems in the electronic convergence, the user should:'
    If (simulation_data%dft%smear%type=='fix_occupancy') Then
       Write (messages(2), '(1x,a)') '- check if the system is a metal instead of an insulator.&
                                    & If so, change the option for "smearing"'
      Call info(messages, 2)
    Else
      Write (messages(2), '(1x,a)')  '- increase the value of NEXTRA_BANDS using "bands"'
      Write (messages(3), '(1x,a)')  '- modify the SMEARING_SCHEME using "smearing"'
      Write (messages(4), '(1x,2a)') '- change the method for ELECTRONIC_MINIMIZER (CG by default) ', Trim(in_extra)
      Write (messages(5), '(1x,a)')  '- increase SMEARING_WIDTH using "width_smear"'
      Write (messages(6), '(1x,a)')  '- increase the value CUT_OFF_ENERGY with "energy_cutoff"'
      If (simulation_data%dft%edft%fread) Then
        Write (messages(7), '(1x,2a)') '- increase NUM_OCC_CYCLES ', Trim(in_extra)
        Write (messages(8), '(1x,a)')  ' '
        Write (messages(9), '(1x,a)') 'If none of the above works, set "edft" to .False.'
      Else
        Write (messages(7), '(1x,2a)') '- set a value of MIX_CHARGE_AMP (default is 0.8) ', Trim(in_extra)
        If (simulation_data%dft%mixing%fread) Then
          Write (messages(8), '(1x,a)') '- change the option for directive "mixing_scheme"'
        Else
          Write (messages(8), '(1x,a)') '- set directive "mixing_scheme" (Broyden-2nd by default)'
        End If        
        Write (messages(9), '(1x,a)') 'If none of the above works, set "edft" to .True.'
      End If
      Write (messages(10), '(1x,a)') 'IMPORTANT: the system might be non-metallic. If so, set "smearing" to "fix_occupancy"'
        Call info(messages, 10)
        Call info(' ', 1)
    End If

    Write (messages(1), '(1x,3a)') 'The amount of generated information during run can be controlled (',&
                                  & Trim(in_extra), ') by specification of IPRINT (0,1-default,2,3).'
    Call info(messages, 1)

    If (simulation_data%dft%hubbard_info%fread  .Or. &
       simulation_data%net_charge%fread  .Or. &    
       simulation_data%dft%vdw%fread) Then
       warning=.True.
    End If
    If (Trim(simulation_data%simulation%type) == 'md') Then
      If (Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
        If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover' .Or. &
           Trim(simulation_data%motion%thermostat%type) == 'langevin' ) Then
          warning=.True.
        End If
      Else If (Trim(simulation_data%motion%ensemble%type) == 'npt' .Or. (simulation_data%motion%ensemble%type) == 'nph') Then
        warning=.True.  
      End If
    End If

    If (simulation_data%motion%change_cell_volume%stat) Then
      warning=.True.  
    End If

    If (warning) Then
      Call info(' ', 1)
      Write (header, '(1x,a)')  '***IMPORTANT*** From the requested settings of "&block_simulation_settings", it is&
                                    & RECOMMENDED to consider:'
      ! Total charge
      If (simulation_data%net_charge%fread) Then
        If (Abs(simulation_data%net_charge%value) > epsilon(1.0_wp)) Then
          Write (messages(1), '(1x,a,f8.3,a)')  ' - the meaning/correctness of having set a net charge of ',&
                                  & simulation_data%net_charge%value, &
                                  & ' for the modelled system(s).'
          Call print_warnings(header, print_header, messages, 1)
        End If
      End If

      ! Hubbard-related parameters
      If (simulation_data%dft%hubbard_info%fread) Then
        Write (messages(1), '(1x,2a)')  ' - further optimization of electronic minimization parameters (if there are&
                                  & problems from the inclusion of Hubbard corrections)'
        Call print_warnings(header, print_header, messages, 1)
      End If

     ! Pulay
     If (simulation_data%motion%change_cell_volume%stat) Then
       Write (messages(1), '(1x,a)')  ' - increasing directive "energy_cutoff" to minimise the Pulay stress&
                                      & from changing the cell volume. Finite basis set corrections are applied by default.'
       Call print_warnings(header, print_header, messages, 1)
     End If


      ! MD-related parameters
      If (Trim(simulation_data%simulation%type) == 'md') Then
        If (Trim(simulation_data%motion%ensemble%type) == 'nvt'.Or. &
            Trim(simulation_data%motion%ensemble%type) == 'npt') Then
            Write (messages(1), '(1x,a)')   ' - changing MD_ION_T using "relax_time_thermostat"'
            Call print_warnings(header, print_header, messages, 1)
          If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover') Then
             Write (messages(1), '(1x,2a)') ' - change the number of Nose-Hoover chains (MD_NHC_LENGTH) ', Trim(in_extra)
             Call print_warnings(header, print_header, messages, 1) 
          End If
          If (Trim(simulation_data%motion%ensemble%type) == 'npt') Then
            Write (messages(1), '(1x,a)')  ' - changing MD_CELL_T using "relax_time_barostat"'
            Call print_warnings(header, print_header, messages, 1)
          End If  
        Else If (Trim(simulation_data%motion%ensemble%type) == 'nph') Then
          Write (messages(1), '(1x,a)')  ' - changing MD_CELL_T using "relax_time_barostat"'
          Call print_warnings(header, print_header, messages, 1)
        End If
      End If

      ! vdW related parameters
      If (simulation_data%dft%vdw%fread) Then
        If (Trim(simulation_data%dft%vdw%type) == 'g06') Then
          error=.False.
          Do i=1, simulation_data%total_tags
            If (simulation_data%component(i)%atomic_number > 54) Then
              error=.True.
            End If
          End Do
          If (error) Then 
            Write (messages(1),'(1x,a)')  ' - revision of the requested g06 for vdW corrections:&
                             & the defaults for C06 and R0& 
                             & are defined only for elements in the first five rows of periodic table (i.e. H-Xe).'
            Write (messages(2),'(1x,2a)') '   WARNING: at least one of the defined species are beyond this range and the user&
                                      & must define the correct parameters ', Trim(in_extra) 
            Call print_warnings(header, print_header, messages,2) 
          End If
        End If

        If (Trim(simulation_data%dft%vdw%type) == 'jchs') Then 
          Write (messages(1),'(1x,a)')  ' - revision of the requested JCHS-vdW correction.'
        End If

        If (Trim(simulation_data%dft%vdw%type) == 'obs') Then 
          Write (messages(1),'(1x,a)')  ' - revision of the requested OBS-vdW correction.'
        End If

        If (Trim(simulation_data%dft%vdw%type) == 'ts'  .Or. &
           Trim(simulation_data%dft%vdw%type) == 'mbd' ) Then
          error=.False.
          Do i=1, simulation_data%total_tags
            If (simulation_data%component(i)%atomic_number > 86 .Or. &
              (simulation_data%component(i)%atomic_number > 56 .And. simulation_data%component(i)%atomic_number < 72)) Then
              error=.True.
            End If
          End Do
          If (error) Then 
            Write (messages(1),'(1x,3a)')  ' - revision of the requested "', Trim(simulation_data%dft%vdw%type), '" vdW correction:&
                             & the input reference data for non-interacting atoms is available only for elements of the first&
                             & six rows of the periodic table except of lanthanides.'
            Write (messages(2),'(1x,2a)') '   WARNING: at least one of the defined species are beyond this range. The user&
                                      & must reconsider using this vdW correction. See CASTEP manual at ', Trim(web_castep)
            Call print_warnings(header, print_header, messages,2) 
          End If
            If (Trim(simulation_data%dft%vdw%type) == 'ts') Then
              Write (messages(1),'(1x,3a)')  ' - revision of the requested "', Trim(simulation_data%dft%vdw%type),&
                             & '" vdW correction: this functional is not convenient for small simulation cells.'
              Write (messages(2),'(1x,a)')   '   WARNING: use this functional with caution.' 
              Call print_warnings(header, print_header, messages, 2)
            End If 
          End If

        End If

    End If

  End Subroutine warnings_castep

  Subroutine check_recpot_castep(internal, simulation_data, i)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if PP file is consistent or not with the .recpot format for CASTEP 
    !
    ! author    - i.scivetti September 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   )   :: internal 
    Type(simul_type), Intent(In   )   :: simulation_data
    Integer(Kind=wi), Intent(In   )   :: i 
    
    Character(Len=256) :: message, messages(2)
    Character(Len=256) :: opium, pp_name, xc, element, end_file

    Integer(Kind=wi)   :: io
    Logical            :: loop, is_recpot, fopium
    
    Write (pp_name, '(a)')  Trim(simulation_data%dft%pseudo_pot(i)%file_name)
    Write (end_file,'(1x,2a)') '***ERROR: incomplete (or empty) pseudo potential file ', Trim(pp_name)
          loop=.True.    
    
    fopium=.True.
    is_recpot=.False.
    ! Find the finger print "OPIUM"
    Do While (loop)
      Read (internal, Fmt='(a)', iostat=io) opium
      If (Index(opium,'opium') /= 0 .Or. Index(opium,'OPIUM') /= 0) Then
        loop=.False.
        is_recpot=.True.
      End If
      If (is_iostat_end(io)) Then
        loop=.False.
      End If
    End Do
    If (is_recpot) Then
      Rewind internal
      loop=.True.
      ! read XC correlation
      Do While (loop)
        Read (internal, Fmt=*, iostat=io) xc
        If (is_iostat_end(io)) Then
           Call error_stop(end_file)
        End If
        If (Trim(xc) == '[XC]') Then
          loop=.False.
        End If 
      End Do
      Read (internal, Fmt=*, iostat=io) xc
      If (is_iostat_end(io)) Then
         Call error_stop(end_file)
      End If
      Call capital_to_lower_case(xc)
      If ( Trim(xc) /= 'lda' .And. Trim(xc) /= 'gga') Then
        Write (message, '(1x,2a)') '***ERROR: Unrecognizable XC level for the pseudo potential file ', Trim(pp_name)
        Call error_stop(message)
      End If 
      ! read element
      Rewind internal
      loop=.True.
      Do While (loop)
        Read (internal, *, iostat=io) element 
        If (is_iostat_end(io)) Then
           Call error_stop(end_file)
        End If
        Call capital_to_lower_case(element)
        If (Trim(element) == '[atom]' .Or. Trim(element) == '[Atom]') Then
          loop=.False.
        End If
      End Do
      Read (internal, Fmt=*, iostat=io) element
      If (is_iostat_end(io)) Then
        Call error_stop(end_file)
      End If
      If (Trim(xc) /= Trim(simulation_data%dft%xc_level%type)) Then
        Write (message, '(1x,3a)') '***ERROR: Inconsistency between the XC level of pseudo potential file ', Trim(pp_name),& 
                                & ' and the specification of directive "XC_level".'
        Call error_stop(message) 
      End if 
    Else
      Write (messages(1),'(1x,5a)') '***WARNING: file ', Trim(pp_name), ' for the pseudo potential of species "',&
                          & Trim(simulation_data%dft%pseudo_pot(i)%tag), '" is not recognised as generated by OPIUM&
                          & (http://opium.sourceforge.net/index.html).'
      Write (messages(2),'(1x,a)') '   ALC_EQCM will continue but there is no certainity the&
                                  & pseudpotentials setup will be correct for ONETEP.' 
      Call info(messages, 2)
      fopium=.False.
    End If

    ! More checks
    If ((Trim(element) /= Trim(simulation_data%dft%pseudo_pot(i)%element)) .And. fopium) Then
      Write (message, '(1x,5a)') '***ERROR: pseudo potential file ', Trim(pp_name), ' does not correspond to element "', &
                                & Trim(simulation_data%dft%pseudo_pot(i)%element), '" as specified by the tag&
                                & in sub-block &pseudo_potentials'
      Call error_stop(message) 
    End if 

  End Subroutine check_recpot_castep

  Subroutine check_usp_castep(internal, simulation_data, i)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if PP file is consistent or not with the usp format for CASTEP
    !
    ! author    - i.scivetti September 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    Integer(Kind=wi), Intent(In   )   :: internal 
    Type(simul_type), Intent(In   )   :: simulation_data
    Integer(Kind=wi), Intent(In   )   :: i 
    
    Character(Len=256) :: message, messages(2), string
    Character(Len=256) :: pp_name, xc, element, data_line(2)

    Integer(Kind=wi)   :: io, j
    Logical            :: loop, is_usp
    
    Write (pp_name, '(a)')  Trim(simulation_data%dft%pseudo_pot(i)%file_name)
    loop=.True. 
    is_usp=.False.    
    
    ! Find the line where to read the data
    Do While (loop)
      Read (internal, Fmt='(a)', iostat=io) string
      If (is_iostat_end(io)) Then
        loop=.False.
      End If
      string=Trim(Adjustl(string))
      If (string(1:1) == '|')Then 
        Read (string, Fmt=*) data_line(1), data_line(2)
        If (Trim(data_line(1)) == '|' .And. Trim(data_line(2)) == 'Element:') Then
          loop=.False.
          is_usp=.True.
        End If
        If (is_iostat_end(io)) Then
          loop=.False.
        End If
      End If
    End Do
    If (is_usp) Then
      Backspace internal
      ! read element and XC 
      Read (internal, Fmt=*, iostat=io) (string, j=1,2), element, (string, j=1,6), xc
      Call capital_to_lower_case(xc)
      If (Trim(xc) /= Trim(simulation_data%dft%xc_version%type)) Then
        Write (message, '(1x,3a)') '***ERROR: Inconsistency between the XC version of pseudo potential file ', Trim(pp_name),& 
                                & ' and the specification of directive "XC_version".'
        Call error_stop(message) 
      End if 

    Else
      Write (messages(1),'(1x,5a)') '***ERROR: file ', Trim(pp_name), ' for the pseudo potential of species "',&
                          & Trim(simulation_data%dft%pseudo_pot(i)%tag), '" 1) is empty/incomplete or 2) has a structure that&
                          & is not compatible with to the .usp format valid for CASTEP.'
      Write (messages(2),'(1x,a)') '   If the user wants ultra-soft potentials, please remove &pseuso_potentials. This will&
                                  & instruct CASTEP to generate the pseudo potentials on the fly.' 
      Call info(messages, 2)
      Call error_stop(' ') 
    End If

    ! More checks
    If (Trim(element) /= Trim(simulation_data%dft%pseudo_pot(i)%element)) Then
      Write (message, '(1x,5a)') '***ERROR: pseudo potential file ', Trim(pp_name), ' does not correspond to element "', &
                                & Trim(simulation_data%dft%pseudo_pot(i)%element), '" as specified by the tag&
                                & in sub-block &pseudo_potentials'
      Call error_stop(message) 
    End if 

   
  End Subroutine check_usp_castep 

End Module dft_castep
