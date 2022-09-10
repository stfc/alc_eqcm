!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that automatically sets up input files for geometry relaxation 
! and molecular dynamics of the generated atomistic models 
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author        - i.scivetti April-Oct 2021
! Refactoring   - i.scivetti Oct       2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module simulation_setup

  Use constants,        Only : date_RELEASE

  Use dft_castep,       Only : define_castep_settings,&
                               print_castep_settings, &
                               warnings_castep
  Use dft_cp2k,         Only : define_cp2k_settings,&
                               print_cp2k_settings, &
                               warnings_cp2k
  Use dft_onetep,       Only : define_onetep_settings,&
                               print_onetep_settings, &
                               warnings_onetep
  Use dft_vasp,         Only : define_vasp_settings,&
                               print_vasp_settings,&
                               warnings_vasp 
  Use fileset,          Only : file_type, &
                               FOLDER_DFT,&  
                               FILE_SET_EQCM
  Use numprec,          Only : wi, &
                               wp
  Use references,       Only : bib_ca, bib_hl, bib_pz, bib_wigner, &
                               bib_vwn, bib_pade, bib_pw92, bib_pw91, bib_am05, bib_pbe, bib_rp, bib_revpbe, &
                               bib_pbesol, bib_blyp, bib_wc, bib_xlyp, bib_scan, bib_rpw86pbe, &
                               bib_g06, bib_obs, bib_jchs, bib_dftd2, bib_dftd3, bib_dftd3bj, bib_ts, bib_tsh, bib_mbd,&
                               bib_ddsc, bib_vdwdf, bib_optpbe, bib_optb88, bib_optb86b, bib_vdwdf2, bib_vdwdf2b86r,&
                               bib_SCANrVV10, bib_VV10, bib_AVV10S, bib_tunega, bib_rpw86, bib_fisher, web_d3bJ, &
                               web_d3bJ, web_vasp, web_onetep, web_castep, web_cp2k

  Use simulation_tools, Only : simul_type
  Use stoichiometry,    Only : stoich_type
  Use unit_output,      Only : error_stop,&
                               info 

  Implicit None
  Private

  ! Error in the masses for each species
  Real(Kind=wp), Parameter, Public :: species_mass_error = 0.001_wp
  
  Public :: check_simulation_settings
  Public :: summary_simulation_settings
  Public :: warning_simulation_settings  

Contains

  Subroutine check_simulation_settings(files, stoich_data, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check directives for the generation of input files for 
    ! atomitic level simulations. Assign default values according to the 
    ! requested output format
    !
    ! author    - i.scivetti May-June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(stoich_type),  Intent(In   ) :: stoich_data
    Type(simul_type),   Intent(InOut) :: simulation_data
 
    Integer(Kind=wi)   :: ifolder
    Character(Len=256) :: messages(6)
    Character(Len=256) :: message
    Character(Len=256) :: error_block

    error_block = '***ERROR in &block_simulation_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
    Call info(' ', 1)
    
    Write (messages(1),'(1x,a)')  '***IMPORTANT: by specification of "&block_simulation_settings", the user has'
    Write (messages(2),'(1x,a)')  'requested to generate input files for atomistic level simulations.'
    Call info(messages, 2)

    If (Trim(simulation_data%code_format) /= 'xyz' .And. Trim(simulation_data%code_format) /= 'cif') Then
      Write (messages(1),'(1x,a)') 'From the specification of directive "output_model_format", the generated'
      Write (messages(2),'(1x,3a)')  'files will be consistent with the code "', Trim(simulation_data%code_format), '".'
      Call info(messages, 2)
    Else
      Write (messages(1),'(1x,3a)')  '***ERROR: format ', Trim(simulation_data%code_format),&
                                  & ' (output_model_format) is only valid to generate atomistic&
                                  & structures, but it is not valid as input format for atomistic simulations.'
      Write (messages(2),'(1x,a)')  '          If the user still wants to generate simulation input files,&
                                  & change the option of "output_model_format".'
      Call info(messages, 2)
      Call error_stop(' ')
    End If

    If (Trim(simulation_data%code_format) /= 'vasp' .And. &
       Trim(simulation_data%code_format) /= 'cp2k'  .And. &
       Trim(simulation_data%code_format) /= 'onetep'  .And. &
       Trim(simulation_data%code_format) /= 'castep') Then
      Write (messages(1),'(1x,3a)') '***ERROR: the generation of input files for simulation compatible with the "',&
                                    & Trim(simulation_data%code_format), '" code is NOT implemented (yet).'
      Write (messages(2),'(1x, a)') '          Available options are:'
      Write (messages(3),'(1x, a)') '          - VASP'
      Write (messages(4),'(1x, a)') '          - CP2K'
      Write (messages(5),'(1x, a)') '          - CASTEP'
      Write (messages(6),'(1x, a)') '          - ONETEP'
      Call info(messages, 6)
      Call error_stop(message)
    End If

    ! Type of simulation (compulsory, options: MD or relax_geometry)
    If (simulation_data%simulation%fread) Then
      If (simulation_data%simulation%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) settings for "simulation_type" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%simulation%type) /= 'relax_geometry' .And.&
           Trim(simulation_data%simulation%type) /= 'md'         ) Then
          Write (message,'(1x,4a)') Trim(error_block), &
                                    &' Invalid specification for directive "simulation_type": ',  &
                                    & Trim(simulation_data%simulation%type),& 
                                    &'. Have you missed the specification? Valid options: "relaxation" and "MD"'
          Call error_stop(message)
        End If
      End If
    Else 
      Write (message,'(2(1x,a))') Trim(error_block),&
                               & 'The user must specify directive "simulation_type" (either relaxation or MD)'
      Call error_stop(message)
    End If

    ! Atomic interactions (compulsory, options: DFT for now)
    If (simulation_data%theory_level%fread) Then
      If (simulation_data%theory_level%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) settings for "atomic_interactions" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%theory_level%type) /= 'dft') Then
          Write (message,'(1x,4a)') Trim(error_block), &
                                    &' Invalid specification for directive "atomic_interactions": ',  &
                                    & Trim(simulation_data%theory_level%type),& 
                                    &'. Have you missed the specification? Valid option: "DFT"'
          Call error_stop(message)
        End If
      End If
    Else 
      Write (message,'(2(1x,a))') Trim(error_block),&
                               & 'The user must specify directive "atomic_interactions" (DFT for now)'
      Call error_stop(message)
    End If

    ! Total system charge (optional)
    If (simulation_data%net_charge%fread) Then
      If (simulation_data%net_charge%fail) Then
        Write (message,'(2(1x,a))') Trim(error_block), 'Wrong (or missing) settings for "net_charge" directive.&
                                  & Value must be a real number'
        Call error_stop(message)
      Else
      ! Net charge
        If (Trim(simulation_data%code_format) == 'cp2k' .Or. Trim(simulation_data%code_format) == 'onetep') Then
          If (Abs(simulation_data%net_charge%value-NINT(simulation_data%net_charge%value))> epsilon(1.0_wp)) Then
            If (Trim(simulation_data%code_format) == 'cp2k') Then
              Write (messages(1),'(1x,2a)') Trim(error_block), ' Directive "net_charge" for CP2K simulations must&
                                        & be a number with zero decimals.'
            Else If (Trim(simulation_data%code_format) == 'onetep') Then
              Write (messages(1),'(1x,2a)') Trim(error_block), ' Directive "net_charge" for ONETEP simulations must&
                                        & be a number with zero decimals.'                             
            End If
            Call info(messages, 1)
            Call error_stop(' ')
          End If
        End If
      End If
    End If

    ! Check if blocks have been properly defined
    If (Trim(simulation_data%theory_level%type) == 'dft') Then
      If (simulation_data%dft%pp_info%stat) Then
        ! Check if folder DFT exists 
        Call execute_command_line('[ -d '//Trim(FOLDER_DFT)//' ]', exitstat=ifolder)
        If (ifolder/=0) Then
          Call info(' ', 1)
          Write (messages(1), '(1x,3a)') '***ERROR: folder ', Trim(FOLDER_DFT), ' cannot be found.'
          Write (messages(2), '(1x,a)') 'This folder must contain folder PPs (for pseudo potentials)&
                                       & and BASIS_SET file (for the CP2K code)'
          Write (messages(3), '(1x,a)') 'The requested analysis cannot be conducted. Pleased create the folder&
                                       & and add the required information.'
          Call info(messages, 3)
          Call error_stop(' ')
        End If
      End If

      ! Check DFT settings
      If (simulation_data%dft%generate) Then
        Call check_dft_settings(files, simulation_data)
      Else
        Write (message,'(2(1x,a))') Trim(error_block),&
                               &'From the option of directive "atomic_interactions", the user must define&
                               & sub-block "&dft_settings". See manual for correct syntax.'
        Call error_stop(message)  
      End If
    
    End If 
  
    If (simulation_data%motion%generate) Then 
      ! Ckeck motion settings (for ions)
      Call check_motion_settings(files, stoich_data, simulation_data)
    Else
      Write (message,'(2(1x,a))') Trim(error_block),&
                               & 'The user must define sub-block "&motion_settings" compatible with the choice&
                               & for "simulation_type". See manual for syntax.'
      Call error_stop(message)  
    End If 

    If (simulation_data%extra_info%stat) Then
      If (simulation_data%extra_directives%N0 == 0) Then
        Write (messages(1),'(2(1x,a))') Trim(error_block),&
                                    & 'block &extra_directives is defined but there are no user-defined&
                                    & directives. Either define directives or remove the block'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    Else
      simulation_data%extra_directives%N0 = 0      
    End If

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Assign default values and check compability against the requested format (code)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    If (Trim(simulation_data%code_format) == 'vasp') Then
      Call define_vasp_settings(files, simulation_data)
    Else If (Trim(simulation_data%code_format) == 'cp2k') Then
      Call define_cp2k_settings(files, simulation_data)
    Else If (Trim(simulation_data%code_format) == 'onetep') Then
      Call define_onetep_settings(files, simulation_data)
    Else If (Trim(simulation_data%code_format) == 'castep') Then
      Call define_castep_settings(files, simulation_data)
    End If    

    If (Trim(simulation_data%theory_level%type) == 'dft') Then
      ! Obtain references for the XC_version
      Call obtain_xc_reference(simulation_data%dft%xc_version%type, simulation_data%dft%xc_ref)
      If (simulation_data%dft%vdw%fread) Then
        ! Obtain references for the vdW_correction
        Call obtain_vdw_reference(simulation_data%dft%vdw%type, simulation_data%dft%vdw_ref) 
      End If
    End If 
    
    
  End Subroutine check_simulation_settings
  
  Subroutine obtain_xc_reference(xc_version, xc_reference)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to obtain the reference for the XC functional 
    !
    ! author    - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: xc_version
    Character(Len=*), Intent(  Out) :: xc_reference

    If (Trim(xc_version)     == 'ca'    ) Then 
      xc_reference= 'Ceperley-Alder (CA) '//Trim(bib_ca)
    Else If (Trim(xc_version) == 'hl'    ) Then
      xc_reference= 'Hedin-Lundqvist (HL)'//Trim(bib_hl)       
    Else If (Trim(xc_version) == 'pz'    ) Then
      xc_reference='Perdew-Zunger (PZ) '//Trim(bib_pz)       
    Else If (Trim(xc_version) == 'wigner') Then
      xc_reference= 'Wigner '//Trim(bib_wigner)       
    Else If (Trim(xc_version) == 'vwn'   ) Then
      xc_reference= 'Vosko-Wilk-Nusair (VWN) '//Trim(bib_vwn)       
    Else If (Trim(xc_version) == 'pade'  ) Then
      xc_reference= 'PADE '// Trim(bib_pade)       
    Else If (Trim(xc_version) == 'am05'  ) Then
      xc_reference= 'Armiento-Mattsson (AM05) '//Trim(bib_am05)       
    Else If (Trim(xc_version) == 'pw91'  ) Then
      xc_reference= 'Perdew-Wang 91 (PW91)'//Trim(bib_pw91)       
    Else If (Trim(xc_version) == 'pw92'  ) Then
      xc_reference= 'Perdew-Wang 92 (PW92) '//Trim(bib_pw92)       
    Else If (Trim(xc_version) == 'pbe'   ) Then
      xc_reference= 'Perdew-Burke-Ernzerhof (PBE) '//Trim(bib_pbe)       
    Else If (Trim(xc_version) == 'rp'    ) Then
      xc_reference= 'Hammer-Hansen-Norskov (RP/RPBE) '//Trim(bib_rp)       
    Else If (Trim(xc_version) == 'revpbe') Then
      xc_reference= 'revPBE '//Trim(bib_revpbe)       
    Else If (Trim(xc_version) == 'pbesol') Then
      xc_reference= 'PBE for solids (PBEsol) '//Trim(bib_pbesol)       
    Else If (Trim(xc_version) == 'wc'    ) Then
      xc_reference= 'Wu-Cohen (WC) '//Trim(bib_wc)       
    Else If (Trim(xc_version) == 'blyp'  ) Then
      xc_reference= 'Becke-Lee-Young-Parr (BLYP) '//Trim(bib_blyp)       
    Else If (Trim(xc_version) == 'xlyp'  ) Then
      xc_reference= 'Xu-Goddard (XLYP) '//Trim(bib_xlyp)       
    Else If (Trim(xc_version) == 'or'  ) Then
      xc_reference= 'XC term of the optPBE functional'
    Else If (Trim(xc_version) == 'bo'  ) Then
      xc_reference= 'XC term of the optB88 functional'
    Else If (Trim(xc_version) == 'mk'  ) Then
      xc_reference= 'XC term of the optb86b functional'
    Else If (Trim(xc_version) == 'ml'  ) Then
      xc_reference= 'XC term of the vdW-DF2-b86r functional'
    Else If (Trim(xc_version) == 'scan'  ) Then
      xc_reference= 'meta-GGA functional (SCAN) '//Trim(bib_scan)
    Else If (Trim(xc_version) == 'rpw86pbe'  ) Then
      xc_reference= 'rPW86PBE functional '//Trim(bib_rpw86pbe)
    End If
 
  End Subroutine obtain_xc_reference

  Subroutine obtain_vdw_reference(vdw_type, vdw_reference)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to obtain the refrrence for the selected vdW correction 
    !
    ! author    - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: vdw_type
    Character(Len=*), Intent(  Out) :: vdw_reference

    If (Trim(vdw_type)      == 'dft-d2'      ) Then
      vdw_reference = 'Grimme DFT-D2 '//Trim(bib_dftd2) 
    Else If (Trim(vdw_type)  == 'g06'         ) Then
      vdw_reference = 'Grimme 2006 (g06) '//Trim(bib_g06)
    Else If (Trim(vdw_type)  == 'obs'         ) Then
      vdw_reference = 'Ortmann-Bechstedt-Schmidt (OBS) '//Trim(bib_obs)
    Else If (Trim(vdw_type)  == 'jchs'        ) Then
      vdw_reference = 'Jurecka-Cerny-Hobza-Salahub (JCHS) '//Trim(bib_jchs)
    Else If (Trim(vdw_type)  == 'dft-d3'      ) Then
      vdw_reference = 'Grimme DFT-D3 with no damping '//Trim(bib_dftd3)
    Else If (Trim(vdw_type)  == 'dft-d3-bj'   ) Then
      vdw_reference = 'Grimme D3 with Becke-Jonson damping (DTF-D3-BJ) '//Trim(bib_dftd3bj)
    Else If (Trim(vdw_type)  == 'ts'          ) Then
      vdw_reference = 'Tkatchenko-Scheffler (TS) method '//Trim(bib_ts)
    Else If (Trim(vdw_type)  == 'tsh'         ) Then
      vdw_reference = 'Tkatchenko-Scheffler method with Hirshfeld partitioning (TSH) '//Trim(bib_tsh)
    Else If (Trim(vdw_type)  == 'mbd'         ) Then
      vdw_reference = 'Many-body dispersion energy method '//Trim(bib_mbd)
    Else If (Trim(vdw_type)  == 'ddsc'        ) Then
      vdw_reference = 'DFT-DDsC '//Trim(bib_ddsc)
    Else If (Trim(vdw_type)  == 'vdw-df'      ) Then
      vdw_reference = 'non-local vdW-DF method '//Trim(bib_vdwdf)
    Else If (Trim(vdw_type)  == 'optpbe'      ) Then
      vdw_reference = 'non-local optPBE method '//Trim(bib_optpbe)
    Else If (Trim(vdw_type)  == 'optb88'      ) Then
      vdw_reference = 'non-local optB88 method '//Trim(bib_optb88)
    Else If (Trim(vdw_type)  == 'optb86b'     ) Then
      vdw_reference = 'non-local optB86b method '//Trim(bib_optb86b) 
    Else If (Trim(vdw_type)  == 'vdw-df2'     ) Then
      vdw_reference = 'non-local vdW-DF2 method '//Trim(bib_vdwdf2)
    Else If (Trim(vdw_type)  == 'vdw-df2-b86r') Then
      vdw_reference = 'non-local vdW-DF2-B86R method '//Trim(bib_vdwdf2b86r)
    Else If (Trim(vdw_type)  == 'scan+rvv10'  ) Then
      vdw_reference = 'non-local SCAN+rVV10 method '//Trim(bib_scanrvv10)
    Else If (Trim(vdw_type)  == 'avv10s'      ) Then
      vdw_reference = 'non-local VV10 method '//Trim(bib_vv10)
    Else If (Trim(vdw_type)  == 'vv10'        ) Then
      vdw_reference = 'non-local AVV10S method '//Trim(bib_AVV10s)
    End If

  End Subroutine obtain_vdw_reference


  Subroutine check_dft_settings(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check DFT directives for the generation of input files for 
    ! atomisitic level simulations. 
    !
    ! author    - i.scivetti May-June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256)  :: messages(21)
    Character(Len=256)  :: message
    Character(Len=256)  :: error_dft
    Character(Len=256)  :: to_file, pp_path
    Logical             :: error, safe
    Integer(Kind=wi)    :: i, j, ic, inorm

    error_dft = '***ERROR in &dft_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
    pp_path   = Trim(FOLDER_DFT)//'/PPs/'    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! parallelization (only for VASP)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! parallelization of bands 
    If (simulation_data%dft%npar%fread) Then
      If (simulation_data%dft%npar%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong (or missing) settings for "npar" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%code_format) /= 'vasp') Then
          Write (message,'(1x,4a)') Trim(error_dft), ' Specification of "npar" is incompatible with code "', &
                                    & Trim(simulation_data%code_format), '". Please remove this directive.' 
          Call error_stop(message)
        Else
          If (simulation_data%dft%npar%value <= 0) Then
            Write (message,'(1x,2a)') Trim(error_dft), ' Value for "npar" must be larger than zero. Please correct.'
            Call error_stop(message)
          End If
        End If
      End If
    Else
      If (Trim(simulation_data%code_format) == 'vasp') Then
        Write (message,'(1x,2a)') Trim(error_dft), ' The user must specify directive "npar" for VASP simulations.'
        Call error_stop(message)
      End If 
    End If

    ! parallelization of kpoints 
    If (simulation_data%dft%kpar%fread) Then
      If (simulation_data%dft%kpar%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong (or missing) settings for "kpar" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%code_format) /= 'vasp') Then
          Write (message,'(1x,4a)') Trim(error_dft), ' Specification of "kpar" is incompatible with code "', &
                                    & Trim(simulation_data%code_format), '". Please remove this directive.' 
          Call error_stop(message)
        Else
          If (simulation_data%dft%kpar%value <= 0) Then
            Write (message,'(1x,2a)') Trim(error_dft), ' Value for "kpar" must be larger than zero. Please correct.'
            Call error_stop(message)
          End If
        End If
      End If
    Else
      If (Trim(simulation_data%code_format) == 'vasp') Then
        Write (messages(1),'(1x,2a)') Trim(error_dft), ' The user must specify directive "kpar" for VASP simulations.'
        Write (messages(2),'(1x,a)')  'For large systems, the user is advised to use the Gamma point and "kpar" should&
                                      & be set to 1.'
        Call info(messages, 2)
        Call error_stop(' ')
      End If 
    End If

    !!!!!!!!!!!!!!!!!!!!!!
    ! Electronic structure
    !!!!!!!!!!!!!!!!!!!!!!

    ! XC level
    If (simulation_data%dft%xc_level%fread) Then
      If (simulation_data%dft%xc_level%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong settings for "XC_level" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%dft%xc_level%type) /= 'lda'  .And.&
           Trim(simulation_data%dft%xc_level%type) /= 'gga' ) Then
          Write (message,'(2(1x,a))') Trim(error_dft), &
                                    &'Invalid specification for directive "XC_level". Options are LDA or GGA'
          Call error_stop(message)
        End If
      End If
    Else
      Write (message,'(2(1x,a))') Trim(error_dft), 'The user must specify directive "XC_level" for&
                                & exchange-correlation. Options are LDA or GGA.'
      Call error_stop(message)
    End If

    ! XC version
    Write (messages(2),'(1x,a)')  '==== LDA-level =================='
    Write (messages(3),'(1x,2a)')  '- CA      Ceperley-Alder         ', Trim(bib_ca)
    Write (messages(4),'(1x,2a)')  '- HL      Hedin-Lundqvist        ', Trim(bib_hl)
    Write (messages(5),'(1x,2a)')  '- PZ      Perdew-Zunger          ', Trim(bib_pz)
    Write (messages(6),'(1x,2a)')  '- Wigner  Wigner                 ', Trim(bib_wigner)
    Write (messages(7),'(1x,2a)')  '- VWN     Vosko-Wilk-Nusair      ', Trim(bib_vwn)
    Write (messages(8),'(1x,2a)')  '- PADE    PADE functional        ', Trim(bib_pade)  
    Write (messages(9),'(1x,2a)')  '- PW92    Perdew-Wang 92         ', Trim(bib_pw92)
    Write (messages(10),'(1x,a)')  '==== GGA-level =================='
    Write (messages(11),'(1x,2a)') '- PW91    Perdew-Wang 91         ', Trim(bib_pw91)
    Write (messages(12),'(1x,2a)') '- AM05    Armiento-Mattsson      ', Trim(bib_am05)
    Write (messages(13),'(1x,2a)') '- PBE     Perdew-Burke-Ernzerhof ', Trim(bib_pbe)
    Write (messages(14),'(1x,2a)') '- RP      Hammer-Hansen-Norskov  ', Trim(bib_rp)
    Write (messages(15),'(1x,2a)') '- revPBE  revPBE                 ', Trim(bib_revpbe)
    Write (messages(16),'(1x,2a)') '- PBEsol  PBE for solids         ', Trim(bib_pbesol)
    Write (messages(17),'(1x,2a)') '- BLYP    Becke-Lee-Young-Parr   ', Trim(bib_blyp)
    Write (messages(18),'(1x,2a)') '- WC      Wu-Cohen               ', Trim(bib_wc)
    Write (messages(19),'(1x,2a)') '- XLYP    Xu-Goddard             ', Trim(bib_xlyp)
    Write (messages(20),'(1x,2a)') '================================='
    Write (messages(21),'(1x,2a)') 'IMPORTANT: not all the above XC functionals are implemented for all&
                                  & DFT codes. Please refer to the ALC_EQCM manual for details.'
    If (simulation_data%dft%xc_version%fread) Then
      If (simulation_data%dft%xc_version%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft),&
                                      & 'Wrong settings for "XC_version" directive. Implemented options are:'
        Call info(messages, 21)
        Call error_stop(' ')
      Else
        If (Trim(simulation_data%dft%xc_version%type) /= 'ca'     .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'hl'     .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'pz'     .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'wigner' .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'vwn'    .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'pade'   .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'am05'   .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'pw91'   .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'pw92'   .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'pbe'    .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'rp'     .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'revpbe' .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'pbesol' .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'wc'     .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'blyp'   .And.&
           Trim(simulation_data%dft%xc_version%type) /= 'xlyp' ) Then
          Write (messages(1),'(2(1x,a))') Trim(error_dft), 'Unrecognised specification for directive&
                                         & "XC_version". Implemented options are:'
          Call info(messages, 21)
          Call error_stop(' ')
        End If
          If (Trim(simulation_data%dft%xc_level%type) == 'lda') Then
             If (Trim(simulation_data%dft%xc_version%type) /= 'ca' .And.&
                Trim(simulation_data%dft%xc_version%type) /= 'hl' .And.&
                Trim(simulation_data%dft%xc_version%type) /= 'pz' .And.&
                Trim(simulation_data%dft%xc_version%type) /= 'wigner' .And.&
                Trim(simulation_data%dft%xc_version%type) /= 'pade' .And.&
                Trim(simulation_data%dft%xc_version%type) /= 'pw92' .And.&
                Trim(simulation_data%dft%xc_version%type) /= 'vwn' ) Then
               Write (messages(1),'(1x,4a)') Trim(error_dft), &
                          & ' Directive "XC_version" has been set to "', Trim(Trim(simulation_data%dft%xc_version%type)),&
                          & '", which is incompatible with LDA. Please review the XC settings'
               Call info(messages, 1)
               Call error_stop(' ')
             End If 
          Else If (Trim(simulation_data%dft%xc_level%type) == 'gga' ) Then
            If (Trim(simulation_data%dft%xc_version%type) /= 'am05'   .And.&
               Trim(simulation_data%dft%xc_version%type)  /= 'pw91'   .And.&
               Trim(simulation_data%dft%xc_version%type)  /= 'pbe'    .And.&
               Trim(simulation_data%dft%xc_version%type)  /= 'rp'     .And.&
               Trim(simulation_data%dft%xc_version%type)  /= 'revpbe' .And.&
               Trim(simulation_data%dft%xc_version%type)  /= 'pbesol' .And.&
               Trim(simulation_data%dft%xc_version%type)  /= 'wc'     .And.&
               Trim(simulation_data%dft%xc_version%type)  /= 'xlyp'   .And.&
               Trim(simulation_data%dft%xc_version%type)  /= 'blyp' ) Then
               Write (messages(1),'(1x,4a)') Trim(error_dft), &
                          &' Directive "XC_version" has been set to "', Trim(Trim(simulation_data%dft%xc_version%type)),&
                          & '", which is incompatible with GGA. Please review the XC settings'
               Call info(messages, 1)
               Call error_stop(' ')
            End If
          End If
      End If
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_dft), 'The user must specify directive "XC_version" for&
                               & exchange-correlation. Implemented options are:'
      Call info(messages, 21)
      Call error_stop(' ')
    End If

    ! Energy cutoff (compulsory)
    If (simulation_data%dft%encut%fread) Then
      If (simulation_data%dft%encut%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong (or missing) settings for "energy_cutoff" directive.'
        Call error_stop(message)
      Else
        If (simulation_data%dft%encut%value < epsilon(simulation_data%dft%encut%value)) Then
          Write (message,'(2(1x,a))') Trim(error_dft), &
                                    &'Input value for "energy_cutoff" MUST be larger than zero'
          Call error_stop(message)
        End If
        If (Trim(simulation_data%dft%encut%units) /= 'ev' .And. &
           Trim(simulation_data%dft%encut%units) /= 'ry') Then
           Write (message,'(1x,4a)') Trim(error_dft), &
                                    &' Invalid specification for the units of "energy_cutoff": ', &
                                    &  Trim(simulation_data%dft%encut%units),&
                                    &'. Units MUST BE in eV or Ry. Have you missed data in the specification?'
          Call error_stop(message)                          
        End If    
      End If
    Else  
      Write (message,'(2(1x,a))') Trim(error_dft), 'The user must specify directive "energy_cutoff"&
                               & (value and units, see manual)'
      Call error_stop(message)
    End If

    ! Smearing 
    If (simulation_data%dft%smear%fread) Then
      If (simulation_data%dft%smear%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong (or missing) settings for "smearing" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%dft%smear%type) /= 'gaussian'      .And.&
           Trim(simulation_data%dft%smear%type)  /= 'fermi'         .And.&
           Trim(simulation_data%dft%smear%type)  /= 'mp'            .And.&
           Trim(simulation_data%dft%smear%type)  /= 'window'        .And.&
           Trim(simulation_data%dft%smear%type)  /= 'fix_occupancy' .And.&
           Trim(simulation_data%dft%smear%type)  /= 'tetrahedron'  ) Then
          Write (messages(1),'(4(1x,a))') Trim(error_dft), &
                                    &'Invalid specification of directive "smearing":', Trim(simulation_data%dft%smear%type), &
                                    &'. Options are:'
          Write (messages(2),'(1x,a)') '- Gaussian       (Gaussian distribution)' 
          Write (messages(3),'(1x,a)') '- Fermi          (Fermi-Dirac distribution)'
          Write (messages(4),'(1x,a)') '- Tetrahedron    (Tetrahedron method with Blochl corrections)'
          Write (messages(5),'(1x,a)') '- Window         (Energy window centered at the Fermi level)' 
          Write (messages(6),'(1x,a)') '- Fix-Occupancy  (Fix the occupancies of the bands)' 
          Write (messages(7),'(1x,a)') '- MP             (Methfessel-Paxton method)'
          Call info(messages, 7)
          Call error_stop(' ')
        End If
      End If
    End If

    ! Basis sets (compulsory for basis set codes)
    If (Trim(simulation_data%code_format) == 'cp2k') Then
      If (.Not. simulation_data%dft%basis_info%fread) Then
         Write (messages(1),'(1x,a)') ' '
         Write (messages(2),'(1x,a)') '**********************************************************************'
         Write (messages(3),'(1x,a)') '*** WARNING: the "&basis_set" sub-block has not been defined.      ***'
         Write (messages(4),'(1x,a)') '***          ALC_EQCM will not generate files with the basis sets! ***'
         Write (messages(5),'(1x,a)') '**********************************************************************'
         Call info(messages, 5)
      Else
        ! Check if user has included all the tags
        Do i=1, simulation_data%total_tags
          error=.True.
          Do j=1, simulation_data%total_tags
            If (Trim(simulation_data%dft%basis_set(i)%tag)==Trim(simulation_data%component(j)%tag)) Then
              simulation_data%dft%basis_set(i)%element=simulation_data%component(j)%element
              error=.False.
            End If
          End Do
          If (error) Then
            Write (message,'(1x,a,1x,3a)') Trim(error_dft), 'Atomic tag "', Trim(simulation_data%dft%pseudo_pot(i)%tag),&
                                     & '" declared in "&basis_set" has not been defined in&
                                     & "&block_input_composition". Please check'
            Call error_stop(message)
          End If
          ! Check if it has the right type
          If (Trim(simulation_data%dft%basis_set(i)%type) /= 'sz'  .And. & 
             Trim(simulation_data%dft%basis_set(i)%type) /= 'dz'  .And. &
             Trim(simulation_data%dft%basis_set(i)%type) /= 'szp' .And. &
             Trim(simulation_data%dft%basis_set(i)%type) /= 'dzp' .And. &
             Trim(simulation_data%dft%basis_set(i)%type) /= 'tzp' .And. &
             Trim(simulation_data%dft%basis_set(i)%type) /= 'tz2p' ) Then
             Write (messages(1),'(1x,4a)') Trim(error_dft), ' Invalid basis set specification for tag "', &
                                      & Trim(simulation_data%dft%basis_set(i)%tag), '". Valid options:'
             Write (messages(2),'(1x,a)') 'SZ   (Single Zeta)              '
             Write (messages(3),'(1x,a)') 'DZ   (Double Zeta)              ' 
             Write (messages(4),'(1x,a)') 'SZP  (Single Zeta Polarizable)  '
             Write (messages(5),'(1x,a)') 'DZP  (Double Zeta Polarizable)  '
             Write (messages(6),'(1x,a)') 'TZP  (Triple Zeta Polarizable)  '
             Write (messages(7),'(1x,a)') 'TZ2P (Triple Zeta 2-Polarizable)'
             Call info(messages, 7)
             Call error_stop(' ')
          End If
        End Do
      End If
    End If


    ! Pseudopotentials (compulsory for onetep)
    If (.Not. simulation_data%dft%pp_info%stat) Then
      Write (messages(1),'(1x,a)') ' '
      Write (messages(2),'(1x,a)') '*************************************************************************'
      Write (messages(3),'(1x,a)') '*** WARNING: the "&pseudo_potentials" sub-block has not been defined. ***'
      Write (messages(4),'(1x,a)') '***          ALC_EQCM will not generate pseudopotential files!        ***'
      Write (messages(5),'(1x,a)') '*************************************************************************'
      Call info(messages, 5)
    Else
      ! Check if user has included all the tags 
      Do i=1, simulation_data%total_tags 
        error=.True.
        Do j=1, simulation_data%total_tags
          If (Trim(simulation_data%dft%pseudo_pot(i)%tag)==Trim(simulation_data%component(j)%tag)) Then
            simulation_data%dft%pseudo_pot(i)%element=simulation_data%component(j)%element
            error=.False.  
          End If
        End Do
        If (error) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_dft), 'Atomic tag "', Trim(simulation_data%dft%pseudo_pot(i)%tag),&
                                   & '" declared in "&pseudo_potentials" has not been defined in&
                                   & "&block_input_composition". Please check'
          Call error_stop(message)                                      
        End If
      End Do
      ! Check if all pseudopotential files exist
      Do i=1, simulation_data%total_tags
        safe=.False.
        to_file=Trim(pp_path)//Trim(simulation_data%dft%pseudo_pot(i)%file_name)
        Inquire(File=Trim(to_file), Exist=safe)  
        If (.not. safe) Then
          Write (message,'(1x,6a)') '***ERROR: File ', Trim(pp_path), Trim(simulation_data%dft%pseudo_pot(i)%file_name),&
                                   & ' not found. Please check the defined pseudo potentials in "&pseudo_potentials"&
                                   & and files in folder ', Trim(pp_path),& 
                                   & '. Does PPs subfolder exist? If not, create it and copy the relevant files into.'
          Call error_stop(message)
        End If
      End Do
    End If

    ! Width smearing (optional)
    If (simulation_data%dft%width_smear%fread) Then
      If (simulation_data%dft%width_smear%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong (or missing) settings for "width_smear" directive.'
        Call error_stop(message)
      Else
        If (simulation_data%dft%width_smear%value < epsilon(simulation_data%dft%width_smear%value)) Then
          Write (message,'(2(1x,a))') Trim(error_dft), &
                                    &'Input value for "width_smear" MUST be larger than zero'
          Call error_stop(message)
        End If  
        If (Trim(simulation_data%dft%width_smear%units) /= 'ev') Then
           Write (message,'(4a)')  Trim(error_dft), ' Invalid units of directive "width_smear": ', &
                                  Trim(simulation_data%dft%width_smear%units), '. Units must be eV'
           Call error_stop(message)
        End If
      End If
    End If

    ! SCF energy tolerance (optional) 
    If (simulation_data%dft%delta_e%fread) Then
      If (simulation_data%dft%delta_e%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong settings for "SCF_energy_tolerance" directive.'
        Call error_stop(message)
      Else
        If (simulation_data%dft%delta_e%value < epsilon(simulation_data%dft%delta_e%value)) Then
          Write (message,'(2(1x,a))') Trim(error_dft), &
                                    &'Input value for "SCF_energy_tolerance" MUST be larger than zero'
          Call error_stop(message)
        End If
        If (Trim(simulation_data%dft%delta_e%units) /= 'ev' .And. &
           Trim(simulation_data%dft%delta_e%units) /= 'hartree' ) Then
           Write (message,'(4a)')  Trim(error_dft), ' Invalid units of directive "SCF_energy_tolerance": ', &
                                  Trim(simulation_data%dft%delta_e%units), '. Check manual for the valid units'
           Call error_stop(message)
        End If   
      End If
    Else
      Write (message,'(2(1x,a))') Trim(error_dft), 'The user must specify directive "SCF_energy_tolerance" to&
                                & define the level of convergence for the electronic problem.'
      Call error_stop(message)
    End If

    ! SCF steps for electronic convergence (optional)
    If (simulation_data%dft%scf_steps%fread) Then
      If (simulation_data%dft%scf_steps%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong (or missing) settings for "scf_steps" directive'
        Call error_stop(message)
      Else
        If (simulation_data%dft%scf_steps%value < 0) Then
          Write (message,'(1x,2a)') Trim(error_dft), ' Number of "scf_steps" must be positive and&
                                                      & sufficiently large to reach electronic convergence&
                                                      & (this is responsibility of the user)'
          Call error_stop(message)
        End If
      End If
    Else
      Write (message,'(2(1x,a))') Trim(error_dft), 'The user must specify directive "scf_steps" to&
                                & define the maximum amount of iterations for the electronic convergence.'
      Call error_stop(message)
    End If

    ! NGWF information compulsory only for ONETEP
    If (simulation_data%dft%ngwf_info%fread) Then
      If (Trim(simulation_data%code_format) /= 'onetep') Then
        Write (message,'(1x,a,1x,3a)')  Trim(error_dft), 'Sub-block &ngwf is irrelevant for the requested code format "', &
                                     & Trim(simulation_data%code_format), '" and must be removed. Specitication of &ngwf&
                                     & is only compulsory for ONETEP.'
        Call error_stop(message)
      Else
        ! Check if user has included all the tags
        Do i=1, simulation_data%total_tags
          error=.True.
          Do j=1, simulation_data%total_tags
            If (Trim(simulation_data%dft%ngwf(i)%tag)==Trim(simulation_data%component(j)%tag)) Then
              simulation_data%dft%ngwf(i)%element=simulation_data%component(j)%element
              error=.False.
            End If
          End Do
          If (error) Then
            Write (message,'(1x,a,1x,3a)') Trim(error_dft),&
                                     & 'Atomic tag "', Trim(simulation_data%dft%ngwf(i)%tag), '" declared in&
                                     & "&ngwf" has not been defined in "&block_input_composition". Please check.'
            Call error_stop(message)
          End If
        End Do

        ! Check values for number and radius
        Do i=1, simulation_data%total_tags
          If (simulation_data%dft%ngwf(i)%ni == 0 .Or. simulation_data%dft%ngwf(i)%ni < -1) Then
            Write (message,'(1x,a,1x,3a)') Trim(error_dft),&
                                     & 'Number of NGWFs for atomic tag "', Trim(simulation_data%dft%ngwf(i)%tag),&
                                     & ' must be an integer larger than 0. The user can also set "-1" for default values&
                                     & (check block &ngwf)'
            Call error_stop(message)
          End If
          If (simulation_data%dft%ngwf(i)%radius <= 0.0_wp) Then
            Write (message,'(1x,a,1x,3a)') Trim(error_dft),&
                                     & 'The radius for the NGWF of atomic tag "', Trim(simulation_data%dft%hubbard(i)%tag),&
                                     & '" must be positive (check block &ngwf)'
            Call error_stop(message)
          End If
          If (simulation_data%dft%ngwf(i)%radius <= 1.0_wp) Then
            Write (message,'(1x,3a,f8.2,a)') '***WARNING: The radius for the NGWF of atomic tag "',&
                                     & Trim(simulation_data%dft%hubbard(i)%tag),&
                                     & '" has been set to ', simulation_data%dft%ngwf(i)%radius, &
                                     &' Angstrom, which is a rather small value (check block &ngwf).'
            Call info(message, 1)
          End If
        End Do
      End If

    Else
      If (Trim(simulation_data%code_format) == 'onetep') Then
        Write (message,'(1x,a,1x,3a)')  Trim(error_dft), 'For the generation of input files for simulations with the code "', &
                                     & Trim(simulation_data%code_format), '", it is compulsory the user defines&
                                     & sub-block &ngwf. Please refer to the manual for correct format and syntax.'
        Call error_stop(message)
      End If
    End If

    ! Orbital transformation (optional vs smearing, only for CP2K)
    If (simulation_data%dft%ot%fread) Then
      If (simulation_data%dft%ot%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong settings for "OT" directive&
                                 & (Orbital Transformation). Set either .True. or .False.'
        Call error_stop(message)
      End If
    Else
      ! By default, there is no OT 
      simulation_data%dft%ot%stat=.False.
    End If

    ! Ensemble DFT (optional vs density mixing, only for CASTEP and ONETEP)
    If (simulation_data%dft%edft%fread) Then
      If (simulation_data%dft%edft%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong settings for "EDFT" directive&
                                 & (Ensemble-DFT). Set either .True. or .False.'
        Call error_stop(message)
      End If
    Else
      ! By default, there is no EDFT 
      simulation_data%dft%edft%stat=.False.
    End If

    ! Check the code vs EDFT
    If (simulation_data%dft%edft%stat) Then
      If (Trim(simulation_data%code_format) /= 'onetep' .And. &
          Trim(simulation_data%code_format) /= 'castep' ) Then
        Write (message,'(1x,a,1x,3a)')  Trim(error_dft), 'eDFT simulations are not possible for the requested code format "', &
                               & Trim(simulation_data%code_format), '". Please remove the "edft" directive or set it to .False.'
        Call error_stop(message)
      End If            
    End If 

    ! GC-DFT functionality
    If (simulation_data%dft%gc%activate%stat) Then
      If (Trim(simulation_data%code_format) /= 'onetep' ) Then
        Write (message,'(1x,a,1x,3a)')  Trim(error_dft), 'GC-DFT simulations are not possible for the requested code format "', &
                & Trim(simulation_data%code_format), '". Please remove the &gcdft sub-block.'
        Call error_stop(message)
      End If
      If (Trim(simulation_data%process)/= 'electrodeposition') Then
        Write (message,'(1x,a,1x,3a)')  Trim(error_dft), 'GC-DFT simulations are only possible for surfaces systems and&
                                        & electrodeposition processes.'
        Call error_stop(message)
      End If        
    End If        

    ! precision (only compulsory for VASP)
    If (simulation_data%dft%precision%fread) Then
      If (simulation_data%dft%precision%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong settings for "precision" directive.'
        Call error_stop(message)
      End If
    End If

    ! vdW settings (optional)
    If (simulation_data%dft%vdw%fread) Then
      If (simulation_data%dft%vdw%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong settings for "vdW" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%dft%vdw%type) /= 'dft-d2'    .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'g06'       .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'obs'       .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'jchs'      .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'dft-d3'    .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'dft-d3-bj' .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'ts'        .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'tsh'       .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'mbd'       .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'ddsc'      .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'vdw-df'    .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'optpbe'   .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'optb88'   .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'optb86b'  .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'vdw-df2'      .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'vdw-df2-b86r' .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'scan+rvv10'   .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'avv10s'   .And.&
          Trim(simulation_data%dft%vdw%type)  /= 'vv10'   ) Then
          Write (messages(1),'(1x,4a)') Trim(error_dft), &
                                  & ' Invalid specification of directive "vdW": ', Trim(simulation_data%dft%vdw%type),&
                                  &'. Valid options are:'
          Write (messages(2),'(1x,2a)')  '- G06            Grimme 2006 ', Trim(bib_g06)
          Write (messages(3),'(1x,2a)')  '- OBS            Ortmann-Bechstedt-Schmidt ', Trim(bib_obs)
          Write (messages(4),'(1x,2a)')  '- JCHS           Jurecka-Cerny-Hobza-Salahub ', Trim(bib_jchs)
          Write (messages(5),'(1x,2a)')  '- DFT-D2         Grimme D2 ', Trim(bib_dftd2)
          Write (messages(6),'(1x,2a)')  '- DFT-D3         Grimme D3 with no damping ', Trim(bib_dftd3)
          Write (messages(7),'(1x,2a)')  '- DFT-D3-BJ      Grimme D3 with Becke-Jonson damping ', Trim(bib_dftd3bj)
          Write (messages(8),'(1x,2a)')  '- TS             Tkatchenko-Scheffler method ', Trim(bib_ts)
          Write (messages(9),'(1x,2a)')  '- TSH            TS method with Hirshfeld partitioning ', Trim(bib_tsh)
          Write (messages(10),'(1x,2a)') '- MBD            Many-body dispersion energy method ', Trim(bib_mbd)
          Write (messages(11),'(1x,2a)') '- dDsC           DFT-DDsC ', Trim(bib_ddsc)
          Write (messages(12),'(1x,2a)') '- vdW-DF         X (revPBE), C (LDA), vdW (vdW-DF) ', Trim(bib_vdwdf)  
          Write (messages(13),'(1x,2a)') '- optPBE         X (OPTPBE), C (LDA), vdW (vdW-DF) ', Trim(bib_optpbe)    
          Write (messages(14),'(1x,2a)') '- optB88         X (OPTB88), C (LDA), vdW (vdW-DF) ', Trim(bib_optb88)
          Write (messages(15),'(1x,2a)') '- optB86b        Optimized Becke86b ', Trim(bib_optb86b)
          Write (messages(16),'(1x,2a)') '- vdW-DF2        X (rPW86), C (LDA), vdW (vdW-DF 2)  ', Trim(bib_vdwdf2)
          Write (messages(17),'(1x,2a)') '- vdW-DF2-B86R   Hamada version of vdW-DF2 ', Trim(bib_vdwdf2b86r)
          Write (messages(18),'(1x,2a)') '- SCAN+rVV10     SCAN + non-local correlation part of the rVV10 ', Trim(bib_scanrvv10)
          Write (messages(19),'(1x,2a)') '- VV10           X (rPW86), C (PBE), vdW (rVV10) ', Trim(bib_vv10)
          Write (messages(20),'(1x,2a)') '- AVV10S         X (AM05), C (AM05), vdW (rVV10-sol) ', Trim(bib_AVV10s)
          Write (messages(21),'(1x,2a)') 'IMPORTANT: not all the above vdW functionals are implemented for all&
                                  & DFT codes. Please refer to the ALC_EQCM manual for details.'
          Call info(messages, 20)
          Call error_stop(' ')
        End If
      End If
    End If

    ! Prevent a couple of combinations
    If (simulation_data%dft%vdw%fread) Then
      If (Trim(simulation_data%dft%vdw%type) == 'dft-d2' .Or. &
          Trim(simulation_data%dft%vdw%type) == 'dft-d3' .Or. &
          Trim(simulation_data%dft%vdw%type) == 'dft-d3' ) Then
          If (Trim(simulation_data%dft%xc_version%type) == 'am05') Then
            Write (messages(1),'(1x,4a)')  Trim(error_dft), ' Until ', Trim(date_RELEASE), ', there was no evidence of previous&
                                        & DFT-D2/DFT-D3/DFT-D3-BJ parametrization/simulations with the "AM05" XC functional.&
                                        & The user is advised to consider a different XC functional.'
            Call info(messages, 1)
            Call error_stop(' ') 
          End If
      End If    

      If (Trim(simulation_data%dft%vdw%type) == 'dft-d2') Then
          If (Trim(simulation_data%dft%xc_version%type) == 'pbesol') Then
            Write (messages(1),'(1x,2a)')  Trim(error_dft), ' The DFT-D2 parametrization with the "PBEsol" XC&   
                                & functional is not implemented. The user is advised to consider a different XC functional.'
            Call info(messages, 1)
            Call error_stop(' ') 
          End If
      End If
    End If
    
    ! Bands for electronic structure calculation              
    If (simulation_data%dft%bands%fread) Then
      If (simulation_data%dft%bands%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong (or missing) settings for "bands" directive'
        Call error_stop(message)
      Else
        If (simulation_data%dft%bands%value < 1) Then
          Write (message,'(1x,2a)') Trim(error_dft), ' Value for "bands" must be >= 1. Please change'
          Call error_stop(message)
        End If
      End If
    End If

    ! Maximum l_orbital (optional, only for VASP)
    If (simulation_data%dft%max_l_orbital%fread) Then
      If (simulation_data%dft%max_l_orbital%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong (or missing) settings for "max_l_orbital" directive'
        Call error_stop(message)
      Else
        If (simulation_data%dft%max_l_orbital%value < 0 .Or. simulation_data%dft%max_l_orbital%value > 3 ) Then
          Write (message,'(1x,2a)') Trim(error_dft), ' Wrong value for max_l_orbital. Allowed values: 0(s), 1(p), 2(d) and 3(f)'
          Call error_stop(message)
        End If
      End If
    End If

    ! Spin polarised (optional)
    If (simulation_data%dft%spin_polarised%fread) Then
      If (simulation_data%dft%spin_polarised%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong settings for "spin_polarised" directive.&
                                 & Either .True. or .False.'
        Call error_stop(message)
      End If
    Else
      ! By default, settings are non spin-polarised
      simulation_data%dft%spin_polarised%stat=.False.
    End If

    !Magnetization (optional)
    If (simulation_data%dft%mag_info%fread) Then
      If (.Not. simulation_data%dft%spin_polarised%stat) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Inclusion of magnetization for the atoms requires to set&
                                     & "spin_polarised" to .True.'
        Call error_stop(message)
      End If
      ! Check if user has included all the tags
      Do i=1, simulation_data%total_tags
        error=.True.
        Do j=1, simulation_data%total_tags
          If (Trim(simulation_data%dft%magnetization(i)%tag)==Trim(simulation_data%component(j)%tag)) Then
            simulation_data%dft%magnetization(i)%element=simulation_data%component(j)%element
            error=.False.
          End If
        End Do
        If (error) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_dft),&
                                   & 'Atomic tag "', Trim(simulation_data%dft%magnetization(i)%tag), '" declared in&
                                   & "&magnetization" has not been defined in "&block_input_composition". Please check'
          Call error_stop(message)
        End If
      End Do
      ! Check if initial all magnetizations are zero, in which case ALC_EQCM aborts 
      error=.True.
      Do i=1, simulation_data%total_tags
        If (Abs(simulation_data%dft%magnetization(i)%value) > epsilon(1.0_wp)) Then
           error=.False.
        End If
      End Do
      If (error) Then
          Write (message,'(1x,2a)') Trim(error_dft), ' All initial magnetic moments in block &magnetization are set to zero.&
                                   & Please check. If no magnetization is required, remove the block.'
          Call error_stop(message)
       End If
    End If

    ! total_magnetization (optional)
    If (simulation_data%dft%total_magnetization%fread) Then
      If (simulation_data%dft%total_magnetization%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong (or missing) settings for "total_magnetization" directive'
        Call error_stop(message)
      Else
        If (.Not. simulation_data%dft%spin_polarised%stat) Then
          Write (message,'(2(1x,a))') Trim(error_dft), 'Specification of "total_magnetization" requires of a spin polarised&
                                     & simulation. Please set "spin_polarised" to .True.'
          Call error_stop(message)
        End If
        If (.Not. simulation_data%dft%mag_info%fread) Then
          Write (message,'(2(1x,a))') Trim(error_dft), 'Specification of "total_magnetization" requires the definition of&
                                     & initial magnetic moments(spins) via sub-block &magnetization.'
          Call error_stop(message)
        End If
       ! Net magnetization
        If (Trim(simulation_data%code_format) == 'onetep') Then
          If (Abs(simulation_data%dft%total_magnetization%value- &
            & NINT(simulation_data%dft%total_magnetization%value))> epsilon(1.0_wp)) Then
            If (.Not. simulation_data%dft%edft%stat) Then
              Write (messages(1),'(1x,2a)') Trim(error_dft), ' In ONETEP, unless EDFT is selected, directive&
                                        & "total_magnetization" must be a number with zero decimals.'
              Call info(messages, 1)
              Call error_stop(' ')
            End If
          End If
        End If
      End If
    End If

    !Hubbard (optional)
    If (simulation_data%dft%hubbard_info%fread) Then
      If (.Not. simulation_data%dft%mag_info%fread) Then
        Write (message,'(1x,a,1x,a)') Trim(error_dft), 'Inclusion of Hubbard corrections for atoms requires of initial&
                                     & magnetic moments via "&magnetization"'
        Call error_stop(message)
      End If  
      ! Check if user has included all the tags
      Do i=1, simulation_data%total_tags
        error=.True.
        Do j=1, simulation_data%total_tags
          If (Trim(simulation_data%dft%hubbard(i)%tag)==Trim(simulation_data%component(j)%tag)) Then
            simulation_data%dft%hubbard(i)%element=simulation_data%component(j)%element
            error=.False.
          End If
        End Do
        If (error) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_dft),&
                                   & 'Atomic tag "', Trim(simulation_data%dft%hubbard(i)%tag), '" declared in&
                                   & "&hubbard" has not been defined in "&block_input_composition". Please check'
          Call error_stop(message)
        End If
      End Do

      ic=0
      simulation_data%dft%hubbard_all_U_zero=.False.
      ! Check values for l_orbital, J and U
      Do i=1, simulation_data%total_tags
        If (simulation_data%dft%hubbard(i)%l_orbital < 0 .And. simulation_data%dft%hubbard(i)%l_orbital > 3) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_dft),&
                                   & 'Value of l_orbital for atomic tag "', Trim(simulation_data%dft%hubbard(i)%tag), '" must&
                                   & be either 0, 1, 2 and 3 (check block &hubbard)'
          Call error_stop(message)
        End If
        If (simulation_data%dft%hubbard(i)%U < 0.0_wp) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_dft),&
                                   & 'Value of U for atomic tag "', Trim(simulation_data%dft%hubbard(i)%tag), '" must&
                                   & be positive or zero (check block &hubbard)'
          Call error_stop(message)
        End If

        If (Abs(simulation_data%dft%hubbard(i)%U) < epsilon(1.0_wp)) Then
          ic=ic+1
        End If

        If (simulation_data%dft%hubbard(i)%J < 0.0_wp) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_dft),&
                                   & 'Value of J for atomic tag "', Trim(simulation_data%dft%hubbard(i)%tag), '" must&
                                   & be positive or zero (check block &hubbard)'
          Call error_stop(message)
        End If

        If ((simulation_data%dft%hubbard(i)%U - simulation_data%dft%hubbard(i)%J) < 0.0_wp) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_dft),&
                                   & 'The value of U for atomic tag "', Trim(simulation_data%dft%hubbard(i)%tag), '" must&
                                   & be larger than J (check block &hubbard)'
          Call error_stop(message)
        End If

      End Do
      
      If (ic==simulation_data%total_tags) Then
         simulation_data%dft%hubbard_all_U_zero=.True.
         If (Trim(simulation_data%code_format) /= 'onetep') Then
           Write (message,'(1x,2a)') Trim(error_dft),&
                                   & ' All values for U in &hubbard sub-block are set to zero! Please check'
           Call error_stop(message)
         End If
      End If
    End If

    ! kpoint sampling (optional)
    If (simulation_data%dft%kpoints%fread) Then
      If (simulation_data%dft%kpoints%fail) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Wrong (or missing) settings for "kpoints" directive (see manual)'
        Call error_stop(message)
      Else
        Do i =1, 3
          If (simulation_data%dft%kpoints%value(i) < 1) Then
            Write (message,'(1x,2a,i2,a)') Trim(error_dft), ' Number of kpoints along dimension ', i, &
                                         &' must be larger than zero!'
            Call error_stop(message)
          End If
        End Do 
        If (Trim(simulation_data%dft%kpoints%tag) /= 'mpack' .And. &
           Trim(simulation_data%dft%kpoints%tag) /= 'automatic' ) Then
           Write (messages(1),'(2a)')  Trim(error_dft), 'Invalid type of mesh for "kpoints". Options:'
           Write (messages(2),'(a)')   '- MPack (Monkhorst-Pack)'
           Write (messages(3),'(a)')   '- Automatic'
          Call info(messages, 3)
          Call error_stop(' ')
        End If
        ! Total number of kpoints
        simulation_data%dft%total_kpoints=1
        Do i = 1,3
          simulation_data%dft%total_kpoints= simulation_data%dft%total_kpoints * simulation_data%dft%kpoints%value(i)
        End Do
        error=.False.
        If (Trim(simulation_data%process)=='electrodeposition') Then
          If (Trim(simulation_data%normal_vector)=='c1') Then 
            If (simulation_data%dft%kpoints%value(1)/=1) Then
              error=.True.
              inorm=1
            End If
          Else If (Trim(simulation_data%normal_vector)=='c2') Then 
            If (simulation_data%dft%kpoints%value(2)/=1) Then
              error=.True.
              inorm=2
            End If
          Else If (Trim(simulation_data%normal_vector)=='c3') Then 
            If (simulation_data%dft%kpoints%value(3)/=1) Then
              error=.True.
              inorm=3
            End If
          End If
        End If
        If (error) Then
          Write (messages(1),'(2(1x,a),2a)') Trim(error_dft), 'The user has specified the that surface model for&
                                   & electrodeposition is perpendicular to the cell vector "',&
                                   & Trim(simulation_data%normal_vector), '".'
          Write (messages(2),'(1x,a,i2, a)') 'The assigned number of kpoints associated with this vector&
                                     & (in the reciprocal space) is ', simulation_data%dft%kpoints%value(inorm),&
                                     & ' but it must be equal to 1 !!! Please change the&
                                     & corresponding value in directive "kpoints".'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End If
    Else
      ! k-points set to Gamma 
      Do i=1, 3
        simulation_data%dft%kpoints%value(i)=1
      End Do
      simulation_data%dft%total_kpoints = 1
      simulation_data%dft%kpoints%tag = 'mpack'
    End If

  End Subroutine check_dft_settings

  Subroutine check_motion_settings(files, stoich_data, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check motion related directive for the generation of input
    ! files for atomistic level simulations
    !
    ! author    - i.scivetti May-June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(stoich_type),  Intent(In   ) :: stoich_data
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256)  :: messages(5)
    Character(Len=256)  :: message
    Character(Len=256)  :: error_motion
    Integer(Kind=wi)    :: i, j, k
    Logical             :: error, loop 

    Integer(Kind=wi)    :: fail
    Real(Kind=wp),     Allocatable :: species_mass(:)

    error_motion = '***ERROR in &motion_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
    
    ! Relaxation method (compulsory for relaxation) 
    If (simulation_data%motion%relax_method%fread) Then
      If (simulation_data%motion%relax_method%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong settings for "relax_method" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%simulation%type) == 'md') Then
          Write (message,'(2(1x,a))') Trim(error_motion),& 
                                    &'Specification of "relax_method" is not meaningful for MD settings. Please remove it.'
          Call error_stop(message)
        End If        
      End If
    Else
      If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'The user must specify directive "relax_method" for relaxation'
        Call error_stop(message)        
      End If
    End If

    ! Condition for masses
    If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      If (simulation_data%motion%mass_info%stat) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Definition of block &masses is not needed for geometry relaxation.&
                                   & Please remove it.'
        Call error_stop(message)        
      End If
    End If

    ! Ensemble (compulsory for MD)
    If (simulation_data%motion%ensemble%fread) Then
      If (simulation_data%motion%ensemble%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong settings for "ensemble" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%motion%ensemble%type) /= 'nve'  .And.&
           Trim(simulation_data%motion%ensemble%type) /= 'nvt'  .And.&
           Trim(simulation_data%motion%ensemble%type) /= 'npt'  .And.&
           Trim(simulation_data%motion%ensemble%type) /= 'nph'  ) Then
          Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                    &'Invalid specification for directive "ensemble". Options are:'
          Write (messages(2),'(1x,a)') '- NVE (Microcanonical ensemble)'
          Write (messages(3),'(1x,a)') '- NVT (Canonical ensemble)'
          Write (messages(4),'(1x,a)') '- NpT (Isothermal-Isobaric ensemble)'
          Write (messages(5),'(1x,a)') '- NpH (Isoenthalpic-Isobaric ensemble)'
          Call info(messages, 5)
          Call error_stop(' ')
        Else
          If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
             Write (message,'(2(1x,a))') Trim(error_motion), 'Directive "ensemble" is irrelevant for geometry relaxation.&
                                       & Please remove it'
             Call error_stop(message)
          End If
        End If
      End If
    Else
      If (Trim(simulation_data%simulation%type) == 'md') Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'The user must specify directive "ensemble" for MD'
        Call error_stop(message)
      End If
    End If

    ! Force tolerance (optional)
    If (simulation_data%motion%delta_f%fread) Then
      If (simulation_data%motion%delta_f%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong (or missing) settings for "force_tolerance" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%simulation%type) == 'md') Then
          Write (message,'(2(1x,a))') Trim(error_motion),& 
                                    &'Specification of "force_tolerance" is not meaningful for MD settings. Please remove it.'
          Call error_stop(message)
        End If        
        If (simulation_data%motion%delta_f%value(1) < epsilon(simulation_data%motion%delta_f%value(1))) Then
          Write (message,'(2(1x,a))') Trim(error_motion), &
                                    &'Input value for "force_tolerance" MUST be larger than zero'
          Call error_stop(message)
        End If
      End If
    End If

    ! Ionic steps (compulsory)
    If (simulation_data%motion%ion_steps%fread) Then
      If (simulation_data%motion%ion_steps%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong (or missing) settings for "ion_steps" directive'
        Call error_stop(message)
      Else
        If (simulation_data%motion%ion_steps%value <= 0) Then
          Write (message,'(1x,2a)') Trim(error_motion), ' Number of "ion_step" must be larger than 0 (at least)'
          Call error_stop(message)
        End If
      End If
    Else
      Write (message,'(2(1x,a))') Trim(error_motion), 'The user must specify directive "ion_steps"'
      Call error_stop(message)
    End If

    ! Timestep (optional)
    If (simulation_data%motion%timestep%fread) Then
      If (simulation_data%motion%timestep%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong (or missing) settings for "timestep" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
          Write (message,'(2(1x,a))') Trim(error_motion),& 
                                    &'Specification of "timestep" is not meaningful for "relaxation". Please remove it.'
          Call error_stop(message)
        End If        
        
        If (simulation_data%motion%timestep%value < epsilon(simulation_data%motion%timestep%value)) Then
          Write (message,'(2(1x,a))') Trim(error_motion), &
                                    &'Input value for "timestep" MUST be larger than zero'
          Call error_stop(message)
        End If
        If (Trim(simulation_data%motion%timestep%units) /= 'fs' .And. &
           Trim(simulation_data%motion%timestep%units) /= 'fsec') Then
           Write (message,'(2a)')  Trim(error_motion), 'Units for directive "timestep" must be "fs" or "fsec" (fento-seconds) '
          Call info(message, 1)
          Call error_stop(' ')
        End If
      End If
    End If

    ! Modify simulation cell volume (optional)
    If (simulation_data%motion%change_cell_volume%fread) Then
      If (simulation_data%motion%change_cell_volume%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong settings for "change_cell_volume" directive.&
                                 & Either .True. or .False.'
        Call error_stop(message)
      End If
    Else
      simulation_data%motion%change_cell_volume%stat=.False.
    End If

    ! Modify simulation cell shape (optional)
    If (simulation_data%motion%change_cell_shape%fread) Then
      If (simulation_data%motion%change_cell_shape%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong settings for "change_cell_shape" directive.&
                                 & Either .True. or .False.'
        Call error_stop(message)
      End If
    Else
      simulation_data%motion%change_cell_shape%stat=.False.
    End If

    ! Temperature (compulsory if MD)
    If (simulation_data%motion%temperature%fread) Then
      If (simulation_data%motion%temperature%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong (or missing) settings for "temperature" directive.&
                                & Both value and units are needed.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
          Write (message,'(2(1x,a))') Trim(error_motion),& 
                                    &'Specification of "temperature" is not meaningful for "relaxation". Please remove it.'
          Call error_stop(message)
        End If        
        If (simulation_data%motion%temperature%value <= 0.0_wp) Then
          Write (message,'(2(1x,a))') Trim(error_motion), &
                                    &'Input value for "temperature" MUST be larger than zero!!'
          Call error_stop(message)
        End If
        If (Trim(simulation_data%motion%temperature%units) /= 'k') Then
          Write (message,'(2(1x,a))') Trim(error_motion), &
                                    &'Wrong units for directive "temperature". Units must be in K'
          Call error_stop(message)
        End If
      End If
    Else
      If (Trim(simulation_data%simulation%type) == 'md') Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'The user must specify directive "temperature" for MD'
        Call error_stop(message)
      End If
    End If

    ! Thermostat (compulsory if ensemble is NVT of NPT)
    If (simulation_data%motion%thermostat%fread) Then
      If (simulation_data%motion%thermostat%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong (or missing) settings for "thermostat" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
          Write (message,'(2(1x,a))') Trim(error_motion),& 
                                    &'Specification of "thermostat" is not meaningful for "relaxation". Please remove it.'
          Call error_stop(message)
        End If        
        If (Trim(simulation_data%motion%ensemble%type) == 'nve' .Or. &
           Trim(simulation_data%motion%ensemble%type) == 'nph') Then
           Write (message,'(2(1x,a))') Trim(error_motion), 'Specification of "thermostat" is incompatible for the&
                                    & chosen ensemble. Please remove it and rerun'
           Call error_stop(message)
        End If
      End If
    Else
      If (Trim(simulation_data%simulation%type) == 'md') Then
        If (Trim(simulation_data%motion%ensemble%type) /= 'nve' .And. &
           Trim(simulation_data%motion%ensemble%type) /= 'nph')then
          Write (message,'(2(1x,a))') Trim(error_motion), 'The user must specify directive "thermostat" for MD'
          Call error_stop(message)        
        End If
      End If
    End If

    ! Relaxation time for the thermostat (compulsory if there is a thermostat)
    If (simulation_data%motion%relax_time_thermostat%fread) Then
      If (simulation_data%motion%relax_time_thermostat%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong (or missing) settings for&
                                  & "relax_time_thermostat" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
          Write (message,'(2(1x,a))') Trim(error_motion),& 
                                    &'Specification of "relax_time_thermostat" is not meaningful for "relaxation".&
                                    & Please remove it.'
          Call error_stop(message)
        End If        

        If (Trim(simulation_data%motion%ensemble%type) == 'nve' .Or. &
           Trim(simulation_data%motion%ensemble%type) == 'nph') Then
           Write (message,'(2(1x,a))') Trim(error_motion), 'Specification of "relax_time_thermostat" is incompatible for the&
                                    & chosen ensemble. Please remove it and rerun.'
           Call error_stop(message)
        End If
 
        If (simulation_data%motion%relax_time_thermostat%value < &
         & epsilon(simulation_data%motion%relax_time_thermostat%value)) Then
          Write (message,'(2(1x,a))') Trim(error_motion), &
                                    &'Input value for "relax_time_thermostat" MUST be larger than zero'
          Call error_stop(message)
        End If
        If (Trim(simulation_data%motion%relax_time_thermostat%units) /= 'fs' .And. &
           Trim(simulation_data%motion%relax_time_thermostat%units) /= 'fsec') Then
           Write (message,'(2a)')  Trim(error_motion), 'Units for directive "relax_time_thermostat" must be&
                                 & "fs" or "fsec" (fento-seconds) '
          Call info(message, 1)
          Call error_stop(' ')
        End If
      End If
    End If

    ! Pressure (compulsory if NPT or NHH)
    If (simulation_data%motion%pressure%fread) Then
      If (simulation_data%motion%pressure%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong (or missing) settings for "pressure" directive.&
                                & Both value and units are required.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%simulation%type) == 'md') Then
          If (Trim(simulation_data%motion%ensemble%type) == 'nve' .Or. &
             Trim(simulation_data%motion%ensemble%type) == 'nvt' ) Then
             Write (message,'(1x,4a)') Trim(error_motion),& 
                                      &' Specification of "pressure" is not meaningful for ensemble "', &
                                      & Trim(simulation_data%motion%ensemble%type), '". Please remove it.'
            Call error_stop(message)
          End If     
        End If     
        If (Trim(simulation_data%motion%pressure%units) /= 'kb' .And. &
           Trim(simulation_data%motion%pressure%units) /= 'kbar' ) Then
          Write (message,'(2(1x,a))') Trim(error_motion), &
                                    &'Invalid units for "pressure". Please check manual'
          Call error_stop(message)
        End If
      End If
    Else
      If (Trim(simulation_data%simulation%type) == 'md') Then
        If (Trim(simulation_data%motion%ensemble%type) == 'npt'  .Or. & 
           Trim(simulation_data%motion%ensemble%type) == 'nph') Then
           Write (message,'(1x,4a)') Trim(error_motion), ' The user must specify directive "pressure" for MD simulations&
                                    & with the "', Trim(simulation_data%motion%ensemble%type), '" ensemble.'
          Call error_stop(message)
        End If
      End If
      ! Assign pressure of 0 atm
      simulation_data%motion%pressure%value=0.0_wp
      simulation_data%motion%pressure%units='kb' 
    End If
       
    ! Barostat (no needed for VASP and CP2K)
    If (simulation_data%motion%barostat%fread) Then
      If (simulation_data%motion%barostat%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong (or missing) settings for "barostat" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
          Write (message,'(2(1x,a))') Trim(error_motion),& 
                                    &'Specification of "barostat" is not meaningful for "relaxation". Please remove it.'
          Call error_stop(message)
        End If        
        If (Trim(simulation_data%motion%ensemble%type) == 'nve' .Or. &
           Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
           Write (message,'(2(1x,a))') Trim(error_motion), 'Specification of "barostat" is incompatible for the&
                                    & selected ensemble. Please remove it and rerun'
           Call error_stop(message)
        Else
          If (Trim(simulation_data%code_format) == 'vasp' .Or. &
             Trim(simulation_data%code_format) == 'cp2k') Then
             Write (messages(1),'(2(1x,a))') Trim(error_motion),&
                                         & 'Specification of "barostat" is not needed for the selected code format.'
             Write (messages(2),'(1x,a)') 'Barostast specifications are determined from "ensemble" and/or "relax_time_barostat".&
                                         & Please remove it and rerun.'
             Call info(messages, 2)
             Call error_stop(' ')
          End If 
        End If
      End If
    Else
      If (Trim(simulation_data%simulation%type) == 'md') Then
        If (Trim(simulation_data%motion%ensemble%type) /= 'nve' .And. &
           Trim(simulation_data%motion%ensemble%type) /= 'nvt')then
          If (simulation_data%code_format /= 'vasp'   .And. &
              simulation_data%code_format /= 'onetep' .And. &
              simulation_data%code_format /= 'cp2k') Then
             Write (message,'(2(1x,a))') Trim(error_motion), 'The user must specify directive "barostat" for the selected emsemble'
             Call error_stop(message)        
          End If
        End If
      End If
    End If

    ! Relaxation time for the barostat (compulsory for NPT and NPH)
    If (simulation_data%motion%relax_time_barostat%fread) Then
      If (simulation_data%motion%relax_time_barostat%fail) Then
        Write (message,'(2(1x,a))') Trim(error_motion), 'Wrong (or missing) settings for&
                                  & "relax_time_barostat" directive.'
        Call error_stop(message)
      Else
        If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
          Write (message,'(2(1x,a))') Trim(error_motion),& 
                                    &'Specification of "relax_time_barostat" is not meaningful for "relaxation".&
                                    & Please remove it.'
          Call error_stop(message)
        End If        

        If (Trim(simulation_data%motion%ensemble%type) == 'nve' .Or. &
           Trim(simulation_data%motion%ensemble%type) == 'nvt' ) Then
           Write (message,'(2(1x,a))') Trim(error_motion), 'Specification of "relax_time_barostat" is incompatible for the&
                                    & selected ensemble. Please remove it and rerun'
           Call error_stop(message)
        End If

        If (simulation_data%motion%relax_time_barostat%value < epsilon(simulation_data%motion%relax_time_barostat%value)) Then
          Write (message,'(2(1x,a))') Trim(error_motion), &
                                    &'Input value for "relax_time_barostat" MUST be larger than zero'
          Call error_stop(message)
        End If
        If (Trim(simulation_data%motion%relax_time_barostat%units) /= 'fs' .And. &
           Trim(simulation_data%motion%relax_time_barostat%units) /= 'fsec') Then
           Write (message,'(2a)')  Trim(error_motion), 'Units for directive "relax_time_barostat" must be&
                                 & "fs" or "fsec" (fento-seconds) '
          Call info(message, 1)
          Call error_stop(' ')
        End If
      End If
    End If

    ! About simulation cell  
    If (Trim(simulation_data%simulation%type) == 'md') Then
      If (Trim(simulation_data%motion%ensemble%type) == 'nve' .Or. Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
        If (simulation_data%motion%change_cell_volume%stat .Or. simulation_data%motion%change_cell_shape%stat) Then
          Write (message,'(1x,4a)') Trim(error_motion), ' Ensemble "', Trim(simulation_data%motion%ensemble%type), '" must keep the&
                                 & simulation cell fixed. Please set "change_cell_volume" and "change_cell_shape" for&
                                 & .False. or remove them.' 
          Call error_stop(message)        
        End If
      ElseIf (Trim(simulation_data%motion%ensemble%type) == 'npt' .Or.&
            & Trim(simulation_data%motion%ensemble%type) == 'nph') Then
        If ( (.Not. simulation_data%motion%change_cell_volume%stat) .Or. &
           & (.Not. simulation_data%motion%change_cell_shape%stat)) Then 
          Write (message,'(1x,4a)') Trim(error_motion), ' For ensemble "', Trim(simulation_data%motion%ensemble%type), & 
                                 & '" the volume and shape of the simulation cell must be allowed to change.&
                                 & Set both directives "change_cell_volume" and "change_cell_shape" to .True. and rerun.' 
          Call error_stop(message)        
        End If
      End If
    Else If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      If ( (.Not. simulation_data%motion%change_cell_volume%stat) .And. (.Not. simulation_data%motion%change_cell_shape%stat)) Then 
        If (simulation_data%motion%pressure%fread) Then
          Write (messages(1),'(1x,2a)') Trim(error_motion), ' Definition of "pressure" for a system with fixed volume and shape is&
                                      & meaningless.' 
          Write (messages(2),'(1x,a)') 'Set one or both directives "change_cell_volume" and "change_cell_shape"&
                                      & to .True. and rerun.' 
          Call info(messages, 2)
          Call error_stop(' ')        
        End If 
      End If
    End If

    ! Check settings of block &masses
    If (simulation_data%motion%mass_info%stat) Then
      ! Check if user has included all the tags
      Do i=1, simulation_data%total_tags
        error=.True.
        Do j=1, simulation_data%total_tags
          If (Trim(simulation_data%motion%mass(i)%tag)==Trim(simulation_data%component(j)%tag)) Then
            simulation_data%motion%mass(i)%element=simulation_data%component(j)%element
            error=.False.
          End If
        End Do
        If (error) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_motion),&
                                   & 'Atomic tag "', Trim(simulation_data%motion%mass(i)%tag), '" declared in&
                                   & "&masses" has not been defined in "&block_input_composition". Please check'
          Call error_stop(message)
        End If
      End Do

      ! Check there is no negative not zero value for mass 
      Do i=1, simulation_data%total_tags
        If (simulation_data%motion%mass(i)%value < 0.0_wp .Or. &
            Abs(simulation_data%motion%mass(i)%value) < epsilon(1.0_wp)) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_motion),&
                                   & 'Value of mass assigned to atomic tag "', Trim(simulation_data%motion%mass(i)%tag),&
                                   & '" must be positive and different from zero (check block &masses)'
          Call error_stop(message)
        End If
      End Do
      ! Corroborate if masses defined in &block_species are compatible with the information in &masses
      Allocate(species_mass(stoich_data%num_species%value), Stat=fail)
      If (fail > 0) Then
        Write (message,'(1x,1a)') '***ERROR: unsuccessful allocations of arrays in check_motion_settings'
        Call error_stop(message)
      End If
     
      Do i=1, stoich_data%num_species%value
        species_mass(i)=0.0_wp
        Do j=1, stoich_data%species(i)%num_components 
          k=1
          loop=.True.
          Do While (k <= simulation_data%total_tags .And. loop)
            If (Trim(simulation_data%motion%mass(k)%tag) == Trim(stoich_data%species(i)%component%tag(j))) Then
              species_mass(i)=species_mass(i)+stoich_data%species(i)%component%N0(j)*simulation_data%motion%mass(k)%value
              loop=.False.
            End If
            k=k+1
          End Do
        End Do
          If(Abs(species_mass(i)-stoich_data%species(i)%mass) > species_mass_error)Then
            Write (message,'(1x,a,1x,3a)') Trim(error_motion),&
                         & 'Inconsistencies between the mass of species "', Trim(stoich_data%species(i)%tag),&
                         & '" defined in &block_species and the masses for its constituents (see &masses and&
                         & &Block_Species_Components).'
            Call error_stop(message)
          End If
      End Do 

      Deallocate(species_mass)

    End If

  End Subroutine check_motion_settings

  Subroutine warning_simulation_settings(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about: 
    ! 1) the need of any possible further check to fully set the simulation.
    ! 2) problems the user may find when executing the jobs for the requested
    !    choice of settings
    !  
    ! A message is printed to OUT_EQCM file, following the generation of 
    ! files for simulation. So far, the following codes are implementend
    ! - VASP (DFT)
    ! - CP2K (DFT)
    !
    ! author    - i.scivetti May-June 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(7), website, code

    If (Trim(simulation_data%code_format)=='vasp') Then
      code='VASP'
      website=Trim(web_vasp)
    Else If (Trim(simulation_data%code_format)=='cp2k') Then
      code='CP2K'
      website=Trim(web_cp2k)
    Else If (Trim(simulation_data%code_format)=='castep') Then
      code='CASTEP'
      website=Trim(web_castep)
    Else If (Trim(simulation_data%code_format)=='onetep') Then
      code='ONETEP'
      website=Trim(web_onetep)
    End If

    Call info(' ', 1)
    Write (messages(1), '(1x,a)')  '*******************************************************************'
    Write (messages(2), '(1x,2a)') 'Aspects to take into consideration for simulations with ', Trim(code)
    Write (messages(3), '(1x,a)')  '*******************************************************************'
    Write (messages(4), '(1x,a)')  'Specification of directives for the generated files might need to be adjusted&
                                 & depending on the system under study.'
    Write (messages(5), '(1x,3a)') 'Implementation of directives has been tested and validated using version ',&
                              & Trim(simulation_data%code_version), ' of the code.'
    Write (messages(6), '(1x,3a)')  'The user is responsible to check if the defined settings are compatible with&
                                 & the version of the ', Trim(code),' code used for the simulations.'
    Write (messages(7), '(1x,3a)') 'For more information visit ', Trim(website), ', and read the corresponding section&
                                  & of the ALC_EQCM manual.'  
    Call info(messages, 7)


    If (.Not. simulation_data%dft%pp_info%stat) Then
      Call info(' ', 1)
      Call info(' *** WARNING *************************************************************************', 1)
      If (Trim(simulation_data%code_format) /= 'cp2k') Then
        Write (messages(1), '(1x,a)') '*** If pseudopotential files have been already generated, please ignore this message:'
        Call info(messages, 1)
      End If
      Write (messages(1), '(1x,a)') '    ALC_EQCM has NOT generated pseudopotential files because the &pseudo_potentials' 
      If (Trim(simulation_data%code_format) /= 'castep') Then
        Write (messages(2), '(1x,a)') '    sub-block has not been defined. Thus, the generated files are NOT SUFFICIENT to'
        Write (messages(3), '(1x,a)') '    execute the DFT simulations. The user MUST set the &pseudo_potentials sub-block within'
        Write (messages(4), '(1x,a)') '    &dft_settings (and relevant files within the DFT folder) for automatic generation'
        Write (messages(5), '(1x,a)') '    of pseudopotentials. It is NOT RECOMMENDED to set pseudopotential info manually'
      Else
        Write (messages(2), '(1x,a)') '    sub-block has not been defined. Still, the generated files are sufficient to'
        Write (messages(3), '(1x,a)') '    execute the DFT simulations since CASTEP will generate ultra-soft pseudopotentials'
        Write (messages(4), '(1x,a)') '    "on-the-fly". If other type of pseudopotentials is needed, the user must set the'
        Write (messages(5), '(1x,a)') '    &pseudo_potentials sub-block for automatic generation of files'  
      End If        
      Call info(messages, 5)
      Call info(' *************************************************************************************', 1)
    End If
     
    If (.Not. simulation_data%dft%basis_info%stat) Then
      If (Trim(simulation_data%code_format) == 'cp2k') Then
        Call info(' ', 1)
        Call info(' *** WARNING *************************************************************************', 1)
        Write (messages(1), '(1x,a)') '    ALC_EQCM has NOT generated basis set files because the &basis_set sub-block' 
        Write (messages(2), '(1x,a)') '    has not been defined. Thus, the generated files are NOT SUFFICIENT to'
        Write (messages(3), '(1x,a)') '    execute the DFT simulations. The user MUST set the &basis_set sub-block'
        Write (messages(4), '(1x,a)') '    within &dft_settings for automatic generation of basis set input files.'
        Write (messages(5), '(1x,a)') '    It is NOT RECOMMENDED to set information for the basis sets manually'
        Call info(messages, 5)
        Call info(' *************************************************************************************', 1)
      End If        
    End If

    If (simulation_data%extra_info%fread) Then
      Call info(' ', 1)
      Call info(' *** WARNING *************************************************************************', 1)
      Write (messages(1), '(1x,a)') '  - ALC_EQCM only checked that the information added in "&extra_directives" has not'
      Write (messages(2), '(1x,a)') '    been already defined from the settings provided in "&block_simulation_settings"'   
      Write (messages(3), '(1x,a)') '  - specification in "&extra_directives" for functionalities that are not implemented'
      Write (messages(4), '(1x,a)') '    in ALC_EQCM were printed, but their correctness is full responsibility of the user'
      Call info(messages, 4)
      Call info(' *************************************************************************************', 1)
    End If

    If (Trim(simulation_data%code_format)=='cp2k') Then
      Call info(' ', 1)
      Call info(' *** WARNING ****************************************************************************', 1)
      Write (messages(1), '(1x,a)') '   Due to the complex block structure of the input.cp2k file, the use of'
      Write (messages(2), '(1x,a)') '   "&extra_directives" is not allowed for the generation of files for CP2K simulations.' 
      Write (messages(3), '(1x,a)') '   Although most directives are defined from &block_simulation_settings, there are'   
      Write (messages(4), '(1x,a)') '   keywords that have been set arbitrarily (see ALC_EQCM manual) based on our previous'
      Write (messages(5), '(1x,a)') '   experience with CP2K. Unfortunately, changes to such directives must be set manually.' 
      Call info(messages, 5)
      Call info(' ****************************************************************************************', 1)
    End If

    If (Trim(simulation_data%code_format)=='vasp') Then
      Call warnings_vasp(simulation_data)
    Else If (Trim(simulation_data%code_format)=='cp2k') Then
      Call warnings_cp2k(simulation_data)
    Else If (Trim(simulation_data%code_format)=='castep') Then
      Call warnings_castep(simulation_data)
    Else If (Trim(simulation_data%code_format)=='onetep') Then
      Call warnings_onetep(simulation_data)
    End If 

  End Subroutine warning_simulation_settings

  Subroutine summary_dft_settings(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to summarise DFT settings from the information 
    ! provided by the user 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),  Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(5)
    Logical :: loop
    Integer(Kind=wi) :: i

    If (Trim(simulation_data%code_format)=='vasp') Then
      Write (messages(1),'(1x,a)') 'Input files have been generated for execution with the VASP code.'
      Write (messages(3),'(1x,a)') 'POSCAR is a copy of the generated atomic structure (SAMPLE.vasp).'  
      If (simulation_data%dft%pp_info%stat) Then
        Write (messages(2),'(1x,a)') 'Files INCAR, KPOINTS, POTCAR and POSCAR are located at each sub-folder.'  
        Call info(messages, 3)
      Else
        Write (messages(2),'(1x,a)') 'Files INCAR, KPOINTS and POSCAR are located at each sub-folder.'  
        Write (messages(4),'(1x,a)') 'WARNING: Pseudopotential file (POTCAR) has not been set due to&
                                    & the missing sub-block &pseudo_potentials.'
        Call info(messages, 4)
      End If         
    Else If (Trim(simulation_data%code_format)=='cp2k') Then    
      Write (messages(1),'(1x,a)') 'Input files have been generated for execution with the CP2K code.'
      Write (messages(2),'(1x,a)') 'File input.cp2k with DFT settings is located at each sub-folder.'  
      If (simulation_data%dft%pp_info%stat) Then
        Write (messages(3),'(1x,2a)') 'Pseudopotentials for the atomic species have also been copied to each sub-folder.&
                                    & See file ', Trim(simulation_data%dft%pseudo_pot(1)%file_name) 
      Else
        Write (messages(3),'(1x,a)') 'WARNING: Pseudopotentials files have not been set due to the missing&
                                    & sub-block &pseudo_potentials.'
      End If         
      Call info(messages, 3)
    Else If (Trim(simulation_data%code_format)=='onetep') Then    
      Write (messages(1),'(1x,a)') 'Input files have been generated for execution with the ONETEP code.'
      Write (messages(2),'(1x,a)') 'File model.dat with DFT settings is located at each sub-folder.'  
      If (simulation_data%dft%pp_info%stat) Then
        Write (messages(3),'(1x,a)') 'Pseudopotentials for the atomic species have also been copied to each sub-folder.'  
      Else
        Write (messages(3),'(1x,a)') 'WARNING: Pseudopotentials files have not been set due to the missing&
                                    & sub-block &pseudo_potentials.'
      End If         
      Call info(messages, 3)
    Else If (Trim(simulation_data%code_format)=='castep') Then    
      Write (messages(1),'(1x,a)') 'Input files have been generated for execution with the CASTEP code.'
      Write (messages(2),'(1x,a)') 'Files model.param and model.cell with DFT settings are located at each sub-folder.'  
      If (simulation_data%dft%pp_info%stat) Then
        Write (messages(3),'(1x,a)') 'Pseudopotentials for the atomic species have also been copied to each sub-folder.'  
        Call info(messages, 3)
      Else
        Write (messages(3),'(1x,a)') 'WARNING: Pseudopotentials files have not been set due to the missing&
                                    & sub-block &pseudo_potentials.'
      End If         
      Call info(messages, 3)
    End If

    Call info(' === DFT settings:', 1)
    ! Grand Canonical
    If (simulation_data%dft%gc%activate%fread) Then
       Write (messages(1),'(1x,a)')  '- Grand-Canonical approximation for electrochemical conditions'
       Call info(messages, 1)
    End If        

    ! Spin polarization 
    If (simulation_data%dft%spin_polarised%fread) Then
       Write (messages(1),'(1x,a)')         '- spin-polarised calculation'
    Else
       Write (messages(1),'(1x,a)')         '- non-spin-polarised calculation'
    End If
    Call info(messages,1)

    ! XC and vdW
     Write (messages(1),'(1x,2a)') '- XC level:   ', Trim(simulation_data%dft%xc_level%type)
     Write (messages(2),'(1x,2a)') '- XC version: ', Trim(simulation_data%dft%xc_ref)
     Call info(messages,2)
     If (simulation_data%dft%vdw%fread) Then
       Write (messages(1),'(1x,2a)') '- vdW corrections: ', Trim(simulation_data%dft%vdw_ref)
     Else
       Write (messages(1),'(1x,a)')  '- vdW corrections are NOT included'         
     End If
     Call info(messages,1)

    ! vdW kernel
    If (simulation_data%dft%need_vdw_kernel) Then
      Write (messages(1),'(1x,3a)') '- each sub-folder also contains the supporting file "',&
                                 & Trim(simulation_data%dft%vdw_kernel_file), '" needed to compute vdW corrections'
      Call info(messages,1)
    End If

    ! Max SCF steps 
    Write (messages(1),'(1x,a,i4)')      '- maximum SCF steps for electronic convergence: ', simulation_data%dft%scf_steps%value
    ! Energy cutoff
    Write (messages(2),'(1x,a,f8.2,1x,a)') '- energy cutoff: ', simulation_data%dft%encut%value, &
                                            & Trim(simulation_data%dft%encut%units) 
    ! Energy tolerance 
    Write (messages(3),'(1x,a,e12.4,1x,a)') '- energy tolerance: ', simulation_data%dft%delta_e%value, &
                                      & simulation_data%dft%delta_e%units
 
    Call info(messages,3)
 
    If (Trim(simulation_data%code_format)=='vasp') Then
      If (simulation_data%dft%mixing%fread) Then      
        Write (messages(1),'(1x,2a)')     '- mixing scheme: ', Trim(simulation_data%dft%mixing%type)
      Else 
        Write (messages(1),'(1x,a)')      '- Pulay mixing scheme is set by default'      
      End If  
      Call info(messages,1)
      Write (messages(1),'(1x,2a)')       '- smearing method: ', Trim(simulation_data%dft%smear%type)
      Call info(messages,1)
      If (Trim(simulation_data%dft%smear%type) /= 'tetrahedron'  ) Then
        Write (messages(1),'(1x,a,f6.3,1x,a)')  '- smearing width: ', simulation_data%dft%width_smear%value, &
                                          &  Trim(simulation_data%dft%width_smear%units)
        Call info(messages,1)
      End If
    Else If (Trim(simulation_data%code_format)=='cp2k') Then
      If (simulation_data%dft%ot%stat) Then
        Write (messages(1),'(1x,a)') '- Orbital Transformation (OT) method for wavefunction optimisation' 
        Call info(messages,1)
      Else
        Write (messages(1),'(1x,a)')             '- Standard LAPACK diagonalization for the KS equations' 
        Write (messages(2),'(1x,2a)')            '- mixing scheme: ', Trim(simulation_data%dft%mixing%type)
        Write (messages(3),'(1x,2a)')            '- smearing method: ', Trim(simulation_data%dft%smear%type)
        Write (messages(4),'(1x,a,f12.3,1x,a)')  '- smearing width: ', simulation_data%dft%width_smear%value, &
                                                  &  Trim(simulation_data%dft%width_smear%units)
        Call info(messages,4)
      End If
    Else If (Trim(simulation_data%code_format)=='castep') Then  
      If (Trim(simulation_data%dft%smear%type) == 'fix_occupancy' ) Then
        Write (messages(1),'(1x,a)') '- fix occupancy for the electronic states'
        Call info(messages,1)
      Else
        If (simulation_data%dft%edft%fread) Then
          Write (messages(1),'(1x,a,f6.2,1x,a)')  '- Ensemble DFT with an electronic temperature of ', &
                                                 & simulation_data%dft%width_smear%value, '(in eV)' 
        Else
          If (simulation_data%dft%mixing%fread) Then 
            If (Trim(simulation_data%dft%mixing%type)    == 'broyden-2nd') Then 
              Write (messages(1),'(1x,a)')      '- mixing scheme:   Broyden-2nd (it is named as Broyden)'
            Else        
              Write (messages(1),'(1x,2a)')     '- mixing scheme: ', Trim(simulation_data%dft%mixing%type)
            End If
          Else 
            Write (messages(1),'(1x,a)')      '- Broyden-2nd mixing scheme is set by default (named as Broyden in CASTEP)'      
          End If  
        End If
        Write (messages(2),'(1x,2a)')           '- smearing method: ', Trim(simulation_data%dft%smear%type)
        Write (messages(3),'(1x,a,f6.3,1x,a)')  '- smearing width: ', simulation_data%dft%width_smear%value, &
                                  &  Trim(simulation_data%dft%width_smear%units)
        Call info(messages,3)
      End If
    Else If (Trim(simulation_data%code_format)=='onetep') Then 
      If (simulation_data%dft%edft%fread) Then
         Write (messages(1),'(1x,a)')            '- Ensemble DFT computation using the Fermi-Dirac distribution'
         Write (messages(2),'(1x,2a)')           '- mixing scheme: ', Trim(simulation_data%dft%mixing%type)
         Write (messages(3),'(1x,a,f6.3,1x,a)')  '- smearing width: ', simulation_data%dft%width_smear%value, &
                                                 &  Trim(simulation_data%dft%width_smear%units)
        Call info(messages,3)          
      Else
         Write (messages(1),'(1x,a)')        '- standard ONETEP computation for systems with band gap'
         Write (messages(2),'(1x,2a)')       '- mixing scheme: ', Trim(simulation_data%dft%mixing%type)
         Call info(messages, 2)
      End If
    End If

    If (simulation_data%dft%basis_info%fread) Then  
      If (Trim(simulation_data%code_format)=='cp2k') Then
         Write (messages(1),'(1x,a)')        '- the atomic basis set for the participating species is defined&
                                             & in sub-block &basis_set'
         Call info(messages,1)
      End If 
    End If 
 
    ! number of bands 
    If (simulation_data%dft%bands%fread) Then
      If (Trim(simulation_data%code_format)=='vasp') Then
        Write (messages(1),'(1x,a,i5)')        '- total number of bands is ', simulation_data%dft%bands%value
      Else 
        Write (messages(1),'(1x,a,i5)')        '- number of extra bands is ', simulation_data%dft%bands%value
      End If
      Call info(messages,1)
    End If

    ! precision (VASP)
    If (Trim(simulation_data%code_format)=='vasp') Then
      Write (messages(1),'(1x,2a)') '- precision:   ', Trim(simulation_data%dft%precision%type)
      Call info(messages,1)
    End If

    ! k-points info 
    If (simulation_data%dft%total_kpoints>1) Then
      Write (messages(1),'(1x,a,i3,3a)') '- a set of ', simulation_data%dft%total_kpoints, ' k-points to sample the reciprocal&
                                      & space using the "', Trim(simulation_data%dft%kpoints%tag), '" scheme (see manual)'
    Else
      Write (messages(1),'(1x,a)')        '- only the Gamma point is used for the reciprocal space' 
    End If
    Call info(messages,1)

    ! Initial magnetization 
    If (simulation_data%dft%mag_info%fread) Then
      Write (messages(1),'(1x,a)')       '- an initial magnetization is assigned to each atomic site'
      Call info(messages,1)
    End If

    ! Total magnetization    
    If (simulation_data%dft%total_magnetization%fread) Then
      Write (messages(1),'(1x,a,f8.3)')  '- the total magnetization is set to ', simulation_data%dft%total_magnetization%value
      Call info(messages,1)
    End If

    ! Hubbard corrections
    If (simulation_data%dft%hubbard_info%fread) Then
      If (.Not. simulation_data%dft%hubbard_all_U_zero) Then
        If (Trim(simulation_data%code_format)=='vasp' .Or. Trim(simulation_data%code_format)=='onetep') Then
          loop=.False.
          Do i=1, simulation_data%total_tags
            If (Abs(simulation_data%dft%hubbard(i)%J)>epsilon(1.0_wp)) Then
              loop=.True.
            End If
          End Do
          If (loop) Then
            Write (messages(1),'(1x,a)')  '- anisotropic (U-J) Hubbard corrections are imposed to correct&
                                         & for the occupancy of selected atomic sites'
          Else
            Write (messages(1),'(1x,a)')  '- isotropic (U-J) Hubbard corrections are imposed to correct&
                                         & for the occupancy of selected atomic sites'
          End If
          Call info(messages,1)
        Else
           Write (messages(1),'(1x,a)')   '- isotropic (U-J) Hubbard corrections are imposed to correct&
                                         & for the occupancy of selected atomic sites'
           Call info(messages,1)
        End If
      End If
    End If

  End Subroutine summary_dft_settings

  Subroutine summary_motion_settings(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to summarise motion settings from the information 
    ! provided by the user 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),  Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(5)

    If (Trim(simulation_data%simulation%type) == 'relax_geometry' .Or.&
    Trim(simulation_data%simulation%type) == 'md'         ) Then
      Call info(' === Ion-related settings:', 1)
      Write (messages(1),'(1x,a,i5)')     '- number of ionic steps: ', simulation_data%motion%ion_steps%value
      If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
        Write (messages(2),'(1x,2a)')    '- relaxation method: ', Trim(simulation_data%motion%relax_method%type)
        Write (messages(3),'(1x,a,f8.5,2(1x,a))')  '- force tolerance: ', simulation_data%motion%delta_f%value(1), &
                                           & Trim(simulation_data%motion%delta_f%units(1)),&
                                           & Trim(simulation_data%motion%delta_f%units(2))
        If ((.Not. simulation_data%motion%change_cell_volume%stat) .And. (.Not. simulation_data%motion%change_cell_shape%stat)) Then
          Write (messages(4),'(1x,a)')     '- the supercell is kept fixed'
        Else If (simulation_data%motion%change_cell_volume%stat .And. (.Not. simulation_data%motion%change_cell_shape%stat)) Then
          Write (messages(4),'(1x,a)')     '- the supercell volume is allowed to change but the shape is kept fixed'
        Else If ((.Not. simulation_data%motion%change_cell_volume%stat) .And. simulation_data%motion%change_cell_shape%stat) Then
          Write (messages(4),'(1x,a)')     '- the supercell shape is allowed to change but the volume is kept fixed'
        Else If (simulation_data%motion%change_cell_shape%stat .And. simulation_data%motion%change_cell_shape%stat) Then
          Write (messages(4),'(1x,a)')     '- both volume and shape of supercell are allowed to change'
        End If
        Call info(messages,3)
      Else If (Trim(simulation_data%simulation%type) == 'md') Then
        Write (messages(2),'(1x,2a)') '- MD ensemble: ', Trim(simulation_data%motion%ensemble%type)
        Write (messages(3),'(1x,a,f6.2,1x,a)') '- timestep: ', simulation_data%motion%timestep%value, &
                                                & Trim(simulation_data%motion%timestep%units)
        Call info(messages, 3)
        If (simulation_data%motion%mass_info%stat) Then
          Write (messages(1),'(1x,a)') '- masses for atomic species are defined in sub-block &masses' 
          Call info(messages, 1)
        End If
        If (simulation_data%motion%temperature%fread) Then
          Write (messages(1),'(1x,a,f6.2,1x,a)') '- temperature: ', simulation_data%motion%temperature%value,&
                                                & Trim(simulation_data%motion%temperature%units)  
          Call info(messages, 1)
        End If
        If (simulation_data%motion%thermostat%fread) Then
          Write (messages(1),'(1x,2a)') '- thermostat type: ', Trim(simulation_data%motion%thermostat%type)
          Call info(messages, 1)
        End If
        If (simulation_data%motion%relax_time_thermostat%fread) Then
          Write (messages(1),'(1x,a,f6.2,1x,a)') '- relaxation time for thermostat: ',&
                                                & simulation_data%motion%relax_time_thermostat%value,&
                                                & Trim(simulation_data%motion%relax_time_thermostat%units)  
          Call info(messages, 1)
        End If
        If (simulation_data%motion%barostat%fread) Then
          Write (messages(1),'(1x,2a)') '- barostat type: ', Trim(simulation_data%motion%barostat%type)
          Call info(messages, 1)
        End If
        If (simulation_data%motion%relax_time_barostat%fread) Then
          Write (messages(1),'(1x,a,f10.2,1x,a)') '- relaxation time for barostat: ', &
                                                & simulation_data%motion%relax_time_barostat%value,&
                                                & Trim(simulation_data%motion%relax_time_barostat%units)  
          Call info(messages, 1)
        End If
      End If 
      If (simulation_data%motion%pressure%fread) Then
        Write (messages(1),'(1x,a,f16.2,1x,a)') '- external pressure: ', simulation_data%motion%pressure%value,&
                                              & Trim(simulation_data%motion%pressure%units)  
        Call info(messages, 1)
      End If
    End If 

    Call info(' ', 1)

  End Subroutine summary_motion_settings

  Subroutine summary_simulation_settings(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to summarise simulation settings 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),  Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(2)

    Call info(' === General settings:', 1)
    ! Type of simulation
    If (Trim(simulation_data%simulation%type) == 'md') Then
      Write (messages(1),'(1x,a)')  '- type of simulation requested: Molecular Dynamics (MD)'
    Else If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      Write (messages(1),'(1x,a)')  '- type of simulation requested: geometry relaxation'
    Else If (Trim(simulation_data%simulation%type) == 'singlepoint') Then
      Write (messages(1),'(1x,a)')  '- type of simulation set: single point'
    End If
    Call info(messages, 1)
    ! Level of theory for interatomic interactions
    If (Trim(simulation_data%theory_level%type) == 'dft') Then
      Write (messages(1),'(1x,a)')  '- level of theory for interatomic interactions: DFT'
    End If
    Call info(messages, 1)

    ! Total system charge
    If (Abs(simulation_data%net_charge%value)>epsilon(1.0_wp)) Then
      If (Trim(simulation_data%code_format)/='cp2k' .And. Trim(simulation_data%code_format)/='onetep') Then
         Write (messages(1),'(1x,a,f8.3)') '- net charge: ', simulation_data%net_charge%value 
      Else
         Write (messages(1),'(1x,a,i4)')   '- net charge: ', Nint(simulation_data%net_charge%value) 
      End If
    Else
      Write (messages(1),'(1x,a)')         '- neutral supercell'
    End If
    Call info(messages,1)

    ! Print summary
    Call summary_dft_settings(simulation_data)
    Call summary_motion_settings(simulation_data)

  End Subroutine summary_simulation_settings

End Module simulation_setup
