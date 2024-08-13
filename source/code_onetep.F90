!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to check, define and print directives for simulations
! with ONETEP. This module also warns the user about aspects to take
! into consideration when performing simulations
!
! Copyright - Ada Lovelace Centre
!
! Author    - i.scivetti Oct 2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module code_onetep
  
  Use constants,        Only : date_RELEASE, &
                               max_components

  Use fileset,          Only : file_type, &
                               FILE_SET_EQCM, &
                               FILE_SET_SIMULATION,&
                               FOLDER_DFT 
  Use numprec,          Only : wi, &
                               wp
  Use process_data,     Only : capital_to_lower_case, remove_symbols
  Use references,       Only : bib_AVV10s, bib_blyp, bib_dftd2, bib_optb88, bib_optpbe, bib_pbe, bib_pbesol, &
                               bib_pw91, bib_pw92, bib_pz, bib_revpbe, bib_rp, bib_vdwdf, bib_vdwdf2, bib_vv10, &
                               bib_vwn, bib_xlyp, bib_fg, bib_fisicaro, bib_andreussi, bib_pbeq_onetep, &
                               bib_neutral_onetep, bib_gcdft_onetep 
  Use simulation_setup, Only : simul_type,&
                               max_directives, &
                               type_extra,&
                               type_ref_data
  Use simulation_tools, Only : check_extra_directives, &
                               check_initial_magnetization, &
                               check_settings_set_extra_directives, &
                               check_settings_single_extra_directive,&
                               print_extra_directives, &
                               print_warnings, &
                               record_directive, &
                               scan_extra_directive, &
                               set_reference_database                               
  Use stoichiometry,    Only : stoich_type
  Use unit_output,      Only : error_stop,&
                               info 

  Implicit None
  Private

  Public :: define_onetep_settings, print_onetep_settings, advise_onetep
  Public :: summary_solvation_onetep, summary_electrolyte_onetep 

Contains

  Subroutine define_onetep_settings(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define (and check) settings for ONETEP directives (set defaults values)
    !
    ! author    - i.scivetti July  2021
    ! refact    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data
    
    Call define_onetep_dft(files, simulation_data)
    Call define_onetep_motion(files, simulation_data)
    ! Solvation	
    If (simulation_data%solvation%info%stat) Then
      Call define_solvation_onetep(files, simulation_data)
      ! Boltzman ions
      If (simulation_data%electrolyte%info%stat) Then
        Call define_electrolyte_onetep(files, simulation_data)
      End If
    End If  
    
  End Subroutine define_onetep_settings

  Subroutine print_onetep_settings(files, net_elements, list_tag, list_number, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print settings from ONETEP directives
    !
    ! author    - i.scivetti May   2021
    ! refac     - i.scivetti Sept  2022 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),  Intent(InOut) :: files(:)
    Integer(Kind=wi), Intent(In   ) :: net_elements
    Character(Len=8), Intent(In   ) :: list_tag(max_components) 
    Integer(Kind=wi), Intent(In   ) :: list_number(max_components)
    Type(simul_type), Intent(InOut) :: simulation_data

    Integer(Kind=wi)   :: iunit
    Character(Len=256) :: message

    Integer(Kind=wi)   :: ic

    ic=1
    ! Open FILE_SET_SIMULATION file
    Open(Newunit=files(FILE_SET_SIMULATION)%unit_no, File=files(FILE_SET_SIMULATION)%filename,Status='Replace')
    iunit=files(FILE_SET_SIMULATION)%unit_no 

    Write (iunit,'(a)')  '##############################'
    Write (iunit,'(a)')  '# File generated with ALC_EQCM'
    Write (iunit,'(a)')  '##############################'
    Write (iunit,'(a)') ' '
 
    Write (iunit,'(a)') '##### Type of calculation'
    If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      Write (message,'(a)') 'task :  geometryoptimization' 
    Else If (Trim(simulation_data%simulation%type) == 'md') Then
      Write (message,'(a)') 'task :  moleculardynamics' 
    Else If (Trim(simulation_data%simulation%type) == 'singlepoint') Then
      Write (message,'(a)') 'task :  singlepoint'
    End If
    Call record_directive(iunit, message, 'task', simulation_data%set_directives%array(ic), ic)

    Call print_onetep_dft(iunit, ic, net_elements, list_tag, list_number, simulation_data)
    Call print_onetep_motion(iunit, ic, simulation_data)

    If (simulation_data%solvation%info%stat) Then
      Call print_onetep_solvation(iunit, ic, simulation_data)
    End If

    If (simulation_data%electrolyte%info%stat) Then
      Call print_onetep_electrolyte(iunit, ic, simulation_data)
    End If
    
   ! Total number of set directives
   simulation_data%set_directives%N0=ic-1

   If (simulation_data%extra_info%stat) Then
     Call print_extra_directives(iunit, simulation_data%extra_directives, simulation_data%set_directives, &
                            & simulation_data%code_format, simulation_data%simulation%type)
   End If

   Write (iunit,'(a)') ' '
   Write (iunit,'(a)') '#### Simulation cell and atomic coordinates'
   Write (iunit,'(a)') '#=========================================='

   Close(iunit)

  End Subroutine print_onetep_settings
    
  Subroutine advise_onetep(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about ONETEP settings 
    !
    ! author    - i.scivetti July 2021
    ! refact    - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    ! DFT
    Call advise_dft_onetep(simulation_data)
    ! Motion
    Call advise_motion_onetep(simulation_data)

    If (simulation_data%solvation%info%stat) Then
      Call advise_solvation_onetep(simulation_data)

      If (simulation_data%electrolyte%info%stat) Then
        Call advise_electrolyte_onetep(simulation_data)
      End If

      Call advise_multigrid_onetep(simulation_data)
      
    End If
    
    Call warnings_onetep_dft(simulation_data)
    
  End Subroutine advise_onetep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! DFT  
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine define_onetep_dft(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define (and check) DFT settings for ONETEP directives (set defaults values)
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256) :: message, messages(15)
    Character(Len=256) :: error_dft, error_gcdft
    Character(Len=256) :: path

    Integer(Kind=wi)   :: i, ifile 
  
    simulation_data%dft%onetep_paw = .False. 
    error_dft    = '***ERROR in &dft_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
    error_gcdft  = '***ERROR in &gcdft (file '//Trim(files(FILE_SET_EQCM)%filename)//'):' 

    ! latest version of the code
    If (simulation_data%solvation%info%stat) Then
      simulation_data%code_version= '6.1.3.6'
    Else
      simulation_data%code_version= '6.0'
    End If        

    ! Check that atomic tags do not have more than 4 characters
    Do i=1, simulation_data%total_tags
      If(Len_Trim(simulation_data%component(i)%tag)>4)Then
        Write (messages(1),'(1x,a)') '***ERROR: In ONETEP, tags for atomic species must not contain&
                                    & more than 4 characters. Please rename atomic tags'
      End If
    End Do

    ! Check XC_version
    If (Trim(simulation_data%dft%xc_version%type) /= 'pz'     .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'vwn'    .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'pw92'   .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'pw91'   .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'pbe'    .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'rp'     .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'revpbe' .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'pbesol' .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'xlyp'   .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'blyp' ) Then
      Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                  &'Invalid specification for directive "XC_version" for ONETEP.&
                                  & Implemented options for ONETEP are:'
      Write (messages(2),'(1x,a)')   '==== LDA-level =================='
      Write (messages(3),'(1x,2a)')  '- PZ      Perdew-Zunger          ', Trim(bib_pz)  
      Write (messages(4),'(1x,2a)')  '- VWN     Vosko-Wilk-Nusair      ', Trim(bib_vwn)
      Write (messages(5),'(1x,2a)')  '- PW92    Perdew-Wang 92         ', Trim(bib_pw92)
      Write (messages(6),'(1x,a)')   '==== GGA-level =================='
      Write (messages(7),'(1x,2a)')  '- PW91    Perdew-Wang 91         ', Trim(bib_pw91)
      Write (messages(8),'(1x,2a)')  '- PBE     Perdew-Burke-Ernzerhof ', Trim(bib_pbe)
      Write (messages(9),'(1x,2a)')  '- RP      Hammer-Hansen-Norskov  ', Trim(bib_rp)
      Write (messages(10),'(1x,2a)') '- revPBE  revPBE                 ', Trim(bib_revpbe)
      Write (messages(11),'(1x,2a)') '- PBEsol  PBE for solids         ', Trim(bib_pbesol)
      Write (messages(12),'(1x,2a)') '- BLYP    Becke-Lee-Young-Parr   ', Trim(bib_blyp)
      Write (messages(13),'(1x,2a)') '- XLYP    Xu-Goddard             ', Trim(bib_xlyp)
      Write (messages(14),'(1x,a)')  '================================='
      Call info(messages, 14)
      Call error_stop(' ')
    End If

    ! XC base
    If (Trim(simulation_data%dft%xc_version%type) == 'pz'    .Or.&
      Trim(simulation_data%dft%xc_version%type)   == 'pw92'  .Or.&
      Trim(simulation_data%dft%xc_version%type)   == 'vwn'    ) Then
      simulation_data%dft%xc_base='lda'
    Else If (Trim(simulation_data%dft%xc_version%type) == 'pw91'   .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'pbe'    .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'rp'     .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'revpbe' .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'pbesol' .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'xlyp'   .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'blyp') Then
      simulation_data%dft%xc_base='pbe'
    End If

    ! Pseudopotentials
    If (simulation_data%dft%pp_info%stat) Then
      Call check_pseudo_potentials_onetep(simulation_data, error_dft)
    End If 

    ! vdW settings 
    simulation_data%dft%need_vdw_kernel=.False.
    If (simulation_data%dft%vdw%fread) Then
      If (Trim(simulation_data%dft%vdw%type) /= 'dft-d2'  .And.&
         Trim(simulation_data%dft%vdw%type) /= 'vdw-df'   .And.&
         Trim(simulation_data%dft%vdw%type) /= 'optb88'   .And.&
         Trim(simulation_data%dft%vdw%type) /= 'optpbe'   .And.&
         Trim(simulation_data%dft%vdw%type) /= 'vdw-df2'  .And.&
         Trim(simulation_data%dft%vdw%type) /= 'avv10s'   .And.&
         Trim(simulation_data%dft%vdw%type) /= 'vv10'   ) Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                  &'Invalid specification of directive "vdW" for ONETEP. Valid options are:'
        Write (messages(2),'(1x,2a)')  '- DFT-D2   Grimme D2 ', Trim(bib_dftd2)
        Write (messages(3),'(1x,2a)')  '- vdW-DF   X (revPBE), C (LDA), vdW (vdW-DF) ', Trim(bib_vdwdf)  
        Write (messages(4),'(1x,2a)')  '- optPBE   X (OPTPBE), C (LDA), vdW (vdW-DF) ', Trim(bib_optpbe)
        Write (messages(5),'(1x,2a)')  '- optB88   X (OPTB88), C (LDA), vdW (vdW-DF) ', Trim(bib_optb88)
        Write (messages(6),'(1x,2a)')  '- vdW-DF2  X (rPW86), C (LDA), vdW (vdW-DF 2)  ', Trim(bib_vdwdf2)
        Write (messages(7),'(1x,2a)')  '- AVV10S   X (AM05), C (AM05), vdW (rVV10-sol) ', Trim(bib_AVV10s)
        Write (messages(8),'(1x,2a)')  '- VV10     X (rPW86), C (PBE), vdW (rVV10) ', Trim(bib_vv10)
        Call info(messages, 8)
        Call error_stop(' ')
      End If

      If (simulation_data%dft%spin_polarised%stat .And. Trim(simulation_data%dft%vdw%type) /= 'dft-d2') Then
         Write (message,'(1x,2a)') Trim(error_dft), ' Spin-polarisation with non-local vdW corrections is not supported&
                                & in ONETEP. For spin-polarised-vdW use the DFT-D2 method of Grimme.'
        Call error_stop(message)
      End If  
      
      If (Trim(simulation_data%dft%xc_level%type) /= 'gga') Then
        Write (message,'(1x,4a)') Trim(error_dft), &
                                &' Dispersion correction type "', Trim(simulation_data%dft%vdw%type), '" requires the&
                                &  option for directive XC_level. Please change it.'
        Call error_stop(message)
      End If

      ! Find if kernel exists
      If (Trim(simulation_data%dft%vdw%type) == 'vdw-df'   .Or.  &
          Trim(simulation_data%dft%vdw%type) == 'optb88'   .Or.  &
          Trim(simulation_data%dft%vdw%type) == 'optpbe'   .Or.  &
          Trim(simulation_data%dft%vdw%type) == 'vdw-df2'  .Or.  &
          Trim(simulation_data%dft%vdw%type) == 'avv10s'   .Or.  &
          Trim(simulation_data%dft%vdw%type) == 'vv10'   ) Then
          If (Trim(simulation_data%dft%vdw%type) == 'vdw-df'   .And.&
              Trim(simulation_data%dft%vdw%type) == 'optb88'   .And.&
              Trim(simulation_data%dft%vdw%type) == 'optpbe'   .And.&
              Trim(simulation_data%dft%vdw%type) == 'vdw-df2') Then
             simulation_data%dft%vdw_kernel_file='vdW_df_kernel' 
          Else
            simulation_data%dft%vdw_kernel_file='vdW_vv10_kernel' 
          End If 
          Call execute_command_line('[ -f '//Trim(path)//Trim(simulation_data%dft%vdw_kernel_file)//' ]', exitstat=ifile)
          If (ifile/=0) Then 
            Write (messages(1),'(1x,5a)') '***WARNING: Kernel file "', Trim(simulation_data%dft%vdw_kernel_file), '" needed for&
                                     & vdW corrections "', Trim(simulation_data%dft%vdw%type), '" is not found wihtin folder DFT.'
            Write (messages(2),'(1x, a)') '   ONETEP will generate this kernel file during the calculation.'   
            Call info(messages,2)
          Else
            simulation_data%dft%need_vdw_kernel=.True.
          End If
      End If

      ! Check consistency between XC version and vdW non-local corrections
      If (Trim(simulation_data%dft%vdw%type) == 'vdw-df') Then
        If (Trim(simulation_data%dft%xc_version%type)  /= 'revpbe') Then
          Write (message,'(1x,4a)') Trim(error_dft), ' XC_version "', Trim(simulation_data%dft%xc_version%type), &
                           & '" is incompatible with the vdW-DF correction. Change XC_version to revPBE.'
          Call error_stop(message)
        End If
      End If

      If (Trim(simulation_data%dft%vdw%type) == 'optpbe'   .Or.&
         Trim(simulation_data%dft%vdw%type) == 'optb88'   .Or.&
         Trim(simulation_data%dft%vdw%type) == 'vdw-df2'      .Or.&
         Trim(simulation_data%dft%vdw%type) == 'vv10' .Or.&
         Trim(simulation_data%dft%vdw%type) == 'avv10s'   ) Then
         If (Trim(simulation_data%dft%vdw%type) == 'optpbe') Then 
           simulation_data%dft%xc_version%type ='or'
         Else If (Trim(simulation_data%dft%vdw%type) == 'optb88') Then
           simulation_data%dft%xc_version%type ='bo'
         Else If (Trim(simulation_data%dft%vdw%type) == 'vdw-df2') Then
           simulation_data%dft%xc_version%type ='ml'
         Else If (Trim(simulation_data%dft%vdw%type) == 'avv10s') Then
           simulation_data%dft%xc_version%type ='am05'
         Else If (Trim(simulation_data%dft%vdw%type) == 'vv10') Then
           simulation_data%dft%xc_version%type ='rpw86pbe'
         End If
         Call info(' ', 1)
         Write (messages(1),'(1x,5a)') '***WARNING: XC_version will be changed to "', Trim(simulation_data%dft%xc_version%type),&
                                 & '" to include set the requested "',  Trim(simulation_data%dft%vdw%type),&
                                 & '" type of dispersion corrections'   
         Call info(messages,1)
      End If
      
    End If

    !Check magnetization has the Hubbard block defined
    If (simulation_data%dft%mag_info%fread) Then
      If (.Not. simulation_data%dft%hubbard_info%fread) Then
        Write (messages(1),'(1x,2a)') Trim(error_dft), ' In ONETEP the block &hubbard is required to set the initial magnetization&
                                      & for the atomic sites, even if the system does not need Hubbard corrections.'
        Write (messages(2),'(1x,a)') 'If the system does not need Hubbard corrections, the user must set U=J=0 and specify the&
                                    & "l_orbital" (see manual for the specification of block &hubbard).'
        Write (messages(3),'(1x,a)') 'The values for the initial magnetizations provided in block &magnetization are used to apply&
                                    & a spin-splitting (sigma) to the corresponding subspace.'
        Call info(messages, 3)
        Call error_stop(' ')  
      End If
    End If


    ! Prevent Orbital transformation
    If (simulation_data%dft%ot%stat) Then
      Write (message,'(2(1x,a))') Trim(error_dft), 'Requested Orbital Transformation via directive "OT"&
                                 & is not possible for ONETEP simulations. Please remove it'
      Call error_stop(message)
    End If

    ! GAPW 
    If (simulation_data%dft%gapw%stat) Then
      Write (message,'(2(1x,a))') Trim(error_dft), 'Requested "Gaussian Augmented Plane Waves" method via the "gapw"&
                                 & diretitve is not possible for ONETEP simulations. Please remove it'
      Call error_stop(message)
    End If

    ! Energy cutoff
    If (Trim(simulation_data%dft%encut%units)/='ev') Then
       Write (message,'(2(1x,a))') Trim(error_dft), &
                                   &'Units for directive "energy_cutoff" for ONETEP simulations must be in eV'
       Call error_stop(message)
    End If
    simulation_data%dft%encut%units='eV'

   ! precision
    If (simulation_data%dft%precision%fread) Then
      Write (message,'(2(1x,a))') Trim(error_dft), 'For ONETEP, "precision" directive is not needed. Please remove it'
      Call error_stop(message)
    End If
   
    ! SCF energy tolerance 
    If (simulation_data%dft%delta_e%fread) Then
      If (Trim(simulation_data%dft%delta_e%units) /= 'ev' ) Then
         Write (message,'(2a)')  Trim(error_dft), ' Units for directive "SCF_energy_tolerance" in ONETEP must be eV'
         Call info(message, 1)
         Call error_stop(' ')
      End If
    End If

    ! There must be only one k-point (Gamma) for ONETEP
    If (simulation_data%dft%total_kpoints /= 1) Then
       Write (messages(1),'(2(1x,a))') Trim(error_dft), 'In ONETEP, only one k-point (Gamma) is allowed. The user can&
                                      & simply remove directive "kpoints" from the SET_EQCM and rerun.'
       Call info(messages, 1)
       Call error_stop(' ')
    End If

    ! EDFT for metallic sytems
    If (simulation_data%dft%edft%stat) Then
      If (simulation_data%dft%smear%fread) Then
        If (Trim(simulation_data%dft%smear%type) /= 'fermi') Then
          Write (messages(1),'(2(1x,a))') Trim(error_dft), 'The only allowed smearing scheme for EDFT in ONETEP is&
                                          & "fermi". Please change the setting of "smearing" accordingly.'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      Else
        Write (messages(1),'(2(1x,a))') Trim(error_dft), 'The user must set option "fermi" for directive "smearing"&
                                       & for EDFT simulations with ONETEP.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
      ! Width of smearing
      If (.Not. simulation_data%dft%width_smear%fread) Then
        simulation_data%dft%width_smear%value= 0.20_wp
        simulation_data%dft%width_smear%units= 'eV'
      End If 
    Else
      ! Both "smearing" and "width_smear" are  incompatible if the simulation is not EDFT
      If (simulation_data%dft%smear%fread) Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft), 'In ONETEP, the definition of "smearing" is meaningless&
                                        & if the simulation is not EDFT. Please review the settings.' 
        Call info(messages, 1)
        Call error_stop(' ')
      End If

      If (simulation_data%dft%width_smear%fread) Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft), 'In ONETEP, the definition of "width_smear" is meaningless&
                                        & if the simulation is not EDFT. Please review the settings.' 
        Call info(messages, 1)
        Call error_stop(' ')
      End If

      If (simulation_data%dft%bands%fread) Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft), 'In ONETEP, the definition of "bands" is meaningless&
                                        & if the simulation is not EDFT. Please review the settings.' 
        Call info(messages, 1)
        Call error_stop(' ')
      End If

    End If

   ! GC-DFT functionality
   If (simulation_data%dft%gc%activate%stat) Then
     If (.Not. simulation_data%dft%edft%stat) Then
       Write (message,'(2(1x,a))')  Trim(error_gcdft), 'In ONETEP, requested GC-DFT simulations need Ensemble DFT treatment.&
                                   & Set directive "edft" to .True. and rerun.'
       Call error_stop(message)
     End If
     If (simulation_data%dft%gc%reference_potential%fread) Then
       If (simulation_data%dft%gc%reference_potential%fail) Then
          Write (message,'(2(1x,a))') Trim(error_gcdft), 'Wrong (or missing) settings for "reference_potential" directive.'
          Call error_stop(message)
       Else
         If (Trim(simulation_data%dft%gc%reference_potential%units) /= 'ev') Then
            Write (message,'(3a)')  Trim(error_gcdft), ' Invalid units of directive "reference_potential". Units must be "eV"'
            Call error_stop(message)
         End If
       End If
     Else
       Write (message,'(2(1x,a))')  Trim(error_gcdft), 'Requested GC-DFT simulation needs the definition of&
                                   & the "reference_potential" directive in the sub-block &gcdft.'
       Call error_stop(message)
     End If

     If (simulation_data%dft%gc%electrode_potential%fread) Then
       If (simulation_data%dft%gc%electrode_potential%fail) Then
          Write (message,'(2(1x,a))') Trim(error_gcdft), 'Wrong (or missing) settings for "electrode_potential" directive.'
          Call error_stop(message)
       Else
         If (Trim(simulation_data%dft%gc%electrode_potential%units) /= 'v') Then
            Write (message,'(3a)')  Trim(error_gcdft), ' Invalid units of directive "electrode_potential". Units must be "V"'
            Call error_stop(message)
         End If
       End If
     Else
       Write (message,'(2(1x,a))')  Trim(error_gcdft), 'Requested GC-DFT simulation needs the definition of&
                                   & the "electrode_potential" directive in the sub-block &gcdft'
       Call error_stop(message)
     End If

     If (simulation_data%dft%gc%electron_threshold%fread) Then
       If (simulation_data%dft%gc%electron_threshold%fail) Then
          Write (message,'(2(1x,a))') Trim(error_gcdft), 'Wrong (or missing) settings for "electron_threshold" directive.'
          Call error_stop(message)
       Else
         If (simulation_data%dft%gc%electron_threshold%value < epsilon(1.0_wp)) Then
            Write (message,'(3a)')  Trim(error_gcdft), ' Input value for "electron_threshold" MUST be larger than zero'
            Call error_stop(message)
         End If
       End If
     Else
       Write (message,'(2(1x,a))')  Trim(error_gcdft), 'Requested GC-DFT simulation needs the definition of&
                                   & the "electron_threshold" directive in the sub-block &gcdft.'
       Call error_stop(message)
     End If
   End If

    ! Mixing
    If (simulation_data%dft%mixing%fread) Then
      If (simulation_data%dft%edft%stat) Then
        If (Trim(simulation_data%dft%mixing%type) /= 'damp_fixpoint'  .And.&
           Trim(simulation_data%dft%mixing%type)  /= 'pulay')   Then
           Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                        &'Invalid specification of directive "mixing_scheme" in ONETEP. Options for EDFT are:'
           Write (messages(2),'(1x,a)') '- Damp_fixpoint'
           Write (messages(3),'(1x,a)') '- Pulay'
           Call info(messages, 3)
           Call error_stop(' ')
        End If
      Else
        If (Trim(simulation_data%dft%mixing%type) /= 'lnv'  .And.&
            Trim(simulation_data%dft%mixing%type) /= 'dnk_linear'  .And.&
            Trim(simulation_data%dft%mixing%type) /= 'ham_linear'  .And.&
            Trim(simulation_data%dft%mixing%type) /= 'dnk_pulay'  .And.&
            Trim(simulation_data%dft%mixing%type) /= 'ham_pulay'  .And.&
            Trim(simulation_data%dft%mixing%type) /= 'dnk_listi'  .And.&
            Trim(simulation_data%dft%mixing%type) /= 'ham_listi'  .And.&
            Trim(simulation_data%dft%mixing%type) /= 'dnk_listb'  .And.&
            Trim(simulation_data%dft%mixing%type) /= 'ham_listb')   Then
           Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                        &'Invalid specification of directive "mixing_scheme" in ONETEP.&
                                        & Options for band gap systems (i.e. no EDFT) are:'
           Write (messages(2),'(1x,a)') '- lnv'
           Write (messages(2),'(1x,a)') '- dnk_linear'
           Write (messages(3),'(1x,a)') '- ham_linear'
           Write (messages(4),'(1x,a)') '- dnk_pulay'
           Write (messages(5),'(1x,a)') '- ham_pulay'
           Write (messages(6),'(1x,a)') '- dnk_listi'
           Write (messages(7),'(1x,a)') '- ham_listi'
           Write (messages(8),'(1x,a)') '- dnk_listb'
           Write (messages(9),'(1x,a)') '- ham_listb'
           Call info(messages, 9)
           Call error_stop(' ')
        End If
      End If
    Else
      If (simulation_data%dft%edft%stat) Then
        simulation_data%dft%mixing%type='damp_fixpoint'
      Else
        simulation_data%dft%mixing%type='lnv'
      End If        
    End If

    ! Check if basis set was defined, complain and abort
    If (simulation_data%dft%basis_info%fread) Then
      Write (message,'(1x,2a)') Trim(error_dft), &
                        &' Definition of basis sets is not required for ONETEP (NGWF are determined self-consistently).&
                        & Please remove sub-block &basis_set and rerun.'
      Call error_stop(message)
    End If

    If (simulation_data%motion%ion_steps%fread) Then
      If (simulation_data%motion%ion_steps%value == 1) Then
        simulation_data%simulation%type='singlepoint'
        Call info(' ***WARNING: since the number of ionic steps was set to 1, the simulation was changed to SinglePoint', 1)
      End If
    End If

  End Subroutine define_onetep_dft    

  Subroutine print_onetep_dft(iunit, ic, net_elements, list_tag, list_number, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print DFT directives for ONETEP
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Integer(Kind=wi), Intent(InOut) :: ic    
    Integer(Kind=wi), Intent(In   ) :: net_elements
    Character(Len=8), Intent(In   ) :: list_tag(max_components) 
    Integer(Kind=wi), Intent(In   ) :: list_number(max_components)
    Type(simul_type), Intent(InOut) :: simulation_data

    Integer(Kind=wi)   :: i, j, k
    Integer(Kind=wi)   :: atomic_number
    Logical            :: loop, loop2
    Real(Kind=wp)      :: mag_ini(max_components)
    Character(Len=256) :: message

    Integer(Kind=wi)   :: ici

    Write (iunit,'(a)')    '    '
    Write (iunit,'(a)') '##### Electronic structure'
    Write (iunit,'(a)') '#========================='

    If (.Not. simulation_data%dft%vdw%fread) Then 
      Write (iunit,'(a)') '#==== Exchange and correlation'
      If (Trim(simulation_data%dft%xc_version%type) == 'pz') Then
        Write (iunit,'(2a)') '# Perdew-Zunger (PZ) functional ', Trim(bib_pz) 
        Write (message,'(a)')  'xc_functional :  CAPZ'
      Else If (Trim(simulation_data%dft%xc_version%type) == 'vwn') Then 
        Write (iunit,'(2a)') '# Vosko-Wilk-Nusair (VWN) functional ', Trim(bib_vwn)
        Write (message,'(a)')  'xc_functional :  VWN'
      Else If (Trim(simulation_data%dft%xc_version%type) == 'pw92') Then 
        Write (iunit,'(2a)') '# Perdew-Wang 92 (PW92) functional ', Trim(bib_pw92)
        Write (message,'(a)')  'xc_functional :  PW92'
      Else If (Trim(simulation_data%dft%xc_version%type) == 'pw91') Then 
        Write (iunit,'(2a)') '# Perdew-Wang 91 (PW91) XC functional ', Trim(bib_pw91) 
        Write (message,'(a)')  'xc_functional :  PW91'
      Else If (Trim(simulation_data%dft%xc_version%type) == 'pbe') Then
        Write (iunit,'(2a)') '# Perdew-Burke-Ernzerhof (PBE) XC functional ', Trim(bib_pbe)
        Write (message,'(a)')  'xc_functional :  PBE'
      Else If (Trim(simulation_data%dft%xc_version%type) == 'rp') Then 
        Write (iunit,'(2a)') '# Hammer-Hansen-Norskov (RPEB) XC functional ', Trim(bib_rp)
        Write (message,'(a)')  'xc_functional :  RPBE'
      Else If (Trim(simulation_data%dft%xc_version%type) == 'revpbe') Then 
        Write (iunit,'(2a)') '# revPBE XC functional ', Trim(bib_revpbe)
        Write (message,'(a)')  'xc_functional :  REVPBE'
      Else If (Trim(simulation_data%dft%xc_version%type) == 'pbesol') Then 
        Write (iunit,'(2a)') '# PBE for solids (PBEsol) XC functional ', Trim(bib_pbesol) 
        Write (message,'(a)')  'xc_functional :  PBESOL'
      Else If (Trim(simulation_data%dft%xc_version%type) == 'xlyp') Then 
        Write (iunit,'(2a)') '# Xu-Goddard (XLYP) XC functional ', Trim(bib_xlyp)
        Write (message,'(a)')  'xc_functional :  XLYP'
      Else If (Trim(simulation_data%dft%xc_version%type) == 'blyp') Then
        Write (iunit,'(2a)') '# Becke-Lee-Young-Parr (BLYP) XC functional ', Trim(bib_blyp)
        Write (message,'(a)')  'xc_functional :  BLYP'
      End If 
      Call record_directive(iunit, message, 'xc_functional', simulation_data%set_directives%array(ic), ic)

    Else ! vdW corrections
      Write (iunit,'(3a)') '#==== Exchange and correlation + ', Trim(simulation_data%dft%vdw%type),&
                           & ' dispersion corrections'
      If (Trim(simulation_data%dft%vdw%type) == 'dft-d2') Then
        Write (iunit,'(2a)') '# Perdew-Burke-Ernzerhof (PBE) XC functional ', Trim(bib_pbe)
        Write (message,'(a)')  'xc_functional :  PBE'
        Call record_directive(iunit, message, 'xc_functional', simulation_data%set_directives%array(ic), ic)
        Write (iunit,'(2a)') '# Damping correction of Grimme DFT-D2 ', Trim(bib_dftd2)
        Write (message,'(a)')  'dispersion    :    4'
        Call record_directive(iunit, message, 'dispresion', simulation_data%set_directives%array(ic), ic)
      Else If (Trim(simulation_data%dft%vdw%type) == 'optpbe') Then
        Write (iunit,'(2a)') '# vdW-optPBE non-local corrections ', Trim(bib_optpbe)
        Write (message,'(a)')  'xc_functional :  optPBE'
        Call record_directive(iunit, message, 'xc_functional', simulation_data%set_directives%array(ic), ic)
      Else If (Trim(simulation_data%dft%vdw%type) == 'optb88') Then 
        Write (iunit,'(2a)') '# vdW-optB88 non-local corrections ', Trim(bib_optb88)
        Write (message,'(a)')  'xc_functional :  optB88'
        Call record_directive(iunit, message, 'xc_functional', simulation_data%set_directives%array(ic), ic)
      Else If (Trim(simulation_data%dft%vdw%type) == 'vdw-df') Then 
        Write (iunit,'(2a)') '# vdW-DF non-local corrections ', Trim(bib_vdwdf)
        Write (message,'(a)')  'xc_functional :  vdWDF'
        Call record_directive(iunit, message, 'xc_functional', simulation_data%set_directives%array(ic), ic)
      Else If (Trim(simulation_data%dft%vdw%type) == 'vdw-df2') Then 
        Write (iunit,'(2a)') '# vdW-DF2 non-local corrections ', Trim(bib_vdwdf2)
        Write (message,'(a)')  'xc_functional :  vdWDF2'
        Call record_directive(iunit, message, 'xc_functional', simulation_data%set_directives%array(ic), ic)
      Else If (Trim(simulation_data%dft%vdw%type) == 'vv10') Then 
        Write (iunit,'(2a)') '# vdW-VV10 non-local corrections ', Trim(bib_vv10)
        Write (message,'(a)')  'xc_functional :  VV10'
        Call record_directive(iunit, message, 'xc_functional', simulation_data%set_directives%array(ic), ic)
      Else If (Trim(simulation_data%dft%vdw%type) == 'avv10s') Then 
        Write (iunit,'(2a)') '# vdW-AVV10s non-local corrections ', Trim(bib_AVV10s)
        Write (message,'(a)')  'xc_functional :  AVV10S'
        Call record_directive(iunit, message, 'xc_functional', simulation_data%set_directives%array(ic), ic)
      End If 
    End If
    
    Write (iunit,'(a)') '#==== Convergence parameters'
    Write (message,'(a,i4,a)') 'maxit_ngwf_cg : ', simulation_data%dft%scf_steps%value,&
                           & ' # maximum number of iterations for the NGWF conjugate gradients optimization'
    Call record_directive(iunit, message, 'maxit_ngwf_cg', simulation_data%set_directives%array(ic), ic)
    Write (message,'(a,f8.2,1x,2a)') 'cutoff_energy : '  , simulation_data%dft%encut%value, Trim(simulation_data%dft%encut%units), &
                              & '  # Energy cutoff'
    Call record_directive(iunit, message, 'cutoff_energy', simulation_data%set_directives%array(ic), ic)
    Write (message,'(a,e10.3,1x,2a)') 'elec_energy_tol :'  , simulation_data%dft%delta_e%value, &
                                   & Trim(simulation_data%dft%delta_e%units), ' # Energy tolerance'
    Call record_directive(iunit, message, 'elec_energy_tol', simulation_data%set_directives%array(ic), ic)
    Write (message,'(a)')  'kernel_update :  T  # Update the density kernel when taking a trial step for NGWF optimization'
    Call record_directive(iunit, message, 'kernel_update', simulation_data%set_directives%array(ic), ic)

    ! net charge                               
    If (simulation_data%net_charge%fread) Then
      Write (iunit,'(a)') '#==== Net charge'
      Write (message,'(a,i6)') 'charge : ', Nint(simulation_data%net_charge%value)
      Call record_directive(iunit, message, 'charge', simulation_data%set_directives%array(ic), ic)
    End If
    
    ! Spin polarised
    If (simulation_data%dft%spin_polarised%stat) Then
      Write (iunit,'(a)') '#==== Spin polarised calculation'
      Write (message,'(a)')      'spin_polarised :   T  ' 
      Call record_directive(iunit, message, 'spin_polarised', simulation_data%set_directives%array(ic), ic)
    End If

    ! EDFT
    If (simulation_data%dft%edft%stat) Then
      Write (iunit,'(a)')           ' '
      Write (iunit,'(a)') '#==== Ensemble DFT'
      Write (message,'(a)')           'edft  :  T'
      Call record_directive(iunit, message, 'edft', simulation_data%set_directives%array(ic), ic)
      Write (message,'(a,f8.2,2x,a)') 'edft_smearing_width : ', simulation_data%dft%width_smear%value, ' eV  # smearing width'  
      Call record_directive(iunit, message, 'edft_smearing_width', simulation_data%set_directives%array(ic), ic)
      ! Mixing
      Write (message,'(2a)') 'edft_update_scheme  :  ', Trim(simulation_data%dft%mixing%type)
      Call record_directive(iunit, message, 'edft_update_scheme', simulation_data%set_directives%array(ic), ic)

      If (simulation_data%dft%bands%fread) Then
          Write (message,'(a,i4,a)') 'edft_extra_bands : ', simulation_data%dft%bands%value, &
                                   & ' # add extra bands to reach/improve convergence'
      Else 
          Write (message,'(a)')     'edft_extra_bands :  -1 # bands is equal to the total number of NGWFs' 
      End If
      Call record_directive(iunit, message, 'edft_extra_bands', simulation_data%set_directives%array(ic), ic)
    Else
      ! DIIS Mixing
      Write (message,'(2a)') 'kernel_diis_scheme  :  ', Trim(simulation_data%dft%mixing%type)
      Call record_directive(iunit, message, 'kernel_diis_scheme', simulation_data%set_directives%array(ic), ic)
    End If

    ! GC-DFT
    If (simulation_data%dft%gc%activate%stat) Then
      Write (iunit,'(a)')           ' '
      Write (iunit,'(a)') '#==== Grand Canonical (GC)'
      Write (message,'(a)')           'edft_grand_canonical  :  T'
      Call record_directive(iunit, message, 'edft_grand_canonical', simulation_data%set_directives%array(ic), ic)
      Write (message,'(a,f8.2,2x,a)') 'edft_reference_potential : ', simulation_data%dft%gc%reference_potential%value, ' eV'  
      Call record_directive(iunit, message, 'dft_reference_potential', simulation_data%set_directives%array(ic), ic)
      Write (message,'(a,f8.2,2x,a)') 'edft_electrode_potential : ', simulation_data%dft%gc%electrode_potential%value, ' V'  
      Call record_directive(iunit, message, 'edft_electrode_potential', simulation_data%set_directives%array(ic), ic)
      Write (message,'(a,e10.3)')      'edft_nelec_thres : ', simulation_data%dft%gc%electron_threshold%value  
      Call record_directive(iunit, message, 'edft_nelec_thres', simulation_data%set_directives%array(ic), ic)
    End If

    
   ! Magnetization 
   If (simulation_data%dft%mag_info%fread) Then
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
     ! Fix total magnetization 
     If (simulation_data%dft%total_magnetization%fread) Then
       Write (iunit,'(a)') ' '
       Write (iunit,'(a)') '#==== Total magnetization'
       If (simulation_data%dft%edft%stat) Then
         Write (message,'(a,f7.3)') 'spin  : ', simulation_data%dft%total_magnetization%value
       Else
         Write (message,'(a,i4)') 'spin  : ', Int(simulation_data%dft%total_magnetization%value)
       End If 
       Call record_directive(iunit, message, 'spin', simulation_data%set_directives%array(ic), ic)

       Call check_initial_magnetization(net_elements, list_tag, list_number, mag_ini,&
                                      & simulation_data%dft%total_magnetization%value)
       If (simulation_data%dft%edft%stat) Then
         Write (message,'(a)') 'edft_spin_fix  :   -1' 
         Call record_directive(iunit, message, 'edft_spin_fix', simulation_data%set_directives%array(ic), ic)
       End If

     End If
     ! Set initial magnetization over atomic sites
       Write (iunit,'(a)') ' '
       If (simulation_data%dft%hubbard_info%fread) Then
          Write (iunit,'(a)') '#==== Hubbard corrections + initial values of sigma'     
       Else
          Write (iunit,'(a)') '#==== Initial magnetization (use of &hubbard block with no&
                             & corrections but sigma values)'     
       End If
       Write (iunit,'(a)') '%block hubbard'
       Do i=1, net_elements
         j=1
         loop=.True.
         Do While (j <= simulation_data%total_tags .And. loop)
           If (Trim(list_tag(i))==Trim(simulation_data%dft%hubbard(j)%tag)) Then
             If (Abs(simulation_data%dft%hubbard(i)%U) > epsilon(1.0_wp) .Or. &
                 Abs(mag_ini(i)) > epsilon(1.0_wp)) Then 
               Write (iunit,'(1x,a5,2x,i1,(2(2x,f5.2)),2x,a,(2(2x,f5.2)))') Trim(list_tag(i)),&
                                                  &simulation_data%dft%hubbard(j)%l_orbital,&
                                                  & simulation_data%dft%hubbard(j)%U,& 
                                                  & simulation_data%dft%hubbard(j)%J,&
                                                  & '-10', 0.0_wp, mag_ini(i) 
               loop=.False.
             End If
           End If
           j=j+1
         End Do
       End Do
       Write (iunit,'(a)') '%endblock hubbard'
   End If

   ! Print pseudo_potentials
   If (simulation_data%dft%pp_info%stat) Then
     ! Pseudo potentials
     Write (iunit,'(a)') '  '
     Write (iunit,'(a)') '#==== Pseudopotentials'
     If (simulation_data%dft%onetep_paw) Then
       Write (iunit,'(a)') 'paw : T'
       Call record_directive(iunit, message, 'paw', simulation_data%set_directives%array(ic), ic) 
     End If
     Write (iunit,'(a)') '%block species_pot'
     Do i=1, net_elements
       j=1
       loop=.True.
       Do While (j <= simulation_data%total_tags .And. loop)
         If (Trim(list_tag(i))==Trim(simulation_data%dft%pseudo_pot(j)%tag)) Then
           Write (iunit,'(1x,a5,2x,a)') Trim(list_tag(i)), Trim(simulation_data%dft%pseudo_pot(j)%file_name)
           loop=.False.
         End If
         j=j+1
       End Do
     End Do   
     Write (iunit,'(a)') '%endblock species_pot'
     ! Initial Pseudo-atomic orbitals 
     Write (iunit,'(a)') '  '
     Write (iunit,'(a)') '#==== Initial pseudo-atomic orbital set'
     Write (iunit,'(a)') '%block species_atomic_set'    
     Do i=1, net_elements
       Write (iunit,'(1x,a5,2x,a)') Trim(list_tag(i)), 'SOLVE'
     End Do   
     Write (iunit,'(a)') '%endblock species_atomic_set'    
   End If
    
   ! block species
   Write (iunit,'(a)') '  '
   Write (iunit,'(a)') '#==== Block species'
   Write (iunit,'(a)') '%block species'
   Write (iunit,'(a)') 'ang'
   Do i=1, net_elements
     j=1
     loop=.True.
     Do While (j <= simulation_data%total_tags .And. loop)
       If (Trim(list_tag(i))==Trim(simulation_data%dft%ngwf(j)%tag)) Then
         loop=.False.
         loop2=.True.
         k=1
         Do While (k <= simulation_data%total_tags .And. loop2)
           If (Trim(list_tag(i))==Trim(simulation_data%component(k)%tag)) Then
             atomic_number=simulation_data%component(k)%atomic_number
             loop2=.False.
             ! Check if the size of the simulation cell is adequate for the radii of the NGWF defined	
             Do ici = 1, 3
               If (2*simulation_data%dft%ngwf(j)%radius >= simulation_data%cell_length(ici)) Then
                 Write (message,'(1x, 3a,f8.2,a)') '***PROBLEMS: the NGWF diameter of species ', &
                                     & Trim(simulation_data%dft%ngwf(j)%tag), ' is ',&
                                     & 2*simulation_data%dft%ngwf(j)%radius, ' Angstrom, which is larger than the size&
                                     & of the supercell at least in one of the three-directions. Please check and&
                                     & enlarge the size of the model with directive "repeat_input_model".'
                 Call error_stop(message)
               End If
             End Do
             
             Write (iunit,'(1x,2a5,2i5,f6.2)') Trim(simulation_data%dft%ngwf(j)%tag),  &
                                               & Trim(simulation_data%dft%ngwf(j)%element),&
                                               & atomic_number, &
                                               & simulation_data%dft%ngwf(j)%ni,       &
                                               & simulation_data%dft%ngwf(j)%radius
           End If
           k=k+1
         End Do
       End If
       j=j+1
     End Do
   End Do   
   Write (iunit,'(a)') '%endblock species'   

  End Subroutine print_onetep_dft

  Subroutine advise_dft_onetep(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about DFT settings for ONETEP 
    !
    ! author    - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data
    
    Character(Len=256)  :: messages(9)
    Character(Len=256)  :: in_extra
    Integer(Kind=wi)    :: num_ref_data

    ! Reference data to compare against &extra_directives
    Type(type_ref_data) :: ref_extra_data(max_directives)
    Character(Len=256)  :: exceptions(10)
    Logical             :: extradir_header 

    in_extra='using the &extra_directives block'
    
    Call info(' ', 1)
    Write (messages(1), '(1x,a)')  'In case of problems in the electronic convergence, the user should:'
    Write (messages(2), '(1x,a)')  ' - check the number and radii of the NGWF for each species via block &ngwf' 
    Write (messages(3), '(1x,a)')  ' - increase the value of "maxit_ngwf_cg" via ALC_EQCM directive "scf_steps"'
    Write (messages(4), '(1x,2a)') ' - optimise the initialisation of the density kernel via "maxit_palser_mano"&
                                    & and "maxit_pen" ', Trim(in_extra)

    If (simulation_data%dft%edft%fread) Then
      Write (messages(5), '(1x,a)')  '===== Hints for EDFT'
      Write (messages(6), '(1x,a)')  ' - increase the value of directive "bands"'
      Write (messages(7), '(1x,a)')  ' - change the value of directive "width_smearing"'
      Call info(messages, 7)
      Call set_reference_database(ref_extra_data, num_ref_data,'ONETEP', 'EDFT')
      exceptions(1)='edft_ham_diis_size'
      extradir_header=.False.      
      Call check_settings_set_extra_directives(ref_extra_data, num_ref_data, simulation_data%extra_directives,&
                                         & exceptions, 1, extradir_header)
      If (Trim(simulation_data%dft%mixing%type) == 'damp_fixpoint') Then
        Call check_settings_single_extra_directive('edft_ham_diis_size', ref_extra_data, num_ref_data,&
                                             & simulation_data%extra_directives, extradir_header, 'mixing_scheme')      
        Write (messages(1), '(1x,a)')  'If the above does not work, change directive "mixing_scheme" to "Pulay"'
      Else If (Trim(simulation_data%dft%mixing%type) == 'pulay') Then
        Call check_settings_single_extra_directive('edft_ham_diis_size', ref_extra_data, num_ref_data,&
                                             & simulation_data%extra_directives, extradir_header)      
        Write (messages(1), '(1x,a)')  'If the above does not work, change directive "mixing_scheme" to "Damp_Fixpoint"'
      End If        
      Call info(messages, 1)
    Else
      If (Trim(simulation_data%dft%mixing%type) /= 'lnv') Then
        Write (messages(5), '(1x,3a)') ' - define/change LNV related keywords ("minit_lnv", "maxit_lnv", etc) ',&
                                       & Trim(in_extra), ' (see ONETEP manual). Alternatively, change the option for&
                                       & directive "mixing_scheme"' 
      Else
        Write (messages(5), '(1x,3a)') ' - define/change KERNEL_DIIS keywords for control, tolerance and/or level-shift ',&
                               & Trim(in_extra), ' (see ONETEP manual). Alternatively, change the option for directive&
                               & "mixing_scheme"' 
      End If        
      Write (messages(6), '(1x,a)')  ' - IMPORTANT: check if the system tends to metallise, in which case set "EDFT" to .True.'
      Call info(messages, 6)

      Write (messages(1), '(1x,2a)') 'For linear scaling DFT the user should optimise "kernel_cutoff" ', Trim(in_extra)
      Call info(messages, 1)
    End If   

  End Subroutine advise_dft_onetep

  Subroutine warnings_onetep_dft(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to recommed the user about DFT settings for ONETEP 
    !
    ! author    - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data
    
    Character(Len=256)  :: messages(9), header
    Character(Len=256)  :: in_extra
    Logical             :: error
    Integer(Kind=wi)    :: i

    Logical             :: warning, print_header

    in_extra='using the &extra_directives block'

    warning=.False.
    print_header=.True.

    If (simulation_data%dft%hubbard_info%fread  .Or. simulation_data%net_charge%fread  .Or. simulation_data%dft%vdw%fread) Then
       warning=.True.
    End If

    
    If (warning) Then
      print_header=.True.
      Call info(' ', 1)
      Write (header, '(1x,a)')  '***IMPORTANT*** From the requested settings of "&block_simulation_settings", it is&
                                    & RECOMMENDED to consider:'
      ! Total charge
      If (simulation_data%net_charge%fread) Then
        If (Abs(simulation_data%net_charge%value) > epsilon(1.0_wp)) Then
          Write (messages(1), '(1x,a,f8.3,a)')  ' - the meaning/correctness of having set a net charge of ',&
                                  & simulation_data%net_charge%value, ' for the modelled system(s).'
          Write (messages(2), '(1x,2a)')        ' - truncating the Coulomb interactions (if the charged system includes vacuum)&
                                                & via directives "coulomb_cutoff_type" and "coulomb_cutoff_radius" ',&
                                                & Trim(in_extra)
          Call print_warnings(header, print_header, messages, 2)
        End If
      End If

      If (simulation_data%dft%mag_info%fread) Then
        Write (messages(1), '(1x,a)')  ' - checking the convergence of the magnetic solution'
        Call print_warnings(header, print_header, messages, 1)
      End If
      
   
      ! Hubbard-related parameters
      If (simulation_data%dft%hubbard_info%fread) Then
        If (.Not. simulation_data%dft%hubbard_all_U_zero) Then
          Write (messages(1), '(1x,2a)')  ' - further optimization of electronic minimization parameters&
                                         & (if there are problems from the inclusion of Hubbard corrections)'
          Call print_warnings(header, print_header, messages, 1)
        End If
      End If

      ! vdW related parameters
      If (simulation_data%dft%vdw%fread) Then
        If (Trim(simulation_data%dft%vdw%type) == 'dft-d2') Then
          error=.False.
          Do i=1, simulation_data%total_tags
            If (simulation_data%component(i)%atomic_number > 54) Then
              error=.True.
            End If
          End Do
          If (error) Then 
            Write (messages(1),'(1x,a)')  ' - revision of the requested DFT-D2 vdW correction:&
                             & defaults parameters are defined only for elements in the first&
                             & five rows of periodic table (i.e. H-Xe).'
            Write (messages(2),'(1x,2a)') '   WARNING: at least one of the defined species are beyond this range and the user&
                                      & must define the correct parameters (via VDW_PARAMS) ', Trim(in_extra) 
            Call print_warnings(header, print_header, messages,2) 
          End If

          If (Trim(simulation_data%dft%xc_version%type) /= 'pbe') Then
             Write (messages(1),'(1x,a)')  ' - revision of the requested DFT-D2 vdW correction:&
                             & the user should manually change the settings for VDW_DCOEFF and VDW_PARAMS'
            Call print_warnings(header, print_header, messages,2) 
          End If

        End If

        If (Trim(simulation_data%dft%vdw%type) == 'vdw-df'  .Or.&
           Trim(simulation_data%dft%vdw%type) == 'optpbe'   .Or.&
           Trim(simulation_data%dft%vdw%type) == 'optb88'   .Or.&
           Trim(simulation_data%dft%vdw%type) == 'vdw-df2'  .Or.&
           Trim(simulation_data%dft%vdw%type) == 'aavv10s' .Or.&
           Trim(simulation_data%dft%vdw%type) == 'vv10'   ) Then
           Write (messages(1),'(1x,3a)')  ' - for the requested "', Trim(simulation_data%dft%vdw%type),&
                                       & '" dispersion correction, "cutoff_energy" should be optimised for accuracy and efficiency'
           Call print_warnings(header, print_header, messages, 1)
        End If 
      End If    
    End If  

    Call info(' ', 1)
    Write (messages(1), '(1x,3a)') 'The efficiency in the parallelization can be optimised ', &
                                Trim(in_extra), ' with the following directives:'
    Write (messages(2), '(1x, a)') ' - threads_max          (number of OpenMP threads in outer loops)'
    Write (messages(3), '(1x, a)') ' - threads_num_fftboxes (number of threads to use in OpenMP-parallel FFTs)'
    Write (messages(4), '(1x, a)') ' - threads_per_fftbox   (number of nested threads used for FFT box operations)'
    Write (messages(5), '(1x, a)') ' - threads_per_cellfft  (number of threads to use in OpenMP-parallel FFTs on simulation cell)'
    Write (messages(6), '(1x, a)') ' - threads_num_mkl      (number of threads to use in MKL routines)'
    Write (messages(7), '(1x, a)') 'WARNING:'
    Write (messages(8), '(1x, a)') ' - "threads_max" must be equal to "threads_num_fftboxes" and both must be consistent with&
                                   & the value of given to "export OMP_NUM_THREADS" in the job submission script'
    Write (messages(9), '(1x, a)') ' - "threads_num_mkl" must only be defined is MKL is used for compilation'
    Call info(messages, 9)

    ! SCALAPACK
    Write (messages(1), '(1x,a)')  ' '
    Write (messages(2), '(1x,a)')  'If ONETEP is interfaced to SCALAPACK, diretives "eigensolver_orfac" and&
                                   & "eigensolver_abstol" can be used.'
    Write (messages(3), '(1x,3a)') 'I/O can be controlled ', Trim(in_extra), ' with the following directives (see ONETEP manual):'
    Write (messages(4), '(1x,a)')  ' - write_denskern, write_tightbox_ngwfs, write_converged_dk_ngwfs, write_hamiltonian'
    Write (messages(5), '(1x,a)')  ' - read_denskern, read_tightbox_ngwfs, read_hamiltonian'
    Write (messages(6), '(1x,a)')  ' - output_detail (to specify the level of detail for the generated output)'
    Call info(messages, 6)

    If (Trim(simulation_data%simulation%type) == 'md' .Or. Trim(simulation_data%simulation%type) == 'relax_geometry') Then 
      Call info(' ', 1)
      Write (messages(1), '(1x,3a)') 'By setting "write_xyz : T" (', Trim(in_extra), ') the atomic coordinates are printed&
                                    & as a .xyz file.'
      Call info(messages, 1)
    End If

  End Subroutine warnings_onetep_dft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Motion  
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine define_onetep_motion(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define (and check) settings for ONETEP directives (set defaults values)
    !
    ! author    - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256) :: messages(15)
    Character(Len=256) :: error_motion

    Integer(Kind=wi)   :: i 
    Logical            :: error
  
    error_motion = '***ERROR in &motion_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
 
    !Relaxation method
    If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      If (Trim(simulation_data%motion%relax_method%type) /= 'bfgs'   .And. &
        Trim(simulation_data%motion%relax_method%type)  /= 'lbfgs'  ) Then
        Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                &'Invalid specification of directive "relax_method" for ONETEP. Implemented options are:'
        Write (messages(2),'(1x,a)') '- BFGS  (Broyden-Fletcher-Goldfarb-Shanno)'
        Write (messages(3),'(1x,a)') '- LBFGS (Linear Broyden-Fletcher-Goldfarb-Shanno)'
        Call info(messages, 3)
        Call error_stop(' ')
      End If

      If (simulation_data%motion%ion_steps%fread) Then
         If (simulation_data%motion%ion_steps%value == 2) Then
          Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                &' In ONETEP, "ion_steps" for geometry relaxation must be larger than 2. Please change.'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If

      ! Only fixed simulation cells 
      If (simulation_data%motion%ion_steps%value > 2) Then
        If (simulation_data%motion%change_cell_volume%stat .Or. simulation_data%motion%change_cell_shape%stat) Then
          Write (messages(1),'(1x,4a)') Trim(error_motion), ' Up to ', Trim(date_RELEASE), ', ONETEP does not&
                                    & allow to perform geometry relaxation where the shape/size of the simulation cell changes.' 
          Write (messages(2),'(1x,a)') 'Therefore, both directives "change_cell_volume" and "change_cell_shape" must be&
                                    & set to .False. or simply removed.'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End If
      
    End If
    
    ! Force tolerance
    error=.False.
    If (simulation_data%motion%delta_f%fread) Then
      If (Trim(simulation_data%motion%delta_f%units(1)) /= 'ev') Then
        error=.True.
      End If
      If (Trim(simulation_data%motion%delta_f%units(2)) /= 'angstrom-1' ) Then
        error=.True.
      End If
    Else
      simulation_data%motion%delta_f%units(1)='eV' 
      simulation_data%motion%delta_f%units(2)='Angstrom-1'
      simulation_data%motion%delta_f%value(1)= 0.01 
    End If
 
    If (error) Then
      Write (messages(1),'(2a)')  Trim(error_motion), ' Invalid units of directive "force_tolerance" for ONETEP.&
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
          Trim(simulation_data%motion%ensemble%type) /= 'nvt') Then
          Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                    &'Invalid specification of "ensemble" for ONETEP. Available options are:'
          Write (messages(2),'(1x,a)') '- NVE (Microcanonical ensemble)'
          Write (messages(3),'(1x,a)') '- NVT (Canonical ensemble)'
          Call info(messages, 3)
          Call error_stop(' ')
        End If
    End If

   ! Thermostat
    If (Trim(simulation_data%simulation%type) == 'md') Then    
      If (simulation_data%motion%thermostat%fread) Then
        If (Trim(simulation_data%motion%thermostat%type) /= 'langevin'  .And. &
           Trim(simulation_data%motion%thermostat%type) /= 'andersen'  .And. &
           Trim(simulation_data%motion%thermostat%type) /= 'nose-hoover'  ) Then
          Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                  &'Specification for "thermostat" is not supported by ONETEP. Options are:'
          Write (messages(2),'(1x,a)') '- Andersen'
          Write (messages(3),'(1x,a)') '- Langevin'
          Write (messages(4),'(1x,a)') '- Nose-Hoover'
          Call info(messages, 4)
          Call error_stop(' ')
        End If
      End If
    End If

    ! Relaxation time for the thermostat
    If (Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
      If (.Not. simulation_data%motion%relax_time_thermostat%fread) Then
          Write (messages(1),'(1x,4a)') Trim(error_motion), ' In ONETEP, thermostat "', &
                                       Trim(simulation_data%motion%thermostat%type), '" requires the specification&
                                       & of "relax_time_thermostat", which is missing.'
          Call info(messages, 1)
          Call error_stop(' ')
      End If
    End If

    If (simulation_data%extra_info%stat) Then
      ! Check if user defined directives contain only symbol ":"
      Do i = 1, simulation_data%extra_directives%N0
        Call check_extra_directives(simulation_data%extra_directives%array(i), &
                                    simulation_data%extra_directives%key(i),   &
                                    simulation_data%extra_directives%set(i), ':', 'ONETEP')
      End Do
    End If

  End Subroutine define_onetep_motion  

  Subroutine print_onetep_motion(iunit, ic, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print motion directives for ONETEP
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Integer(Kind=wi), Intent(InOut) :: ic    
    Type(simul_type), Intent(InOut) :: simulation_data

    Character(Len=256) :: thermo
    Character(Len=256) :: message

   Write (iunit,'(a)') ' '
   If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
     Write (iunit,'(a)') '#### Geometry relaxation'
     Write (iunit,'(a)') '#======================='
     ! atoms      
     Write (message,'(a,i4,a)') 'geom_max_iter : ', simulation_data%motion%ion_steps%value, ' # Number of ionic steps'
     Call record_directive(iunit, message, 'geom_max_iter', simulation_data%set_directives%array(ic), ic) 
     If (Trim(simulation_data%motion%relax_method%type) == 'bfgs') Then
       Write (message, '(a)') 'geom_method :  CARTESIAN # geometry relaxation with Broyden-Fletcher-Goldfarb-Shanno method'
       Call record_directive(iunit, message, 'geom_method', simulation_data%set_directives%array(ic), ic) 
     Else If (Trim(simulation_data%motion%relax_method%type) == 'lbfgs') Then
       Write (message, '(a)') 'geom_lbfgs  :  T  # geometry relaxation with the linear Broyden-Fletcher-Goldfarb-Shanno method'
       Call record_directive(iunit, message, 'geom_lbfgs', simulation_data%set_directives%array(ic), ic) 
     End If

     Write (message,'(a,f6.2,a)') 'geom_force_tol : ', simulation_data%motion%delta_f%value(1), ' ev/ang'   
     Call record_directive(iunit, message, 'geom_lbfgs', simulation_data%set_directives%array(ic), ic) 

   Else If (Trim(simulation_data%simulation%type) == 'md') Then
     Write (iunit,'(a)') '#### Molecular dynamics'
     Write (iunit,'(a)') '#======================'
     Write (message,'(a,f6.2,a)') 'md_delta_t  :    ', simulation_data%motion%timestep%value,    ' fs'
     Call record_directive(iunit, message, 'md_delta_t', simulation_data%set_directives%array(ic), ic) 
     Write (message,'(a,i5,a)')   'md_num_iter    : ', simulation_data%motion%ion_steps%value, ' # Number of MD steps'
     Call record_directive(iunit, message, 'md_num_iter', simulation_data%set_directives%array(ic), ic) 

     If (Trim(simulation_data%motion%ensemble%type) == 'nve') Then
       Write (iunit, '(a)') '#==== Thermostat "None" for NVE ensemble'
       Write (iunit, '(a)') '%block thermostat'
       Write (iunit,'(2x,a,2x,i5,2x,a,f10.2,a)')  '0', simulation_data%motion%ion_steps%value, &
                                                & 'None', simulation_data%motion%temperature%value, ' K'
       Write (iunit, '(a)') '%endblock thermostat'   
     Else If (Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
       Write (iunit, '(a)') '#==== Thermostat for the NVT ensemble'
       Write (iunit, '(a)') '%block thermostat'
       If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover') Then
         thermo = 'nosehoover'
       Else If (Trim(simulation_data%motion%thermostat%type) == 'langevin') Then
         thermo = 'langevin'
       Else If (Trim(simulation_data%motion%thermostat%type) == 'andersen') Then
         thermo = 'andersen'
       End If
       Write (iunit,'(2x,a,2x,i5,2x,a,f10.2,a)')  '0', simulation_data%motion%ion_steps%value, &
                                                & Trim(thermo), simulation_data%motion%temperature%value, ' K'
       Write (iunit, '(2x,a,f6.2,a)') 'tau  =  ', simulation_data%motion%relax_time_thermostat%value, ' fs'
       If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover') Then
         Write (iunit, '(2x,a)') 'nchain  =  3    # number of thermostats in the Nose-Hoover chain'  
         Write (iunit, '(2x,a)') 'nstep   =  20   # number of substeps used to integrate the equation of motion&
                                          & of the Nose-Hoover coordinates'
       Else If (Trim(simulation_data%motion%thermostat%type) == 'langevin') Then
         Write (iunit, '(2x,a)') 'damp  =  0.2    # Langevin dumping parameter'
       Else If (Trim(simulation_data%motion%thermostat%type) == 'andersen') Then
         Write (iunit, '(2x,a)') 'mix  =  1.0     # collision amplitude of the Andersen thermostat'
       End If
       Write (iunit, '(a)') '%endblock thermostat'   
     End If
   End If
   
  End Subroutine print_onetep_motion

  Subroutine advise_motion_onetep(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about motion settings
    ! for ONETEP
    !
    ! author    - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data
    
    Character(Len=256)  :: messages(9), header
    Character(Len=256)  :: in_extra

    Logical :: print_header

    in_extra='using the &extra_directives block'
    print_header=.True.
    

    If (Trim(simulation_data%simulation%type) == 'md') Then
      Write (header, '(1x,a)')  'Regarding the MD settings, the user should consider:'
      If (Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
          Call info(' ', 1)
          Write (messages(1), '(1x,a)')  ' - changing "tau" using ALC_EQCM directive "relax_time_thermostat"'
          Call print_warnings(header, print_header, messages, 1)  
          Write (messages(1), '(2x,a)')  'Unfortunately, due to the structure of "%block thermostat", changes to the following&
                                   & directives, if required, must be applied manually (block &extra_directives cannot be used)'
        If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover') Then
          Write (messages(2), '(1x,a)')  ' - "nchain" (number of thermostats in the Nose-Hoover chain)'
          Write (messages(3), '(1x,a)')  ' - "nstep"  (number of substeps used to integrate the equation of motion&
                                        & of the Nose-Hoover coordinates)'
          Call print_warnings(header, print_header, messages, 3)
        Else If (Trim(simulation_data%motion%thermostat%type) == 'langevin' ) Then
          Write (messages(2), '(1x,a)')  ' - "damp" (Langevin dumping parameter)'
          Call print_warnings(header, print_header, messages, 2)
        Else If (Trim(simulation_data%motion%thermostat%type) == 'andersen' ) Then
          Write (messages(2), '(1x,a)')  ' - "mix" (collision amplitude of the Andersen thermostat)'
          Call print_warnings(header, print_header, messages, 2)
        End If
      End If
    End If
    
  End Subroutine advise_motion_onetep
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Solvation  
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine define_solvation_onetep(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check solvation related directives for the generation of input
    ! files for atomistic level simulations with implicit solvent
    !
    ! author    - i.scivetti September 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256)  :: messages(6)
    Character(Len=256)  :: message
    Character(Len=256)  :: error_sol
    Logical             :: error 
    Integer(Kind=wi)    :: i, j         
    Real(Kind=wp)       :: mini,maxi   

    error=.False.
    error_sol = '***ERROR in &solvation (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'

    If (simulation_data%motion%ion_steps%value /= 1) Then
       Write (message,'(2(1x,a))') 'ERROR: Forces including the solvent are not yet supported for PBCs. Set&
                                  & "ion_steps" to 1 for a single point simulation.'
       Call info(message, 1)
       Call error_stop(' ')
    End If
    
    ! Automatization (compulsory)
    If (simulation_data%solvation%in_vacuum_first%fread) Then
      If (simulation_data%solvation%in_vacuum_first%fail) Then
        Write (message,'(2(1x,a))') Trim(error_sol), 'Wrong settings for "in_vacuum_first" directive. Either .True. or .False.&
                                  & (.True. is recommended to avoid the burden of reseting the calculation)'
        Call error_stop(message)
      End If
    Else
        Write (message,'(2(1x,a))') Trim(error_sol), 'The user must specify directive "in_vacuum_first" (.True. is recommended)'
        Call error_stop(message)
    End If

    ! Relaxation dielectric model (compulsory)
    Write (messages(2),'(1x,a)')  'Implemented options:'
    Write (messages(3),'(1x,a)')  '- Fixed            (Fixed cavity model, always recommended) '
    Write (messages(4),'(1x,a)')  '- Self_consistent  (Self-consistent cavity model, computational and numerical demanding)'
    If (simulation_data%solvation%cavity_model%fread) Then
      If (simulation_data%solvation%cavity_model%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Wrong settings for "cavity_model" directive.'
        error=.True.
      Else
        If (Trim(simulation_data%solvation%cavity_model%type) /= 'fixed' .And. &
            Trim(simulation_data%solvation%cavity_model%type) /= 'self_consistent') Then
          Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "cavity_model" is not valid'
          error=.True.
        End If
      End If
    Else
        Write (messages(1),'(2(1x,a))') Trim(error_sol), 'The user must specify directive "cavity_model"'
        error=.True.
    End If

    If (error) Then
      Call info(messages,4)
      Call error_stop(' ')
    End If
    
    ! Dielectric function (compulsory)
    Write (messages(2),'(1x,a)')  'Implemented options:'
    Write (messages(3),'(1x,a)')  '- Fattebert-Gygi (Fattebert-Gygi model) '
    Write (messages(4),'(1x,a)')  '- Andreussi      (Andreussi model)'
    Write (messages(5),'(1x,a)')  '- Soft_sphere    (atomic radii based model)'
    If (simulation_data%solvation%dielectric_function%fread) Then
      If (simulation_data%solvation%dielectric_function%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Wrong settings for "dielectric_function" directive.'
        error=.True.
      Else
        If (Trim(simulation_data%solvation%dielectric_function%type) /= 'fattebert-gygi' .And. &
          Trim(simulation_data%solvation%dielectric_function%type) /= 'andreussi'       .And. &
          Trim(simulation_data%solvation%dielectric_function%type) /= 'soft_sphere') Then
          Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "dielectric_function" is not valid.'
          error=.True.     
        End If
      End If
    Else
        Write (messages(1),'(2(1x,a))') Trim(error_sol), 'The user must specify directive "dielectric_function"'
        error=.True.
    End If

    If (error) Then
      Call info(messages,5)
      Call error_stop(' ')
    End If

    ! Set the biblio
    If (Trim(simulation_data%solvation%dielectric_function%type) == 'fattebert-gygi') Then
      simulation_data%solvation%bib_epsilon=bib_fg
    Else If(Trim(simulation_data%solvation%dielectric_function%type) == 'andreussi') Then
      simulation_data%solvation%bib_epsilon=bib_andreussi
    Else If(Trim(simulation_data%solvation%dielectric_function%type) == 'soft_sphere') Then
      simulation_data%solvation%bib_epsilon=bib_fisicaro
    End If
     
    ! Important condition 
    If(Trim(simulation_data%solvation%cavity_model%type) == 'self_consistent' .And. &
       Trim(simulation_data%solvation%dielectric_function%type) /= 'soft_sphere') Then
       Write (messages(1),'(2(1x,a))') Trim(error_sol), '"Soft_sphere" approximation is incompatible with the& 
                                      & self_consistent model. Please review the settings.'
       Call info(messages, 1)
       Call error_stop(' ')
    End If   

    ! Permittivity of bulk
    If (simulation_data%solvation%permittivity_bulk%fread) Then
      If (simulation_data%solvation%permittivity_bulk%fail) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "permittivity_bulk" is not valid'
         Call info(messages, 1)
         Call error_stop(' ')
      Else
         If(simulation_data%solvation%permittivity_bulk%value < 1.0_wp) Then
           Write (messages(1),'(2(1x,a))') Trim(error_sol), '"permittivity_bulk" for the implicit solvent must be larger than 1.0'
           Call info(messages, 1)
           Call error_stop(' ')
         End If
      End If
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Directive "permittivity_bulk" must be defined'
      Call info(messages, 1)
      Call error_stop(' ')
    End If

    ! Smear_ion_width
    If (simulation_data%solvation%smear_ion_width%fread) Then
      If (simulation_data%solvation%smear_ion_width%fail) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "smear_ion_width" is not valid'
         Call info(messages, 1)
         Call error_stop(' ')
      Else
         If((Abs(simulation_data%solvation%smear_ion_width%value) < Abs(epsilon(1.0_wp))) .Or. &
            (simulation_data%solvation%smear_ion_width%value < epsilon(1.0_wp))) Then
           Write (messages(1),'(2(1x,a))') Trim(error_sol), '"smear_ion_width" must be a positive real'
           Call info(messages, 1)
           Call error_stop(' ')
         End If
         If (Trim(simulation_data%solvation%smear_ion_width%units) /= 'bohr') Then
           Write (message,'(2(1x,a))') Trim(error_sol), &
                                    &'Units for "smear_ion_width" must be in Bohr. Please change'
           Call error_stop(message)
         End If
      End If
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Directive "smear_ion_width" must be defined&
                                    & (recommended value is 0.8 Bohr)'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
    ! Density threshold
    If (simulation_data%solvation%density_threshold%fread) Then
      If (simulation_data%solvation%dielectric_function%type /= 'fattebert-gygi') Then
          Write (messages(1),'(1x,4a)') Trim(error_sol), ' Directive "density_threshold" is NOT needed when&  
                                   & "dielectric_function" is set to "', Trim(simulation_data%solvation%dielectric_function%type), &
                                   & '". Please review the settings.'
         Call info(messages, 1)
         Call error_stop(' ')

      End If

      If (simulation_data%solvation%density_threshold%fail) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "density_threshold" is not valid'
         Call info(messages, 1)
         Call error_stop(' ')
      Else
         If((Abs(simulation_data%solvation%density_threshold%value) < Abs(epsilon(1.0_wp))) .Or. &
            (simulation_data%solvation%density_threshold%value < epsilon(1.0_wp))) Then
           Write (messages(1),'(2(1x,a))') Trim(error_sol), '"density_threshold" must be a positive real'
           Call info(messages, 1)
           Call error_stop(' ')
         End If
      End If
    Else
       If (simulation_data%solvation%dielectric_function%type == 'fattebert-gygi') Then
         simulation_data%solvation%density_threshold%value=0.00035_wp
       End If
    End If

    ! Beta
    If (simulation_data%solvation%beta%fread) Then
      If (simulation_data%solvation%dielectric_function%type == 'fattebert-gygi') Then
          Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Directive "beta" is only needed when&  
                                        & "dielectric_function" is set to "Fattebert-Gygi". Please review the settings.'
          Call info(messages, 1)
          Call error_stop(' ')
      End If

      If (simulation_data%solvation%beta%fail) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "beta" is not valid'
         Call info(messages, 1)
         Call error_stop(' ')
      Else
         If((Abs(simulation_data%solvation%beta%value) < Abs(epsilon(1.0_wp))) .Or. &
            (simulation_data%solvation%beta%value < epsilon(1.0_wp))) Then
           Write (messages(1),'(2(1x,a))') Trim(error_sol), '"beta" must be a positive real'
           Call info(messages, 1)
           Call error_stop(' ')
         End If
      End If
    Else
       simulation_data%solvation%beta%value=1.3_wp
    End If
    
    ! Minimum density threshold
    If (simulation_data%solvation%density_min_threshold%fread) Then
      If (simulation_data%solvation%dielectric_function%type /= 'andreussi') Then
          Write (messages(1),'(1x,4a)') Trim(error_sol), ' Directive "density_min_threshold" is NOT needed when&  
                                   & "dielectric_function" is set to "', Trim(simulation_data%solvation%dielectric_function%type), &
                                   & '". Please review the settings.'
          Call info(messages, 1)
          Call error_stop(' ')
      End If

      If (simulation_data%solvation%density_min_threshold%fail) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "density_min_threshold" is not valid'
         Call info(messages, 1)
         Call error_stop(' ')
      Else
         If((Abs(simulation_data%solvation%density_min_threshold%value) < Abs(epsilon(1.0_wp))) .Or. &
            (simulation_data%solvation%density_min_threshold%value < epsilon(1.0_wp))) Then
           Write (messages(1),'(2(1x,a))') Trim(error_sol), '"density_min_threshold" must be a positive real'
           Call info(messages, 1)
           Call error_stop(' ')
         End If
      End If
    Else
      If (simulation_data%solvation%dielectric_function%type == 'andreussi') Then 
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Directive "density_min_threshold" must be defined for the "Andreussi"&
                                       & method.'
         Call info(messages, 1)
         Call error_stop(' ')
      End If
    End If
        
    ! Maximum density threshold
    If (simulation_data%solvation%density_max_threshold%fread) Then
      If (simulation_data%solvation%dielectric_function%type /= 'andreussi') Then
          Write (messages(1),'(1x,4a)') Trim(error_sol), ' Directive "density_max_threshold" is NOT needed when&  
                                   & "dielectric_function" is set to "', Trim(simulation_data%solvation%dielectric_function%type), &
                                   & '". Please review the settings.'
          Call info(messages, 1)
          Call error_stop(' ')
      End If

      If (simulation_data%solvation%density_max_threshold%fail) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "density_max_threshold" is not valid'
         Call info(messages, 1)
         Call error_stop(' ')
      Else
         If((Abs(simulation_data%solvation%density_max_threshold%value) < Abs(epsilon(1.0_wp))) .Or. &
            (simulation_data%solvation%density_max_threshold%value < epsilon(1.0_wp))) Then
           Write (messages(1),'(2(1x,a))') Trim(error_sol), '"density_max_threshold" must be a positive real'
           Call info(messages, 1)
           Call error_stop(' ')
         End If
      End If
    Else
      If (simulation_data%solvation%dielectric_function%type == 'andreussi') Then 
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Directive "density_max_threshold" must be defined for the "Andreussi"&
                                       & method.'
         Call info(messages, 1)
         Call error_stop(' ')
      End If
    End If

    mini=simulation_data%solvation%density_min_threshold%value
    maxi=simulation_data%solvation%density_max_threshold%value

    If (Trim(simulation_data%solvation%dielectric_function%type) == 'andreussi') Then 
      If (maxi < mini .Or. Abs(mini-maxi) < epsilon(1.0_wp)) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), '"density_max_threshold" must be larger than "density_min_threshold"'
         Call info(messages, 1)
         Call error_stop(' ')
      End If
    End If

    ! Test for "&soft_sphere_radii
    If (simulation_data%solvation%soft_radii_info%stat .And. &
        Trim(simulation_data%solvation%dielectric_function%type) /= 'soft_sphere') Then
        Write (messages(1),'(1x,4a)') Trim(error_sol), ' Block "&soft_sphere_radii" MUST NOT be defined when&  
                                   & "dielectric_function" is set to "', Trim(simulation_data%solvation%dielectric_function%type), &
                                   & '". Please review the settings.'
        Call info(messages, 1)
        Call error_stop(' ')
    End If

    ! Check settings of block &soft_sphere_radii
    If (simulation_data%solvation%soft_radii_info%stat) Then
      ! Check if user has included all the tags
      Do i=1, simulation_data%total_tags
        error=.True.
        Do j=1, simulation_data%total_tags
          If (Trim(simulation_data%solvation%soft_radii(i)%tag)==Trim(simulation_data%component(j)%tag)) Then
            simulation_data%solvation%soft_radii(i)%element=simulation_data%component(j)%element
            error=.False.
          End If
        End Do
        If (error) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_sol),&
                                   & 'Atomic tag "', Trim(simulation_data%solvation%soft_radii(i)%tag), '" declared in&
                                   & "&soft_sphere_radii" has not been defined in "&block_input_composition". Please check'
          Call error_stop(message)
        End If
      End Do

      ! Check there is no negative not zero value
      Do i=1, simulation_data%total_tags
        If (simulation_data%solvation%soft_radii(i)%value < 0.0_wp .Or. &
            Abs(simulation_data%solvation%soft_radii(i)%value) < epsilon(1.0_wp)) Then
          Write (message,'(1x,a,1x,3a)') Trim(error_sol),&
                                   & 'Value of  soft_sphere_radius to atomic tag "', Trim(simulation_data%motion%mass(i)%tag),&
                                   & '" must be positive and different from zero (check block &soft_sphere_radii)'
          Call error_stop(message)
        End If
      End Do

    End If    
    
    ! Soft sphere scale
    If (simulation_data%solvation%soft_sphere_scale%fread) Then
      If (simulation_data%solvation%dielectric_function%type /= 'soft_sphere') Then
          Write (messages(1),'(1x,4a)') Trim(error_sol), ' Directive "soft_sphere_scale" is NOT needed when&  
                                   & "dielectric_function" is set to "', Trim(simulation_data%solvation%dielectric_function%type), &
                                   & '". Please review the settings.'
         Call info(messages, 1)
         Call error_stop(' ')
      Else 
         If (simulation_data%solvation%soft_radii_info%stat) Then
          Write (messages(1),'(1x,3a)') Trim(error_sol), ' Directive "soft_sphere_scale" does NOT apply when&  
                                   & sphere radii are defined with sub-block &soft_sphere_radii and MUST be removed.&
                                   & Please review the settings.'
         Call info(messages, 1)
         Call error_stop(' ')
         End If
      End If

      If (simulation_data%solvation%soft_sphere_scale%fail) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "soft_sphere_scale" is not valid'
         Call info(messages, 1)
         Call error_stop(' ')
      Else
         If((Abs(simulation_data%solvation%soft_sphere_scale%value) < Abs(epsilon(1.0_wp))) .Or. &
            (simulation_data%solvation%soft_sphere_scale%value < epsilon(1.0_wp))) Then
           Write (messages(1),'(2(1x,a))') Trim(error_sol), '"soft_sphere_scale" must be a positive real'
           Call info(messages, 1)
           Call error_stop(' ')
         End If
      End If
    Else
       If (simulation_data%solvation%dielectric_function%type == 'soft_sphere') Then
         simulation_data%solvation%soft_sphere_scale%value=1.33_wp
       End If
    End If

    ! Soft sphere radii delta
    If (simulation_data%solvation%soft_sphere_delta%fread) Then
      If (simulation_data%solvation%dielectric_function%type /= 'soft_sphere') Then
          Write (messages(1),'(1x,4a)') Trim(error_sol), ' Directive "soft_sphere_delta" is NOT needed when&  
                                   & "dielectric_function" is set to "', Trim(simulation_data%solvation%dielectric_function%type), &
                                   & '" and must be removed. Please review the settings.'
         Call info(messages, 1)
         Call error_stop(' ')
      End If

      If (simulation_data%solvation%soft_sphere_delta%fail) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "soft_sphere_delta" is not valid'
         Call info(messages, 1)
         Call error_stop(' ')
      Else
         If((Abs(simulation_data%solvation%soft_sphere_delta%value) < Abs(epsilon(1.0_wp))) .Or. &
            (simulation_data%solvation%soft_sphere_delta%value < epsilon(1.0_wp))) Then
           Write (messages(1),'(2(1x,a))') Trim(error_sol), '"soft_sphere_delta" must be a positive real'
           Call info(messages, 1)
           Call error_stop(' ')
         End If
      End If
    Else
       If (simulation_data%solvation%dielectric_function%type == 'soft_sphere') Then
         simulation_data%solvation%soft_sphere_delta%value=0.5_wp
       End If
    End If
  
    ! Apolar terms (compulsory)
    Write (messages(2),'(1x,a)')  'Implemented options:'
    Write (messages(3),'(1x,a)')  '- SASA              (Solvent-Accessible Surface-Area approximation) '
    Write (messages(4),'(1x,a)')  '- SAV               (Surface-Accessible Volume approximation)'
    Write (messages(5),'(1x,a)')  '- Only_cavitation   (only cavitation terms)'
    Write (messages(6),'(1x,a)')  '- None              (No apolar corrections)'
    If (simulation_data%solvation%apolar_terms%fread) Then
      If (simulation_data%solvation%apolar_terms%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Wrong settings for "apolar_terms" directive.'
        error=.True.
      Else
          If (Trim(simulation_data%solvation%apolar_terms%type) /= 'sasa' .And. &
              Trim(simulation_data%solvation%apolar_terms%type) /= 'sav'       .And. &
              Trim(simulation_data%solvation%apolar_terms%type) /= 'only_cavitation'       .And. &
              Trim(simulation_data%solvation%apolar_terms%type) /= 'none') Then
              Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "apolar_terms" is not valid'
            error=.True.
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))') Trim(error_sol), 'The user must specify directive "apolar_terms"'
       error=.True.
    End If

    If (error) Then
      Call info(messages, 6)
      Call error_stop(' ')
    End If
    
    ! SASA definiton
    Write (messages(2),'(1x,a)')  'Implemented options:'
    Write (messages(3),'(1x,a)')  '- Density           '
    Write (messages(4),'(1x,a)')  '- Isosurface        '
    If (simulation_data%solvation%sasa_definition%fread) Then
      If (simulation_data%solvation%sasa_definition%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Wrong settings for "sasa_definition" directive.'
        error=.True.
      Else
          If (Trim(simulation_data%solvation%sasa_definition%type) /= 'density' .And. &
              Trim(simulation_data%solvation%sasa_definition%type) /= 'isosurface') Then
              Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "sasa_definition" is not valid'
            error=.True.
          End If
      End If
    Else
       If (Trim(simulation_data%solvation%dielectric_function%type) /= 'soft_sphere') Then
         If (Trim(simulation_data%solvation%apolar_terms%type) /= 'none') Then
           Write (messages(1),'(2(1x,a))') Trim(error_sol), 'The user must specify directive "sasa_definition".&
                                         & Please refer to the ONETEP manual'
           error=.True.
         End If
       End If
    End If

    If (error) Then
      Call info(messages, 4)
      Call error_stop(' ')
    End If

    If (Trim(simulation_data%solvation%sasa_definition%type) == 'isosurface' .And. &
        Trim(simulation_data%solvation%dielectric_function%type) == 'andreussi' ) Then
        Write (messages(1),'(2(1x,a))') Trim(error_sol), 'the only possible option of "sasa_definition" for the "andreussi"&
                                       & method is "density". Please review the settings.'
        Call info(messages, 1)
        Call error_stop(' ')
    End If


    ! Apolar scaling factor
    If (simulation_data%solvation%apolar_scaling%fread) Then
      If (Trim(simulation_data%solvation%apolar_terms%type) == 'none') Then
          Write (messages(1),'(1x,2a)') Trim(error_sol), ' Directive "apolar_scaling" is NOT needed when&  
                                   & "apolar_terms" is set to "None" and must be removed. Please review the settings.'
         Call info(messages, 1)
         Call error_stop(' ')
      End If

      If (simulation_data%solvation%apolar_scaling%fail) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'Specification of "apolar_scaling" is not valid'
         Call info(messages, 1)
         Call error_stop(' ')
      Else
         If((Abs(simulation_data%solvation%apolar_scaling%value) < Abs(epsilon(1.0_wp))) .Or. &
            (simulation_data%solvation%apolar_scaling%value < epsilon(1.0_wp))) Then
           Write (messages(1),'(2(1x,a))') Trim(error_sol), '"apolar_scaling" must be a positive real'
           Call info(messages, 1)
           Call error_stop(' ')
         End If
      End If
    Else
      If (Trim(simulation_data%solvation%apolar_terms%type) /= 'none') Then
          Write (messages(1),'(2(1x,a))') Trim(error_sol), 'The user must define "apolar_scaling" (positive real). It must be&
                                        & equal to 1.0 if the "only_cavitation" option is set for "apolar_terms". See manual'
          Call info(messages, 1)
          Call error_stop(' ')
      End If
    End If

    If (Trim(simulation_data%solvation%apolar_terms%type) == 'only_cavitation') Then
       If (Abs(1.0_wp-simulation_data%solvation%apolar_scaling%value) > Abs(epsilon(1.0_wp))) Then
         Write (messages(1),'(2(1x,a))') Trim(error_sol), 'For the "only_cavitation" option, "apolar_scaling" must be equal to 1.0'
         Call info(messages, 1)
         Call error_stop(' ')
       End If
    End If

    ! Solvent Pressure
    If (simulation_data%solvation%solvent_pressure%fread) Then
      If (simulation_data%solvation%solvent_pressure%fail) Then
        Write (message,'(2(1x,a))') Trim(error_sol), 'Wrong (or missing) settings for "solvent_pressure" directive.&
                                & Both value and units are required. See manual' 
        Call error_stop(message)
      Else
        If (Trim(simulation_data%solvation%apolar_terms%type) /= 'sav') Then
          Write (message,'(1x,4a)') Trim(error_sol), ' Definition of "solvent_pressure" is incompatible with the "',&
                                & Trim(simulation_data%solvation%apolar_terms%type), '" option for "apolar_terms".&
                                & Please review settings.' 
          Call error_stop(message)
        End If
        If (Trim(simulation_data%solvation%solvent_pressure%units) /= 'gpa') Then
          Write (message,'(2(1x,a))') Trim(error_sol), &
                                    &'Units for "solvet_pressure" must be in GPa. Please change'
          Call error_stop(message)
        End If
      End If
    Else
       If (Trim(simulation_data%solvation%apolar_terms%type) == 'sav') Then
          Write (message,'(2(1x,a))') Trim(error_sol), 'The user must define "solvent_pressure", needed for the SAV method.' 
          Call error_stop(message)
       End If
    End If    
    
    ! Solvent surface tension
    If (simulation_data%solvation%surf_tension%fread) Then
      If (simulation_data%solvation%surf_tension%fail) Then
        Write (message,'(2(1x,a))') Trim(error_sol), 'Wrong (or missing) settings for "surf_tension" directive.&
                                & Both value and units are required. See manual' 
        Call error_stop(message)
      Else
        If (Trim(simulation_data%solvation%apolar_terms%type) == 'none') Then
          Write (message,'(1x,2a)') Trim(error_sol), ' Definition of "surf_tension" is not required when "apolar_terms"&
                                & is set to "None". Please review settings and change' 
          Call error_stop(message)
        End If
        If (Trim(simulation_data%solvation%surf_tension%units(1)) /= 'n' .Or. &
            Trim(simulation_data%solvation%surf_tension%units(2)) /= 'm-1') Then
          Write (message,'(2(1x,a))') Trim(error_sol), &
                                    &'Units for "surf_tension" must be in "N m-1". Please change this input to comply with&
                                    & the required units.'
          Call error_stop(message)
        End If
      End If
    Else
       If (Trim(simulation_data%solvation%apolar_terms%type) /= 'none') Then
          Write (message,'(2(1x,a))') Trim(error_sol), 'The user must define "surf_tension" to set the computation of&
                                     & the cavitation term.' 
          Call error_stop(message)
       End If
    End If    
    
  End Subroutine define_solvation_onetep

  Subroutine print_onetep_solvation(iunit, ic, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print solvation directives for ONETEP
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Integer(Kind=wi), Intent(InOut) :: ic    
    Type(simul_type), Intent(InOut) :: simulation_data

    Character(Len=256) :: message
    Integer(Kind=wi)   :: j
    
     Write (iunit,'(a)') ' '
     Write (iunit,'(a)') '#==== Solvation with implicit solvent '
     Write (message,'(a)') 'is_implicit_solvent :  T  # Turns on the implicit solvent'
     Call record_directive(iunit, message, 'is_implicit_solvent', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a)') 'is_smeared_ion_rep :  T  # Turns on the smear-ion representation, consistent with solvation'
     Call record_directive(iunit, message, 'is_smeared_ion_rep', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,f5.2)') 'is_smeared_ion_width : ', simulation_data%solvation%smear_ion_width%value
     Call record_directive(iunit, message, 'is_smeared_ion_width', simulation_data%set_directives%array(ic), ic)

     If (simulation_data%solvation%in_vacuum_first%stat) Then
       Write (message,'(a)') 'is_auto_solvation :    T   # Solvation will follow an initial simulation in vacuum'
     Else
       Write (message,'(a)') 'is_auto_solvation :    F   # Solvation will start from scratch'
     End If
     Call record_directive(iunit, message, 'is_auto_solvation', simulation_data%set_directives%array(ic), ic)
     If (Trim(simulation_data%solvation%cavity_model%type) == 'fixed') Then
       Write (message,'(a)') 'is_dielectric_model :  fix_initial   # the cavity will be fixed' 
     Else If (Trim(simulation_data%solvation%cavity_model%type) == 'self_consistent') Then
       Write (message,'(a)') 'is_dielectric_model :  self_consistent # the cavity will change its shape'
     End If
     Call record_directive(iunit, message, 'is_dielectric_model', simulation_data%set_directives%array(ic), ic)
     Write (iunit,'(a)') '# Details for the dielectric model'
     If (Trim(simulation_data%solvation%dielectric_function%type) == 'fattebert-gygi') Then
       Write (message,'(2a)') 'is_dielectric_function :  FGF    # ', Trim(simulation_data%solvation%bib_epsilon)
     Else 
       Write (message,'(2a,2(1x,a))') 'is_dielectric_function : ', Trim(simulation_data%solvation%dielectric_function%type), '#',&
                                  &  Trim(simulation_data%solvation%bib_epsilon)
     End If        
     Call record_directive(iunit, message, 'is_dielectric_function', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,f7.3)') 'is_bulk_permittivity : ', simulation_data%solvation%permittivity_bulk%value
     Call record_directive(iunit, message, 'is_bulk_permittivity', simulation_data%set_directives%array(ic), ic)

     If (Trim(simulation_data%solvation%dielectric_function%type) /= 'soft_sphere')Then
       If (Trim(simulation_data%solvation%dielectric_function%type) == 'fattebert-gygi' .Or. &
           Trim(simulation_data%solvation%dielectric_function%type) == 'andreussi') Then
         Write (message,'(a,f7.3)') 'is_solvation_beta : ', simulation_data%solvation%beta%value
         Call record_directive(iunit, message, 'is_solvation_beta', simulation_data%set_directives%array(ic), ic)
       End If  
       If (Trim(simulation_data%solvation%dielectric_function%type) == 'andreussi') Then
         Write (message,'(a,f10.6)') 'is_density_min_threshold : ', simulation_data%solvation%density_min_threshold%value
         Call record_directive(iunit, message, 'is_density_min_threshold', simulation_data%set_directives%array(ic), ic)
         Write (message,'(a,f10.6)') 'is_density_max_threshold : ', simulation_data%solvation%density_max_threshold%value
         Call record_directive(iunit, message, 'is_density_max_threshold', simulation_data%set_directives%array(ic), ic)
       Else If (Trim(simulation_data%solvation%dielectric_function%type) == 'fattebert-gygi') Then 
         Write (message,'(a,f10.6)') 'is_density_threshold : ', simulation_data%solvation%density_threshold%value
         Call record_directive(iunit, message, 'is_density_threshold', simulation_data%set_directives%array(ic), ic)
       End If
     Else
       Write (message,'(a,f10.6)') 'is_soft_sphere_delta : ', simulation_data%solvation%soft_sphere_delta%value
       Call record_directive(iunit, message, 'is_soft_sphere_delta', simulation_data%set_directives%array(ic), ic)
       If (simulation_data%solvation%soft_radii_info%stat) Then
         Write (iunit,'(a)') '%block is_soft_sphere_radii'
         Do j=1, simulation_data%total_tags
               Write (iunit,'(1x,a3,4x,f8.3)') Trim(simulation_data%solvation%soft_radii(j)%tag), &
                                            & simulation_data%solvation%soft_radii(j)%value
         End Do
         Write (iunit,'(a)') '%endblock is_soft_sphere_radii'
         Write (iunit,'(a)') ' '
       Else
         Write (message,'(a,f10.6)') 'is_soft_sphere_scale : ', simulation_data%solvation%soft_sphere_scale%value
         Call record_directive(iunit, message, 'is_soft_sphere_scale', simulation_data%set_directives%array(ic), ic)
       End If
     End If  
     
     If (Trim(simulation_data%solvation%apolar_terms%type) == 'none') Then
       Write (iunit,'(a)') '# Apolar terms will be ommited'
       Write (message,'(a)') 'is_include_apolar :  F'
     Else
       Write (iunit,'(a)') '# Settings for apolar corrections'
       Write (message,'(a)') 'is_include_apolar :  T'
     End If
     Call record_directive(iunit, message, 'is_include_apolar', simulation_data%set_directives%array(ic), ic)
     
     If (Trim(simulation_data%solvation%apolar_terms%type) /= 'none') Then
       If (Trim(simulation_data%solvation%apolar_terms%type) == 'only_cavitation') Then
          Write (message,'(a)')        'is_apolar_scaling_factor : 1.0 (only cavitation)'
          Call record_directive(iunit, message, 'is_apolar_scaling_factor', simulation_data%set_directives%array(ic), ic)
       Else
          Write (message,'(2a)') 'is_apolar_method : ', Trim(simulation_data%solvation%apolar_terms%type)
          Call record_directive(iunit, message, 'is_apolar_method', simulation_data%set_directives%array(ic), ic)
          Write (message,'(a,f10.6)') 'is_apolar_scaling_factor : ', simulation_data%solvation%apolar_scaling%value
          Call record_directive(iunit, message, 'is_apolar_scaling_factor', simulation_data%set_directives%array(ic), ic)
       End If

       If (Trim(simulation_data%solvation%dielectric_function%type) /= 'soft_sphere')Then 
         Write (message,'(2a)') 'is_apolar_sasa_definition : ', Trim(simulation_data%solvation%sasa_definition%type)
         Call record_directive(iunit, message, 'is_apolar_sasa_definition', simulation_data%set_directives%array(ic), ic) 
       End If

       If (Trim(simulation_data%solvation%apolar_terms%type) == 'sasa' .Or. &
           Trim(simulation_data%solvation%apolar_terms%type) == 'sav') Then
            Write (message,'(a,f8.4,a)') 'is_solvent_surf_tension : ', simulation_data%solvation%surf_tension%value, '   N/m' 
            Call record_directive(iunit, message, 'is_solvent_surf_tension', simulation_data%set_directives%array(ic), ic)
          If (Trim(simulation_data%solvation%apolar_terms%type) == 'sav') Then
            Write (message,'(a,f8.2,a)') 'is_solvent_pressure : ', simulation_data%solvation%solvent_pressure%value, '   GPa' 
            Call record_directive(iunit, message, 'is_solvent_pressure', simulation_data%set_directives%array(ic), ic)
          End If
       End If
     End If

     Write (iunit,'(a)') '# PBCs for solvation'
     Write (message,'(a)') 'multigrid_bc :   P  P  P'
     Call record_directive(iunit, message, 'multigrid_bc', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a)') 'ion_ion_bc :     P  P  P'
     Call record_directive(iunit, message, 'ion_ion_bc', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a)') 'pspot_bc :       P  P  P'
     Call record_directive(iunit, message, 'pspot_bc', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a)') 'smeared_ion_bc : P  P  P'
     Call record_directive(iunit, message, 'smeared_ion_bc', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a)') 'vdw_bc :         P  P  P'
     Call record_directive(iunit, message, 'vdw_bc', simulation_data%set_directives%array(ic), ic)
     
  End Subroutine print_onetep_solvation
  
  Subroutine summary_solvation_onetep(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to summarise solvation settings from the information
    ! provided by the user via block &solvation
    !
    ! author    - i.scivetti September 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(7)

    If (simulation_data%solvation%in_vacuum_first%fread) Then
      Write (messages(1),'(1x,a)')  ' - a calculation in vacuum will be run before including the solvent' 
    Else
     Write (messages(1),'(1x,a)')  '- the simulation including the solvent will be executed without any&
                                   & previous calculation in vacuum' 
    end If
    Write (messages(2),'(2(1x,a))')  ' - type of solvent cavity: ', Trim(simulation_data%solvation%cavity_model%type)
    Write (messages(3),'(1x,a)')     ' - PBCs will be set for the solvent' 
    Write (messages(4),'(1x,a)')     '= Details for the dielectric model'
    Write (messages(5),'(2(1x,a))')     ' - method: ', Trim(simulation_data%solvation%dielectric_function%type)
    Write (messages(6),'(1x,a,f10.6)')  ' - permittivity of the bulk solvent: ', simulation_data%solvation%permittivity_bulk%value
    Write (messages(7),'(1x,a,f5.2,a)') ' - smear width for ions:  ', simulation_data%solvation%smear_ion_width%value, '  Bohr' 
    Call info(messages, 7) 

    If (Trim(simulation_data%solvation%dielectric_function%type) /= 'soft_sphere')Then
        Write (messages(1),'(1x,a,f10.6)') ' - beta factor: ', simulation_data%solvation%beta%value
       Call info(messages, 1)

       If (Trim(simulation_data%solvation%dielectric_function%type) == 'fattebert-gygi') Then
         Write (messages(1),'(1x,a,f10.6)') ' - density threshold: ', simulation_data%solvation%density_threshold%value
         Call info(messages, 1)
       Else If (Trim(simulation_data%solvation%dielectric_function%type) == 'andreussi') Then
         Write (messages(1),'(1x,a,f10.6)') ' - minimum density threshold: ', simulation_data%solvation%density_min_threshold%value
         Write (messages(2),'(1x,a,f10.6)') ' - maximum density threshold: ', simulation_data%solvation%density_max_threshold%value
         Call info(messages, 2)
       End If
    Else
      Write (messages(1),'(1x,a,f10.6)') ' - steepness of soft sphere transition : ', &
                                       & simulation_data%solvation%soft_sphere_delta%value
      If (simulation_data%solvation%soft_radii_info%stat) Then
        Write (messages(2),'(1x,a)') ' - soft sphere radii are defined in sub-block &soft_sphere_radii'
        Call info(messages, 2)
      Else
        Write (messages(2),'(1x,a)') ' - soft sphere radii will be set equal to the vdW radii scaled by a factor (see next line)'
        Write (messages(3),'(1x,a,f10.6)') ' - scaling of soft sphere radii: ', &
                                         & simulation_data%solvation%soft_sphere_scale%value
        Call info(messages, 3)
      End If
    End If  

    If (Trim(simulation_data%solvation%apolar_terms%type) /= 'none') Then
      Write (messages(1),'(1x,a)')     '= Apolar corrections'
      Write (messages(2),'(1x,2a)')    ' - approximation: ', Trim(simulation_data%solvation%apolar_terms%type)
      Call info(messages, 2)

      If (simulation_data%solvation%sasa_definition%fread) Then
          Write (messages(1),'(1x,3a)') ' - method to compute the solvent accessible surface area: ', &
                                       & Trim(simulation_data%solvation%sasa_definition%type) 
          Call info(messages, 1) 
      End If 

      Write (messages(1),'(1x,a,f8.4,a)') ' - surface tension for solvent: ',&
                                         & simulation_data%solvation%surf_tension%value, '   N/m' 
      Call info(messages, 1)

      If (Trim(simulation_data%solvation%apolar_terms%type) == 'sav') Then
        Write (messages(1),'(1x,a,f8.4,a)') ' - solvent pressure: ', simulation_data%solvation%solvent_pressure%value, '   GPa&
                                         & (not a physical pressure but a volumetric correction to solvation with the SAV method)'
        Call info(messages, 1)
      End If
    Else
      Write (messages(1),'(1x,a)')  ' - ommiting apolar corrections to the solvation energy' 
      Call info(messages, 1) 
    End If

  End Subroutine summary_solvation_onetep

  
  Subroutine advise_solvation_onetep(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about ONETEP settings
    ! for solvation
    !
    ! author    - i.scivetti August 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(9)

    Write (messages(1), '(1x,a)') '====== Hints for implicit solvent'
    Call info(messages, 1)
    
    If (.Not. simulation_data%solvation%in_vacuum_first%stat) Then
      Write (messages(1), '(1x,a)') '- solvation simulation from scratch (without a previous simulation in vacuum)&
                                   & is NOT recommended. It is ALWAYS SAFER to set keyword "in_vacuum_first" to .True.'
      Call info(messages, 1)                             
    End If
    If (Trim(simulation_data%solvation%cavity_model%type) == 'self_consistent') Then
       Write (messages(1), '(1x,a)') '- simulations with the self consistent model for the cavity are often numerical and&
                                    & computational very challenging. For large systems, it is recommended to set the&
                                    & "fixed" option for "cavity_model"'
       Call info(messages, 1)                             
    End If

    If (.Not. simulation_data%solvation%permittivity_bulk%fread) Then
      Write (messages(1), '(1x,a)') '- the value of the "permittivity_bulk" for the solvent is set by default.&
                                   & Make sure to set the correct value'
       Call info(messages, 1)                             
    End If

    If (Trim(simulation_data%solvation%dielectric_function%type) == 'fattebert-gygi') Then
      If (.Not. simulation_data%solvation%density_threshold%fread) Then
        
         Write (messages(1), '(1x,a)') '- the value of the "density_threshold" for the Fattebert-Gygi function is set by default.&
                                   & Verify the value is appropriate to the purpose of the simulation'
         Call info(messages, 1)
      End If   
      If (.Not. simulation_data%solvation%beta%fread) Then
        Write (messages(1), '(1x,a)') '- the value of the "beta" is set by default.&
                                   & Check the value is appropriate to the purpose of the simulation'
        Call info(messages, 1)                             
      End If

    End If

    If (Trim(simulation_data%solvation%dielectric_function%type) == 'andreussi') Then
      Write (messages(1), '(1x,a)') '- make sure to optimise the values of "density_min_threshold" and&
                                   & "density_max_threshold" for the solvent under consideration'
      Call info(messages, 1)                             
    End If
    
    
    If (Trim(simulation_data%solvation%dielectric_function%type) == 'soft_sphere')Then
      If (simulation_data%solvation%soft_radii_info%stat) Then
        Write (messages(1), '(1x,a)') '- check if the radii defined in &soft_sphere_radii for the involved species&
                                     & are correct. Alternatively, the user can remove the sub-block and use vdW parameters'
      Else
        If (simulation_data%solvation%soft_sphere_scale%fread) Then
          Write (messages(1), '(1x,a)') '- the value of the "soft_sphere_scale" to scale the vdW radii for the "soft_sphere"&
                                        & method is by default'
          Call info(messages, 1)                             
        End If
        Write (messages(1), '(1x,a)') '- the use of scaled vdW radii for the "soft_sphere" method needs the optimization of&
                                    & "soft_sphere_scale". Alternatively, the user can define the &soft_sphere_radii sub-block'
      End If
      Call info(messages, 1)
      
      If (simulation_data%solvation%soft_sphere_delta%fread) Then
        Write (messages(1), '(1x,a)') '- the value of the "soft_sphere_delta" to control the steepness of the dielectric&
                                      & function is set by default'
        Call info(messages, 1)                             
      End If
    End If

    If (Trim(simulation_data%solvation%apolar_terms%type) == 'none') Then
      Write (messages(1), '(1x,a)') '- results from the omission of apolar terms must be interpreted with special care'
      Call info(messages, 1)
    End If

    If (Trim(simulation_data%solvation%cavity_model%type) == 'self_consistent') Then
      Write (messages(1), '(1x,a)') '- for the self_consistent cavity model, set an optimal value for&
                                    & "fine_grid_scale" (equal to 2.0 by default). Use the &extra_directives block'
      Call info(messages, 1)
    End If

    Write (messages(1), '(1x,a)') 'IMPORTANT'
    Write (messages(2), '(1x,a)') 'Computation of ion forces including the implicit solvent are not implemented to date.'
    Write (messages(3), '(1x,a)') 'Advanced solvation directives can be set with the &extra_directives block&
                                  & (only if the user knows what to do, see ONETEP notes).'
    Write (messages(4), '(1x,a)') 'Exclusion of regions from solvent filling is not possible with PBCs.'
                                  
    Call info(messages, 4)
   If (.Not. simulation_data%solvation%both_surfaces) Then
     Call info(' ', 1)
     Write (messages(1), '(1x,a)') '****************************************************************************************'
     Write (messages(2), '(1x,a)') 'ATTENTION!!!'
     Write (messages(3), '(1x,a)') 'Deposited species are located at one side of the slab only. The solvation model will act'
     Write (messages(4), '(1x,a)') 'unevenly. We advise setting the "both_surfaces" directive to .True. (see manual).'    
     Write (messages(5), '(1x,a)') '****************************************************************************************'
     Call info(messages, 5)
     Call info(' ', 1)
   End If
    
  End Subroutine advise_solvation_onetep

  Subroutine advise_multigrid_onetep(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Advise about the DL_MG solver
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(2)

    ! Reference data to compare against &extra_directives
    Type(type_ref_data) :: ref_extra_data(max_directives)
    Character(Len=256)  :: exceptions(10)
    Logical             :: extradir_header 
    Integer(Kind=wi)    :: num_ref_data
    
    Write (messages(1), '(a)')    ' '
    Write (messages(2), '(1x,a)') 'Adjust relevant settings of the DL_MG solver.'
    Call info(messages, 2)
    Call set_reference_database(ref_extra_data, num_ref_data, 'ONETEP', 'DL_MG')
    extradir_header=.False.
    Call check_settings_set_extra_directives(ref_extra_data, num_ref_data, simulation_data%extra_directives,&
                                               & exceptions, 0, extradir_header)
    
  End Subroutine advise_multigrid_onetep
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Electrolyte  
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine define_electrolyte_onetep(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check electrolyte related directives for the generation of input
    ! files for atomistic level simulations with implicit solvent+ electrolyte
    !
    ! author    - i.scivetti September 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256)  :: messages(9)
    Character(Len=256)  :: message
    Character(Len=256)  :: error_elect
    Logical             :: error 
    Integer(Kind=wi)    :: i, j           

    error=.False.
    error_elect = '***ERROR in &electrolyte (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'


    ! method for the solution of the Poisson-Boltzmann equation (compulsory)
    Write (messages(2),'(1x,a)')  'Implemented options for the treatment of the Poisson-Boltzmann equation:'
    Write (messages(3),'(1x,a)')  '- Linearised     (linear-simplified approximation)'
    Write (messages(4),'(1x,a)')  '- Full           (full, non-linear equation)'
    If (simulation_data%electrolyte%solver%fread) Then
      If (simulation_data%electrolyte%solver%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_elect), 'Wrong settings for "solver" directive.'
        error=.True.
      Else
        If (Trim(simulation_data%electrolyte%solver%type) /= 'linearised' .And. &
            Trim(simulation_data%electrolyte%solver%type) /= 'full') Then
          Write (messages(1),'(2(1x,a))') Trim(error_elect), 'Specification of "solver" is not valid'
          error=.True.
        End If
      End If
    Else
        Write (messages(1),'(2(1x,a))') Trim(error_elect), 'The user must specify directive "solver"'
        error=.True.
    End If

    If (error) Then
      Call info(messages,4)
      Call error_stop(' ')
    End If
    
    ! Temperature of the Boltzmann ions
    If (simulation_data%electrolyte%boltzmann_temp%fread) Then
      If (simulation_data%electrolyte%boltzmann_temp%fail) Then
        Write (message,'(2(1x,a))') Trim(error_elect), 'Wrong (or missing) settings for "boltzmann_temperature" directive.&
                                & Both value and units are needed.'
        Call error_stop(message)
      Else
        If (simulation_data%electrolyte%boltzmann_temp%value <= 0.0_wp) Then
          Write (message,'(2(1x,a))') Trim(error_elect), &
                                    &'Input value for "boltzmann_temperature" MUST be larger than zero!!'
          Call error_stop(message)
        End If
        If (Trim(simulation_data%electrolyte%boltzmann_temp%units) /= 'k') Then
          Write (message,'(2(1x,a))') Trim(error_elect), &
                                    &'Wrong units for directive "boltzmann_temperature". Units must be in K'
          Call error_stop(message)
        End If
      End If
    Else
        Write (message,'(2(1x,a))') Trim(error_elect), 'The user must specify the "boltzmann_temperature" directive.'
        Call error_stop(message)
    End If

    ! Method to neutralise the charge (compulsory)
    Write (messages(2),'(1x,2a)') 'Implemented options for charge neutralization. See Ref. ', Trim(bib_neutral_onetep)
    Write (messages(3),'(1x,a)')  '- Jellium                 (Jellium method)'
    Write (messages(4),'(1x,a)')  '- Accessible_Jellium      (Modified jellium method)'
    Write (messages(5),'(1x,a)')  '- Counterions_Auto        (Use concentration with optimal shift parameters)'
    Write (messages(6),'(1x,a)')  '- Counterions_Auto_Linear (Linear version of Counterions_Auto)'
    Write (messages(7),'(1x,a)')  '- Counterions_Fixed       (Use concentration with shift parameters from %boltzmann_ions)'
    If (simulation_data%electrolyte%neutral_scheme%fread) Then
      If (simulation_data%electrolyte%neutral_scheme%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_elect), 'Wrong settings for "neutral_scheme" directive.'
        error=.True.
      Else
        If (Trim(simulation_data%electrolyte%neutral_scheme%type) /= 'jellium' .And. &
            Trim(simulation_data%electrolyte%neutral_scheme%type) /= 'accessible_jellium'      .And. &
            Trim(simulation_data%electrolyte%neutral_scheme%type) /= 'counterions_auto'        .And. &
            Trim(simulation_data%electrolyte%neutral_scheme%type) /= 'counterions_auto_linear' .And. &
            Trim(simulation_data%electrolyte%neutral_scheme%type) /= 'counterions_fixed') Then
          Write (messages(1),'(2(1x,a))') Trim(error_elect), 'Specification of "neutral_scheme" is not valid'
          error=.True.
        End If
      End If
    Else
        Write (messages(1),'(2(1x,a))') Trim(error_elect), 'The user must specify directive "neutral_scheme"'
        error=.True.
    End If

    If (error) Then
      Call info(' ', 1) 
      Call info(messages, 7)
      Call error_stop(' ')
    End If

    
    ! Prevent no neutralization scheme if GC-DFT is requested
    If (simulation_data%dft%gc%activate%stat .And. &
        Trim(simulation_data%electrolyte%neutral_scheme%type) == 'none') Then
      Write (message,'(2(1x,a))') Trim(error_elect), 'Option "none" for "neutral_scheme" is no compatible with a GC-DFT&
                                & simulation, where the number of electrons will be changed during the simulation.&
                                & Please select a neutralization scheme.'
        Call error_stop(message)
    End If

    ! Definition of &boltzmann_ions
    If (.Not. simulation_data%electrolyte%boltzmann_ions_info%stat) Then
       Write (message,'(2(1x,a))') Trim(error_elect), 'The user must define &boltzmann_ions'
       Call error_stop(message)
    End If

    ! Check the necs_shift
    If (Trim(simulation_data%electrolyte%neutral_scheme%type) == 'counterions_fixed' .And. &
       (.Not. simulation_data%electrolyte%set_necs_shift)) Then
         Write (message,'(2(1x,a))') Trim(error_elect), 'The option "counterions_fixed" for "neutral_scheme" needs the&
                                & definition of field "necs_shift" in &boltzmann_ions'
        Call error_stop(message)
    End If

    If (Trim(simulation_data%electrolyte%neutral_scheme%type) /= 'counterions_fixed' .And. &
       (simulation_data%electrolyte%set_necs_shift)) Then
         Write (message,'(2(1x,a))') Trim(error_elect), 'The selected option for "neutral_scheme" does NOT need the&
                                & definition of field "necs_shift" in &boltzmann_ions. Please remove the "necs_shift"&
                                & field.'
        Call error_stop(message)
    End If
    
    ! Definition of &boltzmann_ions
    If (.Not. simulation_data%electrolyte%solvent_radii_info%stat) Then
      Write (message,'(2(1x,a))') Trim(error_elect), 'The user must define &solvent_radii (for each atomic species)'
      Call error_stop(message)
    End If
    
    ! check solvent_radii
    ! Check if user has included all the tags
    If (simulation_data%electrolyte%solvent_radii_info%stat) Then
      Do i=1, simulation_data%total_tags
          error=.True.
          Do j=1, simulation_data%total_tags
            If (Trim(simulation_data%electrolyte%solvent_radii(i)%tag)==Trim(simulation_data%component(j)%tag)) Then
              simulation_data%electrolyte%solvent_radii(i)%element=simulation_data%component(j)%element
              error=.False.
            End If
          End Do
          If (error) Then
            Write (message,'(1x,a,1x,3a)') Trim(error_elect),&
                                     & 'Atomic tag "', Trim(simulation_data%electrolyte%solvent_radii(i)%tag), '" declared in&
                                     & "&solvent_radii" has not been defined in "&block_input_composition". Please check'
            Call error_stop(message)
          End If
      End Do
    End If

    ! Method to neutralise the charge (compulsory)
    Write (messages(2),'(1x,2a)') 'Implemented options for steric potential. See Ref. ', Trim(bib_gcdft_onetep)
    Write (messages(3),'(1x,a)')  '- Hard-core               (Not recommended due to numerical problems)'
    Write (messages(4),'(1x,a)')  '- Smoothed-Hard-core      (RECOMMENDED)'
    If (simulation_data%electrolyte%steric_potential%fread) Then
      If (simulation_data%electrolyte%steric_potential%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_elect), 'Wrong settings for "steric_potential" directive.'
        error=.True.
      Else
        If (Trim(simulation_data%electrolyte%steric_potential%type) /= 'hard-core' .And. &
            Trim(simulation_data%electrolyte%steric_potential%type) /= 'smoothed-hard-core') Then
          Write (messages(1),'(2(1x,a))') Trim(error_elect), 'Specification of "steric_potential" is not valid'
          error=.True.
        End If
      End If
    Else
        Write (messages(1),'(2(1x,a))') Trim(error_elect), 'The user must specify directive "steric_potential"'
        error=.True.
    End If

    If (error) Then
      Call info(' ', 1) 
      Call info(messages, 4)
      Call error_stop(' ')
    End If
    
    ! Steric isodensity
    If (simulation_data%electrolyte%steric_isodensity%fread) Then
      If (simulation_data%electrolyte%steric_isodensity%fail) Then
        Write (message,'(2(1x,a))') Trim(error_elect), 'Wrong (or missing) settings for "steric_isodensity" directive.'
        Call error_stop(message)
      Else
        If (simulation_data%electrolyte%steric_isodensity%value <= 0.0_wp) Then
          Write (message,'(2(1x,a))') Trim(error_elect), &
                                    &'Input value for "steric_isodensity" MUST be larger than zero!!'
          Call error_stop(message)
        End If
      End If
    Else
      simulation_data%electrolyte%steric_isodensity%value=0.003_wp 
    End If

    
    ! Steric smearing
    If (simulation_data%electrolyte%steric_smearing%fread) Then
      If (simulation_data%electrolyte%steric_smearing%fail) Then
        Write (message,'(2(1x,a))') Trim(error_elect), 'Wrong (or missing) settings for "steric_smearing" directive.&
                                & Both value and units are needed.'
        Call error_stop(message)
      Else
        If (simulation_data%electrolyte%steric_smearing%value <= 0.0_wp) Then
          Write (message,'(2(1x,a))') Trim(error_elect), &
                                    &'Input value for "steric_smearing" MUST be larger than zero!!'
          Call error_stop(message)
        End If
        If (Trim(simulation_data%electrolyte%steric_smearing%units) /= 'k') Then
          Write (message,'(2(1x,a))') Trim(error_elect), &
                                    &'Wrong units for directive "steric_smearing". Units must be in Bohr'
          Call error_stop(message)
        End If
      End If
    Else
      If (Trim(simulation_data%electrolyte%steric_potential%type) == 'smoothed-hard-core') Then
         simulation_data%electrolyte%steric_smearing%value=0.4_wp 
      End If
    End If

    If (simulation_data%electrolyte%steric_smearing%fread .And. &
        Trim(simulation_data%electrolyte%steric_potential%type) /= 'smoothed-hard-core') Then
        Write (message,'(2(1x,a))') Trim(error_elect), &
                                    & 'Definition of "steric_smearing" is only relevant if the smoothed-hard-core option&
                                    & is selected for the "steric_potential" directive.  Please review the settings.'
          Call error_stop(message)
    End If    
    
    ! Steric isodensity
    If (simulation_data%electrolyte%capping%fread) Then
      If (simulation_data%electrolyte%capping%fail) Then
        Write (message,'(2(1x,a))') Trim(error_elect), 'Wrong (or missing) settings for "capping" directive.'
        Call error_stop(message)
      Else
        If (simulation_data%electrolyte%capping%value < 0.0_wp) Then
          Write (message,'(2(1x,a))') Trim(error_elect), &
                                    &'Input value for "capping" MUST be larger than zero!!'
          Call error_stop(message)
        End If
      End If
    Else
      If (Trim(simulation_data%electrolyte%solver%type) == 'full') Then 
        simulation_data%electrolyte%capping%value=0.0_wp 
      End if
    End If

    
    If (simulation_data%electrolyte%capping%fread .And. &
        Trim(simulation_data%electrolyte%solver%type) /= 'full') Then
        Write (message,'(2(1x,a))') Trim(error_elect), &
                                    & 'Definition of "capping" is irrelevant if "solver" is not set to "full".&
                                    & Please review the settings.'
        Call error_stop(message)
    End If    

    
  End Subroutine define_electrolyte_onetep
  
  Subroutine print_onetep_electrolyte(iunit, ic, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print electrolyte directives for ONETEP
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Integer(Kind=wi), Intent(InOut) :: ic    
    Type(simul_type), Intent(InOut) :: simulation_data

    Integer(Kind=wi)   :: i
    Character(Len=256) :: message
        
     Write (iunit,'(2a)') '#==== Electrolyte settings (Poisson-Boltzmann Equation). See ', Trim(bib_pbeq_onetep)
     Write (message,'(2a)') 'is_pbe : ', Trim(simulation_data%electrolyte%solver%type) 
     Call record_directive(iunit, message, 'is_pbe', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,f8.2,a)') 'is_pbe_temperature :       ', simulation_data%electrolyte%boltzmann_temp%value, '   K' 
     Call record_directive(iunit, message, 'is_pbe_temperature', simulation_data%set_directives%array(ic), ic)
     Write (message,'(2a)') 'is_pbe_neutralisation_scheme :   ', Trim(simulation_data%electrolyte%neutral_scheme%type) 
     Call record_directive(iunit, message, 'is_pbe_neutralisation_scheme', simulation_data%set_directives%array(ic), ic)

     Write (iunit,'(a)') '# Steric potential'
     If (Trim(simulation_data%electrolyte%steric_potential%type) == 'smoothed-hard-core') Then
       Write (message,'(a)') 'is_steric_pot_type :   M   # Smoothed-hard-core'
     Else If (Trim(simulation_data%electrolyte%steric_potential%type) == 'hard-core') Then
       Write (message,'(a)') 'is_steric_pot_type :   H   # Hard-core (Not recommended)'
     End If
     Call record_directive(iunit, message, 'is_steric_pot_type', simulation_data%set_directives%array(ic), ic)

     Write (iunit,'(a)') '# Size of solvation shell for explicit atomic species'         
     Write (iunit,'(a)') '%block species_solvent_radius'
     Do i=1, simulation_data%total_tags
       Write (iunit,'(3x,a,3x,f6.2)') Trim(simulation_data%electrolyte%solvent_radii(i)%tag),&
                                  simulation_data%electrolyte%solvent_radii(i)%value  
     End Do
     Write (iunit,'(a)') '%endblock species_solvent_radius'

     Write (message,'(a,f9.5)') 'is_hc_steric_dens_isovalue :   ', simulation_data%electrolyte%steric_isodensity%value
     Call record_directive(iunit, message, 'is_hc_steric_dens_isovalue', simulation_data%set_directives%array(ic), ic)
     If (Trim(simulation_data%electrolyte%steric_potential%type) == 'smoothed-hard-core') Then
       Write (message,'(a,f7.3,a)') 'is_hc_steric_smearing :     ', simulation_data%electrolyte%steric_smearing%value, '  Bohr' 
       Call record_directive(iunit, message, 'is_hc_steric_smearing', simulation_data%set_directives%array(ic), ic)
     End If

     Write (iunit,'(a)') '# Information for Boltzmann ions'
     Write (iunit,'(a)') '%block sol_ions'
     If (Trim(simulation_data%electrolyte%neutral_scheme%type) == 'counterions_fixed')Then
       Write (iunit,'(3x,a)') '#Species   Charge          Conc.          NECS_shift'
       Do i=1, simulation_data%electrolyte%number_boltzmann_ions
         Write (iunit,'(3x,a,3(8x,f7.3))') Trim(simulation_data%electrolyte%boltzmann_ions(i)%tag),&
                                           simulation_data%electrolyte%boltzmann_ions(i)%charge,   &
                                           simulation_data%electrolyte%boltzmann_ions(i)%conc,     &
                                           simulation_data%electrolyte%boltzmann_ions(i)%necs_shift
       End Do
     Else
       Write (iunit,'(3x,a)') '#Species   Charge          Conc.'
       Do i=1, simulation_data%electrolyte%number_boltzmann_ions
         Write (iunit,'(3x,a,2(8x,f7.3))') Trim(simulation_data%electrolyte%boltzmann_ions(i)%tag),&
                                         simulation_data%electrolyte%boltzmann_ions(i)%charge,   &
                                         simulation_data%electrolyte%boltzmann_ions(i)%conc
       End Do                                  
     End If
     Write (iunit,'(a)') '%endblock sol_ions'
  
  End Subroutine print_onetep_electrolyte
  
  
  Subroutine summary_electrolyte_onetep(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to summarise electrolyte settings from the information
    ! provided by the user via block &electrolyte
    !
    ! author    - i.scivetti September 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(7)

    Write (messages(1),'(2(1x,a))')  ' - level of approximation for the Poisson-Boltzmann equation: ', &
                                     & Trim(simulation_data%electrolyte%solver%type)
    Write (messages(2),'(1x,a,f6.2,1x,a)') ' - temperature of the Boltzmann ions: ', &
                                           & simulation_data%electrolyte%boltzmann_temp%value, ' K'
    Write (messages(3),'(2(1x,a))')  ' - scheme to neutralise the charge: ', Trim(simulation_data%electrolyte%neutral_scheme%type)
    Write (messages(4),'(1x,a)')     ' - solvation shell for explicit atomic species is defined in the&
                                    & "species_solvent_radius" block'
    Write (messages(5),'(2(1x,a))')     ' - Steric potential type: ', Trim(simulation_data%electrolyte%steric_potential%type)
    Write (messages(6),'(1x,a)')       ' - information for the Boltzmann ions is defined in the "sol_ions" block'
    Call info(messages, 6) 

  End Subroutine summary_electrolyte_onetep

  
  Subroutine advise_electrolyte_onetep(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about ONETEP settings
    ! for electrolyte
    !
    ! author    - i.scivetti September 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(9)

    Write (messages(1), '(1x,a)') '====== Hints for the electrolyte'
    Call info(messages, 1)
    
    If (Trim(simulation_data%electrolyte%steric_potential%type) /= 'smoothed-hard-core') Then
      Write (messages(1), '(1x,a)') '- the recommended option for the steric_potential is "smoothed-hard-core".&
                                    & The user must consider changing this setting'
      Call info(messages, 1)                             
    End If
    
    Write (messages(1), '(1x,a)') '- In case of troubleshooting, please consult the ONETEP manual for possible causes&
                                  & and solutions.'
    Call info(messages, 1)                             

  End Subroutine advise_electrolyte_onetep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Pseudo-potentials  
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine check_pseudo_potentials_onetep(simulation_data, error_dft) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check PPs for setting input files  
    !
    ! author    - i.scivetti September 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(InOut) :: simulation_data
    Character(Len=256), Intent(In   ) :: error_dft
    
    Character(Len=256) :: pp_path, pp_file, pp_name
    Character(Len=256) :: ref_pp_extension, pp_extension
    
    Character(Len=256) :: message

    Integer(Kind=wi)   :: i, internal
    
    pp_path   = Trim(FOLDER_DFT)//'/PPs/'    
    
    ! Check if all PP files correspond to the same format
      Do i =1, simulation_data%total_tags
        ! check the extension of the file (recpot vs. abinit vs. unknown)
        Write (pp_name, '(a)')  Trim(simulation_data%dft%pseudo_pot(i)%file_name)
        If (Index(pp_name,'.abinit') /= 0 .Or. Index(pp_name,'.recpot') /= 0) Then

          If (Index(pp_name,'.recpot') /= 0) Then
            pp_extension='.recpot'
            If (i==1) Then
              ref_pp_extension='.recpot'
            End If         

          Else If (Index(pp_name,'.abinit') /= 0) Then
            pp_extension='.abinit'
            !---------------------------------------------
            ! This should be removed when ONETEP+PAW works
            !---------------------------------------------
            Write (message,'(1x,4a)') Trim(error_dft), ' Up to ', Trim(date_RELEASE), ', the PAW method in ONETEP is under&
                                     & revision. The user is advised to use norm-conserving pseudopotentials (files with&
                                     & extension .recpot).'
            Call error_stop(message)
            !--------------------------------------------
            !--------------------------------------------
            If (i==1) Then
              simulation_data%dft%onetep_paw = .True.
              ref_pp_extension='.abinit'
            End If

          End If

          If (Trim(pp_extension) /= Trim(ref_pp_extension)) Then
             Write (message, '(1x,3a)') '***ERROR: pseudo potential files defined in &pseudo_potentials contain&
                                       & different extensions. For ONETEP, all files must be either ".recpot"&
                                       & or ".abinit" type.'
             Call error_stop(message)
          End If

        Else
          Write (message,'(1x,5a)') '***ERROR: The extension of file ',  Trim(pp_name), ' for the pseudo potential of species "',&
                                  & Trim(simulation_data%dft%pseudo_pot(i)%tag), '" is not recognised by ONETEP.&
                                  & Valid extensions are ".abinit" and ".recpot".'
          Call error_stop(message)
        End If 

      End Do 

      ! Check consistency between pseudpotentials and XC directive
      Do i=1, simulation_data%total_tags
        pp_file=Trim(pp_path)//Trim(simulation_data%dft%pseudo_pot(i)%file_name)
        ! Open PP file
        Open(Newunit=internal, File=Trim(pp_file), Status='old')
        If (Trim(ref_pp_extension) =='.recpot') Then
          Call check_recpot_onetep(internal, simulation_data, i)  
        Else If (Trim(ref_pp_extension) == '.abinit') Then
          Call check_abinit(internal, simulation_data, i) 
        End If
        ! Close PP file
        Close(internal)
      End Do

  End Subroutine check_pseudo_potentials_onetep    
      
  Subroutine check_recpot_onetep(internal, simulation_data, i)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if PP file is consistent or not with the .recpot format for ONETEP 
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
      If (Trim(xc) /= 'lda' .And. Trim(xc) /= 'gga') Then
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

  End Subroutine check_recpot_onetep

  Subroutine check_abinit(internal, simulation_data, i)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if PP file is consistent or not with the .abinit format for onetep
    !
    ! author    - i.scivetti September 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   )   :: internal 
    Type(simul_type), Intent(In   )   :: simulation_data
    Integer(Kind=wi), Intent(In   )   :: i 
    
    Character(Len=256) :: message
    Character(Len=256) :: pp_name, element, end_file
    Integer(Kind=wi)   :: io, j
 
    Write (pp_name, '(a)')  Trim(simulation_data%dft%pseudo_pot(i)%file_name)
    Write (end_file,'(1x,2a)') '***ERROR: incomplete (or empty) pseudo potential file ', Trim(pp_name)

    Read (internal, Fmt=*, iostat=io) (element, j=1, 6)
    If (is_iostat_end(io)) Then
      Call error_stop(end_file)
    End If
    If (Index(pp_name,'PBE') /= 0 .Or. Index(pp_name,'pbe') /= 0) Then
      If (Trim(simulation_data%dft%xc_base) /= 'pbe' ) Then
        Write (message, '(1x,3a)') '***ERROR: Inconsistency between the version of pseudo potential file "', Trim(pp_name),& 
                                & '" (PBE) and the specification of directive "XC_version".'
        Call error_stop(message) 
      End if 
    Else If (Index(pp_name,'LDA') /= 0 .Or. Index(pp_name,'lda') /= 0) Then
      If (Trim(simulation_data%dft%xc_base) /= 'lda') Then
        Write (message, '(1x,3a)') '***ERROR: Inconsistency between the level of pseudo potential file "', Trim(pp_name),& 
                                & '" (LDA) and the specification of directive "XC_level".'
        Call error_stop(message) 
      End if
 
    Else 
      Write (message, '(1x,2a)') '***ERROR: Unrecognizable XC level/version for pseudo potential file ', Trim(pp_name)
      Call error_stop(message)
    End If 

    ! More checks
    If ((Trim(element) /= Trim(simulation_data%dft%pseudo_pot(i)%element))) Then
      Write (message, '(1x,5a)') '***ERROR: pseudo potential file ', Trim(pp_name), ' does not correspond to element "', &
                                & Trim(simulation_data%dft%pseudo_pot(i)%element), '" as specified by the tag&
                                & in sub-block &pseudo_potentials'
      Call error_stop(message) 
    End if 

  End Subroutine check_abinit 

End Module code_onetep
