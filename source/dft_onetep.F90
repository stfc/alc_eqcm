!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to check, define and print DFT directives for simulations
! with ONETEP. This module also warns the user about aspects to take
! into consideration when performing simulations
!
! Copyright - Ada Lovelace Centre
!
! Author    - i.scivetti Oct 2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module dft_onetep
  
  Use constants,        Only : date_RELEASE, &
                               max_components

  Use fileset,          Only : file_type, &
                               FILE_SET_EQCM, &
                               FILE_SET_SIMULATION,&
                               FOLDER_DFT 
  Use numprec,          Only : wi, &
                               wp
  Use process_data,     Only : capital_to_lower_case
  Use references,       Only : bib_AVV10s, bib_blyp, bib_dftd2, bib_optb88, bib_optpbe, bib_pbe, bib_pbesol, &
                               bib_pw91, bib_pw92, bib_pz, bib_revpbe, bib_rp, bib_vdwdf, bib_vdwdf2, bib_vv10, &
                               bib_vwn, bib_xlyp
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

  Public :: define_onetep_settings, print_onetep_settings, warnings_onetep
  
Contains

  Subroutine define_onetep_settings(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define (and check) settings for ONETEP directives (set defaults values)
    !
    ! author    - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256) :: message, messages(15)
    Character(Len=256) :: error_dft, error_motion
    Character(Len=256) :: path
    Character(Len=256) :: pp_path, pp_file, pp_name
    Character(Len=256) :: ref_pp_extension, pp_extension

    Integer(Kind=wi)   :: i, ifile, internal 
    Logical            :: error
  
    simulation_data%dft%onetep_paw = .False. 
    error_dft    = '***ERROR in &dft_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
    error_motion = '***ERROR in &motion_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'

    pp_path   = Trim(FOLDER_DFT)//'/PPs/'    

    ! latest vrrsion of the code
    simulation_data%code_version= '6.0'

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
    !!!!!!!!!!!!!!!!!!
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
          Write (message,'(1x,4a)') Trim(error_dft), ' Up to ', Trim(date_RELEASE), ', the PAW method in ONETEP is under revision.&
                                   & The user is advised to use norm-conserving pseudopotentials (files with extension .recpot).'
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

      If (simulation_data%dft%spin_polarised%fread .And. Trim(simulation_data%dft%vdw%type) /= 'dft-d2') Then
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

    ! Ions related settings
    !!!!!!!!!!!!!!!!!!!!!!!
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
        Call check_extra_directives(simulation_data%extra_directives%array(i), ':', 'ONETEP', i)
      End Do
    End If

    ! Prevent definition of masses for ONETEP
    If (simulation_data%motion%mass_info%stat) Then
      Write (message,'(2(1x,a))') Trim(error_motion), 'Explicit definition of atomic masses is not implemented ONETEP.&
                                 & Please remove the &masses block.'
      Call error_stop(message)
    End If
      
  End Subroutine define_onetep_settings  

  Subroutine print_onetep_settings(files, net_elements, list_tag, list_number, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print settings from ONETEP directives
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),  Intent(InOut) :: files(:)
    Integer(Kind=wi), Intent(In   ) :: net_elements
    Character(Len=8), Intent(In   ) :: list_tag(max_components) 
    Integer(Kind=wi), Intent(In   ) :: list_number(max_components)
    Type(simul_type), Intent(InOut) :: simulation_data

    Integer(Kind=wi)   :: i, j, k
    Integer(Kind=wi)   :: iunit, atomic_number
    Logical            :: loop, loop2
    Real(Kind=wp)      :: mag_ini(max_components)
    Character(Len=256) :: thermo, word, word2
    Character(Len=256) :: message, messages(2), tag

    Integer(Kind=wi)   :: ic
    Logical            :: found

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
      Write (iunit,'(3a)') '#==== Ensemble DFT'
      Write (message,'(a)')           'edft  :  T'
      Call record_directive(iunit, message, 'edft', simulation_data%set_directives%array(ic), ic)
      Write (message,'(a,f8.2,2x,a)') 'edft_smearing_width : ', simulation_data%dft%width_smear%value, ' eV  # smearing width'  
      Call record_directive(iunit, message, 'edft_smearing_width', simulation_data%set_directives%array(ic), ic)

      If (simulation_data%dft%bands%fread) Then
          Write (message,'(a,i4,a)') 'edft_extra_bands : ', simulation_data%dft%bands%value,&
                                & ' # add extra bands to reach/improve convergence'
      Else 
          Write (message,'(a)')       'edft_extra_bands :  -1 # bands is equal to the total number of NGWFs' 
      End If
      Call record_directive(iunit, message, 'edft_extra_bands', simulation_data%set_directives%array(ic), ic)
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
             Do ic = 1, 3
               If (2*simulation_data%dft%ngwf(j)%radius >= simulation_data%cell_length(ic)) Then
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

   ! Total number of set directives
   simulation_data%set_directives%N0=ic-1

   If (simulation_data%extra_info%stat) Then
     Write (iunit,'(a)') ' '
     Write (iunit,'(a)') '##### User defined directives'
     Write (iunit,'(a)') '#============================'
     found=.False.
     Do i=1, simulation_data%extra_directives%N0
       Write (iunit,'(a)') Trim(Adjustl(simulation_data%extra_directives%array(i)))
       If (Index(Trim(Adjustl(simulation_data%extra_directives%array(i))), '#') /= 1 ) Then
         Read(simulation_data%extra_directives%array(i), Fmt=*) word, word2
         Call capital_to_lower_case(word)
         Call capital_to_lower_case(word)
           If (Trim(word)=='&block' .And. Trim(word2)=='thermostat') Then
             If (Trim(simulation_data%simulation%type) == 'md') Then
               Write (message, '(1x,a)') '***ERROR in sub-block &extra_directives: "&block thermostat" CANNOT be defined as&
                                        & part of &extra_directives. Please remove it.'
               Call info(message, 1)                         
               Call error_stop(' ')                         
             End If  
           Else
           Call scan_extra_directive(simulation_data%extra_directives%array(i), simulation_data%set_directives, &
                                 & ':', found, tag)         
           If (found)Then
             Close(iunit)
             Call info(' ', 1)
             Write (messages(1), '(1x,3a)') '***ERROR in sub-bock &extra_directives: directive "', Trim(tag),&
                                     & '" has already been defined.'
             Write (messages(2), '(1x,3a)') 'Please check temporary file ', Trim(files(FILE_SET_SIMULATION)%filename), &
                                     & '. The user must review the settings of &extra_directives.&
                                     & Directives must not be duplicated.'
             Call info(messages, 2)
             Call error_stop(' ')
           Else
             simulation_data%set_directives%N0=simulation_data%set_directives%N0+1
             simulation_data%set_directives%array(simulation_data%set_directives%N0)=Trim(tag)
           End If
         End If 
       End If
     End Do
   End If

   Write (iunit,'(a)') ' '
   Write (iunit,'(a)') '#### Simulation cell and atomic coordinates'
   Write (iunit,'(a)') '#=========================================='

   Close(iunit)

  End Subroutine print_onetep_settings
  
  Subroutine warnings_onetep(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about ONETEP settings 
    !
    ! author    - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(9), header
    Character(Len=256)  :: in_extra
    Logical             :: warning, error, print_header
    Integer(Kind=wi)    :: i

    warning=.False.
    print_header=.True.
    in_extra='using the &extra_directives block'
   
    If (simulation_data%extra_info%fread) Then
      Call info(' ', 1)
      Call info(' ***WARNING****************', 1)
      Write (messages(1), '(1x,a)') ' - ALC_EQCM only checked that the information added in "&extra_directives" has not been'
      Write (messages(2), '(1x,a)') '   already defined from the settings provided in "&block_simulation_settings"'   
      Write (messages(3), '(1x,a)') ' - specification of directives in "&extra_directives" for functionalities that are not'
      Write (messages(4), '(1x,a)') '   implemented in ALC_EQCM have been printed, but its correctness is full responsibility&
                                    & of the user'
      Call info(messages, 4)
      Call info(' **************************', 1)
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
                                   & the value of given to "export OMP_NUM_THREADS" en el job submission script'
    Write (messages(9), '(1x, a)') ' - "threads_num_mkl" must only be defined is MKL is used for compilation'
    Call info(messages, 9)
    Call info(' ', 1)
    Write (messages(1), '(1x,a)')  'In case of problems in the electronic convergence, the user should:'
    Write (messages(2), '(1x,a)')  ' - check the number and radii of the NGWF for each species via block &ngwf' 
    Write (messages(3), '(1x,2a)') ' - adjust the values for "minit_lnv" and "maxit_lnv" ', Trim(in_extra)
    Write (messages(4), '(1x,a)')  ' - increase the value of "maxit_ngwf_cg" via ALC_EQCM directive "scf_steps"'
    Write (messages(5), '(1x,2a)') ' - optimise the initialisation of the density kernel via "maxit_pulser_nano"&
                                    &and/or "maxit_pen" ', Trim(in_extra)
    If (simulation_data%dft%edft%fread) Then
      Write (messages(6), '(1x,a)')  ' - increase the value of "edft_extra_bands" using the ALC_EQCM directive "bands"'
      Write (messages(7), '(1x,a)')  ' - change the value of "edft_smearing_width" via the ALC_EQCM directive "width_smearing"'
      Write (messages(8), '(1x,2a)') ' - increase the value of "edft_init_maxit" ', Trim(in_extra)
      Call info(messages, 8)
    Else
      Write (messages(6), '(1x,a)')  ' - check if the system tends to metallise, in which case the user must ser EDFT'
      Call info(messages, 6)
    End If   

    If (.Not. simulation_data%dft%edft%fread) Then
      Write (messages(1), '(1x,2a)') 'For linear scaling DFT the user should optimise "kernel_cutoff" ', Trim(in_extra)
      Call info(messages, 1)
    End If   

    Write (messages(1), '(1x,3a)') 'I/O can be controlled ', Trim(in_extra), ' with the following directives (see ONETEP manual):'
    Write (messages(2), '(1x,a)')  ' - write_denskern, write_tightbox_ngwfs, write_converged_dk_ngwfs'
    Write (messages(3), '(1x,a)')  ' - read_denskern, read_tightbox_ngwfs, read_hamiltonian (EDFT)'
    Write (messages(4), '(1x,a)')  ' - output_detail (to specify the level of detail for the generated output)'
    Call info(messages, 4)
   
    Write (messages(1), '(1x,3a)') 'Atomic forces can be printed by setting "write_forces : T" ', Trim(in_extra), '.'
    Call info(messages, 1)
    If (Trim(simulation_data%simulation%type) == 'md' .Or. Trim(simulation_data%simulation%type) == 'relax_geometry') Then 
      Write (messages(1), '(1x,3a)') 'By setting "write_xyz : T" (', Trim(in_extra), ') the atomic coordinates are printed&
                                    & as a .xyz file.'
      Call info(messages, 1)
    End If

    If (simulation_data%dft%hubbard_info%fread  .Or. simulation_data%net_charge%fread  .Or. simulation_data%dft%vdw%fread) Then
       warning=.True.
    End If
    
    If (Trim(simulation_data%simulation%type) == 'md') Then
      If (Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
        If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover' .Or. &
            Trim(simulation_data%motion%thermostat%type) == 'andersen'    .Or. &
            Trim(simulation_data%motion%thermostat%type) == 'langevin' ) Then
            warning=.True.
        End If
      End If
    End If

    If (warning) Then
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
      
          ! MD-related parameters
      If (Trim(simulation_data%simulation%type) == 'md') Then
        If (Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
            Write (messages(1), '(1x,a)')  ' - changing "tau" using ALC_EQCM directive "relax_time_thermostat"'
            Call info(messages, 1)  
            Write (messages(1), '(1x,a)')  'Unfortunately, due to the structure of "%block thermostat", changes to the following&
                                       & directives, if required, must be applied manually (&extra_directives cannot be used)'
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
      
    End If

  End Subroutine warnings_onetep

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

End Module dft_onetep
