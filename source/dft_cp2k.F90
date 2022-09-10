!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to check, define and print DFT directives for simulations
! with CP2K. This module also warns the user about aspects to take
! into consideration when performing simulations
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author    - i.scivetti Oct 2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module dft_cp2k
  
  Use constants,        Only : date_RELEASE , &
                               max_components, &
                               Ha_to_eV,       &
                               K_to_eV

  Use fileset,          Only : file_type, &
                               FILE_SET_EQCM, &
                               FILE_SET_SIMULATION,&
                               FOLDER_DFT 
  Use numprec,          Only : wi, &
                               wp
  Use process_data,     Only : get_word_length
  Use references,       Only : bib_xlyp, bib_wigner, bib_wc, bib_vwn, bib_vdwdf2b86r, bib_vdwdf2, bib_vdwdf, bib_scanrvv10, &
                               bib_rp, bib_revpbe, bib_pz, bib_pw91, bib_pw92, bib_pbesol, bib_pbe, bib_hl, bib_dftd2,      &
                               bib_dftd3, bib_dftd3bj, bib_optb86b, bib_optb88, bib_am05, bib_blyp, bib_optpbe, bib_pade,   &
                               bib_slater, bib_tunega, bib_rpw86, bib_scan
  Use simulation_tools, Only : simul_type,&
                               check_initial_magnetization, &
                               print_warnings 
  Use stoichiometry,    Only : stoich_type
  Use unit_output,      Only : error_stop,&
                               info 

  Implicit None
  Private

  Public :: define_cp2k_settings, print_cp2k_settings, warnings_cp2k
  
Contains

  Subroutine define_cp2k_settings(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define (and check) settings for CP2K directives (set defaults values)
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256) :: message, messages(19)
    Character(Len=256) :: error_dft, error_motion
    Character(Len=256) :: exec_grep, path, pp_path
    Character(Len=256) :: word, root, basis, line, potential

    Integer(Kind=wi)   :: i, j, k, l, io, maxcol 
    Integer(Kind=wi)   :: root_length
    Integer(Kind=wi)   :: internal, nlines 
    Logical            :: error, safe, loop
   
    maxcol=1000 
    error_dft    = '***ERROR in &dft_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
    error_motion = '***ERROR in &motion_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'

    pp_path   = Trim(FOLDER_DFT)//'/PPs/'    

    ! latest vrrsion of the code
    simulation_data%code_version= '8.1' 

    ! Check XC_version
    If (Trim(simulation_data%dft%xc_version%type)  /= 'hl'     .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'pz'     .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'wigner' .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'vwn'    .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'pade'   .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'pw92'   .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'pw91'   .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'wc'     .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'am05'   .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'pbe'    .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'rp'     .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'revpbe' .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'pbesol' .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'blyp'   .And.&
      Trim(simulation_data%dft%xc_version%type)  /= 'xlyp' ) Then
      Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                  &'Invalid specification for directive "XC_version" for CP2K.&
                                  & Implemented options for CP2K are:'
      Write (messages(2),'(1x,a)')  '==== LDA-level =================='
      Write (messages(3),'(1x,2a)')  '- HL      Hedin-Lundqvist        ', Trim(bib_hl)  
      Write (messages(4),'(1x,2a)')  '- PZ      Perdew-Zunger          ', Trim(bib_pz)
      Write (messages(5),'(1x,2a)')  '- Wigner  Wigner                 ', Trim(bib_wigner)
      Write (messages(6),'(1x,2a)')  '- VWN     Vosko-Wilk-Nusair      ', Trim(bib_vwn)
      Write (messages(7),'(1x,2a)')  '- PADE    Pade functional        ', Trim(bib_pade)
      Write (messages(8),'(1x,2a)')  '- PW92    Perdew-Wang 92         ', Trim(bib_pw92)
      Write (messages(9),'(1x,a)')   '==== GGA-level =================='
      Write (messages(10),'(1x,2a)') '- PW91    Perdew-Wang 91         ', Trim(bib_pw91)
      Write (messages(11),'(1x,2a)') '- WC      Wu-Cohen               ', Trim(bib_wc)
      Write (messages(12),'(1x,2a)') '- AM05    Armiento-Mattsson      ', Trim(bib_am05)
      Write (messages(13),'(1x,2a)') '- PBE     Perdew-Burke-Ernzerhof ', Trim(bib_pbe)
      Write (messages(14),'(1x,2a)') '- RP      Hammer-Hansen-Norskov  ', Trim(bib_rp)
      Write (messages(15),'(1x,2a)') '- revPBE  revPBE                 ', Trim(bib_revpbe)
      Write (messages(16),'(1x,2a)') '- PBEsol  PBE for solids         ', Trim(bib_pbesol)
      Write (messages(17),'(1x,2a)') '- BLYP    Becke-Lee-Young-Parr   ', Trim(bib_blyp)
      Write (messages(18),'(1x,2a)') '- XLYP    Xu-Goddard             ', Trim(bib_xlyp)
      Write (messages(19),'(1x,a)') '=================================='
      Call info(messages, 19)
      Call error_stop(' ')
    End If

    ! XC base
    If (Trim(simulation_data%dft%xc_version%type)  == 'hl'   .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'pz'     .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'wigner' .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'pade'   .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'pw92'   .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'vwn'    ) Then
      simulation_data%dft%xc_base='PADE'
    Else If (Trim(simulation_data%dft%xc_version%type) == 'pw91' .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'wc'     .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'am05'   .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'pbe'    .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'rp'     .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'revpbe' .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'xlyp' .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'pbesol') Then  
      simulation_data%dft%xc_base='PBE'
    Else If (Trim(simulation_data%dft%xc_version%type)  == 'blyp' ) Then
      simulation_data%dft%xc_base='BLYP'
    End If

    ! vdW settings 
    simulation_data%dft%need_vdw_kernel=.False.
    If (simulation_data%dft%vdw%fread) Then
      If (Trim(simulation_data%dft%vdw%type) /= 'dft-d2'   .And.&
         Trim(simulation_data%dft%vdw%type) /= 'dft-d3'    .And.&
         Trim(simulation_data%dft%vdw%type) /= 'dft-d3-bj' .And.&
         Trim(simulation_data%dft%vdw%type) /= 'vdw-df'    .And.&
         Trim(simulation_data%dft%vdw%type) /= 'optb88'    .And.&
         Trim(simulation_data%dft%vdw%type) /= 'optb86b'   .And.&
         Trim(simulation_data%dft%vdw%type) /= 'optpbe'    .And.&
         Trim(simulation_data%dft%vdw%type) /= 'vdw-df2'   .And.&
         Trim(simulation_data%dft%vdw%type) /= 'vdw-df2-b86r' .And.&
         Trim(simulation_data%dft%vdw%type) /= 'scan+rvv10'   ) Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                  &'Invalid specification of directive "vdW" for CP2K. Valid options are:'
        Write (messages(2),'(1x,2a)')  '- DFT-D2         Grimme D2 ', Trim(bib_dftd2)
        Write (messages(3),'(1x,2a)')  '- DFT-D3         Grimme D3 with no damping ', Trim(bib_dftd3)
        Write (messages(4),'(1x,2a)')  '- DFT-D3-BJ      Grimme D3 with Becke-Jonson damping ', Trim(bib_dftd3bj)
        Write (messages(5),'(1x,2a)')  '- vdW-DF         X (revPBE), C (LDA), vdW (vdW-DF) ', Trim(bib_vdwdf)  
        Write (messages(6),'(1x,2a)')  '- optPBE         X (OPTPBE), C (LDA), vdW (vdW-DF) ', Trim(bib_optpbe)    
        Write (messages(7),'(1x,2a)')  '- optB88         X (OPTB88), C (LDA), vdW (vdW-DF) ', Trim(bib_optb88)
        Write (messages(8),'(1x,2a)')  '- optB86b        Optimized Becke86b ', Trim(bib_optb86b)
        Write (messages(9),'(1x,2a)')  '- vdW-DF2        X (rPW86), C (LDA), vdW (vdW-DF 2)  ', Trim(bib_vdwdf2)
        Write (messages(10),'(1x,2a)') '- vdW-DF2-B86R   Hamada version of vdW-DF2 ', Trim(bib_vdwdf2b86r)
        Write (messages(11),'(1x,2a)') '- SCAN+rVV10     SCAN + non-local correlation part of the rVV10 ', Trim(bib_scanrvv10)
        Call info(messages, 11)
        Call error_stop(' ')
      End If

      If (Trim(simulation_data%dft%xc_level%type) /= 'gga') Then
        Write (message,'(1x,4a)') Trim(error_dft), &
                                &' Dispersion correction type "', Trim(simulation_data%dft%vdw%type), '" requires the GGA&
                                & option for directive XC_level. Please change'
        Call error_stop(message)
      End If

     ! Check if the kernel exists 
      If (Trim(simulation_data%dft%vdw%type) == 'dft-d3' .Or. &
         Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then
         If (Trim(simulation_data%dft%xc_base)=='PADE') Then
            Write (message,'(1x,4a)') Trim(error_dft), &
                                &' Dispersion correction type "', Trim(simulation_data%dft%vdw%type), '" is not compatible&
                                & with the choice of "', Trim(simulation_data%dft%xc_version%type),&
                                &'" for the XC. Please review'
           Call error_stop(message)
         End If
         simulation_data%dft%need_vdw_kernel=.True.
         simulation_data%dft%vdw_kernel_file='dftd3.dat'
         Inquire(File=Trim(FOLDER_DFT)//'/'//Trim(simulation_data%dft%vdw_kernel_file), Exist=safe)
         If (.not.safe) Then
           Write (message,'(1x,5a)') '***ERROR: File ',&
                                   & Trim(FOLDER_DFT)//'/'//Trim(simulation_data%dft%vdw_kernel_file),&
                                   & ', needed for "', Trim(simulation_data%dft%vdw%type), '" dispersion corrections,&
                                   & does not exist. Please copy the file from the CP2K repository and rerun.'
           Call error_stop(message)
         End If
      End If
    
      If (Trim(simulation_data%dft%vdw%type) == 'vdw-df'       .Or.&
         Trim(simulation_data%dft%vdw%type) == 'optb88'   .Or.&
         Trim(simulation_data%dft%vdw%type) == 'optpbe'   .Or.&
         Trim(simulation_data%dft%vdw%type) == 'optb86b'  .Or.&
         Trim(simulation_data%dft%vdw%type) == 'vdw-df2'      .Or.&
         Trim(simulation_data%dft%vdw%type) == 'vdw-df2-b86r' .Or.&
         Trim(simulation_data%dft%vdw%type) == 'scan+rvv10'   ) Then        
         simulation_data%dft%need_vdw_kernel=.True.
         If (Trim(simulation_data%dft%vdw%type) == 'vdw-df'       .Or.&
           Trim(simulation_data%dft%vdw%type) == 'optb88'   .Or.&
           Trim(simulation_data%dft%vdw%type) == 'optpbe'   .Or.&
           Trim(simulation_data%dft%vdw%type) == 'optb86b'  .Or.&
           Trim(simulation_data%dft%vdw%type) == 'vdw-df2-b86r' .Or.&
           Trim(simulation_data%dft%vdw%type) == 'vdw-df2') Then
           simulation_data%dft%vdw_kernel_file= 'vdW_kernel_table.dat' 
         Else If (Trim(simulation_data%dft%vdw%type) == 'scan+rvv10') Then
           simulation_data%dft%vdw_kernel_file= 'rVV10_kernel_table.dat'
         End If
         Inquire(File=Trim(FOLDER_DFT)//'/'//Trim(simulation_data%dft%vdw_kernel_file), Exist=safe)
         If (.not.safe) Then
           Write (message,'(1x,5a)') '***ERROR: Kernel file ', &
                                   & Trim(FOLDER_DFT)//'/'//Trim(simulation_data%dft%vdw_kernel_file),&
                                   & ', needed for "', Trim(simulation_data%dft%vdw%type),&
                                   & '" dispersion corrections, does not exist. Please copy the file and rerun.'
           Call error_stop(message)
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
         Trim(simulation_data%dft%vdw%type) == 'optb86b'  .Or.&
         Trim(simulation_data%dft%vdw%type) == 'vdw-df2'      .Or.&
         Trim(simulation_data%dft%vdw%type) == 'vdw-df2-b86r' .Or.&
         Trim(simulation_data%dft%vdw%type) == 'scan+rvv10'   ) Then
         If (Trim(simulation_data%dft%vdw%type) == 'optpbe') Then 
           simulation_data%dft%xc_version%type ='or'
         Else If (Trim(simulation_data%dft%vdw%type) == 'optb88') Then
           simulation_data%dft%xc_version%type ='bo'
         Else If (Trim(simulation_data%dft%vdw%type) == 'optb86b'  .Or. &
                 Trim(simulation_data%dft%vdw%type) == 'vdw-df2-b86r' ) Then
           simulation_data%dft%xc_version%type ='mk'
         Else If (Trim(simulation_data%dft%vdw%type) == 'vdw-df2') Then
           simulation_data%dft%xc_version%type ='ml'
         Else If (Trim(simulation_data%dft%vdw%type) == 'scan+rvv10') Then
           simulation_data%dft%xc_version%type ='scan'
         End If
         Call info(' ', 1)
         Write (messages(1),'(1x,5a)') '***WARNING: XC_version will be changed to "', Trim(simulation_data%dft%xc_version%type),&
                                 & '" to include set the requested "',  Trim(simulation_data%dft%vdw%type),&
                                 & '" type of dispersion corrections'   
         Call info(messages,1)
      End If
       
       
     End If

    ! Energy cutoff
    If (Trim(simulation_data%dft%encut%units)/='ry') Then
       Write (message,'(2(1x,a))') Trim(error_dft), &
                                   &'Units for directive "energy_cutoff" for CP2K simulations must be in Ry'
       Call error_stop(message)
    End If
    simulation_data%dft%encut%units='Ry'

    ! Orbital transformation
    If (simulation_data%dft%ot%stat) Then
      If (simulation_data%dft%smear%fread) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Setting "smearing" is inconsistent with the&
                                 & requested Orbital Transformation (OT). Please remove it'
        Call error_stop(message)
      End If 
      If (simulation_data%dft%width_smear%fread) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Setting "width_smear" is inconsistent with the&
                                 & requested Orbital Transformation (OT). Please remove it'
        Call error_stop(message)
      End If
      If (simulation_data%dft%bands%fread) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Setting "bands" is inconsistent with the&
                                 & requested Orbital Transformation (OT). Please remove it'
        Call error_stop(message)
      End If 
      If (simulation_data%dft%mixing%fread) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'Setting "mixing_scheme" is inconsistent with the&
                                 & requested Orbital Transformation (OT). Please remove it'
        Call error_stop(message)
      End If 
    Else ! Smearing
      ! Mixing
      If (simulation_data%dft%mixing%fread) Then
        If (Trim(simulation_data%dft%mixing%type)   /= 'kerker'       .And.&
           Trim(simulation_data%dft%mixing%type)    /= 'linear'       .And.&
           Trim(simulation_data%dft%mixing%type)    /= 'broyden'      .And.&
           Trim(simulation_data%dft%mixing%type)    /= 'broyden-2nd'  .And.&
           Trim(simulation_data%dft%mixing%type)    /= 'multisecant'  .And.&
           Trim(simulation_data%dft%mixing%type)    /= 'pulay')   Then
           Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                        &'Invalid specification of directive "mixing_scheme" for VASP. Options are:'
           Write (messages(2),'(1x,a)') '- Kerker'
           Write (messages(3),'(1x,a)') '- Linear'
           Write (messages(4),'(1x,a)') '- Broyden'
           Write (messages(5),'(1x,a)') '- Broyden-2nd'
           Write (messages(6),'(1x,a)') '- Multisecant'
           Write (messages(7),'(1x,a)') '- Pulay'
           Call info(messages, 7)
           Call error_stop(' ')
        End If
      Else
        simulation_data%dft%mixing%type='linear'
      End If

      ! Smearing           
      If (simulation_data%dft%smear%fread) Then 
        If (Trim(simulation_data%dft%smear%type) /= 'window'    .And.&
          Trim(simulation_data%dft%smear%type) /= 'fermi') Then
          Write (messages(1),'(2(1x,a))') Trim(error_dft), 'The required specification of directive&
                                      & "smearing" is not valid in CP2K. Implemented options are:'
          Write (messages(2),'(1x,a)') '- Fermi       (Fermi-Dirac distribution)'
          Write (messages(3),'(1x,a)') '- Window      (Energy window centered at the Fermi level)' 
          Write (messages(4),'(1x,a)') 'IMPORTANT: method "List" is not implemented in ALC_EQCM'
          Call info(messages, 4)
          Call error_stop(' ')
        End If
      Else
         Write (messages(1),'(2(1x,a))') Trim(error_dft), 'By default (or by setting "OT .False.") the electronic&
                              & problem in CP2K is solved via diagonalization.'
         Write (messages(2),'(1x,a)') 'The user must define directive "smearing".' 
         Call info(messages, 2)
         Call error_stop(' ')
      End If
      ! Width of smearing
      If (.Not. simulation_data%dft%width_smear%fread) Then
        simulation_data%dft%width_smear%value= 0.20_wp
        simulation_data%dft%width_smear%units= 'eV'
      End If 
    
      If (Trim(simulation_data%dft%smear%type) == 'fermi') Then
        simulation_data%dft%width_smear%value= simulation_data%dft%width_smear%value/K_to_eV
        simulation_data%dft%width_smear%units= 'K'
      Else If (Trim(simulation_data%dft%smear%type) == 'window') Then
        simulation_data%dft%width_smear%value= simulation_data%dft%width_smear%value/Ha_to_eV
        simulation_data%dft%width_smear%units= 'Hartree'
      End If
      If (.Not. simulation_data%dft%bands%fread) Then 
         Write (messages(1),'(2(1x,a))') Trim(error_dft), 'For simulations in CP2K with diagonalization and smearing, it is&
                                        & necessary to add a reasonable number of extra bands with directive "bands".'
         Call info(messages, 1)
         Call error_stop(' ')
      End If
    End If

   ! precision
    If (simulation_data%dft%precision%fread) Then
      Write (message,'(2(1x,a))') Trim(error_dft), 'For CP2K, "precision" directive is not needed. Please remove it'
      Call error_stop(message)
    End If
   
   ! SCF energy tolerance 
    If (simulation_data%dft%delta_e%fread) Then
      If (Trim(simulation_data%dft%delta_e%units) == 'ev') Then
        simulation_data%dft%delta_e%value=simulation_data%dft%delta_e%value/Ha_to_eV
        simulation_data%dft%delta_e%units='Hartree'
      End If
    End If
    
    ! Checkings related to basis sets
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    If (simulation_data%dft%basis_info%fread) Then
      ! Check BASIS_SET file exists
      Inquire(File=Trim(FOLDER_DFT)//'/BASIS_SET', Exist=safe)
      If (.not.safe) Then
        Write (messages(1),'(1x,3a)') '***ERROR: File BASIS_SET cannot be found in folder ', Trim(FOLDER_DFT),&
                               & '. This file is needed to set the basis parameters for CP2K simulations.'
        Write (messages(2),'(1x,a)') '   The user must set the content of this file using the various&
                                   & basis set from the CP2K repository,'
        Write (messages(3),'(1x,a)') '   depending on the specification for "basis_set" and the participating atomic species'
        Call info(messages, 3)
        Call error_stop(' ')
      End If

      ! Rename the type of basis
      Do i=1, simulation_data%total_tags
        If (Trim(simulation_data%dft%basis_set(i)%type) /= 'sz'  .And. &
           Trim(simulation_data%dft%basis_set(i)%type) /= 'dz'  .And. &
           Trim(simulation_data%dft%basis_set(i)%type) /= 'szp' .And. &
           Trim(simulation_data%dft%basis_set(i)%type) /= 'dzp' .And. &
           Trim(simulation_data%dft%basis_set(i)%type) /= 'tzp' .And. &
           Trim(simulation_data%dft%basis_set(i)%type) /= 'tz2p' ) Then
           Write (messages(1),'(1x,4a)') Trim(error_dft), ' Invalid CP2K basis set specification for tag "', &
                                    & Trim(simulation_data%dft%basis_set(i)%tag), '". Valid options:'
           Write (messages(2),'(1x,a)') 'SZ    (Single Zeta)'
           Write (messages(3),'(1x,a)') 'DZ    (Double Zeta)'
           Write (messages(4),'(1x,a)') 'SZP   (Single Zeta Polarizable)'  
           Write (messages(5),'(1x,a)') 'DZP   (Double Zeta Polarizable)'   
           Write (messages(6),'(1x,a)') 'TZP   (Triple Zeta Polarizable)' 
           Write (messages(7),'(1x,a)') 'TZ2P  (Triple Zeta 2-Polarizable)'
           Call info(messages, 7)
           Call error_stop(' ')
        End If
        ! Transform to the CP2K lab`lling
        If (Trim(simulation_data%dft%basis_set(i)%type) == 'sz')  simulation_data%dft%basis_set(i)%type='SZV'
        If (Trim(simulation_data%dft%basis_set(i)%type) == 'dz')  simulation_data%dft%basis_set(i)%type='DZV'
        If (Trim(simulation_data%dft%basis_set(i)%type) == 'szp') simulation_data%dft%basis_set(i)%type='SZVP'
        If (Trim(simulation_data%dft%basis_set(i)%type) == 'dzp') simulation_data%dft%basis_set(i)%type='DZVP'
        If (Trim(simulation_data%dft%basis_set(i)%type) == 'tzp') simulation_data%dft%basis_set(i)%type='TZVP'
        If (Trim(simulation_data%dft%basis_set(i)%type) == 'tz2p')simulation_data%dft%basis_set(i)%type='TZV2P'
      End Do

      ! Check if element and basis are defined in the BASIS_SET file
      Do i=1, simulation_data%total_tags
        path=Trim(FOLDER_DFT)//'/BASIS_SET'
        If (Trim(simulation_data%dft%xc_base)=='PBE') Then
          basis= Trim(simulation_data%dft%basis_set(i)%type)//'-MOLOPT-SR-GTH '
        Else 
          basis= Trim(simulation_data%dft%basis_set(i)%type)//'-GTH-'//Trim(simulation_data%dft%xc_base)
        End If
        ! Error message, just in case-----------------------------
        Write (messages(1),'(1x,9a)') '*** ERROR in BASIS_SET file: Basis ', Trim(basis), &
                         & ' cannot be found for atomic species "', Trim(simulation_data%dft%basis_set(i)%element),&
                         & '" (for tag "', Trim(simulation_data%dft%basis_set(i)%tag), &
                         & '"), which is required to compute the DFT problem using the "',&
                         & Trim(simulation_data%dft%xc_version%type),'" XC energy functional.'
        Write (messages(3),'(1x,5a)') '    If basis ', Trim(basis), ' for element "',&
                                    & Trim(simulation_data%dft%basis_set(i)%element),& 
                                    &'" is not available, the user should either consider:'
        Write (messages(4),'(1x,3a)') '      i) changing the basis_set for tag "', Trim(simulation_data%dft%basis_set(i)%tag),&
                                             & '" in block &basis_set'
        Write (messages(5),'(1x,a)')  '     ii) changing "XC_level" and "XC_version" directives in SET_EQCM'
        Write (messages(6),'(1x,a)')  '    iii) contacting the CP2K forum for assistance at https://groups.google.com/g/cp2k'
        Write (messages(7),'(1x,a)')  '     iv) generating a new basis set for the element (only if the user knows what to do)'
        
        If (Trim(simulation_data%dft%xc_base)=='PBE') Then
          Write (messages(2),'(1x,3a)') '    Please refer to files BASIS_MOLOPT, BASIS_MOLOPT_UCL&
                                      & and BASIS_MOLOPT_LnPP1 of the CP2K repository and add&
                                      & the basis to ', Trim(FOLDER_DFT), '/BASIS_SET.'
        Else 
          Write (messages(2),'(1x,3a)') '    Please refer to the BASIS_SET file of the CP2K repository&
                                      & and add the basis to ', Trim(FOLDER_DFT), '/BASIS_SET.'
        End If
        ! -------------------------------------------------------
        word=Trim(simulation_data%component(i)%element) 
        exec_grep='grep "'//Trim(word)//' " '//Trim(path)//' > xc.dat'
        Call execute_command_line(exec_grep)
        Call execute_command_line('wc -l xc.dat > nlines.dat')
        Open(Newunit=internal, File='nlines.dat' ,Status='old')
        Read (internal, Fmt=*, iostat=io) nlines
        If (nlines==0) Then
          Close(internal)
          Call execute_command_line('rm  nlines.dat xc.dat')
          Call info(messages, 7)
          Call error_stop(' ')
        End If
        Close(internal)
 
        root=repeat(' ', 256)
        root=Trim(basis)
        loop=.True.
        Open(Newunit=internal, File='xc.dat' ,Status='old')
        j=1
        Do While (j<=nlines .And. loop)
          Read (internal, Fmt='(a)') line
          Do k=1, maxcol 
            word=repeat(' ', 256)
            Read (line,Fmt=*, iostat=io ) (word, l=1,k)
            If (io /= 0 ) Then
              exit
            Else
              If (root==word) Then
                simulation_data%dft%basis_set(i)%basis=root
                loop=.False.
              End If 
            End If
          End Do
          j=j+1 
        End Do
        If (loop) Then
          Call execute_command_line('rm  nlines.dat xc.dat')
          Call info(messages, 7)
          Call error_stop(' ')
        End If
        Close(internal)
        Call execute_command_line('rm  nlines.dat xc.dat')
      End Do
    End If

    !Pseudo potentials
    !!!!!!!!!!!!!!!!!!! 
    ! Check if element is defined in the POTENTIAL file
    If (simulation_data%dft%pp_info%stat) Then 
      Do i=1, simulation_data%total_tags
        path=Trim(pp_path)//Trim(simulation_data%dft%pseudo_pot(i)%file_name)
        potential= 'GTH-'//Trim(simulation_data%dft%xc_base)
        ! Error message, just in case-----------------------------
        Write (messages(1),'(1x,9a)') '*** ERROR in file ', Trim(simulation_data%dft%pseudo_pot(i)%file_name), &
                         & ': Potential ', Trim(potential), &
                         & ' cannot be found for atomic species "', Trim(simulation_data%component(i)%element),&
                         & '", which is required to compute the DFT problem using the "',&
                         & Trim(simulation_data%dft%xc_version%type),&
                         & '" XC energy functional.' 
        Write (messages(2),'(1x,3a)') '    Please refer to the potential files of the CP2K repository.&
                                    & If potential ', Trim(potential), ' is not available, the user should& 
                                    & consider changing the settings for "XC_level" and "XC_version"'
        ! -------------------------------------------------------
        word=Trim(simulation_data%component(i)%element) 
        exec_grep='grep "'//Trim(word)//' " '//Trim(path)//' > xc.dat'
        Call execute_command_line(exec_grep)
        Call execute_command_line('wc -l xc.dat > nlines.dat')
        Open(Newunit=internal, File='nlines.dat' ,Status='old')
        Read (internal, Fmt=*, iostat=io) nlines
        If (nlines==0) Then
          Close(internal)
          Call execute_command_line('rm  nlines.dat xc.dat')
          Call info(messages, 2)
          Call error_stop(' ')
        End If
        Close(internal)

        root=potential
        loop=.True.
        Call get_word_length(root,root_length)
        Open(Newunit=internal, File='xc.dat' ,Status='old')
        j=1
        Do While (j<=nlines .And. loop)
          Read (internal, Fmt='(a)') line
          Do k=1, maxcol
            Read (line,Fmt=*, iostat=io ) (word, l=1,k)
            If (io /= 0 ) Then
              exit
            Else
              If (root(1:root_length)==word(1:root_length)) Then
                loop=.False.
              End If
            End If
          End Do
          j=j+1
        End Do
        If (loop) Then
          Close(internal)
          Call execute_command_line('rm  nlines.dat xc.dat')
          Call info(messages, 2)
          Call error_stop(' ')
        End If
        Close(internal)
        simulation_data%dft%pseudo_pot(i)%potential=potential
        Call execute_command_line('rm  nlines.dat xc.dat')
      End Do
    End If

    ! max_l_orbital   
    If (simulation_data%dft%max_l_orbital%fread) Then
      Write (message,'(2(1x,a))') Trim(error_dft), 'For CP2K, "max_l_orbital" directive is not needed.&
                                & Please remove it'
      Call error_stop(message)
    End If

    ! Total magnetization
    If (simulation_data%dft%total_magnetization%fread) Then
      If (simulation_data%dft%ot%stat) Then
        Write (messages(1),'(1x,2a)') Trim(error_dft), ' Directive "total_magnetization" is not compatible&
                                    & with the requested OT (Orbital Transformation)'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If

   ! Transform U and J from eV to Hartee   
   If (simulation_data%dft%hubbard_info%fread) Then
     Do i=1, simulation_data%total_tags
       simulation_data%dft%hubbard(i)%U=simulation_data%dft%hubbard(i)%U/Ha_to_eV
       simulation_data%dft%hubbard(i)%J=simulation_data%dft%hubbard(i)%J/Ha_to_eV
     End Do
   End If
    
    ! Ions related settings
    !!!!!!!!!!!!!!!!!!!!!!!
    !Relaxation method
    If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      If (Trim(simulation_data%motion%relax_method%type) /= 'cg'    .And. &
        Trim(simulation_data%motion%relax_method%type)  /= 'bfgs'   .And. &
        Trim(simulation_data%motion%relax_method%type)  /= 'lbfgs'  ) Then
        Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                &'Invalid specification of directive "relax_method" for CP2K. Implemented options are:'
        Write (messages(2),'(1x,a)') '- CG    (Conjugate Gradient)'
        Write (messages(3),'(1x,a)') '- BFGS  (Broyden-Fletcher-Goldfarb-Shanno)'
        Write (messages(4),'(1x,a)') '- LBFGS (Linear Broyden-Fletcher-Goldfarb-Shanno)'
        Call info(messages, 4)
        Call error_stop(' ')
      End If
      If ((.Not. simulation_data%motion%change_cell_volume%stat) .And. simulation_data%motion%change_cell_shape%stat) Then
        Write (messages(1),'(1x,4a)') Trim(error_motion), ' Up to ', Trim(date_RELEASE), ', CP2K does not&
                                  & allow to perform optimisations by keeing the volumen fixed and&
                                  & changing shape of the simulation cell.' 
        Write (messages(2),'(1x,a)') 'For cell optimisation, the volume of the cell must always be set to&
                                  & change, independently if the cell shape is allowed to change or not.' 
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End If
     
    ! Force tolerance
    error=.False.
    If (simulation_data%motion%delta_f%fread) Then
      If (Trim(simulation_data%motion%delta_f%units(1)) /= 'hartree' ) Then
        error=.True.
      End If
      If (Trim(simulation_data%motion%delta_f%units(2)) /= 'bohr-1' ) Then
        error=.True.
      End If
    Else
      simulation_data%motion%delta_f%units(1)='Hartree' 
      simulation_data%motion%delta_f%units(2)='Bohr-1'
      simulation_data%motion%delta_f%value(1)= 4.50000000e-4 
    End If
 
    If (error) Then
      Write (messages(1),'(2a)')  Trim(error_motion), ' Invalid units of directive "force_tolerance" for CP2K.&
                                & Units must be "Hartree Bohr-1"'
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
           Trim(simulation_data%motion%ensemble%type) /= 'nph') Then
           Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                    &'Invalid specification of "ensemble" for CP2K. Options are:'
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
        If (Trim(simulation_data%motion%thermostat%type) /= 'ad_langevin'  .And. &
           Trim(simulation_data%motion%thermostat%type) /= 'csvr'  .And. &
           Trim(simulation_data%motion%thermostat%type) /= 'gle'  .And. &
           Trim(simulation_data%motion%thermostat%type) /= 'nose-hoover'  ) Then
          Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                  &'Specification for "thermostat" is not supported by CP2K. Options are:'
          Write (messages(2),'(1x,a)') '- CSVR        (Canonical Sampling through Velocity Rescaling)'
          Write (messages(3),'(1x,a)') '- GLE         (Generalised Langevin Equation)'
          Write (messages(4),'(1x,a)') '- Ad_Langevin (Adaptative Langevin)'
          Write (messages(5),'(1x,a)') '- Nose-Hoover'
          Call info(messages, 5)
          Call error_stop(' ')
        End If
      End If
    End If

    ! Relaxation time for the thermostat
    If (Trim(simulation_data%motion%ensemble%type) == 'nvt' .Or. &
       Trim(simulation_data%motion%ensemble%type) == 'npt' ) Then
      If (.Not. simulation_data%motion%relax_time_thermostat%fread) Then
        If (Trim(simulation_data%motion%thermostat%type) /= 'gle') Then
          Write (messages(1),'(1x,4a)') Trim(error_motion), ' For CP2K, thermostat "', &
                                       Trim(simulation_data%motion%thermostat%type), '" requires the specification&
                                       & of "relax_time_thermostat", which is missing.'
          Call info(messages, 1)
          Call error_stop(' ')
        End If 
      End If
    End If

    ! Relaxation time for the barostat
    If (Trim(simulation_data%motion%ensemble%type) == 'npt' .Or. &
       Trim(simulation_data%motion%ensemble%type) == 'nph') Then 
      If (.Not. simulation_data%motion%relax_time_barostat%fread) Then
        Write (messages(1),'(1x,4a)') Trim(error_motion), ' For "', &
                                     Trim(simulation_data%motion%ensemble%type), '" simulations in CP2K, the&
                                     & user must specify "relax_time_barostat", which is missing.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If
      
    ! Convert pressure to bars
    If (simulation_data%motion%pressure%fread) Then
      simulation_data%motion%pressure%value= 1000.0_wp *simulation_data%motion%pressure%value
      simulation_data%motion%pressure%units='bar'
    End If

    ! &extra_directives are not allowed to generate simulation files for CP2K
    If (simulation_data%extra_info%stat) Then
      Write (messages(1),'(1x,a)') '***ERROR: definition of extra directives is not possible for CP2K.&
                                  & This is due to the block structure of the input-cp2k.dat file.&
                                  & Please remove sub-block &extra_directives.'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
  End Subroutine define_cp2k_settings  
 
  Subroutine print_cp2k_settings(files, net_elements, list_element, list_tag, list_number, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define settings for CP2K directives
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),  Intent(InOut) :: files(:)
    Integer(Kind=wi), Intent(In   ) :: net_elements
    Character(Len=2), Intent(In   ) :: list_element(max_components) 
    Character(Len=8), Intent(In   ) :: list_tag(max_components) 
    Integer(Kind=wi), Intent(In   ) :: list_number(max_components)
    Type(simul_type), Intent(In   ) :: simulation_data

    Integer(Kind=wi)   :: iunit, i, j, k
    Character(Len=256) :: kpoint_sampling, periodic
    Logical            :: loop
    Real(Kind=wp)      :: mag_ini(max_components)

    ! Open FILE_SET_SIMULATION file
    Open(Newunit=files(FILE_SET_SIMULATION)%unit_no, File=files(FILE_SET_SIMULATION)%filename,Status='Replace')
    iunit=files(FILE_SET_SIMULATION)%unit_no 

    Write (iunit,'(a)')  '##############################'
    Write (iunit,'(a)')  '# File generated with ALC_EQCM'
    Write (iunit,'(a)')  '##############################'
    Write (iunit,'(a)') ' '
 
    ! Global definitions
    Write (iunit,'(a)') '##### Global definitions'
    Write (iunit,'(a)') '&GLOBAL'
    Write (iunit,'(2x,a)') 'PROJECT   model'
    If (Trim(simulation_data%simulation%type) == 'md') Then
      Write (iunit,'(2x,a)') 'RUN_TYPE   MD' 
    Else If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      Write (iunit,'(2x,a)') 'RUN_TYPE   GEO_OPT'
    End If
    Write (iunit,'(2x,a)')   'PRINT_LEVEL  MEDIUM' 
    Write (iunit,'(2x,a)')   'WALLTIME     300000' 
    Write (iunit,'(a)')   '&END GLOBAL'
    Write (iunit,'(a)')   ' '
    ! Electronic structure settings
    Write (iunit,'(a)') '##### Settings for DFT, SYSTEM and IONS directives'
    Write (iunit,'(a)') '&FORCE_EVAL'
    Write (iunit,'(2x,a)') 'METHOD  Quickstep'
  
    If (simulation_data%motion%change_cell_volume%stat .Or. &
       simulation_data%motion%change_cell_shape%stat) Then
       Write (iunit,'(2x,a)') 'STRESS_TENSOR   NUMERICAL'
    End If

    Write (iunit,'(2x,a)') ' '
    Write (iunit,'(2x,a)') '######################'
    Write (iunit,'(2x,a)') '### DFT directives ###'
    Write (iunit,'(2x,a)') '######################'
    Write (iunit,'(2x,a)') '&DFT'
    If (simulation_data%dft%spin_polarised%stat) Then
      Write (iunit,'(4x,a)') 'UKS  .True.  # Spin-polarised calculation'
    End If
    If (simulation_data%dft%hubbard_info%fread) Then
      Write (iunit,'(4x,a)') 'PLUS_U_METHOD   MULLIKEN    # DFT+U method with Mulliken charges'
    End If
    If (simulation_data%net_charge%fread) Then
      Write (iunit,'(4x,a,i6)') 'CHARGE  ', Nint(simulation_data%net_charge%value)
    End If
    Write (iunit,'(a)')   ' '
    Write (iunit,'(4x,a)') '#==== Basis set filename' 
    Write (iunit,'(4x,a)') 'BASIS_SET_FILE_NAME  BASIS_SET'
    If (simulation_data%dft%pp_info%stat) Then 
      Write (iunit,'(4x,a)') '#==== Pseudopotential filename' 
      Write (iunit,'(4x,2a)') 'POTENTIAL_FILE_NAME  ', Trim(simulation_data%dft%pseudo_pot(1)%file_name) 
    End If
    Write (iunit,'(a)')   ' '
    Write (iunit,'(4x,a)') '#==== Self-consistency'
    Write (iunit,'(4x,a)') '&QS'
    Write (iunit,'(6x,a,e10.3,a)') 'EPS_DEFAULT', simulation_data%dft%delta_e%value/10000.0_wp
    Write (iunit,'(6x,a)')         'METHOD  GPW'
    If (simulation_data%dft%total_kpoints > 1) Then
      Write (iunit,'(6x,a)')       'EXTRAPOLATION  USE_GUESS  ! required for K-Point sampling'
    End If
    Write (iunit,'(4x,a)') '&END QS'
    Write (iunit,'(4x,a)') ' '
   
    If (.Not. simulation_data%dft%vdw%fread) Then 
      Write (iunit,'(4x,a)') '#==== Exchange and correlation'
      Write (iunit,'(4x,a)') '&XC'
      If (Trim(simulation_data%dft%xc_version%type)  == 'hl'   .Or.&
        Trim(simulation_data%dft%xc_version%type)  == 'pz'     .Or.&
        Trim(simulation_data%dft%xc_version%type)  == 'wigner' .Or.&
        Trim(simulation_data%dft%xc_version%type)  == 'pw92'   .Or.&
        Trim(simulation_data%dft%xc_version%type)  == 'vwn' ) Then
        Write (iunit,'(6x,a)') '&XC_FUNCTIONAL'
        Write (iunit,'(8x,2a)') '# Use the Dirac-Slater exchange functional ', Trim(bib_slater)
        Write (iunit,'(8x,a)')  '&LDA_X'
        Write (iunit,'(8x,a)')  '&LDA_X'
        Write (iunit,'(8x,a)')  '&END LDA_X'
        If (Trim(simulation_data%dft%xc_version%type)  == 'hl') Then
          Write (iunit,'(8x,2a)') '# Hedin-Lundqvist (HL) correlation ', Trim(bib_hl) 
          Write (iunit,'(8x,a)')  '&LDA_C_HL'
          Write (iunit,'(8x,a)')  '&END LDA_C_HL'
        Else If (Trim(simulation_data%dft%xc_version%type)  == 'pz') Then
          Write (iunit,'(8x,2a)') '# Perdew-Zunger (PZ) correlation ', Trim(bib_pz) 
          Write (iunit,'(8x,a)')  '&LDA_C_PZ'
          Write (iunit,'(8x,a)')  '&END LDA_C_PZ'
        Else If (Trim(simulation_data%dft%xc_version%type)  == 'wigner') Then
          Write (iunit,'(8x,2a)') '# Wigner correlation ', Trim(bib_wigner) 
          Write (iunit,'(8x,a)')  '&LDA_C_WIGNER'
          Write (iunit,'(8x,a)')  '&END LDA_C_WIGNER'
        Else If (Trim(simulation_data%dft%xc_version%type)  == 'vwn') Then
          Write (iunit,'(8x,2a)') '# Vosko-Wilk-Nusair (VWN) correlation ', Trim(bib_vwn) 
          Write (iunit,'(8x,a)')  '&LDA_C_VWN'
          Write (iunit,'(8x,a)')  '&END LDA_C_VWN'
        Else If (Trim(simulation_data%dft%xc_version%type)  == 'pw92') Then
          Write (iunit,'(8x,2a)') '# Perdew-Wang 92 (PW92) correlation ', Trim(bib_pw92) 
          Write (iunit,'(8x,a)')  '&PW92'
          Write (iunit,'(8x,a)')  '&END PW92'
        End If
      End If
      If (Trim(simulation_data%dft%xc_version%type)  == 'pade'  ) Then
        Write (iunit,'(6x,2a)') '# LDA-Pade (PADE) XC functional ', Trim(bib_pade) 
        Write (iunit,'(6x,a)')  '&XC_FUNCTIONAL'
        Write (iunit,'(8x,a)')  '&PADE'
        Write (iunit,'(8x,a)')  '&END PADE'
      End If
      If (Trim(simulation_data%dft%xc_version%type)  == 'pw91'  ) Then
        Write (iunit,'(6x,2a)') '# Perdew-Wang 91 (PW91) XC functional ', Trim(bib_pw91) 
        Write (iunit,'(6x,a)')  '&XC_FUNCTIONAL'
        Write (iunit,'(8x,a)')  '&GGA_X_PW91'
        Write (iunit,'(8x,a)')  '&END GGA_X_PW91'
        Write (iunit,'(8x,a)')  '&GGA_C_PW91'
        Write (iunit,'(8x,a)')  '&END GGA_C_PW91'
      End If 
      If (Trim(simulation_data%dft%xc_version%type)  == 'wc'  ) Then
        Write (iunit,'(6x,2a)') '# Wu-Cohen (WC) exchange (no correlation) ', Trim(bib_wc) 
        Write (iunit,'(6x,a)')  '&XC_FUNCTIONAL'
        Write (iunit,'(8x,a)')  '&GGA_X_WC'
        Write (iunit,'(8x,a)')  '&END GGA_X_WC'
      End If 
      If (Trim(simulation_data%dft%xc_version%type)  == 'am05'  ) Then
        Write (iunit,'(6x,2a)') '# Armiento-Mattsson (AM05) XC functional ', Trim(bib_am05)       
        Write (iunit,'(6x,a)')  '&XC_FUNCTIONAL'
        Write (iunit,'(8x,a)')  '&GGA_X_AM05'
        Write (iunit,'(8x,a)')  '&END GGA_X_AM05'
        Write (iunit,'(8x,a)')  '&GGA_C_AM05'
        Write (iunit,'(8x,a)')  '&END GGA_C_AM05'
      End If 
      If (Trim(simulation_data%dft%xc_version%type)  == 'rp') Then
        Write (iunit,'(6x,2a)') '# Hammer-Hansen-Norskov (RP) XC functional ', Trim(bib_rp) 
        Write (iunit,'(6x,a)')  '&XC_FUNCTIONAL'
        Write (iunit,'(8x,a)')  '&GGA_X_RPBE'
        Write (iunit,'(8x,a)')  '&END GGA_X_RPBE'
        Write (iunit,'(8x,a)')  '&GGA_C_PBE'
        Write (iunit,'(8x,a)')  '&END GGA_C_PBE'
      End If
      If (Trim(simulation_data%dft%xc_version%type)  == 'blyp'  ) Then
        Write (iunit,'(6x,2a)') '# Becke-Lee-Young-Parr (BLYP) XC functional ', Trim(bib_blyp)
        Write (iunit,'(6x,a)')  '&XC_FUNCTIONAL'
        Write (iunit,'(8x,a)')  '&BLYP'
        Write (iunit,'(8x,a)')  '&END BLYP'
      End If
      If (Trim(simulation_data%dft%xc_version%type)  == 'xlyp'  ) Then
        Write (iunit,'(6x,2a)') '# Xu-Goddard (XLYP) XC functional ', Trim(bib_xlyp)
        Write (iunit,'(6x,a)')  '&XC_FUNCTIONAL'
        Write (iunit,'(8x,a)')  '&GGA_XC_XLYP'
        Write (iunit,'(8x,a)')  '&END GGA_XC_XLYP'
      End If
      If (Trim(simulation_data%dft%xc_version%type) == 'pbe'     .Or.&
        Trim(simulation_data%dft%xc_version%type)  == 'revpbe' .Or.&
        Trim(simulation_data%dft%xc_version%type)  == 'pbesol' ) Then
          Write (iunit,'(6x,a)') '&XC_FUNCTIONAL'
          Write (iunit,'(8x,a)') '&PBE'
        If (Trim(simulation_data%dft%xc_version%type)  == 'pbe') Then
          Write (iunit,'(10x,2a)') '# Perdew-Burke-Ernzerhof (PBE) XC functional ', Trim(bib_pbe)
          Write (iunit,'(10x,a)')  'PARAMETRIZATION  ORIG'
        Else If (Trim(simulation_data%dft%xc_version%type)  == 'revpbe') Then
          Write (iunit,'(10x,2a)') '# revPBE XC functional ', Trim(bib_revpbe)
          Write (iunit,'(10x,a)')  'PARAMETRIZATION  REVPBE'
        Else If (Trim(simulation_data%dft%xc_version%type)  == 'pbesol' ) Then
          Write (iunit,'(10x,2a)') '# PBE for solids (PBEsol) XC functional ', Trim(bib_pbesol)        
          Write (iunit,'(10x,a)')  'PARAMETRIZATION  PBESOL'
        End If
        Write (iunit,'(8x,a)') '&END PBE'
      End If
      Write (iunit,'(6x,2a)') '&END XC_FUNCTIONAL' 
      Write (iunit,'(4x,a)')  '&END XC'

    ! vdW corrections
    Else
      Write (iunit,'(4x,3a)') '#==== Exchange and correlation + "', Trim(simulation_data%dft%vdw%type),&
                           & '" dispersion corrections'
      Write (iunit,'(4x,a)') '&XC'
      Write (iunit,'(6x,a)') '&XC_FUNCTIONAL'
      ! Pair potentials
      If (Trim(simulation_data%dft%vdw%type) == 'dft-d2' .Or. &
         Trim(simulation_data%dft%vdw%type) == 'dft-d3' .Or. & 
         Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj' ) Then
        If (Trim(simulation_data%dft%xc_version%type) == 'pbe'     .Or.&
          Trim(simulation_data%dft%xc_version%type)  == 'revpbe' .Or.&
          Trim(simulation_data%dft%xc_version%type)  == 'pbesol' ) Then
            Write (iunit,'(8x,a)') '&PBE'
          If (Trim(simulation_data%dft%xc_version%type)  == 'pbe') Then
            Write (iunit,'(10x,2a)') '# Perdew-Burke-Ernzerhof (PBE) XC functional ', Trim(bib_pbe)
            Write (iunit,'(10x,a)')  'PARAMETRIZATION  ORIG'
          Else If (Trim(simulation_data%dft%xc_version%type)  == 'revpbe') Then
            Write (iunit,'(10x,2a)') '# revPBE XC functional ', Trim(bib_revpbe)
            Write (iunit,'(10x,a)') 'PARAMETRIZATION  REVPBE'
          Else If (Trim(simulation_data%dft%xc_version%type)  == 'pbesol' ) Then
            Write (iunit,'(10x,2a)') '# PBE for solids (PBEsol) XC functional ', Trim(bib_pbesol)
            Write (iunit,'(10x,a)') 'PARAMETRIZATION  PBESOL'
          End If
          Write (iunit,'(8x,a)') '&END PBE'
        End If
        If (Trim(simulation_data%dft%xc_version%type)  == 'blyp'  ) Then
          Write (iunit,'(8x,2a)') '# Becke-Lee-Young-Parr (BLYP) XC functional ', Trim(bib_blyp)
          Write (iunit,'(8x,a)')  '&BLYP'
          Write (iunit,'(8x,a)')  '&END BLYP'
        End If
        If (Trim(simulation_data%dft%xc_version%type)  == 'am05'  ) Then
          Write (iunit,'(8x,2a)') '# Armiento-Mattsson (AM05) XC functional ', Trim(bib_am05) 
          Write (iunit,'(8x,a)')  '&GGA_X_AM05'
          Write (iunit,'(8x,a)')  '&END GGA_X_AM05'
          Write (iunit,'(8x,a)')  '&GGA_C_AM05'
          Write (iunit,'(8x,a)')  '&END GGA_C_AM05'
        End If 
        If (Trim(simulation_data%dft%xc_version%type)  == 'rp') Then
          Write (iunit,'(8x,2a)') '# 1) Hammer-Hansen-Norskov (RPBE), only exchange part ', Trim(bib_rp) 
          Write (iunit,'(8x,a)')  '&GGA_X_RPBE'
          Write (iunit,'(8x,a)')  '&END GGA_X_RPBE'
          Write (iunit,'(8x,2a)') '# 2) Perdew-Burke-Ernzerhof (PBE), only the correlation part ', Trim(bib_pbe)
          Write (iunit,'(8x,a)')  '&GGA_C_PBE'
          Write (iunit,'(8x,a)')  '&END GGA_C_PBE'
        End If
        Write (iunit,'(6x,2a)') '&END XC_FUNCTIONAL'

        ! vdW
        Write (iunit,'(6x,2a)') '&VDW_POTENTIAL'
        Write (iunit,'(8x,a)')    'DISPERSION_FUNCTIONAL  PAIR_POTENTIAL'
        Write (iunit,'(8x,a)')    '&PAIR_POTENTIAL'
        If (Trim(simulation_data%dft%vdw%type) == 'dft-d2') Then
          Write (iunit,'(10x,2a)')  '# Grimme DFT-D2 vdW correction ', Trim(bib_dftd2)  
          Write (iunit,'(10x,a)')   'TYPE DFTD2'
          If (Trim(simulation_data%dft%xc_version%type)  == 'blyp') Then
            Write (iunit,'(10x,a)')   'REFERENCE_FUNCTIONAL  BLYP'
          Else
            Write (iunit,'(10x,a)')   'REFERENCE_FUNCTIONAL  PBE'
            If (Trim(simulation_data%dft%xc_version%type)  == 'pbe') Then
              Write (iunit,'(10x,2a)')  'SCALING =  0.75 # Scaling factor for PBE ', Trim(bib_dftd2) 
            Else If (Trim(simulation_data%dft%xc_version%type)  == 'revpbe') Then
              Write (iunit,'(10x,2a)')  'SCALING =  1.25 # Scaling factor for revPBE ', Trim(bib_dftd3)
            Else If (Trim(simulation_data%dft%xc_version%type)  == 'rp') Then
              Write (iunit,'(10x,2a)')  'SCALING =  1.25  # Scaling factor for RPBE(RP) ',  Trim(bib_tunega)
            Else If (Trim(simulation_data%dft%xc_version%type)  == 'pbesol') Then
              Write (iunit,'(10x,a)')  '#SCALING =  '
            Else If (Trim(simulation_data%dft%xc_version%type)  == 'am05') Then
              Write (iunit,'(10x,a)')  '#SCALING =  '
            End If
          End If
        Else If (Trim(simulation_data%dft%vdw%type) == 'dft-d3' .Or. &
                Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then
           If (Trim(simulation_data%dft%vdw%type) == 'dft-d3') Then
             Write (iunit,'(10x,2a)')  '# Grimme D3 vdW correction ', Trim(bib_dftd3)
             Write (iunit,'(10x,2a)')  'TYPE   DFTD3'
           Else If (Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then
             Write (iunit,'(10x,2a)')  '# Grimme D3 with Becke-Jonson damping ', Trim(bib_dftd3bj)
             Write (iunit,'(10x,2a)')  'TYPE   DFTD3(BJ)'
           End If
           If (Trim(simulation_data%dft%xc_version%type)  == 'pbe') Then
             Write (iunit,'(10x,a)')   'REFERENCE_FUNCTIONAL  PBE'
           Else If (Trim(simulation_data%dft%xc_version%type)  == 'blyp') Then
             Write (iunit,'(10x,a)')   'REFERENCE_FUNCTIONAL  BLYP'
           Else If (Trim(simulation_data%dft%xc_version%type)  == 'revpbe') Then
             Write (iunit,'(10x,a)')   'REFERENCE_FUNCTIONAL  revPBE'
           Else If (Trim(simulation_data%dft%xc_version%type)  == 'pbesol') Then
             Write (iunit,'(10x,a)')   'REFERENCE_FUNCTIONAL  PBEsol'
           Else If (Trim(simulation_data%dft%xc_version%type)  == 'rp') Then
             If (Trim(simulation_data%dft%vdw%type) == 'dft-d3') Then
               Write (iunit,'(10x,a)')  '# D3 scaling parameters for the RPBE (RP) functional (obtained from&
                                         & the execution of VASP)'
               Write (iunit,'(10x,2a)')  'D3_SCALING     1.000   0.8720   0.5140 '
             Else If (Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then
               Write (iunit,'(10x,a)')  '# D3 scaling parameters for the RPBE (RP) functional (obtained from&
                                         & the execution of VASP)'
               Write (iunit,'(10x,2a)')  'D3BJ_SCALING   1.0000  0.1820   0.8318   4.0094 '
             End If
           Else If (Trim(simulation_data%dft%xc_version%type)  == 'am05') Then
             If (Trim(simulation_data%dft%vdw%type) == 'dft-d3') Then
               Write (iunit,'(10x,a)')  '# DFT-D3 vdW parameters are not available for AM05 (to March 2022)'
               Write (iunit,'(10x,a)')  '#D3_SCALING'
             Else If (Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then
               Write (iunit,'(10x,a)')  '# DFT-D3-BJ vdW parameters are not available for AM05 (to March 2022)'
               Write (iunit,'(10x,2a)')  '#D3BJ_SCALING'
             End If
           End If
           Write (iunit,'(10x,a)')  'PARAMETER_FILE_NAME   dftd3.dat'
        End If
        Write (iunit,'(8x,a)')  '&END PAIR_POTENTIAL'
      
      ! Non-local potentials
      Else If (Trim(simulation_data%dft%vdw%type) == 'vdw-df'      .Or. &
              Trim(simulation_data%dft%vdw%type) == 'optpbe'  .Or. &
              Trim(simulation_data%dft%vdw%type) == 'optb86b' .Or.&
              Trim(simulation_data%dft%vdw%type) == 'optb88') Then
        If (Trim(simulation_data%dft%vdw%type) == 'vdw-df') Then
          Write (iunit,'(8x,2a)') '### vdW-DF non-local corrections ', Trim(bib_vdwdf)
          Write (iunit,'(8x,2a)') '# 1) Perdew-Burke-Ernzerhof (PBE) functional (only exchange) ', Trim(bib_pbe)
          Write (iunit,'(8x,a)') '&PBE'
          Write (iunit,'(10x,a)')   'PARAMETRIZATION revPBE'
          Write (iunit,'(10x,a)')   'SCALE_C 0.0'
          Write (iunit,'(8x,a)') '&END PBE'
        Else If (Trim(simulation_data%dft%vdw%type) == 'optpbe') Then
          Write (iunit,'(8x,2a)') '### vdW-optPBE non-local corrections ', Trim(bib_optpbe)
          Write (iunit,'(8x,a)') '# 1) Optimised PBE functional (optPBE, only exchange part)' 
          Write (iunit,'(8x,a)') '&GGA_X_OPTPBE_VDW'
          Write (iunit,'(8x,a)') '&END GGA_X_OPTPBE_VDW' 
        Else If (Trim(simulation_data%dft%vdw%type) == 'optb88') Then
          Write (iunit,'(8x,2a)') '### vdW-optB88 non-local corrections ', Trim(bib_optb88)
          Write (iunit,'(8x,a)')  '# 1) Optimised B88 functional (optB88, only exchange part)'
          Write (iunit,'(8x,a)')  '&GGA_X_OPTB88_VDW'
          Write (iunit,'(8x,a)')  '&END GGA_X_OPTB88_VDW' 
        Else If (Trim(simulation_data%dft%vdw%type) == 'optb86b') Then
          Write (iunit,'(8x,2a)') '### vdW-optB86b non-local corrections ', Trim(bib_optb86b)
          Write (iunit,'(8x,2a)') '# 1) Optimised B86b functional (optB86b, only exchange part)'
          Write (iunit,'(8x,a)') '&GGA_X_OPTB86B_VDW'
          Write (iunit,'(8x,a)') '&END GGA_X_OPTB86B_VDW' 
        End If
        Write (iunit,'(8x,2a)') '# 2) Perdew-Wang 92 (PW92), only the correlation part ', Trim(bib_pw92) 
        Write (iunit,'(8x,a)') '&PW92'
        Write (iunit,'(8x,a)') '&END PW92'
        Write (iunit,'(6x,a)')'&END XC_FUNCTIONAL'

        Write (iunit,'(6x,a)') '&VDW_POTENTIAL'
        Write (iunit,'(8x,a)')   'DISPERSION_FUNCTIONAL  NON_LOCAL'
        Write (iunit,'(8x,a)')   '&NON_LOCAL'
        Write (iunit,'(10x,a)')     'TYPE DRSLL'
        Write (iunit,'(10x,2a)')    'KERNEL_FILE_NAME  ', Trim(simulation_data%dft%vdw_kernel_file)   
        Write (iunit,'(10x,a)')     'CUTOFF 50'
        Write (iunit,'(8x,a)')    '&END NON_LOCAL'
      Else If (Trim(simulation_data%dft%vdw%type) == 'vdw-df2' .Or. & 
              Trim(simulation_data%dft%vdw%type) == 'vdw-df2-b86r') Then
        If (Trim(simulation_data%dft%vdw%type) == 'vdw-df2') Then
          Write (iunit,'(8x,2a)') '### vdW-DF2 non-local corrections ', Trim(bib_vdwdf2)
          Write (iunit,'(8x,2a)') '# 1) Refitted Perdew-Wang 86 (only exchange) ', Trim(bib_rpw86)
          Write (iunit,'(8x,a)')  '&GGA_X_RPW86'
          Write (iunit,'(8x,a)')  '&END GGA_X_RPW86'
        Else If (Trim(simulation_data%dft%vdw%type) == 'vdw-df2-b86r') Then
          Write (iunit,'(8x,2a)') '### vdW-DF2-b86r non-local corrections ', Trim(bib_vdwdf2b86r)
          Write (iunit,'(8x, a)') '# 1) Revised Becke 86 with modified gradient correction (for exchange)'
          Write (iunit,'(8x,a)') '&GGA_X_B86_R'
          Write (iunit,'(8x,a)') '&END GGA_X_B86_R'
        End If
        Write (iunit,'(8x,2a)') '# 2) Perdew-Wang 92 (PW92), only the correlation part ', Trim(bib_pw92) 
        Write (iunit,'(8x,a)') '&PW92'
        Write (iunit,'(8x,a)') '&END PW92'
        Write (iunit,'(6x,a)')'&END XC_FUNCTIONAL'

        Write (iunit,'(6x,a)') '&VDW_POTENTIAL'
        Write (iunit,'(8x,a)')   'DISPERSION_FUNCTIONAL  NON_LOCAL'
        Write (iunit,'(8x,a)')   '&NON_LOCAL'
        Write (iunit,'(10x,a)')     'TYPE LMKLL'
        Write (iunit,'(10x,2a)')    'KERNEL_FILE_NAME  ', Trim(simulation_data%dft%vdw_kernel_file)   
        Write (iunit,'(10x,a)')     'CUTOFF 50'
        Write (iunit,'(8x,a)')    '&END NON_LOCAL'

      Else If (Trim(simulation_data%dft%vdw%type) == 'scan+rvv10') Then
        Write (iunit,'(8x,2a)') '### SCAN+rVV10 non-local corrections ', Trim(bib_scanrvv10)
        Write (iunit,'(8x,2a)') '# SCAN exchange and correlation ', Trim(bib_scan)
        Write (iunit,'(8x,a)') '&MGGA_X_SCAN'
        Write (iunit,'(8x,a)') '&END MGGA_X_SCAN'
        Write (iunit,'(8x,a)') '&MGGA_C_SCAN'
        Write (iunit,'(8x,a)') '&END MGGA_C_SCAN'

        Write (iunit,'(6x,a)') '&VDW_POTENTIAL'
        Write (iunit,'(8x,a)')   'DISPERSION_FUNCTIONAL  NON_LOCAL'
        Write (iunit,'(8x,a)')   '&NON_LOCAL'
        Write (iunit,'(10x,a)')     'TYPE RVV10'
        Write (iunit,'(10x,2a)')    'KERNEL_FILE_NAME  ', Trim(simulation_data%dft%vdw_kernel_file)   
        Write (iunit,'(10x,a)')     'CUTOFF 50'
        Write (iunit,'(8x,a)')    '&END NON_LOCAL'
      End If
      ! Close vdW block
      Write (iunit,'(6x,2a)') '&END VDW_POTENTIAL'
      Write (iunit,'(4x,a)')  '&END XC'
    End If

    ! XC_GRID block (commented)
    Write (iunit,'(4x,a)')  '#Tricks to speed up the computation of the XC potential&
                          & (only if you know what you are doing)'
    Write (iunit,'(4x,a)')  '#&XC_GRID'
    Write (iunit,'(4x,a)')  '#&END XC_GRID'

    Write (iunit,'(4x,a)')  ' '
    Write (iunit,'(4x,a)') '#==== Settings for the SCF process'
    Write (iunit,'(4x,a)') '&SCF'
    Write (iunit,'(6x,a)')       'SCF_GUESS   ATOMIC'
    Write (iunit,'(6x,a,e10.3)') 'EPS_SCF  ',  simulation_data%dft%delta_e%value/100.0_wp
    Write (iunit,'(6x,a,i4)')    'MAX_SCF  ', simulation_data%dft%scf_steps%value 

    If (simulation_data%dft%ot%stat) Then
      Write (iunit,'(6x,a)') '#== Orbital Transformation (OT) settings'
      Write (iunit,'(6x,a)') '&OT'
      Write (iunit,'(8x,a)')   'PRECONDITIONER  FULL_SINGLE_INVERSE'
      Write (iunit,'(8x,a)')   'MINIMIZER  DIIS'
      Write (iunit,'(8x,a)')   'N_DIIS 7'
      Write (iunit,'(6x,a)') '&END OT'
      Write (iunit,'(6x,a)') '&OUTER_SCF'
      Write (iunit,'(8x,a)')   'MAX_SCF 10'
      Write (iunit,'(6x,a)') '&END OUTER_SCF'
    Else
      Write (iunit,'(6x,a)') '#== Diagonalisation'
      Write (iunit,'(6x,a,i4)') 'ADDED_MOS ', simulation_data%dft%bands%value
      Write (iunit,'(6x,a)') '&DIAGONALIZATION'
      Write (iunit,'(8x,a)') 'ALGORITHM  STANDARD'
      Write (iunit,'(6x,a)') '&END DIAGONALIZATION'
      Write (iunit,'(6x,a)') '#== Smearing'
      Write (iunit,'(6x,a)') '&SMEAR'
      If (simulation_data%dft%smear%type=='window') Then
        Write (iunit,'(8x,a)') 'METHOD  ENERGY_WINDOW' 
        Write (iunit,'(8x,a,f12.5)') 'WINDOW_SIZE', simulation_data%dft%width_smear%value
      Else If (simulation_data%dft%smear%type=='fermi') Then
        Write (iunit,'(8x,a)') 'METHOD  FERMI_DIRAC'
        Write (iunit,'(8x,a,f12.5)') 'ELECTRONIC_TEMPERATURE ', simulation_data%dft%width_smear%value
      End If
      If (simulation_data%dft%total_magnetization%fread) Then 
        Write (iunit,'(8x,a,f12.5)') 'FIXED_MAGNETIC_MOMENT  ', simulation_data%dft%total_magnetization%value
      End If
      Write (iunit,'(6x,a)') '&END SMEAR'
      Write (iunit,'(6x,a)') '#== Mixing'
      Write (iunit,'(6x,a)') '&MIXING'
      If (Trim(simulation_data%dft%mixing%type)        == 'kerker') Then     
        Write (iunit,'(8x,a)') 'METHOD  KERKER_MIXING'
      Else If (Trim(simulation_data%dft%mixing%type)    == 'linear') Then      
        Write (iunit,'(8x,a)') 'METHOD DIRECT_P_MIXING'
      Else If (Trim(simulation_data%dft%mixing%type)    == 'broyden') Then     
        Write (iunit,'(8x,a)') 'METHOD BROYDEN_MIXING'
      Else If (Trim(simulation_data%dft%mixing%type)    == 'broyden-2nd') Then 
        Write (iunit,'(8x,a)') 'METHOD BROYDEN_MIXING_NEW'
      Else If (Trim(simulation_data%dft%mixing%type)    == 'multisecant') Then 
        Write (iunit,'(8x,a)') 'METHOD MULTISECANT_MIXING'
      Else If (Trim(simulation_data%dft%mixing%type)    == 'pulay') Then
        Write (iunit,'(8x,a)') 'METHOD PULAY_MIXING'
      End If
      Write (iunit,'(6x,a)') '&END MIXING'
    End If
    Write (iunit,'(4x,a)') '&END SCF'

    !k-points
    If (.Not. simulation_data%dft%ot%stat) Then
      If (simulation_data%dft%total_kpoints > 1) Then
        Write (iunit,'(a)') '  '
        Write (iunit,'(4x,a)') '#== k-point sampling'
        Write (iunit,'(4x,a)') '&KPOINTS'
        If (Trim(simulation_data%dft%kpoints%tag)=='mpack') Then
          kpoint_sampling='MONKHORST-PACK'
        Else If (Trim(simulation_data%dft%kpoints%tag)=='automatic') Then
          kpoint_sampling='GENERAL'
        End If
        Write (iunit,'(6x,2a,3i3)') 'SCHEME  ', Trim(kpoint_sampling),&
                                 & (simulation_data%dft%kpoints%value(i), i= 1, 3)
        Write (iunit,'(6x,a)') 'SYMMETRY  TRUE'
        Write (iunit,'(6x,a)') 'WAVEFUNCTIONS  REAL' 
        Write (iunit,'(6x,a)') 'FULL_GRID  .TRUE.'
        Write (iunit,'(4x,a)') '&END KPOINTS' 
      End If
    End If

    Write (iunit,'(4x,a)') ' '
    Write (iunit,'(4x,a)') '#==== Settings for multigrid information'
    Write (iunit,'(4x,a)') '&MGRID'
    Write (iunit,'(6x,a)') 'NGRIDS    4'
    Write (iunit,'(6x,a,f12.3)') 'CUTOFF', simulation_data%dft%encut%value
    Write (iunit,'(6x,a,f12.3)') 'REL_CUTOFF', simulation_data%dft%encut%value/6.0
    Write (iunit,'(4x,a)') ' &END MGRID'

    Write (iunit,'(4x,a)') ' '
    Write (iunit,'(4x,a)') '#==== Poisson resolutor'
    Write (iunit,'(4x,a)') '&POISSON'
    If (Trim(simulation_data%process)=='electrodeposition') Then
      Write (iunit,'(6x,a)') 'POISSON_SOLVER  ANALYTIC'
      If (Trim(simulation_data%normal_vector)=='c1') Then
        periodic='YZ'
      Else If (Trim(simulation_data%normal_vector)=='c2') Then
        periodic= 'XZ'
      Else If (Trim(simulation_data%normal_vector)=='c3') Then
        periodic= 'XY'
      End If
      Write (iunit,'(6x,2a)') 'PERIODIC  ', Trim(periodic)
    Else
      periodic= 'XYZ' 
      Write (iunit,'(6x,a)') 'POISSON_SOLVER  PERIODIC'
      Write (iunit,'(6x,2a)') 'PERIODIC ', Trim(periodic)
    End If
    Write (iunit,'(4x,a)') '&END POISSON'
    Write (iunit,'(a)')    ' '
    Write (iunit,'(2x,a)') '&END DFT'
    Write (iunit,'(2x,a)') ' '
    Write (iunit,'(2x,a)') '#########################'
    Write (iunit,'(2x,a)') '### SYSTEM directives ###'
    Write (iunit,'(2x,a)') '#########################'
    Write (iunit,'(2x,a)') '&SUBSYS'
    Write (iunit,'(2x,a)') ' '
    Write (iunit,'(4x,a)') '#==== Simulation cell'
    Write (iunit,'(4x,a)') '&CELL'
    Write (iunit,'(6x,2a)')  'PERIODIC ', Trim(periodic) 
    Write (iunit,'(6x,a)')   'SYMMETRY  NONE'
    Write (iunit,'(6x,a,3f12.6)')   'A ', (simulation_data%cell(1,i), i= 1,3)
    Write (iunit,'(6x,a,3f12.6)')   'B ', (simulation_data%cell(2,i), i= 1,3)
    Write (iunit,'(6x,a,3f12.6)')   'C ', (simulation_data%cell(3,i), i= 1,3)
    Write (iunit,'(4x,a)') '&END CELL'
    Write (iunit,'(4x,a)') ' '
    Write (iunit,'(4x,a)') '#==== Topology'
    Write (iunit,'(4x,a)') '&TOPOLOGY'
    Write (iunit,'(6x,a)') 'COORD_FILE_NAME    SAMPLE.cp2k'
    Write (iunit,'(6x,a)') 'COORD_FILE_FORMAT  XYZ'
    Write (iunit,'(4x,a)') '&END TOPOLOGY'
    Write (iunit,'(4x,a)') ' '
    Write (iunit,'(4x,a)') '#==== Description for each atomic species'

    Do k=1, net_elements
      Write (iunit,'(4x,2a)') '&KIND  ', Trim(list_tag(k))
      Write (iunit,'(6x,2a)')  'ELEMENT ', Trim(list_element(k))

      j=1
      loop=.True.
      If (simulation_data%dft%basis_info%fread) Then
        Do While (j <= simulation_data%total_tags .And. loop)
          If (Trim(simulation_data%dft%basis_set(j)%tag)==Trim(list_tag(k))) Then
             Write (iunit,'(6x,2a)')  'BASIS_SET ', Trim(simulation_data%dft%basis_set(j)%basis)
             loop=.False.
          End If
          j=j+1 
        End Do
      End If

      If (simulation_data%dft%pp_info%stat) Then
        j=1
        loop=.True.
        Do While (j <= simulation_data%total_tags .And. loop)
          If (Trim(simulation_data%dft%pseudo_pot(j)%tag)==Trim(list_tag(k))) Then
            Write (iunit,'(6x,2a)')  'POTENTIAL ', Trim(simulation_data%dft%pseudo_pot(j)%potential)
            loop=.False.
          End If
          j=j+1 
        End Do
      End If

      If (simulation_data%dft%mag_info%fread) Then
        j=1
        loop=.True.
        Do While (j <= simulation_data%total_tags .And. loop)
          If (Trim(simulation_data%dft%magnetization(j)%tag)==Trim(list_tag(k))) Then
            Write (iunit,'(6x,a,f8.3)') 'MAGNETIZATION ', simulation_data%dft%magnetization(j)%value 
            loop=.False.
          End If 
          j=j+1 
        End Do
      End If

      If (simulation_data%motion%mass_info%fread) Then
        j=1
        loop=.True.
        Do While (j <= simulation_data%total_tags .And. loop)
          If (Trim(simulation_data%motion%mass(j)%tag)==Trim(list_tag(k))) Then
            Write (iunit,'(6x,a,f8.3)') 'MASS ', simulation_data%motion%mass(j)%value
            loop=.False.
          End If
          j=j+1
        End Do
      End If

      If (simulation_data%dft%hubbard_info%fread) Then
        j=1
        loop=.True.
        Do While (j <= simulation_data%total_tags .And. loop)
          If (Trim(simulation_data%dft%hubbard(j)%tag)==Trim(list_tag(k))) Then
            If (Abs(simulation_data%dft%hubbard(j)%U) > epsilon(1.0_wp) .Or. &
               Abs(simulation_data%dft%hubbard(j)%J) > epsilon(1.0_wp)) Then
              Write (iunit,'(6x,a)') '&DFT_PLUS_U'
              Write (iunit,'(8x,a,i2)')   'L ', simulation_data%dft%hubbard(j)%l_orbital 
              Write (iunit,'(8x,a,f8.3,a)') 'U_minus_J ', & 
                                      & simulation_data%dft%hubbard(j)%U-simulation_data%dft%hubbard(j)%J , &
                                      & ' # (in Hartree)'
              loop=.False.
            End If 
          End If 
          j=j+1 
        End Do
        If (.Not. loop) Then
          Write (iunit,'(6x,a)') '&END DFT_PLUS_U'
        End If
      End If
      Write (iunit,'(4x,a)')  '&END KIND'
    End Do

    If (simulation_data%dft%total_magnetization%fread) Then
       Call check_initial_magnetization(net_elements, list_tag, list_number, mag_ini,&
                                      & simulation_data%dft%total_magnetization%value)
    End If

    Write (iunit,'(a)')    ' '
    Write (iunit,'(2x,a)') '&END SUBSYS'
    Write (iunit,'(a)')    ' '
    Write (iunit,'(a)') '&END FORCE_EVAL'

    Write (iunit,'(a)')    ' '
    Write (iunit,'(a)') '#######################'
    Write (iunit,'(a)') '### IONS directives ###'
    Write (iunit,'(a)') '#######################'
    Write (iunit,'(a)') '&MOTION'
    If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      If (simulation_data%motion%change_cell_volume%stat .Or. &
        simulation_data%motion%change_cell_shape%stat) Then
        Write (iunit,'(2x,a)') '  '
        Write (iunit,'(2x,a)') '#==== Cell and geometry optimization'
        Write (iunit,'(2x,a)') '&CELL_OPT'
        Write (iunit,'(4x,a)')        'TYPE   DIRECT_CELL_OPT'
        Write (iunit,'(4x,a,i4)')     'MAX_ITER ', simulation_data%motion%ion_steps%value
        Write (iunit,'(4x,2a)')       'OPTIMIZER  ', Trim(simulation_data%motion%relax_method%type)
        Write (iunit,'(4x,a,e10.3)')  'MAX_FORCE  ', simulation_data%motion%delta_f%value(1)
        Write (iunit,'(4x,a,9e10.2)') 'EXTERNAL_PRESSURE  ',&
                                      & simulation_data%motion%pressure%value, 0.0_wp, 0.0_wp,& 
                                      & 0.0_wp, simulation_data%motion%pressure%value, 0.0_wp,& 
                                      & 0.0_wp, 0.0_wp, simulation_data%motion%pressure%value 
        Write (iunit,'(4x,a)')        'PRESSURE_TOLERANCE  1.00E+2'
        If (.Not. simulation_data%motion%change_cell_shape%stat) Then
           Write (iunit,'(4x,a)')    'KEEP_ANGLES'
        End If
        Write (iunit,'(2x,a)') '&END CELL_OPT'
      Else
        Write (iunit,'(2x,a)') '  '
        Write (iunit,'(2x,a)') '#==== Geometry optimization'
        Write (iunit,'(2x,a)') '&GEO_OPT'
        Write (iunit,'(4x,a,i7)')    'MAX_ITER ', simulation_data%motion%ion_steps%value
        Write (iunit,'(4x,2a)')      'OPTIMIZER  ', Trim(simulation_data%motion%relax_method%type)
        Write (iunit,'(4x,a,e10.3)') 'MAX_FORCE  ', simulation_data%motion%delta_f%value(1)
        Write (iunit,'(2x,a)') '&END GEO_OPT'
      End If
    ElseIf (Trim(simulation_data%simulation%type) == 'md') Then  
        Write (iunit,'(2x,a)') '  '
        Write (iunit,'(2x,a)') '#==== MD settings'
        Write (iunit,'(2x,a)') '&MD'
        If (Trim(simulation_data%motion%ensemble%type)=='npt') Then
          Write (iunit,'(4x,a)')      'ENSEMBLE    NPT_F'
        Else If (Trim(simulation_data%motion%ensemble%type)=='nph') Then
          Write (iunit,'(4x,a)')      'ENSEMBLE    NPE_F'
        Else
          Write (iunit,'(4x,2a)')     'ENSEMBLE    ', Trim(simulation_data%motion%ensemble%type)
        End If
        Write (iunit,'(4x,a,i7)')   'STEPS     ', simulation_data%motion%ion_steps%value
        Write (iunit,'(4x,a,f7.2)') 'TIMESTEP     ', simulation_data%motion%timestep%value
        Write (iunit,'(4x,a,f7.2)') 'TEMPERATURE  ', simulation_data%motion%temperature%value
        If (Trim(simulation_data%motion%ensemble%type) == 'nvt' .Or. &
           Trim(simulation_data%motion%ensemble%type) == 'npt') Then
          Write (iunit,'(4x,a)') '&THERMOSTAT'
          Write (iunit,'(6x,a)') 'REGION  MASSIVE'
          If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover') Then
            Write (iunit,'(6x,2a)') 'TYPE   NOSE'
            Write (iunit,'(6x,a)') '&NOSE'
            Write (iunit,'(8x,a)') 'LENGTH   3'
            Write (iunit,'(8x,a)') 'YOSHIDA  3'
            Write (iunit,'(8x,a)') 'MTS      2'
            Write (iunit,'(8x,a,f8.2)') 'TIMECON  ', simulation_data%motion%relax_time_thermostat%value
            Write (iunit,'(6x,a)') '&END NOSE'
          ElseIf (Trim(simulation_data%motion%thermostat%type) == 'gle') Then
            Write (iunit,'(6x,2a)') 'TYPE  GLE'
            Write (iunit,'(6x,a)') '&GLE'
            Write (iunit,'(6x,a)') '&END GLE'
          ElseIf (Trim(simulation_data%motion%thermostat%type) == 'ad_langevin') Then
            Write (iunit,'(6x,2a)') 'TYPE  D_LANGEVIN'
            Write (iunit,'(6x,a)') '&AD_LANGEVIN'
            Write (iunit,'(8x,a,f8.2)') 'TIMECON_LANGEVIN  ', simulation_data%motion%relax_time_thermostat%value
            Write (iunit,'(6x,a)') '&END AD_LANGEVIN'
          ElseIf (Trim(simulation_data%motion%thermostat%type) == 'csvr') Then
            Write (iunit,'(6x,2a)') 'TYPE  CSVR'
            Write (iunit,'(6x,a)') '&CSVR'
            Write (iunit,'(8x,a,f8.2)') 'TIMECON  ', simulation_data%motion%relax_time_thermostat%value
            Write (iunit,'(6x,a)') '&END CSVR'
          End If
          Write (iunit,'(4x,a)') '&END THERMOSTAT'
          If (Trim(simulation_data%motion%ensemble%type) == 'npt') Then
            Write (iunit,'(4x,a)') '&BAROSTAT'
            Write (iunit,'(6x,a,f8.2)') 'PRESSURE [bar] ', simulation_data%motion%pressure%value
            Write (iunit,'(6x,a,f10.2)') 'TIMECON  [fs]  ', simulation_data%motion%relax_time_barostat%value
            Write (iunit,'(4x,a)') '&END BAROSTAT'
          End If
        Else If (Trim(simulation_data%motion%ensemble%type) == 'nph') Then
            Write (iunit,'(4x,a)') '&BAROSTAT'
            Write (iunit,'(6x,a,f8.2)') 'PRESSURE [bar] ', simulation_data%motion%pressure%value
            Write (iunit,'(6x,a,f10.2)')'TIMECON  [fs]  ', simulation_data%motion%relax_time_barostat%value
            Write (iunit,'(4x,a)') '&END BAROSTAT'
        End If
        Write (iunit,'(2x,a)') '&END MD'
    End If

    Write (iunit,'(a)')    ' '
    Write (iunit,'(a)') '&END MOTION'

   Close(iunit)

  End Subroutine print_cp2k_settings

  Subroutine warnings_cp2k(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about CP2K settings 
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256)  :: messages(8), header
    Character(Len=256)  :: in_file
    Logical             :: warning, error, print_header
    Integer(Kind=wi)    :: i

    warning=.False.
    print_header=.True.

    Call info(' ', 1)
    Write (messages(1), '(1x,a)') 'In case of problems, or to speed up the electronic convergence,&
                                & the user could try changing:'
    Write (messages(2), '(1x,a)') ' - EPS_DEFAULT and EPS_SCF, which are controlled with directive "scf_energy_tolerance"'
    Write (messages(3), '(1x,a)') ' - CUTOFF and REL_CUTOFF, set to 1.0 and 1/6 the value of "energy_cutoff" (in Ry), respectively'
    Write (messages(4), '(1x,a)') ' - the value of NGRIDS (manually)'  
    ! OT vs magnetization
    If (simulation_data%dft%ot%stat) Then
      Write (messages(5), '(1x,a)')  ' - the parameters in the &OT block (PRECONDITIONER/MINIMIZER)&
                                     & as well as the &OUTER_SCF block (manually)'
      Call info(messages, 5)                               
    Else
      Write (messages(5), '(1x,a)')  ' - the option of directive "mixing_scheme". Add/adjust manually related&
                                     & mixing parameters in block &MIXING'
      Write (messages(6), '(1x,a)')  ' - the settings in block &DIAGONALIZATION (manually)'
      Write (messages(7), '(1x,a)')  ' - the value of ADDED_MOS with directive "bands"' 
      Write (messages(8), '(1x,a)')  ' - the ELECTRONIC_TEMPERATURE using directive "width_smear"'                            
      Call info(messages, 8)                               
    End If

    ! XC_GRID
    Write (messages(1), '(1x,a)') 'Expert users might also use tricks to speed up the calculation&
                                 & of the XC potential by activating block &XC_GRID.'
    Call info(messages, 1)

    Write (messages(1), '(1x,a)') 'The user can manually insert the sub-block "&PRINT" to print different physical quantities&
                                 & either within blocks "&DFT" and/or "&MOTION".'    
    Call info(messages, 1)

    If (.Not. simulation_data%dft%ot%stat) Then
      If (simulation_data%dft%total_kpoints > 1) Then
        Write (messages(1), '(1x,2a)')  '- specifications of the k-point sampling'
        Call info(messages, 1)
      End If
    End If

    If (simulation_data%dft%mag_info%fread  .Or. &
       simulation_data%net_charge%fread  .Or. &    
       simulation_data%dft%vdw%fread) Then
       warning=.True.
    End If
    If (Trim(simulation_data%simulation%type) == 'md') Then
      If (Trim(simulation_data%motion%ensemble%type) == 'nvt' .Or. &
         Trim(simulation_data%motion%ensemble%type) == 'npt' .Or. &
         Trim(simulation_data%motion%ensemble%type) == 'nph' ) Then
         warning=.True.
      End If
    End If

    ! Pulay
    If (simulation_data%motion%change_cell_volume%stat) Then
       warning=.True.
    End If

    If (warning) Then
      in_file='in the model.inp file.'
      Call info(' ', 1)
      Write (header, '(1x,a)')  '***IMPORTANT*** From the requested settings of "&block_simulation_settings", it is&
                                    & RECOMMENDED to consider:'
      ! Total charge
      If (simulation_data%net_charge%fread) Then
        If (Abs(simulation_data%net_charge%value) > epsilon(1.0_wp)) Then
          Write (messages(1), '(1x,a,f8.2,a)')  ' - the meaning/correctness of having set a net charge of ',&
                                  & simulation_data%net_charge%value, ' for the modelled system(s).'
          Call print_warnings(header, print_header, messages, 1)
        End If
      End If

      ! magnetization-related parameters
      If (simulation_data%dft%mag_info%fread) Then
        If (simulation_data%dft%ot%stat) Then
          Write (messages(1), '(1x,a)')  ' - checking the convergence of the magnetic solution with OT,&
                                        & particularly for Ferro-magnetic ordering'
        Else
          Write (messages(1), '(1x,a)')  ' - checking the convergence of the magnetic solution'
        End If
        Call print_warnings(header, print_header, messages, 1)
      End If

      ! Pulay
      If (simulation_data%motion%change_cell_volume%stat) Then
         Write (messages(1), '(1x,a)')  ' - increasing the value for directive "energy_cutoff" to minimise the Pulay stress&
                                        & from changing the cell volume'
         Call print_warnings(header, print_header, messages, 1)
      End If

      ! MD-related parameters
      If (Trim(simulation_data%simulation%type) == 'md') Then
        If (Trim(simulation_data%motion%ensemble%type) == 'nvt' .Or. Trim(simulation_data%motion%ensemble%type) == 'npt') Then
          Write (messages(1), '(1x,a)')  ' - changing (manually) the region the thermostat is attached to (REGION)' 
          Call print_warnings(header, print_header, messages, 1)
          If (Trim(simulation_data%motion%thermostat%type) /= 'gle') Then
            If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover' .Or. &
               Trim(simulation_data%motion%thermostat%type) == 'csvr' ) Then
               Write (messages(1), '(1x,a)')  ' - optimising the setting for "relax_time_thermostat" (TIMECON)'
            Else If (Trim(simulation_data%motion%thermostat%type) == 'ad_langevin') Then
               Write (messages(1), '(1x,a)')  ' - optimising the setting for "relax_time_thermostat" (TIMECON_NH)' 
            End If   
            Call print_warnings(header, print_header, messages, 1)
          End If
          
          If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover') Then
             Write (messages(1), '(1x,a)')  ' - adjusting (manually) the values of LENGTH, YOSHIDA and MTS of the &NOSE block&
                                            & (Nose-Hoover thermostat)' 
             Call print_warnings(header, print_header, messages, 1)
          End If
          If (Trim(simulation_data%motion%ensemble%type) == 'npt') Then
            Write (messages(1), '(1x,a)')  ' - optimising the setting for "relax_time_barostat" (TIMECON in block &BAROSTAT)' 
            Call print_warnings(header, print_header, messages, 1)
          End If
        Else If (Trim(simulation_data%motion%ensemble%type) == 'nph') Then
          Write (messages(1), '(1x,a)')  ' - optimising the setting for "relax_time_barostat" (TIMECON in block &BAROSTAT)'
          Write (messages(2), '(1x,a)')  ' - adding a thermostat to the barostat to control cell fluctuation (only for experts)' 
          Call print_warnings(header, print_header, messages, 2)
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
                                      & must define the correct parameters ', Trim(in_file) 
            Call print_warnings(header, print_header, messages,2) 
          End If

          If (Trim(simulation_data%dft%xc_version%type) /= 'pbe'  .And. &
             Trim(simulation_data%dft%xc_version%type) /= 'blyp' ) Then
            Write (messages(1),'(1x,a)')  ' - revision of the requested DFT-D2 vdW correction:&
                             & please double check the value assigned to SCALING.'
            If (Trim(simulation_data%dft%xc_version%type)  == 'am05') Then
              Write (messages(2),'(1x,3a)')  '   WARNING: To ', Trim(date_RELEASE), ', there is no evidence of previous&
                                        & DFT-D2 simulations with the "AM05" XC functional.&
                                        & IS THE USER CONVINCED ABOUT ADDING VDW CORRECTIONS TO THIS XC FUNCTIONAL?&
                                        & If so, define the correct values in the &PAIR_POTENTIAL block'
            Else If (Trim(simulation_data%dft%xc_version%type)  == 'pbesol') Then
              Write (messages(2),'(1x,3a)')  '   For the requested XC functional "',&
                                        & Trim(simulation_data%dft%xc_version%type),&
                                        & '", the user must define the settings in the &PAIR_POTENTIAL block'
            Else
              Write (messages(2),'(1x,3a)')  '   For the requested XC functional "',&
                                        & Trim(simulation_data%dft%xc_version%type),&
                                        & '", ALC_EQCM has set the corresponding value SCALING. This should be sufficient.&
                                        & Still, the user might want to adjust settings in the &PAIR_POTENTIAL block'
            End If
            Call print_warnings(header, print_header, messages,2) 
          End If

        End If

        If (Trim(simulation_data%dft%vdw%type) == 'dft-d3'   .Or. &
           Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then
          error=.False.
          Do i=1, simulation_data%total_tags
            If (simulation_data%component(i)%atomic_number > 94) Then
              error=.True.
            End If
          End Do
          If (error) Then 
            Write (messages(1),'(1x,3a)')  ' - revision of the requested "', Trim(simulation_data%dft%vdw%type), '" vdW correction:&
                             & the defaults parameters are defined only for elements between H and Pu.'
            Write (messages(2),'(1x,3a)') '   WARNING: at least one of the defined species are beyond this range. The user&
                                      & must define the correct parameters ', Trim(in_file), &
                                      & ' Visit http://www.thch.uni-bonn.de/tc/dftd3 for details'
            Call print_warnings(header, print_header, messages,2) 
          End If

          If (Trim(simulation_data%dft%xc_version%type) /='pbe'    .And. &
             Trim(simulation_data%dft%xc_version%type) /='blyp'   .And. &
             Trim(simulation_data%dft%xc_version%type) /='revpbe' .And. &
             Trim(simulation_data%dft%xc_version%type) /='pbesol' ) Then
             Write (messages(1),'(1x,3a)')  ' - revision of the requested "', Trim(simulation_data%dft%vdw%type),&
                                          & '" vdW correction: the defaults for damping parameters are only available&
                                          & for "PBE", "BLYP", "revPBE" and "PBEsol".'
             If (Trim(simulation_data%dft%xc_version%type) == 'am05') Then
               Write (messages(2),'(1x,3a)')   '   WARNING: To ', Trim(date_RELEASE), ', there is no evidence of previous&
                                        & DFT-D3/DFT-D3-BJ simulations with the "AM05" XC functional.&
                                        & IS THE USER CONVINCED ABOUT ADDING VDW CORRECTIONS TO THIS XC FUNCTIONAL?'
               If (Trim(simulation_data%dft%vdw%type) == 'dft-d3') Then
                 Write (messages(3),'(1x,2a)') '   If so, define the values for D3_SCALING.'
               Else If (Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then
                 Write (messages(3),'(1x,2a)') '   If so, define the values for D3BJ_SCALING.'
               End If

             Else If (Trim(simulation_data%dft%xc_version%type) == 'rp') Then
               Write (messages(2),'(1x,3a)')   '   For the requested XC functional "',&
                                         & Trim(simulation_data%dft%xc_version%type),&
                                         & '", ALC_EQCM has set the corresponding vdW parameters.&
                                         & This is often sufficient.'
               If (Trim(simulation_data%dft%vdw%type) == 'dft-d3') Then
                  Write (messages(3),'(1x,a)') '   Still, the user might want to adjust settings for D3_SCALING.'
               Else If (Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then
                  Write (messages(3),'(1x,a)') '   Still, the user might want to adjust settings for D3BJ_SCALING.'
               End If
             End If
             Write (messages(4),'(1x,a)')   '   Visit http://www.thch.uni-bonn.de/tc/dftd3 for details'
             
             Call print_warnings(header, print_header, messages,4)
          End If
          If (Trim(simulation_data%dft%vdw%type) == 'dft-d3' ) Then
            Write (messages(1),'(1x,a)') ' - vdW-C9 elements are NOT be computed (expensive for large systems). However,&
                                        & the user can change this by'
            Write (messages(2),'(1x,a)') '   setting "REFERENCE_C9_TERM .True." in the block &PAIR_POTENTIAL of the&
                                        & input.cp2k file' 

            Call info(messages,2)
          End If
        End If

        If (Trim(simulation_data%dft%vdw%type) == 'vdw-df'       .Or.&
           Trim(simulation_data%dft%vdw%type) == 'optpbe'   .Or.&
           Trim(simulation_data%dft%vdw%type) == 'optb88'   .Or.&
           Trim(simulation_data%dft%vdw%type) == 'optb86b'  .Or.&
           Trim(simulation_data%dft%vdw%type) == 'vdw-df2'      .Or.&
           Trim(simulation_data%dft%vdw%type) == 'vdw-df2-b86r' .Or.&
           Trim(simulation_data%dft%vdw%type) == 'scan+rvv10'   ) Then
           Write (messages(1),'(1x,3a)')  ' - the following points for the requested "', Trim(simulation_data%dft%vdw%type),&
                                       & '" dispersion correction:'
           Write (messages(2),'(1x,a)')   '   1) CUTOFF directive within &NON_LOCAL block should be optimised&
                                              & for accuracy and efficiency'
           Write (messages(3),'(1x,a)')   '   2) LDA (PADE) potentials could be used in principle (not recommended)'
           Write (messages(4),'(1x,a)')   '   3) this vdW approximation is not defined for spin-polarised systems,&
                                        & but it is still possible to perform spin-polarised simulations'
           Call print_warnings(header, print_header, messages,4)
        End If 
      End If
    End If

  End Subroutine warnings_cp2k

End Module dft_cp2k
