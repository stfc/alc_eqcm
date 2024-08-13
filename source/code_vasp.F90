!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to check, define and print DFT directives for simulations
! with VASP. This module also warns the user about aspects to take
! into consideration when performing simulations
!
! Copyright: 2022-2024 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author    - i.scivetti Oct  2021
! Refact    - i.scivetti Sept 2022
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module code_vasp
  
  Use constants,        Only : date_RELEASE, &
                               max_components

  Use fileset,          Only : file_type, &
                               FILE_SET_EQCM, &
                               FILE_SET_SIMULATION,&
                               FILE_KPOINTS, &
                               FOLDER_DFT 
  Use numprec,          Only : wi, &
                               wp
  Use process_data,     Only : capital_to_lower_case
  Use references,       Only : bib_am05, bib_blyp, bib_ca, bib_ddsc, bib_dftd2, bib_dftd3, bib_dftd3bj, bib_fisher, bib_hl, &
                               bib_mbd, bib_optb86b, bib_optpbe, bib_pbesol, bib_pw91, bib_pz, bib_revpbe, bib_rp, bib_scan, &
                               bib_scanrvv10, bib_ts, bib_tsh, bib_tunega, bib_vdwdf, bib_optb88, bib_pbe, bib_vdwdf2, bib_vwn, &
                               bib_vdwdf2, web_D3BJ, bib_wigner, bib_vdwdf2b86r
  Use simulation_setup, Only : simul_type
  Use simulation_tools, Only : check_extra_directives,&
                               print_extra_directives, &
                               check_initial_magnetization, &                               
                               print_warnings, &
                               record_directive, &
                               scan_extra_directive                                
  Use stoichiometry,    Only : stoich_type
  Use unit_output,      Only : error_stop,&
                               info 

  Implicit None
  Private

  Public :: define_vasp_settings, print_vasp_settings, advise_vasp
  
Contains
  Subroutine define_vasp_settings(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define (and check) settings for VASP directives (set defaults values)
    !
    ! author    - i.scivetti May  2021
    ! refact    - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data


    ! latest vrrsion of the code
    simulation_data%code_version= '5.4.1' 
    
    ! DFT 
    Call define_vasp_dft(files, simulation_data)
    ! motion
    Call define_vasp_motion(files, simulation_data)

  End Subroutine define_vasp_settings

  Subroutine print_vasp_settings(files, net_elements, list_element, list_tag, list_number, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print VASP settings for simulation 
    !
    ! author    - i.scivetti May   2021
    ! refact    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),  Intent(InOut) :: files(:)
    Integer(Kind=wi), Intent(In   ) :: net_elements 
    Character(Len=2), Intent(In   ) :: list_element(max_components)
    Character(Len=8), Intent(In   ) :: list_tag(max_components) 
    Integer(Kind=wi), Intent(In   ) :: list_number(max_components)
    Type(simul_type), Intent(InOut) :: simulation_data

    If (simulation_data%dft%pp_info%stat) Then
      Call print_vasp_potcar(net_elements, list_element, simulation_data)
    End If

    Call print_vasp_kpoints(files, simulation_data)

    Call print_vasp_incar(files, net_elements, list_element, list_tag, list_number, simulation_data)
    
  End Subroutine print_vasp_settings

  Subroutine print_vasp_incar(files, net_elements, list_element, list_tag, list_number, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print VASP settings for simulation 
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),  Intent(InOut) :: files(:)
    Integer(Kind=wi), Intent(In   ) :: net_elements 
    Character(Len=2), Intent(In   ) :: list_element(max_components)
    Character(Len=8), Intent(In   ) :: list_tag(max_components) 
    Integer(Kind=wi), Intent(In   ) :: list_number(max_components)
    Type(simul_type), Intent(InOut) :: simulation_data

    Integer(Kind=wi)  ::  iunit
    Integer(Kind=wi)   :: ic

    ic=1
    
    ! Open FILE_SET_SIMULATION file
    Open(Newunit=files(FILE_SET_SIMULATION)%unit_no, File=files(FILE_SET_SIMULATION)%filename,Status='Replace')
    iunit=files(FILE_SET_SIMULATION)%unit_no

    Write (iunit,'(a)')  '##############################'
    Write (iunit,'(a)')  '# File generated with ALC_EQCM'
    Write (iunit,'(a)')  '##############################'
    Write (iunit,'(a)') ' '
 
    ! DFT part of INCAR
    Call print_vasp_incar_dft(iunit, ic, net_elements, list_element, list_tag, list_number, simulation_data)
    ! Motion part of INCAR
    Call print_vasp_incar_motion(iunit, ic, net_elements, list_tag, simulation_data)

    ! Total number of set directives
    simulation_data%set_directives%N0=ic-1

    If (simulation_data%extra_info%stat) Then
      Call print_extra_directives(iunit, simulation_data%extra_directives, simulation_data%set_directives, &
                            & simulation_data%code_format, simulation_data%simulation%type)
    End If

    Close(iunit)
    
  End Subroutine print_vasp_incar  
    
  Subroutine advise_vasp(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about VASP settings 
    !
    ! author    - i.scivetti May  2021
    ! refact    - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    ! DFT advise
    Call advise_vasp_dft(simulation_data)
    ! motion
    Call advise_vasp_motion(simulation_data)
    ! DFT warnings 
    Call warnings_vasp_dft(simulation_data)
    
  End Subroutine advise_vasp  
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! DFT  
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine define_vasp_dft(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define (and check) DFT settings for VASP directives (set defaults values)
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256)  :: message, messages(16)
    Character(Len=256)  :: error_dft

    Logical             :: safe

    error_dft    = '***ERROR in &dft_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'
    
    ! latest version of the code
    simulation_data%code_version= '5.4.1' 

    ! Check XC_version
    If (Trim(simulation_data%dft%xc_version%type) /= 'ca'     .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'hl'     .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'pz'     .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'wigner' .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'vwn'    .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'am05'   .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'pw91'   .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'pbe'    .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'rp'     .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'revpbe' .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'pbesol' .And.&
      Trim(simulation_data%dft%xc_version%type)   /= 'blyp' ) Then
      Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                  &'Invalid specification for directive "XC_version" for VASP.&
                                  & Implemented options for VASP are:'
      Write (messages(2),'(1x,a)')   '==== LDA-level =================='
      Write (messages(3),'(1x,2a)')  '- CA      Ceperley-Alder         ', Trim(bib_ca)  
      Write (messages(4),'(1x,2a)')  '- HL      Hedin-Lundqvist        ', Trim(bib_hl)
      Write (messages(5),'(1x,2a)')  '- PZ      Perdew-Zunger          ', Trim(bib_pz)
      Write (messages(6),'(1x,2a)')  '- Wigner  Wigner                 ', Trim(bib_wigner)
      Write (messages(7),'(1x,2a)')  '- VWN     Vosko-Wilk-Nusair      ', Trim(bib_vwn)
      Write (messages(8),'(1x,a)')   '==== GGA-level =================='
      Write (messages(9),'(1x,2a)')  '- AM05    Armiento-Mattsson      ', Trim(bib_am05)
      Write (messages(10),'(1x,2a)') '- PW91    Perdew-Wang 91         ', Trim(bib_pw91)
      Write (messages(11),'(1x,2a)') '- PBE     Perdew-Burke-Ernzerhof ', Trim(bib_pbe)
      Write (messages(12),'(1x,2a)') '- RP      Hammer-Hansen-Norskov  ', Trim(bib_rp)
      Write (messages(13),'(1x,2a)') '- revPBE  revPBE                 ', Trim(bib_revpbe)
      Write (messages(14),'(1x,2a)') '- PBEsol  PBE for solids         ', Trim(bib_pbesol)
      Write (messages(15),'(1x,2a)') '- BLYP    Becke-Lee-Young-Parr   ', Trim(bib_blyp)
      Write (messages(16),'(1x,a)')  '================================='
      Call info(messages, 16)
      Call error_stop(' ')
    End If

    
    If (Trim(simulation_data%dft%xc_version%type) == 'ca'     .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'hl'     .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'pz'     .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'wigner' .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'vwn'    ) Then
      simulation_data%dft%xc_base='ca'
    Else If (Trim(simulation_data%dft%xc_version%type)  == 'pw91') Then
      simulation_data%dft%xc_base='pw91'
    Else If (Trim(simulation_data%dft%xc_version%type) == 'am05'   .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'pbe'    .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'rp'     .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'revpbe' .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'pbesol' .Or.&
      Trim(simulation_data%dft%xc_version%type)  == 'blyp' ) Then
      simulation_data%dft%xc_base='pbe'
    End If

    ! Check if basis set was defined, complain and abort
    If (simulation_data%dft%basis_info%fread) Then
      Write (message,'(1x,2a)') Trim(error_dft), &
                        &' Definition of basis sets is not required for VASP.&
                        & Please remove &block_basis_set and rerun.'
      Call error_stop(message)
    End If

    ! Pseudopotentials
    If (simulation_data%dft%pp_info%stat) Then 
      Call check_pseudo_potentials_vasp(simulation_data)
    End If

    ! vdW settings
    !!!!!!!!!!!!!!!!!!!!!!
    simulation_data%dft%need_vdw_kernel=.False.
    If (simulation_data%dft%vdw%fread) Then
      If (Trim(simulation_data%dft%vdw%type) /= 'dft-d2'   .And.&
         Trim(simulation_data%dft%vdw%type) /= 'dft-d3' .And.&
         Trim(simulation_data%dft%vdw%type) /= 'dft-d3-bj' .And.&
         Trim(simulation_data%dft%vdw%type) /= 'ts' .And.&
         Trim(simulation_data%dft%vdw%type) /= 'tsh' .And.&
         Trim(simulation_data%dft%vdw%type) /= 'mbd'       .And.&
         Trim(simulation_data%dft%vdw%type) /= 'ddsc'      .And.&
         Trim(simulation_data%dft%vdw%type) /= 'vdw-df'    .And.&
         Trim(simulation_data%dft%vdw%type) /= 'optpbe'   .And.&
         Trim(simulation_data%dft%vdw%type) /= 'optb88'   .And.&
         Trim(simulation_data%dft%vdw%type) /= 'optb86b'  .And.&
         Trim(simulation_data%dft%vdw%type) /= 'vdw-df2'      .And.&
         Trim(simulation_data%dft%vdw%type) /= 'vdw-df2-b86r' .And.&
         Trim(simulation_data%dft%vdw%type) /= 'scan+rvv10'   ) Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                  &'Invalid specification of directive "vdW" for VASP. Valid options are:'
        Write (messages(2),'(1x,a)')  '- DFT-D2 (method of Grimme)'
        Write (messages(3),'(1x,a)')  '- DFT-D3 (method of Grimme zero damping)'
        Write (messages(4),'(1x,a)')  '- DFT-D3-BJ (method of Grimme, Becke-Jonson damping)'
        Write (messages(5),'(1x,a)')  '- TS   (Tkatchenko-Scheffler method)'
        Write (messages(6),'(1x,a)')  '- TSH  (Tkatchenko-Scheffler method with iterative Hirshfeld partitioning)'
        Write (messages(7),'(1x,a)')  '- MBD  (Many-body dispersion energy method)'
        Write (messages(8),'(1x,a)')  '- dDsC (dispersion correction)'
        Write (messages(9),'(1x,a)')  '- vdW-DF'
        Write (messages(10),'(1x,a)') '- optPBE'
        Write (messages(11),'(1x,a)') '- optB88'
        Write (messages(12),'(1x,a)') '- optB86b'
        Write (messages(13),'(1x,a)') '- vdW-DF2'
        Write (messages(14),'(1x,a)') '- vdW-DF2-B86R'
        Write (messages(15),'(1x,a)') '- SCAN+rVV10'
        Call info(messages, 15)
        Call error_stop(' ')
      End If

      ! vdW only for GGA type
      If (Trim(simulation_data%dft%xc_level%type) /= 'gga') Then
        Write (message,'(1x,4a)') Trim(error_dft), &
                                &' Dispersion correction type "', Trim(simulation_data%dft%vdw%type), '" requires of GGA&
                                & option for directive XC_level. Please change'
        Call error_stop(message)
      End If

      ! DFT-D3 and DFT-D3(BJ) checking for only PBE functionals
      If (Trim(simulation_data%dft%vdw%type) == 'dft-d3'  .Or. &
         Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then
         If (Trim(simulation_data%dft%xc_base)=='pw91') Then
           Write (message,'(1x,4a)') Trim(error_dft), &
                                  &' Dispersion correction type "', Trim(simulation_data%dft%vdw%type), '" requires of PBE&
                                  & pseudo-potentials. Please check the provided pseudo-potential files.'
           Call error_stop(message)
         End If
      End If

      ! Check if the kernel exists 
      If (Trim(simulation_data%dft%vdw%type) == 'vdw-df'       .Or.&
         Trim(simulation_data%dft%vdw%type) == 'optpbe'   .Or.&
         Trim(simulation_data%dft%vdw%type) == 'optb88'   .Or.&
         Trim(simulation_data%dft%vdw%type) == 'optb86b'  .Or.&
         Trim(simulation_data%dft%vdw%type) == 'vdw-df2'      .Or.&
         Trim(simulation_data%dft%vdw%type) == 'vdw-df2-b86r' .Or.&
         Trim(simulation_data%dft%vdw%type) == 'scan+rvv10'   ) Then

        If (Trim(simulation_data%dft%xc_base)=='pw91') Then
          Write (message,'(1x,4a)') Trim(error_dft), &
                                  &' Dispersion correction type "', Trim(simulation_data%dft%vdw%type), '" requires of PBE&
                                  & pseudo-potentials.'
          Call error_stop(message)
        End If
        simulation_data%dft%need_vdw_kernel=.True.
        simulation_data%dft%vdw_kernel_file= 'vdw_kernel.bindat' 
        Inquire(File=Trim(FOLDER_DFT)//'/'//Trim(simulation_data%dft%vdw_kernel_file), Exist=safe)
        If (.not.safe) Then
          Write (message,'(1x,5a)') '***ERROR: Kernel file ',&
                                  & Trim(FOLDER_DFT)//'/'//Trim(simulation_data%dft%vdw_kernel_file), &
                                  & ', needed for "',  Trim(simulation_data%dft%vdw%type), '" dispersion corrections,&
                                  & does not exist. Please copy the file and rerun.'
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

    ! max_l_orbital   
    If (.Not. simulation_data%dft%max_l_orbital%fread) Then
      If (simulation_data%dft%spin_polarised%stat) Then
        Write (message,'(2(1x,a))') Trim(error_dft), 'For VASP settings, spin polarised calculations requires&
                                   & the specification of "max_l_orbital"'
        Call error_stop(message)
      End If
    End If

    ! Orbital transformation
    If (simulation_data%dft%ot%stat) Then
      Write (message,'(2(1x,a))') Trim(error_dft), 'Requested Orbital Transformation via directive "OT"&
                                 & is not possible for VASP simulations. Please remove it'
      Call error_stop(message)
    End If

    ! GAPW 
    If (simulation_data%dft%gapw%stat) Then
      Write (message,'(2(1x,a))') Trim(error_dft), 'Requested "Gaussian Augmented Plane Waves" method via the "gapw"&
                                 & directive is not possible for VASP simulations. Please remove it'
      Call error_stop(message)
    End If

   ! Energy cutoff 
    If (Trim(simulation_data%dft%encut%units)/='ev') Then
       Write (message,'(2(1x,a))') Trim(error_dft), &
                                   &'Units for directive "energy_cutoff" for VASP simulations must be in eV'
       Call error_stop(message)
    End If
      
    ! precision
    If (simulation_data%dft%precision%fread) Then
      If (Trim(simulation_data%dft%precision%type) /= 'low'      .And.&
         Trim(simulation_data%dft%precision%type) /= 'normal'   .And.&
         Trim(simulation_data%dft%precision%type) /= 'accurate' ) Then
         Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                    &'Invalid specification of directive "precision" for VASP simulations. Options are:'
         Write (messages(2),'(1x,a)') '- Low'
         Write (messages(3),'(1x,a)') '- Normal'
         Write (messages(4),'(1x,a)') '- Accurate'
         Call info(messages, 4)
         Call error_stop(' ')
       End If 
       If (Trim(simulation_data%dft%precision%type) =='low' ) Then
         If (simulation_data%dft%vdw%fread) Then
           If (Trim(simulation_data%dft%vdw%type) == 'ts'   .Or.  &
              Trim(simulation_data%dft%vdw%type) == 'tsh'   .Or. &
              Trim(simulation_data%dft%vdw%type) == 'mbd'   .Or. &
             Trim(simulation_data%dft%vdw%type) == 'ddsc') Then    
             Write (messages(1),'(1x,4a)') Trim(error_dft), &
                             & ' requested "', Trim(simulation_data%dft%vdw%type),&
                             & '" dispersion correction necessarily requires of option "Accurate"&
                             & for directive "precision". Please change and rerun.'
             Call info(messages,1)
             Call error_stop(' ')
           End If
         End If
       End If
    Else       
      Write (message,'(2(1x,a))') Trim(error_dft), 'The user must specify directive "precision" for VASP simulation'
      Call error_stop(message)
    End If

   ! Mixing
   If (simulation_data%dft%mixing%fread) Then
     If (Trim(simulation_data%dft%mixing%type)   /= 'kerker'       .And.&
        Trim(simulation_data%dft%mixing%type)    /= 'tchebycheff'  .And.&
        Trim(simulation_data%dft%mixing%type)    /= 'broyden-2nd'  .And.&
        Trim(simulation_data%dft%mixing%type)    /= 'pulay')   Then
        Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                     &'Invalid specification of directive "mixing_scheme" for VASP. Options are:'
        Write (messages(2),'(1x,a)') '- Kerker'
        Write (messages(3),'(1x,a)') '- Tchebycheff'
        Write (messages(4),'(1x,a)') '- Broyden-2nd'
        Write (messages(5),'(1x,a)') '- Pulay'
        Call info(messages, 5)
        Call error_stop(' ')
     End If
   Else
     simulation_data%dft%mixing%type='pulay'
   End If

    ! Smear
    If (simulation_data%dft%smear%fread) Then
      If (Trim(simulation_data%dft%smear%type) /= 'gaussian'    .And.&
         Trim(simulation_data%dft%smear%type) /= 'fermi'       .And.&
         Trim(simulation_data%dft%smear%type) /= 'mp'          .And.&
         Trim(simulation_data%dft%smear%type) /= 'tetrahedron'  ) Then
         Write (messages(1),'(2(1x,a))') Trim(error_dft), &
                                      &'Invalid specification of directive "smearing" for VASP. Options are:'
         Write (messages(2),'(1x,a)') '- Gaussian       (Gaussian distribution)' 
         Write (messages(3),'(1x,a)') '- Fermi          (Fermi-Dirac distribution)'
         Write (messages(4),'(1x,a)') '- Tetrahedron    (Tetrahedron method with Blochl corrections)'
         Write (messages(5),'(1x,a)') '- MP             (Methfessel-Paxton method)'
         Call info(messages, 5)
         Call error_stop(' ')
      End If
    Else
      Write (message,'(2(1x,a))') Trim(error_dft), 'The user must specify directive "smearing" for VASP simulations'
      Call error_stop(message)
    End If 
        
        
    ! Width smearing
    If (Trim(simulation_data%dft%smear%type) /= 'tetrahedron') Then
      If (.Not. simulation_data%dft%width_smear%fread) Then
        simulation_data%dft%width_smear%value=0.20_wp
        simulation_data%dft%width_smear%units='eV'
      Else
        If (Trim(simulation_data%dft%width_smear%units) /= 'ev') Then
           Write (message,'(2a)')  Trim(error_dft), ' Units of directive "width_smear" for VASP must be in eV'
          Call error_stop(message)
        End If
      End If
    End If
    
    ! SCF energy tolerance 
    If (simulation_data%dft%delta_e%fread) Then
      If (Trim(simulation_data%dft%delta_e%units) /= 'ev' ) Then
         Write (message,'(2a)')  Trim(error_dft), ' Units for directive "SCF_energy_tolerance" in VASP must be eV'
         Call info(message, 1)
         Call error_stop(' ')
      End If
    End If

  End Subroutine define_vasp_dft

  Subroutine print_vasp_kpoints(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print KPOINTS file 
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),  Intent(InOut) :: files(:)
    Type(simul_type), Intent(InOut) :: simulation_data

    Integer(Kind=wi)  ::  iunit, i

    ! Open FILE_KPOINTS file
    Open(Newunit=files(FILE_KPOINTS)%unit_no, File=files(FILE_KPOINTS)%filename,Status='Replace')
    iunit=files(FILE_KPOINTS)%unit_no

    Write (iunit,'(a)') 'KPOINTS'
    Write (iunit,'(a)') '0'
    If (Trim(simulation_data%dft%kpoints%tag)=='mpack') Then
      Write (iunit,'(a)') 'Monkhorst-Pack'
    Else If (Trim(simulation_data%dft%kpoints%tag)=='automatic') Then
      Write (iunit,'(a)') 'Automatic'
    End If
    Write (iunit,'(3i3)') (simulation_data%dft%kpoints%value(i), i= 1, 3) 
    Write (iunit,'(3i3)') (                               0, i =1, 3)

    Close(iunit)

  End Subroutine print_vasp_kpoints
  

  Subroutine print_vasp_potcar(net_elements, list_element, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print POTCAR file 
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: net_elements 
    Character(Len=2), Intent(In   ) :: list_element(max_components)
    Type(simul_type), Intent(InOut) :: simulation_data

    Integer(Kind=wi)   :: i, j
    Character(Len=256) :: exec_cat, pp_path
    Character(Len=256) :: pseudo_list(max_components), pseudo_final
    Logical            :: loop
  
    pp_path   = Trim(FOLDER_DFT)//'/PPs/'    
     ! Pseudo potentials
     Do i=1, net_elements
       j=1
       loop=.True.
       Do While (j <= simulation_data%total_tags .And. loop)
         If (Trim(list_element(i))==Trim(simulation_data%dft%pseudo_pot(j)%element)) Then
           pseudo_list(i)= Trim(pp_path)//Trim(simulation_data%dft%pseudo_pot(j)%file_name)
           loop=.False.
         End If
         j=j+1
       End Do
     End Do
     Write (pseudo_final,'(*(1x,a))') (Trim(pseudo_list(i)), i=1, net_elements)
     exec_cat= 'cat '//Trim(pseudo_final)//' > POTCAR'
     Call execute_command_line(exec_cat)
    
  End Subroutine print_vasp_potcar
    
  Subroutine print_vasp_incar_dft(iunit, ic, net_elements, list_element, list_tag, list_number, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print the DFT part of the INCAR file 
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Integer(Kind=wi), Intent(InOut) :: ic 
    Integer(Kind=wi), Intent(In   ) :: net_elements 
    Character(Len=2), Intent(In   ) :: list_element(max_components)
    Character(Len=8), Intent(In   ) :: list_tag(max_components) 
    Integer(Kind=wi), Intent(In   ) :: list_number(max_components)
    Type(simul_type), Intent(InOut) :: simulation_data

    Integer(Kind=wi)  ::  i, j
    Integer(Kind=wi), Dimension(max_components)  ::  l_orbital
    Real(Kind=wp), Dimension(max_components) :: uc, jc  
    Logical  ::  loop
    Character(Len=256) :: mag, nel, mag_list(max_components), mag_final
    Character(Len=256) :: message 
    Real(Kind=wp)      :: number_electrons
    Real(Kind=wp)      :: mag_ini(max_components)

    Write (iunit,'(a)') '##### Parallelization'
    Write (message,'(a,i4,a)') 'NPAR = ', simulation_data%dft%npar%value,  ' #  parellelization of bands'
    Call record_directive(iunit, message, 'NPAR', simulation_data%set_directives%array(ic), ic)
    Write (message,'(a,i4,a)') 'KPAR = ', simulation_data%dft%kpar%value,  ' #  parellelization of k-points'
    Call record_directive(iunit, message, 'KPAR', simulation_data%set_directives%array(ic), ic)
    Write (iunit,'(a)')    '    '
    Write (iunit,'(a)') '##### Electronic structure'
    Write (iunit,'(a)') '#========================='
    Write (message,'(a,i4,a)') 'NELM = ', simulation_data%dft%scf_steps%value, ' # maximum number of SC steps'
    Call record_directive(iunit, message, 'NELM', simulation_data%set_directives%array(ic), ic)
    Write (message,'(a,f6.2,a)') 'ENCUT = '  , simulation_data%dft%encut%value, ' # Energy cutoff'
    Call record_directive(iunit, message, 'ENCUT', simulation_data%set_directives%array(ic), ic)
    Write (message,'(3a)')       'PREC = ', Trim(simulation_data%dft%precision%type),    ' # Precision for calculation'
    Call record_directive(iunit, message, 'PREC', simulation_data%set_directives%array(ic), ic)
    Write (message,'(a,e10.3,a)') 'EDIFF = '  , simulation_data%dft%delta_e%value, ' # Energy tolerance'
    Call record_directive(iunit, message, 'EDIFF', simulation_data%set_directives%array(ic), ic)
    If (Trim(simulation_data%dft%smear%type) == 'gaussian') Then
      Write (message,'(a)')      'ISMEAR = 0   # Gaussian mearing'
    Else If (Trim(simulation_data%dft%smear%type) == 'fermi') Then
      Write (message,'(a)')      'ISMEAR = -1  # Fermi smearing'
    Else If (Trim(simulation_data%dft%smear%type) == 'mp') Then
      Write (message,'(a)')      'ISMEAR =  2  # MP (Methfessel-Paxton order 2)'  
    Else If (Trim(simulation_data%dft%smear%type) == 'tetrahedron') Then
      Write (message,'(a)')      'ISMEAR =  -5  # Tetrahedron'  
    End If
    Call record_directive(iunit, message, 'ISMEAR', simulation_data%set_directives%array(ic), ic)
    If (Trim(simulation_data%dft%smear%type) /= 'tetrahedron') Then
      Write (message,'(a,f6.2,a)') 'SIGMA = '  , simulation_data%dft%width_smear%value, ' # Smearing width in eV'
      Call record_directive(iunit, message, 'SIGMA', simulation_data%set_directives%array(ic), ic)
    End If
    If (simulation_data%dft%spin_polarised%stat) Then
      Write (message,'(a)')      'ISPIN =  2  # Spin-polarised' 
    Else
      Write (message,'(a)')      'ISPIN =  1  # Non spin-polarised' 
    End If
    Call record_directive(iunit, message, 'ISPIN', simulation_data%set_directives%array(ic), ic)

    ! Mixing
    If (Trim(simulation_data%dft%mixing%type)         == 'kerker') Then
       Write (message,'(a)')      'IMIX =  1  # Kerker mixing' 
    Else If (Trim(simulation_data%dft%mixing%type)    == 'tchebycheff') Then 
       Write (message,'(a)')      'IMIX =  2  # Tchebycheff mixing' 
    Else If (Trim(simulation_data%dft%mixing%type)    == 'broyden-2nd') Then
       Write (message,'(a)')      'IMIX =  4 ; WC=0   # Broyden-2nd mixing' 
    Else If (Trim(simulation_data%dft%mixing%type)    == 'pulay')   Then
       Write (message,'(a)')      'IMIX =  4  # Pulay mixing' 
    End If
    Call record_directive(iunit, message, 'IMIX', simulation_data%set_directives%array(ic), ic)

   ! Define the XC part
    If (Trim(simulation_data%dft%xc_version%type) == 'ca') Then
      Write (iunit,'(2a)')     '### Ceperley-Alder (CA) functional ', Trim(bib_ca) 
      Write (message,'(a)')      'GGA =  CA'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
    Else If (Trim(simulation_data%dft%xc_version%type) == 'hl') Then
      Write (iunit,'(2a)')     '### Hedin-Lundqvist (HL) functional ', Trim(bib_hl) 
      Write (message,'(a)')      'GGA =  HL'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
    Else If (Trim(simulation_data%dft%xc_version%type) == 'pz') Then
      Write (iunit,'(2a)')     '### Perdew-Zunger (PZ) functional ', Trim(bib_pz)
      Write (message,'(a)')      'GGA =  PZ'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
    Else If (Trim(simulation_data%dft%xc_version%type) == 'vwn') Then
      Write (iunit,'(2a)')     '### Vosko-Wilk-Nusair functional ', Trim(bib_vwn)
      Write (message,'(a)')      'GGA =  VW'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
    Else If (Trim(simulation_data%dft%xc_version%type) == 'pw91') Then 
      Write (iunit,'(2a)')     '### Perdew-Wang 91 (PW91) XC functional ', Trim(bib_pw91) 
      Write (message,'(a)')      'GGA =  91'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
    Else If (Trim(simulation_data%dft%xc_version%type) == 'am05') Then 
      Write (iunit,'(2a)')     '### Armiento-Mattsson (AM05) functional ', Trim(bib_am05)
      Write (message,'(a)')      'GGA =  AM'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
    Else If (Trim(simulation_data%dft%xc_version%type) == 'pbe') Then
      Write (iunit,'(2a)')     '### Perdew-Burke-Ernzerhof (PBE) XC functional ', Trim(bib_pbe)   
      Write (message,'(a)')      'GGA =  PE'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
    Else If (Trim(simulation_data%dft%xc_version%type) == 'rp') Then
      Write (iunit,'(2a)')     '### Hammer-Hansen-Norskov functional ', Trim(bib_rp) 
      Write (message,'(a)')      'GGA =  RP'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
    Else If (Trim(simulation_data%dft%xc_version%type) == 'revpbe') Then
      Write (iunit,'(2a)')     '### revPBE functional ', Trim(bib_revpbe)
      Write (message,'(a)')      'GGA =  RE'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
    Else If (Trim(simulation_data%dft%xc_version%type) == 'pbesol') Then 
      Write (iunit,'(2a)')     '### PBE for solids (PBEsol) functional ', Trim(bib_pbesol)
      Write (message,'(a)')      'GGA =  PS'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
    Else If (Trim(simulation_data%dft%xc_version%type) == 'blyp') Then
      Write (iunit,'(2a)')     '### Becke-Lee-Young-Parr (BLYP) functional ', Trim(bib_blyp)
      Write (message,'(a)')      'GGA = B5'
      Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
      Write (message,'(a)')      'ALDAX = 1.00'
      Call record_directive(iunit, message, 'ALDAX', simulation_data%set_directives%array(ic), ic)
      Write (message,'(a)')      'AGGAX = 1.00'
      Call record_directive(iunit, message, 'ALDAC', simulation_data%set_directives%array(ic), ic)
      Write (message,'(a)')      'AGGAC = 1.00'
      Call record_directive(iunit, message, 'AGGAX', simulation_data%set_directives%array(ic), ic)
      Write (message,'(a)')      'ALDAC = 0.00'
      Call record_directive(iunit, message, 'AGGAC', simulation_data%set_directives%array(ic), ic)
    End If 
    Write (iunit,'(a)')     '###'
 
   If (Trim(simulation_data%dft%xc_level%type) == 'gga') Then
     Write (message,'(a)')      'GGA_COMPACT =  .FALSE.  # For GGA symmetry' 
     Call record_directive(iunit, message, 'GGA_COMPACT', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a)')      'LASPH       =  .TRUE.   # non-sperical contributions' 
     Call record_directive(iunit, message, 'LASPH', simulation_data%set_directives%array(ic), ic)
   End If
   
   If (simulation_data%large_cell) Then
     Write (message,'(a)') 'AMIN  = 0.01    # Large supercell'
     Call record_directive(iunit, message, 'AMIN', simulation_data%set_directives%array(ic), ic)
   End If   

   If (simulation_data%dft%vdw%fread) Then
     Write (iunit,'(a)') ' '
     Write (iunit,'(a)') '#==== vdW (please check option against VASP version)'
     If (Trim(simulation_data%dft%vdw%type)     == 'dft-d2'   ) Then
       Write (iunit,'(2a)')   '# Method DFT-D2 of Grimme ', Trim(bib_dftd2)
       Write (message,'(a)')   'IVDW = 10'    
       Call record_directive(iunit, message, 'IVDW', simulation_data%set_directives%array(ic), ic)
       If (Trim(simulation_data%dft%xc_version%type) == 'pbe') Then
         Write (message,'(a)')  'VDW_S6 = 0.75   # Scaling factor S6 for PBE'
         Call record_directive(iunit, message, 'VDW_S6', simulation_data%set_directives%array(ic), ic)
       Else If (Trim(simulation_data%dft%xc_version%type) == 'pw91') Then
         Write (message,'(2a)') 'VDW_S6 = 0.7    # Scaling factor S6 for PW91 ', Trim(bib_dftd3)
         Call record_directive(iunit, message, 'VDW_S6', simulation_data%set_directives%array(ic), ic)
       Else If (Trim(simulation_data%dft%xc_version%type) == 'blyp') Then
         Write (message,'(a)')  'VDW_S6 = 1.2    # Scaling factor S6 for BLYP'
         Call record_directive(iunit, message, 'VDW_S6', simulation_data%set_directives%array(ic), ic)
       Else If (Trim(simulation_data%dft%xc_version%type) == 'revpbe') Then
         Write (message,'(2a)') 'VDW_S6 = 1.25   # Scaling factor for revPBE ', Trim(bib_dftd3)
         Call record_directive(iunit, message, 'VDW_S6', simulation_data%set_directives%array(ic), ic)
       Else If (Trim(simulation_data%dft%xc_version%type) == 'pbesol') Then
         Write (message,'(a)') 'VDW_S6 = 1.00   # Scaling factor S6 for PBEsol ', Trim(bib_fisher)
         Call record_directive(iunit, message, 'VDW_S6', simulation_data%set_directives%array(ic), ic)
         Write (message,'(a)') 'VDW_SR = 1.42   # Scaling factor SR'
         Call record_directive(iunit, message, 'VDW_SR', simulation_data%set_directives%array(ic), ic)
       Else If (Trim(simulation_data%dft%xc_version%type) == 'rp') Then 
         Write (message,'(a)') 'VDW_S6 = 1.25   # Scaling factor S6 for RPBE(RP) ', Trim(bib_tunega)
         Call record_directive(iunit, message, 'VDW_S6', simulation_data%set_directives%array(ic), ic)
       End If
     Else If (Trim(simulation_data%dft%vdw%type) == 'dft-d3'   ) Then
       Write (iunit,'(2a)')   '# Method DFT-D3 of Grimme with no damping ', Trim(bib_dftd3)
       Write (message,'(a)')   'IVDW = 11       # zero damping DFT-D3 method of Grimme'
       Call record_directive(iunit, message, 'IVDW', simulation_data%set_directives%array(ic), ic)
       If (Trim(simulation_data%dft%xc_version%type) == 'blyp') Then
         Write (message,'(a)') 'VDW_S6 = 1.000  # DFT-D3 parameters for BLYP'
         Call record_directive(iunit, message, 'VDW_S6', simulation_data%set_directives%array(ic), ic)
         Write (message,'(a)') 'VDW_S8 = 1.682'
         Call record_directive(iunit, message, 'VDW_S8', simulation_data%set_directives%array(ic), ic)
         Write (message,'(a)') 'VDW_SR = 1.094'
         Call record_directive(iunit, message, 'VDW_SR', simulation_data%set_directives%array(ic), ic)
       End If
     Else If (Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then 
       Write (iunit,'(2a)')   '# Method DFT-D3-BJ of Grimme with Becke-Jonson damping ', Trim(bib_dftd3bj)
       Write (message,'(a)')   'IVDW = 12'
       Call record_directive(iunit, message, 'IVDW', simulation_data%set_directives%array(ic), ic)
       If (Trim(simulation_data%dft%xc_version%type) == 'blyp') Then
         Write (message,'(2a)') 'VDW_S6 = 1.0000  # Parameters for BLYP ', Trim(web_D3BJ)
         Call record_directive(iunit, message, 'VDW_S6', simulation_data%set_directives%array(ic), ic)
         Write (message,'(a)')  'VDW_S8 = 2.6996'
         Call record_directive(iunit, message, 'VDW_S8', simulation_data%set_directives%array(ic), ic)
         Write (message,'(a)')  'VDW_A1 = 0.4298'
         Call record_directive(iunit, message, 'VDW_A1', simulation_data%set_directives%array(ic), ic)
         Write (message,'(a)')  'VDW_A2 = 4.2359'
         Call record_directive(iunit, message, 'VDW_A2', simulation_data%set_directives%array(ic), ic)
       End If
     Else If (Trim(simulation_data%dft%vdw%type) == 'ts'       ) Then
       Write (iunit,'(2a)') '# Tkatchenko-Scheffler method (TS) ', Trim(bib_ts)
       Write (message,'(a)')   'IVDW = 20'
       Call record_directive(iunit, message, 'IVDW', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%dft%vdw%type) == 'tsh'      ) Then
       Write (iunit,'(2a)') '# Tkatchenko-Scheffler method with iterative Hirshfeld partitioning ', Trim(bib_tsh)
       Write (message,'(a)')   'IVDW = 21'
       Call record_directive(iunit, message, 'IVDW', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%dft%vdw%type) == 'mbd'      ) Then
       Write (iunit,'(2a)')  '# Many-body dispersion energy method MBD ', Trim(bib_mbd) 
       Write (message,'(a)')   'IVDW = 202'
       Call record_directive(iunit, message, 'IVDW', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%dft%vdw%type) == 'ddsc'     ) Then
       Write (iunit,'(2a)')  '# DDsC dispersion correction method ', Trim(bib_ddsc) 
       Write (message,'(a)')     'IVDW = 4 '
       Call record_directive(iunit, message, 'IVDW', simulation_data%set_directives%array(ic), ic)
       If (Trim(simulation_data%dft%xc_version%type) /= 'pbe'   .And.&
          Trim(simulation_data%dft%xc_version%type) /= 'revpbe') Then
         If (Trim(simulation_data%dft%xc_version%type) == 'blyp') Then
           Write (iunit,'(2a)') '#---BLYP parameters', Trim(bib_ddsc)
           Write (message,'(a)')  'VDW_S6 = 3.6   # a0'
           Call record_directive(iunit, message, 'VDW_S6', simulation_data%set_directives%array(ic), ic)
           Write (message,'(a)')  'VDW_SR = 1.79  # b0'
           Call record_directive(iunit, message, 'VDW_SR', simulation_data%set_directives%array(ic), ic)
         End If
       End If
     Else If (Trim(simulation_data%dft%vdw%type) == 'vdw-df') Then 
        Write (iunit,'(2a)') '# vdW-DF non-local corrections (GGA = RE is defined above) ', Trim(bib_vdwdf)
        Write (message,'(a)')  'LUSE_VDW = .TRUE.'
        Call record_directive(iunit, message, 'LUSE_VDW', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'AGGAC = 0.0000' 
        Call record_directive(iunit, message, 'AGGAC', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%dft%vdw%type) == 'optpbe') Then   
        Write (iunit,'(2a)') '# vdW-optPBE non-local corrections ', Trim(bib_optpbe)
        Write (message,'(a)')  'GGA = OR' 
        Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'LUSE_VDW = .TRUE.' 
        Call record_directive(iunit, message, 'LUSE_VDW', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'AGGAC = 0.0000' 
        Call record_directive(iunit, message, 'AGGAC', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%dft%vdw%type) == 'optb88') Then  
        Write (iunit,'(2a)') '# vdW-optB88 non-local corrections ', Trim(bib_optb88)
        Write (message,'(a)')  'GGA = BO' 
        Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'PARAM1 = 0.1833333333' 
        Call record_directive(iunit, message, 'PARAM1', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'PARAM2 = 0.2200000000' 
        Call record_directive(iunit, message, 'PARAM2', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'LUSE_VDW = .TRUE.' 
        Call record_directive(iunit, message, 'LUSE_VDW', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'AGGAC = 0.0000' 
        Call record_directive(iunit, message, 'AGGAC', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%dft%vdw%type) == 'optb86b') Then 
        Write (iunit,'(2a)') '# vdW-optB86b non-local corrections ', Trim(bib_optb86b)
        Write (message,'(a)')  'GGA = MK' 
        Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'PARAM1 = 0.1234' 
        Call record_directive(iunit, message, 'PARAM1', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'PARAM2 = 1.0000' 
        Call record_directive(iunit, message, 'PARAM2', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'LUSE_VDW = .TRUE.' 
        Call record_directive(iunit, message, 'LUSE_VDW', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'AGGAC = 0.0000' 
        Call record_directive(iunit, message, 'AGGAC', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%dft%vdw%type) == 'vdw-df2') Then    
        Write (iunit,'(2a)') '# vdW-DF2 non-local corrections ', Trim(bib_vdwdf2)
        Write (message,'(a)')  'GGA = ML' 
        Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'LUSE_VDW = .TRUE.' 
        Call record_directive(iunit, message, 'LUSE_VDW', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'Zab_vdW = -1.8867' 
        Call record_directive(iunit, message, 'Zab_vdW', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'AGGAC = 0.0000' 
        Call record_directive(iunit, message, 'AGGAC', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%dft%vdw%type) == 'vdw-df2-b86r') Then
        Write (iunit,'(2a)') '# vdW-DF2-b86r non-local corrections ', Trim(bib_vdwdf2b86r)
        Write (message,'(a)')  'GGA      = MK' 
        Call record_directive(iunit, message, 'GGA', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'LUSE_VDW = .TRUE' 
        Call record_directive(iunit, message, 'LUSE_VDW', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'PARAM1   = 0.1234' 
        Call record_directive(iunit, message, 'PARAM1', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'PARAM2   = 0.711357' 
        Call record_directive(iunit, message, 'PARAM2', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'Zab_vdW  = -1.8867' 
        Call record_directive(iunit, message, 'Zab_vdW', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'AGGAC    = 0.0000' 
        Call record_directive(iunit, message, 'AGGAC', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%dft%vdw%type) == 'scan+rvv10') Then
        Write (iunit,'(2a)') '# SCAN+rVV10 non-local corrections ', Trim(bib_scanrvv10)
        Write (message,'(2a)') 'METAGGA  = SCAN  ', Trim(bib_scan) 
        Call record_directive(iunit, message, 'METAGGA', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'LUSE_VDW = .TRUE' 
        Call record_directive(iunit, message, 'LUSE_VDW', simulation_data%set_directives%array(ic), ic)
        Write (message,'(a)')  'BPARAM = 15.7' 
        Call record_directive(iunit, message, 'METAGGA', simulation_data%set_directives%array(ic), ic)
     End If
     Write (iunit,'(a)') ' '
   End If

   If (simulation_data%dft%bands%fread) Then
     Write (message,'(a,i4,a)') 'NBANDS = ', simulation_data%dft%bands%value,&
                             & ' # total bands to improve convergence'
     Call record_directive(iunit, message, 'NBANDS', simulation_data%set_directives%array(ic), ic)
   End If

   If (simulation_data%net_charge%fread) Then
     Write (iunit,'(a,f8.3)') '#==== Net system charge equal to ', simulation_data%net_charge%value
     number_electrons=0.0_wp
     Do i=1, net_elements
       number_electrons=number_electrons+list_number(i)*simulation_data%component(i)%valence_electrons
     End Do
     If (number_electrons < 1) Then
       Write (message, '(1x,a)') '***ERROR: it appears that the value assigned for "net_charge" is unrealistic.&
                                & Is the user sure about this setting? Please check'
       Call error_stop(message)
     End If
     Write (message,'(a,f20.4)') 'NELECT= ', number_electrons-simulation_data%net_charge%value
     Call record_directive(iunit, message, 'NELECT', simulation_data%set_directives%array(ic), ic)
     Write (iunit,'(a)') ' '
   End If

   If (simulation_data%dft%spin_polarised%stat) Then
     Write (iunit,'(a)') '#==== Additional spin-polarised directives'  
     If (simulation_data%dft%max_l_orbital%fread) Then
       If (simulation_data%dft%max_l_orbital%value<2) Then
         Write (message,'(a)') 'LMAXMIX = 2'
         Call record_directive(iunit, message, 'LMAXMIX', simulation_data%set_directives%array(ic), ic)
       ElseIf (simulation_data%dft%max_l_orbital%value==2) Then
         Write (message,'(a)') 'LMAXMIX = 4'
         Call record_directive(iunit, message, 'LMAXMIX', simulation_data%set_directives%array(ic), ic)
       ElseIf (simulation_data%dft%max_l_orbital%value==3) Then
         Write (message,'(a)') 'LMAXMIX = 6'
         Call record_directive(iunit, message, 'LMAXMIX', simulation_data%set_directives%array(ic), ic)
       End If
     Else
       Write (message,'(a)') 'LMAXMIX = 2'
       Call record_directive(iunit, message, 'LMAXMIX', simulation_data%set_directives%array(ic), ic)
     End If

     If (Trim(simulation_data%dft%xc_version%type) == 'pw91') Then
       Write (message,'(a)') 'VOSKOWN = 1'
       Call record_directive(iunit, message, 'VOSKOWN', simulation_data%set_directives%array(ic), ic)
     End If
     Write (iunit,'(a)') ' '
   End If

   If (simulation_data%dft%mag_info%fread) Then
     Write (iunit,'(a)') '#==== Magnetization'
     If (simulation_data%dft%total_magnetization%fread) Then
       Write (message,'(a,f6.2)') 'NUPDOWN = ', simulation_data%dft%total_magnetization%value
       Call record_directive(iunit, message, 'NUPDOWN', simulation_data%set_directives%array(ic), ic)
     Else
       Write (message,'(a)') 'NUPDOWN = -1'
       Call record_directive(iunit, message, 'NUPDOWN', simulation_data%set_directives%array(ic), ic)
     End If
     mag_list=repeat(' ', 256)
     Do i=1, net_elements 
       j=1
       loop=.True.
       Do While (j <= simulation_data%total_tags .And. loop)
         If (Trim(list_tag(i))==Trim(simulation_data%dft%magnetization(j)%tag)) Then
           mag_ini(i)=simulation_data%dft%magnetization(j)%value
           Write (nel,*) list_number(i)
           Write (mag,'(f5.2)') simulation_data%dft%magnetization(j)%value
           mag_list(i)= Trim(nel)//'*'//Adjustl(Trim(mag))//'  '
           loop=.False.
         End If
         j=j+1
       End Do
     End Do
     Write (mag_final,'(*(1x,a10))') (Adjustl(Trim(mag_list(i))), i=1, net_elements)
     Write (message,'(2a)') 'MAGMOM = ', Trim(mag_final)
     Call record_directive(iunit, message, 'MAGMOM', simulation_data%set_directives%array(ic), ic)
     Write (iunit,'(a)') ' '

     If (simulation_data%dft%total_magnetization%fread) Then
       Call check_initial_magnetization(net_elements, list_tag, list_number, mag_ini,&
                                      & simulation_data%dft%total_magnetization%value)
     End If 

   End If

   If (simulation_data%dft%hubbard_info%fread) Then
     Write (iunit,'(a)') '#==== Hubbard corrections'
     Write (message,'(a)') 'LDAU = .TRUE.'
     Call record_directive(iunit, message, 'LDAU', simulation_data%set_directives%array(ic), ic)
     loop=.False.
     Do i=1, simulation_data%total_tags
       If (Abs(simulation_data%dft%hubbard(i)%J)>epsilon(1.0_wp)) Then 
         loop=.True.
       End If
     End Do
     If (loop) Then
       Write (message,'(a)') 'LDAUTYPE = 1'
       Call record_directive(iunit, message, 'LDAUTYPE', simulation_data%set_directives%array(ic), ic)
     Else
       Write (message,'(a)') 'LDAUTYPE = 2'
       Call record_directive(iunit, message, 'LDAUTYPE', simulation_data%set_directives%array(ic), ic)
     End If 

     Do i=1, net_elements
       j=1
       loop=.True.
       Do While (j <= simulation_data%total_tags .And. loop)
         If (Trim(list_element(i))==Trim(simulation_data%dft%hubbard(j)%element)) Then
           If (Abs(simulation_data%dft%hubbard(j)%U)<epsilon(1.0_wp) .And.&
               Abs(simulation_data%dft%hubbard(j)%J)<epsilon(1.0_wp)) Then
             l_orbital(i)= -1
           Else
             l_orbital(i)= simulation_data%dft%hubbard(j)%l_orbital 
           End If
           uc(i)=simulation_data%dft%hubbard(j)%U 
           jc(i)=simulation_data%dft%hubbard(j)%j 
           loop=.False.
         End If
         j=j+1
       End Do      
     End Do
     Write (message,'(a,*(i6))')   'LDAUL = ', (l_orbital(i), i=1,net_elements)
     Call record_directive(iunit, message, 'LDAUL', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,*(f6.2))') 'LDAUU = ', (uc(i),        i=1,net_elements)
     Call record_directive(iunit, message, 'LDAUU', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,*(f6.2))') 'LDAUJ = ', (jc(i),        i=1,net_elements)
     Call record_directive(iunit, message, 'LDAUJ', simulation_data%set_directives%array(ic), ic)
   End If

  End Subroutine print_vasp_incar_dft   

  Subroutine advise_vasp_dft(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about DFT settings 
    !
    ! author   - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256) :: messages(8)
    Character(Len=256) :: in_extra

    in_extra='using the &extra_directives block'

    Call info(' ', 1)
    Write (messages(1), '(1x,a)')  'The efficiency in the parallelization can be optimised by:'
    Write (messages(2), '(1x,a)')  ' - adjusting NPAR via ALC_EQCM directive "npar".'
    Write (messages(3), '(1x,2a)') ' - setting directives LREAL and/or LPLANE ', Trim(in_extra) 
    Write (messages(4), '(1x,a)')  'In case of problems in the electronic convergence, the user should try:'
    Write (messages(5), '(1x,2a)') ' - changing the mixing parameters (AMIX and BMIX) ', Trim(in_extra)
    Write (messages(6), '(1x,a)')  ' - increasing the value of NBANDS via ALC_EQCM directive "bands"'
    Write (messages(7), '(1x,2a)') ' - changing the settings of MAXMIX, NELMIN and NELMDL ', Trim(in_extra)
    Call info(messages, 7)
    If (Trim(simulation_data%dft%mixing%type)  == 'pulay' .Or. Trim(simulation_data%dft%mixing%type)  == 'broyden-2nd') Then
      Write (messages(1), '(1x,3a)')  ' - changing default values for MIXPRE and INIMIX ', Trim(in_extra), ' (see manual)' 
      Call info(messages, 1)
    End If    
    If (Trim(simulation_data%dft%smear%type) /= 'tetrahedron') Then
      Write (messages(1), '(1x,a)')  ' - increasing the value of SIGMA via ALC_EQCM directive "width_smear"'
      Call info(messages, 1)
    End If 
    Write (messages(1), '(1x,a)')  ' - selecting a different option for "mixing_scheme"'
    Write (messages(2), '(1x,2a)') 'If problems persist, try setting "IALGO = 48" ', Trim(in_extra)
    Call info(messages, 2)

    If (simulation_data%dft%total_kpoints > 1) Then
      If (simulation_data%dft%kpar%value==1) Then
        Write (messages(1), '(1x,a)') 'Since the  Brillouin zone is set to be sampled with more than one k-points,&
                                    & the user might consider changing KPAR changed via directive "kpar".'
        Call info(messages, 1)
      End If
    End If    

    Write (messages(1), '(1x,3a)') 'I/O can be controlled ', Trim(in_extra), ' with the following directives (see VASP manual):'
    Write (messages(2), '(1x,a)')  ' - LWAVE, LCHARG, LVTOT and LVHAR (depending on the system and memory requirements)'
    Write (messages(3), '(1x,a)')  ' - LORBIT (for DOS analysis and printing of magnetic moments)'
    Call info(messages, 3)


  End Subroutine advise_vasp_dft

  Subroutine warnings_vasp_dft(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to recommend the user about CASTEP settings 
    !
    ! author    - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data
    
    Character(Len=256) :: message
    Character(Len=256) :: messages(10), header
    Character(Len=256) :: in_extra
    Logical            :: error
    Integer(Kind=wi)   :: i

    Logical            :: warning, print_header

    warning=.False.
    print_header=.True.
    
    in_extra='using the &extra_directives block'
    
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
                                  & simulation_data%net_charge%value, &
                                  & ' for the modelled system(s).'
          Call print_warnings(header, print_header, messages, 1)
        End If
      End If

      ! Hubbard-related parameters
      If (simulation_data%dft%hubbard_info%fread) Then
        Write (messages(1), '(1x,3a)')  ' - changing the settings for AMIX_MAG and/or BMIX_MAG directives ', Trim(in_extra),&
                                  & ' (if electronic convergence fails from the inclusion of Hubbard corrections)'
        Call print_warnings(header, print_header, messages, 1)
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
                             & the defaults for VDW_C6 and VDW_R0& 
                             & are defined only for elements in the first five rows of periodic table (i.e. H-Xe).'
            Write (messages(2),'(1x,2a)') '   WARNING: at least one of the defined species are beyond this range and the user&
                                      & must define the correct parameters ', Trim(in_extra) 
            Call print_warnings(header, print_header, messages,2) 
          End If
          If (Trim(simulation_data%dft%xc_version%type)  /= 'pbe') Then
            Write (messages(1),'(1x,a)')  ' - revision of the requested DFT-D2 vdW correction:&
                             & the defaults for parameters controlling the damping function (VDW_S6, VDW_SR, VDW_D)&
                             & are available only for the PBE functional.'
            If (Trim(simulation_data%dft%xc_version%type)  == 'am05') Then
              Write (messages(2),'(1x,4a)')  '   WARNING: To ', Trim(date_RELEASE), ', there is no evidence of previous&
                                        & DFT-D2 simulations with the "AM05" XC functional.&
                                        & IS THE USER CONVINCED ABOUT ADDING VDW CORRECTIONS TO THIS XC FUNCTIONAL?&
                                        & If so, define the values for VDW_S6, VDW_SR, VDW_C6 and VDW_R0 ', Trim(in_extra)
            Else
              Write (messages(2),'(1x,5a)')  '   For the requested XC functional "',&
                                        & Trim(simulation_data%dft%xc_version%type),&
                                        & '", ALC_EQCM has set the corresponding value of VDW_S6 ', Trim(in_extra),&
                                        & ' This is often sufficient.'
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
                                      & must define the correct parameters ', Trim(in_extra), &
                                      & ' Visit http://www.thch.uni-bonn.de/tc/dftd3 for details'
            Call print_warnings(header, print_header, messages,2) 
          End If
          If (Trim(simulation_data%dft%xc_version%type) /='pbe'    .And. &
             Trim(simulation_data%dft%xc_version%type) /='rp'     .And. &
             Trim(simulation_data%dft%xc_version%type) /='revpbe' .And. &
             Trim(simulation_data%dft%xc_version%type) /='pbesol' ) Then
             Write (messages(1),'(1x,3a)')  ' - revision of the requested "', Trim(simulation_data%dft%vdw%type),&
                                            '" vdW correction: the defaults for damping parameters are only available&
                                           & for "PBE", "RP", "revPBE" and "PBEsol".'
             If (Trim(simulation_data%dft%xc_version%type) == 'am05') Then
               Write (messages(2),'(1x,3a)')   '   WARNING: To ', Trim(date_RELEASE), ', there is no evidence of previous&
                                        & DFT-D3/DFT-D3-BJ simulations with the "AM05" XC functional.&
                                        & IS THE USER CONVINCED ABOUT ADDING VDW CORRECTIONS TO THIS XC FUNCTIONAL?'
               If (Trim(simulation_data%dft%vdw%type) == 'dft-d3') Then
                 Write (messages(3),'(1x,2a)') '   If so, define the values for VDW_S6, VDW_S8 and VDW_SR ', Trim(in_extra)  
               Else If (Trim(simulation_data%dft%vdw%type) == 'dft-d3-bj') Then
                 Write (messages(3),'(1x,2a)') '   If so, define the values for VDW_S6, VDW_S8, VDW_A1 and VDW_A2 ',&
                                               & Trim(in_extra)  
               End If
               
             Else If (Trim(simulation_data%dft%xc_version%type) == 'blyp') Then
               Write (messages(2),'(1x,3a)')   '   For the requested XC functional "',&
                                         & Trim(simulation_data%dft%xc_version%type),&
                                         & '", ALC_EQCM has set the corresponding vdW parameters in the INCAR file.'
             End If
             Call print_warnings(header, print_header, messages,3)
          End If 
        End If

        If (Trim(simulation_data%dft%vdw%type) == 'ts'   .Or. &
           Trim(simulation_data%dft%vdw%type) == 'tsh'  .Or. &
           Trim(simulation_data%dft%vdw%type) == 'mbd'  .Or. &
           Trim(simulation_data%dft%vdw%type) == 'ddsc' ) Then
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
                                      & must define the correct parameters for all atoms in the system via VDW_ALPHA,&
                                      & VDW_C6 and VDW_R0 ', Trim(in_extra)
            Call print_warnings(header, print_header, messages,2) 
          End If
        End If

        Write (messages(1),'(1x,3a)')  ' - revision of the requested "', Trim(simulation_data%dft%vdw%type), '" vdW correction:'
        Write (messages(2),'(1x,a)')   '   1) this method requires the use of POTCAR files from the PAW dataset&
                                       & version 5.2 or later'
        Write (messages(3),'(1x,a)')   '   2) bear in mind that the charge-density dependence of gradients is neglected'

        If (Trim(simulation_data%dft%vdw%type) == 'ddsc') Then
          Write (messages(4),'(1x,a)') '   3) the input reference data for non-interacting atoms is available only&
                                      & for elements of the first six rows of the periodic table except of lanthanide.'
          Call print_warnings(header, print_header, messages,4)
          If (Trim(simulation_data%dft%xc_version%type) /= 'pbe'    .And. &
             Trim(simulation_data%dft%xc_version%type) /= 'revpbe' .And. &
             Trim(simulation_data%dft%xc_version%type) /= 'blyp'   ) Then
             Write (message,'(1x,4a)')     '   4) WARNING: damping parameters&
                                      & are available only for PBE, revPBE and BLYP. For the requested "', & 
                                      & Trim(simulation_data%dft%xc_version%type), '" XC-functional, VDW_S6 (a0) and VDW_SR (b0)&
                                      & must be specified ', Trim(in_extra)
             Call print_warnings(header, print_header, message,1)
          End If
        Else If (Trim(simulation_data%dft%vdw%type) == 'ts'  .Or. &
                Trim(simulation_data%dft%vdw%type) == 'tsh' .Or. &
                Trim(simulation_data%dft%vdw%type) == 'mbd') Then
          Write (messages(4),'(1x,a)')   '   3) the method is NOT compatible with the setting ADDGRID=.TRUE.'
          If (Trim(simulation_data%dft%vdw%type) == 'ts' .Or. Trim(simulation_data%dft%vdw%type) == 'tsh') Then
             Write (messages(5),'(1x,a)')   '   4) modification of directives LVDW_EWALD, VDW_C6 and VDW_R0 is available&
                                      & from VASP 5.3.4'
             Write (messages(6),'(1x,a)')   '   5) the input reference data for non-interacting atoms is available only&
                                      & for elements of the first six rows of the periodic table except of lanthanide.'
             Call print_warnings(header, print_header, messages,6)

             If (Trim(simulation_data%dft%xc_version%type) /= 'pbe') Then     
               Write (message,'(1x,4a)')    '   6) WARNING: Defaults parameters controlling the damping function&
                                      & are available only for the PBE functional. For the requested "', & 
                                      & Trim(simulation_data%dft%xc_version%type), '" XC-functional, VDW_SR must be specified ',&
                                      & Trim(in_extra)
               Call print_warnings(header, print_header, message,1)
             End If
          Else If (Trim(simulation_data%dft%vdw%type) == 'mbd') Then
             Write (messages(5),'(1x,a)')   '   4) modification of directives VDW_C6 and VDW_R0 is available&
                                      & from VASP 5.3.4'
             Write (messages(6),'(1x,a)')   '   5) the method has sometimes numerical problems if highly polarizable atoms&
                                          & are located at short distances'
             Write (messages(7),'(1x,a)')   '   6) due to the long-range nature of dispersion interactions, the convergence&
                                          & of energy with respect to the number of k-points should be carefully examined'
             Call print_warnings(header, print_header, messages,7)
             If (Trim(simulation_data%dft%xc_version%type) /= 'pbe') Then     
               Write (message,'(1x,4a)')    '   7) WARNING: Defaults parameters controlling the damping function&
                                      & are available only for the PBE functional. For the requested "', & 
                                      & Trim(simulation_data%dft%xc_version%type), '" XC-functional, VDW_SR must be specified ',&
                                      & Trim(in_extra)
               Call print_warnings(header, print_header, message,1)
             End If
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
           Write (messages(2),'(1x,a)')   '   1) LDA (CA) pseudo potentials (PPs) could be used in principle,&
                                       & we have only allowed the use of PBE-PPs'
           Write (messages(3),'(1x,a)')   '   2) this vdW approximation is not defined for spin-polarised systems,&
                                        & but it is still possible to perform spin-polarised simulations'
           Call print_warnings(header, print_header, messages,3)
        End If 
      End If
    End If

  End Subroutine warnings_vasp_dft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Motion  
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine define_vasp_motion(files, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define (and check) motion settings for VASP directives (set defaults values)
    !
    ! author    - i.scivetti May 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256)  :: message, messages(16)
    Character(Len=256)  :: error_motion

    Integer(Kind=wi)    :: i
    Logical             :: error

    error_motion = '***ERROR in &motion_settings (file '//Trim(files(FILE_SET_EQCM)%filename)//'):'

    ! latest version of the code
    simulation_data%code_version= '5.4.1' 

    ! Ions related settings
    !!!!!!!!!!!!!!!!!!!!!!!
    !Relaxation method
    If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
      If (Trim(simulation_data%motion%relax_method%type) /= 'cg'    .And.&
        Trim(simulation_data%motion%relax_method%type) /= 'qn'  ) Then
        Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                &'Invalid specification of directive "relax_method" for VASP. Implemented options are:'
        Write (messages(2),'(1x,a)') '- CG (Conjugate Gradient)'
        Write (messages(3),'(1x,a)') '- QN (Quasi-Newton)'
        Call info(messages, 3)
        Call error_stop(' ')
      End If
    End If

    If (simulation_data%motion%change_cell_volume%stat .And. (.Not. simulation_data%motion%change_cell_shape%stat)) Then
       Write (messages(1),'(2(1x,a))') Trim(error_motion), 'Simulations that change the cell volume and keep the shape&
                                      & are not allowed in VASP. Please change option.'
       Call info(messages, 1)
       Call error_stop(' ')
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
      Write (messages(1),'(2a)')  Trim(error_motion), 'Invalid units of directive "force_tolerance" for VASP.&
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
                                    &'Invalid specification of "ensemble" for VASP. Options are:'
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
        If (Trim(simulation_data%motion%thermostat%type) /= 'andersen'  .And. &
           Trim(simulation_data%motion%thermostat%type) /= 'multi_andersen'  .And. &
           Trim(simulation_data%motion%thermostat%type) /= 'langevin'  .And. &
           Trim(simulation_data%motion%thermostat%type) /= 'nose-hoover'  ) Then
          Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                                  &'Specification for "thermostat" is not supported by VASP. Options are:'
          Write (messages(2),'(1x,a)') '- Andersen'
          Write (messages(3),'(1x,a)') '- Multi_Andersen'
          Write (messages(4),'(1x,a)') '- Langevin'
          Write (messages(5),'(1x,a)') '- Nose-Hoover'
          Call info(messages, 5)
          Call error_stop(' ')
        End If
      End If
    End If

    ! conditions
    If (simulation_data%motion%ensemble%fread) Then
      If (Trim(simulation_data%motion%ensemble%type) == 'npt') Then
        If (Trim(simulation_data%motion%thermostat%type)  == 'andersen'  .Or. &
           Trim(simulation_data%motion%thermostat%type)  == 'multi_andersen'  .Or. &
           Trim(simulation_data%motion%thermostat%type)  == 'nose-hoover'  ) Then
          Write (message,'(1x,6a)') Trim(error_motion), ' Ensemble "', Trim(simulation_data%motion%ensemble%type),&
                               & '" with thermostat "', Trim(simulation_data%motion%thermostat%type),&
                               & '" is not supported by VASP. Only Langevin is supported for NPT'
          Call error_stop(message)        
        End If
      End If
    End If
    
    ! Relax time thermostat
    If (simulation_data%motion%relax_time_thermostat%fread) Then
      If (Trim(simulation_data%motion%thermostat%type) /= 'langevin') Then
        Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                        &'In VASP, specification for "relax_time_thermostat" is only meaningful for the Langevin&
                        & thermostat.'
        Call info(messages, 1)       
        Call error_stop(' ')
      End If
    Else
      If (Trim(simulation_data%motion%thermostat%type) == 'langevin') Then 
        Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                        &'Missing specification "relax_time_thermostat" for the Langevin thermostat.'
        Write (messages(2),'(1x,a)')  '   Atomic friction values are set as the inverse of the&
                                    & "relax_time_thermostat" value in units of ps-1.' 
        Call info(messages, 2)       
        Call error_stop(' ')
 
      End If
    End If

    ! Relax time barostat
    If (Trim(simulation_data%motion%ensemble%type) == 'npt') Then
      If (.Not. simulation_data%motion%relax_time_barostat%fread) Then
        Write (messages(1),'(2(1x,a))') Trim(error_motion), &
                       &'In VASP, NPT simulations require the specification for "relax_time_barostat", which is missing.'
        Write (messages(2),'(1x,a)')  '   Friction value for the lattice  will be taken as the inverse of the&
                                    & "relax_time_barostat" value in units of ps-1.' 
        Call info(messages, 2)
        Call error_stop(' ')     
      End If
    End If
   
    If (simulation_data%extra_info%stat) Then
      ! Check if user defined directives contain only symbol "="
      Do i = 1, simulation_data%extra_directives%N0
        Call check_extra_directives(simulation_data%extra_directives%array(i), &
                                    simulation_data%extra_directives%key(i),   &
                                    simulation_data%extra_directives%set(i), '=', 'VASP')
      End Do
    End If
      
  End Subroutine define_vasp_motion
  
  Subroutine print_vasp_incar_motion(iunit, ic, net_elements, list_tag, simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print the motion part of the INCAR file 
    !
    ! author    - i.scivetti Sept  2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Integer(Kind=wi), Intent(InOut) :: ic 
    Integer(Kind=wi), Intent(In   ) :: net_elements 
    Character(Len=8), Intent(In   ) :: list_tag(max_components) 
    Type(simul_type), Intent(InOut) :: simulation_data

    Integer(Kind=wi)  ::  i, j
    Logical  ::  loop
    Character(Len=256) :: mass_list(max_components), mass_final
    Character(Len=256) :: message 

   Write (iunit,'(a)') ' '
   If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
     Write (iunit,'(a)') '#### Geometry relaxation'
     Write (iunit,'(a)') '#======================='
     ! atoms      
     Write (message,'(a,i4,a)') 'NSW = ', simulation_data%motion%ion_steps%value, ' # Number of ionic steps'
     Call record_directive(iunit, message, 'NSW', simulation_data%set_directives%array(ic), ic)
     If (Trim(simulation_data%motion%relax_method%type) == 'cg') Then
       Write (message, '(a)') 'IBRION = 2  # geometry relaxation with conjugate gradients'
       Call record_directive(iunit, message, 'IBRION', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%motion%relax_method%type) == 'qn') Then
       Write (message, '(a)') 'IBRION = 1  # geometry relaxation with the Quasi-Newton method'
       Call record_directive(iunit, message, 'IBRION', simulation_data%set_directives%array(ic), ic)
     End If
     Write (message,'(a,f6.2,a)') 'EDIFFG = ', -simulation_data%motion%delta_f%value(1), ' # Force tolerance eV/Angstrom' 
     Call record_directive(iunit, message, 'EDIFFG', simulation_data%set_directives%array(ic), ic)

     ! Cell
     If (simulation_data%motion%change_cell_volume%stat .And. simulation_data%motion%change_cell_shape%stat) Then
       Write (message,'(a)')      'ISIF =  3  # Also relax cell (change shape and volume)'  
       Call record_directive(iunit, message, 'ISIF', simulation_data%set_directives%array(ic), ic)
     ElseIf ((.Not. simulation_data%motion%change_cell_volume%stat) .And. simulation_data%motion%change_cell_shape%stat) Then
       Write (message,'(a)')      'ISIF =  4  # Also relax cell (change shape but keep volume fixed)'  
       Call record_directive(iunit, message, 'ISIF', simulation_data%set_directives%array(ic), ic)
     ElseIf ((.Not.simulation_data%motion%change_cell_volume%stat) .And. &
             (.Not. simulation_data%motion%change_cell_shape%stat)) Then
       Write (message,'(a)')      'ISIF =  2  # Keep cell fixed'  
       Call record_directive(iunit, message, 'ISIF', simulation_data%set_directives%array(ic), ic)
     End If


     If (simulation_data%motion%change_cell_volume%stat .Or. simulation_data%motion%change_cell_shape%stat) Then
       Write (message,'(a,f6.2,a)') 'PSTRESS = '  , simulation_data%motion%pressure%value, ' # External pressure (in kB)'
       Call record_directive(iunit, message, 'PSTRESS', simulation_data%set_directives%array(ic), ic)
     End If

   Else If (Trim(simulation_data%simulation%type) == 'md') Then
     Write (iunit,'(a)') '#### Molecular dynamics'
     Write (message, '(a)')       'IBRION = 0'
     Call record_directive(iunit, message, 'IBRION', simulation_data%set_directives%array(ic), ic)
     Write (message, '(a)')       'ISYM = 0'
     Call record_directive(iunit, message, 'ISYM', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,f6.2,a)') 'POTIM = '  , simulation_data%motion%timestep%value,    ' # Time step in fs'
     Call record_directive(iunit, message, 'POTIM', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,f6.2,a)') 'TEBEG = '  , simulation_data%motion%temperature%value, ' # Temperature in Kelvin'
     Call record_directive(iunit, message, 'TEBEG', simulation_data%set_directives%array(ic), ic)
     Write (message,'(a,i4,a)')   'NSW = ', simulation_data%motion%ion_steps%value, ' # Number of ionic steps'
     Call record_directive(iunit, message, 'NSW', simulation_data%set_directives%array(ic), ic)

     If (simulation_data%motion%mass_info%stat) Then
       Do i=1, net_elements
         j=1
         loop=.True.
         Do While (j <= simulation_data%total_tags .And. loop)
           If (Trim(list_tag(i))==Trim(simulation_data%motion%mass(j)%tag)) Then
             Write (mass_list(i),'(f8.3)') simulation_data%motion%mass(j)%value
             loop=.False.
           End If
           j=j+1
         End Do
       End Do
       Write (mass_final,'(*(1x,a10))') (Adjustl(Trim(mass_list(i))), i=1, net_elements)
       Write (iunit,'(a)') '# ==== Atomic masses' 
       Write (message,'(2a)') 'POMASS = ', Trim(mass_final)
       Call record_directive(iunit, message, 'POMASS', simulation_data%set_directives%array(ic), ic)
       Write (iunit,'(a)') ' '
     End If 

     If (Trim(simulation_data%motion%ensemble%type) == 'nve') Then
       Write (iunit, '(a)') '#==== NVE ensemble'
       Write (message, '(a)') 'MDALGO = 0'
       Call record_directive(iunit, message, 'MDALGO', simulation_data%set_directives%array(ic), ic)
       Write (message, '(a)') 'SMASS = -3'
       Call record_directive(iunit, message, 'SMASS', simulation_data%set_directives%array(ic), ic)
       Write (message, '(a)') 'ANDERSEN_PROB = 0.0'
       Call record_directive(iunit, message, 'ANDERSEN_PROB', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
       If (Trim(simulation_data%motion%thermostat%type) == 'andersen') Then
         Write (iunit, '(a)') '#==== NVT ensemble, Andersen thermostat'
         Write (message, '(a)') 'MDALGO = 1'
         Call record_directive(iunit, message, 'MDALGO', simulation_data%set_directives%array(ic), ic)
         Write (message, '(a)') 'ISIF   = 2'
         Call record_directive(iunit, message, 'ISIF', simulation_data%set_directives%array(ic), ic)
       Else If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover') Then
         Write (iunit, '(a)') '#==== NVT ensemble, Nose-Hoover thermostat'
         Write (message, '(a)') 'MDALGO = 2'
         Call record_directive(iunit, message, 'MDALGO', simulation_data%set_directives%array(ic), ic)
         Write (message, '(a)') 'ISIF   = 2'
         Call record_directive(iunit, message, 'ISIF', simulation_data%set_directives%array(ic), ic)
       Else If (Trim(simulation_data%motion%thermostat%type) == 'langevin') Then
         Write (iunit, '(a)') '#==== NVT ensemble, Langevin thermostat'
         Write (message, '(a)') 'MDALGO = 2'
         Call record_directive(iunit, message, 'MDALGO', simulation_data%set_directives%array(ic), ic)
         Write (message, '(a)') 'ISIF   = 2'
         Call record_directive(iunit, message, 'ISIF', simulation_data%set_directives%array(ic), ic)
         Write (message, '(a,*(f8.2))') 'LANGEVIN_GAMMA = ',&
                                    &( 1000.0_wp/simulation_data%motion%relax_time_thermostat%value, i=1,net_elements) 
         Call record_directive(iunit, message, 'LANGEVIN_GAMMA', simulation_data%set_directives%array(ic), ic)
       Else If (Trim(simulation_data%motion%thermostat%type) == 'multi_andersen') Then
         Write (iunit, '(a)') '#==== NVT ensemble, Multiple Andersen thermostat'
         Write (message, '(a)') 'MDALGO = 13'
         Call record_directive(iunit, message, 'MDALGO', simulation_data%set_directives%array(ic), ic)
         Write (message, '(a)') 'ISIF   = 2'
         Call record_directive(iunit, message, 'ISIF', simulation_data%set_directives%array(ic), ic)
       End If
     Else If (Trim(simulation_data%motion%ensemble%type) == 'npt') Then
       Write (iunit, '(a)') '#==== NPT ensemble, Langevin thermostat'
       Write (message,'(a)')         'MDALGO = 3'
       Call record_directive(iunit, message, 'MDALGO', simulation_data%set_directives%array(ic), ic)
       Write (message,'(a)')         'ISIF   = 3'
       Call record_directive(iunit, message, 'ISIF', simulation_data%set_directives%array(ic), ic)
       Write (message,'(a,f6.2,a)')  'PSTRESS = '  , simulation_data%motion%pressure%value, ' # External pressure (in kB)'
       Call record_directive(iunit, message, 'PSTRESS', simulation_data%set_directives%array(ic), ic)
       Write (message,'(a,*(f8.2))') 'LANGEVIN_GAMMA = ',&
                                    &( 1000.0_wp/simulation_data%motion%relax_time_thermostat%value, i=1,net_elements) 
       Call record_directive(iunit, message, 'LANGEVIN_GAMMA', simulation_data%set_directives%array(ic), ic)
       Write (message, '(a,f10.2)')   'LANGEVIN_GAMMA_L = ', 1000.0_wp/simulation_data%motion%relax_time_barostat%value
       Call record_directive(iunit, message, 'LANGEVIN_GAMMA_L', simulation_data%set_directives%array(ic), ic)
     Else If (Trim(simulation_data%motion%ensemble%type) == 'nph') Then
       Write (iunit, '(a)') '#==== NPH ensemble'
       Write (message, '(a)')       'MDALGO = 3'
       Call record_directive(iunit, message, 'MDALGO', simulation_data%set_directives%array(ic), ic)
       Write (message, '(a)')       'ISIF   = 3'
       Call record_directive(iunit, message, 'ISIF', simulation_data%set_directives%array(ic), ic)
       Write (message,'(a,f6.2,a)') 'PSTRESS = '  , simulation_data%motion%pressure%value, ' # External pressure (in kB)'
       Call record_directive(iunit, message, 'PSTRESS', simulation_data%set_directives%array(ic), ic)
       Write (message, '(a)')       'LANGEVIN_GAMMA_L = 0.0'
       Call record_directive(iunit, message, 'LANGEVIN_GAMMA_L', simulation_data%set_directives%array(ic), ic)
     End If
   End If

  End Subroutine print_vasp_incar_motion


  Subroutine advise_vasp_motion(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to instruct the user about motions settings 
    !
    ! author   - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(In   ) :: simulation_data

    Character(Len=256) :: messages(8), header
    Character(Len=256) :: in_extra
    Logical            :: print_header

    in_extra='using the &extra_directives block'
    
    print_header=.True.      
    
    ! MD-related parameters
    If (Trim(simulation_data%simulation%type) == 'md') Then
      Write (header, '(1x,a)')  'Regarding the MD convergence, the user should consider:'
      If (Trim(simulation_data%motion%ensemble%type) == 'nvt') Then
        If (Trim(simulation_data%motion%thermostat%type) == 'nose-hoover') Then
          Write (messages(1), '(1x,2a)')  ' - changing SMASS (must be > 0) ', Trim(in_extra)
          Call print_warnings(header, print_header, messages, 1)
        Else If (Trim(simulation_data%motion%thermostat%type) == 'langevin' ) Then
          Write (messages(1), '(1x,a)')  ' - optimizing the values for LANGEVIN_GAMMA using ALC_EQCM directive&
                                        & "relax_time_thermostat"'
          Call print_warnings(header, print_header, messages, 1)
        Else If (Trim(simulation_data%motion%thermostat%type) == 'andersen' .Or. &
                 Trim(simulation_data%motion%thermostat%type) == 'multi-andersen' ) Then
          Write (messages(1), '(1x,2a)')  ' - changing the value for ANDERSEN_PROB ', Trim(in_extra)
          Call print_warnings(header, print_header, messages, 1)
        End If
      Else If (Trim(simulation_data%motion%ensemble%type) == 'npt') Then
        Write (messages(1), '(1x,2a)')  ' - modifying the value for PMASS ', Trim(in_extra)
        Write (messages(2), '(1x,a)')   ' - optimising LANGEVIN_GAMMA using ALC_EQCM directive "relax_time_thermostat"'
        Write (messages(3), '(1x,2a)')  ' - optimising LANGEVIN_GAMMA_L using ALC_EQCM directive "relax_time_barostat"'
        Call print_warnings(header, print_header, messages, 3)
      Else If (Trim(simulation_data%motion%ensemble%type) == 'nph') Then
        Write (messages(1), '(1x,2a)')  ' - changing the value for PMASS ', Trim(in_extra)
        Write (messages(2), '(1x,2a)')  ' - optimising LANGEVIN_GAMMA_L using ALC_EQCM directive "relax_time_barostat"'
        Call print_warnings(header, print_header, messages, 2)
      End If
    End If

    ! Pulay
    If (simulation_data%motion%change_cell_volume%stat) Then
      If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
        Write (header, '(1x,a)')  'Regarding the relaxation of the simulation cell, the user should consider:'
      End If
      Write (messages(1), '(1x,a)')  ' - increasing the ENCUT via ALC_EQCM directive "enery_cutoff" to minimise the&
                                      & Pulay stress from changing the cell volume'
      Call print_warnings(header, print_header, messages, 1)
      If (Trim(simulation_data%simulation%type) == 'relax_geometry') Then
         Write (messages(1), '(1x,a)')  ' - setting directive "precision" to "Accurate", if possible,&
                                        & for accurate cell relaxation'
         Call print_warnings(header, print_header, messages, 1)
      End If 
    End If

  End Subroutine advise_vasp_motion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Pseudopotentials  
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine check_pseudo_potentials_vasp(simulation_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check PPs for VASP
    !
    ! author    - i.scivetti Sept 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(simul_type),   Intent(InOut) :: simulation_data

    Character(Len=256)  :: message
    Character(Len=256)  :: exec_grep, path
    Character(Len=256)  :: word, word2, word3
    Character(Len=256)  :: xc, pp_path

    Integer(Kind=wi)    :: i, k, io, internal, ipos

    pp_path   = Trim(FOLDER_DFT)//'/PPs/'    

    ! Check consistency between pseudpotentials and XC directive
    Do i=1, simulation_data%total_tags
      exec_grep='grep "LEXC" '//Trim(pp_path)//Trim(simulation_data%dft%pseudo_pot(i)%file_name)//' > xc.dat'
      Call execute_command_line(exec_grep)
      Open(Newunit=internal, File='xc.dat' ,Status='old')
      Read (internal, Fmt=*, iostat=io) word, word, word
      If (Trim(word)=='PE') Then
        xc='pbe'
      Else If (Trim(word)=='91') Then
        xc='pw91'
      Else If (Trim(word)=='CA') Then
        xc='ca'
      End If
      If (Trim(xc) /= Trim(simulation_data%dft%xc_base)) Then
        Write (message,'(1x,5a)') '***ERROR in ',  Trim(pp_path),' folder: File ',&
                           & Trim(simulation_data%dft%pseudo_pot(i)%file_name), &
                           & ' corresponds to a exchange-correlation type of functional that is different from the option&
                           & specified in directive "XC_version". Please check if this PP is consistent with VASP.'
        Call error_stop(message)
      End If
      Close(internal)
      Call execute_command_line('rm xc.dat')
    End Do

    ! Check consistency between pseudpotentials and elements
    Do i=1, simulation_data%total_tags
      path=Trim(pp_path)//Trim(simulation_data%dft%pseudo_pot(i)%file_name)
      Open(Newunit=internal, File=path ,Status='old')
      Read (internal, Fmt=*, iostat=io) word, word2, word3 
      If (io /= 0 ) Then
        Write (message,'(1x,10a)') '***ERROR in ', Trim(pp_path), ' folder: File ', &
                                   Trim(simulation_data%dft%pseudo_pot(i)%file_name), ' appears to have the wrong headings.&
                                 & Please check.'
        Call error_stop(message)
      End If
      ipos= Index(word2,'_') 
      If (ipos /= 0) Then
        word2(ipos:256)=' ' 
      End If 
      If (Trim(word2) /= Trim(simulation_data%dft%pseudo_pot(i)%element)) Then
        Write (message,'(1x,10a)') '***ERROR in ', Trim(pp_path), ' folder: File ', &
                           & Trim(simulation_data%dft%pseudo_pot(i)%file_name), &
                           & ' corresponds to element "', Trim(word2), '", which does not agree with with the element&
                           & associated with tag "', Trim(simulation_data%dft%pseudo_pot(i)%tag), '" as set in&
                           & "&pseudo_potentials". Please check if this PP is consistent with VASP.'
        Call error_stop(message)
      End If
      Do k=1, simulation_data%total_tags
        If (Trim(simulation_data%component(k)%tag)==Trim(simulation_data%dft%pseudo_pot(i)%tag)) Then 
          Read (internal, Fmt=*, iostat=io) simulation_data%component(k)%valence_electrons
        End If
      End Do
      Close(internal)
    End Do
  
  End Subroutine check_pseudo_potentials_vasp
  
End Module code_vasp
