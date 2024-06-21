!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module with subroutines tools for simulation
!
! Copyright: 2022-2024 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
!
! Author        - i.scivetti April-Oct 2021
! Refactoring   - i.scivetti Oct       2021
! Refactoring   - i.scivetti Sep       2022
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module simulation_tools

  Use constants,        Only : max_components

  Use numprec,          Only : wi, &
                               wp

  Use process_data,     Only : capital_to_lower_case, &
                               remove_symbols  
                             
  Use references,       Only : bib_ca, bib_hl, bib_pz, bib_wigner, &
                               bib_vwn, bib_pade, bib_pw92, bib_pw91, bib_am05, bib_pbe, bib_rp, bib_revpbe, &
                               bib_pbesol, bib_blyp, bib_wc, bib_xlyp, bib_scan, bib_rpw86pbe, &
                               bib_g06, bib_obs, bib_jchs, bib_dftd2, bib_dftd3, bib_dftd3bj, bib_ts, bib_tsh, bib_mbd,&
                               bib_ddsc, bib_vdwdf, bib_optpbe, bib_optb88, bib_optb86b, bib_vdwdf2, bib_vdwdf2b86r,&
                               bib_SCANrVV10, bib_VV10, bib_AVV10S, bib_tunega, bib_rpw86, bib_fisher, web_d3bJ, &
                               web_d3bJ, web_vasp, web_onetep, web_castep, web_cp2k
                             
  Use simulation_setup, Only:  type_extra, &
                               type_ref_data                             
  Use unit_output,      Only : error_stop,&
                               info 

  ! Error in the initialization of magnetization                            
  Real(Kind=wp), Parameter, Public :: error_mag = 0.00001_wp
                               
  Public ::  obtain_xc_reference, obtain_vdw_reference
  Public ::  check_extra_directives, check_settings_set_extra_directives, extra_directive_setting
  Public ::  check_settings_single_extra_directive
  Public ::  scan_extra_directive, print_extra_directives
  Public ::  print_warnings
  Public ::  set_reference_database
  Public ::  record_directive 
  Public ::  check_initial_magnetization

Contains

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

  Subroutine check_extra_directives(sentence, key, set, symbol, code)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check the structure of sub-block &extra_directives 
    !
    ! author        - i.scivetti Jul 2021
    ! modification  - i.scivetti Aug  2022 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=256), Intent(In   ) :: sentence
    Character(Len=256), Intent(  Out) :: key      
    Character(Len=256), Intent(  Out) :: set      
    Character(Len=*),   Intent(In   ) :: symbol 
    Character(Len=*),   Intent(In   ) :: code

    Character(Len=256)  :: word 
    Character(Len=256)  :: messages(4)
    Integer(Kind=wi)    :: io  
    Logical             :: error

    error=.False.
    key=' '
    set=' ' 

    Write (messages(1),'(1x,a)') '***ERROR: in block &extra_directives: user-defined directive'
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
   
    If (Index(Trim(Adjustl(sentence)), '#') == 0) Then
      Read(sentence, Fmt='(a)') word      
      Call remove_symbols(word,Trim(symbol))
      Read(word, Fmt=*, iostat=io) key, set 
      If (is_iostat_end(io)) Then
        Call info(' ', 1)       
        Write (messages(1), '(1x,3a)') '***ERROR in sub-bock &extra_directives: problems in the definiton of directive "',&
                                & Trim(key),'". Please check format, syntax and any missing information.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
      Call capital_to_lower_case(key)
      Call capital_to_lower_case(set)    
    End If

  End Subroutine check_extra_directives  
 
  Subroutine check_settings_set_extra_directives(ref_data, num_ref_data, extra_directives, exception_keys, &
                                               & number_exceptions, extradir_header)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if reference extra directives (ref_data) are defined or not within the
    ! &extra_directives block. If defined, check if the definition is correct.
    ! If there are not defined, print a message advising the user to consider for 
    ! the directive. Number_exceptions is the number of keywords that will not 
    ! be checked, and these keyword are defined in exception_keys
    !
    ! author    - i.scivetti August 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(type_ref_data), Intent(In   ) :: ref_data(:)
    Integer(Kind=wi),    Intent(In   ) :: num_ref_data
    Type(type_extra),    Intent(In   ) :: extra_directives
    Character(Len=*),    Intent(In   ) :: exception_keys(:)
    Integer(Kind=wi),    Intent(In   ) :: number_exceptions
    Logical,             Intent(InOut) :: extradir_header

    Logical             :: found, print_msn
    Character(Len=256)  :: set, msn
    Character(Len=256)  :: messages(50)
    Integer(Kind=wi)    :: j, k, indx, iex

    indx=0

    Do k= 1, num_ref_data
      found=.False.
      iex=1 

      Do While (iex <= number_exceptions .And. (.Not. found)) 
        If (Trim(ref_data(k)%key) == Trim(exception_keys(iex))) Then
          found=.True.
        End If
        iex=iex+1
      End Do

      If (.Not. found) Then

        j=1
        Do While (j <= extra_directives%N0 .And. (.Not. found))
          If (Trim(ref_data(k)%key) == Trim(extra_directives%key(j))) Then
            found=.True.
            set=extra_directives%set(j)
          End If
          j=j+1
        End Do

        Call extra_directive_setting(ref_data(k), set, found, print_msn, msn)
        If (print_msn) Then
           indx=indx+1
           messages(indx)=msn     
        End If  

      End If
    End Do
    
    If (indx > 0) Then
      If (.Not. extradir_header) Then
        Call info(' The following extra keywords can be defined in the &extra_directives sub-block:', 1)
        extradir_header=.True.
      End If
      Call info(messages, indx)
    End If  

  End Subroutine check_settings_set_extra_directives

  Subroutine extra_directive_setting(ref_data, set, found, print_msn, msn, condition)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check what to print or not for a given directive, with name ref_data%key  
    !
    ! author    - i.scivetti August 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(type_ref_data), Intent(In   ) :: ref_data
    Character(Len=*),    Intent(In   ) :: set
    Logical,             Intent(In   ) :: found
    Logical,             Intent(  Out) :: print_msn
    Character(Len=*),    Intent(  Out) :: msn
    Character(Len=*),  Optional,  Intent(In   ) :: condition

    Integer(Kind=wi)    :: extra_int, io
    Logical             :: extra_logic    
    Real(Kind=wp)       :: extra_real
    Character(Len=256)  :: action_dir

    print_msn=.False.

    If (Present(condition))then
      action_dir='complain'
    Else
      action_dir='allow'
    End If    
 
    If (found) Then
      If (Trim(ref_data%keytype) == 'integer') Then 
        Read(set, Fmt=*, iostat=io) extra_int
      Else If (Trim(ref_data%keytype) == 'logical') Then
        Read(set, Fmt=*, iostat=io) extra_logic
      Else If (Trim(ref_data%keytype) == 'real') Then
        Read(set, Fmt=*, iostat=io) extra_real
      End If
      If (io/=0) Then
        print_msn=.True.
        Write (msn, '(1x,5a)') '  *** PROBLEMS *** keyword "', Trim(ref_data%key),&
                                  & '" (defined in &extra_direcitves) must be ', Trim(ref_data%keytype),&
                                  & '. PLEASE FIX THIS KEYWORD. To keep the already generated models,&
                                  & change the "analysis" directive to "hpc_simulation_files" and rerun.'
      End If
      If (Trim(action_dir) == 'complain') Then
        print_msn=.True.
        Write (msn, '(1x,5a)') '  *** PROBLEMS *** keyword "', Trim(ref_data%key), '" (defined in&
                                  & &extra_direcitves) is not compatible with the option of "', Trim(condition),&
                                  & '". PLEASE FIX. To keep the already generated models,&
                                  & use the "hpc_simulation_files" option for "analysis" and rerun.'
      End If
    Else
      If (Trim(action_dir) == 'allow') Then
        print_msn=.True.
        Write (msn, '(1x,5a,1x,a)') '  * ', Trim(ref_data%key), Trim(ref_data%msn),&
                               & ' Default is ', Trim(ref_data%set_default), Trim(ref_data%units)
      End If
    End If

  End Subroutine extra_directive_setting
  

  Subroutine check_settings_single_extra_directive(single_dir, ref_data, num_ref_data, extra_directives,&
                                                 & extradir_header, condition)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The same as check_settings_set_extra_directives but for a single selected
    ! keyword, given by single_dir
    !
    ! author    - i.scivetti August 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*),    Intent(In   ) :: single_dir
    Type(type_ref_data), Intent(In   ) :: ref_data(:)
    Integer(Kind=wi),    Intent(In   ) :: num_ref_data
    Type(type_extra),    Intent(In   ) :: extra_directives
    Logical,             Intent(InOut) :: extradir_header
    Character(Len=*), Optional,  Intent(In   ) :: condition

    Logical             :: found, print_msn
    Type(type_ref_data) :: single_ref_data
    Character(Len=256)  :: set
    Character(Len=256)  :: message
    Integer(Kind=wi)    :: j

    found=.False.
    j=1
    Do While (j <= num_ref_data .And. (.Not. found))
      If (Trim(single_dir) == Trim(ref_data(j)%key)) Then
        found=.True.
        single_ref_data=ref_data(j)
      End If
      j=j+1
    End Do

    found=.False.
    j=1
    Do While (j <= extra_directives%N0 .And. (.Not. found))
      If (Trim(single_dir) == Trim(extra_directives%key(j))) Then
        found=.True.
        set=extra_directives%set(j)
      End If
      j=j+1
    End Do

    If (Present(condition)) Then
      Call extra_directive_setting(single_ref_data, set, found, print_msn, message, condition)
    Else
      Call extra_directive_setting(single_ref_data, set, found, print_msn, message)
    End If        
   
    If (print_msn) Then
      If (.Not. extradir_header) Then
        Call info(' The following extra keywords can be defined in the &extra_directives sub-block:', 1)
        extradir_header=.True.
      End If
      Call info(message, 1)
    End If  

  End Subroutine check_settings_single_extra_directive

  Subroutine print_extra_directives(iunit, extra_directives, set_directives, code, simul_type)
    Integer(Kind=wi),    Intent(In   ) :: iunit
    Type(type_extra),    Intent(In   ) :: extra_directives
    Type(type_extra),    Intent(InOut) :: set_directives          
    Character(Len=*),    Intent(In   ) :: code
    Character(Len=*),    Intent(In   ) :: simul_type
  
    Character(Len=256)  :: messages(10)
    Integer(Kind=wi)    :: i
    Logical :: found
    
    found=.False.
  
     Write (iunit,'(a)') ' '
     Write (iunit,'(a)') '##### Extra directives'
     Write (iunit,'(a)') '#====================='
     Do i=1, extra_directives%N0
       Write (iunit,'(a)') Trim(Adjustl(extra_directives%array(i)))
       If (Index(Trim(Adjustl(extra_directives%array(i))), '#') /= 1 ) Then
         If (Trim(code)=='onetep') Then
           If (Trim(extra_directives%key(i)) =='&block'     .And. &
              Trim(extra_directives%set(i)) =='thermostat') Then
             If (Trim(simul_type) == 'md') Then
               Write (messages(1), '(1x,a)') '***ERROR in sub-block &extra_directives: "&block thermostat" CANNOT be defined as&
                                      & part of &extra_directives. Please remove it.'
               Call info(messages, 1)                         
               Call error_stop(' ')                         
             End If  
           End If
         End If
         Call scan_extra_directive(extra_directives%key(i), set_directives, found)
         If (found)Then
           Call info(' ', 1)
           Write (messages(1), '(1x,a)')  '*************************************************************************************'
           Write (messages(2), '(1x,a)')  '** WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING **'
           Write (messages(3), '(1x,a)')  '*************************************************************************************'
           Write (messages(4), '(1x,3a)') '  PROBLEMS in sub-bock &extra_directives: directive "', &
                                          & Trim(extra_directives%key(i)), '" has already'
           Write (messages(5), '(1x, a)') '  been set from the definition of ALC_EQCM directives. The user must review the'
           Write (messages(6), '(1x, a)') '  settings of &extra_directives. THERE MUST NOT BE DUPLICATION OF DIRECTIVES.' 
           Write (messages(7), '(1x,3a)') '  To keep the already generated models, remove "', &
                                         Trim(extra_directives%key(i)), '" from the &extra_directives' 
           Write (messages(8), '(1x, a)') '  sub-block, change the "analysis" directive to "hpc_simulation_files" and rerun.'
           Write (messages(9), '(1x, a)') '*************************************************************************************'
           Write (messages(10), '(1x,a)') '*************************************************************************************'
           Call info(messages, 10)
           Call info(' ', 1)
         End If
         set_directives%N0=set_directives%N0+1
         set_directives%array(set_directives%N0)=Trim(extra_directives%key(i))
       End If
     End Do
     
  End Subroutine print_extra_directives 
  
  Subroutine scan_extra_directive(sentence, set_directives, found)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to scan the directives defined in block &extra_directives against
    ! the directives set from the definitions of &block_simulation_settings 
    ! 
    ! author        - i.scivetti July 2021
    ! modification  - i.scivetti Aug  2022 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: sentence 
    Type(type_extra), Intent(In   ) :: set_directives
    Logical,          Intent(  Out) :: found
   
    Character(Len=256) :: word
    Integer(Kind=wi)   :: j

    found=.False.

    j=1
    Do While (j <= set_directives%N0 .And. (.Not. found))
      word=set_directives%array(j)
      Call capital_to_lower_case(word)
      If (Trim(sentence) == Trim(word)) Then 
        found=.True.
      End If
      j=j+1
    End Do

  End Subroutine scan_extra_directive

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
  
  Subroutine set_reference_database(ref_database, num_ref_data, code, functionality)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set reference database for keywords compatible with the requested
    ! code and functionality. Info read in &extra_directives will be
    ! compared against the following settings
    !
    ! author    - i.scivetti August 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(type_ref_data), Intent(  Out) :: ref_database(:)
    Integer(Kind=wi),    Intent(  Out) :: num_ref_data
    Character(Len=*),    Intent(In   ) :: code
    Character(Len=*),    Intent(In   ) :: functionality

    If (Trim(code) == 'ONETEP') Then
      If (Trim(functionality) == 'EDFT') Then

        ref_database(1)%key='edft_maxit'
        ref_database(1)%keytype='integer'
        ref_database(1)%msn=': maximum number of inner loop iterations.'
        ref_database(1)%set_default='10'
        ref_database(1)%units=' '

        ref_database(2)%key='edft_write_occ'
        ref_database(2)%keytype='logical'
        ref_database(2)%msn=': writes the occupancies and the energy levels to file with .occ extension.'
        ref_database(2)%set_default='False'
        ref_database(2)%units=' '

        ref_database(3)%key='edft_trial_step'
        ref_database(3)%keytype='real'
        ref_database(3)%msn=': sets the value of lambda, fixing the step size of the inner loop&
                           & and switching off the line search for optimum. If negative the normal line search is used.'
        ref_database(3)%set_default='-1'
        ref_database(3)%units=' '

        ref_database(4)%key='edft_free_energy_thres'
        ref_database(4)%keytype='real'
        ref_database(4)%msn=': maximum difference in the Helmholtz free energy functional per atom between two&
                               & consecutive iterations.'
        ref_database(4)%set_default='1.0E-6'
        ref_database(4)%units='Hartree'

        ref_database(5)%key='edft_energy_thres'
        ref_database(5)%keytype='real'
        ref_database(5)%msn=': maximum difference in the energy functional per atom between two consecutive iterations.'
        ref_database(5)%set_default='1.0E-6'
        ref_database(5)%units='Hartree'

        ref_database(6)%key='edft_entropy_thres'
        ref_database(6)%keytype='real'
        ref_database(6)%msn=': maximum difference in the entropy per atom between two consecutive iterations.'
        ref_database(6)%set_default='1.0E-6'
        ref_database(6)%units='Hartree'

        ref_database(7)%key='edft_rms_gradient_thres'
        ref_database(7)%keytype='real'
        ref_database(7)%msn=': maximum RMSgradient.'
        ref_database(7)%set_default='1.0E-4'
        ref_database(7)%units='Hartree'

        ref_database(8)%key='edft_commutator_thres'
        ref_database(8)%keytype='real'
        ref_database(8)%msn=': maximum value of the Hamiltonian-Kernel commutator.'
        ref_database(8)%set_default='1.0E-5'
        ref_database(8)%units='Hartree'
        
        ref_database(9)%key='edft_fermi_thres'
        ref_database(9)%keytype='real'
        ref_database(9)%msn=': maximum change in the Fermi energy between two consecutive iterations.'
        ref_database(9)%set_default='1.0E-3'
        ref_database(9)%units='Hartree'
        
        ref_database(10)%key='edft_round_evals'
        ref_database(10)%keytype='integer'
        ref_database(10)%msn=': when set to n>0, the occupancies that result from the Fermi-Dirac distribution&
                                  & are rounded to n significant figures.'
        ref_database(10)%set_default='-1'
        ref_database(10)%units=' '

        ref_database(11)%key='edft_max_step'
        ref_database(11)%keytype='real'
        ref_database(11)%msn=': maximum step during line search.'
        ref_database(11)%set_default='1.0'
        ref_database(11)%units=' '

        ref_database(12)%key='ngwf_cg_rotate'
        ref_database(12)%keytype='logical'
        ref_database(12)%msn=': rotation of eigenvectors to the new NGWF representation once these are updated.'
        ref_database(12)%set_default='False'
        ref_database(12)%units=' '

        ref_database(13)%key='edft_ham_diis_size'
        ref_database(13)%keytype='integer'
        ref_database(13)%msn=': maximum number of Hamiltonians used from previous iterations to&
                               & generate the new guess through Pulay mixing.'
        ref_database(13)%set_default='10'
        ref_database(13)%units=' '

        ref_database(14)%key='write_hamiltonian'
        ref_database(14)%keytype='logical'
        ref_database(14)%msn=': write the Hamiltonian matrix on a .ham file.'
        ref_database(14)%set_default='False'
        ref_database(14)%units=' '

        ref_database(15)%key='read_hamiltonian'
        ref_database(15)%keytype='logical'
        ref_database(15)%msn=': read the Hamiltonian matrix from a .ham file.'
        ref_database(15)%set_default='False'
        ref_database(15)%units=' '

        num_ref_data=15
      End If

     If (Trim(functionality) == 'DL_MG') Then

       ref_database(1)%key='mg_use_cg'
       ref_database(1)%keytype='logical'
       ref_database(1)%msn=': if .True., it turns on the conjugate gradient solver to increase&
                           & stability at expenses of reducing the performance.'
       ref_database(1)%set_default='False'
       ref_database(1)%units=' '

       ref_database(2)%key='mg_max_iters_vcycle'
       ref_database(2)%keytype='integer'
       ref_database(2)%msn=': maximum number of multigrid V-cycle iterations.'
       ref_database(2)%set_default='50'
       ref_database(2)%units=' '

       ref_database(3)%key='mg_vcyc_smoother_iter_pre'
       ref_database(3)%keytype='integer'
       ref_database(3)%msn=': sets the number of V-cycle smoother iterations pre-smoothing.&
                           & Set it to 4 and 8 for difficult systems.'
       ref_database(3)%set_default='2'
       ref_database(3)%units=' '

       ref_database(4)%key='mg_vcyc_smoother_iter_post'
       ref_database(4)%keytype='integer'
       ref_database(4)%msn=': sets the number of V-cycle smoother iterations post-smoothing.&
                          & Set it to 4 and 8 for difficult systems.'
       ref_database(4)%set_default='1'
       ref_database(4)%units=' '

       ref_database(5)%key='mg_defco_fd_order'
       ref_database(5)%keytype='integer'
       ref_database(5)%msn=': sets the discretization order used when solving the P(B)E.&
                                & Value must be even an integer. Recommended value is 8.'
       ref_database(5)%set_default='8'
       ref_database(5)%units=' '

       ref_database(6)%key='mg_max_iters_newton'
       ref_database(6)%keytype='integer'
       ref_database(6)%msn=': sets the maximum number of Newton method iterations. Only relevant if&
                            & solver is set to "full" (&electrolyte). Increase it for difficult systems.'
       ref_database(6)%set_default='30'
       ref_database(6)%units=' '
     
        num_ref_data=6
      End If

    End If

  End Subroutine set_reference_database
  
  
End Module simulation_tools
