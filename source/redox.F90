!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for redox related variables and subroutines
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! author: i.scivetti August 2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module redox 

  Use electrode,    Only : electrode_type 
  Use eqcm,         Only : eqcm_type
  Use fileset,      Only : file_type, &
                           FILE_CHARACT, &
                           FOLDER_ANALYSIS 
  Use numprec,      Only : wi,wp 
  Use unit_output,  Only : error_stop,&
                           info

  Implicit None
  Private

  ! Redox variables for Electrochemical Characterization
  Type, Public :: redox_type
    Private

    Integer(Kind=wi), Public     :: ncycles
    Integer(Kind=wi), Public     :: maxpoints
    Integer(Kind=wi), Public     :: limit_cycles
    Integer(Kind=wi), Public     :: legs
    Real(Kind=wp),    Public     :: min_voltage
    Real(Kind=wp),    Public     :: max_voltage
    Logical,          Public     :: l_voltage_range

    ! Charge variables
    Real(Kind=wp), Allocatable, Public :: DQ_ox(:)
    Real(Kind=wp), Allocatable, Public :: Q_ox(:,:)
    Real(Kind=wp), Allocatable, Public :: DQ_red(:)
    Real(Kind=wp), Allocatable, Public :: Q_red(:,:)

    ! Mass variables
    Real(Kind=wp), Allocatable, Public :: mass(:,:)
    Real(Kind=wp), Allocatable, Public :: M_ox(:,:)
    Real(Kind=wp), Allocatable, Public :: M_red(:,:)
    Real(Kind=wp), Allocatable, Public :: DM_ox(:)
    Real(Kind=wp), Allocatable, Public :: DM_red(:)
    Real(Kind=wp), Allocatable, Public :: DM_half(:)
    Real(Kind=wp), Allocatable, Public :: DM_residual(:)

    ! Difference variables
    Real(Kind=wp), Allocatable, Public :: Q_diff(:)

    ! Other Variables
    Integer(Kind=wi), Allocatable,  Public :: points(:)
    Real(Kind=wp), Allocatable,     Public :: time(:,:)
    Real(Kind=wp), Allocatable,     Public :: voltage(:,:)
    Real(Kind=wp), Allocatable,     Public :: current(:,:)
    Character(Len=16), Allocatable, Public :: label_leg(:,:)

    ! Limits
    Integer(Kind=wi), Allocatable, Public  :: voltage_index(:,:,:)

  Contains
    Private
    Procedure :: init_variables => allocate_redox_variables
    Procedure :: init_mass      => allocate_redox_mass_arrays
    Procedure :: init_charge    => allocate_redox_charge_arrays
    Final     :: cleanup
  End Type redox_type


  Public :: redox_characterization 

Contains

  Subroutine allocate_redox_variables(T, flag_time)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocation of essential redox arrays 
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(redox_type), Intent(InOut)  :: T
    Logical          , Intent(In   )  :: flag_time 

    Character(Len=256)  :: message
    Integer(Kind=wi) :: fail(4)

    Allocate(T%voltage(T%maxpoints,T%ncycles) , Stat=fail(1))
    Allocate(T%current(T%maxpoints,T%ncycles)   , Stat=fail(2))
    Allocate(T%points(T%ncycles)                , Stat=fail(3))
    Allocate(T%label_leg(2,T%ncycles)           , Stat=fail(4))

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for general arrays&
                                & (subroutine allocate_redox_variables)'
      Call error_stop(message)
    End If

    !Initialise
    T%voltage=0.0_wp
    T%current=0.0_wp
    T%points=0
    T%label_leg= ' '

    If (flag_time) Then
      Allocate(T%time(T%maxpoints,T%ncycles) , Stat=fail(1))
      If (fail(1) > 0) Then
        Write (message,'(1x,1a)') '***ERROR: Allocation problems for elapsed time&
                                  & (subroutine allocate_redox_variables)'
        Call error_stop(message)
      End If
      T%time=0.0_wp
    End If

  End Subroutine allocate_redox_variables

  Subroutine allocate_redox_mass_arrays(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocation of redox arrays for mass changes. Variables allocated depend 
    ! on the data read from experiments
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(redox_type), Intent(InOut)  :: T

    Character(Len=256)  :: message
    Integer(Kind=wi) :: fail(7)

    Allocate(T%DM_ox(T%ncycles),             Stat=fail(1))
    Allocate(T%DM_red(T%ncycles),            Stat=fail(2))
    Allocate(T%DM_residual(T%ncycles),       Stat=fail(3))
    Allocate(T%DM_half(T%ncycles),           Stat=fail(4))
    Allocate(T%mass(T%maxpoints,T%ncycles),  Stat=fail(5))
    Allocate(T%M_ox(T%maxpoints,T%ncycles),  Stat=fail(6))
    Allocate(T%M_red(T%maxpoints,T%ncycles), Stat=fail(7))

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of redox mass related arrays&
                                & (subroutine allocate_redox_mass_arrays)'
      Call error_stop(message)
    End If

    !Initialise
    T%mass=0.0_wp
    T%DM_ox=0.0_wp
    T%DM_red=0.0_wp
    T%DM_residual=0.0_wp
    T%DM_half=0.0_wp

  End Subroutine allocate_redox_mass_arrays

  Subroutine allocate_redox_charge_arrays(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocation of redox arrays for charge changes. Variables allocated 
    ! depend on the data read from experiments
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(redox_type), Intent(InOut)  :: T

    Character(Len=256) :: message
    Integer(Kind=wi)   :: fail(5)

    Allocate(T%DQ_ox(T%ncycles)         , Stat=fail(1))
    Allocate(T%DQ_red(T%ncycles)        , Stat=fail(2))
    Allocate(T%Q_diff(T%ncycles)        , Stat=fail(3))
    Allocate(T%Q_ox(T%maxpoints,T%ncycles)  , Stat=fail(4))
    Allocate(T%Q_red(T%maxpoints,T%ncycles) , Stat=fail(5))

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems of redox charge arrays&
                                & (subroutine allocate_redox_mass_arrays)'
      Call error_stop(message)
    End If

    !Initialise
    T%DQ_ox=0.0_wp
    T%DQ_red=0.0_wp
    T%Q_diff=0.0_wp
    T%Q_ox=0.0_wp
    T%Q_red=0.0_wp

  End Subroutine allocate_redox_charge_arrays

  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocation for the final process
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(redox_type) :: T

    If (Allocated(T%DQ_ox)) Then
      Deallocate(T%DQ_ox)
    End If

    If (Allocated(T%DQ_red)) Then
      Deallocate(T%DQ_red)
    End If

    If (Allocated(T%Q_ox)) Then
      Deallocate(T%Q_ox)
    End If

    If (Allocated(T%Q_red)) Then
      Deallocate(T%Q_red)
    End If

    If (Allocated(T%DM_ox)) Then
      Deallocate(T%DM_ox)
    End If

    If (Allocated(T%DM_red)) Then
      Deallocate(T%DM_red)
    End If

    If (Allocated(T%mass)) Then
      Deallocate(T%mass)
    End If

    If (Allocated(T%DM_residual)) Then
      Deallocate(T%DM_residual)
    End If

    If (Allocated(T%DM_half)) Then
      Deallocate(T%DM_half)
    End If

  End Subroutine cleanup
 
  Subroutine redox_characterization(files, redox_data, eqcm_data, electrode_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to quantify change of mass and total charge during CV cycling
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut) :: files(:)
    Type(redox_type),     Intent(InOut) :: redox_data
    Type(eqcm_type),      Intent(InOut) :: eqcm_data
    Type(electrode_type), Intent(InOut) :: electrode_data 

    Character(Len=256)   :: messages(3)
    
    redox_data%l_voltage_range=.False.

    Call info(' ',1)
    Write (messages(1),'(1x,a)') 'Characterization analysis'
    Write (messages(2),'(1x,a)') '========================='
    Call info(messages, 2)

    Call redox_cycles(eqcm_data, redox_data)

    !Extract values and integrate quantities
    If (eqcm_data%current%fread) Then
      Call integrate_redox_current(eqcm_data, redox_data)
    End If   

    !Extract values and integrate quantities
    If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
      Call extract_Dmass(eqcm_data, redox_data, electrode_data)
    End If

    Write (messages(1),'(1x,a)') 'Relevant computed quantities from the reported data'
    Call info(messages, 1)
    Call integrated_quantitites_summary(redox_data, eqcm_data) 
    Call print_redox_characterization(files, redox_data, eqcm_data)

    If (eqcm_data%voltage_range%fread) Then 
      ! Determine if the requested voltage range includes of not the whole CV domain
      Call check_voltage_range_domain(eqcm_data, redox_data)
      Call info(' ',1)

      If (redox_data%l_voltage_range) Then
        If (eqcm_data%current%fread) Then
          Call integrate_redox_current(eqcm_data, redox_data)
        End If   
        Write (messages(1),'(1x,a,2(f8.2,a))') 'Computed quantities within the requested voltage range: ', &
                                                eqcm_data%voltage_range%value(1),  ' V /', &
                                                eqcm_data%voltage_range%value(2),  ' V'
        Call info(messages, 1)
        If (eqcm_data%current%fread) Then
          Call integrate_redox_current(eqcm_data, redox_data)
        End If   
        If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
          Call extract_Dmass(eqcm_data, redox_data, electrode_data)
        End If

        Call integrated_quantitites_summary(redox_data, eqcm_data) 
      Else
        Write (messages(1),'(1x,a)') '**** WARNING: The set of recorded EQCM is within the requested voltage range.'
        Write (messages(2),'(1x,a)') '     There is NO need to recompute the quantities of the table above.'
        Call info(messages, 2)
      End If  
    End If  

  End Subroutine redox_characterization

  Subroutine check_voltage_range_domain(eqcm_data, redox_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to assign the eqcm values of CV cycles to redox cycles
    ! 
    ! author    - i.scivetti July 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),    Intent(InOut) :: eqcm_data
    Type(redox_type),   Intent(InOut) :: redox_data

    Integer(Kind=wi)   :: i, k, ic
    Real(Kind=wp)      :: average, DV
    Character(Len=256) :: message

    redox_data%min_voltage= Huge(1.0_wp)
    redox_data%max_voltage=-Huge(1.0_wp)

    ic=0
    average=0.0_wp

    DV=Abs(eqcm_data%voltage_range%value(2)-eqcm_data%voltage_range%value(1))

    Do i = 1, redox_data%limit_cycles
      Do k=1, redox_data%points(i)
         If (k < redox_data%points(i))Then
           ic=ic+1
           average=average+Abs(redox_data%voltage(k+1,i)-redox_data%voltage(k,i))           
         End If        
         If(redox_data%min_voltage > redox_data%voltage(k,i)) Then
           redox_data%min_voltage= redox_data%voltage(k,i)      
         End If        
         If(redox_data%max_voltage < redox_data%voltage(k,i)) Then
           redox_data%max_voltage= redox_data%voltage(k,i)      
         End If        
      End Do 
    End Do

    average=average/ic
    If (average > DV) Then     
       Write(message, '(a)') '****ERROR: the voltage range defined with the "voltage_range" directive is smaller than&
                            & the average discretiation voltage used for the EQCM experiment! Please review the settings'
       Call error_stop(message)               
    End If 

    If ((eqcm_data%voltage_range%value(1) > redox_data%min_voltage) .Or. &
        (eqcm_data%voltage_range%value(2) < redox_data%max_voltage)) Then
        redox_data%l_voltage_range=.True.
    End If    

  End Subroutine check_voltage_range_domain

  Subroutine redox_cycles(eqcm_data, redox_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define the number of redox cycles and the number of point per each 
    ! redox cycle
    ! 
    ! author    - i.scivetti July 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),    Intent(In   ) :: eqcm_data
    Type(redox_type),   Intent(InOut) :: redox_data

    Integer(Kind=wi)  :: i, j, l, k, m
    Integer(Kind=wi)  :: icycle, n_tot, n_ini
    Logical           :: flag, flag_neg, flag_pos, flag_ini
    Integer(Kind=wi)  :: points(eqcm_data%ncycles)

    Integer(Kind=wi), Allocatable :: i1(:)
    Integer(Kind=wi), Allocatable :: i2(:)
    Integer(Kind=wi), Allocatable :: i3(:)

    ! Obtain the total number of points
    n_tot=0
    Do i = 1, eqcm_data%ncycles
      Do j= 1, 2
        Do k=1, eqcm_data%segment_points(j,i)
          n_tot=n_tot+1 
        End Do
      End Do
    End Do 

    ! Allocate temporary arrays 
    Allocate(i1(n_tot)) 
    Allocate(i2(n_tot)) 
    Allocate(i3(n_tot)) 

    ! define index arrays
    m=0
    Do i = 1, eqcm_data%ncycles
      Do j= 1, 2
        Do k=1, eqcm_data%segment_points(j,i)
          m=m+1 
          i1(m)=k
          i2(m)=j
          i3(m)=i
        End Do
      End Do
    End Do 

    !Initialise variables
    l=0 
    icycle=1
    flag =.True.
    flag_neg=.False.
    flag_pos=.False.
    flag_ini=.True.
    points=0

    ! Set oxidation and reduction parts of the redox cycles
    Do m=1, n_tot 
      If (eqcm_data%scan_sweep_ref=='negative' .And. eqcm_data%current%value(i1(m),i2(m),i3(m)) < 0.0_wp) Then
        flag_neg=.True.
        If (flag_ini) Then
          flag_ini=.False.
          n_ini=m      
        End If       
      End If          

      If (flag_neg) Then
        If (i3(m)/=1 .And. eqcm_data%label_leg(i2(m),i3(m))=='negative' .And. &
            eqcm_data%current%value(i1(m),i2(m),i3(m))<0.0_wp .And. flag) Then
          points(icycle)=l
          icycle=icycle+1
          l=0
          flag=.False.
        End If
        l=l+1
        If (eqcm_data%label_leg(i2(m),i3(m))=='positive' .And. (.Not. flag)) Then
           flag=.True.
        End If   
      End If
  
      If (eqcm_data%scan_sweep_ref=='positive' .And. eqcm_data%current%value(i1(m),i2(m),i3(m)) >= 0.0_wp) Then
        flag_pos=.True. 
        If (flag_ini) Then
          flag_ini=.False.
          n_ini=m      
        End If       
      End If          

      If (flag_pos) Then
        If (i3(m)/=1 .And. eqcm_data%label_leg(i2(m),i3(m))=='positive' .And. &
          eqcm_data%current%value(i1(m),i2(m),i3(m)) > 0.0_wp .And. flag) Then
          points(icycle)=l
          icycle=icycle+1
          l=0
          flag=.False.
        End If
        l=l+1
        If (eqcm_data%label_leg(i2(m),i3(m))=='negative' .And. (.Not. flag)) Then
           flag=.True.
        End If
      End If
    End Do 
    points(icycle)=l

    redox_data%ncycles   = icycle
    redox_data%maxpoints = maxval(points)

    ! Initialise redox related arrays
    Call redox_data%init_variables(eqcm_data%time%fread)
    ! label
    Call label_redox_processes(eqcm_data, redox_data)

    ! Allocate charge variables
    If (eqcm_data%current%fread) Then
      Call redox_data%init_charge()
    End If

    ! Allocate mass variables
    If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
      Call redox_data%init_mass()
    End If
 
    ! Define the amount of cycles
    Do i=1, redox_data%ncycles
      redox_data%points(i)=points(i)
    End Do

    l=0 
    icycle=1

    ! Assign EQCM data to redox 
    Do m=n_ini, n_tot
      If(l+1 > redox_data%points(icycle))Then
        l=1
        icycle=icycle+1
      Else
        l=l+1
      End If
      
      redox_data%voltage(l,icycle) = eqcm_data%voltage%value(i1(m),i2(m),i3(m))
      redox_data%current(l,icycle) = eqcm_data%current%value(i1(m),i2(m),i3(m))

      If (eqcm_data%time%fread) Then
        redox_data%time(l,icycle)=eqcm_data%time%value(i1(m),i2(m),i3(m))
      End If

      If (eqcm_data%mass%fread) Then
          redox_data%mass(l,icycle)=eqcm_data%mass%value(i1(m),i2(m),i3(m))
      Else 
        If (eqcm_data%mass_frequency%fread) Then
          redox_data%mass(l,icycle)=eqcm_data%mass_frequency%value(i1(m),i2(m),i3(m))
        End If
      End If
    End Do 
    
    ! Deallocate
    Deallocate(i3) 
    Deallocate(i2) 
    Deallocate(i1) 

  End Subroutine redox_cycles   

  Subroutine integrated_quantitites_summary(redox_data, eqcm_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Summarise the information found for integrated charge and/or mass for
    ! each cycle
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(redox_type),     Intent(InOut) :: redox_data
    Type(eqcm_type),      Intent(In   ) :: eqcm_data

    Integer(Kind=wi)   :: i, j
    Character(Len=256) :: message, line
    Character(Len= 32) :: fmt1
    Character(Len= 32) :: fmt2
    Character(Len= 32) :: fmt3
    Character(Len= 32) :: fmt4
    Character(Len= 32) :: fmt5

    Logical                      :: lmlch, lch       

    lch=(.Not. eqcm_data%mass_frequency%fread) .And.  (eqcm_data%current%fread)
    lmlch=(eqcm_data%mass_frequency%fread) .And.  (eqcm_data%current%fread)

    If (lmlch) Then
      fmt1 = '(1x,4(a,5x))'
      fmt2 = '(1x,i3, 5x,a,2(9x,f12.3))'
      fmt3 = '(1x,    8x,a,2(9x,f12.3))'
      fmt4 = '(1x,i3, 5x,a)'
      fmt5 = '(1x,    8x,a)'
    ElseIf (lch) Then
      fmt1 = '(1x,3(a,5x))'
      fmt2 = '(1x,i3, 6x,a,9x,f12.3)'
      fmt3 = '(1x,    9x,a,9x,f12.3)'
      fmt4 = '(1x,i3, a)'
      fmt5 = '(1x,    a)'
    End If

    If (lmlch) Write (line,'(1x,a)') '-----------------------------------------------------------------'
    If (lch)   Write (line,'(1x,a)') '--------------------------------------------'
    Call info(line, 1)
    If (lmlch) Write (message, fmt1) 'Cycle', 'Process', 'Total charge [mC]', 'Mass change  [ng]'
    If (lch)   Write (message, fmt1) 'Cycle', 'Process', 'Total charge [mC]'
    Call info(message, 1)
    If (lmlch) Write (line,'(1x,a)') '-----------------------------------------------------------------'
    If (lch)   Write (line,'(1x,a)') '--------------------------------------------'
    Call info(line, 1)

    Do i = eqcm_data%range_cycles%value(1), redox_data%limit_cycles
      Do j= 1, 2
        If (j==1) Then
          If (redox_data%label_leg(j,i) /= 'undefined') Then
            If (redox_data%label_leg(j,i)=='oxidation') Then
              If (lmlch) Write (message, fmt2) i, Trim(redox_data%label_leg(j,i)), redox_data%DQ_ox(i), redox_data%DM_ox(i)
              If (lch)   Write (message, fmt2) i, Trim(redox_data%label_leg(j,i)), redox_data%DQ_ox(i)
            Else If (redox_data%label_leg(j,i)=='reduction') Then
              If (lmlch) Write (message, fmt2) i, Trim(redox_data%label_leg(j,i)), redox_data%DQ_red(i), redox_data%DM_red(i)
              If (lch)   Write (message, fmt2) i, Trim(redox_data%label_leg(j,i)), redox_data%DQ_red(i)
            End If
          Else
            Write (message, fmt4) i, ' ' 
          End If
        Else
          If (redox_data%label_leg(j,i) /= 'undefined') Then
            If (redox_data%label_leg(j,i)=='oxidation') Then
              If (lmlch) Write (message, fmt3) Trim(redox_data%label_leg(j,i)), redox_data%DQ_ox(i), redox_data%DM_ox(i)
              If (lch)   Write (message, fmt3) Trim(redox_data%label_leg(j,i)), redox_data%DQ_ox(i)
            Else If (redox_data%label_leg(j,i)=='reduction') Then
              If (lmlch) Write (message, fmt3) Trim(redox_data%label_leg(j,i)), redox_data%DQ_red(i), redox_data%DM_red(i)
              If (lch)   Write (message, fmt3) Trim(redox_data%label_leg(j,i)), redox_data%DQ_red(i)
            End If
          Else
            Write (message, fmt5) ' ' 
          End If
        End If
         Call info(message, 1)
      End Do
    End Do

    Call info(line, 1)

  End Subroutine integrated_quantitites_summary   

  Subroutine print_redox_characterization(files, redox_data, eqcm_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print computed information to file CHARACTERIZATION for characterization 
    ! of the redox cycles
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(redox_type),  Intent(InOut) :: redox_data
    Type(eqcm_type),   Intent(In   ) :: eqcm_data
 
    Integer(Kind=wi)  :: iunit, i
    Real(Kind=wp)     :: dMdQ_red, dMdQ_ox, dQdQ 
    Character(Len=64) :: lchg_diff, lchg_ratio

    Call info(' ', 1)
    Call info(' Print results to file '//Trim(FOLDER_ANALYSIS)//'/'//Trim(files(FILE_CHARACT)%filename), 1)

    Open(Newunit=files(FILE_CHARACT)%unit_no, File=files(FILE_CHARACT)%filename,Status='Replace')
    iunit=files(FILE_CHARACT)%unit_no

    If (eqcm_data%current%fread) Then
      If (eqcm_data%scan_sweep_ref=='negative') Then
        lchg_diff ='|Q_red+Q_ox| (mC)'
        lchg_ratio='|Q_red/Q_ox|'
      ElseIf (eqcm_data%scan_sweep_ref=='positive') Then
        lchg_diff ='|Q_ox+Q_red| (mC)'
        lchg_ratio='|Q_ox/Q_red|'
      End If
    End If

    If (eqcm_data%mass_frequency%fread .And. (eqcm_data%current%fread)) Then
      Write (iunit,'(1x,11(a,4x))') '#Cycle', 'Dmass_oxidation [ng]', 'Dmass_reduction [ng]', &
                                       'Dmass_residual [ng]', 'Dmass_half_redox [ng]', &
                                       'Q_oxidation [mC]', 'Q_reduction [mC]',&
                                        Trim(lchg_diff), Trim(lchg_ratio), '|Dmass/Q|_red [ng/mC]', &
                                       '|Dmass/Q|_ox [ng/mC]' 
      Do i=eqcm_data%range_cycles%value(1), redox_data%limit_cycles
        If (Abs(redox_data%DQ_ox(i))<epsilon(redox_data%DQ_ox(i))) Then
          dMdQ_ox=0.0_wp
        Else
          dMdQ_ox=Abs(redox_data%DM_ox(i)/redox_data%DQ_ox(i))
        End If
        If (Abs(redox_data%DQ_red(i))<epsilon(redox_data%DQ_red(i))) Then
          dMdQ_red=0.0_wp
        Else
          dMdQ_red=Abs(redox_data%DM_red(i)/redox_data%DQ_red(i))
        End If

        If (Abs(redox_data%DQ_ox(i))<epsilon(redox_data%DQ_ox(i)) .Or. &
           Abs(redox_data%DQ_red(i))<epsilon(redox_data%DQ_red(i)) ) Then
            dQdQ=0.0_wp 
        Else
          If (lchg_ratio=='|Q_ox/Q_red|') Then
            dQdQ=Abs(redox_data%DQ_ox(i)/redox_data%DQ_red(i))
          ElseIf (lchg_ratio=='|Q_red/Q_ox|') Then
            dQdQ=Abs(redox_data%DQ_red(i)/redox_data%DQ_ox(i))
          End If
        End If       

        Write (iunit,'(i3,10(10x,f12.3))')  i, redox_data%DM_ox(i), redox_data%DM_red(i), &
                                           redox_data%DM_residual(i), redox_data%DM_half(i),&
                                           redox_data%DQ_ox(i), redox_data%DQ_red(i), &
                                           redox_data%Q_diff(i), dQdQ, dMdQ_red, dMdQ_ox
      End Do

    Else If (.Not. eqcm_data%mass_frequency%fread .And. eqcm_data%current%fread) Then
      Write (iunit,'(1x,a,8x,4(a,10x))') '#Cycle','Q_oxidation [mC]', 'Q_reduction [mC]',&
                                        Trim(lchg_diff), Trim(lchg_ratio)
      Do i= eqcm_data%range_cycles%value(1), redox_data%limit_cycles

        If (Abs(redox_data%DQ_ox(i))<epsilon(redox_data%DQ_ox(i)) .Or. &
           Abs(redox_data%DQ_red(i))<epsilon(redox_data%DQ_red(i)) ) Then
            dQdQ=0.0_wp 
        Else
          If (lchg_ratio=='|Q_ox/Q_red|') Then
            dQdQ=Abs(redox_data%DQ_ox(i)/redox_data%DQ_red(i))
          ElseIf (lchg_ratio=='|Q_red/Q_ox|') Then
            dQdQ=Abs(redox_data%DQ_red(i)/redox_data%DQ_ox(i))
          End If
        End If       

        Write (iunit,'(i3,4(13x,f12.3))')  i, redox_data%DQ_ox(i), redox_data%DQ_red(i), &
                                          redox_data%Q_diff(i), dQdQ
      End Do

    End If  

    ! Close file
    Close(iunit)
    ! mv files
    Call execute_command_line('mv '//Trim(files(FILE_CHARACT)%filename)//' '//Trim(FOLDER_ANALYSIS))

  End Subroutine print_redox_characterization  

  Subroutine label_redox_processes(eqcm_data, redox_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Identify redox fragment is reduction, oxidation or undefined
    ! 
    ! author    - i.scivetti July 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),    Intent(In   ) :: eqcm_data
    Type(redox_type),   Intent(InOut) :: redox_data

    Logical            :: flag
    Integer(Kind=wi)   :: i, j
    Character(Len=256) :: message
    Character(Len=256) :: messages(2)

    If (redox_data%ncycles < eqcm_data%range_cycles%value(2)) Then
      redox_data%limit_cycles=redox_data%ncycles
      Write (message,'(1x,a,i3)') '***IMPORTANT: Identified REDOX cycles from the input CV data (file DATA_EQCM): ', &
                                 redox_data%ncycles
      Call info(message,1)
      Call info(' ', 1)
    Else
      redox_data%limit_cycles=eqcm_data%range_cycles%value(2)
    End If
   
    flag=.True.

    Do i = 1, redox_data%ncycles
      Do j= 1, 2
        If (eqcm_data%label_leg(j,i)=='negative') Then
          redox_data%label_leg(j,i)= 'reduction'
        Else If (eqcm_data%label_leg(j,i)=='positive') Then
          redox_data%label_leg(j,i)= 'oxidation'
        Else If (eqcm_data%label_leg(j,i)=='undefined') Then
          redox_data%label_leg(j,i)= 'undefined'
        End If
      End Do
      If (Any(redox_data%label_leg(:,i)=='undefined')) Then
        If (flag) Then
          If (redox_data%limit_cycles > (i-2)) Then
            Write (messages(1),'(1x,2(a,i3),a)') '***ERROR: Cycle ', i, ' is not a proper CV cycle, and the upper limit of ',&
                                            redox_data%limit_cycles,' cycles is not suitable for characterization analysis.'
            If (i-2<=0) Then
              Write (message,'(1x,a,i3,a)') 'Fix problems with cycle ', i, ', or eliminate incomplete cycles from DATA_EQCM'
            Else
              Write (message,'(1x,(a,i3))') 'The only possible characterization with the given CV data is to set&
                                           & the upper limit of directive cycle equal or lower than ', i-2 
            End If    
            Call info(messages,1)
            Call error_stop(message)
          End If
          flag=.False.
        End If
      End If

    End Do

  End Subroutine label_redox_processes   

  Subroutine integrate_redox_current(eqcm_data, redox_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Integrate the current within the voltage/time range for oxidation/reduction
    ! The operation is computed within the voltage range 
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),   Intent(InOut) :: eqcm_data
    Type(redox_type),  Intent(InOut) :: redox_data

    Integer(Kind=wi) :: i, k
    Real(Kind=wp)    :: Dchg, Dt

    Do i = 1, redox_data%limit_cycles 
      redox_data%DQ_ox(i)=0.0_wp
      redox_data%DQ_red(i)=0.0_wp
      Dchg=0.0_wp
      If (redox_data%points(i) >2) Then
        Do k=2, redox_data%points(i)
          If (eqcm_data%time%fread) Then
            Dt=Abs(redox_data%time(k,i)-redox_data%time(k-1,i))
          Else
            Dt=Abs(redox_data%voltage(k,i)-redox_data%voltage(k-1,i))/eqcm_data%scan_rate%value(1)
          End If 
          Dchg=Dt*(redox_data%current(k,i)+redox_data%current(k-1,i))/2.0_wp

          If (redox_data%current(k,i) >= 0.0_wp) Then
            If (.Not. redox_data%l_voltage_range) Then      
              redox_data%DQ_ox(i)=redox_data%DQ_ox(i)+Dchg
              redox_data%Q_ox(k,i)=redox_data%Q_ox(k-1,i)+Dchg
            Else
              If ((eqcm_data%voltage_range%value(1) <= redox_data%voltage(k,i)) .And. &
                 (eqcm_data%voltage_range%value(2) >= redox_data%voltage(k,i))) Then   
                  redox_data%DQ_ox(i)=redox_data%DQ_ox(i)+Dchg
                  redox_data%Q_ox(k,i)=redox_data%Q_ox(k-1,i)+Dchg
              End If        
            End If        
          Else
            If (.Not. redox_data%l_voltage_range) Then      
              redox_data%DQ_red(i)=redox_data%DQ_red(i)+Dchg
              redox_data%Q_red(k,i)=redox_data%Q_red(k-1,i)+Dchg
            Else
              If ((eqcm_data%voltage_range%value(1) <= redox_data%voltage(k,i)) .And. &
                 (eqcm_data%voltage_range%value(2) >= redox_data%voltage(k,i))) Then   
                  redox_data%DQ_red(i)=redox_data%DQ_red(i)+Dchg
                  redox_data%Q_red(k,i)=redox_data%Q_red(k-1,i)+Dchg
              End If
            End If
          End If
        End Do
      End If 
      redox_data%Q_diff(i)=Abs(redox_data%DQ_ox(i)+redox_data%DQ_red(i))
    End Do    

  End Subroutine integrate_redox_current

  Subroutine extract_Dmass(eqcm_data, redox_data, electrode_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Obtain mass variation for oxidation/reduction. 
    ! Mass densities are multiplied by the electrode area
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(redox_type),     Intent(InOut) :: redox_data
    Type(eqcm_type),      Intent(InOut) :: eqcm_data
    Type(electrode_type), Intent(In   ) :: electrode_data 

    Integer(Kind=wi)   :: i, k
    Real(Kind=wp)      :: dV, dVfw, mass_limit
    Logical            :: fpass
    Character(Len=256) :: message

    Real(Kind=wp)      :: mox_f(2)
    Real(Kind=wp)      :: mox_b(2)
    Real(Kind=wp)      :: mred_f(2)
    Real(Kind=wp)      :: mred_b(2)

    Logical            :: fox_f(2)
    Logical            :: fox_b(2)
    Logical            :: fred_f(2)
    Logical            :: fred_b(2)

    If (.Not. redox_data%l_voltage_range) Then
      ! Only if mass is not read from DATA_EQCM
      If (.Not. eqcm_data%mass%fread) Then 
        redox_data%mass=electrode_data%area_geom%value*redox_data%mass  
      End If

      Do i = 1, redox_data%limit_cycles
        mass_limit=0.0_wp
        k=1
        fpass=.True.

        Do While (k <= redox_data%points(i)-1 .And. fpass)
          dV=redox_data%voltage(k+1,i)-redox_data%voltage(k,i)
          If (eqcm_data%scan_sweep_ref=='negative' .And. dV>0 ) Then
            If (redox_data%current(k,i)>0) Then
              Mass_limit = redox_data%mass(k-1,i)
              fpass=.False.
            End If
          End If
   
          If (eqcm_data%scan_sweep_ref=='positive' .And. dV<0 ) Then
            If (redox_data%current(k,i)<0) Then
              Mass_limit = redox_data%mass(k-1,i)
              fpass=.False.
            End If
          End If
          k=k+1
        End Do

        If (Abs(mass_limit)<epsilon(mass_limit)) Then
          Write (message,'(1x,a,i3)') 'Problems to determine the change of mass for cycle ', i
          Call error_stop(message)   
        End If 
  
        If (eqcm_data%scan_sweep_ref=='negative') Then
          redox_data%DM_red(i)=      mass_limit-redox_data%mass(1,i)  
          redox_data%DM_ox(i) =      redox_data%mass(redox_data%points(i),i)-mass_limit  
          redox_data%DM_residual(i)= redox_data%DM_red(i)+redox_data%DM_ox(i)
          If (i== 1) Then
            redox_data%DM_half(i)=redox_data%DM_red(i) 
          Else
            redox_data%DM_half(i)=redox_data%DM_red(i)+redox_data%DM_residual(i-1)
            redox_data%DM_residual(i)=redox_data%DM_residual(i)+redox_data%DM_residual(i-1) 
          End If
        Else If (eqcm_data%scan_sweep_ref=='positive') Then
          redox_data%DM_red(i)= redox_data%mass(redox_data%points(i),i)-mass_limit 
          redox_data%DM_ox(i) = mass_limit-redox_data%mass(1,i)
          redox_data%DM_residual(i)=redox_data%DM_ox(i)+redox_data%DM_red(i)
          If (i== 1) Then
            redox_data%DM_half(i)=redox_data%DM_ox(i) 
          Else
            redox_data%DM_half(i)=redox_data%DM_ox(i)+redox_data%DM_residual(i-1)
            redox_data%DM_residual(i)=redox_data%DM_residual(i)+redox_data%DM_residual(i-1) 
          End If
        End If
      End Do
  
    Else   
      Do i = 1, redox_data%limit_cycles
        k=1
        ! Flags to determine voltage domains
        fpass=.True.
        fox_f(:)=.True.
        fox_b(:)=.True.
        fred_f(:)=.True.
        fred_b(:)=.True.

        ! Mass values 
        mred_f(:)=0.0_wp
        mred_b(:)=0.0_wp
        mox_f(:)=0.0_wp
        mox_b(:)=0.0_wp


        Do While (k <= redox_data%points(i)-1 .And. fpass)
          dV=redox_data%voltage(k+1,i)-redox_data%voltage(k,i)
          If (k<redox_data%points(i)-1) Then
            dVfw=redox_data%voltage(k+2,i)-redox_data%voltage(k+1,i)
          End If
          If (redox_data%current(k,i) >= 0.0_wp) Then
            If (dV>0) Then
              If ((eqcm_data%voltage_range%value(1) <= redox_data%voltage(k,i)) .And. fox_f(1)) Then
                mox_f(1)=redox_data%mass(k,i)
                fox_f(1)=.False.      
              End If
              If ((eqcm_data%voltage_range%value(2) < redox_data%voltage(k,i)) .And. fox_f(2)) Then
                mox_f(2)=redox_data%mass(k-1,i)
                fox_f(2)=.False.      
              End If
              If ((dVfw<0.0_wp .Or. (k+1 == redox_data%points(i))) .And. fox_f(2)) Then
                mox_f(2)=redox_data%mass(k+1,i)
                fox_f(2)=.False.      
              End If         
              ! Red
              If (fred_f(1) .And. (.Not. fred_b(1))) Then
                 mred_f(1)=redox_data%mass(k,i)
                 fred_f(1)=.False.  
               End If
               If (fred_f(2) .And. (.Not. fred_b(2))) Then
                 mred_f(2)=redox_data%mass(k,i)
                 fred_f(2)=.False.  
              End If
            Else If (dV<0) Then 
              If ((eqcm_data%voltage_range%value(2) >= redox_data%voltage(k,i)) .And. fox_b(2) .And. (.Not. fox_f(2))) Then
                mox_b(2)=redox_data%mass(k,i)
                fox_b(2)=.False.  
              End If
              If ((eqcm_data%voltage_range%value(1) > redox_data%voltage(k,i)) .And. fox_b(1) .And. (.Not. fox_f(1))) Then
                mox_b(1)=redox_data%mass(k-1,i)
                fox_b(1)=.False.  
              End If

              If (k+1==redox_data%points(i)) Then
                If (fox_b(2)) Then
                  mox_b(2)=redox_data%mass(k,i)
                End If
                If (fox_b(1)) Then
                  mox_b(1)=redox_data%mass(k,i)
                End If
              End If

            End If
          Else
            If (dV <0) Then
              If (fox_b(1) .And. (.Not. fox_f(1))) Then
                mox_b(1)=redox_data%mass(k,i)
                fox_b(1)=.False.  
              End If
              If (fox_b(2) .And. (.Not. fox_f(2))) Then
                mox_b(2)=redox_data%mass(k,i)
                fox_b(2)=.False.  
              End If

              ! Red
              If ((eqcm_data%voltage_range%value(2) >= redox_data%voltage(k,i)) .And. fred_b(1)) Then
                mred_b(1)=redox_data%mass(k,i)
                fred_b(1)=.False.      
              End If
              If ((eqcm_data%voltage_range%value(1) > redox_data%voltage(k,i)) .And. fred_b(2)) Then
                mred_b(2)=redox_data%mass(k-1,i)
                fred_b(2)=.False.      
              End If
              If ((dVfw>0.0_wp .Or. (k+1 == redox_data%points(i))) .And. fred_b(2)) Then
                mred_b(2)=redox_data%mass(k+1,i)
                fred_b(2)=.False.      
              End If         

            Else If (dV>0) Then
              If ((eqcm_data%voltage_range%value(1) <= redox_data%voltage(k,i)) .And. fred_f(1) .And. (.Not. fred_b(1))) Then
                mred_f(1)=redox_data%mass(k,i)
                fred_f(1)=.False.  
              End If
              If ((eqcm_data%voltage_range%value(2) < redox_data%voltage(k,i)) .And. fred_f(2) .And. (.Not. fred_b(2))) Then
                mred_f(2)=redox_data%mass(k-1,i)
                fred_f(2)=.False.  
              End If
              If (k+1==redox_data%points(i)) Then
                If (fred_f(2)) Then
                  mred_f(2)=redox_data%mass(k,i)
                End If
                If (fred_f(1)) Then
                  mred_f(1)=redox_data%mass(k,i)
                End If        
              End If
            End If
          End If
          k=k+1
        End Do

        redox_data%DM_ox(i)= (mox_f(2)-mox_f(1))+(mox_b(2)-mox_b(1))
        redox_data%DM_red(i)= (mred_f(2)-mred_f(1))+(mred_b(2)-mred_b(1))

      End Do

    End If

  End Subroutine extract_Dmass

End Module redox
