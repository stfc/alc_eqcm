!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for:
!  - reading of the DATA_EQCM file
!  - definition of eqcm variables
!  - data analysis
!  - mass calibration
!  - massograms
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:        i.scivetti April    2020
! Contribution:  i.scivetti October  2020 (mass calibration)
!                i.scivetti November 2020 (massogram)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module eqcm

  Use fileset,      Only : file_type, &
                           FILE_CALIBRATION, &
                           FILE_DATA_EQCM, &
                           FILE_FILTERED_CURRENT, &
                           FILE_FILTERED_MASS, &
                           FILE_MASSOGRAM, &
                           FILE_RAW_CURRENT, &
                           FILE_RAW_MASS, &
                           FILE_SET_EQCM, &
                           FILE_SPEC_MASS,&
                           FILE_SPEC_CURRENT,&
                           FOLDER_ANALYSIS, & 
                           refresh_out_eqcm

  Use filtering,    Only : filter_type, &
                           filter_data  

  Use fft,          Only : fft_type,   &
                           fft_points, &
                           fft_setup,  &
                           fft_spectrum, &
                           fftnmax, &
                           print_spectrum   
 
  Use input_types,  Only : in_integer_array, &
                           in_logic,   &
                           in_scalar,  &
                           in_param,   & 
                           in_param_array, & 
                           in_string

  Use numprec,      Only : wi,&
                           wp 
  
  Use process_data, Only : capital_to_lower_case, &
                           detect_rubbish,        &
                           remove_symbols,        &
                           remove_front_tabs
  Use unit_output,  Only : error_stop, &
                           info

  Implicit None
  Private

  ! Maximum number of EQCM variables available
  Integer(Kind=wi) :: max_num_var = 9

  Type, Private :: exp_data
    Real(Kind=wp),    Allocatable  :: value(:,:,:)
    Real(Kind=wp)     :: convert
    Character(Len=16) :: units= repeat(' ',16)
    Integer(Kind=wi)  :: col = 0
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Logical           :: warn = .False. 
  End Type

  ! Type for eqcm data and analysis
  Type, Public :: eqcm_type
    Private
    ! Type of electrochemical process: electrodeposition or intercalation
    Type(in_string), Public       :: process

    ! Efficiency of the intercalation
    Type(in_scalar), Public       ::  efficiency

    !Offset to EQCM current?
    Type(in_logic),    Public     :: current_offset

    ! EQCM type software
    Type(in_string),   Public     :: software

    ! EQCM type analysis
    Type(in_string),   Public     :: analysis

    ! Conversion factors for units
    Real(Kind=wp), Public         :: fact_current 

    ! Maximum of data points among all segments of all cycles
    Integer(Kind=wi), Public      :: max_points

    ! Measured points per segment of CV
    Integer(Kind=wi), Allocatable, Public  :: segment_points(:,:)

    ! Label the Voltage sweeps as positive and negative 
    Character(Len=16), Allocatable, Public :: label_leg(:,:)

    ! Reference scan sweep
    Character(Len=16), Public              :: scan_sweep_ref

    ! Input EQCM data
    Real(Kind=wp), Allocatable, Public     :: whole_set(:,:,:,:)

    ! EQCM arrays of interest (for now)
    Type(exp_data),   Public  :: current
    Type(exp_data),   Public  :: charge
    Type(exp_data),   Public  :: time
    Type(exp_data),   Public  :: voltage 
    Type(exp_data),   Public  :: mass_frequency
    Type(exp_data),   Public  :: mass 
    Type(exp_data),   Public  :: resistance

    ! Sensitivity factor (Volts to frequency)
    Type(in_scalar),  Public  :: V_to_Hz  

    ! Range of EQMC cycles be indentified from EQCM data
    Integer(Kind=wi), Public  :: ncycles 

    ! Number of EQMC cycles to be analysed (from SET_EQCM) 
    Type(in_integer_array), Public :: range_cycles

    ! Voltage range 
    Type(in_param_array),   Public :: voltage_range

    ! Scan rate
    Type(in_param_array),   Public :: scan_rate

    ! Number of colums of EQMC files
    Integer(Kind=wi),       Public :: columns

  Contains
      Private
      Procedure         :: init_arrays => allocate_eqcm_data_arrays
      Procedure, Public :: init_input_variables => allocate_eqcm_input_variables
      Final             :: cleanup
  End Type eqcm_type


  Public :: read_eqcm_data, eqcm_spectra, eqcm_filter, print_eqcm_data
  Public :: eqcm_check_voltage_range, eqcm_cycling_summary
  Public :: eqcm_mass_calibration, eqcm_massogram, get_eqcm_charge

Contains

  Subroutine allocate_eqcm_input_variables(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate essential eqcm input variables to be read from SET_EQCM 
    ! 
    ! author    - i.scivetti Nov 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(eqcm_type), Intent(InOut)  :: T

    Integer(Kind=wi)   :: fail(5)
    Character(Len=256) :: message

    Allocate(T%range_cycles%value(2),  Stat=fail(1)) 
    Allocate(T%voltage_range%value(2), Stat=fail(2))
    Allocate(T%voltage_range%units(1), Stat=fail(3))   
    Allocate(T%scan_rate%value(1),     Stat=fail(4))
    Allocate(T%scan_rate%units(2),     Stat=fail(5))   

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: allocation problems for essential input&
                                & EQCM variables (subroutine allocate_eqcm_input_variables)'
      Call error_stop(message)
    End If

    T%range_cycles%value(:) = 0
    T%voltage_range%value(1)= -huge(0.0_wp)
    T%voltage_range%value(2)=  huge(0.0_wp)

  End Subroutine allocate_eqcm_input_variables

  Subroutine allocate_eqcm_data_arrays(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate eqcm arrays depending on what it was read in DATA_EQCM 
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(eqcm_type), Intent(InOut)  :: T

    Integer(Kind=wi)   :: ncycles 
    Integer(Kind=wi)   :: maxpoints 
    Integer(Kind=wi)   :: col 
    Character(Len=256) :: message
    Character(Len=256) :: error

    Integer(Kind=wi) :: fail(4), fail2

    ncycles=T%ncycles
    maxpoints=T%max_points
    col=T%columns 

    error='***ERROR: Allocation problems for'

    Allocate(T%segment_points(2, ncycles)           , Stat=fail(1))
    Allocate(T%whole_set(col,maxpoints, 2, ncycles) , Stat=fail(2))
    Allocate(T%voltage%value(maxpoints, 2, ncycles)  , Stat=fail(3))
    Allocate(T%label_leg(2, ncycles)                , Stat=fail(4))

    If (Any(fail > 0)) Then
      Write (message,'(2(1x,a))') Trim(error), 'arrays to extract EQCM data& 
                                  & (subroutine allocate_eqcm_data_arrays)'
      Call error_stop(message)
    End If

    T%whole_set       = 0.0_wp
    T%voltage%value = 0.0_wp

    If (T%current%fread) Then
      Allocate(T%current%value(maxpoints, 2,ncycles)     , Stat=fail(1))
      Allocate(T%charge%value(maxpoints, 2,ncycles)      , Stat=fail(2))
      If (Any(fail > 0)) Then
        Write (message,'(2(1x,a))') Trim(error), 'EQCM current and charge arrays&
                                  & (subroutine allocate_eqcm_data_arrays)'
        Call error_stop(message)
      End If
      T%current%value = 0.0_wp
      T%charge%value = 0.0_wp
    Else
      If (T%charge%fread) Then
        Allocate(T%charge%value(maxpoints, 2,ncycles)      , Stat=fail2)
        If (Any(fail > 0)) Then
          Write (message,'(2(1x,a))') Trim(error), 'EQCM charge array&
                                    & (subroutine allocate_eqcm_data_arrays)'
          Call error_stop(message)
        End If
        T%charge%value = 0.0_wp
      End IF 
    End If

    If (T%mass_frequency%fread) Then
      Allocate(T%mass_frequency%value(maxpoints, 2,ncycles) , Stat=fail2)
      If (fail2 > 0) Then
        Write (message,'(2(1x,a))') Trim(error), 'EQCM mass_frequency array&
                                  & (subroutine allocate_eqcm_data_arrays)'
        Call error_stop(message)
      End If
      T%mass_frequency%value = 0.0_wp
    End If

    If (T%mass%fread) Then
      Allocate(T%mass%value(maxpoints, 2,ncycles) , Stat=fail2)
      If (fail2 > 0) Then
        Write (message,'(2(1x,a))') Trim(error), 'EQCM mass array&
                                  & (subroutine allocate_eqcm_data_arrays)'
        Call error_stop(message)
      End If
      T%mass%value = 0.0_wp
    End If

    If (T%time%fread) Then
      Allocate(T%time%value(maxpoints, 2,ncycles)   , Stat=fail2)
      If (fail2 > 0) Then
        Write (message,'(2(1x,a))') Trim(error), 'time array&
                                  & (subroutine allocate_eqcm_data_arrays)'
        Call error_stop(message)
      End If
      T%time%value = 0.0_wp
    End If

  End Subroutine allocate_eqcm_data_arrays


  Subroutine cleanup(T)
    Type(eqcm_type) :: T

    If (Allocated(T%segment_points)) Then
      Deallocate(T%segment_points)
    End If

    If (Allocated(T%whole_set)) Then
      Deallocate(T%whole_set)
    End If

    If (Allocated(T%voltage%value)) Then
      Deallocate(T%voltage%value)
    End If

    If (Allocated(T%time%value)) Then
      Deallocate(T%time%value)
    End If

    If (Allocated(T%mass_frequency%value)) Then
      Deallocate(T%mass_frequency%value)
    End If

    If (Allocated(T%mass%value)) Then
      Deallocate(T%mass%value)
    End If

    If (Allocated(T%current%value)) Then
      Deallocate(T%current%value)
    End If

    If (Allocated(T%range_cycles%value)) Then
      Deallocate(T%range_cycles%value)
    End If

    If (Allocated(T%voltage_range%value)) Then
      Deallocate(T%voltage_range%value)
    End If

    If (Allocated(T%voltage_range%units)) Then
      Deallocate(T%voltage_range%units)
    End If

    If (Allocated(T%scan_rate%value)) Then
      Deallocate(T%scan_rate%value)
    End If

    If (Allocated(T%scan_rate%units)) Then
      Deallocate(T%scan_rate%units)
    End If

  End Subroutine cleanup
  

  Subroutine read_eqcm_data(files, eqcm_data, maxdim)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set of subroutines to extract information from the experimental EQCM data
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(eqcm_type),   Intent(InOut) :: eqcm_data  
    Integer(Kind=wi),  Intent(In   ) :: maxdim

    Integer(Kind=wi)   :: iunit
    Logical            :: safe
    Character(Len=256) :: unit_name
    Character(Len=256) :: message

    ! Open the SET_EQCM file with EQCM settings
    Inquire(File=files(FILE_DATA_EQCM)%filename, Exist=safe)

    ! Check if file exists
    If (.not.safe) Then
      Write (message,'(3(1x,a))') '***ERROR - File', Trim(files(FILE_DATA_EQCM)%filename), 'not found'
      Call error_stop(message)
    Else
      Open(Newunit=files(FILE_DATA_EQCM)%unit_no, File=files(FILE_DATA_EQCM)%filename,Status='old')
      ! Assign unit to inuit for clearer coding, as using files(FILE_DATA_EQCM)%unit_no is cumbersome
      iunit=files(FILE_DATA_EQCM)%unit_no
      unit_name=Trim(files(FILE_DATA_EQCM)%filename)
    End If

    Call eqcm_variables_units(iunit, unit_name, eqcm_data)

    Call eqcm_cycles_points(iunit, unit_name, eqcm_data, maxdim)

    Call eqcm_values(iunit, unit_name, eqcm_data)

    Call eqcm_name_process(eqcm_data)

    ! Close file
    close(files(FILE_DATA_EQCM)%unit_no)

  
  End Subroutine read_eqcm_data


  Subroutine eqcm_variables_units(iunit, unit_name, eqcm_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read involved variables and units from the headers of DATA_EQCM
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   )   :: iunit
    Character(Len=64),  Intent(In   )   :: unit_name
    Type(eqcm_type),    Intent(InOut)   :: eqcm_data
    

    Character(Len=256)               :: message
    Character(Len=256)               :: messages(6) 
    Character(Len=256)               :: header, first_line
    Character(Len= 32)               :: error_set
    Character(Len= 32)               :: check_header

    Integer(Kind=wi) :: i, j, ic 
    Integer(Kind=wi) :: io, io2

    Logical          :: cont
    Real(Kind=wp)    :: var

    Character(Len= 16)  :: units(max_num_var)
    Character(Len= 64)  :: variables(max_num_var)
    Character(Len= 64)  :: heads(max_num_var)
    Real(Kind=wp)       :: element(max_num_var)
    Logical             :: column_fail(max_num_var)

    Character(Len=*), Parameter :: fmt1 = '(5(1x,a))'

    error_set = '***ERROR in file '//Trim(unit_name)//' -'
    check_header='Check labelling for header' 

    ic=0
    cont=.True.
    Rewind iunit

    Do While (cont)
      Read (iunit,'(a)', iostat=io) header
      Call remove_front_tabs(header)
      header=Trim(Adjustl(header))
      If (.Not. is_iostat_end(io)) Then
        If (header(1:1)==' ') Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), 'Heading MUST be in the first line (no previous empty lines)'
          Call info(messages, 1)
          Call error_stop(' ')
        End If 
        cont=.False.
      Else If (is_iostat_end(io)) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'Empty file?'
        Call error_stop(message)
      End If 
    End Do

    cont=.True.

    Do While (cont)
      Read (iunit,'(a)', iostat=io) first_line
      Call detect_rubbish(first_line, unit_name)
      Call remove_front_tabs(first_line)
      first_line=Trim(Adjustl(first_line))
      If (.Not. is_iostat_end(io)) Then
         If ( (first_line(1:1) == ' ')) Then
         Else
          Read (first_line, * , iostat=io2) var 
          If (io2 == 0) Then
            cont=.False. 
          End If
         End If
      Else If (is_iostat_end(io)) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'Empty file? Also check the use commas&
                                  & to represent values with decimals. If so, commas should be replaced with points'
        Call error_stop(message)
      End If 
    End Do
 
    ic=0
    Do 
      Read (first_line, *, iostat=io) (element(i), i=1,ic+1)
      If (io/=0) Then
        exit
      Else
        Backspace iunit
        ic=ic+1
        If (ic == max_num_var)Then
          Write (message,'(2(1x,a))') Trim(error_set), 'Problems detected in the first line that contains numerical& 
                                     & data. Please check the file.'
          Call error_stop(message)
        End If
      End IF
    End Do
    eqcm_data%columns=ic

    Call info(' ', 1) 
    Write (messages(1),'(1x,a)') 'Structure and details of EQCM data'
    Write (messages(2),'(1x,a)') '=================================='
    Call info(messages,2)
    Write (messages(1),'(1x,a,i2,a)') 'The code has detected ',  eqcm_data%columns, ' EQCM variables (columns)& 
                                     & in file DATA_EQCM.'
    Call info(messages, 1)
    If (eqcm_data%columns==0) Then
       Call info(' ',1)
       Write (message,'(2(1x,a))') Trim(error_set), 'Inconsistencies between the heading and data structure.'
       Call error_stop(message)
    End If

    column_fail=.False.
    variables=' '
    units=' ' 

    ! Remove separators (\/) and parenthesis "{[()]}" (if present)
    Call remove_symbols(header,'/|\#')
    Call remove_symbols(header,'{[()]}')

    Read (header, Fmt=*, iostat=io) (variables(i), units(i),  i=1,eqcm_data%columns) 
    Do i =1, eqcm_data%columns
      Call capital_to_lower_case(variables(i))
      Call capital_to_lower_case(units(i))
    End Do 

    ic=0

    Do i=1, eqcm_data%columns
      ! Check there is no missing names for variables
      If (variables(i)==' ') Then
        Write (message,'(1x,2(a,1x),i2,(1x,a))')  Trim(error_set), 'Missing variable name for the header of column',i,&
                                                  & '. Please check the structure of the header' 
        Call error_stop(message) 
      End If

      ! Check there is no missing names for units
      If (units(i)==' ') Then
         Write (message,'(3(1x,a))') Trim(error_set), 'Missing units for one of the variables!', check_header
        Call error_stop(message)
     End If

      If (Index(variables(i),'time') /= 0) Then
        heads(i)= 'time'
        eqcm_data%time%col  = i
        eqcm_data%time%fread = .True.
        eqcm_data%time%units = units(i)
        If (eqcm_data%time%units(1:1) /= 's') Then
          Write (message,fmt1) Trim(error_set), 'Specified units "', Trim(eqcm_data%time%units), &
                              & '" for time is different from "s" (seconds).', check_header
          Call error_stop(message)
        End If
        ic=ic+1
   
      Else If (Index(variables(i),'potential') /= 0 .Or. Index(variables(i),'voltage') /= 0) Then
        heads(i)= 'voltage'
        eqcm_data%voltage%col  = i
        eqcm_data%voltage%fread = .True.
        eqcm_data%voltage%units = units(i)
        If (eqcm_data%voltage%units /= 'v') Then
          Write (message,fmt1) Trim(error_set), 'Specified units "', Trim(eqcm_data%voltage%units), &
                                   & '" for potential/voltage is different from "V" (Volts).', check_header 
          Call error_stop(message)
        End If
        ic=ic+1

      Else If (Index(variables(i),'current')   /= 0) Then
        heads(i)= 'current'
        eqcm_data%current%col    = i
        eqcm_data%current%fread   = .True.
        eqcm_data%current%units   = units(i)
        If (eqcm_data%current%units  == 'a') Then
          eqcm_data%current%convert = 1000.0_wp
        Else If (eqcm_data%current%units  == 'ma') Then
          eqcm_data%current%convert = 1.0_wp
        Else If (eqcm_data%current%units  == 'ua') Then
          eqcm_data%current%convert = 0.001_wp
        Else If (eqcm_data%current%units  == 'na') Then
          eqcm_data%current%convert = 0.000001_wp
        Else
          Write (message,fmt1) Trim(error_set), 'Specified units "', Trim(eqcm_data%current%units), &
                                    & '" for current is not valid. Options: A, mA (mili), uA (micro) or nA (nano).',&
                                    & check_header 
          Call error_stop(message)
        End If
        eqcm_data%current%units   = 'mA'
        ic=ic+1

      Else If (Index(variables(i),'charge') /= 0) Then
        heads(i)= 'charge'
        eqcm_data%charge%col     = i
        eqcm_data%charge%fread    = .True.
        eqcm_data%charge%units    = units(i)
        If (eqcm_data%charge%units == 'c') Then
          eqcm_data%charge%convert  = 1000.0_wp
        Else If (eqcm_data%charge%units   == 'mc') Then
          eqcm_data%charge%convert  = 1.0_wp
        Else If (eqcm_data%charge%units   == 'uc') Then
          eqcm_data%charge%convert  = 0.001_wp
        Else If (eqcm_data%charge%units   == 'nc') Then
          eqcm_data%charge%convert  = 0.000001_wp
        Else
          Write (message,fmt1) Trim(error_set), 'Specified units "', Trim(eqcm_data%charge%units),&
                                    & '" for charge is not valid. Options: C, mC (mili), uC (micro) or nC (nano). ',&
                                    & check_header 
          Call error_stop(message)
        End If
        eqcm_data%charge%units    = 'mC'
        ic=ic+1

      Else If (Index(variables(i),'mass') /= 0) Then
        heads(i)= 'mass'
        eqcm_data%mass%col     = i
        eqcm_data%mass%fread    = .True.
        eqcm_data%mass%units    = units(i)
        If (eqcm_data%mass%units == 'mg') Then
          eqcm_data%mass%convert  = 1000000.0_wp
        Else If (eqcm_data%mass%units   == 'ug') Then
          eqcm_data%mass%convert  = 1000.0_wp
        Else If (eqcm_data%mass%units   == 'ng') Then
          eqcm_data%mass%convert  = 1.0_wp
        Else If (eqcm_data%mass%units   == 'g') Then
          eqcm_data%mass%convert  = 1.0E+9_wp
        Else
          Write (message,fmt1) Trim(error_set), 'Specified units "', Trim(eqcm_data%mass%units),&
                                    & '" for "mass" is not valid. Options: g(grams), mg (miligrams), ug (micrograms)&
                                    & or ng (nanograms). ', check_header 
          Call error_stop(message)
        End If
        eqcm_data%charge%units    = 'ng'
        ic=ic+1

      Else If (Index(variables(i),'frequency') /= 0) Then 
        heads(i)= 'mass-frequency'
        eqcm_data%mass_frequency%col  = i
        eqcm_data%mass_frequency%fread = .True.
        eqcm_data%mass_frequency%units = units(i)
        If (eqcm_data%mass_frequency%units   == 'hz') Then
          eqcm_data%mass_frequency%convert  = 1.0_wp
        Else If (eqcm_data%mass_frequency%units   == 'v') Then
          If (eqcm_data%V_to_Hz%value <= epsilon(eqcm_data%V_to_Hz%value)) Then
            Write (message,'(1x,a)') 'Units of mass-frequency are given in Volts,&
                                    & and this requires of the V_to_Hz conversion factor.'
            Call info(message,1) 
            Write (message,'(1x,a)') '***ERROR - No specified value of V_to_Hz in file SET_EQCM.' 
            Call error_stop(message)
          Else 
            eqcm_data%mass_frequency%convert  = eqcm_data%V_to_Hz%value
          End If 
        Else
          Write (message,fmt1) Trim(error_set), 'Specified units "', Trim(eqcm_data%mass_frequency%units), &
                                    & '" for mass-frequency is not valid. Options: V or Hz.', check_header
          Call error_stop(message)
        End If
        eqcm_data%mass_frequency%units = 'Hz'
        ic=ic+1

      Else If (Index(variables(i),'resistance') /= 0) Then
        heads(i)= 'resistance'
        eqcm_data%resistance%col     = i
        eqcm_data%resistance%fread    = .True.
        eqcm_data%resistance%units    = units(i)
        
        If (eqcm_data%resistance%units /= 'v' .And. eqcm_data%resistance%units /= 'ohm') Then
          Write (message,'(1x,3a)') Trim(error_set), 'Specified units for resistance is different from "V" (Volts)&
                                  & and "Ohm". ', check_header
          Call error_stop(message)
        End If
        ic=ic+1

      Else
        Write (message,'(1x,2(a,1x),i2,3a)')  Trim(error_set), 'Name of variable for the header of column', i, &
                                                 & ' is NOT recognised: "', Trim(variables(i)), &
                                                 & '". Please check the labels for the columns.'
        column_fail(i)= .True.
        Call error_stop(message) 
      End If

    End Do

    If (ic==0) Then
      Call info(' ',1)
      Write (message,'(2(1x,a))') Trim(error_set), 'Please check if header exists or the structure of header is correct'
      Call error_stop(message)
    End If


    Write (messages(1),'(1x,a)') 'Description of the EQCM information provided by each column:'
    Call info(messages, 1)
    ! Print information identified in each column 
    Write (message,'(1x,a)') '----------------------------------'
    Call info(message, 1)
    Write (message,'(1x,2(a,3x),a)') 'Column', '|', 'Variable'
    Call info(message, 1)
    Write (message,'(1x,a)') '----------------------------------'
    Call info(message, 1)
    Do i=1, eqcm_data%columns
      If (.Not.column_fail(i)) Then
        Write (message,'(1x,i2,7x,a,3x,a)') i, '|' , heads(i)
        Call info(message, 1)
      End If
    End Do
    Write (message,'(1x,a)') '----------------------------------'
    Call info(message, 1)

   ! If potential (or voltage) has not been found, abort
   If (.Not. eqcm_data%voltage%fread) Then
     Write (message,'(4(1x,a))') Trim(error_set), 'No column for potential or voltage has been identified.',& 
                               & 'EQCM analysis is not possible without voltage data.', check_header
     Call error_stop(message) 
   End If
   
   !Checking columns are not duplicated
    Do i = 1, eqcm_data%columns-1
      Do j = i+1, eqcm_data%columns
        If (Trim(heads(i)) == Trim(heads(j)) .And. &
          (.Not. column_fail(i)) .And. (.Not. column_fail(j)) ) Then 
          Write (message,'(1x,a,2(a,i2),a)') Trim(error_set), 'Columns ', i, ' and ', j, ' have header labels &
                                         &that suggest they correspond to the same variable. Please check and correct!'
          Call error_stop(message)
        End If
      End Do
    End Do

    If (eqcm_data%mass_frequency%fread .And. eqcm_data%mass%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)')      '*** IMPORTANT **************************************'
      Write (messages(2),'(1x,a)')      '*** Mass-frequency and mass have both been recorded.'
      If (eqcm_data%analysis%type == 'mass_calibration' ) Then
        Write (messages(3),'(1x,a)')    '*** The requested analysis will be executed using mass-frequency'
      Else        
        Write (messages(3),'(1x,a)')    '*** The requested analysis will be executed using mass'
      End If
      Write (messages(4),'(1x,a)')      '****************************************************'
      Call info(messages, 4)
    End If

  End Subroutine eqcm_variables_units


  Subroutine eqcm_cycles_points(iunit, unit_name, eqcm_data, maxdim)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Identify the number of cylces and the number of points per cycle
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Character(Len=256), Intent(In   ) :: unit_name
    Type(eqcm_type),    Intent(InOut) :: eqcm_data
    Integer(Kind=wi),            Intent(In   ) :: maxdim

    ! Local variables
    Character(Len=256) :: message
    Character(Len=256) :: messages(4), line(3)
    Character(Len=256) :: first_line, header 
    Integer(Kind=wi)   :: i, j, k, m
    Integer(Kind=wi)   :: ioh, ioh2, iocol, iclc
    Integer(Kind=wi)   :: io(3), iostr(3) 
    Integer(Kind=wi)   :: icycles, isteps
    Integer(Kind=wi)   :: cvpoints(2,fftnmax)
    Integer(Kind=wi)   :: linesback, ihead 

    Real(Kind=wp)  :: element(3,max_num_var)
    Real(Kind=wp)  :: v(3), var, vchange
    Real(Kind=wp)  :: vstore(maxdim)

    Character(Len= 64) :: error_set
    Logical :: lv(3)
    Logical :: lsteps, lcycles, lcvloop
    Logical :: lendset, cont, read_header
    Logical :: fpass

    lcycles  =.True.
    icycles  = 0
    cvpoints = 0
    read_header = .True.

    error_set = '***ERROR in file '//Trim(unit_name)//' -' 

    Write (messages(1),'(2(1x,a))') Trim(error_set), 'It seems there are wrong settings with the file....'
    Write (messages(2),'((1x,a))')  '1) Remove all non-ascii characters from file before running the code.'  
    Write (messages(3),'(3(1x,a))') 'In Linux, this can be done by copying', Trim(unit_name), &
                                  &'to file input_file.dat, for example, and then execute:'
    Write (messages(4),'(2(1x,a))') "tr -cd '\11\12\40-\176' < input_file.dat >", Trim(unit_name)                      
 
    Rewind iunit
    Read (iunit,'(a)') header
    If (header=="") Then
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'Heading MUST be in the first line (no previous empty lines)'
      Call info(messages, 1)
      Call error_stop(' ') 
    End If
    Rewind iunit

    Do While (lcycles) 
      lcvloop=.True.
      j=1
      Do While ((j <= 2) .And. lcvloop)
        isteps = 0
        lsteps= .True.         
        cont=.True.
        lendset=.False. 
        ihead=0    

        Do While (cont)
          Read (iunit,'(a)', iostat=ioh) first_line
          Call remove_front_tabs(first_line)
          first_line=Trim(Adjustl(first_line))
          If (.Not. is_iostat_end(ioh)) Then
              If ( (first_line(1:1) == ' ')) Then
              Else
              Read (first_line, * , iostat=ioh2) var
              If (ioh2 == 0) Then
                cont=.False. 
                Backspace iunit
              Else If (ioh2 > 0) Then
                ihead=ihead+1
                fpass=.False.
                If (icycles==0) Then
                  If (j==1 .And. ihead > 2) Then
                    fpass=.True.
                  ElseIf (j==2 .And. ihead > 1) Then
                    fpass=.True.
                  End If
                Else
                  If (ihead > 1) Then
                    fpass=.True.
                  End If 
                End If

                If (fpass) Then
                  lsteps=.False.
                  cont=.False.
                   Backspace iunit
                End If 
              End If
            End If
          Else If (is_iostat_end(ioh)) Then
            Write (message,'(2(1x,a))') Trim(error_set), 'End of file?'
            Call error_stop(message)
          End If 
        End Do
        iclc=1
        vstore=0.0_wp
 
        Do While (lsteps)
            cont=.True.
            lv=.False.
            Do While (cont)
              Read (iunit, '(a)', iostat=io(1)) line(1)
              Call detect_rubbish(line(1), unit_name)
              Call remove_front_tabs(line(1))
              line(1)=Trim(Adjustl(line(1)))
              If (.Not. is_iostat_end(io(1))) Then
                If (line(1)(1:1) == ' ') Then
                Else
                  Read (line(1), *, iostat=iostr(1)) (element(1,i),  i=1,eqcm_data%columns)
                  If (iostr(1) == 0) Then
                    v(1)= element(1,eqcm_data%voltage%col)
                    Do m= iclc-1, 1, -1
                      If ( Abs(v(1)-vstore(m)) < epsilon(0.0_wp) ) Then 
                        Call info(messages,4)
                        Call error_stop(' ')  
                      End If
                    End Do
                    vstore(iclc)=v(1)
                    iclc=iclc+1 
                    lv(1)=.True. 
                    isteps=isteps+1
                  Else
                    Read (line(1), *, iostat=iocol) (element(1,i),  i=1, 2)
                    If (iocol == 0) Then
                       Write (messages(1),'(2(1x,a),i1,a,i3)') Trim(error_set), 'Missing values somewhere in segment ', j, &
                                                              &' of cycle ', icycles+1
                       Write (messages(2),'(1x,a,i2,a)')        'According to the heading, there should be',&
                                                              & eqcm_data%columns, ' values per line containing numerical data.'  
                       Write (messages(3),'(1x,a)')             'Please check file.'
                       Call info(messages,3)
                       Call error_stop(' ')           
                    Else
                      If (line(1)==header) Then
                        Call info(messages,4)
                        Write (messages(1),'((1x,a))') '2) If the problem still persists, look at the headers that separate'
                        Write (messages(2),'((1x,a))') 'each CV set of data and remove them all (except the first header)'
                        Write (messages(3),'((1x,a))') 'as recommended in section "Setting input files" of the manual.'
                        Call info(messages,3)
                        Call error_stop(' ')  
                      End If
                      Backspace iunit
                      lsteps=.False.
                    End If
                  End If
                  cont=.False. 
                End If
              Else If (is_iostat_end(io(1))) Then
                lcycles=.False.
                lsteps=.False.
                lcvloop=.False.
                cont=.False.
                lendset=.True.  
              End If
            End Do
          
            If ( (.Not. lendset) .And. lsteps) Then
              k=2 
              linesback=0
              Do While (k <= 3)
                cont=.True.
                Do While (cont)
                  Read (iunit, '(a)', iostat=io(k)) line(k)
                  Call detect_rubbish(line(k), unit_name)
                  Call remove_front_tabs(line(k))
                  line(k)=Trim(Adjustl(line(k)))
                  linesback=linesback+1
                  If (.Not. is_iostat_end(io(k))) Then
                    If ( (line(k)(1:1) == ' ') ) Then
                    Else
                      Read (line(k), *, iostat=iostr(k)) (element(k,i),  i=1,eqcm_data%columns)
                      If (iostr(k) == 0) Then
                        v(k)= element(k,eqcm_data%voltage%col)
                        lv(k)=.True. 
                      Else
                        lendset=.True.  
                      End If
                      cont=.False. 
                    End If
                  Else If (is_iostat_end(io(k))) Then
                    cont=.False.
                    lendset=.True.  
                    k=3
                  End If
                End Do
                k=k+1
              End Do
              Do k=1,linesback
                Backspace iunit
              End Do

              If (All(lv)) Then
                vchange = (v(3)-v(2))*(v(2)-v(1))
                If ( vchange < 0.0_wp ) Then
                  lsteps=.False.
                Else If (vchange <= epsilon(vchange)) Then
                  Write (message,'(2(1x,a))') Trim(error_set), 'CV data has two consecutive values with the same voltage. ' 
                  Call error_stop(message)
                End If 
              End If

            End If

        End Do

        cvpoints(j,icycles+1)=isteps
        j=j+1

      End Do

      If (j>2) Then
        icycles=icycles+1
      End If

    End Do

    If (icycles == 0 ) Then
      Write (message,'(1x,a)') '***WARNING - CV data does even complete a cycle. ' 
      Call info(message,1)
    End If

    If (Any(cvpoints(:,icycles+1)>0)) Then
      eqcm_data%ncycles =icycles+1
    Else
      eqcm_data%ncycles =icycles
    End If

    ! Find the maximum number of points meassured among all CV segments
    eqcm_data%max_points = maxval(cvpoints)

    ! Initialise eqcm arrays
    Call eqcm_data%init_arrays() 

    ! Assign points to each CV segment
    Do i =1, eqcm_data%ncycles
      eqcm_data%segment_points(:,i)=cvpoints(:,i)
    End Do


  End Subroutine eqcm_cycles_points

  Subroutine eqcm_values(iunit, unit_name, eqcm_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read values for each variable and cycle. Assign to corresponding arrays
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),            Intent(In   ) :: iunit
    Character(Len=256), Intent(In   ) :: unit_name
    Type(eqcm_type),    Intent(InOut) :: eqcm_data

    Character(Len=256) :: message
    Character(Len=256) :: first_line
    Integer(Kind=wi)   :: i, j, k, l
    Integer(Kind=wi)   :: ioh, ioh2
   
    Real(Kind=wp)      :: var
    Logical            :: cont 

    Character(Len= 64)                :: error_set 

    error_set = '***ERROR in file '//Trim(unit_name)//' -' 

    Rewind iunit

    Do i = 1, eqcm_data%ncycles 
      Do j= 1, 2
        Do k = 1, eqcm_data%segment_points(j,i)
          cont=.True.
          Do While (cont)
            Read (iunit,'(a)', iostat=ioh) first_line
            Call remove_front_tabs(first_line)
            first_line=Trim(Adjustl(first_line))
            If (.Not. is_iostat_end(ioh)) Then
              If ( (first_line(1:1) == ' ') ) Then
              Else
                Read (first_line,*, iostat=ioh2) var
                If (ioh2 == 0) Then

                  Read (first_line,*) (eqcm_data%whole_set(l,k,j,i),  l=1, eqcm_data%columns)
                  eqcm_data%voltage%value(k,j,i) = eqcm_data%whole_set(eqcm_data%voltage%col, k, j, i)

                  If (eqcm_data%charge%fread) Then
                    eqcm_data%charge%value(k,j,i) = eqcm_data%charge%convert*eqcm_data%whole_set(eqcm_data%charge%col, k, j, i)
                  End If
                  
                  If (eqcm_data%current%fread) Then
                    eqcm_data%current%value(k,j,i) = eqcm_data%current%convert * eqcm_data%whole_set(eqcm_data%current%col, k, j, i)
                  End If
                  If (eqcm_data%time%fread) Then
                    eqcm_data%time%value(k,j,i) = eqcm_data%whole_set(eqcm_data%time%col, k, j, i)
                  End If
                  If (eqcm_data%mass_frequency%fread) Then
                    eqcm_data%mass_frequency%value(k,j,i) = eqcm_data%mass_frequency%convert*&
                                                    & eqcm_data%whole_set(eqcm_data%mass_frequency%col, k, j, i)
                  End If
                  If (eqcm_data%mass%fread) Then
                    eqcm_data%mass%value(k,j,i) = eqcm_data%mass%convert*&
                                                 & eqcm_data%whole_set(eqcm_data%mass%col, k, j, i)
                  End If

                  cont=.False.

                Else If (ioh2 > 0) Then
                End If
              End If
            Else If (is_iostat_end(ioh)) Then
              Write (message,'(2(1x,a))') Trim(error_set), 'End of file'
              Call error_stop(message)
            End If 
          End Do
        End Do
      End Do
    End Do

  End Subroutine eqcm_values

  Subroutine eqcm_check_voltage_range(eqcm_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if voltage range is consistent with the EQCM data 
    !
    ! author    - i.scivetti July 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),      Intent(InOut) :: eqcm_data

    Integer(Kind=wi)                :: i, j, k, fail(1)
    Integer(Kind=wi), Allocatable   :: range_points(:,:)
    Character(Len=256)              :: message

    Allocate(range_points(2,eqcm_data%ncycles), Stat=fail(1))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for "range_points" array&
                               & in subroutine eqcm_check_voltage_range'
      Call error_stop(message)
    End If

    range_points=0

    Do i = eqcm_data%range_cycles%value(1), eqcm_data%range_cycles%value(2)
      Do j = 1, 2
        Do k = 1, eqcm_data%segment_points(j,i)
          If ((eqcm_data%voltage_range%value(1) < eqcm_data%voltage%value(k,j,i)) .And. &
             (eqcm_data%voltage_range%value(2) > eqcm_data%voltage%value(k,j,i)) ) Then
            range_points(j,i)=range_points(j,i)+1
          End If
        End Do
      End Do
    End Do

    If (.Not. Any(range_points>0)) Then
      Write (message,'(1x,a)') '***ERROR - No CV values found within the domain defined by voltage_range'
      Call error_stop(message)     
    End If 

    Do i = eqcm_data%range_cycles%value(1), eqcm_data%range_cycles%value(2)
      Do j = 1, 2
        If (eqcm_data%segment_points(j,i) > 2 .And. range_points(j,i)==0) Then
          Write (message,'(1x,a,i3,a,1x,a)') '***WARNING: Cycle ', i, ' ('//Trim(eqcm_data%label_leg(j,i))//')-', &
                                           'No CV values found within the domain defined by voltage_range' 
          Call info(message,1)
        End If
      End Do
    End Do

    Deallocate(range_points)

  End Subroutine eqcm_check_voltage_range


  Subroutine eqcm_name_process(eqcm_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Identify if the voltage sweep is poitive or negative (or undefined)
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),      Intent(InOut) :: eqcm_data

    Integer(Kind=wi)    :: i, j
    Real(Kind=wp)       :: dv
    Logical             :: fpass

   fpass=.True. 

    Do i = 1, eqcm_data%ncycles
      Do j= 1, 2
        If (eqcm_data%segment_points(j,i) < 2) Then
          eqcm_data%label_leg(j,i)='undefined'
        Else
          dv=eqcm_data%voltage%value(2,j,i)-eqcm_data%voltage%value(1,j,i)
          If (dv>0.0_wp) Then
            eqcm_data%label_leg(j,i)= 'positive'
          ElseIf (dv<0.0_wp) Then
            eqcm_data%label_leg(j,i)= 'negative'
          End If 
        End If
      End Do
    End Do

    i=1
    eqcm_data%scan_sweep_ref='undefined'

    Do While (i<=eqcm_data%ncycles .And. fpass)
      If (eqcm_data%label_leg(1,i) == 'undefined') Then
        i=i+1
      Else
        eqcm_data%scan_sweep_ref= eqcm_data%label_leg(1,i)
        fpass=.False.
      End If
    End Do


  End Subroutine eqcm_name_process


  Subroutine eqcm_cycling_summary(files, eqcm_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Summarise the information found for each of the chosen cycles by 
    ! directive "cycles" in SET_EQCM 
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(eqcm_type),    Intent(InOut) :: eqcm_data

    Integer(Kind=wi) :: i, j
    Character(Len=256) :: append, remove_ascii
    Character(Len=256) :: message
    Character(Len=*), Parameter :: fmt1 = '(1x,4(a,4x))'
    Character(Len=*), Parameter :: fmt2 = '(1x,i3,7x,a,7x,i4,6x,f7.3,a,f7.3)'    
    Character(Len=*), Parameter :: fmt3 = '(1x,  10x,a,7x,i4,6x,f7.3,a,f7.3)'    
    Character(Len=*), Parameter :: fmt4 = '(1x,i3,7x,a)'    
    Character(Len=*), Parameter :: fmt5 = '(1x,  10x,a)'    

    append= '. Provided EQCM data is NOT consistent with CV cycling!'
    remove_ascii=' ***IMPORTANT: File DATA_EQCM might contain non-ASCII characters.&
                   & Please make sure to remove them.'   
 
    Call info(' ', 1)
    Write (message,'(1x,a,i3,a)') 'A total of ', eqcm_data%ncycles,&
                                  ' CV cycles have been identified from the EQCM data'
    Call info(message, 1)
    Write (message,'(1x,a)') '-----------------------------------------------------------------'
    Call info(message, 1)
    Write (message, fmt1 ) 'Cycle', 'Voltage sweep', 'Points', 'Voltage range [V]'
    Call info(message, 1)
    Write (message,'(1x,a)') '-----------------------------------------------------------------'
    Call info(message, 1)
    Do i = 1, eqcm_data%ncycles 
      Do j= 1, 2
        If (j==1) Then
          If (eqcm_data%label_leg(j,i) /= 'undefined') Then
            Write (message, fmt2) i, Trim(eqcm_data%label_leg(j,i)), eqcm_data%segment_points(j,i), &
                               eqcm_data%voltage%value(1,j,i), ' --> ' , &
                               eqcm_data%voltage%value(eqcm_data%segment_points(j,i),j,i)
          Else
            Write (message, fmt4) i, Trim(eqcm_data%label_leg(j,i))  
          End If
        Else
          If (eqcm_data%label_leg(j,i) /= 'undefined') Then
            Write (message, fmt3)    Trim(eqcm_data%label_leg(j,i)), eqcm_data%segment_points(j,i), &
                               eqcm_data%voltage%value(1,j,i), ' --> ' , &
                               eqcm_data%voltage%value(eqcm_data%segment_points(j,i),j,i)
          Else
            Write (message, fmt5) Trim(eqcm_data%label_leg(j,i))  
          End If
        End If
         Call info(message, 1)
      End Do
    End Do
    Write (message,'(1x,a)') '-----------------------------------------------------------------'
    Call info(message, 1)

    Do i= 1, eqcm_data%ncycles
      Do j=1,2
        If (eqcm_data%label_leg(j,i) == 'undefined') Then
          Write (message,'(1x,a,1x,i3,1x,a)') '***WARNING - Cycle ', i, ' is NOT complete'
          Call info(message, 1)
        End If         
      End Do
    End Do

    Do i= 1, eqcm_data%ncycles
        If (i<eqcm_data%range_cycles%value(2)) Then
          If ( (eqcm_data%label_leg(1,i) == eqcm_data%label_leg(2,i)) .And. &
              (eqcm_data%label_leg(1,i) /= 'undefined') ) Then
            Write (message,'(3(1x,a),i3,a)') '***ERROR- Two consecutive', Trim(eqcm_data%label_leg(1,i)), &
                                           &'Voltage sweeps are found for cycle', i, Trim(append)
            Call info(message, 1)
            Call error_stop(remove_ascii)
          End If
          If ((eqcm_data%label_leg(2,i) == eqcm_data%label_leg(1,i+1)) .And. &
                  (eqcm_data%label_leg(2,i) /= 'undefined') ) Then
            Write (message,'(2(1x,a),2(a,i3),a)') '***ERROR- Two consecutive', Trim(eqcm_data%label_leg(1,i)), &
                                              & 'Voltage sweeps between cycles ', i, ' and ', i+1, Trim(append)
            Call info(message, 1)
            Call error_stop(remove_ascii)
          End If         
        Else
          If ( (eqcm_data%label_leg(1,i) == eqcm_data%label_leg(2,i)) .And. &
              (eqcm_data%label_leg(1,i) /= 'undefined') ) Then
            Write (message,'(3(1x,a),i3)') '***ERROR- Two consecutive', Trim(eqcm_data%label_leg(1,i)),&
                                         &' Voltage sweeps for cycle', i
            Call info(message, 1)
            Call error_stop(remove_ascii)
          End If
        End If
    End Do


    ! Refresh out_eqcm
    Call refresh_out_eqcm(files)

  End Subroutine eqcm_cycling_summary


  Subroutine eqcm_spectra(eqcm_data,fft_data, files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Perform spectra analysis for the range of cycles selected in SET_EQCM.
    ! If EQCM mass/mass-frequency has been recorded in DATA_EQCM, results are 
    ! printed to SPEC_MASS.
    ! If EQCM current has been recorded in DATA_EQCM, results are printed 
    ! to SPEC_CURRENT.
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),    Intent(InOut)  :: eqcm_data
    Type(fft_type),     Intent(InOut)  :: fft_data
    Type(file_type),    Intent(InOut)  :: files(:)

    Character(Len=256) :: messages(3)
    Character(Len=256) :: message
    Integer(Kind=wi)   :: i, j, npoints
    Real(Kind=wp)      :: domain 
    Integer(Kind=wi)   :: unit_mass_freq, unit_current
    Integer(Kind=wi)   :: fail(3)

    Call info(' ', 1)
    Call info(' Spectra analysis', 1)
    Call info(' ================', 1)

    ! Set further endpoints specification for FFT
    Write (message,'(1x,a)') 'FFT mesh settings:'
    Call info(message,1)
    If (eqcm_data%current%fread) Then
      Write (message,'(1x,a,i4)') '- data points to set padding (current)        ', &
                               & fft_data%end_current%value
      Call info(message,1)
    End If

    If (eqcm_data%mass%fread) Then
      Write (message,'(1x,a,i4)') '- data points to set padding (mass)           ', &
                               & fft_data%end_mass%value
      Call info(message,1)
    Else    
      If (eqcm_data%mass_frequency%fread) Then
        Write (message,'(1x,a,i4)') '- data points to set padding (mass-frequency) ', &
                                 & fft_data%end_mass%value
        Call info(message,1)
      End If
    End If

    fft_data%endpoints=max(fft_data%end_mass%value, fft_data%end_current%value)

    Write (messages(1),'(a)') '  '
    If (eqcm_data%range_cycles%value(1)==eqcm_data%range_cycles%value(2)) Then
      Write (messages(2),'(1x,a,i3)') 'Analysis for cycle ', eqcm_data%range_cycles%value(1)
    Else
      Write (messages(2),'(1x,(2(a,i3)))') 'Analysis for cycle ', & 
                                       & eqcm_data%range_cycles%value(1), ' to cycle ', eqcm_data%range_cycles%value(2)
    End If
    Call info(messages,2)
 

    If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
      ! Open file SPEC_MASS_FREQUENCY
      Open(Newunit=files(FILE_SPEC_MASS)%unit_no, File=files(FILE_SPEC_MASS)%filename, Status='Replace')
      unit_mass_freq=files(FILE_SPEC_MASS)%unit_no
    End If
    
    If (eqcm_data%current%fread) Then
      ! Open file SPEC_CURRENT
      Open(Newunit=files(FILE_SPEC_CURRENT)%unit_no, File=files(FILE_SPEC_CURRENT)%filename, Status='Replace')
      unit_current=files(FILE_SPEC_CURRENT)%unit_no
    End If

    If (eqcm_data%time%fread) Then
      If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
        Write (unit_mass_freq,'(a,8x,a)') '# Frequency [Hz]', 'Magnitude [a.u.]'
      End If
      If (eqcm_data%current%fread) Then
        Write (unit_current,  '(a,8x,a)') '# Frequency [Hz]', 'Magnitude [a.u.]' 
      End If
    Else
      If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
        Write (unit_mass_freq,'(a,8x,a)') '# Inverse of voltage [1/V]', 'Magnitude [a.u.]' 
      End If
      
      If (eqcm_data%current%fread) Then
        Write (unit_current,'(a,8x,a)')   '# Inverse of voltage [1/V]', 'Magnitude [a.u.]'
      End If
    End If


    Do i = eqcm_data%range_cycles%value(1), eqcm_data%range_cycles%value(2) 
      Do j= 1, 2
        If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
          Write (unit_mass_freq,'(a,i3,3x,2a)') '# Cycle ', i,  'Voltage sweep: ', Trim(eqcm_data%label_leg(j,i))
        End If
        If (eqcm_data%current%fread) Then
          Write (unit_current,  '(a,i3,3x,2a)') '# Cycle ', i,  'Voltage sweep: ', Trim(eqcm_data%label_leg(j,i))
        End If

        If (eqcm_data%segment_points(j,i) > 2*fft_data%endpoints) Then
          If (eqcm_data%time%fread) Then
            domain=Abs(eqcm_data%time%value(1,j,i)-eqcm_data%time%value(eqcm_data%segment_points(j,i),j,i))
          Else
            domain=Abs(eqcm_data%voltage%value(1,j,i)-eqcm_data%voltage%value(eqcm_data%segment_points(j,i),j,i))
          End If
          npoints= eqcm_data%segment_points(j,i)
          Call fft_points(fft_data, npoints)
          ! Allocate fft arrays
          Allocate(fft_data%array(fft_data%ntot), Stat=fail(1))
          Allocate(fft_data%rec_domain(fft_data%ntot/2+1), Stat=fail(2)) 
          Allocate(fft_data%magnitude(fft_data%ntot/2+1),  Stat=fail(3)) 
          If (Any(fail > 0)) Then
            Write (message,'(1x,1a)') '***ERROR: Allocation problems for fft_data arrays in subroutine "eqcm_spectra"'
            Call error_stop(message)
          End If            
          
          If (eqcm_data%current%fread) Then
            Call fft_setup(fft_data, npoints, eqcm_data%current%value(:,j,i), domain, fft_data%end_current%value)
            Call fft_spectrum(fft_data)
            Call print_spectrum(unit_current,fft_data)
          End If
          If (eqcm_data%mass%fread) Then
            Call fft_setup(fft_data, npoints, eqcm_data%mass%value(:,j,i), domain, fft_data%end_mass%value)
            Call fft_spectrum(fft_data)
            Call print_spectrum(unit_mass_freq,fft_data)
          Else        
            If (eqcm_data%mass_frequency%fread) Then
              Call fft_setup(fft_data, npoints, eqcm_data%mass_frequency%value(:,j,i), domain, fft_data%end_mass%value)
              Call fft_spectrum(fft_data)
              Call print_spectrum(unit_mass_freq,fft_data)
            End If
          End If

          ! Deallocate
          Deallocate(fft_data%magnitude)
          Deallocate(fft_data%rec_domain)
          Deallocate(fft_data%array)
        Else
          If (eqcm_data%mass%fread) Then
          Else  
            If (eqcm_data%mass_frequency%fread) Then
              Write (message, '(a,i3,a)') ' *** WARNING for cycle ', i, ': FFT problems with mass-frequency.'
              Call info(message, 1)
              Write (unit_mass_freq,'(a)')   '# Either no data or number of the EQCM points is not sufficient for the '   
              Write (unit_mass_freq,'(2a)')  '# choice of endpoints_mass in ', Trim(files(FILE_SET_EQCM)%filename)
              Write (unit_mass_freq,'(a)')   '  '
            End If
          End If

          If (eqcm_data%current%fread) Then
            Write (message, '(a,i3,a)') ' *** WARNING for cycle ', i, ': FFT problems with current.'
            Call info(message, 1)
            Write (unit_current,'(a)')  '# Either no data or number of the EQCM points is not sufficient for the '   
            Write (unit_current,'(2a)') '# choice of endpoints_current in ',  Trim(files(FILE_SET_EQCM)%filename)
            Write (unit_current,'(a)')  '  '
          End If
        End If 
      End Do
    End Do
    
    Call info(' ', 1)
    If ((eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) .And. eqcm_data%current%fread) Then
      Write (messages(1),'(4(1x,a))') 'Print results to files', Trim(FOLDER_ANALYSIS)//'/'//Trim(files(FILE_SPEC_MASS)%filename),&
                                    & 'and', Trim(FOLDER_ANALYSIS)//'/'//Trim(files(FILE_SPEC_CURRENT)%filename) 
    Else If ((eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) .And. (.Not. eqcm_data%current%fread)) Then
      Write (messages(1),'(2(1x,a))') 'Print results to file', Trim(FOLDER_ANALYSIS)//'/'//Trim(files(FILE_SPEC_MASS)%filename)
    Else If ((.Not. (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread)) .And. eqcm_data%current%fread) Then
      Write (messages(1),'(2(1x,a))') 'Print results to files', Trim(FOLDER_ANALYSIS)//'/'//Trim(files(FILE_SPEC_CURRENT)%filename) 
    End If
    Call info(messages,1)

    ! Close files 
    If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
      ! Close file SPEC_MASS_FREQUENCY
      Close(files(FILE_SPEC_MASS)%unit_no)
      ! Move file
      Call execute_command_line('mv '//Trim(files(FILE_SPEC_MASS)%filename)//' '//Trim(FOLDER_ANALYSIS))
    End If
    
    If (eqcm_data%current%fread) Then
      ! Close file SPEC_CURRENT
      Close(files(FILE_SPEC_CURRENT)%unit_no)
      ! Move file
      Call execute_command_line('mv '//Trim(files(FILE_SPEC_CURRENT)%filename)//' '//Trim(FOLDER_ANALYSIS))
    End If

  End Subroutine eqcm_spectra

  Subroutine eqcm_filter(eqcm_data,fft_data,filter)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Perform filtering of raw EQCM data via a low-pass Gaussian filter with
    ! cutoff frequency defined in SET_EQCM. Filtering is done via a FFT.
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),   Intent(InOut)    :: eqcm_data
    Type(fft_type),    Intent(InOut)    :: fft_data
    Type(filter_type), Intent(InOut)    :: filter

    Character(Len=256) :: messages(3)
    Character(Len=256) :: message
    Integer(Kind=wi)   :: i, j, npoints
    Real(Kind=wp)      :: domain 
    Integer(Kind=wi)   :: fail(2)

    Call info(' ', 1)
    Call info(' Filtering data', 1)
    Call info(' ==============', 1)

    ! Set further endpoints specification for FFT
    Write (message,'(1x,a)') 'FFT mesh for filtering:'
    Call info(message,1)
    If (eqcm_data%current%fread) Then
      Write (message,'(1x,a,i3)') 'Current points to average padding:        ', &
                               & fft_data%end_current%value
      Call info(message,1)
    End If

    If (eqcm_data%mass%fread) Then
      Write (message,'(1x,a,i3)')   'Mass points to average padding:           ', &
                               & fft_data%end_mass%value
      Call info(message,1)
    Else    
      If (eqcm_data%mass_frequency%fread) Then
        Write (message,'(1x,a,i3)') 'Mass-frequency points to average padding: ', &
                                 & fft_data%end_mass%value
        Call info(message,1)
      End If
    End If
    
    fft_data%endpoints=max(fft_data%end_mass%value, fft_data%end_current%value)

    Write (messages(1),'(a)') '  '
    If (eqcm_data%range_cycles%value(1)==eqcm_data%range_cycles%value(2)) Then
      Write (messages(2),'(1x,a,i3)') 'Process data for cycle ', eqcm_data%range_cycles%value(1)
    Else
      Write (messages(2),'(1x,(2(a,i3)))') 'Process data for range between cycle ', & 
                                        & eqcm_data%range_cycles%value(1), ' and cycle ', eqcm_data%range_cycles%value(2)
    End If
    Call info(messages,2)

    Do i = eqcm_data%range_cycles%value(1), eqcm_data%range_cycles%value(2) 
      Do j= 1, 2
        If (eqcm_data%segment_points(j,i) > 2*fft_data%endpoints) Then 
          If (eqcm_data%time%fread) Then
            domain=Abs(eqcm_data%time%value(1,j,i)-eqcm_data%time%value(eqcm_data%segment_points(j,i),j,i))
          Else
            domain=Abs(eqcm_data%voltage%value(1,j,i)-eqcm_data%voltage%value(eqcm_data%segment_points(j,i),j,i))
          End If
          npoints= eqcm_data%segment_points(j,i)
          Call fft_points(fft_data, npoints)
          ! Allocate fft arrays
          Allocate(fft_data%array(fft_data%ntot), Stat=fail(1))
          Allocate(fft_data%rec_domain(fft_data%ntot/2+1), Stat=fail(2)) 
          If (Any(fail > 0)) Then
            Write (message,'(1x,1a)') '***ERROR: Allocation problems for fft_data arrays in subroutine "eqcm_filter"'
            Call error_stop(message)
          End If            

          If (eqcm_data%current%fread) Then
            Call fft_setup(fft_data, npoints, eqcm_data%current%value(:,j,i), domain, fft_data%end_current%value)
            Call filter_data(fft_data, filter, npoints, eqcm_data%current%value(:,j,i))
          End If

          If (eqcm_data%mass%fread) Then
              Call fft_setup(fft_data, npoints, eqcm_data%mass%value(:,j,i), domain, fft_data%end_mass%value)
              Call filter_data(fft_data, filter, npoints, eqcm_data%mass%value(:,j,i))
          Else        
            If (eqcm_data%mass_frequency%fread) Then
              Call fft_setup(fft_data, npoints, eqcm_data%mass_frequency%value(:,j,i), domain, fft_data%end_mass%value)
              Call filter_data(fft_data, filter, npoints, eqcm_data%mass_frequency%value(:,j,i))
            End If
          End If

          ! Deallocate
          Deallocate(fft_data%rec_domain) 
          Deallocate(fft_data%array) 
        Else
          Write (messages(1),'(1x,a,i3,2a)') '*** WARNING: Cycle ', i, '- Voltage sweep: ', Trim(eqcm_data%label_leg(j,i))
          Write (messages(2),'(1x,a)')       'Either no data or number of the EQCM points is not sufficient for the&
                                            & choice of endpoints'
          Write (messages(3),'(1x,a)')      'This part of the EQCM data will not be filtered.' 
          Call info(messages,3)
        End If
      End Do
    End Do

  End Subroutine eqcm_filter

  Subroutine print_eqcm_data(eqcm_data, fft_data, files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print EQCM data if either "print_eqcm_raw" or "print_eqcm_filter" option 
    ! was selected. Current and/or mass (or mass-density) will be printed to file:    
    ! - FILTERED_CURRENT if print_eqcm_filter and there is data for current
    ! - FILTERED_MASS    if print_eqcm_filter and there is data for mass-frequency 
    ! - RAW_CURRENT      if print_eqcm_raw and there is data for current
    ! - RAW_MASS         if print_eqcm_raw and there is data for mass-frequency
    ! 
    ! author    - i.scivetti June 2020
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),   Intent(InOut)    :: eqcm_data
    Type(fft_type),    Intent(InOut)    :: fft_data
    Type(file_type),   Intent(InOut)    :: files(:)

    Integer(Kind=wi)  :: i, j, k, l
    Integer(Kind=wi)  :: iunit_mass, iunit_current

    Logical           :: frange

    Character(Len=256) :: messages(3), header
    Character(Len=256) :: unit_name_mass, unit_name_current 

    If (eqcm_data%analysis%type == 'print_eqcm_filter') Then
      ! Open file FILTERED_CURRENT
      If (eqcm_data%current%fread) Then
        Open(Newunit=files(FILE_FILTERED_CURRENT)%unit_no, File=files(FILE_FILTERED_CURRENT)%filename, Status='Replace')
        iunit_current=files(FILE_FILTERED_CURRENT)%unit_no
        unit_name_current=Trim(files(FILE_FILTERED_CURRENT)%filename)
      End If
      ! Open file FILTERED_MASS
      If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
        Open(Newunit=files(FILE_FILTERED_MASS)%unit_no, File=files(FILE_FILTERED_MASS)%filename, Status='Replace')
        iunit_mass=files(FILE_FILTERED_MASS)%unit_no
        unit_name_mass=Trim(files(FILE_FILTERED_MASS)%filename)
      End If
    Else
      ! Open file RAW_CURRENT
      If (eqcm_data%current%fread) Then
        Open(Newunit=files(FILE_RAW_CURRENT)%unit_no, File=files(FILE_RAW_CURRENT)%filename,Status='Replace')
        iunit_current=files(FILE_RAW_CURRENT)%unit_no
        unit_name_current=Trim(files(FILE_RAW_CURRENT)%filename)
      End If

      ! Open file RAW_MASS
      If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
        Open(Newunit=files(FILE_RAW_MASS)%unit_no, File=files(FILE_RAW_MASS)%filename,Status='Replace')
        iunit_mass=files(FILE_RAW_MASS)%unit_no
        unit_name_mass=Trim(files(FILE_RAW_MASS)%filename)
      End If

    End If

    !Call info(messages,1)
    If (eqcm_data%current%fread .And. (.Not. (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread))) Then
      Write (messages(1),'(1x,2a)') 'Print EQCM current to file ', Trim(FOLDER_ANALYSIS)//'/'//Trim(unit_name_current) 
    Else If ((eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) .And. (.Not. eqcm_data%current%fread)) Then
      If (eqcm_data%mass%fread) Then      
        Write (messages(1),'(1x,2a)') 'Print EQCM mass to file ', Trim(FOLDER_ANALYSIS)//'/'//Trim(unit_name_mass) 
      Else
        If (eqcm_data%mass_frequency%fread) Then      
          Write (messages(1),'(1x,2a)') 'Print EQCM mass-density to file ', Trim(FOLDER_ANALYSIS)//'/'//Trim(unit_name_mass) 
        End If        
      End If        
    Else If (eqcm_data%current%fread .And. (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread)) Then  
      If (eqcm_data%mass%fread) Then      
        header='Print EQCM current and mass to files'
      Else
        If (eqcm_data%mass_frequency%fread) Then
          header='Print EQCM current and mass-density to file'
        End If        
      End If        
      Write (messages(1),'(1x,a,1x,4a)') Trim(header), &
                                   & Trim(FOLDER_ANALYSIS)//'/'//Trim(unit_name_current), ' and ',&
                                   & Trim(FOLDER_ANALYSIS)//'/'//Trim(unit_name_mass), ', respectively.' 
    End If
    Call info(' ', 1)
    Call info(messages,1)

    If (eqcm_data%current%fread) Then
      Write (iunit_current,'(a,6x,a)')    '# Potential [V]', 'Current [mA]'
    End If

    If (eqcm_data%mass%fread) Then
      Write (iunit_mass,'(a,6x,a)')    '# Potential [V]', 'Mass [ng]'
    Else         
      If (eqcm_data%mass_frequency%fread) Then
        Write (iunit_mass,'(a,6x,a)')    '# Potential [V]', 'Mass-density [ng/cm2]'
      End If
    End If

    Do i = eqcm_data%range_cycles%value(1), eqcm_data%range_cycles%value(2)
      Do j= 1, 2

        If (eqcm_data%current%fread) Then
          Write (iunit_current,'(a,i3,2a)') '# Cycle ', i, '- Voltage sweep: ', Trim(eqcm_data%label_leg(j,i))
        End If
        If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
          Write (iunit_mass,'(a,i3,2a)')   '# Cycle ', i, '- Voltage sweep: ', Trim(eqcm_data%label_leg(j,i))
        End If 

        If (eqcm_data%segment_points(j,i) > 2*fft_data%endpoints) Then
          l=eqcm_data%segment_points(j,i)
        Else
          l=0 
        Endif   
        Do k =1, l
            frange=.True.
            If (eqcm_data%voltage_range%fread) Then
              frange= (eqcm_data%voltage_range%value(1) < eqcm_data%voltage%value(k,j,i)) .And. &
                      (eqcm_data%voltage_range%value(2) > eqcm_data%voltage%value(k,j,i)) 
            End If            

            If (frange) Then
              If (eqcm_data%current%fread) Then
                Write (iunit_current,'(f12.4,8x,f12.4)')  eqcm_data%voltage%value(k,j,i), eqcm_data%current%value(k,j,i)
              End If             
              If (eqcm_data%mass%fread) Then
                Write (iunit_mass,'(f12.4,8x,f12.4)')  eqcm_data%voltage%value(k,j,i), eqcm_data%mass%value(k,j,i)
              Else        
                If (eqcm_data%mass_frequency%fread) Then
                  Write (iunit_mass,'(f12.4,8x,f12.4)')  eqcm_data%voltage%value(k,j,i), eqcm_data%mass_frequency%value(k,j,i)
                End If             
              End If             
            End If

        End Do
        If (l==0) Then
          If (eqcm_data%current%fread) Then
            Write (iunit_current,'(a)') '# Either no data or number of the EQCM points&
                                       & is not sufficient for the choice of endpoints'
          End If
          If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
            Write (iunit_mass,'(a)') '# Either no data or number of the EQCM points is&
                                    & not sufficient for the choice of endpoints'
          End If
        End If
        If (eqcm_data%current%fread) Write (iunit_current,*) ' '
        If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Write (iunit_mass,*) ' '
      End Do
    End Do

    ! Close file and move it to the corresponding directory 
    If (eqcm_data%current%fread) Then
      Call execute_command_line('mv '//Trim(unit_name_current)//' '//Trim(FOLDER_ANALYSIS))
      Close(iunit_current)
    End If
    If (eqcm_data%mass_frequency%fread .Or. eqcm_data%mass%fread) Then
      Call execute_command_line('mv '//Trim(unit_name_mass)//' '//Trim(FOLDER_ANALYSIS))
      Close(iunit_mass)
    End If

  End Subroutine print_eqcm_data

  Subroutine get_eqcm_charge(eqcm_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to calculate the charge from the current for each CV cycle. 
    ! If time is not recorded, the scan_rate will be used to convert the 
    ! voltage to time and integrate the current. 
    ! 
    ! author    - i.scivetti October 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),   Intent(InOut)    :: eqcm_data

    Integer(Kind=wi)     :: i, j, k
    Real(Kind=wp)        :: chg_acc, Dt, Dchg
    Real(Kind=wp)        :: ltime, lpot, lcurrent
    Character(Len=256)   :: message

    lcurrent=0.0_wp
    ltime=0.0_wp
    lpot=0.0_wp
    chg_acc=0.0_wp

    If (.Not. eqcm_data%time%fread) Then
      If (.Not. eqcm_data%scan_rate%fread) Then
        Write (message, '(a)') '***ERROR: computation of accumulated charges requires the specification of the&
                              & scan_rate directive (missing), which indicates the voltage scan rate used for&
                              & EQCM experiments.'
        Call error_stop(message)
      End If 
    End If 

    If (eqcm_data%current_offset%stat) Then
      eqcm_data%current%value(:,:,:)=eqcm_data%current%value(:,:,:)-eqcm_data%current%value(1,1,1)
    End If

    Do i = 1, eqcm_data%ncycles
      Do j= 1, 2
        If (eqcm_data%segment_points(j,i) > 2) Then
          If (i==1 .And. j==1) Then
            eqcm_data%charge%value(1,j,i)=0.0_wp
          Else
            If (eqcm_data%time%fread) Then
              Dt=Abs(eqcm_data%time%value(1,j,i)-ltime)
            Else
              Dt=Abs(eqcm_data%voltage%value(1,j,i)-lpot)/eqcm_data%scan_rate%value(1) 
            End If 
            Dchg=Dt*(eqcm_data%current%value(1,j,i)+lcurrent)/2.0_wp
            eqcm_data%charge%value(1,j,i)=Dchg+chg_acc
          End If

          Do k= 2, eqcm_data%segment_points(j,i)
            If (eqcm_data%time%fread) Then
              Dt=Abs(eqcm_data%time%value(k,j,i)-eqcm_data%time%value(k-1,j,i))
            Else
              Dt=Abs(eqcm_data%voltage%value(k,j,i)-eqcm_data%voltage%value(k-1,j,i))/eqcm_data%scan_rate%value(1)
            End If 
            Dchg=Dt*(eqcm_data%current%value(k,j,i)+eqcm_data%current%value(k-1,j,i))/2.0_wp
            eqcm_data%charge%value(k,j,i)=Dchg+eqcm_data%charge%value(k-1,j,i)
          End Do
               
          lpot=eqcm_data%voltage%value(k-1,j,i)
          lcurrent=eqcm_data%current%value(k-1,j,i)
          chg_acc=eqcm_data%charge%value(k-1,j,i) 
          If (eqcm_data%time%fread) Then
            ltime=eqcm_data%time%value(k-1,j,i)
          End If
        End If
      End Do
    End Do


  End Subroutine get_eqcm_charge


  Subroutine eqcm_mass_calibration(eqcm_data, files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prints Mass-frequency (Delta f) vs charge variation (Delta Q) for 
    ! calibration purposes. Such exercise is usually conducted using the 
    ! electro-deposition of Ag on Pt, which is known to behave nicely.
    ! 
    ! author    - i.scivetti October 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),   Intent(InOut) :: eqcm_data
    Type(file_type),   Intent(InOut) :: files(:)

    Integer(Kind=wi)   :: i, j, k, l
    Integer(Kind=wi)   :: iunit_calibration
    Character(Len=256) :: messages(3)
    Character(Len=256) :: unit_calibration
 
    ! Open file
    Open(Newunit=files(FILE_CALIBRATION)%unit_no, File=files(FILE_CALIBRATION)%filename, Status='Replace')
    iunit_calibration=files(FILE_CALIBRATION)%unit_no
    unit_calibration =Trim(files(FILE_CALIBRATION)%filename)

    ! Print to OUT_EQCM
    Write (messages(1),'(a)') '  '
    Write (messages(2),'(1x,2a)') 'Print EQCM frequency to mass relation to file ',&
                                & Trim(FOLDER_ANALYSIS)//'/'//Trim(unit_calibration) 
    Call info(messages,2)

    ! Print to file
    Write (iunit_calibration,'(a,6x,a)')    '# Charge [mC]', 'Mass-frequency [Hz]'

    Do i = eqcm_data%range_cycles%value(1), eqcm_data%range_cycles%value(2)
      Do j= 1, 2
       Write (iunit_calibration,'(a,i3,2a)') '# Cycle ', i, '- Voltage sweep: ', Trim(eqcm_data%label_leg(j,i))
        l=eqcm_data%segment_points(j,i)
        Do k =1, l
           Write (iunit_calibration,'(f12.4,8x,f12.4)')  eqcm_data%charge%value(k,j,i), eqcm_data%mass_frequency%value(k,j,i)
        End Do
        If (l==0) Then
          Write (iunit_calibration,'(a)') '# Either no data or number of the EQCM points is not sufficient'
        End If
        Write (iunit_calibration,*) ' '
      End Do
    End Do

    ! Close file 
    Close(iunit_calibration)

    ! Move file to directory 
    Call execute_command_line('mv '//Trim(files(FILE_CALIBRATION)%filename)//' '//Trim(FOLDER_ANALYSIS))

  End Subroutine eqcm_mass_calibration


  Subroutine eqcm_massogram(eqcm_data, files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! computes and prints the massogram from the EQCM mass change 
    ! 
    ! author    - i.scivetti November 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(eqcm_type),   Intent(InOut) :: eqcm_data
    Type(file_type),   Intent(InOut) :: files(:)

    Integer(Kind=wi)   :: i, j, k, l, fail
    Integer(Kind=wi)   :: iunit_massogram
    Character(Len=256) :: unit_massogram
    Character(Len=256) :: messages(3)

    Real(Kind=wp)                  :: Dt, DM 
    Real(Kind=wp), Allocatable     :: massogram(:,:,:)  

    ! Allocate massogram array
    Allocate(massogram(eqcm_data%max_points, 2, eqcm_data%ncycles), Stat=fail)
    If (fail > 0) Then
      Write (messages(1),'(1x,1a)') '***ERROR: Allocation problems for massogram array'
      Call error_stop(messages(1))
    End If 

    massogram=0.0_wp

    Do i = eqcm_data%range_cycles%value(1), eqcm_data%range_cycles%value(2)
      Do j = 1, 2
        If (eqcm_data%segment_points(j,i) >=2) Then
          Do k = 1, eqcm_data%segment_points(j,i)
            ! Computes Dt
            If (k==1) Then
              If (eqcm_data%time%fread) Then
                Dt=Abs(eqcm_data%time%value(k+1,j,i)-eqcm_data%time%value(k,j,i))
              Else
                Dt=Abs(eqcm_data%voltage%value(k+1,j,i)-eqcm_data%voltage%value(k,j,i))/eqcm_data%scan_rate%value(1)
              End If
            ElseIf (k==eqcm_data%segment_points(j,i)) Then
              If (eqcm_data%time%fread) Then 
                Dt=Abs(eqcm_data%time%value(k,j,i)-eqcm_data%time%value(k-1,j,i))
              Else
                Dt=Abs(eqcm_data%voltage%value(k,j,i)-eqcm_data%voltage%value(k-1,j,i))/eqcm_data%scan_rate%value(1) 
              End If
            Else 
              If (eqcm_data%time%fread) Then 
                Dt=Abs(eqcm_data%time%value(k+1,j,i)-eqcm_data%time%value(k-1,j,i))
              Else
                Dt=Abs(eqcm_data%voltage%value(k+1,j,i)-eqcm_data%voltage%value(k-1,j,i))/eqcm_data%scan_rate%value(1)
              End If  
            End If 
            ! Computes DM
            If (eqcm_data%mass%fread) Then
              If (k==1) Then
                DM=eqcm_data%mass%value(k+1,j,i)-eqcm_data%mass%value(k,j,i)
              ElseIf (k==eqcm_data%segment_points(j,i)) Then
                DM=eqcm_data%mass%value(k,j,i)-eqcm_data%mass%value(k-1,j,i)
              Else 
                DM=eqcm_data%mass%value(k+1,j,i)-eqcm_data%mass%value(k-1,j,i)
              End If 
            Else        
              If (eqcm_data%mass_frequency%fread) Then
                If (k==1) Then
                  DM=eqcm_data%mass_frequency%value(k+1,j,i)-eqcm_data%mass_frequency%value(k,j,i)
                ElseIf (k==eqcm_data%segment_points(j,i)) Then
                  DM=eqcm_data%mass_frequency%value(k,j,i)-eqcm_data%mass_frequency%value(k-1,j,i)
                Else 
                  DM=eqcm_data%mass_frequency%value(k+1,j,i)-eqcm_data%mass_frequency%value(k-1,j,i)
                End If 
              End If 
            End If 
            ! Compute DM/dt
            massogram(k,j,i)=DM/Dt
          End Do
        End If
      End Do
    End Do

    ! Open file
    Open(Newunit=files(FILE_MASSOGRAM)%unit_no, File=files(FILE_MASSOGRAM)%filename, Status='Replace')
    iunit_massogram=files(FILE_MASSOGRAM)%unit_no
    unit_massogram =Trim(files(FILE_MASSOGRAM)%filename)

    Write (messages(1),'(a)') '  '
    Write (messages(2),'(1x,2a)') 'Print massogram to file ', Trim(FOLDER_ANALYSIS)//'/'//Trim(unit_massogram)
    Call info(messages,2)

    If (eqcm_data%mass%fread) Then
        Write (iunit_massogram,'(a,6x,a)')    '# Voltage [V]', 'Mass rate [ng/s]'
    Else
      If (eqcm_data%mass_frequency%fread) Then       
        Write (iunit_massogram,'(a,6x,a)')    '# Voltage [V]', 'Mass-density rate [ng/cm2/s]'
      End If        
    End If        

    Do i = eqcm_data%range_cycles%value(1), eqcm_data%range_cycles%value(2)
      Do j= 1, 2
       Write (iunit_massogram,'(a,i3,2a)') '# Cycle ', i, '- Voltage sweep: ', Trim(eqcm_data%label_leg(j,i))
        l=eqcm_data%segment_points(j,i)
        If (l>=2) Then
          Do k =1, l
            Write (iunit_massogram,'(f12.4,8x,f12.4)')  eqcm_data%voltage%value(k,j,i), massogram(k,j,i)
          End Do
        Else
          Write (iunit_massogram,'(a)') '# Either no data or number of the EQCM points is not sufficient'
        End If
        Write (iunit_massogram,*) ' '
      End Do
    End Do

    ! Close file
    Close(iunit_massogram)
     
    ! Deallocate
    Deallocate(massogram)

    ! Move file to directory
    Call execute_command_line('mv '//Trim(files(FILE_MASSOGRAM)%filename)//' '//Trim(FOLDER_ANALYSIS))

  End Subroutine eqcm_massogram

End Module eqcm
