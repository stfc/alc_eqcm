!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Module for input/output files and related subroutines
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author -    i.scivetti  April 2020
!!!!!!!!!!!!!!!!!!!!11!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module fileset

  Use constants,    Only: code_name,&
                          code_VERSION, &
                          date_RELEASE
  Use numprec,      Only: wi
  use unit_output,  Only: info, &
                          set_output_unit 

  Implicit None
  Private

  ! File data
  Type, Public :: file_type
    Private
    ! Filename
    Character(Len=256), Public :: filename
    ! Fortran unit number, set with newunit=T%unit_no
    Integer(Kind=wi), Public   :: unit_no = -2
  Contains
    Procedure, Public :: init => file_type_init
    Procedure, Public :: Close => close_file
  End Type file_type

  ! SET_EQCM file
  Integer(Kind=wi), Parameter, Public :: FILE_SET_EQCM   = 1
  ! OUT_EQC file
  Integer(Kind=wi), Parameter, Public :: FILE_OUT_EQCM   = 2 
  ! DATA_EQCM file
  Integer(Kind=wi), Parameter, Public :: FILE_DATA_EQCM  = 3
  ! Filtered data for current
  Integer(Kind=wi), Parameter, Public :: FILE_FILTERED_CURRENT = 4          
  ! Filtered data for mass 
  Integer(Kind=wi), Parameter, Public :: FILE_FILTERED_MASS = 5          
  ! Raw data for current
  Integer(Kind=wi), Parameter, Public :: FILE_RAW_CURRENT = 6          
  ! Raw data for mass 
  Integer(Kind=wi), Parameter, Public :: FILE_RAW_MASS = 7          
  ! CHARACTERIZATION 
  Integer(Kind=wi), Parameter, Public :: FILE_CHARACT = 8  
  ! CURRENT SPECTRUM file
  Integer(Kind=wi), Parameter, Public :: FILE_SPEC_CURRENT  = 9
  ! DELTA FREQUENCY SPECTRUM file
  Integer(Kind=wi), Parameter, Public :: FILE_SPEC_MASS  = 10
  ! CALIBRATION file
  Integer(Kind=wi), Parameter, Public :: FILE_CALIBRATION = 11 
  ! MASSOGRAM file
  Integer(Kind=wi), Parameter, Public :: FILE_MASSOGRAM = 12 
  ! ELECTRO_DEPOSITION file
  Integer(Kind=wi), Parameter, Public :: FILE_ELECTRO = 13 
  ! INTERCALATION for OXIDATION file
  Integer(Kind=wi), Parameter, Public :: FILE_INTERCALATION_OX  = 14 
  ! INTERCALATION for REDUCTION file
  Integer(Kind=wi), Parameter, Public :: FILE_INTERCALATION_RED = 15
  ! INPUT STRUCTURE
  Integer(Kind=wi), Parameter, Public :: FILE_INPUT_STRUCTURE = 16
  ! OUTPUT VASP STRUCTURE 
  Integer(Kind=wi), Parameter, Public :: FILE_OUTPUT_STRUCTURE = 17
  ! SELECTED SOLUTIONS file
  Integer(Kind=wi), Parameter, Public :: FILE_SELECTED_SOLUTIONS  = 18 
  ! SET_SIMULATION file
  Integer(Kind=wi), Parameter, Public :: FILE_SET_SIMULATION  = 19 
  ! KPOINTS file
  Integer(Kind=wi), Parameter, Public :: FILE_KPOINTS  = 20 
  ! HPC_SETTINGS file
  Integer(Kind=wi), Parameter, Public :: FILE_HPC_SETTINGS  = 21 
  ! RECORD_MODELS file 
  Integer(Kind=wi), Parameter, Public :: FILE_RECORD_MODELS = 22
  ! MODEL_SUMMARY file 
  Integer(Kind=wi), Parameter, Public :: FILE_MODEL_SUMMARY = 23

  ! Size of filename array
  Integer(Kind=wi), Parameter, Public :: NUM_FILES = 23

  ! Folder data
  Character(Len=256), Public :: FOLDER_MODELS     = "ATOMISTIC_MODELS"
  Character(Len=256), Public :: FOLDER_DFT        = "DFT"
  Character(Len=256), Public :: FOLDER_RESTART    = "RESTART"
  Character(Len=256), Public :: FOLDER_INPUT_GEOM = "INPUT_GEOM"
  Character(Len=256), Public :: FOLDER_ANALYSIS   = "ANALYSIS_EQCM"

  Public :: set_system_files, print_header_out_eqcm, wrapping_up, refresh_out_eqcm

Contains

  Subroutine refresh_out_eqcm(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to refresh the output
    !
    ! author    - i.scivetti Oct 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES) 

    Call files(FILE_OUT_EQCM)%close ()
    Open (Newunit=files(FILE_OUT_EQCM)%unit_no, File=files(FILE_OUT_EQCM)%filename, Position='Append')
    Call set_output_unit(files(FILE_OUT_EQCM)%unit_no)

  End Subroutine refresh_out_eqcm

  Subroutine file_type_init(T, filename)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to initialise files
    !
    ! author    - i.scivetti April 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    Class(file_type)                :: T
    Character(Len=*), Intent(In   ) :: filename

    T%filename = Trim(filename)
  End Subroutine file_type_init


  Subroutine set_names_files(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set default names for files
    !
    ! author    - i.scivetti April 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Character(Len=256), Dimension(NUM_FILES)   :: set_names
    Integer(Kind=wi)                           :: file_no

    ! Default file names array
    ! Populate default names array
    set_names(FILE_SET_EQCM)         = "SET_EQCM"
    set_names(FILE_OUT_EQCM)         = "OUT_EQCM"
    set_names(FILE_DATA_EQCM)        = "DATA_EQCM"
    set_names(FILE_RAW_CURRENT)      = "RAW_CURRENT"
    set_names(FILE_RAW_MASS)         = "RAW_MASS"
    set_names(FILE_CHARACT)          = "CHARACTERIZATION"
    set_names(FILE_FILTERED_MASS)    = "FILTERED_MASS"
    set_names(FILE_FILTERED_CURRENT) = "FILTERED_CURRENT"
    set_names(FILE_SPEC_CURRENT)     = "SPEC_CURRENT"
    set_names(FILE_SPEC_MASS)        = "SPEC_MASS"
    set_names(FILE_CALIBRATION)      = "MASS_CALIBRATION"
    set_names(FILE_MASSOGRAM)        = "MASSOGRAM"
    set_names(FILE_ELECTRO)          = "ELECTRO_DEPOSITION"
    set_names(FILE_INTERCALATION_OX) = "INTERCALATION_OX"
    set_names(FILE_INTERCALATION_RED)= "INTERCALATION_RED"
    set_names(FILE_INPUT_STRUCTURE)  = "INPUT_STRUCTURE"
    set_names(FILE_OUTPUT_STRUCTURE) = "OUTPUT_STRUCTURE"
    set_names(FILE_SELECTED_SOLUTIONS) = "SELECTED_SOLUTIONS"
    set_names(FILE_SET_SIMULATION)   = "SET_SIMULATION"
    set_names(FILE_KPOINTS)          = "SET_KPOINTS"
    set_names(FILE_HPC_SETTINGS)     = "HPC_SETTINGS"
    set_names(FILE_RECORD_MODELS)    = "RECORD_MODELS"
    set_names(FILE_MODEL_SUMMARY)    = "MODEL_SUMMARY"

    ! Set default filenames
    Do file_no = 1, NUM_FILES
      Call files(file_no)%init(set_names(file_no))
    End Do
  End Subroutine set_names_files


  Subroutine close_file(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to close files
    !
    ! author    - i.scivetti April 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    Class(file_type) :: T

    Logical :: is_open

    Inquire (T%unit_no, opened=is_open)
    If (is_open) Then
      Close (T%unit_no)
      T%unit_no = -2
    End If

  End Subroutine close_file

  Subroutine set_system_files(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to open OUTPUT file 
    ! 
    ! author    - i.scivetti April 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Call set_names_files(files)   
    Open (Newunit=files(FILE_OUT_EQCM)%unit_no, File=files(FILE_OUT_EQCM)%filename, Status='replace')
    Call set_output_unit(files(FILE_OUT_EQCM)%unit_no)

  End Subroutine set_system_files   

  Subroutine print_header_out_eqcm(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the header to OUTPUT file 
    !  
    ! author        - i.scivetti July 2020
    ! contribution  - i.scivetti Oct  2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Character(Len=*), Parameter :: fmt1 = '(a)'
    Character(Len=*), Parameter :: fmt2 = '(3a)'
    Character(Len=*), Parameter :: fmt3 = '(4a)'
    Character(Len=128)          :: header(18)

    Write (header(1), fmt1)   Repeat("#", 74)
    Write (header(2), fmt2)  "#                      WELCOME TO ", Trim(code_name),  Repeat(" ", 31)//"#"
    Write (header(3), fmt1)  "#  SCD/ALC program for Electrochemical analysis and atomistic modeling   #"
    Write (header(4), fmt1)  "#  using Electrochemical Quartz Crystal Microbalance (EQCM) experiments  #"
    Write (header(5), fmt3)  "#  version:  ", Trim(code_VERSION), Repeat(' ',57),                     "#"
    Write (header(6), fmt3)  "#  release:  ", Trim(date_RELEASE), Repeat(' ',52),                     "#"
    Write (header(7), fmt1)  "#                                                                        #"
    Write (header(8), fmt1)  "#  Copyright:  2022  Ada Lovelace Centre (ALC)                           #"
    Write (header(9), fmt1)  "#              Scientific Computing Department (SCD)                     #"
    Write (header(10), fmt1) "#              Science and Technology Facilities Councils (STFC)         #"
    Write (header(11), fmt1) "#                                                                        #"
    Write (header(12), fmt1) "#  Author:                   Ivan Scivetti (SCD-STFC)                    #"
    Write (header(13), fmt1) "#  Scientific support:       Gilberto Teobaldi (SCD-STFC)                #"
    Write (header(14), fmt1) "#  Experimental partners:    Sofia Diaz-Moreno (Diamond-Spectroscopy)    #"
    Write (header(15), fmt1) "#                            Paul Donaldson    (CLF-ULTRA)               #" 
    Write (header(16), fmt1) "#                            Daniel Bowron (ISIS-Disordered Materials)   #"
    Write (header(17), fmt1) "#                                                                        #"
    Write (header(18), fmt1)  Repeat("#", 74)
    Call info(header, 18)
    
    ! Refresh OUT_EQCM
    Call refresh_out_eqcm(files)

  End Subroutine print_header_out_eqcm

  Subroutine wrapping_up(files)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to print final remarks to OUT_EQCM file 
  ! and close the file 
  !  
  ! author    - i.scivetti July 2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)
 
    Character(Len=*), Parameter :: fmt1 = '(1x,a)'
    Character(Len=*), Parameter :: fmt2 = '(1x,3a)'
    Character(Len=128)          :: appex(7)
     
    Write (appex(1), fmt1)   Repeat(" ", 1)
    Write (appex(2), fmt1)   Repeat("#", 35)
    Write (appex(3), fmt1)  "#                                 #" 
    Write (appex(4), fmt1)  "#  Job has finished successfully  #"
    Write (appex(5), fmt2)  "#  Thanks for using ", Trim(code_name), " !!!  #"
    Write (appex(6), fmt1)  "#                                 #" 
    Write (appex(7), fmt1)   Repeat("#", 35)
    Call info(appex, 7)

    Close(files(FILE_OUT_EQCM)%unit_no)    

  End Subroutine wrapping_up  

End Module fileset
