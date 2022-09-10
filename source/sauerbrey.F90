!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for Sauerbrey correlation
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author: i.scivetti  July 2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module sauerbrey 
   
  Use constants,     Only : code_VERSION
  Use eqcm,          Only : eqcm_type
  Use fileset,       Only : file_type, &
                            FILE_SET_EQCM
  Use input_types,   Only : in_param_array
  Use numprec,       Only : wi, &
                            wp 
  Use system,        Only : system_type
  Use unit_output,   Only : error_stop,&
                            info 
                    
  Implicit None
  Private

  ! Porcentage limit (2%) for the validity of the Sauerbrey approximation
  Real(Kind=wp), Parameter  :: valid_limit = 0.02

  ! Sauerbrey type
  Type, Public :: sauer_type
    Private
    ! Factor to transform mass-frequency to mass/area 
    Type(in_param_array),  Public               :: factor
    ! Maximum delta mass-frequency change
    Real(Kind=wp), Public                 :: max_Dmassfreq
  Contains
      Private
      Procedure, Public :: init_input_variables  => allocate_sauer_input_variables
      Final             :: cleanup
  End Type sauer_type

  Public ::  sauerbrey_correlation, sauerbrey_transformation

Contains

   Subroutine allocate_sauer_input_variables(T)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Allocate essential sauer input variables to be read from SET_EQCM
     !
     ! author    - i.scivetti Nov 2020
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Class(sauer_type), Intent(InOut)  :: T

     Integer(Kind=wi)     :: fail(2)
     Character(Len=256)   :: message

     Allocate(T%factor%value(1),         Stat=fail(1))
     Allocate(T%factor%units(3),         Stat=fail(2))

     If (Any(fail > 0)) Then
       Write (message,'(1x,1a)') '***ERROR: Allocation problems for essential Sauer input variables'
       Call error_stop(message)
     End If

     T%factor%value(1)=0.0_wp

   End Subroutine allocate_sauer_input_variables

   Subroutine cleanup(T)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Deallocation fore the Final procedure
   ! 
   ! author    - i.scivetti Nov 2020
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Type(sauer_type) :: T

     If (Allocated(T%factor%value)) Then
       Deallocate (T%factor%value)
     End If

     If (Allocated(T%factor%units)) Then
       Deallocate (T%factor%units)
     End If

   End Subroutine cleanup

   Subroutine sauerbrey_correlation(files, eqcm_data, system_data, sauer_data)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Uses the mass-frequency from input data, finds the maximum value and 
     ! checks the validity of the Sauerbrey equation 
     ! 
     ! author    - i.scivetti June 2020
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Type(file_type),     Intent(InOut) :: files(:)
     Type(eqcm_type),     Intent(In   ) :: eqcm_data
     Type(system_type),   Intent(In   ) :: system_data     
     Type(sauer_type),    Intent(InOut) :: sauer_data     

     Real(Kind=wp)    :: ratio 
     Character(Len = 256)  :: message, messages(2)

     ! Test sauerbrey condition
     If (.Not. (sauer_data%factor%fread) ) Then
       Call info(' ', 1)
       Write (message,'(1x,3a)') '***ERROR - Computation of the Sauerbrey factor from electrode properties is not&
                               & implemented. The user must define directive "sauerbrey" in file ',&
                               & Trim(files(FILE_SET_EQCM)%filename), ', in units of [Hz ng-1 cm2].'
       Call error_stop(message)
     End If

     ! Test sauerbrey condition
     sauer_data%max_Dmassfreq=maxval(Abs(eqcm_data%mass_frequency%value))
     
     ratio = 100.0_wp * sauer_data%max_Dmassfreq/system_data%quartz_freq%value

     If (ratio > valid_limit) Then
       Write (messages(1),'(1x,a)') '***ERROR - Maximum variation for mass-frequency is larger than 2% of the fundametal&
                              & frequency of the quartz crystal. This does not comply with the condition for the&  
                              & Sauerbrey correlation.'
       Write (messages(2),'(1x,a)') '   Make sure that the recorded data is be indeed analysed with the Sauerbrey model.& 
                      & Otherwise, input mass must be provided by the user. Check also the input for the "quartz_freq" directive.'
       Call info(messages, 2)               
       Call error_stop(' ')
     End If

   End Subroutine sauerbrey_correlation

   Subroutine sauerbrey_transformation(eqcm_data, sauer_data)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Transform mass-frequency to mass-density using the Sauerbrey equation
     ! 
     ! author    - i.scivetti June 2020
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Type(eqcm_type),      Intent(InOut) :: eqcm_data
     Type(sauer_type),     Intent(InOut) :: sauer_data     

     Integer(Kind=wi)      :: i, j, k
     Character(Len = 256)  :: message

     If (Abs(sauer_data%factor%value(1)) <= epsilon(sauer_data%factor%value(1))) Then
       Write (message,'(1x,a)') '***ERROR - Sauerbrey factor has not been specified'
       Call error_stop(message)
     End If 
 
     Do i = 1, eqcm_data%ncycles
       Do j= 1, 2
         Do k = 1, eqcm_data%segment_points(j,i)
            eqcm_data%mass_frequency%value(k,j,i) = eqcm_data%mass_frequency%value(k,j,i)/sauer_data%factor%value(1)
         End Do
       End Do
     End Do

   End Subroutine sauerbrey_transformation

End Module sauerbrey
