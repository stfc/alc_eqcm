!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for definition of variables related to the 
! characteristics of the EQCM device and surroundings 
!
! Copyright: 2022-2024 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author: i.scivetti April 2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module system 
  
  Use input_types,  Only : in_param
  
  Implicit None
  Private

  ! Type to describe the EQCM device/system
  Type, Public :: system_type
    Private
    ! Quartz resonant frequency
    Type(in_param),   Public       :: quartz_freq
  End Type system_type

End Module system
