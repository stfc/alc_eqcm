!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of common types   
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
! Copyright - Ada Lovelace Centre
!
! Author: i.scivetti  Dec 2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module input_types
  Use numprec,    Only: wi,&
                        wp

  Implicit None
  Private

  Type, Public :: in_string
    Character(Len=256) :: type=repeat(' ',256) 
    Logical            :: fread= .False.
    Logical            :: fail = .False.
    Logical            :: warn = .False. 
  End Type

  Type, Public :: in_integer
    Integer(Kind=wi)  :: value
    Character(Len=32) :: tag
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Logical           :: warn = .False. 
  End Type

  Type, Public :: in_integer_array
    Integer(Kind=wi), Allocatable :: value(:)
    Character(Len=32) :: tag
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Logical           :: warn = .False. 
  End Type

  Type, Public :: in_logic
    Logical           :: stat = .False. 
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Logical           :: warn = .False. 
  End Type

  Type, Public :: in_scalar
    Real(Kind=wp)     :: value
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Logical           :: warn = .False. 
  End Type

  Type, Public :: in_param
    Real(Kind=wp)     :: value
    Character(Len=16) :: units
    Real(Kind=wp)     :: convert
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Logical           :: warn = .False. 
  End Type

  Type, Public :: in_param_array
    Real(Kind=wp),     Allocatable :: value(:)
    Character(Len=16), Allocatable :: units(:)
    Real(Kind=wp)     :: convert
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Logical           :: warn = .False. 
  End Type

End Module input_types
