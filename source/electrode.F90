!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for definition of electrode properties and quantities 
!
! Copyright: 2022-2024 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author: i.scivetti June 2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module electrode

  Use input_types,  Only : in_param, &
                           in_scalar

  Implicit None

  Type, Public :: electrode_type
    Private
    ! Total mass of EQCM electrode
    Type(in_param), Public                :: mass
    ! Number of EQCM electrode
    Type(in_scalar),  Public              :: mole
    ! Electrode area
    Type(in_param),   Public              :: area_geom
    ! Scaling factor 
    Type(in_scalar),  Public              :: area_scale

  End Type electrode_type
 
End Module electrode
