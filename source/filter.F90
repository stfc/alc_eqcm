!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Module to filter EQCM data via a convolution
! in reciprocal space with with a Gaussian low pass filter
!
! Copyright: 2022-2024 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! author: i.scivetti July 2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module filtering 

  Use constants,    Only : pi
  Use fft,          Only : fft_imag_criteria,&
                           fft_type, & 
                           direct_fft,&
                           inverse_fft
  Use fileset,      Only : file_type
  Use input_types,  Only : in_param
  Use numprec,      Only : wi,&
                           wp 
  Use unit_output,  Only : error_stop, &
                           info 
  Implicit None

  Private

  ! EQCM types of analysis for Electrochemical Characterization
  Type, Public :: filter_type
    Private
    ! Cutoff value for the reciprocal space 
    Type(in_param),   Public  :: cutoff
  End Type filter_type

  Public :: filter_data 

Contains

  Subroutine filter_data(fft_data, filter, npoints, array)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to remove high frequency components of EQCM values using
    ! a low-pass Gaussian filter
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fft_type),    Intent(InOut)  :: fft_data
    Type(filter_type), Intent(InOut)  :: filter
    Integer(Kind=wi),  Intent(In   )  :: npoints
    Real(Kind=wp),     Intent(InOut)  :: array(npoints)     
    
    Character(Len = 256)  :: messages(3) 
    Integer(Kind=wi)      :: i

    Call direct_fft(fft_data%array)

    Call low_pass_gaussian(fft_data, filter)

    Call inverse_fft(fft_data%array,fft_data%ntot)

    Do i= fft_data%nsides(1), fft_data%nsides(2)
      If (Abs(Aimag(fft_data%array(i))) > fft_imag_criteria) Then
        Write (messages(1),'(1x,a)') '***ERROR - Absolute value of the imaginary part of an element following'
        Write (messages(2),'(1x,a)') '***an inverse FFT is greater than the numerical criteria.'
        Write (messages(3),'(1x,a)') '***If this error manifests something is very wrong with the input EQCM data.'    
        Call info(messages,3)
        Call error_stop(' ')
      Else
        array(i-fft_data%nsides(1)+1)=Real(fft_data%array(i))  
      End If
    End Do

  End Subroutine filter_data

  Subroutine low_pass_gaussian(fft_data,filter_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Applies filtering via convolution
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fft_type),      Intent(InOut)  :: fft_data
    Type(filter_type),   Intent(In   )  :: filter_data

    Integer(Kind=wi) ::  i
    Real(kind=wp)   ::  sigma
    Real(kind=wp)   ::  fact(fft_data%ntot)

    ! Half the Full Width at Half Maximum (FWHM)
    sigma=filter_data%cutoff%value/sqrt(2.0_wp*log(2.0_wp))  

    Do i =1, fft_data%ntot/2+1
      fact(i)=fft_data%rec_domain(i)/sigma
      fact(i)=fact(i)**2
      fact(i)=exp(-fact(i))
    End Do

    Do i =2, fft_data%ntot/2
      fact(fft_data%ntot+2-i)=fact(i)
    End Do

    Do i=1, fft_data%ntot
      fft_data%array(i)=fact(i)*fft_data%array(i)
    End Do 

  End Subroutine low_pass_gaussian


End Module filtering
