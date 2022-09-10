!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for Fast Fourier Transform (FFT)
!
! Copyright - 2022 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author: i.scivetti June 2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module fft 

  Use constants,    Only : pi
  Use fileset,      Only : file_type
  Use input_types,  Only : in_integer
  Use numprec,      Only : wi, &
                           wp 
  Use unit_output,  Only : error_stop, &
                           info
 
  Implicit None

  Private

  ! Maximum points for fft transformation
  Integer(Kind=wi), Public, Parameter     :: fftnmax = 32768
  ! Numerical criteria to neglect contribution of the imaginary part
  Real(Kind=wp), Public, Parameter        :: fft_imag_criteria = 1.0E-10

  Type, Public :: fft_type
    Private
    ! Points for fft transformation
    Integer(Kind=wi), Public              :: ntot 
    ! Indexes for the limits of measured data within the fft domain
    Integer(Kind=wi), Public              :: nsides(2)
    ! Discretization in the reciprocal space 
    Real(Kind=wp), Public                 :: Drec 
    ! Reciprocal space domain 
    Real(kind=wp),    Public, Allocatable :: rec_domain(:) 
    ! FFT array
    Complex(kind=wp), Public, Allocatable :: array(:)
    ! magnitude of FFT elements
    Real(kind=wp),    Public, Allocatable :: magnitude(:) 
    ! number of measured points to be averaged for padding (current)
    Type(in_integer), Public               :: end_current  
    ! number of measured points to be averaged for padding (mass/mass-density)
    Type(in_integer), Public               :: end_mass  
    ! Largest value of endpoints between endpoints_mass and endpoints_current
    Integer(Kind=wi), Public               :: endpoints=0
  End Type fft_type

  Public :: direct_fft, inverse_fft
  Public :: fft_setup, fft_points, fft_spectrum, print_spectrum

Contains

  Subroutine fft_points(fft_data, npoints)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define the number of points for FFT based on the amount of EQCM data
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fft_type),    Intent(InOut)  :: fft_data
    Integer(Kind=wi),           Intent(In   )  :: npoints  

    Integer(Kind=wi)   :: i
    Logical            :: cont
    Character(Len=256) :: message

    ! Find the number of FFT points
    i=0
    cont = .True.
    fft_data%ntot = 1

    Do While (cont)
      If (fft_data%ntot > fftnmax) Then
        Write (message,'(1x,a,i4)') '***ERROR - Number of required points for FFT exceeds the maximum number of ', fftnmax
        Call error_stop(message)
      End If
      If (npoints < fft_data%ntot) Then
        cont = .False.
      Else
        i=i+1
        fft_data%ntot=2**i
      End IF
    End Do

  End Subroutine fft_points


  Subroutine fft_setup(fft_data, npoints, array, domain, points_extreme)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define fft grid and padding
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fft_type),    Intent(InOut)  :: fft_data
    Integer(Kind=wi),  Intent(In   )  :: npoints  
    Real(Kind=wp),     Intent(In   )  :: array(npoints)
    Real(Kind=wp),     Intent(In   )  :: domain
    Integer(Kind=wi),  Intent(In   )  :: points_extreme

    Integer(Kind=wi)  :: i, j
    Character(Len=256)  :: message

    Real(Kind=wp)    :: rarray
    Real(Kind=wp)    :: pad_avg(2)
    Integer(Kind=wi) :: npad(2)
    Integer(Kind=wi) :: nrange(2)


    ! Assing number of points for padding
    If (mod(npoints,2)==0) Then
      npad(1)=(fft_data%ntot-npoints)/2
      npad(2)=npad(1)
    Else
      npad(1)=(fft_data%ntot-npoints-1)/2+1
      npad(2)=npad(1)-1
    End If

    If ( (npad(1)+npad(2)+npoints) /= fft_data%ntot) Then
      Write (message,'(1x,a)') '***ERROR -  Wrong assignment of padding for FFT'
      Call error_stop(message)
    End If

    fft_data%nsides(1)=npad(1)+1
    fft_data%nsides(2)=fft_data%nsides(1)+npoints-1

    pad_avg=0.0_wp

    Do j = 1, 2
      If (j==1) Then
        nrange(1)= 1
        nrange(2)= points_extreme
      Else
        nrange(1)=npoints-points_extreme+1
        nrange(2)=npoints
      End If

      Do i= nrange(1), nrange(2)
         pad_avg(j)=pad_avg(j)+array(i)
      End Do
      pad_avg(j)=pad_avg(j)/points_extreme
    End Do

    Do i = 1, fft_data%ntot
      If (i<fft_data%nsides(1)) Then
        rarray=pad_avg(1)
      ElseIf (i >= fft_data%nsides(1) .And. i <= fft_data%nsides(2)) Then
        rarray=array(i-fft_data%nsides(1)+1)
      ElseIf (i>fft_data%nsides(2)) Then
        rarray=pad_avg(2)
      End If
      fft_data%array(i)=CMPLX(rarray, 0.0_wp, kind=wp)
    End Do

    ! Definition of discretization in the reciprocal space, either the space of frequencies (1/s)  or inverse of voltage (1/V)
    fft_data%Drec=Real(npoints-1,Kind=wp)/(fft_data%ntot-1)/domain

    ! Define grid in the reciprocal space
    Do i =1, fft_data%ntot/2+1
      fft_data%rec_domain(i)=Real(i-1,Kind=wp)*fft_data%Drec
    End Do

  End Subroutine fft_setup

  Subroutine fft_spectrum(fft_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute magnitude of the FFT components
    ! 
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(fft_type),    Intent(InOut)  :: fft_data

    Integer(Kind=wi) :: i

    Call direct_fft(fft_data%array)

    Do i= 1, fft_data%ntot/2+1
      fft_data%magnitude(i)=Abs(fft_data%array(i))
    End Do

  End Subroutine fft_spectrum

  Subroutine print_spectrum(iunit,fft_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Print spectrum in reciprocal space
    !
    ! author    - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   )  :: iunit        
    Type(fft_type),    Intent(InOut)  :: fft_data

    Integer(Kind=wi) :: i

    Do i= 1, fft_data%ntot/2+1
      Write (iunit,'(f12.4,16x,f12.4)')  fft_data%rec_domain(i), fft_data%magnitude(i)
    End Do
    Write (iunit,*) ' '  

  End Subroutine print_spectrum

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutines for direct and inverse FFT based on the Colley-Tukey 
  ! recursive algorithm
  ! 
  ! author    - i.scivetti June 2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Subroutine direct_fft(x)
    Complex(kind=wp), Dimension(:), intent(inout)  :: x

    Call apply_fft(x,1.0_wp)

  End Subroutine direct_fft

  Subroutine inverse_fft(x,N)
    Complex(kind=wp), Dimension(:), Intent(inout)  :: x
    Integer(Kind=wi),               Intent(InOut)  :: N

    Call apply_fft(x,-1.0_wp)
    x=x/N

  End Subroutine inverse_fft

  Recursive Subroutine apply_fft(x,phase)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to execute the Cooley-Tukey FFT
    !
    ! author   - i.scivetti June 2020 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Complex(kind=wp), Dimension(:), Intent(InOut)  :: x
    Real(kind=wp),                  Intent(In   )  :: phase 

    Complex(kind=wp), Dimension(:), Allocatable    :: even, odd

    Complex(kind=wp) :: t
    Integer(Kind=wi) :: N, i
    Integer(Kind=wi) :: fail(2)

    N=size(x)

    If (N <= 1) return
    
    ! Allocate
    Allocate(odd((N+1)/2), Stat=fail(1))
    Allocate(even(N/2),    Stat=fail(2)) 
    If (Any(fail > 0)) Then
      Call error_stop(' ***ERROR: Allocation problems for arrays odd and even in subroutine "apply_fft"')
    End If

    ! divide
    odd =x(1:N:2)
    even=x(2:N:2)

    ! conquer
    Call apply_fft(odd,phase)
    Call apply_fft(even,phase)

    ! combine divide and conquer
    Do i=1,N/2
       t=exp(phase*cmplx(0.0_wp,-2.0_wp*pi*real(i-1,wp)/real(N,wp),kind=wp))*even(i)
       x(i)     = odd(i) + t
       x(i+N/2) = odd(i) - t
    End Do

    ! Deallocate
    Deallocate(odd)
    Deallocate(even)

  end subroutine apply_fft

End Module fft
