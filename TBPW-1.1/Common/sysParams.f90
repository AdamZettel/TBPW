!*************************************************************************
Module SysParams
!*************************************************************************
  implicit none

  save

  integer, parameter         :: double = kind(0.0D0) ! double precision
  real(double), parameter    :: zero  = 0D0     ! whole numbers
  real(double), parameter    :: one   = 1D0      
  real(double), parameter    :: two   = 2D0
  real(double), parameter    :: four  = 4D0
  real(double), parameter    :: half  = 0.50D0

  ! define some mathematical & physical constants

  real(double), parameter    :: epsilon = 1.0d-30
  real(double), parameter    :: twopi = 6.283185307179586476925287D0
  real(double), parameter    :: log2 = 0.301030103010301030103010D0
  real(double), parameter    :: ryd2eV   = 13.60569805D0
  real(double), parameter    :: hartree2eV = two * ryd2eV    ! atomic unit of energy
  real(double), parameter    :: hBar2 = 1.0d0          ! atomic unit
  real(double), parameter    :: electronMass = 1.0d0   ! atomic unit
  complex(double), parameter :: imaginary = (0.0d0, 1.0d0)
  integer, parameter         :: maxDimensions = 3

end Module SysParams
