!*************************************************************************
Module AtomPotentialMod
!*************************************************************************
!
! ***  Function atomPotential(q2,elementLabel) is contained here - It is
! ***  the only part of the plane wave program that needs to be changed 
! ***  for different potentials 
!

contains

!*************************************************************************
  Function atomPotential(q2,elementLabel)
!*************************************************************************
!
!   pseudopotentials for Si, Ga, and As taken from Zhang, Yeh, and Zunger
!   Phys. Rev. B 48 11204 (1993).

    use SysParams,	only : two, zero, double

    implicit none

    character*2, intent(in)  :: elementLabel  ! atomic species
    real(double), intent(in) :: q2            ! squared reciprocal q-value
    real(double)             :: atomPotential     ! pseudopotential V(q)
    real(double)             :: rs
 
    select case (elementLabel)

    case("H")    
       rs = 1.d0  !density parameter in a0, volume = (4pi/3)rs^3
                ! V(q2) = (4pi/vol)/(q2 + k0^2) = (3/rs^3)/(q2 + 2.43/rs)
       atomPotential = - (3.d0/rs**3)/(q2 + 2.43/rs)     !screened potential
       !atomPotential = 12.56637061/(q2)/(1.0+(2.4336/q2))

    case("Si")
      atomPotential = 0.53706*(q2-2.19104)/(2.05716*exp(0.48716*q2)-1)/two

    case("Ga")
      atomPotential = 1.22*(q2-2.45)/(exp(0.54*(q2+2.71))+1)/two

    case("As")
      atomPotential = 0.35*(q2-2.62)/(exp(0.93*(q2-1.57))+1)/two

    case("El")
      atomPotential = zero ! free-electron case

    case default

       stop '***ERROR*** pseudopotential not coded for this element'

    end select

  end Function atomPotential



end Module AtomPotentialMOd
