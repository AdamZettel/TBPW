!*************************************************************************
Module pwHamMod
!*************************************************************************
!
! The generic name of the module "SpecificHamMod" and the subroutines
! "HGenerate", "HamInfoInit", "HamInfoDestroy",  must be used so they are 
! called properly by the general hamiltonian routines
!
! This module contains all the routines needs for seeting up the
! Plane wave hamiltonian h (a 1d array which is upper triangular part
! of the H matrix of size = hSize)  Size is fixed by energy cutoff
! read from input file in the GVectors module
!
! ***  Function atomPotential(q2,elementLabel) is called here - It is
! ***  in module AtomPotentialMod - it is the only subroutine
! ***  needs to be changed for different potentials 
!
! This module uses the following modules:
!    SysParams
!    GVectorsMod
!    Hamiltonian
!    kPointsMod
!    LatticeMod  (only for bLatMetricTensor used to find lengths of k + G)
!
! The module contains the public functions:
!	HGenerate  Generates the plane wave hamiltonian for
!               the given k point K_Points$kPointsLat(:,k).
!	HamInfoInit  Sets up relevent information for the hamiltonian
!	HamInfoDestroy
!
! The module contains the private functions:
!	PotentialInit  allocates space for initializes potentialV.
!   PotentialDestroy deallocates the space for potentialV
!
! The module contains the private variables:
!   potentialV(i) Fourier component of the crystal (pseudo)potential for
!                   G-vector gVectors(i,:)
!*************************************************************************

  use SysParams,	only : double
  use GVectorsMod,      only : gVectorObjectsT
  Use EigenStatesMod,   only : eigenStatesT
  use typeMod,          only : GridT

  implicit none

  save

  type hamInfoT
    integer                    :: hMaxSize
    complex(double),   pointer :: potentialV(:)
    type( gVectorObjectsT ), pointer :: gVectors
  end type hamInfoT
  
  private PotentialInit
  private PotentialDestroy

  contains

!*************************************************************************
  Subroutine HamInfoInit( hamInfo, struct, tagHandler )
!*************************************************************************

    use GVectorsMod,   only : GVectorsInit
    use StructureMod,  only : structureT
    use TagHandlerMod, only : tagHandlerT

    implicit none

    type( hamInfoT ),    pointer :: hamInfo
    type( structureT ),  pointer :: struct
    type( tagHandlerT ), pointer :: tagHandler

    integer :: error

    if( associated( hamInfo )) stop 'Error: SpecificHamiltonian already allocated.'
    allocate( hamInfo, stat = error )
    if( error /= 0 ) stop 'Error: SpecificHamiltonian allocation failed.'

    call GVectorsInit( hamInfo%gVectors, struct, tagHandler )
    call PotentialInit( hamInfo, struct )

  end Subroutine HamInfoInit

!*************************************************************************
  Subroutine HamInfoDestroy( hamInfo )
!*************************************************************************

    use GVectorsMod, only : GVectorsDestroy
    implicit none

    type( hamInfoT ), pointer :: hamInfo

    integer :: error

    if( .not. associated( hamInfo )) stop 'Error: SpecificHamiltonian not allocated.'

    call GVectorsDestroy( hamInfo%gVectors )
    call PotentialDestroy( hamInfo )

    deallocate( hamInfo, stat = error )
    if( error /= 0 ) stop 'Error: specificHamiltonian deallocation failed.'

  end Subroutine HamInfoDestroy

!*************************************************************************
  Subroutine PotentialInit( hamInfo, struct )
!*************************************************************************
!  For all local potentials the Hamiltonian matrix depends on k only through
!  the diagonal terms. The Fourier components of the potential terms 
!  are independnet of k and can be set up only once for all k.

    use SysParams,    only : two
    use StructureMod, only : StructureT
    use AtomPotentialMod

    implicit none

    type( hamInfoT ),   pointer :: hamInfo
    type( StructureT ), pointer :: struct

    integer :: error, i, j


    allocate( hamInfo%potentialV(hamInfo%gVectors%numGVectors), stat=error )
	if( error /= 0 ) stop 'Error allocating potential space.'
                
    hamInfo%potentialV(:) = 0

    do j = 1, struct%atomBasis%numSpecies 
      do i = 1, hamInfo%gVectors%numGVectors

        ! Compute potential at G-vectors
        hamInfo%potentialV(i) = hamInfo%potentialV(i) + &
                & atomPotential(two * hamInfo%gVectors%kineticE(i), struct%atomBasis%elementLabel(j)) * &
                & hamInfo%gVectors%structureFactor(i,j)

      end do
    end do

    ! Explicitely put the potential at the origin to be zero.
    hamInfo%PotentialV(hamInfo%gVectors%gVecIndex(0,0,0)) = cmplx(0.d0,0.d0)

  end Subroutine PotentialInit

!*************************************************************************
  Subroutine PotentialDestroy( hamInfo )
!*************************************************************************

    implicit none

    type( hamInfoT ), pointer :: hamInfo
    integer :: error


    deallocate( hamInfo%potentialV, stat=error )
    if(error /= 0) stop 'Error deallocating potential space.'

  end Subroutine PotentialDestroy


!*************************************************************************
  Subroutine FindhMaxSize( hamInfo, struct, kPoints )
!*************************************************************************
!
!  This subroutine goes through the same loops as for generating the 
!  in order to find the maximum memory needed 
!
!  Outpout is hamInfo%hMaxSize which is defined in module variables and is public
!
    use SysParams,     only : double, two
    use StructureMod,  only : StructureT
    use kPointsMod,    only : kPointsT

    implicit none

    type( StructureT ),   pointer :: struct
    type( kPointsT ),     pointer :: kPoints
    type( hamInfoT ),     pointer :: hamInfo
    integer, dimension(3)          :: gg
    integer                        :: i,j,error,indx
    integer                        :: k, hSize, ktemp 
    real(double)                   :: ek
    real(double), dimension(3)     :: kPlusG

    hamInfo%hMaxSize = 0

    do k = 1, kPoints%numKPoints

    hSize = 0
    kPlusG = 0

! This section is
! inefficient since their are a lot of G vectors in gVectors that
! will never be within the energy cut off radius, because all of the
! G vectors within two times the cut off radius are stored in gVectors
! for the structure factor calculation.

       do i = 1, hamInfo%gVectors%numGVectors

         kPlusG(1:struct%ndim) = kPoints%kPointsLat(:,k) + hamInfo%gVectors%gVectors(:,i)
         ek = dot_product(kPlusG(1:struct%ndim), &
             & matmul(struct%lattice%bLatMetricTensor, kPlusG(1:struct%ndim)))/two

         if(ek <= hamInfo%gVectors%energyCutOff) then
           hSize        = hSize + 1
         end if

       end do


       if(hSize < kPoints%numBands) then
         print *, 'Error number of bands = ', kPoints%numBands, &
             & ' > size of the hamiltonian = ', hSize, &
             & ' for kPoint = ', k
         stop
       end if


       if(hSize > hamInfo%hMaxSize) then
	   hamInfo%hMaxSize = hSize				!hamInfo%hMaxSize set to largest hSize
           ktemp = k
       end if 

    end do
	   
!  print *, 'Maximum size of the hamiltonian = ', hamInfo%hMaxSize, &
!                & ' for k point = ', ktemp, numKPoints
!  print *, 'In recip. lattice basis vectors, this k point is:'
!  print *, kPointsLat(:,ktemp)
	
  end Subroutine FindhMaxSize

!*************************************************************************
  Subroutine KineticGenerate(hamInfo,struct,kPoints,eigenStates,k,Kinetic, & 
       & hSize)
!****************************************************************************

!    This subroutine will fill up the kinetic energy T_{G} array of the 
! hamiltonian. T_{G} = |k+G|^2/2, and the range of G is such that the 
! energy cutoff is satisfied.

    use EigenStatesMod,   only : eigenStatesT
    use SysParams,        only : double, two
    use StructureMod,     only : StructureT
    use kPointsMod,       only : kPointsT

    Implicit None

    type( StructureT ),   pointer :: struct
    type( kPointsT ),     pointer :: kPoints
    type( eigenStatesT ), pointer :: eigenStates
    type( hamInfoT ),     pointer :: hamInfo

    integer, intent(in)           :: k
    integer                       :: hSize,i
    real(double)                  :: ek
    real(double), dimension(3)    :: kPlusG
    real(double), dimension(:)    :: Kinetic


    hSize = 0
    kPlusG = 0

    do i = 1, hamInfo%gVectors%numGVectors
       kPlusG(1:struct%ndim) = kPoints%kPointsLat(:,k) + & 
            & hamInfo%gVectors%gVectors(:,i) 
       ek = dot_product(kPlusG(1:struct%ndim), &
            & matmul(struct%lattice%bLatMetricTensor, & 
            & kPlusG(1:struct%ndim)))/two

       if(ek <= hamInfo%gVectors%energyCutOff) then
          hSize        = hSize + 1
          Kinetic(hSize) = ek
          eigenStates%basisIndex(hSize,k) = i

       end if
    end do

  end Subroutine KineticGenerate

!*****************************************************************************
  Subroutine PotentialRGenerate(hamInfo,Grid)
!*****************************************************************************

!   This subroutine will put the Potential_{G} array on the grid, the backward
! fourier transform into V_{r}

    Type(hamInfoT),           pointer :: hamInfo
    Type(GridT),              pointer :: Grid
    integer                           :: G,i,j,k,L,s
    real(8)                           :: theta,phase

    ! size of the VGrid has to be upto maxGLen
    Grid%V = cmplx(0.d0,0.d0)

    ! Translation
    do i = 1, 2*Grid%Size(1)
       do j = 1, 2*Grid%Size(2)
          do k = 1, 2*Grid%Size(3)
             Grid%V(i,j,k) = hamInfo%PotentialV(Grid%GInvIndex( &
                  & i-1-Grid%Size(1),j-1-Grid%Size(2),k-1-Grid%Size(3)))
          end do
       end do
    end do

    ! This is so that the final grid VPsi_G is centered at the origin.
    do i = 1, 2*Grid%Size(1)
       do j = 1, 2*Grid%Size(2)
          do k = 1, 2*Grid%Size(3)
             s = i+j+k-3
             Grid%V(i,j,k) = Grid%V(i,j,k)*(-1)**s
          end do
       end do
    end do

  end Subroutine PotentialRGenerate

!*************************************************************************
  Subroutine HGenerate( hamInfo, struct, kPoints, eigenStates, &
                      & k, h, sqH, hSize )
!*************************************************************************


    use EigenStatesMod,   only : eigenStatesT
    use SysParams,        only : double, two
    use StructureMod,     only : StructureT
    use kPointsMod,       only : kPointsT

    implicit none

    type( StructureT ),   pointer :: struct
    type( kPointsT ),     pointer :: kPoints
    type( eigenStatesT ), pointer :: eigenStates
    type( hamInfoT ),     pointer :: hamInfo
    integer, intent(in)            :: k
    complex(double), intent(inout) :: h(:)     
    complex(double), intent(inout) :: SqH(:,:)     
    integer, intent(inout)         :: hSize

    integer, dimension(3)          :: gg       
    integer                        :: i,j,error,indx
    real(double)                   :: ek
    real(double), dimension(3)     :: kPlusG

    !!!!!!!!!!!!!!!! CHANGED ddas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: N
    N = 180
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    hSize = 0
    kPlusG = 0

    ! For each G vector calculate the kinetic energy if it is less than the 
    ! cut off energy add a row for it on the hamiltonian insert the kinetic
    ! energy along the diagonal.  BasisIndex is used to store the indices 
    ! of the G vectors that are in the hamiltonian.  This section is
    ! inefficient since their are a lot of G vectors in gVectors that
    ! will never be within the energy cut off radius, because all of the
    ! G vectors within two times the cut off radius are stored in gVectors
    ! for the structure factor calculation.

    do i = 1, hamInfo%gVectors%numGVectors

      kPlusG(1:struct%ndim) = kPoints%kPointsLat(:,k) + hamInfo%gVectors%gVectors(:,i) 
      ek = dot_product(kPlusG(1:struct%ndim), &
	     & matmul(struct%lattice%bLatMetricTensor, kPlusG(1:struct%ndim)))/two

      if(ek <= hamInfo%gVectors%energyCutOff) then

        hSize        = hSize + 1
        h(hSize * (hSize + 1) / 2) = ek

        SqH(hSize,hSize) = ek                       ! For the square matrix

        eigenStates%basisIndex(hSize,k) = i

      end if
    end do

!!$    !!!!!!!!!!!!!!! CHANGED ddas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    do i = 1, N
!!$       SqH(i,i) = cmplx(1.d0,0.d0)    
!!$    end do
!!$
!!$    do i = 1, N
!!$       do j = i+1, N
!!$          SqH(i,j) = cmplx(0.d0,0.d0)
!!$          SqH(j,i) = conjg(SqH(i,j))
!!$       end do
!!$    end do
!!$    !!!!!!!!!!!!!!! CHANGED ddas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    eigenStates%numBasisVectors(k) = hSize

    if(hSize < kPoints%numBands) then
      print *, 'Error number of bands = ', kPoints%numBands, &
             & ' > size of the hamiltonian = ', hSize
      stop
    end if

    ! For each row in the hamiltonian store the kinetic energy along the
    ! diagonal and calculate the off diagonal terms.
    do i = 1, hSize

      do j = i+1, hSize

        gg(1:struct%ndim) = &
      & hamInfo%gVectors%gVectors(1:struct%ndim, eigenStates%basisIndex(i, k))&
    & - hamInfo%gVectors%gVectors(1:struct%ndim, eigenStates%basisIndex(j, k))
        if(struct%ndim < 3) gg(struct%ndim+1:3) = 0
        indx = hamInfo%gVectors%gVecIndex(gg(1),gg(2),gg(3))    ! index for G-G' vector
        h(i+j*(j-1)/2) = hamInfo%potentialV(indx)

        SqH(i,j) = hamInfo%potentialV(indx)
        SqH(j,i) = conjg(SqH(i,j))

      end do
    end do

  end Subroutine HGenerate

!*************************************************************************
  Subroutine CalculateDensity( hamInfo, struct, kPoints, eigenStates, densityArray )
!*************************************************************************
    use GVectorsMod,     only : GVectorScatter
    use EigenStatesMod,  only : eigenStatesT
    use densityArrayMod, only : densityArrayT
    use StructureMod,    only : StructureT
    use kPointsMod,      only : kPointsT

    implicit none

    interface
      subroutine fft(grid, nn, ndim, isign)
        use sysparams, only : double
        implicit none
        integer, intent(in) :: isign, ndim, nn(ndim)
        complex(double), intent(inout) :: grid(:)
      end subroutine fft
    end interface

    type( StructureT ),     pointer :: struct
    type( kPointsT ),       pointer :: kPoints
    type( densityArrayT ), pointer :: densityArray
    type( eigenStatesT ),   pointer :: eigenStates
    type( hamInfoT ),       pointer :: hamInfo
    integer n, k

    densityArray%densityArray = 0

    do k = 1, kPoints%numKPoints
      do n = 1, kPoints%numBands
        call GVectorScatter( hamInfo%gVectors, struct, eigenStates, densityArray, k, n)
        call fft( densityArray%waveFunction, densityArray%gridDimensions, struct%ndim, -1)
        densityArray%waveFunction = densityArray%waveFunction / densityArray%gridSize
        densityArray%densityArray = densityArray%densityArray + densityArray%waveFunction*conjg(densityArray%waveFunction)
      end do
    end do

  end Subroutine CalculateDensity

!*************************************************************************
  Subroutine SpecificOrbitalDensity( hamInfo, struct, kPoints, eigenStates, densityArray, tagHandler )
!*************************************************************************
    use TagHandlerMod,   only : FindTag, tagHandlerT
    use EigenStatesMod,  only : eigenStatesT
    use StructureMod,    only : StructureT
    use kPointsMod,      only : kPointsT
    use GVectorsMod,     only : GVectorScatter
	use densityArrayMod, only : densityArrayT

    implicit none

    interface
      subroutine fft(grid, nn, ndim, isign)
        use sysparams, only : double
        implicit none
        integer, intent(in) :: isign, ndim, nn(ndim)
        complex(double), intent(inout) :: grid(:)
      end subroutine fft
    end interface

    type( StructureT ),     pointer :: struct
    type( kPointsT ),       pointer :: kPoints
    type( densityArrayT ), pointer :: densityArray
    type( eigenStatesT ),   pointer :: eigenStates
    type( tagHandlerT ),    pointer :: tagHandler
    type( hamInfoT ),       pointer :: hamInfo
    integer orbit, k, error

    call FindTag( tagHandler,'DensityForSpecificOrbital', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I3)', iostat = error) orbit
      if(error /= 0) stop 'Error, couldn''t read DensityForSpecificOrbital.'

      print *, 'DensityForSpecificOrbital ', orbit

    else

      stop 'Error, could not find DensityForSpecificOrbital.'

    endif

    if(orbit .gt. kPoints%numBands) stop 'Error, orbit exceeds number of bands.'

    do k = 1, kPoints%numKPoints
      call GVectorScatter( hamInfo%gVectors, struct, eigenStates, densityArray, k, orbit)
      call fft(densityArray%waveFunction, densityArray%gridDimensions, struct%ndim, -1)
      densityArray%waveFunction = densityArray%waveFunction / densityArray%gridSize
      densityArray%densityArray = densityArray%waveFunction*conjg(densityArray%waveFunction)
    end do

  end Subroutine SpecificOrbitalDensity

!******************************************************************************
  Subroutine GridInit(Grid,hamInfo,eigenStates,ndim)
!******************************************************************************

    type(GridT),                               pointer :: Grid
    type(hamInfoT),                            pointer :: hamInfo
    Type(eigenStatesT),                        pointer :: eigenStates    
    integer,                                   intent(in) :: ndim      

    integer                                            :: r

    !  Allocate space for the Grid Type
    allocate(Grid)

    !  Length of the Grid in each direction
    Grid%Size = hamInfo%gVectors%maxGLen

    Grid%fac = dble(8*Grid%Size(1)*Grid%Size(2)*Grid%Size(3))

    !  Number of GVectors
    Grid%NumGVec = hamInfo%gVectors%numGVectors

    !  Allocate space for the VGrid, and initialise it to zero
    allocate(Grid%V(2*Grid%Size(1),2*Grid%Size(2),2*Grid%Size(3)))
    allocate(Grid%Psi(2*Grid%Size(1),2*Grid%Size(2),2*Grid%Size(3)))
    allocate(Grid%VPsi(2*Grid%Size(1),2*Grid%Size(2),2*Grid%Size(3)))

    allocate(Grid%VPsiG(-Grid%Size(1):Grid%Size(1),-Grid%Size(2):Grid%Size(2),&
         &  -Grid%Size(3):Grid%Size(3)))


    Grid%V   = cmplx(0.d0,0.d0)

    !  Allocate and Initialise the Grid Index
    allocate(Grid%GVecIndex(ndim,Grid%NumGVec))
    Grid%GVecIndex = hamInfo%gVectors%gVectors

    allocate(Grid%GInvIndex( &
         & -hamInfo%gVectors%maxGLen(1):hamInfo%gVectors%maxGLen(1), &
         & -hamInfo%gVectors%maxGLen(2):hamInfo%gVectors%maxGLen(2), &
         & -hamInfo%gVectors%maxGLen(3):hamInfo%gVectors%maxGLen(3)))
    Grid%GInvIndex = hamInfo%gVectors%gVecIndex

    allocate(Grid%BasisIndex(EigenStates%maxBasisVectors, & 
         & EigenStates%numKPoints))

  end Subroutine GridInit

!*****************************************************************************

end Module pwHamMod
