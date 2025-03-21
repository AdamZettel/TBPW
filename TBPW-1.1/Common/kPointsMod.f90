!*************************************************************************
Module KPointsMod
!*************************************************************************
! This module contains information about the k-points.
!
!*************************************************************************
! The module contains the public functions:
!    KPointsInit    		reads information on k points and 
!                   	  allocates the space for the dynamic variables
!    KPointsDestroy		deallocates space for the dynamic variables
! 
! The module contains the type:
!   kPointsT which contains the variables:
! numKPoints            total number of k-points used
! kPointsLat(i,j)     i-th component (reciprocal lattice coordinates)
!                     of the j-th k-point
! kPointsCart(i,j)    i-th component (cartesian coordinates * 2pi / a)
!                     of the j-th k-point
! numBands         number of bands to be calculated
! Information for graphs

  use SysParams,	only : double
  use GraphMod,     only : graphT

  implicit none

  save

  type kPointsT
    integer               :: numKPoints       
    real(double), pointer :: kPointsCart(:,:)
    real(double), pointer :: kPointsLat(:,:)
    integer               :: numBands
  end type kPointsT

  integer, parameter           :: scaledByRecipLatVec = 1
  integer, parameter           :: scaledByTwoPi = 2

  private scaledByRecipLatVec
  private scaledByTwoPi
  contains

!*************************************************************************
  Subroutine KPointsInit( kPoints, graphInfo, ndim, lattice, tagHandler )
!*************************************************************************
!    generate k-points around a BZ for plotting band structure plots such that
!    each line segment has a constant density of points to make the graphics
!    look good

    use SysParams,	   only : one, zero, twopi
    use TagHandlerMod, only : FindTag, StringCompare, tagHandlerT
    use LatticeMod,    only : latticeT
    use GraphMod,	   only : graphT


    implicit none

    type( kPointsT ),    pointer :: kPoints
    type( LatticeT ),    pointer :: lattice
    type( graphT ),      pointer :: graphInfo
    type( tagHandlerT ), pointer :: tagHandler
	integer, intent(in)  :: ndim
    real(double), allocatable :: ktemp(:,:)   ! temporary list of k-points

    integer       :: i,j,npt,error, kPointsScale
    character(80) :: buff
    character     :: c
    logical       :: default

    if( associated( kPoints )) stop 'Error, kPoints already allocated.'

    allocate( kPoints, stat = error )
    if( error /= 0 ) stop 'Error, kPoints allocation failed.'

	nullify( graphInfo )  ! used only if associated by another subroutine

    call FindTag( tagHandler, 'NumberOfBands', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) kPoints%numBands
      if(error /= 0) stop 'Error, couldn''t read NumberOfBands.'

      print *, 'NumberOfBands  ', kPoints%numBands

    else

      stop 'Error, could not find NumberOfBands.'

    endif

! Determine whether kpoints are read from input file are in units of ReciprocalLatticeVectors
! or in Cartesian units.  Default is scaledByRecipLatVec
 
    default = .True.

    call FindTag(tagHandler, 'KPointsScale', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(A80)', iostat = error) buff

      if(error .eq. 0) then

        if(StringCompare(buff, 'ReciprocalLatticeVectors')) then

          kPointsScale = scaledByRecipLatVec
          print *, 'KPointsScale ReciprocalLatticeVectors'
          default = .False.

        else if(StringCompare(buff, 'TwoPi/a')) then

          kPointsScale = scaledByTwoPi
          print *, 'KPointsScale TwoPi/a'
          default = .False.

        end if
      end if
    end if

    if(default) then
      kPointsScale = scaledByRecipLatVec
      print *, 'Using default KPointsScale ReciprocalLatticeVectors'
    end if

!   Start Various Possible types of generation of k ponts

    call FindTag( tagHandler, 'KPointsList', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) kPoints%numKPoints
      if(error /= 0) stop 'Error: Could not read number of KPoints after KPointsList.'

      print *, 'Use k point list read from input file '
	  print *, 'Calculation will be done for ' , kPoints%numKPoints, 'k points' 

      allocate(kPoints%kPointsCart(ndim, kPoints%numKPoints), &
        & kPoints%kPointsLat(ndim, kPoints%numKPoints), &
        & ktemp(ndim, kPoints%numKPoints), &
        & stat=error) ! allocate memory for k-points, and temporary ktemp 
      if(error /= 0) stop 'could not allocate memory for k-points'


      read(unit = tagHandler%fileno,fmt = '(A1)', iostat = error) c  !needed to read data from following line

		do j = 1, kPoints%numKPoints
			read(unit = tagHandler%fileno, fmt = *, iostat = error) ktemp(:,j)
			if(error /= 0) stop 'Error reading KPoints.'
			print *, j, ' ', ktemp(:,j)
		enddo

			print *, ' HERE WE ASSUME INPUT K ARE IN LATTICE VECTOR UNITS - NEED TO CAHNGE'

		do j = 1, kPoints%numKPoints
			kPoints%kPointsLat(:,j) = ktemp(:,j)
			print *, j, ' ', ktemp(:,j)
		enddo

    else
      !Checking for next option - read lines and divisions for plotting

      call FindTag( tagHandler, 'NumberOfLines', error)

      if(error .eq. 0) then

	    call KPointsLinesandGraphInit( kPoints, graphInfo, ndim, lattice, tagHandler )

      else
        if(error /= 0) stop 'Error, could not find any tag for kpoints'
	  endif

    endif

 
  end Subroutine KPointsInit

!*************************************************************************
  Subroutine KPointsDestroy( kPoints )
!*************************************************************************

    use GraphMod,	only : GraphDestroy

    implicit none

    type( kPointsT ), pointer :: kPoints
    integer :: error

    if( .not. associated( kPoints )) stop 'Error, kPoints not allocated.'

!    if( associated( kPoints%graph ))  call GraphDestroy( kPoints%graph )

!    deallocate(kPoints%kPointsLat, kPoints%kPointsCart,kPoints%xAxis, stat=error)
	deallocate(kPoints%kPointsLat, kPoints%kPointsCart, stat=error)
    if(error /= 0) stop 'Error deallocating kPointsLat in KPointsDestroy.'

    deallocate( kPoints, stat = error )
    if(error /= 0) stop 'Error deallocating kPoints in KPointsDestroy.'

    nullify( kPoints )

  end Subroutine KPointsDestroy


!*************************************************************************
  Subroutine MonkHorstPackInit( kPoints, ndim, lattice, tagHandler )
!*************************************************************************

    use SysParams,     only : one, half, twopi, double
    use TagHandlerMod,	   only : FindTag, tagHandlerT
    use LatticeMod,    only : LatticeT
!RMM    use DimensionsMod, only : DimensionsT

    implicit none

    type( kPointsT ),    pointer :: kPoints
!RMM    type( DimensionsT ), pointer :: dimensions
    type( LatticeT ),    pointer :: lattice
    type( tagHandlerT ), pointer :: tagHandler
	integer, intent(in)  :: ndim

    integer, allocatable      :: ndiv(:)  ! number of divisions per Dimension
    real(double), allocatable :: initialK(:), deltaK(:)
    character     :: c
    integer       :: i, j, len, start, error
    logical       :: default

    call FindTag( tagHandler,'NumberOfOccupiedBands', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) kPoints%numBands
      if(error /= 0) stop 'Error, couldn''t read NumberOfOccupiedBands.'

      print *, 'NumberOfOccupiedBands  ', kPoints%numBands

    else

      stop 'Error, could not find NumberOfOccupiedBands.'

    endif

    allocate(ndiv(ndim), initialK(ndim), &
           & deltaK(ndim), stat = error)
    if(error /= 0) stop 'Error allocating ndiv in MonkHorstPackInit.'

    call FindTag( tagHandler,'MonkHorstPack', error)
   
    default = .True.

    if(error .eq. 0) then

      print *, 'MonkHorstPack'
      read(unit = tagHandler%fileno, fmt = '(A1)', iostat = error) c
      read(unit = tagHandler%fileno, fmt = *, iostat = error) ndiv(:)
      if(error /= 0) stop 'Error reading MonkHorstPack.'
      print *, ndiv(:)
      default = .False.

    end if

    if(default) then
      print *, 'Using Default MonkHorstPack'
      ndiv = 1
      print *, ndiv(:)
    end if

    kPoints%numKPoints = 1
    do i = 1, ndim 
      kPoints%numKPoints = kPoints%numKPoints * ndiv(i)
    end do


    allocate(kPoints%kPointsCart(ndim, kPoints%numKPoints), &
        & kPoints%kPointsLat(ndim, kPoints%numKPoints), &
        & stat=error) ! allocate memory for k-points
    if(error /= 0) stop 'could not allocate memory for k-points'

    deltaK = one / ndiv
    initialK = half * (deltaK - one)
    len = 1
    start = 2

    kPoints%kPointsLat(:,1) = initialK

    do i = 1, ndim
      do j = 2, ndiv(i)

        kPoints%kPointsLat(:,start:start+len-1) = kPoints%kPointsLat(:,start-len:start-1)
        kPoints%kPointsLat(i,start:start+len-1) = kPoints%kPointsLat(i,start:start+len-1) + deltaK(i)
        start = start + len

      end do

      len = len * ndiv(i)

    end do

    kPoints%kPointsCart = matmul(lattice%bLatVec, kPoints%kPointsLat)  * lattice%latticeConst / twopi

    deallocate(ndiv, initialK, deltaK, stat = error)
    if(error /= 0) stop 'Error deallocating num in MonkHorstPackInit.'

  end Subroutine MonkHorstPackInit

!*************************************************************************
  Subroutine MonkHorstPackDestroy( kPoints )
!*************************************************************************

    implicit none

    type( kPointsT ), pointer :: kPoints
    integer :: error

    if( .not. associated( kPoints )) stop 'Error, kPoints not allocated.'

!    deallocate(kPoints%kPointsLat, kPoints%kPointsCart, kPoints%xAxis, stat=error)
	deallocate(kPoints%kPointsLat, kPoints%kPointsCart, stat=error)
    if(error /= 0) stop 'Error deallocating kPointsLat in KPointsDestroy.'

    deallocate( kPoints, stat = error )
    if(error /= 0) stop 'Error deallocating kPoints in KPointsDestroy.'

    nullify( kPoints )
  
  end Subroutine MonkHorstPackDestroy


!*************************************************************************
  Subroutine KPointsLinesandGraphInit( kPoints, graphInfo, ndim, lattice, tagHandler )
!*************************************************************************
!    generate k-points around a BZ for plotting band structure plots such that
!    each line segment has a constant density of points to make the graphics
!    look good

    use SysParams,	   only : one, zero, twopi
    use TagHandlerMod,    only : FindTag, StringCompare, tagHandlerT
    use LatticeMod,    only : LatticeT
    use GraphMod,	   only : AddLabel, GraphInit, graphT


    implicit none

    type( LatticeT ),    pointer :: lattice
    type( kPointsT ),    pointer :: kPoints
    type( tagHandlerT ), pointer :: tagHandler
	type( graphT ), pointer      :: graphInfo
	integer, intent(in)  :: ndim

    integer                   :: nlines       ! number of lines in k-space to plot
    integer                   :: ndiv         ! number of divisions per line
    integer,      allocatable :: num(:)       ! number of divisions per line
    real(double), allocatable :: length(:)    ! length of line
    real(double), allocatable :: kc(:,:)      ! list of circuit k-points
	real(double), allocatable :: tick(:)      ! list of positions of lables
	character(2), allocatable :: label(:)     ! list of labels for plot
    real(double), allocatable :: dk(:,:)      ! line vector	
    character(2)  :: c
    character(80) :: buff
    integer       :: i,j,npt,error, kPointsScale
    real(double)  :: delta_k
    logical       :: default



    call FindTag( tagHandler, 'NumberOfLines', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) nlines
      if(error /= 0) stop 'Error, couldn''t read NumberOfLines.'

      print *, 'NumberOfLines  ', nlines

    else

      stop 'Error, could not find NumberOfLines.'

    endif

    call FindTag( tagHandler, 'NumberOfDivisions', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) ndiv
      if(error /= 0) stop 'Error, couldn''t read NumberOfDivisions.'

      print *, 'NumberOfDivisions  ', ndiv

    else

      stop 'Error, could not find NumberOfDivisions.'

    endif

    allocate(kc(ndim, nlines + 1), dk(ndim, nlines), &
		   &label(nlines + 1), tick(nlines + 1), &
           & num(nlines), length(nlines), stat = error) 
    if(error /= 0) stop 'Error allocating k-points space in KPointsInit.'

! Determine whether kpoints are read from input file are in units of ReciprocalLatticeVectors
! or in Cartesian units.  Default is scaledByRecipLatVec
 
    default = .True.

    call FindTag(tagHandler, 'KPointsScale', error)

    if(error .eq. 0) then

      read(unit = tagHandler%fileno, fmt = '(A80)', iostat = error) buff

      if(error .eq. 0) then

        if(StringCompare(buff, 'ReciprocalLatticeVectors')) then

          kPointsScale = scaledByRecipLatVec
          print *, 'KPointsScale ReciprocalLatticeVectors'
          default = .False.

        else if(StringCompare(buff, 'TwoPi/a')) then

          kPointsScale = scaledByTwoPi
          print *, 'KPointsScale TwoPi/a'
          default = .False.

        end if
      end if
    end if

    if(default) then
      kPointsScale = scaledByRecipLatVec
      print *, 'Using default KPointsScale ReciprocalLatticeVectors'
    end if


    call FindTag( tagHandler,'KPointsAndLabels', error)
	print *, 'KPointsAndLabels'
    if(error .eq. 0) then

      read(unit = tagHandler%fileno,fmt = '(A1)', iostat = error) c

      read(unit = tagHandler%fileno,fmt =*, iostat = error) kc(:,1), label(1)

      print *, kc(:,1),' ', label(1)
!      call AddLabel( graphInfo, label, zero )
      tick(1) = 0.0d0

      do i = 1, nlines

        ! input circuit information
        read(unit = tagHandler%fileno, fmt = *, iostat = error) kc(:,i+1),label(i+1)
        dk(:,i)    = kc(:,i+1) - kc(:,i)
        length(i)  = SQRT(dot_product(matmul(dk(:,i),lattice%bLatMetricTensor),dk(:,i)))
        tick(i+1) = tick(i) + length(i)
!        call AddLabel( graphInfo, label(i+1), tick(i+1) )
        print *, kc(:,i+1),' ', label(i+1)

      end do

    else

      stop 'Error, could not find KPointsAndLabels.'

    end if

    num = (length/minval(length))*ndiv          ! number of divisions per line

    kPoints%numKPoints  = SUM(num) + 1                        ! total number of k-points

    allocate(kPoints%kPointsCart(ndim, kPoints%numKPoints), &
        & kPoints%kPointsLat(ndim, kPoints%numKPoints), stat=error) ! allocate memory for k-points
    if(error /= 0) stop 'Error, could not allocate memory for k-points.'


 ! generate a continuous set of k-points around the circuit 
    npt = 0
    do i = 1, nlines                             
      delta_k     = one/num(i)                
	  npt         = npt + 1
      kPoints%kPointsLat(:,npt) = kc(:,i)
      do j = 1, num(i) - 1
        npt         = npt + 1
        kPoints%kPointsLat(:,npt) = kc(:,i) + j * delta_k * dk(:,i)
      end do
    end do

    kPoints%kPointsLat(:, npt+1) = kc(:, nlines+1)                ! get the last k-point

 ! Generate information for plotting bands
  	 
   call GraphInit(graphInfo, kpoints%numKpoints, nlines)

      do i = 1, nlines + 1
        call AddLabel( graphInfo, label(i), tick(i) )
      end do

    graphInfo%xAxis = 0

    do i = 2, kPoints%numKPoints

      dk(:,1)        = kPoints%kPointsLat(:,i) - kPoints%kPointsLat(:,i-1)
      graphInfo%xAxis(i) = graphInfo%xAxis(i-1) + &
          & sqrt(dot_product(dk(:,1),matmul(lattice%bLatMetricTensor,dk(:,1))))

    end do

! Calculate k points in Cartesian coordinates for possible use

    if(kPointsScale .eq. scaledByRecipLatVec) then
      kPoints%kPointsCart = matmul(lattice%bLatVec, kPoints%kPointsLat)  &
                          & * lattice%latticeConst / twopi
    else if(kPointsScale .eq. scaledByTwoPi) then
      kPoints%kPointsCart = kPoints%kPointsLat
      kPoints%kPointsLat = matmul(transpose(lattice%aLatVec), kPoints%kPointsCart) / & 
                 & lattice%latticeConst * twopi
    end if

    deallocate(kc, dk, num, length, stat = error)
    if(error /= 0) stop 'Error deallocating kc, etc. in KPointsInit.'

  end Subroutine KPointsLinesandGraphInit

!*************************************************************************
  Subroutine PrintKPoints( kPoints )
!*************************************************************************

    implicit none

    type( kPointsT ), pointer :: kPoints
    integer i

    do i = 1, kPoints%numKPoints
      print *, kPoints%kPointsLat(:,i)
    end do

  end Subroutine PrintKPoints
end Module KPointsMod

