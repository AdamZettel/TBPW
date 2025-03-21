Module densityArrayMod

! Module for density defined on a grid
! densityArrayT is a type containing the wavefunction and density information

!**************************************************************************
!      NEED TO CHENGE TO MAKE MORE UNIVERSAL
!************************************************************************ 
    use SysParams,    only : double

  type densityArrayT
    complex(double), pointer :: densityArray(:), waveFunction(:)
    integer,         pointer :: gridDimensions(:)
    integer                  :: gridSize
  end type densityArrayT

  contains

!*************************************************************************
  Subroutine densityArrayInit( densityArray, ndim, ndiv )
!*************************************************************************
    use SysParams,    only : two, log2
 !RMM   use DimensionsMod, only : DimensionsT
    
    implicit none
    
 !RMM    type( DimensionsT ), pointer :: dimensions
    type( densityArrayT ), pointer :: densityArray
	integer, intent(in) :: ndim
    integer, intent(in) :: ndiv
    integer             :: i, error

    if( associated( densityArray )) stop 'Error, densityArray already allocated.'
    allocate( densityArray, stat=error )
    if( error /= 0 ) stop 'Error, densityArray allocation failed.'

    densityArray%gridSize = 1

    allocate(densityArray%gridDimensions(ndim), stat = error)
    if(error /= 0) stop 'Error allocating gridDimensions.'

    do i = 1, ndim
      densityArray%gridDimensions(i) = int(two**ceiling(log10(dble(ndiv+1))/log2))
      densityArray%gridSize = densityArray%gridSize * densityArray%gridDimensions(i)
    end do

    print *, densityArray%gridDimensions(:)
    allocate( densityArray%densityArray(densityArray%gridSize), &
            & densityArray%waveFunction(densityArray%gridSize), stat = error )
    if(error /= 0) stop 'Error allocating denstiyMatrix in densityArrayInit.'

  end Subroutine densityArrayInit


!*************************************************************************
  Subroutine densityArrayPlot( densityArray, ndim, lattice )
!*************************************************************************

    use SysParams,     only : zero, two, half
    use LatticeMod,    only : LatticeT
 !RMM   use DimensionsMod, only : DimensionsT

    type( densityArrayT ), pointer :: densityArray
 !RMM    type( DimensionsT ),    pointer :: dimensions
    type( LatticeT ),       pointer :: lattice
	integer, intent(in) :: ndim

    integer :: index, i, j, k ,l , error, cartIndex(ndim)
    real :: maxDim(ndim), gridDist(ndim), distance, &
          & maxDist, cartCoordGrid(ndim), invMaxValue
    integer, allocatable :: latCoords(:,:)
    real, allocatable :: cartGridVal(:), minDist(:)

    open(7, file = 'density.txt')
    write(7, '(2I3)') 3, 1
    write(7, '(3I4)') densityArray%gridDimensions(:)
    maxDim = zero

    do i = 1, ndim
      maxDim = maxDim + lattice%aLatVec(:,i)
    end do

    allocate(cartGridVal(densityArray%gridSize), minDist(densityArray%gridSize), &
           & latCoords(ndim, densityArray%gridSize), stat = error)
    if(error /= 0) stop 'Error allocating cartGridVal in densityArrayPlot'

    gridDist = MaxDim/densityArray%gridDimensions
    maxDist = minval(gridDist)
    minDist = maxDist
    cartCoordGrid = zero
    cartGridVal = zero

    do j = 1, densityArray%gridSize

      index = j - 1

      do i = 1, ndim - 1

        latCoords( i, j ) = mod( index, densityArray%gridDimensions( i+1 ))
        index = int( index / densityArray%gridDimensions( i+1 ))

      end do

      latCoords( ndim, j ) = index

    end do

    do j = 1, densityArray%gridSize

      cartCoordGrid = zero

      do i = 1, ndim

        cartCoordGrid(:) = cartCoordGrid(:) + latCoords(i,j) * &
                           & lattice%aLatVec(:,i) / densityArray%gridDimensions(i)

      end do

      do i = 1, ndim

        cartIndex(i) = int(cartCoordGrid(i) / gridDist(i) + half) 

      end do

      index = cartIndex(ndim)

      do i = ndim-1, 1, -1

          index = index * densityArray%gridDimensions(i+1) + cartIndex(i)

      end do  

      index = index + 1

      distance = sqrt(dot_product(cartIndex(:) * gridDist - &
                 & cartCoordGrid(:), cartIndex(:) * gridDist - &
                 & cartCoordGrid(:)))
      if(distance .le. minDist(index)) then
        cartGridVal(index) = real( densityArray%densityArray( j ))
        minDist(index) = distance
      end if
    end do

    invMaxValue = 1 / maxval(cartGridVal)

    write(7, '(6E12.5)') -maxDim(3)/two,  maxDim(3)/two, -maxDim(2)/two, &
                        & maxDim(2)/two, -maxDim(1)/two,  maxDim(1)/two

    do i = 1, densityArray%gridSize
      write(7, '(E11.5)') cartGridVal(i) * invMaxValue
    end do

    deallocate(cartGridVal, minDist, &
           & latCoords, stat = error)
    if(error /= 0) stop 'Error deallocating cartGridVal in densityArrayPlot'

    close(7)

  end Subroutine densityArrayPlot

!*************************************************************************
  Subroutine densityArrayDestroy( densityArray )
!*************************************************************************
    implicit none
    
    type( densityArrayT ), pointer :: densityArray
    integer :: error

    if( .not. associated( densityArray )) stop 'Error, densityArray not allocated.'

    deallocate( densityArray%densityArray, densityArray%waveFunction, &
              & densityArray%gridDimensions, stat = error )
    if(error /= 0) stop 'Error deallocating densityArray.'

    deallocate( densityArray, stat = error )
    if(error /= 0) stop 'Error, densityArray deallocation failed.'

    nullify( densityArray )

  end Subroutine densityArrayDestroy

end Module densityArrayMod

