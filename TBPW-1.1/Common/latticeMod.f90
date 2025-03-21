!*************************************************************************
Module LatticeMod
!*************************************************************************
! This module contains information on the lattice.
! The Lattice portion contains information on the direct and reciprocal 
!   Bravais lattices.
!
!*************************************************************************
! The module contains the public functions:
!    LatticeInit    	allocates arrays; reads lattice constant and vectors;
!                         calculates reciprocal lattice, metric tensors
!    LatticeDestroy 	deallocates arrays for the dynamic variables 
! 
! The module contains the public variables:
!    latticeConst          scale factor for lengths
!    aLatVec(i,j)          i-th component of the (scaled) j-th lattice vector
!    aLatVecLen(j)         length of the j-th lattice vector (a.u.)
!    aLatMetricTensor(i,j) = aLatMetricTensor(j,i): metric tensor in
!                          real space, used to scale the integer multiples
!                          of the reciprocal lattice vectors into real values.
!    bLatVec(i,j)          i-th component of the j-th reciprocal lattice vector
!    bLatVecLen(j)         length of the j-th reciprocal lattice vector
!    bLatMetricTensor(i,j) = bLatMetricTensor(j,i): metric tensor in
!                          reciprocal space, used to scale the integer multiples
!                          of the reciprocal lattice vectors into real values.
!
! The module contains the private functions:
!    InvertMatrix(a)    returns the inverse of a.
  
  use SysParams,    only : double
  implicit none

  type LatticeT
    real(double)             :: latticeConst
  
    real(double),    pointer :: aLatVec(:,:)  
    real(double),    pointer :: aLatVecLen(:)     
    real(double),    pointer :: aLatMetricTensor(:,:) 

    real(double),    pointer :: bLatVec(:,:)  
    real(double),    pointer :: bLatVecLen(:)     
    real(double),    pointer :: bLatMetricTensor(:,:) 
  end type LatticeT

  contains
!*************************************************************************
  Subroutine LatticeInit( lattice, ndim, tagHandler )
!*************************************************************************

    use TagHandlerMod,    only : FindTag, StringCompare, tagHandlerT
    use SysParams,     only : one, twopi

    implicit none

    type( LatticeT ),    pointer :: lattice
    type( tagHandlerT ), pointer :: tagHandler
	integer, intent(in)  :: ndim

    integer :: i, j, k, error
    character :: c
    character(80) :: buff
    logical :: default

!	if( associated( lattice )) stop 'Error: lattice already associated.'

	allocate( lattice, stat = error )

	if( error /= 0 ) stop 'Error allocation lattice.'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate spave for the lattice vectors, the reciprocal lattice
    ! vectors, and the metric tensor.
    allocate( lattice%aLatVec(ndim,ndim), &
            & lattice%bLatVec(ndim,ndim), &
            & lattice%aLatMetricTensor(ndim,ndim), &
            & lattice%bLatMetricTensor(ndim,ndim), &
            & lattice%aLatVecLen(ndim), lattice%bLatVecLen(ndim), &
            & stat = error)

    if(error /= 0) stop 'Error allocating vector space in LatticeInit.'
    ! Allocation succeded.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read in the lattice constant.  If there is no lattice constant in
    ! the input file or if there is an error reading it then the default
    ! value of 1.0d0 will be used.
	
	print *, 'Lattice Constant = scale factor for lengths'
    default = .True.
	
    call FindTag( tagHandler, 'LatticeConstant', error)
    if(error .eq. 0) then
      read(unit = tagHandler%fileno, fmt = '(G15.10)', iostat = error) lattice%latticeConst

      if(error .eq. 0) then
        print *, 'LatticeConstant ', lattice%latticeConst
        default = .False.
      end if
	  
    end if

    if(default) then
      lattice%latticeConst = one
      print *, 'Using default value of LatticeConstant 1.0.'
    end if
    ! latticeConst has been initialized
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read in the lattice vectors and scale them by the lattice constant.
    ! Calculate the lengths of the lattice vectors.
    ! If the lattice vectors are missing or there is an error reading them
    ! the program will halt.
	
    call FindTag( tagHandler, 'LatticeVectors', error)
	print *, 'Lattice Vectors as read in (unscaled) LatticeVectors'

    if(error .eq. 0) then		! read dummy variable c to start 
								! folowing read staement on next line
      read(unit = tagHandler%fileno, fmt ='(A1)') c  

      do i = 1, ndim
        read(unit = tagHandler%fileno, fmt = *, iostat = error) lattice%aLatVec(:,i)
        if(error /= 0) stop 'Error, couldn''t read LatticeVector.'
		print *, lattice%aLatVec(:,i)
      end do

    else
      stop 'Error, could not find lattice%aLatVec.'
    endif

    lattice%aLatVec = lattice%latticeConst * lattice%aLatVec ! rescale lattice vectors

    lattice%aLatMetricTensor = matmul(TRANSPOSE(lattice%aLatVec),lattice%aLatVec)
	
    do i = 1, ndim         ! determine lattice vector lengths
      lattice%aLatVecLen(i) = SQRT(dot_product(lattice%aLatVec(:,i),lattice%aLatVec(:,i)))
    end do
	
    ! lattice%aLatVec and lattice%aLatVecLen has been initialized and scaled.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculate the reciprocal lattice vectors using Invertmatrix.  Then
    ! calculate the metric tensor for the reciprocal lattice vectors.

    ! calculate reciprocal lattice
    lattice%bLatVec=twopi*TRANSPOSE(InvertMatrix(lattice%aLatVec))

    ! reciprocal space metric tensor
    lattice%bLatMetricTensor = matmul(TRANSPOSE(lattice%bLatVec),lattice%bLatVec)
	
    do i = 1, ndim          ! determine lattice vector lengths
      lattice%bLatVecLen(i) = SQRT(dot_product(lattice%bLatVec(:,i),lattice%bLatVec(:,i)))
    end do
	
    !  lattice%bLatVec and lattice%aLatMetricTensor have been initialized.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end Subroutine LatticeInit

!*************************************************************************
  Subroutine LatticeDestroy( lattice )
!*************************************************************************

! !  Deallocate space for the lattice information 

    implicit none

    type( LatticeT ), pointer :: lattice
    integer :: error

    deallocate(lattice%aLatVec, lattice%bLatVec, lattice%aLatMetricTensor, lattice%bLatmetricTensor, &
             & lattice%aLatVecLen, lattice%bLatVecLen, stat = error)
    if(error /= 0) stop 'Error deallocating lattice%aLatVec in LatticeDestroy.'

  end Subroutine LatticeDestroy

!*************************************************************************
  Subroutine LatticePrint( lattice )
!*************************************************************************

    use SysParams,    only : twopi

    implicit none

    type( LatticeT ), pointer :: lattice

    print *, ' '
    print *, 'The scaled lattice vectors are:'
    print *, lattice%aLatVec
    print *, ' '
    print *, 'The reciprocal lattice vectors (in units of 2pi/a) are:'
    print *, lattice%bLatVec/twopi

  end Subroutine LatticePrint

!*************************************************************************
  Function InvertMatrix(a)
!*************************************************************************
!
! This Function is similar to that given in numerical recipes for 
! Fortran 77.  It is use Gauss-Jordan Elimination with full pivoting.
! See Numerical Recipes in Fortran 77, Chapter 2 Section 1.

    use SysParams,   only : double, zero, one

    implicit none

    real(double), intent(in)        :: a(:,:)
    real(double)                    :: InvertMatrix(ubound(a,1),ubound(a,2)) 
	real(double)                    :: biggest, invPivot, pivot, scale
	real, dimension(ubound(a,1))    :: temp
	integer, dimension(ubound(a,1)) :: pivotIndex, rowIndex, colIndex
	integer                         :: i, j, k, l, m, n, row, col


    n = ubound(a,1)                          ! n is the number of rows in a
    InvertMatrix = a
    pivotIndex = 0
    do i = 1, n
      biggest = 0.0d0
      do j = 1, n
        if(pivotIndex(j) .ne. 1) then
          do k = 1, n
            if(pivotIndex(k) .eq. 0) then
              if(abs(InvertMatrix(j,k)) .ge. biggest) then
                biggest = abs(InvertMatrix(j,k))
                row = j
                col = k
              end if
            else if(pivotIndex(k) .gt. 1) then
              stop 'Singular Matrix in InvertMatrix.'
            end if
          end do
        end if
      end do
      pivotIndex(col) = pivotIndex(col) + 1
      if(row .ne. col) then
        temp = InvertMatrix(row,:)
        InvertMatrix(row,:) = InvertMatrix(col,:)
        InvertMatrix(col,:) = temp
      end if
      rowIndex(i) = row
      colIndex(i) = col
      pivot = InvertMatrix(col,col)
      InvertMatrix(col,col) = one
      if(pivot .eq. zero) stop 'Singular Matrix in InvertMatrix.'
      invPivot = one / pivot
      InvertMatrix(col,:) = InvertMatrix(col,:) * invPivot
      do l = 1, n
        if(l .ne. col) then
          scale = InvertMatrix(l,col)
          InvertMatrix(l,col) = zero
          InvertMatrix(l,:) = InvertMatrix(l,:) - InvertMatrix(col,:) * scale
        end if
      end do
    end do
    do l = n, 1, -1
      if(rowIndex(l) .ne. colIndex(l)) then
        temp = InvertMatrix(:,rowIndex(l))
        InvertMatrix(:,rowIndex(l)) = InvertMatrix(:,colIndex(l))
        InvertMatrix(:,colIndex(l)) = temp
      end if
    end do

  end function InvertMatrix

end module LatticeMod
