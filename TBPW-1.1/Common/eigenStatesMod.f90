Module EigenStatesMod

! eigenVectors(i,j,k) the cooefficient of the i-th basis vector of the 
!					normalized eigenvector for the j-th eigenvalue at the 
!                   k-th k-point.
!                   [ For PW's The corresponding G-vector is
!                   gVectors(basisIndex(i),:) ]
!
! eigenValues(j,k)  j-th energy eigenvalue (in ascending order) at the k-th
!                   k point (i=1..numStates) (a.u.)
!
! basisIndex(i,k)   sorting vector defined for each k-point:
!                   BasisVectors(basisIndex(i,k),:) is
!                   a BasisVector for the k-point.

  use SysParams, only : double

  implicit none

  type eigenStatesT
    complex(double), pointer :: eigenVectors(:,:,:), ePsi(:,:)
    real(double),    pointer :: eigenValues(:,:)
    integer,         pointer :: basisIndex(:,:), numBasisVectors(:)
    integer                  :: numStates,hSize,hDIm
    integer                  :: numKPoints
    integer                  :: maxBasisVectors
    real( double )           :: tolerance
    integer                  :: maxIterations
  end type eigenStatesT

  contains

!*************************************************************************
  Subroutine EigenStatesInit( eigenStates, nStates, nKP, nBasis, tol, iter )
!*************************************************************************

    implicit none
    
    type( eigenStatesT ), pointer :: eigenStates
    integer, intent(in)           :: nStates, nKP, nBasis, iter
    real( double )                :: tol
	integer                       :: error
  
    if( associated( eigenStates )) stop 'Error, eigenStates already allocated.'
    allocate( eigenStates, stat = error )
    if( error /= 0 ) stop 'Error, eigenStates allocation failed.'

    eigenStates%numStates = nStates
    eigenStates%numKPoints = nKP
    eigenStates%maxBasisVectors = nBasis
    eigenStates%tolerance = tol
    eigenStates%maxIterations = iter

    allocate( eigenStates%eigenValues( eigenStates%numStates, &
            & eigenStates%numKPoints ), stat = error )
    if( error /= 0 ) stop 'Error allocating eigenValues in EigenStates.'

    allocate( eigenStates%eigenVectors( eigenStates%maxBasisVectors, &
            & eigenStates%numStates, eigenStates%numKPoints ), stat = error )
    if( error /= 0 ) stop 'Error allocating eigenVectors in EigenStates.'

    allocate( eigenStates%ePsi( eigenStates%maxBasisVectors, &
            & eigenStates%numKPoints ), stat = error )
    if( error /= 0 ) stop 'Error allocating eigenVectors in EigenStates.'

    allocate( eigenStates%basisIndex( eigenStates%maxBasisVectors, &
            & eigenStates%numKPoints ), stat = error )
    if( error /= 0 ) stop 'Error allocating basisIndex in EigenStates.'

    allocate( eigenStates%numBasisVectors( eigenStates%numKPoints ), &
            & stat = error)
    if( error /= 0 ) stop 'Error allocating numBasisVectors in EigenStates.'

  end Subroutine EigenStatesInit

!*************************************************************************
  Subroutine EigenStatesDestroy( eigenStates )
!*************************************************************************
 
    implicit none

    type( eigenStatesT ), pointer :: eigenStates
    integer :: error

    if( .not. associated( eigenStates )) stop &
        & 'Error, eigenStates not allocated.'

    deallocate( eigenStates%eigenVectors, eigenStates%eigenValues, &
              & eigenStates%ePsi, eigenStates%basisIndex, &
              & eigenStates%numBasisVectors, stat=error)
    if(error /= 0) stop 'Error eigenValues deallocating in KpointsDestroy.'

    deallocate( eigenStates, stat=error )
    if( error /= 0 ) stop 'Error, eigenState deallocation failed.'

    nullify( eigenStates )
 
  end Subroutine EigenStatesDestroy


!*************************************************************************
  Subroutine PrintEigenvalues(eigenStates, kpt)
!*************************************************************************

    implicit none
   
    type( eigenStatesT ), pointer :: eigenStates
	integer  i,j, kpt

	print *, 'eigenvalues  1 - ', eigenStates%numStates, &
          & ' at k points in lattice vector units '
	i = kpt
	  print *, 'kpoint ', i  ! , kpoints%kpointsLat(:,i)
	    do j = 1, eigenStates%numStates
	       print *, i, eigenStates%eigenValues(j, i)
		end do

  end Subroutine PrintEigenvalues

!*************************************************************************
end Module EigenStatesMod
!*************************************************************************
