!*************************************************************************
Module HamiltonianMod
!*************************************************************************
!	HDiagonalize calls LAPACK rountines

  use SysParams, 	only : double

  implicit none

  save

! hDim        Maximum dimension of the Hamiltonian matrix that gets used.
! hSize       Actual dimension of the Hamiltonian matrix that gets used.
!             The variable name is misleading, but we're stuck with it for now.
! h(i)        I-th element of the Hamiltonian matrix in packed storage
! SqH(i,j)    The complete Hamiltonian matrix 
! work arrays for LAPACK/Conjugate gradient

  type hamiltonianArrayT
    integer                   :: hSize, hDim
    complex(double), pointer  :: h(:)
    complex(double), pointer  :: SqH(:,:)
    complex(double), pointer  :: PotentialR(:)
    real(double),    pointer  :: KineticG(:)
                                            
    complex(double), pointer  :: work(:)
    real(double),    pointer  :: rwork(:)
    integer,         pointer  :: iwork(:)
    integer,         pointer  :: fail(:)
  end type hamiltonianArrayT

! The module contains the public functions:
!	HInit
!	HDestroy
!	HDiagonalize (calls LAPACK rountines)
!	HPrint
! The module contains the public variables:
! The module contains the private functions:
! The module contains the private variables:
 
  contains

!*****************************************************************************
  Subroutine HInit(hamiltonian, maxHDim)
!*****************************************************************************
!   Allocate space for the Hamiltonian and the work matrices.

    implicit none

    type( hamiltonianArrayT ), pointer :: hamiltonian
    integer, intent(in) :: maxHDim

    ! Workspace variables
    integer             :: error

    if( associated( hamiltonian )) stop 'Error, hamiltonian already allocated.'
    allocate( hamiltonian, stat = error )
    if( error /= 0 ) stop 'Error, hamiltonian allocation failed.'

    hamiltonian%hDim  = maxHDim
    
    ! Assume we used the entire Hamiltonian matrix.
    hamiltonian%hSize = maxHDim

    ! The Hamiltonian matrix is stored in packed form. Only the upper triangle
    ! plus half the main diagonal is stored. This is why the allocate statement
    ! looks funny.
    allocate( hamiltonian%h( hamiltonian%hDim * ( hamiltonian%hDim + 1 ) / 2 ), &
            & hamiltonian%sqH( hamiltonian%hDim, hamiltonian%hDim ), &
            & hamiltonian%KineticG(hamiltonian%hDim), &
            & hamiltonian%PotentialR(hamiltonian%hDim), &
            & hamiltonian%fail( hamiltonian%hDim ), stat=error )
    allocate( hamiltonian%work( 2 * hamiltonian%hDim ), &
            & hamiltonian%rwork( 7 * hamiltonian%hDim ), &
            & hamiltonian%iwork( 5 *hamiltonian%hDim ), stat=error )
    if(error /= 0) stop 'Error in allocating Hamiltonian space.'

  end Subroutine HInit

!*****************************************************************************
  Subroutine HDestroy( hamiltonian )
!*****************************************************************************

    implicit none

    type( hamiltonianArrayT ), pointer :: hamiltonian
    integer	:: error

    if( .not. associated( hamiltonian )) stop 'Error, hamiltonian not allocated.'

    deallocate( hamiltonian%h, hamiltonian%sqH, &
              & hamiltonian%work, hamiltonian%rwork, &
              & hamiltonian%iwork, hamiltonian%fail, stat=error )
    if(error /= 0) stop 'Error in deallocating Hamiltonian space.'

    deallocate( hamiltonian, stat = error )
    if(error /= 0) stop 'Error, hamiltonian deallocation failed.'

    nullify( hamiltonian )

  end Subroutine Hdestroy

!*****************************************************************************
  Subroutine HDiagonalize( hamiltonian, eigenStates, k )
!*****************************************************************************
! Hamiltonian matrix is diagonalized by a call to LAPACK routines.

    use SysParams,      only : one, zero
    use EigenStatesMod, only : eigenStatesT

    implicit none

    type( hamiltonianArrayT ), pointer :: hamiltonian
    type( eigenStatesT ),      pointer :: eigenStates

    ! Workspace variables
    integer, intent( in )         :: k
    integer                       :: neig, error

    ! Using macro to make external call more portable between machines
#ifdef AUSLAPACK
    call zhpevx_( 'V','I','U', hamiltonian%hSize, hamiltonian%h, &
               & zero, zero, 1, eigenStates%numStates, -one, &
               & neig, eigenStates%eigenValues(:,k), &
               & eigenStates%eigenVectors(:,:,k), &
               & eigenStates%maxBasisVectors, &
               & hamiltonian%work, hamiltonian%rwork, &
               & hamiltonian%iwork, hamiltonian%fail(1:hamiltonian%hSize), &
               & error )
#else
    call zhpevx( 'V','I','U', hamiltonian%hSize, hamiltonian%h, &
               & zero, zero, 1, eigenStates%numStates, -one, &
               & neig, eigenStates%eigenValues(:,k), &
               & eigenStates%eigenVectors(:,:,k), &
               & eigenStates%maxBasisVectors, &
               & hamiltonian%work, hamiltonian%rwork, &
               & hamiltonian%iwork, hamiltonian%fail(1:hamiltonian%hSize), &
               & error )
#endif

    if(error /= 0) stop 'Error: failed in zhpevx.'

  end Subroutine HDiagonalize




!*****************************************************************************
  Subroutine HPrint( hamiltonian )
!*****************************************************************************

    implicit none

    type( hamiltonianArrayT ), pointer :: hamiltonian
    integer :: i, j

    open(7, file = 'h.dat')

    do i = 1, hamiltonian%hSize
      do j = i, hamiltonian%hSize
        write( unit = 7, fmt = '(F8.4)', advance = 'NO') &
             & real( hamiltonian%h( i + j * ( j - 1 ) / 2 ))
      end do
      write(unit = 7, fmt = '()')
    end do


  end Subroutine HPrint

end Module HamiltonianMod


