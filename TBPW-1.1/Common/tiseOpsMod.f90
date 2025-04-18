!*****************************************************************************
Module tiseOpsMod
!****************************************************************************

  use SysParams, only : double

  Implicit None

  complex( double ), external :: zdotu
  
Contains

!****************************************************************************
  Subroutine EigenValues( eigenStates, minVec, iBand, kPoint)
!*****************************************************************************

    use EigenStatesMod, only : eigenStatesT
    use TypeMod, only : minVecT

    integer, intent(in) :: iBand, kPoint
    type(minVecT),    pointer :: minVec   
    type(eigenStatesT), pointer :: eigenStates



    Call CalcIterEPsi( eigenStates, minVec, iBand, kPoint )

    ! Store the RitzValue in eigenStates
    eigenStates%eigenValues( iBand, kPoint ) = Dot_Product( &
         & eigenStates%eigenVectors(:, iBand, kPoint), &
         & eigenStates%ePsi( :, kPoint ))

  end Subroutine EigenValues

!****************************************************************************
  Subroutine CalcIterEPsi( eigenStates, minVec, iBand, kPoint )
!****************************************************************************

    use EigenStatesMod, only : eigenStatesT
    use TypeMod, only : minVecT

    integer, intent(in)    :: iBand, kPoint
    type(minVecT), pointer :: minVec
    type(eigenStatesT), pointer :: eigenStates

    eigenStates%ePsi( :, kPoint ) = &
        & minVec%cosine * eigenStates%ePsi( :, kPoint ) + &
        & minVec%sine * minVec%EPhi

  end Subroutine CalcIterEPsi

!****************************************************************************
  Subroutine CalculateEPsi( eigenStates, hMatrix,iBand, kPoint)
!****************************************************************************

    use EigenStatesMod, only : eigenStatesT
    use SysParams, only : double

    integer, intent( in ) :: iBand, kPoint
    complex( double ), dimension(:,:) :: hMatrix
    type( eigenStatesT ), pointer :: eigenStates


    integer :: nSize
    complex( double ) :: a, b

    a = cmplx(1.d0,0.d0)
    b = cmplx(0.d0,0.d0)

    nSize = eigenStates%hSize

    !    Evaluate hMatrix * EigenState, this matrix multiplication is done
    ! by a BLAS routine.
    ! Using macro to make external call more portable between machines
#ifdef AUSBLAS
    Call zgemv_('N', nSize, nSize, a, hMatrix, nSize, &
         & eigenStates%eigenVectors( :, iBand, kPoint ), &
         & 1, b, eigenStates%ePsi(:, kPoint ), 1 )
#else
    Call zgemv('N', nSize, nSize, a, hMatrix, nSize, &
         & eigenStates%eigenVectors( :, iBand, kPoint ), &
         & 1, b, eigenStates%ePsi(:, kPoint ), 1 )
#endif

  end Subroutine CalculateEPsi

!****************************************************************************
  Subroutine PreCondition( eigenStates, minVec, hMatrix, iBand, kPoint )
!****************************************************************************

    use EigenStatesMod, only : eigenStatesT
    use SysParams, only : double
    use TypeMod, only : minVecT

    integer, intent(in) :: iBand, kPoint
    complex( double ), dimension(:,:) :: hMatrix
    type(eigenStatesT), pointer :: eigenStates
    type(minVecT), pointer :: minVec
    
    integer :: i, nSize
    real( double ) :: eKinIter, x, k

    nSize = eigenStates%hSize

    !   There is a preconditioning matrix K for every iBand, to calculate
    ! that we have to first find the quantity x, the ratio of the kinetic
    ! energy with the kinetic energy of that iteration. this is different
    ! for every band
    eKinIter = cmplx( 0.d0, 0.d0 )
    do i = 1, nSize
       eKinIter = eKinIter + hMatrix( i, i ) * &
            & conjg( eigenStates%eigenVectors( i, iBand, kPoint )) * &
            & eigenStates%eigenVectors( i, iBand, kPoint )
    end do

    !   Now calculate x and K
    do i = 1, nSize
       x = hMatrix( i, i ) /( eKinIter + 1.d-30 )
       
       k = (27.d0 + 18.d0 * x + 12.d0 * x**2 + 8.d0 * x**3 ) / &
         & (27.d0 + 18.d0 * x + 12.d0 * x**2 + 8.d0 * x**3 + 16.d0 * x**4 )

       ! Update the StDescVector
       minVec%StDesc( i ) = k * minVec%StDesc( i )

    end do

  end Subroutine PreCondition

!***************************************************************************
  Subroutine ExactEnergyMinimize( eigenStates, minVec, hMatrix, iBand, kPoint )
!**************************************************************************

    use EigenStatesMod, only : eigenStatesT
    use SysParams, only : double
    use TypeMod, only : minVecT

    integer, intent(in) :: iBand, kPoint
    complex( double ), dimension(:,:) :: hMatrix
    type(eigenStatesT), pointer :: eigenStates
    type(minVecT), pointer :: minVec

    integer :: i, j, nSize
    real( double ) :: e0, theta, ePhi
    complex( double ) :: phiHPsi, a, b
    complex( double ), dimension(:), allocatable :: hPhi

    nSize = eigenStates%hSize

    a = cmplx(1.d0,0.d0)
    b = cmplx(0.d0,0.d0)

    e0 = cmplx( eigenStates%eigenValues( iBand, kPoint ))

    phiHPsi = Dot_Product( minVec%Phi, eigenStates%ePsi( :, kPoint ))

    allocate( hPhi( nSize ))

    !    Evaluate Phi.hMatrix.Phi, by calling the BLAS routine. First
    ! calulate hMatrix.Phi

    ! Using macro to make external call more portable between machines
#ifdef AUSBLAS
    Call zgemv_('N', nSize, nSize, a, hMatrix, nSize, minVec%Phi, 1, b, hPhi, 1 )
#else
    Call zgemv('N', nSize, nSize, a, hMatrix, nSize, minVec%Phi, 1, b, hPhi, 1 )
#endif
    ePhi = Dot_Product( minVec%Phi, hPhi )

    minVec%EPhi = hPhi
    minVec%EPhi(eigenStates%hsize+1:eigenStates%hDim) = cmplx(0.d0,0.d0)
    deallocate( hPhi )

    theta = 0.5d0 * abs( atan(( 2.d0 * real( phiHPsi )/( e0 - ePhi +10.d-30 ))))
!    print *, 'theta:',theta

    minVec%sine = sin( theta )
    minVec%cosine = cos( theta )

    eigenStates%eigenVectors( :, iBand, kPoint ) = &
         & eigenStates%eigenVectors( :, iBand, kPoint ) * minVec%cosine + &
         & minVec%Phi * minVec%sine

    ! Update the Ritz values
    Call EigenValues( eigenStates, minVec, iBand, kPoint )


  end Subroutine ExactEnergyMinimize

!***************************************************************************
  Subroutine CalcHPsi(size,T_G,vector1,vector2,Grid,iBand,kPoint)
!*****************************************************************************

    ! H in T_G and Grid
    ! Psi in vector1
    ! HPsi in vector2

    use sysparams,              only : double
    use typeMod,                only : GridT
    
    type(GridT),               pointer :: Grid
    integer, intent(in)                :: size
    real(double), dimension(:)         :: T_G
    complex(double), dimension(:)      :: vector1, vector2

    integer, intent(in)                :: iBand,kPoint

    integer :: NSize,r,s
    integer :: i,j,k,G
    real(8) :: fac

    ! Zero the Psi Grid
    Grid%Psi = cmplx(0.d0,0.d0)

    ! Put the  Statevector on the grid (Translation)
    do r = 1, size
       G = Grid%BasisIndex(r,Grid%k)
       Grid%Psi(Grid%GVecIndex(1,G)+Grid%Size(1)+1,Grid%GVecIndex(2,G)+ &
            & Grid%Size(2)+1, Grid%GVecIndex(3,G)+Grid%Size(3)+1) = &
            & vector1(r)
    end do
    
    ! Multiply V_r and Psi_r
    do i = 1, 2*Grid%Size(1)
       do j = 1, 2*Grid%Size(2)
          do k = 1, 2*Grid%Size(3)
!             s = i+j+k-3 ! If the output from potential is not centered.
             Grid%VPsi(i,j,k) = real(Grid%V(i,j,k))*Grid%Psi(i,j,k)!*(-1)**s
          end do
       end do
    end do

    ! Since this is a backward transform, scale the grid
    Grid%VPsi = Grid%VPsi/Grid%fac

    ! Evaluate (T+V)*StateVector
    do r = 1, size
       G = Grid%BasisIndex(r,Grid%k) 
       vector2(r) = T_G(r)*vector1(r) +  &
            & Grid%VPsi(Grid%GVecIndex(1,G)+ &
            & Grid%Size(1)+1,Grid%GVecIndex(2,G)+Grid%Size(2)+1, &
            & Grid%GVecIndex(3,G)+Grid%Size(3)+1)
    end do


  end Subroutine CalcHPsi

!*****************************************************************************
  Subroutine FFTPreCondition(eigenStates,minVec,T_g,iBand,kPoint)
!*****************************************************************************

    use eigenStatesMod,          only : eigenStatesT
    use sysparams,               only : double 
    use typeMod,                 only : minVecT

    type(eigenStatesT),         pointer :: eigenStates
    type(minVecT),              pointer :: minVec

    integer, intent(in)               :: iBand,kPoint
    real(double), dimension(:)        :: T_G

    integer                           :: r
    real(double)                      :: Ekin_Iter,x,K

    !   There is a preconditioning matrix K for every iband, to calculate
    ! that we have to first find the quantity x, the ratio of the kinetic
    ! energy with the kinetic energy of that iteration. this is different
    ! for every band
    Ekin_iter=cmplx(0.d0,0.d0)
    do r = 1, eigenStates%hSize
       Ekin_iter = Ekin_iter + T_g(r) * &
            & conjg(eigenStates%eigenVectors(r,iband,kPoint))* &
            & eigenStates%eigenVectors(r,iband,kPoint)
    end do

    !   Now calculate x and K
    do r = 1, eigenStates%hSize
       x = T_g(r)/(Ekin_iter+1.d-30)
       
       K =    (27.d0+18.d0*x+12.d0*x**2+8.d0*x**3)/ &
            & (27.d0+18.d0*x+12.d0*x**2+8.d0*x**3+16.d0*x**4)

       ! Update the StDescVector
       MinVec%StDesc(r) = K*MinVec%StDesc(r)

    end do


  end Subroutine FFTPreCondition

!******************************************************************
  Subroutine FFTEnergyMinimize(eigenStates,minvec,T_G,Grid,iBand,kPoint)
!******************************************************************

    use EigenStatesMod, only : eigenStatesT
    use SysParams, only : double
    use TypeMod, only : minVecT,GridT

    integer, intent(in) :: iBand, kPoint
    real( double ), dimension(:)  :: T_G
    type(eigenStatesT), pointer :: eigenStates
    type(minVecT), pointer :: minVec
    type(GridT), pointer   :: Grid

    integer            :: nSize
    real(double)       :: theta, ephi, e0
    complex(double)    :: a,b,phiHPsi
    complex( double ), dimension(:), allocatable :: hPhi

    nSize = eigenStates%hSize

    a = cmplx(1.d0,0.d0)
    b = cmplx(0.d0,0.d0)

    e0 = cmplx( eigenStates%eigenValues( iBand, kPoint ))
    phiHPsi = Dot_Product(minVec%Phi(1:nSize), &
         & eigenStates%ePsi(1:nSize,kPoint ))

    allocate(hPhi(nSize))

    !    Evaluate Phi.HMatrix.Phi, by calling the FFT routine. First
    ! calulate HMatrix.Phi
    Call CalcHPSi(nSize,T_G,minVec%Phi,hPhi,Grid,iBand,kPoint)

    ePhi = Dot_Product(MinVec%Phi(1:nSize),hPhi)
    MinVec%EPhi = hPhi
    minVec%EPhi(eigenStates%hsize+1:eigenStates%hDim) = cmplx(0.d0,0.d0)
    deallocate(HPhi)

    theta = 0.5d0*abs(atan((2.d0*real(PhiHPsi)/(e0 - ePhi +10.d-30))))

    MinVec%cosine = cos(theta)
    MinVec%sine = sin(theta)

    eigenStates%eigenVectors( :, iBand, kPoint ) = &
         & eigenStates%eigenVectors( :, iBand, kPoint ) * minVec%cosine + &
         & minVec%Phi * minVec%sine

    ! Update the Ritz values
    Call EigenValues( eigenStates, minVec, iBand, kPoint )

  end Subroutine FFTEnergyMinimize

!*******************************************************************

end Module tiseOpsMod
!******************************************************************







