!*************************************************************************
Module SlaterKosterParamMod
!*************************************************************************
!*** THIS IS THE ONLY PART OF THE PROGRAM THAT CHANGES fOR SPECIFIC MODELS
!*************************************************************************

  use SysParams,			only : half, zero, one, double
  use SysParams,			only : epsilon, hbar2, electronMass
  
  implicit none

  complex(double), external                  :: zdotc  
  
  contains

!!RMM *******************************************************

!*************************************************************************
  Function TBHMatrixElem( vSK, xt, li, mi, lj, mj )
!*************************************************************************
! Calculates the geometric factors for the two center tight binding 
! hamiltonian matrix elements, where xt = vector between atoms.

    implicit none
    complex(double)             :: TBHMatrixElem
    
    real(double), intent(in)    :: vSK(:)    
    integer, intent(in)         :: li, mi, lj, mj
    real(double), intent(in)    :: xt(3)

!  From Koster tables - matrix elements easily written in terms
!  of direction cosines given by the x unit vector
!  See Harrison, "Elec. Str. and Props. of Solids", p. 481 
!  We require that for l = 1, m =1 for px; m=2 for py; m=3 for pz
!  where x,y,z are the 1,2,3 cartesian directions defined in the 
!  structure
     
    select case(li)

      case(0) 
    
        if (lj .eq. 0) TBHMatrixElem = vSK(1)           ! SS 
        if (lj .eq. 1) TBHMatrixElem = vSK(2) * xt(mj)  ! SP sigma
    
      case(1)
        
        if (lj .eq. 0) TBHMatrixElem = - vSK(2) * xt(mi) ! SP sigma
        ! Note minus x(mi) since x is from i to j

        if(lj .eq. 1) then              !PP sigma, pi

          TBHMatrixElem = (vSK(3) - vSK(4)) * xt(mi) * xt(mj)

          if (mi .eq. mj) then

            TBHMatrixElem = TBHMatrixElem + vSK(4)

          end if

        end if
        
    end  select
        
  end Function TBHMatrixElem

!*************************************************************************
  Subroutine TBVInit( vSK, tbModelType, tbParams, speciesI, speciesJ, distance )

!*************************************************************************

!   This subroutine supplies all the information specific to a given
!   material or model.  A case selection is made - new possible forms
!   must be programmed individually
    

    implicit none
 
    integer, intent(in)		   :: tbModelType
    integer, intent(in)		   :: speciesI   
    integer, intent(in)		   :: speciesJ
    real(double), intent(in) 	   :: distance
    real(double), intent(in)       :: tbParams(20)
    real(double), intent(out)      :: vSK(20)
    
    vSK(:) = 0   !resets all Slater Koster parameters to zero
    
    ! tbModelType read in input file with tag 'TightBindingModelType'

        if(tbModelType .eq. 1) then
       	  ! tbModelType = 1 is for Harrison Tight Binding empirical parameters.
          call TBHarrison( vSK, distance )

        else if(tbModelType .ge. 10) then
          ! tbModelType >= 10 is for special cases coded in the subroutine TBVInitSpecial        
          call TBParamsSpecial( vSK, tbModelType, tbParams, speciesI, speciesJ, distance )
        	
        else
        
	  stop 'Model type not supported'
        
        end if
             
  end Subroutine TBVInit

!*************************************************************************
  Subroutine TBHarrison( vSK, distance )
!*************************************************************************
!
! This subroutine supplies all the information specific to a given material
! or model --  Here we have examples of the Harrison matrix elements
! and examples of simple models
! For s and p states the parameters depend only on distance for all elements
! Ref: W. Harrison, "Electronic Structure of Solids".  
! Values are given in the "Solid State Table for the Elements".

! distance assumed to be in Bohr
    
    implicit none
    
    real(double), intent(in) 	   :: distance
    real(double), intent(out)      :: vSK(20)
    ! Harrison Parameters.
    real(double), parameter :: n_ss_0 = -1.40d0 * hBar2/electronMass
    real(double), parameter :: n_sp_0 =  1.84d0 * hBar2/electronMass
    real(double), parameter :: n_pp_0 =  3.24d0 * hBar2/electronMass
    real(double), parameter :: n_pp_1 = -0.81d0 * hBar2/electronMass


    real(double) :: invDist2

    invDist2 = 1.0d0 / distance
    invDist2 = invDist2 * invDist2

    !   V(1) = s orbital coupling
    !   V(2) = s p orbital coupling
    !   V(3) = p p parallel orbital coupling
    !   V(4) = p p perpendicular orbital coupling
    
    vSK(1) = n_ss_0 * invDist2
    vSK(2) = n_sp_0 * invDist2
    vSK(3) = n_pp_0 * invDist2
    vSK(4) = n_pp_1 * invDist2

  end Subroutine TBHarrison

!*************************************************************************
  Subroutine TBParamsSpecial( vSK, tbModelType, tbParams, speciesI, speciesJ, distance  )
!*************************************************************************
!
!   ***  THIS IS THE ONLY SUBROUTINE THAT CHANGES FOR SPECIFIC MOLDELS ***
!
!   This subroutine supplies all the information specific to a given
!   material or model.  A case selection is made - new possible forms
!   must be programmed individually
  
    implicit none
    
    integer, intent(in)		   :: tbModelType    
    integer, intent(in)		   :: speciesI   
    integer, intent(in)		   :: speciesJ   
    real(double), intent(in) 	   :: distance
    real(double), intent(in)       :: tbParams(20)
    real(double), intent(out)      :: vSK(20)
    
    vSK(:) = 0   !resets all Slater Koster parameters to zero
    
    ! tbModelType read in input file with tag 'TightBindingModelType'
    select case(tbModelType)

      !********************************************************************
      case(10)
      ! Lattice of atoms with s orbitals only n.n. hopping
      ! Works for any lattice in any dimension
      ! tbParams(1) is the value of the matrix element assumed to 
      !    apply to all neighbors at a distance < tbParams(2)

        ! Use tbParams read in input file after tag 'TightBindingParameters'
        if(distance .le. tbParams(2)) then
          vSK(1) = tbParams(1)            !t  - first neighbor
        end if

      !********************************************************************
      case(11)
      ! Lattice of atoms with s orbitals with n.n. and n.n.n. hopping
      ! See description of case 10

        ! Use tbParams read in input file after tag 'TightBindingParameters'
        if(distance .le. tbParams(2)) then
          vSK(1) = tbParams(1)            !t  - first neighbor
        else if(distance .le. tbParams(4)) then
          vSK(1) = tbParams(3)            !t' - second neighbor
        end if

     !******************************************************************
     case(12)
     ! Two site 1D lattice with s-p interactions only

        ! Use tbParams read in input file after tag 'TightBindingParameters'
        if(distance .le. tbParams(2)) then
          vSK(1) = tbParams(1)            !t  - first neighbor
        end if

     end select
     
  end Subroutine TBParamsSpecial
  
  
!RMM *************************************************

!*************************************************************************
  Subroutine TBVInit_Gen(V, F, K, xvec, distance, Lmin, Lmax, My )
!*************************************************************************
!
!   ***  THIS IS THE ONLY SUBROUTINE THAT CHANGES fOR SPECIFIC MOLDELS ***
!

!   WARNING: If the indices on the input and output arrays are removed,
!            then this subroutine may not give correct results!

!   NOTE: Do not remove the position of the debug lines.
!    

    implicit none

    integer, intent(in)           :: Lmin, Lmax

    ! Output variables
    complex(double), intent(out)  ::  V(Lmin**2:((Lmax+1)**2-1), &
                                      & Lmin**2:((Lmax+1)**2-1))  

    ! Input variables
    complex(double), intent(in)   :: My(Lmin**2:((Lmax+1)**2-1), &
                                      & Lmin**2:((Lmax+1)**2-1))

    real(double), intent(in)      :: F(Lmin:Lmax,Lmin:Lmax,0:Lmax)                                      
    real(double), intent(in)      :: K(Lmin:Lmax,Lmin:Lmax,0:Lmax)
    real(double), intent(in)      :: xvec(3)

    real(double), intent(in)      :: distance

    ! Workspace variables
    complex(double), pointer      :: Theta(:), Phi(:)
    complex(double), pointer      :: U(:,:), Temp1(:,:), Temp2(:,:)
  
    real(double), pointer         :: H_2c(:,:)
    real(double)                  :: rho
    real(double)                  :: costheta, sintheta
    real(double)                  :: cosphi, sinphi

    integer                       :: minOrb, maxOrb
    integer                       :: Lint,   LintPrime 
    integer                       :: LintSq, LintPrimeSq
    integer                       :: m, error

    ! All the relevant distances and angles you could possibly need
    ! for computing the two-center Slater-Koster parameters
    ! NOTE1: xvec is a unit vector.
    ! NOTE2: epsilon is used to handle the case of two atoms aligned
    !        along the z-axis without the evaluation of an 'if-then'
    !        statement
    rho = SQRT(xvec(1)**2 + xvec(2)**2)
    costheta = xvec(3)
    sintheta = rho
    cosphi = (xvec(1)+epsilon)/(rho+epsilon)
    sinphi = xvec(2)/(rho+epsilon)

    ! Initial and final indices into all the 2D matrices in this subroutine.
    minOrb = Lmin**2
    maxOrb = (Lmax + 1)**2 - 1

    ! Initialize Slater-Koster matrix
    V = CMPLX(zero, zero)

    ! Allocate the H_2c, Phi, Theta, U, Temp1 and Temp2 matrices
    allocate(H_2c(minOrb:maxOrb,minOrb:maxOrb), stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error allocating H_2c matrix.'

    allocate(Phi(minOrb:maxOrb), stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error allocating Phi matrix.'

    allocate(Theta(minOrb:maxOrb), stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error allocating Theta matrix.'
      
    allocate(U(minOrb:maxOrb,minOrb:maxOrb), stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error allocating U matrix.'

    allocate(Temp1(minOrb:maxOrb,minOrb:maxOrb), stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error allocating Temp1 matrix.'
      
    allocate(Temp2(minOrb:maxOrb,minOrb:maxOrb), stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error allocating Temp2 matrix.'

    ! Initialize the H_2c, Phi and Theta matrices
    ! H_2c = F*K
    ! Phi = exp(-i*phi*Lz)
    ! exp(-i*theta*Ly) = My*Theta*My^dag
    ! Temp1 = My*Theta
    ! Temp2 = My^dag*Phi
    H_2c(:,:)   = zero
    Phi(:)      = CMPLX(zero, zero)
    Theta(:)    = CMPLX(zero, zero)
    U(:,:)      = CMPLX(zero, zero)
    Temp1(:,:)  = CMPLX(zero, zero)
    Temp2(:,:)  = CMPLX(zero, zero)

    ! Create the Phi and Theta matrices. Also creates, the main diagonal
    ! of the H_2c matrix which corresponds to the hopping terms for
    ! Lprime = L.
    do Lint = Lmin, Lmax
       LintSq = Lint**2            
       do m = -Lint, Lint 
          H_2c(LintSq + Lint + m, LintSq + Lint + m) = &
          & distance**(F(Lint, Lint, abs(m)))*K(Lint, Lint, abs(m))
          Phi(LintSq + Lint + m) = (CMPLX(cosphi, sinphi))**(-m)
          Theta(LintSq + Lint + m) = (CMPLX(costheta, sintheta))**(-m)
       end do
    end do

    ! Create the off-diagonal parts of H_2c. These are the hopping terms
    ! for the case Lprime \neq = L.
    ! NOTE: The code below may be difficult to understand due to the 
    !       complicated indexing. However, it works and does not involve
    !       any conditional statements which would slow the code down.
    do Lint = Lmin, Lmax
       LintSq = Lint**2
       do LintPrime = Lint+1, Lmax
          LintPrimeSq = LintPrime**2
          do m = -Lint, Lint
             H_2c(LintPrimeSq + LintPrime + m, LintSq + Lint + m) = &
             & distance**(F(Lint, LintPrime, abs(m)))*K(Lint, LintPrime, abs(m))
             H_2c(LintPrimeSq + LintPrime - m, LintSq + Lint + m) = &
             & distance**(F(Lint, LintPrime, abs(m)))*K(Lint, LintPrime, abs(m))
             H_2c(LintSq + Lint + m, LintPrimeSq + LintPrime + m) = &
             & distance**(F(LintPrime, Lint, abs(m)))*K(LintPrime, Lint, abs(m))
             H_2c(LintSq + Lint - m, LintPrimeSq + LintPrime + m) = &
             & distance**(F(LintPrime, Lint, abs(m)))*K(LintPrime, Lint, abs(m))
          end do
       end do
    end do

    ! Debug
    ! write (1,*) "Phi & Theta (from TBVInit_Gen)" 
    ! write (1,*) "costheta = ", costheta
    ! write (1,*) "sintheta = ", sintheta
    ! write (1,*) "cosphi = ", cosphi
    ! write (1,*) "sinphi = ", sinphi

    ! Debug
    ! write (*,*) "H_2c (from TBVInit_Gen)"
    ! write (*,*) "i,  j,  H_2c(i,j)" 
    ! do Lint = minOrb, maxOrb
    !   do LintPrime = minOrb, maxOrb
    !      write(*,*) LintPrime, Lint, H_2c(LintPrime, Lint) 
    !   end do
    ! end do

    ! Debug
    ! write (2,*) "Phi and Theta (from TBVInit_Gen)"
    ! write (2,*) "i  Phi(i) Theta(i)" 
    ! do Lint = minOrb, maxOrb
    !    write(2,*) Lint, Phi(Lint), Theta(Lint)
    ! end do

    ! Create Temp1 = My*Theta and Temp2 = My^dag*Phi
    ! They have been "hand-coded" for efficiency, rather than using
    ! a matrix multiply which would waste a lot of time multiplying
    ! by zeroes.
    do Lint = minOrb, maxOrb
       Temp1(minOrb:maxOrb,Lint) = Theta(Lint)*My(minOrb:maxOrb,Lint)
       Temp2(minOrb:maxOrb,Lint) = Phi(Lint)*CONJG(My(Lint,minOrb:maxOrb))
    end do

    ! Debug
    ! write (1,*) "Temp1 (from TBVInit_Gen)"
    ! write (1,*) "i,  j,  Temp1(i,j)" 
    ! do Lint = minOrb, maxOrb
    !    do LintPrime = minOrb, maxOrb
    !       write (1,1000) LintPrime, Lint, Temp1(LintPrime, Lint)
    !    end do
    ! end do

    ! Debug
    ! write (1,*) "Temp2 (from TBVInit_Gen)"
    ! write (1,*) "i,  j,  Temp2(i,j)" 
    ! do Lint = minOrb, maxOrb
    !    do LintPrime = minOrb, maxOrb
    !       write (1,1000) LintPrime, Lint, Temp2(LintPrime, Lint)
    !    end do
    ! end do

    ! Create U = My*Theta*My^dag*Phi
    U = MATMUL(Temp1, Temp2)

    ! Debug
    ! write (*,*) "U (from TBVInit_Gen)"
    ! write (*,*) "i,  j,  U(i,j)" 
    ! do Lint = minOrb, maxOrb
    !   do LintPrime = minOrb, maxOrb
    !     write (*,1000) LintPrime, Lint, U(LintPrime, Lint)
    !   end do
    ! end do

    ! Create V = U^dag*H_2c*U
    ! NOTE: V is not a Hermitian matrix, because H_2c \neq H_2c\dag.
    V = MATMUL(MATMUL(TRANSPOSE(CONJG(U)),H_2c),U)

    ! Debug
    ! write (1,*) "V (from TBVInit_Gen)"
    ! write (1,*) "i,  j,  V(i,j)" 
    ! do Lint = minOrb, maxOrb
    !    do LintPrime = minOrb, maxOrb
    !       write (1,1000) LintPrime, Lint, V(LintPrime, Lint)
    !    end do
    ! end do
        
    ! Deallocate H_2c, Phi, Theta, U, Temp1 and Temp2 matrices.
    deallocate(H_2c, stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error deallocating H_2c matrix.'

    deallocate(Phi, stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error deallocating Phi matrix.'

    deallocate(Theta, stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error deallocating Theta matrix.'
      
    deallocate(U, stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error deallocating U matrix.'

    deallocate(Temp1, stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error deallocating Temp1 matrix.'
      
    deallocate(Temp2, stat = error)
    if (error /= 0) stop 'TBVInit_Gen: Error deallocating Temp2 matrix.'

  end Subroutine TBVInit_Gen

!*************************************************************************
  Subroutine LyEigenstates(Lmin, Lmax, My)
!*************************************************************************
! Calculates the eigenstates of Ly operator in the Lz basis.
! There are some peculiariatires of the eigenvectors produced by the
! LAPACK ZHBEV routine. Obviously, they are not exact. Do not be surprised
! to find an eigenvector with a component ~ 10^(-16) instead of zero. More
! bizarre is that some of the eigenvectors will have an arbitrary global
! phase factor. This doesn't really matter as long as all the eigenvectors
! are orthonormal, which they are indeed!

    implicit none

    integer, intent(in)           :: Lmin, Lmax
    complex(double), intent(out)  :: My(Lmin**2:((Lmax+1)**2 - 1), &
                                      & Lmin**2:((Lmax+1)**2 - 1))

! Error flag
    integer                       :: error

! Loop indices  
    integer                       :: Lint, m, m1, m2
    integer                       :: LintSq

! Parameters and work space for LAPACK ZHBEV routine
    character, parameter          :: JOBZ = 'V'
    character, parameter          :: UPLO = 'L'
    integer, parameter            :: KD = 1
    integer, parameter            :: LDAB = 2

    integer                       :: INFO, N, LDZ
    real(double)                  :: RWORK(6*Lmax+1), W(2*Lmax+1)
    complex(double)               :: Ly(2,2*Lmax+1), WORK(2*Lmax+1)
    complex(double)               :: Z(2*Lmax+1,2*Lmax+1)

! Temporary holding place for all Ly eigenstates (l=0 thru l=Lmax)
    complex(double)               :: Y(0:((Lmax+1)**2 - 1),0:((Lmax+1)**2 - 1))

! Make sure that Lmin and Lmax are positive.
    if (Lmin < 0) stop 'LyEigenstates: Lmin must be greater than zero'
    if (Lmax < 0) stop 'LyEigenstates: Lmin must be greater than zero'

! Make sure that Lmin < Lmax.
    if (Lmax < Lmin) stop 'LyEigenstates: Lmax must be greater than Lmin'

! Initialize Y.
    Y = CMPLX(zero, zero)

! Explictly do L = 0 case, since it is trivial.
    Y(0,0) = one

! Do other cases up to Lmax
    do Lint = 1, Lmax

! Initialize Ly matrix
	Ly = CMPLX(zero, zero)

! Construct Ly matrix from raising/lowering form
	do m = -Lint, Lint - 1
	   Ly(2,1+m+Lint) = CMPLX(zero, -half*SQRT(DBLE(Lint*(Lint+1)-m*(m+1))))
        end do ! m

! Diagonalize the Ly matrix using LAPACK ZHBEV routine
        LDZ = 2*Lmax + 1
	N   = 2*Lint + 1
	call ZHBEV(JOBZ, UPLO, N, KD, Ly, LDAB, W, Z, LDZ, WORK, RWORK, INFO)

	if (INFO /= 0) stop 'LyEigenstates: ZHBEV convergence error.'

! Put all the Ly eigenvectors in Y
        LintSq = Lint**2
	do m1 = 0, 2*Lint
	   do m2 = 0, 2*Lint
	      Y(LintSq + m1, LintSq + m2) = Z(m1+1, m2+1)
	   end do ! m2
        end do ! m1
    end do ! Lint

! Copy just the Ly eigenvectors of interest (l=Lmin thru l=Lmax)
    My = Y(Lmin**2:((Lmax+1)**2-1), Lmin**2:((Lmax+1)**2-1))

  end Subroutine LyEigenstates

!*************************************************************************
  Subroutine RealOrbitals(Lmin,Lmax,Cb)
!*************************************************************************
! Expands the real orbitals (s,px,py,...) in terms of the spherical
! harmonic. The columns of the matrix B are the expansion coefficients.
! If Lmin = 0 and Lmax = 2, the matrix is organized as follows:
!
! B(:,0) -> s               =              Y_{00}
! B(:.1) -> p_z             =              Y_{10}
! B(:,2) -> p_x             =   2^{-1/2} ( Y_{1-1} - Y_{11} )
! B(:,3) -> p_y             = i*2^{-1/2} ( Y_{1-1} + Y_{11} )
! B(:,4) -> d_{3z^2 - r^2}  =              Y_{20}
! B(:,5) -> d_{xz}          =   2^{-1/2} ( Y_{2-1} - Y_{21} )
! B(:,6) -> d_{yz}          = i*2^{-1/2} ( Y_{2-1} + Y_{21} )
! B(:,7) -> d_{x^2-y^2}     =   2^{-1/2} ( Y_{2-2} + Y_{22} )
! B(:,8) -> d_{xy}          = i*2^{-1/2} ( Y_{2-2} - Y_{22} )
!
! The matrix B for this case is given by
!     [1    0    0    0    0    0    0     0     0  ]
!     [0    0   0.7  0.7i  0    0    0     0     0  ]
!     [0    1    0    0    0    0    0     0     0  ]
!     [0    0  -0.7  0.7i  0    0    0     0     0  ]
! B = [0    0    0    0    0    0    0    0.7   0.7i]
!     [0    0    0    0    0   0.7  0.7i   0     0  ]
!     [0    0    0    0    1    0    0     0     0  ]
!     [0    0    0    0    0  -0.7  0.7i   0     0  ] 
!     [0    0    0    0    0    0    0    0.7  -0.7i]
!
! Need the coefficients    
    implicit none

    integer, intent(in)           :: Lmin, Lmax
    complex(double), intent(out)  :: Cb(Lmin**2:((Lmax+1)**2 - 1),Lmin**2:((Lmax+1)**2 - 1))

! Constants
    real(double), parameter       :: invsqrt2 = 0.70710678118655D0

! Loop indices  
    integer                       :: Lint, LintSq, m, mprime

! Temporary hold place for all the expansion coefficients (l=0 thru l=Lmax)
    complex(double)               :: B(0:((Lmax+1)**2 - 1),0:((Lmax+1)**2 - 1))

! Make sure that Lmin and Lmax are positive.
    if (Lmin < 0) stop 'RealOrbitals: Lmin must be greater than zero'
    if (Lmax < 0) stop 'RealOrbitals: Lmax must be greater than zero'

! Make sure that Lmin < Lmax.
    if (Lmax < Lmin) stop 'RealOrbitals: Lmax must be greater than Lmin'

! Initialize B.
  B = CMPLX(zero, zero)

! There are three cases to take care of.
!
! Case 1: Real orbitals which are given by Y_{l0}.
  do Lint = 0, Lmax
	LintSq = Lint*Lint
	B(LintSq + Lint, LintSq) = one
  end do

! Case 2: Real orbitals which are given by (Y_{l-m} + (-1)^m Y_{lm}^\ast)/sqrt(2),
!         where m is a positive integer.
  do Lint = 1, Lmax
	LintSq = Lint*Lint
	do m = 1, Lint
	   mprime = 2*(m-1) + 1
           B(LintSq + Lint - m, LintSq + mprime) = CMPLX(invsqrt2, zero)
           B(LintSq + Lint + m, LintSq + mprime) = CMPLX(((-1)**m)*invsqrt2, zero)
	end do ! m
  end do ! Lint

! Case 3: Real orbitals which are given by i*(Y_{l-m} - (-1)^m Y_{lm}^\ast)/sqrt(2),
!         where m is a positive integer.
 do Lint = 1, Lmax
	LintSq = Lint*Lint
	do m = 1, Lint
	   mprime = 2*m
           B(LintSq + Lint - m, LintSq + mprime) = CMPLX(zero, invsqrt2)
           B(LintSq + Lint + m, LintSq + mprime) = CMPLX(zero, -((-1)**m)*invsqrt2)
	end do ! m
  end do ! Lint

! Copy just the expansion coefficients of interest (l=Lmin thru l=Lmax).
  Cb = B(Lmin**2:((Lmax+1)**2-1), Lmin**2:((Lmax+1)**2-1))

  end Subroutine RealOrbitals

!*************************************************************************
  Function TBHMatElem_Gen( VSlaterKoster, bra, ket)
!*************************************************************************
! Calculates the geometric factors for the two center tight binding 
! hamiltonian matrix elements, where x = vector between atoms.


    implicit none

    real(double)                :: TBHMatElem_Gen

    complex(double), intent(in) :: VSlaterKoster(:,:)
    complex(double), intent(in) :: bra(:), ket(:)
    
!  Parameters and workspace for BLAS ZGEMV routine
    character, parameter        :: TRANS = 'N'   
    complex(double), parameter  :: ALPHA = (1.d0, 0.d0)
    complex(double), parameter  :: BETA  = (0.d0, 0.d0)
    integer, parameter          :: INCX  = 1
    integer, parameter          :: INCY  = 1

    complex(double), pointer    :: Y(:)
    integer                     :: N, M, LDA
    integer                     :: error

!  Debug
    integer                     :: i
    

!  From Koster tables - matrix elements easily written in terms
!  of direction cosines given by the xvect unit vector
!  See Harrison, "Elec. Str. and Props. of Solids", p. 481 
!  We require that for l = 1, m =1 for px; m=2 for py; m=3 for pz
!  where x,y,z are the 1,2,3 cartesian directions defined in the 
!  structure


    ! Parameters for matrix vector product using BLAS ZGEMV
    ! NOTE: VSlaterKoster is not a Hermitian matrix.
    N = INT(SQRT(REAL(size(VSlaterKoster))))
    M = N
    LDA = N

    ! Allocate the Y vector
    allocate(Y(1:N), stat = error)
    if (error /= 0) stop 'TBHMatElem_Gen: Error allocating Y vector.'

    ! VSlaterKoster|ket>
    call ZGEMV(TRANS, M, N, ALPHA, VSlaterKoster, LDA, ket, INCX, &
             & BETA, Y, INCY)   

    ! Debug
    ! write (*,*) "bra & ket (from TBHMatElem_Gen)"
    ! write (*,*) "i,  bra(i)  ket(i)"
    ! do i = 1, N
    !   write(*,*) i, bra(i), ket(i)
    ! end do
    
    ! <bra|VSlaterKoster|ket>
    ! NOTE: Calling ZDOTC instead of DOT_PRODUCT seems to be 1% faster
    !       This is most likely compiler and machine dependent.
    TBHMatElem_Gen = dot_product(bra,Y)
    ! TBHMatElem_Gen = ZDOTC(N, bra, INCX, Y, INCY)

    deallocate(Y, stat = error)
    if (error /= 0) stop 'TBHMatElem_Gen: Error deallocating Y vector.'

  end Function TBHMatElem_Gen


end Module SlaterKosterParamMod











