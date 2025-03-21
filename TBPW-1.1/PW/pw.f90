!*************************************************************************
Program PW
!*************************************************************************
! PW calculates the band structure along various directions in
! k-space using a plane-wave basis and a fixed "effective" potential.  
! Several choices of empirical pseudopotentials are provided.
! Written by William Mattson and Richard M. Martin of the University 
! of Illinois, based upon program written by K. Glassford and I. Souza.

! Note: atomic units (a.u.) are used throughout the program.

  use sysParams,           only : double, zero
  use TagHandlerMod,       only : tagHandlerT, TagHandlerInit, TagHandlerDestroy 
  use StructureMod,        only : StructureT, StructureInit, StructureDestroy
  use kPointsMod,          only : kPointsT, KPointsInit, KPointsDestroy
  use graphMod,            only : graphT, PlotBands
  use hamiltonianMod,      only : hamiltonianArrayT, HInit, HDestroy, HDiagonalize,  & 
                                  & HPrint
  use eigenStatesMod,      only : eigenStatesT, EigenStatesInit, EigenStatesDestroy, &
                                  & PrintEigenvalues
  use pwHamMod,            only : hamInfoT, HamInfoInit, HamInfoDestroy, &
                                  & HGenerate, FindHMaxSize, GridInit, &
                                  & PotentialRGenerate, KineticGenerate
  use cgParamsMod,         only : cgParamsT, cgParamsInit
  use typeMod,             only : GridT
  use ConjGradMod,         only : cgEigenSystem,cgFFTDiagonalise

  implicit none

  type(tagHandlerT),       pointer :: tagHandler
  type(hamInfoT),          pointer :: hamInfo
  type(StructureT),        pointer :: structure1
  type(kPointsT),          pointer :: kPoints
  type(eigenStatesT),      pointer :: eigenStates
  type(hamiltonianArrayT), pointer :: hamiltonian
  type(graphT),            pointer :: graphinfo
  type(cgParamsT)                  :: cgParams
  type(GridT),             pointer :: Grid

  integer :: kpt, r, ReportNo
  character(80) :: filename

  ReportNo = 110

  print *, ' Input file name' 
  read( 5, '(A80)') filename

  ! Necessary if using conjugate gradient method
  call Random_Seed

  ! Initialize "tag handler" to read the input from file 'filename'.
  call TagHandlerInit( tagHandler, 10, filename)

  ! Read the Conjugate Gradient Parameter
  call cgParamsInit(cgParams,tagHandler)

  ! This will call DimensionsInit, LatticeInit, and AtomsInit.
  call StructureInit( structure1, tagHandler )

  call KPointsInit( kPoints, graphInfo, structure1%ndim, structure1%lattice, tagHandler)

  ! Gvector list and structurefactor
  call HamInfoInit( hamInfo, structure1, tagHandler )

  ! Find hMaxSize = max size of H at any k point
  call FindHMaxSize( hamInfo, structure1, kPoints )
  
  ! Allocates array for H and work arrays
  call HInit( hamiltonian, hamInfo%hMaxSize )

  ! Allocates array for eigenStates array
  call EigenStatesInit( eigenStates, kPoints%numBands, kPoints%numKPoints, &
                      & hamInfo%hMaxSize, -1.0d0, 1)
  eigenStates%eigenVectors = zero

  ! Initialise the FFT Grid and Plans
  Call GridInit(Grid,hamInfo,eigenStates,structure1%ndim)

  open(unit=ReportNo,file='Report.cg',status='replace',action='write', &
       & position='rewind')
  write(ReportNo,*) 'Eigenvalues:'
  write(ReportNo,*) 'The FFT Grid Size +/-',Grid%Size(1)
  write(ReportNo,*) 'The number of GVevtors:',Grid%NumgVec
  write(ReportNo,*) 'Eigenvalues:'

  ! Calculate the Potential and FFT to Vr
  Call PotentialRGenerate(hamInfo,Grid)

  if(cgParams%Switch.eq.1) then 
     do kpt = 1, kPoints%numKPoints 
        call HGenerate( hamInfo, structure1, kPoints, eigenStates, &
             & kpt, hamiltonian%h, hamiltonian%sqH, hamiltonian%hSize )

        write(ReportNo,*) 'k:',kpt

        call cgEigenSystem( hamiltonian,eigenStates,kpt,cgParams%tol, &
             & cgParams%period,ReportNo)

        print *, kpt,'completed out of a total of ',kPoints%numKPoints

     end do ! kpt
  else
     do kpt = 1, kPoints%numKPoints

        call HGenerate( hamInfo, structure1, kPoints, eigenStates, &
             & kpt, hamiltonian%h, hamiltonian%sqH, hamiltonian%hSize )

        write(ReportNo,*) 'k:',kpt
     
        call HDiagonalize( hamiltonian, eigenStates, kpt )

        print *, kpt,'completed out of a total of ',kPoints%numKPoints
        
     end do ! kpt

  end if

    close(unit=ReportNo)

  ! Creates files plot.dat and gnuplot.dat (Display with: gnuplot gnuplot.dat)
  if( associated( graphinfo ) )    call PlotBands( graphinfo, eigenStates )           

  if( .not. associated( graphinfo ) )  call PrintEigenvalues(eigenStates, kPoints%numKPoints)    

  ! Clean up memory
  call EigenStatesDestroy( eigenStates )

  call HDestroy( hamiltonian )
                
  call HamInfoDestroy( hamInfo )

  call KPointsDestroy( kPoints )

  call StructureDestroy( structure1 )

  ! Clean up and close the input file
  call TagHandlerDestroy( tagHandler )

end Program PW






