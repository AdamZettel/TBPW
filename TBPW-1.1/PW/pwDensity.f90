!*************************************************************************
Program PWDensity
!*************************************************************************
! PW calculates the band structure along various directions in
! k-space using a plane-wave basis and a fixed "effective" potential.  
! Several choices of empirical pseudopotentials are provided.
! Written by William Mattson and Richard M. Martin of the University 
! of Illinois, based upon program written by K. Glassford and I. Souza.

! Note: atomic units (a.u.) are used throughout the program.

  use sysParams,           only : double, zero
  use TagHandlerMod,       only : TagHandlerInit, TagHandlerDestroy, tagHandlerT
  use StructureMod,        only : structureT, StructureInit, StructureDestroy
  use KPointsMod,          only : MonkHorstPackinit, MonkHorstPackDestroy, &
                                & kPointsT
  use HamiltonianMod,      only : hamiltonianArrayT, HInit, HDestroy, & 
                                & HDiagonalize, HPrint
  use EigenStatesMod,      only : eigenStatesT, EigenStatesInit, EigenStatesDestroy
  use GVectorsMod,         only : gVectorObjectsT, GVectorsInit, GVectorsDestroy
  use pwHamMod,            only : FindHMaxSize, CalculateDensity, HGenerate, &
                                & HamInfoT, HamInfoInit, HamInfoDestroy
  Use DensityArrayMod,     only : DensityArrayT, &
                                & DensityArrayInit, DensityArrayPlot, &
                                & DensityArrayDestroy
			
  implicit none

  type(HamInfoT),          pointer :: HamInfo
  type(structureT),        pointer :: structure1
  type(kPointsT),          pointer :: kPoints
  type(tagHandlerT),       pointer :: tagHandler
  type(gVectorObjectsT),   pointer :: gVectors
  type(eigenStatesT),      pointer :: eigenStates
  type(hamiltonianArrayT), pointer :: hamiltonian
  type(DensityArrayT),     pointer :: DensityArray


  integer       :: kpt
  character(80) :: filename

  print *, 'Input File Name'
  read(5,'(A80)') filename

  ! Initialize "tag handler" to read the input from file 'filename'.
  call TagHandlerInit( tagHandler, 10, filename)

  ! This will call DimensionsInit, LatticeInit, and AtomsInit.
  call StructureInit( structure1, tagHandler )

  ! Gvector list and structurefactor
  call GVectorsInit( gVectors, structure1, tagHandler )	

  call MonkHorstPackInit( kPoints, structure1%ndim, structure1%lattice, tagHandler)

  Call DensityArrayInit( DensityArray, structure1%ndim, 100)

  Call HamInfoInit( hamInfo, structure1, TagHandler )

  ! Find hMaxSize = max size of H at any k point
  call FindHMaxSize( hamInfo, structure1, kPoints )
  
  ! Allocates array for H and work arrays
  call HInit( hamiltonian, HamInfo%hMaxSize)

  call EigenStatesInit( eigenStates, kPoints%numBands, kPoints%numKPoints, HamInfo%hMaxSize, -1.0d0, 1 )

  eigenStates%eigenVectors = zero

  ! Creates files plot.dat and gnuplot.dat (Display with: gnuplot gnuplot.dat)
  do kpt = 1, kPoints%numKPoints

    call HGenerate( hamInfo, structure1, kPoints, eigenStates, &
         & kpt, hamiltonian%h, hamiltonian%sqH, hamiltonian%hSize )

    call HDiagonalize( hamiltonian, eigenStates, kpt )

        print *, kpt,'completed out of a total of ',kPoints%numKPoints
        
  end do ! kpt
 
  call CalculateDensity( hamInfo, structure1, kPoints, eigenStates, densityArray )

  ! Clean up memory
  call EigenStatesDestroy(eigenStates)

  call HDestroy(hamiltonian)
                
  call MonkhorstPackDestroy( kPoints )

  call GVectorsDestroy( gVectors )

  call StructureDestroy( structure1 )

  ! Clean up and close the input file
  call TagHandlerDestroy( tagHandler )

end Program PWDensity












