!*************************************************************************
Program TB
!*************************************************************************
! TB calculates the band structure along various directions in
! k-space using empirical tight binding 2-center form.   
! Options include a general form of arbitrary angular momenta
! and forms using Harrison's universal matrix elements.
! Several simple models in 1,2 and 3 d are included; others can be added
! by editing one subroutine
!
! Written by Nichols Romero and Richard M. Martin of the Univ. of Illinois, 
! based upon programs written by William Mattson and Y. Kwon.
! 
! Note: atomic units (a.u.) are used throughout the program.

  use sysParams,           only : double, zero
  use tagHandlerMod,       only : tagHandlerT, tagHandlerInit, tagHandlerDestroy 
  use structureMod,        only : structureT, structureInit, structureDestroy
  use kPointsMod,          only : kPointsT, kPointsInit, kPointsDestroy
  use graphMod,            only : graphT, PlotBands
  use hamiltonianMod,      only : hamiltonianArrayT, HInit, HDestroy, &
                                  & HDiagonalize,  HPrint
  use eigenStatesMod,      only : eigenStatesT, eigenStatesInit, & 
                                  & eigenStatesDestroy, PrintEigenvalues
                                  
  !*** Specific to the TB method
  use tbHamMod,            only : hamInfoT, hamInfoInit, hamInfoDestroy, &
  				  & subHamiltonianArrayT, subHInit, subHDestroy, &
                                  & subHGenerate, HGenerate

  implicit none

  type( tagHandlerT ),          pointer :: tagHandler
  type( hamInfoT ),             pointer :: hamInfo
  type( StructureT ),           pointer :: structure1
  type( kPointsT ),             pointer :: kPoints
  type( eigenStatesT ),         pointer :: eigenStates
  type( graphT ),               pointer :: graphinfo
  type( hamiltonianArrayT ),    pointer :: hamiltonian

  !*** Specific to the TB method
  type( subHamiltonianArrayT),  pointer :: subHamiltonian

  ! Workspace variables
  integer         :: kpt
  integer         :: error
  integer         :: ReportNo
  character(80)   :: filename

  print *, ' Input file name' 
  read( 5, '(A80)') filename

  ! Initialize "tag handler" to read the input from file 'filename'.
  call TagHandlerInit( tagHandler, 10, filename)

  ! This will call DimensionsInit, LatticeInit, and AtomsInit.
  call StructureInit( structure1, tagHandler ) 

  call KPointsInit( kPoints, graphInfo, structure1%ndim, structure1%lattice, tagHandler )
  
  !*** Specific to the TB method
  ! Reads on-site energies, sets up neighbors, orbit info, etc.
  call HamInfoInit( hamInfo, structure1, tagHandler )
  
  ! Allocate and generate the sub-Hamiltonian array: subH(orbitJ, orbitI, neigh, atomI)
  ! Array of hamiltonian matrix elements in real space
  call SubHInit( subHamiltonian, maxval(hamInfo%orbits), maxval(hamInfo%numChannels), &
               & maxval(hamInfo%neighbors%neighborList(:)%numNeighbors), &
               & structure1%atomBasis%totalNumAtoms ) 
  call subHGenerate(hamInfo, structure1, subHamiltonian%subH)
  
  ! Allocate array for H and work arrays of dimension numTotalOrbits
  Print *, 'Allocate array for H and work arrays - dimension = ', &
  & hamInfo%numTotalOrbits
  call HInit( hamiltonian, hamInfo%numTotalOrbits)
 
  ! Allocate array for eigenStates array
  call EigenStatesInit( eigenStates, kPoints%numBands, kPoints%numKPoints, &
                      & hamiltonian%hDim, -1.0d0, 1)
  eigenStates%eigenVectors = zero
 
  ! Calculate the Hamiltonian matrix at each k-pt as sum of
  ! phase factors exp(i k xvec) * subH(orbitJ, orbitI, neigh, atomI)
  ! and diagonialize
  do kpt = 1, kPoints%numKPoints
    call HGenerate( hamInfo, structure1, kpt, kPoints, &
                  & subHamiltonian%subH, hamiltonian%h)
    call HPrint( hamiltonian ) 
    call HDiagonalize( hamiltonian, eigenStates, kpt )
  end do ! kpt
  
  Print *, 'Eigenvalues calculated and written to file plot.dat'
  
  ! Create files plot.dat and gnuplot.dat (Display with: gnuplot gnuplot.dat)
  if( associated( graphinfo ) )    call PlotBands( graphinfo, eigenStates )           

  if( .not. associated( graphinfo ) )  call PrintEigenvalues(eigenStates, 1)      

  !  Clean up memory
  call EigenStatesDestroy( eigenStates )
  call SubHDestroy ( subHamiltonian)
  call HDestroy( hamiltonian )                
  call HamInfoDestroy( hamInfo )    
  call KPointsDestroy( kPoints )
  call StructureDestroy( structure1 )

  ! Clean up and close the input file
  call TagHandlerDestroy( tagHandler )

end Program TB
