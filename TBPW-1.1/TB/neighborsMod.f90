!*************************************************************************
Module NeighborsMod
!*************************************************************************
!
! The neighbors module containds information about the neighbors of the
! unit cell atoms.  These neighbors lie within the cut off radius 
! maxDistance.  
!
! The module contains the public functions:
!    NeighborsInit       Reads in maxDistance, calls FindNeighborCells then
!                        calls FindNeighbors.
!    NeighborsDestroy    Deallocates space allocated by FindNeighborCells and
!                        FindNeighbors.
!    FindNeighborCells   Calculates how the unit cell needs to be replicated
!                        to include all neighbors with in maxDistance.  It
!                        then stores the indices of these replicated cells
!                        in space allocated in neighborCells and the number
!                        in numCells.
!    FindNeighbors       For each atom i in the unit cell it finds the
!                        neighbors of that atom first in the unit cell and
!                        then in the replicated cell.  When a neighbor j is
!                        found an entry is made in 
!                        neighborList(i).neighbor(j) containing the distance,
!                        unit vector to, and type of the neighbor j.  For 
!                        symetry considerations only atoms with atom number
!                        greater than or equal to i are stored as neighbors.
!    PrintNeighbors      Prints a list of the neighbors.

  use SysParams, only : double, maxDimensions
  use StructureMod, only : structureT
  implicit none


  ! Type used to store data about a neighbor
  type oneNeighborInfoT
    integer      :: id               ! The number of the corresponding atom
                                     ! in the unit cell.
    real(double) :: x(maxDimensions) ! The unit vector in the direction of 
                                     ! the neighbor.
    real(double) :: distance         ! The distance to the neighbor.
  end type oneNeighborInfoT

  ! Type used to make a list of neighbors for an atom.
  type oneNeighborListT
    integer				:: numNeighbors  ! The number of neighbors.
    type(oneNeighborInfoT), pointer     :: neighbor(:)   ! An array of the neighbors.
  end type oneNeighborListT

  ! Type that describes all neighbors of each atom in the the cell
  type neighborsT
    real(double) :: maxDistance, maxDist2     ! cut off radius and it's square.
    integer, pointer :: atomIDToType(:)       ! Table listing the species number of atom i.
    integer, pointer :: neighborCells(:,:)    ! List of indices for neighboring cells.
    integer              :: numCells          ! The number of neighboring cells.
    type( oneNeighborListT ), pointer :: neighborList(:) ! A set of neighbor lists.
    type( StructureT ),  pointer :: struct
  end type neighborsT



contains

!*************************************************************************
  Subroutine FindNeighborCells( neighbors )
!*************************************************************************
!  FindNeighborCells finds the replicas of the unit cells nessecary to
!  contain all of the atoms with in the distance maxDistance.  It stores
!  the indices of these cells in neighborCells, with a total of numCells.
!  This function simply takes the parallel piped that contains the volume
!  of the sphere of the cut off radius and store all of the cells in that
!  parallel piped.
    use SysParams,    only : twopi

    implicit none

    type( neighborsT ), pointer :: neighbors

    integer :: i, j, k, minCellNum(3), error, cellTrial(3)

    ! Calculate the minimum number of cells along each lattice vector to 
    ! contain all of the neighbors.
    minCellNum = 0
    minCellNum(1:neighbors%struct%ndim) = &
        & int(neighbors%maxDistance *neighbors%struct%lattice%bLatVecLen / twopi) + 1

    ! Allocate the space needed to store the cell indices.
    allocate(neighbors%neighborCells(neighbors%struct%ndim,&
            & product(2*minCellNum+1)), stat = error)
    if(error /= 0) stop 'Error allocating neighborCells in FindNeighborCells.'

    neighbors%numCells = 0
    neighbors%neighborCells = 0

    ! Loop over all of the cells in the parallel piped and store the indices.
    do i = -minCellNum(1), minCellNum(1)
      do j = -minCellNum(2), minCellNum(2)
        do k = -minCellNum(3), minCellNum(3)

            cellTrial = (/i, j, k/)
            neighbors%numCells = neighbors%numCells + 1
            neighbors%neighborCells(:,neighbors%numCells) = &
                & cellTrial(1:neighbors%struct%ndim)

        end do
      end do
    end do

  end Subroutine FindNeighborCells

!*************************************************************************
  Subroutine FindNeighbors( neighbors )
!*************************************************************************
! FindNeighbors finds all of the neighbors with in the maxDistance of the
! atom for each atom in the unit cell. This is done by looping over each
! atom i in the unit cell, and the every other atom j in the unit cell 
! after the atom i(for reasons of symmetry) if the atom j is with in range
! information is stored about it in neighborList(i).neighbors.  The every
! replicated cell in neighborCells is looped over and the same calculation
! for every atom as a possible neighbor in that cell is perfomed.

    implicit none

    type( neighborsT ), pointer :: neighbors
    integer :: error, iElement, iAtom, iCell, &
             & jAtom, iDimension
    real(double) :: x(3), distanceIJ
    type(oneNeighborListT) :: localNeighbor

    ! Allocate temporary space to store neighbors for a single atom.
    allocate(localNeighbor%neighbor(neighbors%struct%atomBasis%totalNumAtoms * neighbors%numCells), stat = error)
    if(error /= 0) stop 'Error allocating localNeighbor in NeighborsInit.'

    ! Allocate space to allocate a neighbor list for each atom.
    allocate(neighbors%neighborList(neighbors%struct%atomBasis%totalNumAtoms), stat = error)
    if(error /= 0) stop 'Error allocating neighborList in NeighborsInit.'

    ! Allocate space for the table.
    allocate(neighbors%atomIDToType(neighbors%struct%atomBasis%totalNumAtoms), stat = error)
    if(error /= 0) stop 'Error allocating localNeighbor in NeighborsInit.'

    jAtom = 0

    ! Create the table that stores the element number for every atom.
    do iElement = 1, neighbors%struct%atomBasis%numSpecies
      do iAtom = 1, neighbors%struct%atomBasis%numAtomsOfSpecies(iElement)

        jAtom = jAtom + 1
    	neighbors%atomIDToType( jAtom ) = iElement

      end do
    end do

    ! Now find neighbors for each atom iAtom in the unit cell.
    do iAtom = 1, neighbors%struct%atomBasis%totalNumAtoms

      localNeighbor%numNeighbors = 0
        
      ! Now loop over all of the neighboring cells to find other neighbors.
      do iCell = 1, neighbors%numCells
        ! Loop over each atom in the neighboring cell with number greater
        ! than or equal to iAtom.
        do jAtom = iAtom, neighbors%struct%atomBasis%totalNumAtoms

          x = 0

          ! Calculate the position of jAtom in the neighboring cell.
          do iDimension = 1, neighbors%struct%ndim

            x = x + neighbors%neighborCells(iDimension,iCell) * &
              & neighbors%struct%lattice%aLatVec(:,iDimension)

          end do

          ! Calculate the distance from iAtom to the replica of jAtom.
          x = neighbors%struct%atomBasis%atomCart(:,iAtom)-neighbors%struct%atomBasis%atomCart(:,jAtom)-x
          distanceIJ = dot_product(x,x)
      
          ! If the distance is less than maxDistance the replica of jAtom
          ! is a neighbor.
          if(distanceIJ .le. neighbors%maxDist2 .and. distanceIJ .ne. 0.0) then

            ! Store information about the neighbor in the local list.
            localNeighbor%numNeighbors = localNeighbor%numNeighbors + 1
            distanceIJ = sqrt(distanceIJ)
            localNeighbor%neighbor(localNeighbor%numNeighbors)%id = jAtom
            localNeighbor%neighbor(localNeighbor%numNeighbors)%x = &
                        & x / distanceIJ
            localNeighbor%neighbor(localNeighbor%numNeighbors)%distance = &
                        & distanceIJ
          end if
        end do
      end do

      ! Allocate space for the neighbor list of iAtom.
      allocate(neighbors%neighborList(iAtom)%neighbor(localNeighbor%numNeighbors), &
             & stat = error)
      if(error /= 0) stop 'Error allocating neighborList in NeighborsInit.'

      ! Copy the information in the local list to the modules neighborList.
      neighbors%neighborList(iAtom)%numNeighbors = localNeighbor%numNeighbors
      neighbors%neighborList(iAtom)%neighbor(:) = &
        & localNeighbor%neighbor(1:localNeighbor%numNeighbors)

    end do
  end Subroutine FindNeighbors

!*************************************************************************
  Subroutine NeighborsInit( neighbors, struct, tagHandler )
!*************************************************************************
! Do not call unless LatticeInit has already been called.
! Neighbors Init reads in the maximumDistance from the input file
! and calls FindNeighborCells and FindNeighbors to find all of the
! neighboring atoms.

    use TagHandlerMod, only : FindTag, tagHandlerT

    implicit none

    type( neighborsT ),  pointer :: neighbors
    type( StructureT ),  pointer :: struct
    type( tagHandlerT ), pointer :: tagHandler

    integer :: error

    if( associated( neighbors )) stop 'Error: neighbors already allocated.'
    allocate( neighbors, stat = error )
    if( error /= 0 ) stop 'Error: neighbors allocation failed.'

    neighbors%struct => struct

    !*************************************************************************
    ! Read in the Maximum distance print it's value, and calculate it's
    ! square for later use.
	if(struct%ndim .gt. maxDimensions .or. struct%ndim .lt. 1) &
      & stop 'Error, neighbors requires three dimensions'

    call FindTag( tagHandler, 'MaximumDistance', error)

	if(error .eq. 0) then

	  read(unit = tagHandler%fileno, fmt = '(G15.10)', iostat = error) neighbors%maxDistance
	  if(error /= 0) stop 'Error, couldn''t read MaximumDistance.'

	else

	  stop 'Error couldn''t find MaximumDistance.'

	end if
    
    print *, 'MaximumDistance ', neighbors%maxDistance

	neighbors%maxDist2 = neighbors%maxDistance**2
    ! maxDistance and maxDist2 are initialized.
    !*************************************************************************

    !*************************************************************************
    ! Find all the neighbors and initialize the neighborList
    call FindNeighborCells( neighbors )
    call FindNeighbors( neighbors )
    ! NeighborList and neighbors%atomIDToType initialized.
    !*************************************************************************

  end Subroutine NeighborsInit

!*************************************************************************
  Subroutine NeighborsDestroy( neighbors )
!*************************************************************************
! Deallocate the neighborList and atomIDToType.

    implicit none

    type( neighborsT ), pointer :: neighbors
    integer :: error, i

    if( .not. associated( neighbors )) stop 'Error: neighbor not allocated.'

    deallocate(neighbors%atomIDToType, neighbors%neighborCells, stat = error)
    if(error /= 0) stop 'Error deallocating atomIDToType in neighborsDestroy.'

    do i = 1, neighbors%struct%atomBasis%totalNumAtoms

      deallocate(neighbors%neighborList(i)%neighbor, stat = error)
      if(error /= 0) stop 'Error deallocating neighbor in neighborsDestroy.'

    end do

    deallocate(neighbors%neighborList, stat = error)
    if(error /= 0) stop 'Error deallocating neighborList in neighborsDestroy.'

    deallocate(neighbors, stat = error)
    if(error /= 0) stop 'Error: deallocating neighbors failed.'

    nullify( neighbors )

  end Subroutine NeighborsDestroy

!*************************************************************************
  Subroutine PrintNumNeighbors( neighbors )
!*************************************************************************
! Print the neighborList for debugging.

    implicit none
    
    type( neighborsT ), pointer :: neighbors
    integer :: i
    
    do i = 1, neighbors%struct%atomBasis%totalNumAtoms
    
      print *, neighbors%struct%atomBasis%elementLabel(neighbors%atomIDToType(i)),' #', i,' has ', &
             & neighbors%neighborList(i)%numNeighbors, ' neighbors.'
    
    end do

  end Subroutine PrintNumNeighbors

end Module NeighborsMod








