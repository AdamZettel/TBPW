!*************************************************************************
Module tbHamMod
!*************************************************************************
!
! The generic name of the module "tbHamMod" and the subroutines
! "HGenerate", "HamInfoInit", "HamInfoDestroy",  must be used so they are 
! called properly by the general hamiltonian routines
!
! The tbHamMod  Module contains the routines and data
! nessecary to generate a tight binding hamiltonian.  All subroutines are
! universal for 2-center tight binding except TBVInit_Gen which determines the
! Slater Koster tight binding interactions between different atoms. Site
! energies are read from the input file.
!
! This module contains the public functions:
!
!    HamInfoInit          Sets up data necesary for the tight
!                         binding hamiltonian.  It uses Lattice,
!                         and KpointsInit, and Hinit.  It reads
!                         in the OrbitsAndEnenrgies and orbit types l,m
!                         for each element.  It finds the total no. orbits
!                         generates an index for the hamiltonian given the
!                         two atom numbers and their orbits.
!
!    HamInfoDestroy       Deallocates space allocated by HamInfoInit
!
!    SubHInit         Allocates space for subH - array of matrix elements
!             of the TB hamiltonian in real space
!
!    SubHDestroy      Deallocates space allocated by SubHInit
!
!    subHGenerate     Generates the matrix elements
!             of the TB hamiltonian in real space subH 
!                         by calling TBVInit_Gen and TBHMatElem_Gen for each atom
!                         and it's neighbors. 
!
!    HGenerate            Generates the hamiltonian for kpoint k, 
!                         by calculating the pahse factors and using the matrix
!             elements precomputed by subHGenerate and stored in the array subH

  use SysParams,    only : double, zero, one
  use NeighborsMod, only : neighborsT

  implicit none

  !******************************************************************************
  type hamInfoT
    !****************************************************************************
    ! Information about the Tight Binding Model
 
    ! NOTE: A '*' denotes that the index on a certain array has been set to 
    ! maximize speed throughput. This mean that if you change the order of the
    ! indices you risk missing the cache at execution time.

    integer          :: tbModelType      ! Model Type  
                         ! 0 = General 2-center Model valid for any Ang. Mom.
                         ! 1 = Harrison universal 2-center
                         ! >=10 = Special models coded in TBParamsSpecial

    integer          :: numTotalOrbits   ! Size of the hamiltonian.

    integer, pointer :: numchannels(:)   ! The number of electron channels for a species
                                         ! Allows different principle quantum numbers
                                         ! Similar the multiple zeta orbitals with same l,m
                                         ! Must be 1 for Harrison model
                                         ! numchannels(species) 

    integer, pointer :: orbits(:,:)      ! The total number of orbits
                                         ! per element per electron channel. (*)
                                         ! orbits(n,species)

    integer, pointer :: orbittype(:,:,:,:)   ! orbittype(1,orbit,n,species)
                                             ! contains the value of ell.(*)
                                             ! orbittype(2,orbit,n,species)
                                             ! contains the value of mstar.(*)
 
    real(double), pointer :: energies(:,:,:) ! These are the orbital energies
                                             ! which are part of the contri-
                                             ! bution to the diagonal terms.(*)
                                             ! energies(orbit,n,species)

    integer, pointer :: hIndex(:,:,:,:,:,:)  ! Stores the index into the 
                                             ! Hamiltonian.(*)
                                             ! hIndex(orbitJ,orbitI,nj,ni,atomJ,atomI)
    !**The following are NOT used by the general routine that works for any angular momentum
     
    real(double) 	:: vSK(20)          ! Array for SK matrix elements
                                             ! used in Harrison and simple TB models

    integer             :: numTBParams       ! Number of tbParams to be read in

    real(double) 	:: tbParams(20)     ! Array for optional parameters that can be
                                             ! read in and used in simple TB models

    !**The following are used only by the general routine that works for any angular momentum
       
    integer                :: numSKparams       ! Number of SK params to be
                                                ! read in
 
    complex(double), pointer :: My(:,:)         ! Used for doing unitary transformations
                                                       
    real(double), pointer :: F(:,:,:,:,:,:,:)   ! This is the matrix of
                                                ! hopping integrals which
                                                ! contains the radial 
                                                ! dependence. (*)
                                                ! F(lI,lJ,m,nI,nJ,speciesI,speciesJ)

    real(double), pointer :: K(:,:,:,:,:,:,:)   ! This is the matrix of 
                                                ! hopping integrals which
                                                ! contains the constant
                                                ! parameters. (*)
                                                ! K(lI,lJ,m,nI,nJ,speciesI,speciesJ)

    integer, pointer :: Lmin(:,:)            ! The minimum value of \ell per
                                             ! species per electron channel. (*)
                                             ! Lmin(n,species)

    integer, pointer :: Lmax(:,:)            ! The maximum value of \ell per
                                             ! species per electron channel. (*)
                                             ! Lmax(n,species)

    type( neighborsT ), pointer :: neighbors
  end type hamInfoT

!******************************************************************************
  type subHamiltonianArrayT
!****************************************************************************
    ! Array for hamiltonian matrix elements in real space 
    ! $H_{m,m^{\prime}}(T)$ with $T \neq 0$ and/or $m \neqm^{\prime}$
    ! The diagonal terms for each orbit are not included

    complex(double), pointer  :: subH(:,:,:,:,:,:)
    integer                   :: dimOrbit, dimChannel, dimNeigh, dimAtom 
  end type subHamiltonianArrayT
 
       real(double) etemp
       integer ltemp, mtemp, numorbtemp
 
  contains

!*************************************************************************
  Subroutine HamInfoInit( hamInfo, struct, tagHandler )
!*************************************************************************
    
    use SysParams,     only : double, hartree2eV
    use TagHandlerMod, only : FindTag, tagHandlerT
    use StructureMod,  only : StructureT
    use NeighborsMod,  only : NeighborsInit

    implicit none

    type( hamInfoT ), pointer    :: hamInfo
    type( StructureT ),  pointer :: struct
    type( tagHandlerT ), pointer :: tagHandler

    ! Workspace variables
    real(double)   :: tempF, tempK
    integer        :: i, j
    integer        :: iorb, index
    integer        :: startn, startorbit
    integer        :: lastn, lastorbit
    integer        :: atomI, atomJ
    integer        :: orbitI, orbitJ
    integer        :: n, ni, nj
    integer        :: l, li, lj
    integer        :: m, mi, mj
    integer        :: mstar, mstarI, mstarJ
    integer        :: species, speciesI, speciesJ
    integer        :: error

    ! abcOrb are indices used when referencing certain sub-blocks of 
    ! They are only present for convenience. The nomenclature may be misleading.
    integer        :: minLmin    ! Minimum value of Lmin.
    integer        :: maxLmax    ! Maximum value of Lmax.
    integer        :: minOrb     ! In general, this number equals (Lmax+1)**2 - 1.
                                 ! Used for indexing purposes only.
    integer        :: maxOrb     ! In general, this number equals Lmin**2.
                                 ! Used for indexing purposes only.
    integer        :: numOrb     ! In general, this number equals maxOrb - minOrb + 1.
                                 ! Used to referece the number of orbitals.
                                 ! Instead of the indices.

    integer        :: maxChannel ! Maximum number of electron channels.
    integer        :: maxNumOrb  ! Maximum value for numOrb
    character      :: c

    !*********************************************************************
    ! Make sure that hamInfo hasn't been already created.
    if( associated( hamInfo )) stop 'HamInfoInit: Error, hamInfo already associated.'
    allocate( hamInfo, stat = error )
    if( error /= 0 ) stop 'HamInfoInit: Error allocating hamInfo.'

   !*********************************************************************
    !Find neighbors of each atom within distance specified by MaximumDistance
    nullify(hamInfo%neighbors) ! This line is important. DO NOT ERASE !
    call NeighborsInit(hamInfo%neighbors, struct, tagHandler)

   !*********************************************************************
    ! "channels" denotes the set of principle quantum numbers defined for each species 
     ! Each channel can specify a set of anugular momenta
     ! Useful to add different orbitals in the basis with the same l,m
     ! numchannels(species) = number for that species
     ! The default is numchannels(species) = 1
    ! Allocate and initialize the channel array.
    allocate(hamInfo%numChannels(struct%atomBasis%numSpecies), stat = error)
    if(error /= 0) stop 'HamInfoInit: Error allocating hamInfo%numChannels.'
    
    hamInfo%numChannels = one	!default
    maxChannel = one		!default
    
    !*********************************************************************  
    ! Read in the  type of Tight Binding Model
    call FindTag( tagHandler,'TightBindingModelType', error)
    if(error /= 0) stop 'Error, couldn''t find TightBindingModelType.'

    read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) hamInfo%tbModelType 
    print *, 'TightBindingModelType = ', hamInfo%tbModelType

    	!*********************************************************************  
    	if (hamInfo%tbModelType .eq. 0) then
    	    	print *, 'TightBindingModelType  0'
      		print *, 'General Tight-binding Model for any Ang. Mom.'
      		print *, 'Written by Nick Romero'
      		print *, 'Must read in', &
            	& 'OrbitsAndEnergies, channels, Lmin, Lmax, and SK matrices'
    	end if
    	
    	!*********************************************************************  
    	if (hamInfo%tbModelType .eq. 1) then
    	      	print *, 'TightBindingModelType  1'
      		print *, 'Harrison Tight-binding Model with universal 2-center terms'
      		print *, 'Limited to only one channel per species - defined by default,&
      		& not read in'
      		print *, 'Must read in', &
            	& 'OrbitsAndEnergies but no other input parameters'
    	end if
         	
    	!*********************************************************************  
    	if (hamInfo%tbModelType .ge. 10) then
    	    	print *, 'TightBindingModelType  Number >= 10'
      		print *, 'Special form that are coded in subroutine TBParamsSpecial'
      		print *, 'Limited to only one channel per species - defined by default,&
      		& not read in'
      		print *, 'Must read in OrbitsAndEnergies and parameters tbParams'
            	print *, 'read in using the keyword TightBindingParameters'		
		    ! First, find the tag.
		    call FindTag( tagHandler, 'TightBindingParameters', error)
		    if(error /= 0) stop 'HamInfoInit: Error, couldn''t find TightBindingParameters'
		    read(unit = tagHandler%fileno, fmt = '(A1)') c ! This should get us to the next line		
		    ! Then, read in the value(s) for the tag.
		      read(unit = tagHandler%fileno, fmt = *, iostat = error) hamInfo%numTBParams		      
		      if (error /= 0) stop 'HamInfoInit: Error reading numTBParams'
		    do i = 1, hamInfo%numTBParams
		      read(unit = tagHandler%fileno, fmt = *, iostat = error) &
		         & hamInfo%tbParams(i)
		      if (error /= 0) stop 'HamInfoInit: Error reading hamInfo%tbParams'
    		    end do
    		    
    		    print *, 'TightBindingParameters'
		    do i = 1, hamInfo%numTBParams
		      print *, i, hamInfo%tbParams(i)
    		    end do
    	  end if
      
    !*********************************************************************  
      	! The following is used only for General Tight-binding Model for any Ang. Mom
      	if (hamInfo%tbModelType .eq. 0) then
   
    ! Read in the number of channels for each species.

    ! First, find the tag.
    call FindTag( tagHandler, 'NumChannelsPerSpecie', error)
    if(error /= 0) stop 'HamInfoInit: Error, couldn''t find channel'

    ! Then, read in the value(s) for the tag.
    read(unit = tagHandler%fileno, fmt = '(A1)') c ! This should get us to the next line
    do i = 1, struct%atomBasis%numSpecies
      read(unit = tagHandler%fileno, fmt = *, iostat = error) &
         & species, hamInfo%numChannels(species)
      if (error /= 0) stop 'HamInfoInit: Error reading hamInfo%numChannels'
    end do

    ! Make sure that the electron channel is greater than or equal to zero
    do species = 1, struct%atomBasis%numSpecies
       if (hamInfo%numChannels(species) < 1) stop 'HamInfoInit: Error, &
         & number of channels for species must be at least zero.'
    end do

    ! Store the maximum number of electron channels in a variable 
    ! used for setting dimensions 
    maxChannel = maxval(hamInfo%numChannels)

      	endif		!if (hamInfo%tbModelType .eq. 0) then
      	

     !********************************************************************* 
     ! RMM 5-19-04
     ! Changed file to calculate Lmin and Lmax instead of reading them in
     ! This is used to fix minOrb and maxOrb for setting dimensions.
     ! XXXX  Must check that this works in all cases
    
       ! Allocate and initialize the Lmin array.
    allocate(hamInfo%Lmin(maxChannel, struct%atomBasis%numSpecies), stat = error)
    if(error /= 0) stop 'HamInfoInit: Error allocating hamInfo%Lmin.'
    hamInfo%Lmin = zero

    ! Allocate and initialize the Lmax array.
    allocate(hamInfo%Lmax(maxChannel, struct%atomBasis%numSpecies), stat = error)
    if(error /= 0) stop 'HamInfoInit: Error allocating hamInfo%Lmax.'
    
     hamInfo%Lmax = zero      
     minLmin = 10 ! Large number - assume actual min always < 10
     maxLmax = one
   
        ! First read of OrbitsAndEnergies to find Lmin and Lmax
        call FindTag( tagHandler,'OrbitsAndEnergies', error)
        if(error /= 0) stop 'HamInfoInit: Error, couldn''t find OrbitsAndEnergies.'
    
        read(unit = tagHandler%fileno, fmt = '(A1)') c
        ! Then, read in the value(s) for the tag.
        ! NOTE: That the OrbitsAndEnergies need NOT be in any particular order after the tag.
        !       However, the code will expect values for all the orbits in every 
        !       channel and species.           
        do i = 1, struct%atomBasis%numSpecies
           do j = 1, hamInfo%numChannels(i)
              read(unit = tagHandler%fileno, fmt = *, iostat = error) species, &
                & n, numorbtemp
              if(error /= 0) stop 'HamInfoInit: Error reading no. of orbits in OrbitsAndEnergies'
     
    	  hamInfo%Lmin(n,species) = 10  ! Large number - assume actual min always < 10
    	  hamInfo%Lmax(n,species) = zero
              ! Read for each orbit and channel (l, mstar, energy)
              do iorb = 1, numorbtemp
                 read(unit = tagHandler%fileno, fmt = *, iostat = error) &
                   & ltemp, mtemp, etemp
                 if(error /= 0) stop 'Error reading types and energies in OrbitsAndEnergies'
                 if(ltemp .lt. hamInfo%Lmin(n,species)) hamInfo%Lmin(n,species) = ltemp
                 if(ltemp .gt. hamInfo%Lmax(n,species)) hamInfo%Lmax(n,species) = ltemp
                 if(ltemp .lt. minLmin) minLmin = ltemp
                 if(ltemp .gt. maxLmax) maxLmax = ltemp
              end do ! iorb         
           end do ! j
        end do ! i
    !RMM
    !********************************************************************* 
  
    !*********************************************************************   
    ! Checks to make sure that Lmin and Lmax are positive and verify that Lmin < Lmax.
    do species = 1, struct%atomBasis%numSpecies
       do n = 1, hamInfo%numChannels(species)
          if (hamInfo%Lmin(n,species) < 0) &
            & stop 'HamInfoInit: Error, Lmin must be greater than zero.'
          if (hamInfo%Lmax(n,species) < 0) &
            & stop 'HamInfoInit: Error, Lmax must be greater than zero.'
          if (hamInfo%Lmax(n,species) < hamInfo%Lmin(n,species)) &
            & stop 'HamInfoInit: Error, Lmax must be greater or equal to Lmin.'
       end do                ! n
    end do                   ! species

    ! Print the values of Lmin and Lmax for each species
    do species = 1, struct%atomBasis%numspecies
       do n = 1, hamInfo%numChannels(species)
          print '(I11,I12,2I6)',species, n, hamInfo%Lmin(n,species), hamInfo%Lmax(n,species)
       end do ! n
    end do    ! species
    
    ! These variables store indices and sizes which are used in determining
    ! the size of the arrays in hamInfo type. The size is determined by the 
    ! largest number of orbitals for any of the species. I'm actually wasting
    ! memory space for convenience
    minLmin   = minval(hamInfo%Lmin) 
    maxLmax   = maxval(hamInfo%Lmax)
    maxNumOrb = (maxLmax + 1)**2 - minLmin**2 
    
    !*********************************************************************
    ! Read in the Orbits and Energies for each element.
    ! NOTE: That the OrbitsAndEnergies need NOT be in any particular order after the tag.
    !       However, the code will expect values for all the orbits in every 
    !       channel and species.
 
    ! Allocate and initialize energy array.
    allocate(hamInfo%energies(maxNumOrb,maxChannel,struct%atomBasis%numSpecies), stat = error)
    if (error /= 0) stop 'HamInfoInit: Error allocating energies matrix.'
    hamInfo%energies = zero

    ! Allocate and initialize orbits array.
    allocate(hamInfo%orbits(maxChannel, struct%atomBasis%numSpecies), stat = error)
    if (error /= 0) stop 'HamInfoInit: Error allocate orbits array.'
    hamInfo%orbits = zero

    ! Allocate and initialize orbitType array.
    allocate(hamInfo%orbitType(2,maxNumOrb,maxChannel,struct%atomBasis%numSpecies), stat = error)
    if (error /= 0) stop 'HamInfoInit: Error allocating orbitType.'
    hamInfo%orbitType = zero

    ! First, find the tag.
    call FindTag( tagHandler,'OrbitsAndEnergies', error)
    if(error /= 0) stop 'HamInfoInit: Error, couldn''t find OrbitsAndEnergies.'

    print *, 'OrbitsAndEnergies'
    read(unit = tagHandler%fileno, fmt = '(A1)') c
    ! Then, read in the value(s) for the tag.

    do i = 1, struct%atomBasis%numSpecies
       do j = 1, hamInfo%numChannels(i)
          read(unit = tagHandler%fileno, fmt = *, iostat = error) species, &
            & n, hamInfo%orbits(n,species)
          if(error /= 0) stop 'HamInfoInit: Error reading no. of orbits in OrbitsAndEnergies'
          print *, 'Element = ', species, 'Channel = ', n

          minOrb = hamInfo%Lmin(n,species)**2
          maxOrb = (hamInfo%Lmax(n,species) + 1)**2 - 1
          numOrb = maxOrb - minOrb + 1

          ! Check for \ell s being consistent with Lmin and Lmax 
          if (n > (hamInfo%numChannels(species)) .OR. (n <= 0)) stop 'HamInfoInit: Electron channel &
            & number is incosistent for this element.'
          if (hamInfo%orbits(n,species) /= numOrb) stop 'HamInfoInit: Number of &
            & orbits for this element and channel are inconsistent with the &
            & values Lmin and Lmax.'
            
          !****************************************
          ! For each orbit and channel read in (l, mstar, energy)
          ! orbittype(1,iorb,n,species) = l     = 0, 1, ...
          ! orbittype(2,iorb,n,species) = mstar = 0 for s; 1 for pz, 2 for px, 3 for py, etc.

          do iorb = 1, hamInfo%orbits(n,species)
             read(unit = tagHandler%fileno, fmt = *, iostat = error) &
               & hamInfo%orbitType(1,iorb,n,species), hamInfo%orbitType(2,iorb,n,species), &
               & hamInfo%energies(iorb,n,species)
             if(error /= 0) stop 'Error reading types and energies in OrbitsAndEnergies'
          end do ! iorb
            
          ! Check for each orbit that the values of (l,mstar) don't violate the rules
          do iorb = 1, hamInfo%orbits(n,species)
             l     = hamInfo%orbitType(1,iorb,n,species)
             mstar = hamInfo%orbitType(2,iorb,n,species)
             if (l < 0) stop 'HamInfoInit: l must be greater than or equal to zero.'
             if (mstar < l**2 .OR. (mstar > (l+1)**2 - 1)) &
               & stop 'HamInfoInit: mstar must be between l**2 and (l+1)**2 - 1'         
          end do ! iorb

         ! For each species, channel, and orbit, print the orbitType and energies.
         print *, ' orbit ', ' ang. mom. ', ' mstar ', '  energy  '
         do iorb = 1, hamInfo%orbits(n,species)
            print '(I7,I11,I7,F10.4)', iorb, hamInfo%orbitType(1,iorb, n, species), &
              & hamInfo%orbitType(2,iorb, n, species), &
              & hamInfo%energies(iorb, n, species)
         end do ! iorb
       end do ! j
    end do ! i
    
    !*************************************************************************  
    ! Energies read in eV must be converted to A.U.
    call FindTag( tagHandler,'InputEnergiesInEV', error)
    if(error == 0) then
       print *, 'InputEnergiesInEV - divide hamInfo%energies by hartree2eV'
      hamInfo%energies = hamInfo%energies / hartree2eV 
    end if
    
    !*********************************************************************
         ! The following is used only for General Tight-binding Model for any Ang. Mom
         if (hamInfo%tbModelType .eq. 0) then
 
    ! Read in the SlaterKoster Two-Center Hopping Integrals
    !
    ! There are two Slater-Koster Matrices:
    ! F(lI, lJ, m, nI, nJ, speciesI, speciesJ) = n
    ! where n is the power of the distance
    ! 
    ! K(lI, lJ, m, nI, nJ, speciesI, speciesJ) = lilj(sigma,pi,delta, etc)
    ! which is distance independent.
    !
    ! Also refer to p. 1451 eq. 8 in
    ! Tight-binding Modelling of Materials
    ! Rep. Prog. Phys. 60 (1997) 1447-1512
    !

    ! Allocate F array
    allocate(hamInfo%F(minLmin:maxLmax, minLmin:maxLmax, 0:maxLmax, &
                     & maxChannel, maxChannel, struct%atomBasis%numSpecies, &
                     & struct%atomBasis%numSpecies), stat = error)
    if (error /= 0) stop 'HamInfoInit: Error allocating F matrix.'

    ! Allocate K array
    allocate(hamInfo%K(minLmin:maxLmax, minLmin:maxLmax, 0:maxLmax, &
                     & maxChannel, maxChannel, struct%atomBasis%numSpecies, &
                     & struct%atomBasis%numSpecies), stat = error)
    if (error /= 0) stop 'HamInfoInit: Error allocating K matrix.'      
 
    ! Initialize SlaterKoster Matrices
    hamInfo%F = zero
    hamInfo%K = zero

    ! First, find the tag.
    call FindTag(tagHandler, 'numSKparams', error)
    if (error /= 0) stop 'HamInfoInit: Error, couldn''t find numSKparams.'

    ! Then, read in the value(s) for the tag.
    read(unit = tagHandler%fileno, fmt = '(I10)', iostat = error) hamInfo%numSKparams 

    ! First, find the tag. 
    call FindTag(tagHandler, 'SKMatrices', error)
    if (error /= 0) stop 'HamInfoInit: Error, couldn''t find SK matrices.'

    print *, 'SKMatrices'
    read(unit = tagHandler%fileno, fmt = '(A1)') c
    print *, ' li ',' lj ',' m ',' ni ',' nj ',' speciesI ',' speciesJ ','   F   ', '   K   '
       
    ! Then, read in the value(s) for the tag.
    do i = 1, hamInfo%numSKparams
       ! Read in the value to a temporary variable.
       read(unit = tagHandler%fileno, fmt = *, iostat = error) &
            & li, lj, m, ni, nj, speciesI, speciesJ, tempF, tempK

       ! By symmetry only abs(m) is relevant
       m = abs(m)

       ! Make sure the indices don't violate any of the rules.
       ! Species 
       if ((speciesI > struct%atomBasis%numSpecies) .OR. (speciesI < 1)) &
         & stop 'HamInfoInit: Error, speciesI does not exist.'
       if ((speciesJ > struct%atomBasis%numSpecies) .OR. (speciesJ < 1)) &
         & stop 'HamInfoInit: Error, speciesJ does not exist.'  

       ! Electron channel
       if ((ni > hamInfo%numChannels(speciesI)) .OR. (ni <= 0)) &
         & stop 'HamInfoInit: Error, ni does not exist for speciesI.'
       if ((nj > hamInfo%numChannels(speciesJ)) .OR. (nj <=0 )) &
         & stop 'HamInfoInit: Error, nj does not exist for speciesJ.'            

       ! Angular Momentum 
       if ((li < hamInfo%Lmin(ni,speciesI)) .OR. (li > hamInfo%Lmax(ni,speciesI))) &
         & stop 'HamInfoInit: Error, li angular momentum out of bounds.'
       if ((lj < hamInfo%Lmin(nj,speciesJ)) .OR. (lj > hamInfo%Lmax(nj,speciesJ))) &
         & stop 'HamInfoInit: Error, lj angular momentum out of bounds.'
       ! Angular momentum z-component
       if (abs(m) > max(li,lj)) stop 'HamInfoInit: Error, m is greater than max(li,lj).'
       if ((abs(m) > min(li,lj)) .AND. (tempK /= zero)) &
         & stop 'HamInfoInit: Error, K must be zero when m is greater than min(li,lj).'

       ! Now, read in the value for the tag since everything checks out.
       hamInfo%F(li, lj, m, ni, nj, speciesI, speciesJ) = tempF
       hamInfo%K(li, lj, m, ni, nj, speciesI, speciesJ) = tempK

       ! Impose symmetrization on SK Matrices
       hamInfo%F(lj, li, m, nj, ni, speciesJ, speciesI) = tempF
       hamInfo%K(lj, li, m, nj, ni, speciesJ, speciesI) = (-1)**(li+lj)*tempK

    end do ! i

    ! For each element and channel pair, print out the SK matrices.
    do speciesI = 1, struct%atomBasis%numSpecies
       do speciesJ = speciesI, struct%atomBasis%numSpecies
          do ni = 1, hamInfo%numChannels(speciesI)
             do nj = 1, hamInfo%numChannels(speciesJ)
                do li = hamInfo%Lmin(ni,speciesI), hamInfo%Lmax(ni,speciesI)
                   do lj = hamInfo%Lmin(nj,speciesJ), hamInfo%Lmax(nj,speciesJ)
                      do m = 0, max(li,lj)
 10                      format (2I4,I3,2I4,2I10,2F7.2)
                         print 10, li,lj,m,ni,nj,speciesI,speciesJ, &
                         & hamInfo%F(li,lj,m,ni,nj,speciesI,speciesJ), &
                         & hamInfo%K(li,lj,m,ni,nj,speciesI,speciesJ)
                      end do ! m
                   end do ! lj
                end do ! li
             end do ! nj
          end do ! ni
       end do ! speciesJ
    end do ! speciesI
   
   endif		!if (hamInfo%tbModelType .eq. 0) then
   
     !*************************************************************************  
     ! Energies read in eV must be converted to A.U.
     call FindTag( tagHandler,'InputEnergiesInEV', error)
     if(error == 0) then
       print *, 'InputEnergiesInEV - divide slaterKoster parameters hamInfo%K by hartree2eV'
       hamInfo%K        = hamInfo%K / hartree2eV
     end if
    
    !*************************************************************************
    ! Allocate and initialize the hamiltonian and the hamiltonian index arrays
    
    allocate(hamInfo%hIndex(maxNumOrb,maxNumOrb,maxChannel,maxChannel, &
                          & struct%atomBasis%totalNumAtoms, &
                          & struct%atomBasis%totalNumAtoms), stat = error)
    if (error /= 0) stop 'HamInfoInit: Error allocating hIndex.'

    ! Initialize the hamiltonian and the hamiltonian index.
    index = 0
    hamInfo%hIndex = 0
    hamInfo%numTotalOrbits = 0

    ! Count all the orbitals by looping over all the atoms and channels.
    do atomI = 1, struct%atomBasis%totalNumAtoms

       ! species of atomI
       speciesI = struct%atomBasis%speciesOfAtom(atomI)

       ! Loop over all the electron channels
       do ni = 1, hamInfo%numChannels(speciesI)

          do orbitI = 1, hamInfo%orbits(ni,speciesI)      

             ! Increment the size of the Hamiltonian
             hamInfo%numTotalOrbits = hamInfo%numTotalOrbits + 1

          end do                ! orbitI
       end do                   ! ni
    end do                      ! atomI
 

    ! Since the Hamiltonian is Hermitian we only need the upper right triangle
    ! of the Hamiltonian, so we create an index into the long vector that 
    ! represents the Hamiltonian.
    do atomJ = 1, struct%atomBasis%totalNumAtoms
       
       ! species of atomJ
       speciesJ = struct%atomBasis%speciesOfAtom(atomJ)

       ! Loop over all electron channels on atomJ
       do nj = 1, hamInfo%numChannels(speciesJ)

          do orbitJ = 1, hamInfo%orbits(nj,speciesJ)
          
             do atomI = 1, atomJ

                ! species of atomJ
                speciesI = struct%atomBasis%speciesOfAtom(atomI)

                if (atomI == atomJ) then
                   lastn = nj
                else
                   lastn = hamInfo%numChannels(speciesI)
                end if

                ! Loop over all electron channels on atomI
                do ni = 1, lastn

                   if ((atomI == atomJ) .AND. (ni == nj)) then
                      lastorbit = orbitJ
                   else
                      lastorbit = hamInfo%orbits(ni,speciesI)
                   end if

                   do orbitI = 1, lastorbit

                      ! Store the index.
                      index = index + 1
                      hamInfo%hIndex(orbitI,orbitJ,ni,nj,atomI,atomJ) = index

                   end do       ! orbitI
                end do          ! ni
             end do             ! atomI
          end do                ! orbitJ
       end do                   ! nj
      end do                    ! atomJ

    print *, 'Total No. of orbitals = size of Hamiltonian = ', &
          & hamInfo%numTotalOrbits

   print *, 'Index = ',index       
  end Subroutine HamInfoInit

!*************************************************************************
  Subroutine HamInfoDestroy( hamInfo )
!*************************************************************************
! HamInfoDestroy deallocates orbits, orbittype, energies, 
!  and hIndex
    
    use NeighborsMod, only : NeighborsDestroy
    implicit none

    type( hamInfoT ),    pointer :: hamInfo
    integer :: error

    if( .not. associated( hamInfo )) stop 'Error: hamInfo not allocated.'

    call NeighborsDestroy( hamInfo%neighbors )
    deallocate( hamInfo%orbits, hamInfo%orbitType, hamInfo%energies, &
              & hamInfo%hIndex, stat = error)
    
    if(error /= 0) stop 'Error deallocating orbits, energies and hIndex.'

    deallocate( hamInfo, stat = error )
    if( error /= 0 ) stop 'Error: hamInfo deallocation failed.'
    nullify( hamInfo )

  end Subroutine HamInfoDestroy


!*****************************************************************************
  Subroutine SubHInit(subHamiltonian, maxOrbit, maxChannel, maxNeigh, maxAtom)
!*****************************************************************************
!   Allocate space for the sub-Hamiltonian

    implicit none

    type( subHamiltonianArrayT), pointer :: subHamiltonian
      
    integer, intent(in)                  :: maxOrbit, maxNeigh
    integer, intent(in)                  :: maxChannel, maxAtom

    ! Workspace Variable(s)
    integer  :: error
 
    if( associated( subHamiltonian )) stop 'Error, subHamiltonian already &
    &  allocated.'
    allocate ( subHamiltonian, stat = error )
    if( error /= 0) stop 'HSubInit: Error in allocating subHamiltonian.'
    
    subHamiltonian%dimOrbit   = maxOrbit
    subHamiltonian%dimChannel = maxChannel
    subHamiltonian%dimNeigh   = maxNeigh
    subHamiltonian%dimAtom    = maxAtom

    ! Allocate subHamiltonian array
    allocate(subHamiltonian%subH(maxOrbit,maxOrbit,maxChannel,maxChannel, &
    & maxNeigh,maxAtom), stat = error)
    if (error /= 0) stop 'HSubInit: Error in allocating subH array.'

  end Subroutine SubHInit
                              

!*****************************************************************************
  Subroutine SubHDestroy(subHamiltonian)
!*****************************************************************************

  implicit none

  type( subHamiltonianArrayT), pointer  :: subHamiltonian

  ! Workspace variable(s)
  integer  :: error

   if( .not. associated( subHamiltonian )) &
   & stop 'SubHDestry: Error, subHamiltonian not allocated.'

   deallocate( subHamiltonian%subH, stat=error )
   if(error /= 0) stop 'SubHDestroy: Error in deallocating subH array.'

   deallocate( subHamiltonian, stat = error )
   if(error /= 0) stop 'SubHDestroy: Error in deallocating subHamiltonian.'

    nullify( subHamiltonian )

  end Subroutine SubHDestroy
  
!*************************************************************************
  Subroutine subHGenerate(hamInfo, struct, subH)
!*************************************************************************
!   Generates the tight binding sub-hamiltonian which is k-point independent.

    use SysParams,            only : double, zero
    use StructureMod,         only : StructureT
    use SlaterKosterParamMod, only : TBVInit_Gen, RealOrbitals, LyEigenstates 
    use SlaterKosterParamMod, only : TBHMatElem_Gen
    use SlaterKosterParamMod, only : TBVInit, TBHMatrixElem
    implicit none

    type( hamInfoT ),      pointer :: hamInfo
    type( StructureT ),    pointer :: struct
  
    ! This is the k-pt independent part of the Hamiltonian.
    ! hsub(orbitJ,orbitI,nj,ni,neigh,atomI)
    ! NOTE: Make sure this agrees with the type declared in 
    ! subHamiltonianMod.f90 otherwise there will be a disaster.

    complex(double), intent(out)   :: subH(:,:,:,:,:,:)

    ! Workspace variables
    complex(double), pointer :: My(:,:)
    complex(double), pointer :: Cb(:,:)
    complex(double), pointer :: VSlaterKoster(:,:)

    complex(double) :: phaseFactor  

    integer         :: index, neigh
    integer         :: startorbit, startn
    integer         :: lastorbit, lastn
    integer     :: atomI, atomJ
    integer         :: speciesI, speciesJ
    integer         :: orbitI, orbitJ
    integer         :: li, mstarI, ni
    integer         :: lj, mstarJ, nj
    integer         :: error
  
    real(double)    :: distance
    real(double)    :: xvec(3) ! 3d vector to neighbor needed in 
                               
    ! minOrb and maxOrb are used for calculating dimensions for tbModelType = 5 
    integer        :: minLmin ! Minimum value of Lmin.
    integer        :: maxLmax ! Maximum value of Lmax.
    integer        :: minOrb  ! In general, this number equals (Lmax+1)**2 - 1.
                              ! Used for indexing purposes only.
    integer        :: maxOrb  ! In general, this number equals Lmin**2.
                              ! Used for indexing purposes only.

    ! Debug
    real(double)   :: MatrixElem

    ! Initialize the Hamiltonian
    subH = zero

    !********************************************************************* 
    !RMM 
    ! The following is used only for the Harrison model and other simple models
    ! It does not use the general rotation method or the general F and K matrices
    if (hamInfo%tbModelType .ne. 0) then

    ! Loop over every atom in the unit cell.
    do atomI = 1, struct%atomBasis%totalNumAtoms

       ! species of atomI
       speciesI = struct%atomBasis%speciesOfAtom(atomI)

       ! Loop over all the neighboring atoms
       do neigh = 1, hamInfo%neighbors%neighborList(atomI)%numNeighbors

          ! number of atomJ in the list of all atoms in the unit cell
          atomJ = hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%id

          ! species of atomJ
          speciesJ = struct%atomBasis%speciesOfAtom(atomJ)

          ! Note: vector to neighbor atomJ is hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%x   
          ! The vector from atomI to atomJ is needed in TBHMatrixElem
          ! distance is magnitude, xvec(3) is unit vector
          distance = hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%distance
          xvec = zero
          xvec(1:struct%ndim) = hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%x
       
          ! Loop over every electron channel on each atom
          do ni = 1, hamInfo%numChannels(speciesI)
        
             if (atomI == atomJ) then
                startn = ni
             else
                startn = 1
             end if

             ! Loop over every electron channel on each atom
             do nj = startn, hamInfo%numChannels(speciesJ)
 
                ! Find vSK(:) Slater Koster matrix elements for the specific model
       		 call TBVInit( hamInfo%vSK, hamInfo%tbModelType, hamInfo%tbParams, &
       		 & speciesI, speciesJ, distance )

                ! Loop over every orbital in each electron channel
                do orbitI = 1, hamInfo%orbits(ni,speciesI)

                   ! Identify l,m for orbitI
                   li = hamInfo%orbitType(1,orbitI,ni,speciesI)
                   mstarI = hamInfo%orbitType(2,orbitI,ni,speciesI)
        
                   if ((atomI == atomJ) .AND. (ni == nj)) then
                      startorbit = orbitI
                   else
                      startorbit = 1
                   end if

                   ! Loop over every orbital in each electron channel
                   do orbitJ = startorbit, hamInfo%orbits(nj,speciesJ)

                      ! Identify l,m for orbitj
                      lj = hamInfo%orbitType(1,orbitJ,nj,speciesJ)
                      mstarJ = hamInfo%orbitType(2,orbitJ,nj,speciesJ)
        
                      subH(orbitI,orbitJ,ni,nj,neigh,atomJ) = &
                      & TBHMatrixElem( hamInfo%vSK, xvec, li, mstarI, lj, mstarJ )
                      
                   end do       ! orbitI
                end do          ! orbitJ
             end do             ! ni
          end do                ! nj
       end do                   ! neigh
      end do                    ! atomJ

	endif    	!end of if tbModelType .ne. 0
	

    !********************************************************************* 
    !RMM 
    ! The following is used only for General Tight-binding Model for any Ang. Mom
    if (hamInfo%tbModelType .eq. 0) then

    ! Allocate workspaces for matrix elements of all orbitals
    ! in one channel on one atom with all orbitals in one channel
    ! on another atoms.  Set maximum dimensions for the largest range of angular
    ! momenta ever required for any matrix elements in the speciic problem 
    ! We will repeatedly use these workspaces.

    minLmin = minval(hamInfo%Lmin)
    maxLmax = maxval(hamInfo%Lmax)
    minOrb = minLmin**2
    maxOrb = (maxLmax + 1)**2 - 1
    ! Allocate My matrix
    allocate(My(minOrb:maxOrb, minOrb:maxOrb), stat = error)
    if (error /= 0) stop 'HSubGenerate: Error allocating My matrix.'
    ! Allocate Cb matrix
    allocate(Cb(minOrb:maxOrb, minOrb:maxOrb), stat = error)
    if (error /= 0) stop 'HSubGenerate: Error allocating Cb matrix.'

        
     ! Calculate the transformation matrices for the maximum range of angular momenta

    ! Calculate the entire My matrix.
    call LyEigenstates(minLmin, maxLmax, My)

! XXXX CHANGE NAME OF Cb? 
    ! Calculate all the directed orbital basis vectors.
    call RealOrbitals(minLmin, maxLmax, Cb)

    ! Loop over every atom in the unit cell.
    do atomI = 1, struct%atomBasis%totalNumAtoms

       ! species of atomI
       speciesI = struct%atomBasis%speciesOfAtom(atomI)

       ! Loop over all the neighboring atoms
       do neigh = 1, hamInfo%neighbors%neighborList(atomI)%numNeighbors

          ! number of atomJ in the list of all atoms in the unit cell
          atomJ = hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%id

          ! species of atomJ
          speciesJ = struct%atomBasis%speciesOfAtom(atomJ)

          ! Note: vector to neighbor atomJ is hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%x   
          ! The vector from atomI to atomJ is needed in TBHMatrixElem_Gen
          ! distance is magnitude, xvec(3) is unit vector
          distance = hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%distance
          xvec = zero
          xvec(1:struct%ndim) = hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%x
       
          ! Loop over every electron channel on each atom
          do ni = 1, hamInfo%numChannels(speciesI)
        
             ! Note that we may be doing extra work
             ! For example
             ! I'm doing (ni,nj) = (1,1) , (1,2), (2,1), (2,2) etc.

             if (atomI == atomJ) then
                startn = ni
             else
                startn = 1
             end if

             ! Loop over every electron channel on each atom
             do nj = startn, hamInfo%numChannels(speciesJ)
              
                ! Calculate the dimensions of VSlaterKoster
                ! Here we find the maximum dimension needed and we allocate
                ! a square matrix even though the non-zero part may not be square
                minLmin = min(hamInfo%Lmin(ni,speciesI), hamInfo%Lmin(nj,speciesJ))
                maxLmax = max(hamInfo%Lmax(ni,speciesI), hamInfo%Lmax(nj,speciesJ))

                minOrb = minLmin**2
                maxOrb = (maxLmax + 1)**2 - 1

                ! Allocate VSlaterKosterMatrix
                allocate(VSlaterKoster(minOrb:maxOrb, minOrb:maxOrb), stat = error)
                if (error /= 0) stop 'HGenerate: Error allocating VSlaterKoster.'

                ! Calculate VSlaterKoster for the specific model
                call TBVInit_Gen(VSlaterKoster, &
                & hamInfo%F(minLmin:maxLmax,minLmin:maxLmax,0:maxLmax,ni,nj,speciesI,speciesJ), &
                & hamInfo%K(minLmin:maxLmax,minLmin:maxLmax,0:maxLmax,ni,nj,speciesI,speciesJ), &
                & xvec, distance, minLmin, maxLmax, My(minOrb:maxOrb,minOrb:maxOrb))


                ! Loop over every orbital in each electron channel
                do orbitI = 1, hamInfo%orbits(ni,speciesI)

                   ! Identify l,m for orbitJ
                   li = hamInfo%orbitType(1,orbitI,ni,speciesI)
                   mstari = hamInfo%orbitType(2,orbitI,ni,speciesI)
        
                   if ((atomI == atomJ) .AND. (ni == nj)) then
                      startorbit = orbitI
                   else
                      startorbit = 1
                   end if

                   ! Loop over every orbital in each electron channel
                   do orbitJ = startorbit, hamInfo%orbits(nj,speciesJ)

                      ! Identify l,m for orbitj
                      lj = hamInfo%orbitType(1,orbitJ,nj,speciesJ)
                      mstarJ = hamInfo%orbitType(2,orbitJ,nj,speciesJ)
        
                      ! Uncomment the line below to use Y_lm basis
                      ! subH(orbitI,orbitJ,ni,nj,neigh,atomJ) = VSlaterKoster(mstarI,mstarJ)

                      ! Uncomment the line below to use directed orbital basis
                      ! Convention for ket and bra are different when we have non-symetric
                      ! angular momentum channels, e.g. s & p on 1 and just s on 2
                      subH(orbitI,orbitJ,ni,nj,neigh,atomJ) = & 
                      & TBHMatElem_Gen(VSlaterKoster, Cb(minOrb:maxOrb,mstarJ), Cb(minOrb:maxOrb,mstarI))
                      
                      ! XXXXX CHECK
                      ! line below doesn't work for case described above
                      ! & TBHMatElem_Gen(VSlaterKoster, Cb(minOrb:maxOrb,mstarI), Cb(minOrb:maxOrb,mstarJ))

                   end do       ! orbitI
                end do          ! orbitJ
        
                ! Deallocate VSlaterKoster
                deallocate(VSlaterKoster, stat = error)
                if (error /= 0) stop 'HGenerate: Error deallocating K matrix.'

             end do             ! ni
          end do                ! nj
       end do                   ! neigh
      end do                    ! atomJ

  ! Deallocate My matrix
  deallocate(My, stat = error)
  if (error /= 0) stop 'HSubGenerate: Error deallocating My matrix.'
  ! Deallocate Cb matrix
  deallocate(Cb, stat = error)
  if (error /= 0) stop 'HSubGenerate: Error deallocating Cb matrix,'

	endif    	!end if tbModelType .eq. 0
	
  end Subroutine subHGenerate

!*************************************************************************
  Subroutine HGenerate( hamInfo, struct, k, kPoints, subH, h)
!*************************************************************************
!   Generates the tight-binding Hamiltonian for each k-point.
!   The generic name of the subroutines  "HGenerate" must be used 
!   so it is  called properly by the general hamiltonian routines

    use SysParams,            only : double, imaginary, zero, twopi
    use kPointsMod,           only : kPointsT
    use StructureMod,         only : StructureT
  
    implicit none

    type( hamInfoT ),      pointer :: hamInfo
    type( StructureT ),    pointer :: struct
    type( kPointsT ),      pointer :: kPoints

    complex(double), intent(in)    :: subH(:,:,:,:,:,:)
    complex(double), intent(inout) :: h(:)
    
    integer, intent(in)            :: k

    ! Workspace variables
    complex(double) :: phaseFactor
  
    integer         :: index, neigh
    integer         :: startn , startorbit
    integer         :: lastn, lastorbit
    integer         :: atomI, atomJ
    integer         :: orbitI, orbitJ
    integer         :: ni, nj
    integer         :: speciesI, speciesJ
    integer         :: error

    real(double)    :: distance
    real(double)    :: xvec(3) 

    ! Initialize the Hamiltonian
    h = zero

    ! Diagonal k-independent terms
    ! These are the orbital-energies terms 
    ! See Martin, Electronic Structure, 2004, Eq. 14.4
    ! Contributions to $H_{m,m}(k)$ from diagonal terms in
    ! $H_{m,m^{\prime}}(T)$ with $m = m^{\prime}$ and $T=0$
    
    do atomI = 1, struct%atomBasis%totalNumAtoms

       ! species of atomI 
       speciesI = struct%atomBasis%speciesOfAtom(atomI)

       ! Loop over every electron channel on each atom
       do ni = 1, hamInfo%numChannels(speciesI)

          ! Loop over every orbital in each electron channel
          do orbitI = 1, hamInfo%orbits(ni,speciesI)

             index = hamInfo%hIndex(orbitI, orbitI, ni, ni, atomI, atomI)
             h(index) = h(index) + hamInfo%energies(orbitI,ni,speciesI)

          end do              ! orbitI
       end do                 ! ni
    end do                    ! atomI

    
    ! Two-center k-dependent terms
    ! See Martin, Electronic Structure, 2004, Eq. 14.4
    ! Contributions to $H_{m,m}(k)$ from terms in
    ! $H_{m,m^{\prime}}(T)$ with $T \neq 0$ and/or $m \neqm^{\prime}$
    ! The phaseFactor is exp(ik xvec) where xvec is the relative position of
    ! the two atoms $xvec = T + tau_m - tau_{m^{\prime}}$ 
    ! This is a different choice from Eq. (14.4) but this only affects the
    ! and has no physical consequence
   
    do atomI = 1, struct%atomBasis%totalNumAtoms

       ! species of atomI
       speciesI = struct%atomBasis%speciesOfAtom(atomI)

       ! Loop over all the neighboring atoms
       do neigh = 1, hamInfo%neighbors%neighborList(atomI)%numNeighbors

          ! Number of atomJ
          atomJ = hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%id

          ! species of atomJ
          speciesJ = struct%atomBasis%speciesOfAtom(atomJ)

          ! Note: vector to neighbor atomJ is hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%x   
          ! This is used in phaseFactor in phaseFactor where the dimension is ndim
          distance = hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%distance
          xvec = zero
          xvec(1:struct%ndim) = hamInfo%neighbors%neighborList(atomI)%neighbor(neigh)%x

          ! Find phase factor for given k-pont - this is only place k enters
          phaseFactor = exp(imaginary * twopi / struct%lattice%latticeConst * &
            & dot_product(kPoints%kPointsCart(:,k), &
            & xvec * distance))

          ! Loop over every electron channel on each atom
          do ni = 1, hamInfo%numChannels(speciesI)

             ! XXXXX CHECK 
             ! I think that the problem with this is that we are over counting.
             ! For example
             ! I'm doing (ni,nj) = (1,1) , (1,2), (2,1), (2,2) etc.

             if (atomI == atomJ) then
                startn = ni
             else
                startn = 1
             end if

             ! Loop over every electron channel on each atom
             do nj = startn, hamInfo%numChannels(speciesJ)

                ! Loop over every orbital in each electron channel
                do orbitI = 1, hamInfo%orbits(ni,speciesI)

                   if ((atomI == atomJ) .AND. (ni == nj)) then
                      startorbit = orbitI
                   else
                      startorbit = 1
                   end if
             
                   ! Loop over every orbital in each electron channel
                   ! Accumulate phaseFactor * matrix element of real space TB hamiltonian 
                   do orbitJ = startorbit, hamInfo%orbits(nj,speciesJ)

                      index = hamInfo%hIndex(orbitI,orbitJ,ni,nj,atomI,atomJ)

                      h(index) = h(index) + phaseFactor*subH(orbitI,orbitJ,ni,nj,neigh,atomJ)

                   end do       ! orbitJ
                end do          ! orbitI
             end do             ! nj
          end do                ! ni
       end do                   ! neigh
    end do                      ! atomI

end Subroutine HGenerate

end Module tbHamMod
