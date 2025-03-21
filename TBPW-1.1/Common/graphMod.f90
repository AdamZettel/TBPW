!*************************************************************************
Module GraphMod
!*************************************************************************
!
! The Graphics module is used  to print a GNU plot with Ticks 
! and Labels and lines from a matrix y, and an axis x.  Multiple y's
! are allowed for each x.  The module can only deal with on graph
! at a time.  Each Graph must be create and destroyed with GraphInit
! and GraphDestroy.  GraphInit must be called before GraphDestroy,
! AddLabel, or CreateGNUPlot.  If creating a second plot, GraphDestroy
! must be called for the first graph before GraphInit can be called for
! the second.
!
! The module contains the public functions:
!	GraphInit(integer :: nl)
!	GraphDestroy
!	AddLabel(character(2) :: label,real(double):: tick)
!	CreateGNUPlot(real(double) :: x(:), real(double) :: y(:,:))
!
! The module contains the public variables:
!	none
!
! The module contains the private functions:
!	none
!
! The module contains the private variables:
!	nlabel	-	The current number of labels on the graph.
!	labelmax-	The maximum number of labels that can be used in the
!				current graph.  It is set by GraphInit(nl) and is equal
!				to nl + 1.
!	xlabel(:)-	The array of labels.  Each label is two characters.
!	xtick(:)-	The x position of each tick(i) for xlabel(i)

  use SysParams,	only : double
  use EigenStatesMod, only : eigenStatesT

  implicit none

  save

  type GraphT
    integer :: nlabel	! current number of labels
    integer :: labelmax	! max number of labels
    integer :: npoints	! number of points on xAxis
    character(2), pointer :: xlabel(:)    ! label marker
    real(double), pointer :: xtick(:)     ! tick-markers
	real(double), pointer :: xAxis(:)     ! points on x axis
  end type GraphT

  contains

!*************************************************************************
  Subroutine GraphInit( graph, npoints, nlines )
!*************************************************************************
!
! Initializes graph data module by creating space 
! for nl+1 labels and ticks.

    type( graphT ), pointer :: graph
    integer, intent(in)	:: npoints, nlines
	integer             :: error

    if( associated( graph )) stop 'Error, graph already allocated.'

    allocate( graph, stat = error )

    if( error /= 0 ) stop 'Error, graph allocation failed.'

    graph%labelmax = nlines + 1
    graph%nlabel = 0
    graph%npoints = npoints   
    allocate(graph%xAxis(npoints), graph%xlabel(graph%labelmax), &
	         & graph%xtick(graph%labelmax), stat=error)
    if(error /= 0) stop 'Error allocating space in GraphInit'

  end Subroutine GraphInit

!*************************************************************************
  Subroutine GraphDestroy( graph )
!*************************************************************************
!
!   Destroys graph data space for nlabel labels and ticks.

    type( graphT ), pointer :: graph
    integer :: error

    if( .not. associated( graph )) stop 'Error, graph not allocated.'

    deallocate( graph%xAxis, graph%xlabel, graph%xtick, stat = error )
    if( error /= 0 ) stop 'Error xlabel deallocating in GraphDestroy.'
    graph%nlabel = 0
    graph%labelmax = 0

    deallocate( graph, stat = error )
    if( error /= 0 ) stop 'Error, graph deallocation failed.'

    nullify( graph )

  end Subroutine GraphDestroy

!*************************************************************************
  Subroutine AddLabel( graph, label, tick )
!*************************************************************************
!
!   Enters graph data labels and ticks.

    type( graphT ), pointer :: graph
    character(2), intent(in)	:: label
    real(double), intent(in)	:: tick

    graph%nlabel = graph%nlabel + 1

    if(graph%nlabel > graph%labelmax) STOP 'Error the number of labels is to larger.'

    graph%xlabel( graph%nlabel ) = label
    graph%xtick( graph%nlabel ) = tick

  end Subroutine Addlabel

!*************************************************************************
  Subroutine PlotBands( graph, eigenStates )
!*************************************************************************
!
!   Plots the bands with labels.

!    use GraphMod,       only : CreateGNUPlot, graphT
    use EigenStatesMod, only : eigenStatesT

    implicit none

    type( graphT ),       pointer :: graph
 !   type( kPointsT ),     pointer :: kPoints
    type( eigenStatesT ), pointer :: eigenStates

! generate a graphics data file
    call CreateGNUPlot( graph, graph%xAxis, eigenStates%eigenValues ) 

  end Subroutine PlotBands

!*************************************************************************
  Subroutine CreateGNUPlot( graph, x, y )
!*************************************************************************
!
! CreateGNUPlot creates the files need to create a two dimensional plot
! with GNU plot.  There are two files create gnuplot.dat, the gnuplot
! command file, and plot.dat, the data file.  To view the graph at the 
! command line gnuplot gnuplot.dat
!
! Inputs:
!	real(double) :: x(i)	The x coordinate of the ith point.
!	real(double) :: y(i,j)	The y coordinate of the ith point for the jth line.
!
! Outputs:
!	file gnuplot.dat
!	file plot.dat
!
    implicit none
 
    type( graphT ), pointer   :: graph
    real(double), intent(in)  :: x(:), y(:,:)
    real(double)              :: ymin,ymax
    integer                   :: i,j,error, ulines, upoints, llines, lpoints

    open(7, file = 'plot.dat')           ! output data for use by gnuplot

    lpoints = lbound(y, 2)
    llines  = lbound(y, 1)
    upoints = ubound(y, 2)
    ulines  = ubound(y, 1)

    if((upoints - lpoints) /= (ubound(x,1) - lbound(x,1))) then
      stop 'Error: unequal number of points in x and y'
    end if

    do i = llines, ulines
      do j = lpoints, upoints

        write(7,700) x(j),y(i,j)

      end do

      write(7,126)   	!Needed to write balank line on IBM
      126 format()

    end do

    close(7)

    ! make a gnuplot data file: at command prompt type the following:
    ! % gnuplot gnuplot.dat
    ! hit any key to end display of plot (program will do this automatically)

    open(7, file = 'gnuplot.dat')
    ymax = maxval(y)
    ymin = minval(y)
    write(7,701) graph%xtick(graph%nlabel),ymin,ymax

    do i = 2, graph%nlabel-1

      write(7,705) graph%xtick(i),ymin,graph%xtick(i),ymax

    end do

    write(7,702,ADVANCE="NO") ( graph%xlabel(i), graph%xtick(i), i=1, graph%nlabel-1)
    write(7,703) graph%xlabel(graph%nlabel), graph%xtick(graph%nlabel)
    write(7,704)

    close(7)

    700 format(8F11.5)
    701 format('set data style lines'/'set nokey'/&
       & 'set xrange [0:',F8.2,']'/'set yrange [',F8.2,' :',F8.2,']')
    702 format('set xtics (',:20('"',A2,'" ',F8.2,','))
    703 format(' ',A2,'" ',F8.2,')')
    704 format('plot "plot.dat"'/'pause -1'/'set output "band.ps"'/&
       & 'set terminal postscript'/&
       & 'replot')
    705 format('set arrow from ',F8.2,',',F8.2,' to ',F8.2,',',F8.2,' nohead')

  end Subroutine CreateGNUPlot

end Module GraphMod
