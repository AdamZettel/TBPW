# This is the master makefile, and make should always be called from this
# directory.  All of the flags and variables for different machines are in this
# file.  Make all changes to these flags and variables here.  If you wish to
# compile for one of the machines listed below uncomment the flags associated
# with that machine and comment out any flags associated with other machines
#
# Make alone or make pw will make the planewave code, make tb will
# make the tight binding code, make clean will remove all 
# compiler generated files, and make save will generate a tar file of all 
# the source and make files.  
#
# If you would like to add a source file to a
# specific project edit the makefile in that directory and add the source
# file name to both the source and object file variables using the same
# format that is used for other files.
# When this makefile is used it changes directories to the target directory
# and starts make in that directory passing the flags and the target.

# Compilation and directory-related variables
MAKE    = make
COMDIR  = Common
PWDIR   = PW
TBDIR   = TB
TAR     = tar cvf
INCLUDE = -I../Common
COMLIB  = libcommon.a

# System specific make file
include arch.make

TIME_STAMP=`date +%b_%d_%H:%M`

# Target-related variables
COMPROJ = common
PWPROJ  = pw
PWDPROJ = pwDensity
TBPROJ  = tb
PROJ    = els

#Plane wave code
$(PWPROJ): $(COMLIB) 
	(cd $(PWDIR); $(MAKE) "FORTRAN=$(FORTRAN)" "FFLAGS=$(FFLAGS)" "INCLUDE=$(INCLUDE)" "SUF=$(SUF)" "COMLIB=$(COMLIB)" "LIBFFTW=$(LIBFFTW)" "LIBS=$(LIBS)" $(PWPROJ))

#Plane wave Density code
$(PWDPROJ): $(COMLIB) 
	(cd $(PWDIR); $(MAKE) "FORTRAN=$(FORTRAN)" "FFLAGS=$(FFLAGS)" "INCLUDE=$(INCLUDE)" "SUF=$(SUF)" "COMLIB=$(COMLIB)" "LIBFFTW=$(LIBFFTW)" "LIBS=$(LIBS)" $(PWDPROJ))

#Tight Binding code
$(TBPROJ): $(COMLIB)
	(cd $(TBDIR); $(MAKE) "FORTRAN=$(FORTRAN)" "FFLAGS=$(FFLAGS)" "INCLUDE=$(INCLUDE)" "SUF=$(SUF)" "COMLIB=$(COMLIB)" "LIBS=$(LIBS)" $(TBPROJ))

#Common library
$(COMLIB): 
	(cd $(COMDIR); $(MAKE) "FORTRAN=$(FORTRAN)" "FFLAGS=$(FFLAGS)" "SUF=$(SUF)" "LIBS=$(LIBS)" $(COMLIB))

save:
	(cd $(COMDIR); $(MAKE) "SUF=$(SUF)" "TAR=$(TAR)" save)
	(cd $(PWDIR); $(MAKE) "SUF=$(SUF)" "TAR=$(TAR)" save)
	(cd $(TBDIR); $(MAKE) "SUF=$(SUF)" "TAR=$(TAR)" save)
	$(TAR) $(PROJ)_$(TIME_STAMP).tar $(COMDIR)/$(COMPROJ).tar \
		$(PWDIR)/$(PWPROJ).tar $(TBDIR)/$(TBPROJ).tar makefile
	rm $(COMDIR)/$(COMPROJ).tar $(PWDIR)/$(PWPROJ).tar \
		$(TBDIR)/$(TBPROJ).tar
	gzip $(PROJ)_*.tar
#	mv $(PROJ)_*.tar.gz $(HOME)/save/.

clean:
	(cd $(COMDIR); $(MAKE) clean)
	(cd $(PWDIR); $(MAKE) clean)
	(cd $(TBDIR); $(MAKE) clean)









