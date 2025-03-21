# This is not the main makefile.  Add only files to be compiled into the
# plane wave code to this file.  See the make file in the parent directory for
# further details.   
PWPROJ  = pw
PWDPROJ = pwDensity
COMDIR  = ../Common

# Source Files
PWSRC = gVectorsMod$(SUF) \
	pwHamMod$(SUF) \
	atomPotentialMod$(SUF) \
	pwDensity$(SUF) \
	cgParamsMod$(SUF) \
	pw$(SUF)

# Object Files
PWCOREMOD = gVectorsMod.o \
	pwHamMod.o \
	atomPotentialMod.o \
	cgParamsMod.o
PWMOD = $(PWCOREMOD) \
	pw.o 
PWDMOD = $(PWCOREMOD) \
	pwDensity.o 

$(PWPROJ): $(PWMOD)
	$(FORTRAN) $(FFLAGS) -o $(PWPROJ) $(PWMOD) \
	$(LIBFFTW) $(LIBS)

$(PWDPROJ): $(PWDMOD)
	$(FORTRAN) $(FFLAGS) -o $(PWDPROJ) $(PWDMOD) \
	$(LIBFFTW) $(LIBS)

# Dependencies
pw.o: pwHamMod.o cgParamsMod.o $(COMDIR)/$(COMLIB)
pwDensity.o: pwHamMod.o $(COMDIR)/$(COMLIB)
pwHamMod.o: gVectorsMod.o atomPotentialMod.o cgParamsMod.o

# Compiler Rules
.SUFFIXES: $(SUF)
$(SUF).o:
	$(FORTRAN) $(FFLAGS) $(INCLUDE) -c $<

save:
	$(TAR) $(PWPROJ).tar $(PWSRC) makefile

clean:
	rm -f  *.o $(PWPROJ) *.mod
