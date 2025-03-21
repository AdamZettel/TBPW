# This is not the main makefile.  Add only files to be compiled into the
# tight binding code to this file.  See the make file in the parent
# directory for further details.   
TBPROJ = tb
COMDIR = ../Common

# Source Files
TBSRC = slaterKosterParamMod$(SUF) \
	neighborsMod$(SUF) \
	tbHamMod$(SUF) \
	tb$(SUF) 

# Object Files
TBMOD = slaterKosterParamMod.o \
	neighborsMod.o \
	tbHamMod.o \
	tb.o 

$(TBPROJ): $(TBMOD)
	$(FORTRAN) $(FFLAGS) -o $(TBPROJ) $(TBMOD) \
	$(LIBS)

# Dependencies
tb.o: 	tbHamMod.o neighborsMod.o slaterKosterParamMod.o\
	$(COMDIR)/$(COMLIB)

tbHamMod.o: neighborsMod.o  slaterKosterParamMod.o

# Compiler Rules
.SUFFIXES: $(SUF)
$(SUF).o:
	$(FORTRAN) $(FFLAGS) $(INCLUDE) -c $<

save:
	$(TAR) $(TBPROJ).tar $(TBSRC) makefile

clean:
	rm -f  *.o $(TBPROJ) *.mod






