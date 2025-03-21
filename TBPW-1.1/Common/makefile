COMLIB = libcommon.a
COMPROJ= common

AR= ar rv
TAR= tar cvf
FFLAGS = -g

COMSRC= atomBasisMod$(SUF) \
	conjGradMod$(SUF) \
	eigenStatesMod$(SUF) \
	densityArrayMod$(SUF) \
	fft$(SUF) \
	graphMod$(SUF) \
	kPointsMod$(SUF) \
	latticeMod$(SUF) \
	structureMod$(SUF) \
	initializeMod$(SUF) \
	minimizeMod$(SUF) \
	hamiltonianMod$(SUF) \
	tiseOpsMod$(SUF) \
	normalizeMod$(SUF) \
	typeMod$(SUF) \
	sysParams$(SUF) \
	tagHandlerMod$(SUF)

COMMOD=  atomBasisMod.o \
	conjGradMod.o \
	eigenStatesMod.o \
	densityArrayMod.o \
	fft.o \
	graphMod.o \
	kPointsMod.o \
	latticeMod.o \
	structureMod.o \
	initializeMod.o \
	minimizeMod.o \
	hamiltonianMod.o \
	tiseOpsMod.o \
	normalizeMod.o \
	typeMod.o \
	sysParams.o \
	tagHandlerMod.o

$(COMLIB): $(COMMOD)
	$(AR) $(COMLIB) $(COMMOD)

latticeMod.o: sysParams.o tagHandlerMod.o
atomBasisMod.o: sysParams.o latticeMod.o tagHandlerMod.o
structureMod.o: sysParams.o latticeMod.o atomBasisMod.o tagHandlerMod.o
fft.o: sysParams.o
graphMod.o: sysParams.o
kPointsMod.o: sysParams.o latticeMod.o graphMod.o tagHandlerMod.o
eigenStatesMod.o: sysParams.o
typeMod.o : sysParams.o
initializeMod.o : typeMod.o sysParams.o eigenStatesMod.o
normalizeMod.o : sysParams.o typeMod.o eigenStatesMod.o
hamiltonianMod.o : sysParams.o typeMod.o eigenStatesMod.o 
tiseOpsMod.o: sysParams.o 
minimizeMod.o : sysParams.o typeMod.o eigenStatesMod.o normalizeMod.o
conjGradMod.o : sysParams.o typeMod.o eigenStatesMod.o normalizeMod.o \
	initializeMod.o minimizeMod.o tiseOpsMod.o hamiltonianMod.o

.SUFFIXES: .f90
$(SUF).o:
	$(FORTRAN) $(FFLAGS) -c $<

save:
	$(TAR) $(COMPROJ).tar $(COMSRC) makefile

clean:
	rm -f  *.o $(COMLIB) *.mod
