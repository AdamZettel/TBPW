# Flags for IBM
FORTRAN = xlf90 -qsuffix=f=f90:cpp=F90
FFLAGS  = -qstrict -O3 
SUF     = .f90
LIBS    = -lblas -L/usr/local/lib -llapack_rs64c  -L../Common -lcommon
