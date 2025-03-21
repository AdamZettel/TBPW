# Flags for Linux using IFC and MKL
# If you are getting undefined reference errors in the external
# subroutine calls try some combinations of the flags below.
# -DAUSLAPACK
# -DAUSBLAS
FORTRAN = f90
FFLAGS  = -O3 -fpp -M../$(COMDIR) 
SUF     = .f90
# Uncomment the lines below if you are using lapack
LIBDIR  = -L/lang/SUNWspro/lib
LAPACK  = -lsunperf
BLAS    = -lsunmath


LIBS    = -L../$(COMDIR) -l$(COMPROJ) \
          $(LIBDIR) $(LAPACK) $(BLAS) 

