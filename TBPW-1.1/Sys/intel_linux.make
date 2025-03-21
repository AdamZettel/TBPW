# Flags for Linux using IFC and MKL
# If you are getting undefined reference errors in the external
# subroutine calls try some combinations of the flags below.
# -DAUSLAPACK
# -DAUSBLAS
FORTRAN = ifort
FFLAGS  = -w -O3 -fpp 
SUF     = .f90
G2C     = -lg2c
# Uncomment the lines below if you are using lapack
LAPACK  = -llapack
BLAS    = -lblas

# Uncomment the lines below if you are using MKL
# LIBDIR  = -L/opt/intel/mkl/lib/32
# LAPACK  = -lmkl_lapack
# BLAS    = -lmkl_p3
# GUIDE   = -lguide
# PTHREAD = -lpthread

LIBS    = -L../$(COMDIR) -l$(COMPROJ) \
          $(LIBDIR) $(LAPACK) $(BLAS) $(GUIDE) $(PTHREAD) $(G2C)

