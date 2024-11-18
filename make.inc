F90	= gfortran
CXX	= /opt/homebrew/bin/g++-14
LAPACK_DIR = /opt/homebrew/opt/lapack/lib
LDFLAGS	= -L${LAPACK_DIR} -llapack

#for cygnus, pegasus
#F90 = ifort
#CXX = icpc
#LDFLAGS = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
#CXXFLAGS = -I/home/EIS00001/gen1127/program/myax/include

FFLAGS =
CXXFLAGS = -I/Users/gen/program/myax/include/


