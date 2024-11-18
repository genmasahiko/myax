F90	= gfortran
CXX	= /opt/homebrew/bin/g++-14
LAPACK_DIR = /opt/homebrew/opt/lapack/lib
LDFLAGS	= -L${LAPACK_DIR} -llapack

FFLAGS =
CXXFLAGS = -I/Users/gen/program/myax/include/

TARGET = \
disp

#for cygnus
#CPP = g++
#LDFLAGS = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
