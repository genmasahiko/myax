
F90	= gfortran
CPP	= /opt/homebrew/bin/g++-13
LAPACK_DIR = /opt/homebrew/opt/lapack/lib
LDFLAGS	= -L${LAPACK_DIR} -llapack

#for cygnus
#CPP = g++
#LDFLAGS = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core

all :
	$(F90) disp.f90 -o disp.x $(LDFLAGS)
	$(F90) vibration.f90 -o vibration.x $(LDFLAGS)
	$(F90) anime.f90 -o anime.x $(LDFLAGS)
	$(CPP) adjustcube.cpp -o adjustcube.x
	$(CPP) potsurf.cpp -o potsurf.x
	$(CPP) h2oangle.cpp -o h2oangle.x
	$(CPP) convcheck.cpp -o convcheck.x

disp:
	$(F90) disp.f90 -o disp.x $(LDFLAGS)

vibration:
	$(F90) vibration.f90 -o vibration.x $(LDFLAGS)

anime:
	$(F90) anime.f90 -o anime.x $(LDFLAGS)

adjust:
	$(CPP) adjustcube.cpp -o adjustcube.x

potsurf:
	$(CPP) potsurf.cpp -o potsurf.x

h2oangle:
	$(CPP) h2oangle.cpp -o h2oangle.x

convcheck:
	$(CPP) convcheck.cpp -o convcheck.x

test:
	$(CPP) -O3 test.cpp -o test.x

clean : 
	rm *.mod
	rm *.o
	rm $(TARGET)
	
