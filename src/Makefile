
F90	= gfortran
LAPACK_DIR = /usr/local/opt/lapack/lib
LDFLAGS	= -L${LAPACK_DIR} -llapack

#for cygnus
#LDFLAGS = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core

all :
	$(F90) disp.f90 -o disp.x $(LDFLAGS)
	$(F90) vibration.f90 -o vibration.x $(LDFLAGS)
	$(F90) anime.f90 -o anime.x $(LDFLAGS)

disp:
	$(F90) disp.f90 -o disp.x $(LDFLAGS)

vibration:
	$(F90) vibration.f90 -o vibration.x $(LDFLAGS)

anime:
	$(F90) anime.f90 -o anime.x $(LDFLAGS)

test:
	$(F90) test.f90 -o test.x $(LDFLAGS)

clean : 
	rm *.mod
	rm *.o
	rm $(TARGET)
	
