
F90	= gfortran
LAPACK_DIR = /usr/local/opt/lapack/lib
LDFLAGS	= -L${LAPACK_DIR} -llapack

all :
	$(F90) disp.f90 -o disp.x $(LDFLAGS)
	$(F90) vibration.f90 -o vibration.x $(LDFLAGS)

disp:
	$(F90) disp.f90 -o disp.x $(LDFLAGS)

vibration:
	$(F90) vibration.f90 -o vibration.x $(LDFLAGS)

clean : 
	rm *.mod
	rm *.o
	rm $(TARGET)
	
