include ../make.inc

.SUFFIXES:
.SUFFIXES: .cpp .f90 .o

BINDIR = ../bin

CXXSRCS = \
	adjustcube.cpp		\
	convcheck.cpp  		\
	h2oangle.cpp   		\
	integrate_norm.cpp	\
	potsurf.cpp    		\
	restart.cpp    		\
	vibration.cpp  		\
	cube.cpp       		\
	data.cpp       		\

FSRCS = \
	anime.f90					\
	disp.f90					\
	vibration.f90				\

TARGETS = \
	adjustcube	\
	convcheck	\
	int			\
	potsurf		\
	restart		\
	h2oangle	\
	anime		\
	disp		\
	vibration	\

OBJS = $(CXXSRCS:.cpp=.o) $(FSRCS:.f90=.o)

all: $(TARGETS)

.f90.o:
	$(F90) $(FFLAGS) -c $< -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

adjustcube: adjustcube.o
	$(CXX) $(CXXFLAGS) -o $@ $^
	ln -sf $(abspath $@) $(BINDIR)/$@

test: test.o
	$(CXX) $(CXXFLAGS) -o $@ $^
	ln -sf $(abspath $@) $(BINDIR)/$@

test: test.o

convcheck: convcheck.o data.o
	$(CXX) $(CXXFLAGS) -o $@ $^
	ln -sf $(abspath $@) $(BINDIR)/$@

disp : disp.o
	$(F90) $(FFLAGS) -o $@ $^
	ln -sf $(abspath $@) $(BINDIR)/$@

vibration : vibration.o
	$(F90) $(FFLAGS) -o $@ $^ $(LDFLAGS)
	ln -sf $(abspath $@) $(BINDIR)/$@

anime : anime.o
	$(F90) $(FFLAGS) -o $@ $^ $(LDFLAGS)
	ln -sf $(abspath $@) $(BINDIR)/$@

potsurf : potsurf.o
	$(CXX) $(CXXFLAGS) -o $@ $^
	ln -sf $(abspath $@) $(BINDIR)/$@

h2oangle : h2oangle.o
	$(CXX) $(CXXFLAGS) -o $@ $^
	ln -sf $(abspath $@) $(BINDIR)/$@

restart : restart.o data.o
	$(CXX) $(CXXFLAGS) -o $@ $^
	ln -sf $(abspath $@) $(BINDIR)/$@

int : integrate_norm.o cube.o
	$(CXX) $(CXXFLAGS) -o $@ $^
	ln -sf $(abspath $@) $(BINDIR)/$@

clean:
	-rm $(BINDIR)/*
	-rm $(TARGETS)
	-rm $(OBJS)
