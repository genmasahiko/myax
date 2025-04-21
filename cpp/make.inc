CXX	= /opt/homebrew/bin/g++-14

#for cygnus, pegasus
#CXX = icpc

CXXFLAGS = -Iinclude

SRCDIR = src
BUILDDIR = build
OBJDIR = build/obj

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

DEPS = $(OBJECTS:.o=.d)
