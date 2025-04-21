include make.inc

TARGETS = $(BUILDDIR)/main

.PHONY = all clean

all: $(TARGETS)

$(TARGETS): $(OBJECTS) | $(BUILDDIR)
	$(CXX) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@ -MF $(@:.o=.d)

$(BUILDDIR) $(OBJDIR):
	mkdir -p $@

clean:
	rm -rf $(BUILDDIR)

