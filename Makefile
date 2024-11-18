include make.inc

SRCDIR = src

.PHONY: all clean $(SRCDIR)

default :
	@echo "Usage: make [TARGET]"
	@echo "TARGETS are as follows:"
	@echo "  all:	 		all programs"
	@echo "  disp:	 		disp"
	@echo "  vibration: 	vibration"
	@echo "  anime: 		anime"
	@echo "  adjust: 		adjustcube"
	@echo "  potsurf: 		potsurf"
	@echo "  h2oangle: 		h2oangle"
	@echo "  convcheck: 	convcheck"
	@echo "  restart: 		restart"
	@echo "  clean: 		clean up"

all: $(SRCDIR)

$(SRCDIR):
	$(MAKE) -C $@

clean:
	$(MAKE) -C $(SRCDIR) clean
	@echo "Cleaned up."
