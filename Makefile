#
# Makefile for BPOP code
# M. Arca Sedda 2025
#


RESULT = BPOP.x
# Use the active conda environment unless an explicit PREFIX is supplied.
PREFIX ?= $(CONDA_PREFIX)
DESTDIR ?=

.PHONY: all $(RESULT) clean install uninstall

all: $(RESULT)

# Always consult the build Makefile so header changes are detected there.
$(RESULT):
	$(MAKE) -C build

clean:
	$(MAKE) -C build clean

install:
	$(MAKE) -C build install PREFIX="$(PREFIX)" DESTDIR="$(DESTDIR)"

uninstall:
	$(MAKE) -C build uninstall PREFIX="$(PREFIX)" DESTDIR="$(DESTDIR)"
