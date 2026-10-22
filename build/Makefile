#
# Makefile for BPOP
# Manuel Arca Sedda 2025
#

# Respect the compiler selected by an activated conda environment.
CXX ?= g++

RESULT   = BPOP.x
CXXFLAGS = -O3 -std=c++17 -I../include
LDFLAGS = -o $(RESULT)

# Install into the active conda environment; never default to /usr/local.
PREFIX ?= $(CONDA_PREFIX)
DESTDIR ?=

VPATH= ../:../src:../include

INC = input_params.h funcGR.h
SOURCE = Main_Bpop_dev.cpp funcGR.cpp

.SUFFIXES : .o .cpp

OBJECTS = $(SOURCE:.cpp=.o)

# The README's clean/build commands and header changes must both trigger a build.
.PHONY: TARGET clean install uninstall

TARGET: $(RESULT) 

$(RESULT): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS)  $(OBJECTS)

%.o: ../src/%.cpp $(addprefix ../include/,$(INC))
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(RESULT)

install: $(RESULT)
	@test -n "$(PREFIX)" || (echo "Activate a conda environment or set PREFIX explicitly." >&2; exit 1)
	install -d "$(DESTDIR)$(PREFIX)/bin" "$(DESTDIR)$(PREFIX)/share/bpop"
	install -m 755 $(RESULT) "$(DESTDIR)$(PREFIX)/bin/$(RESULT)"
	install -m 644 ../include/input_params.h "$(DESTDIR)$(PREFIX)/share/bpop/input_params.h"

uninstall:
	@test -n "$(PREFIX)" || (echo "Activate a conda environment or set PREFIX explicitly." >&2; exit 1)
	rm -f "$(DESTDIR)$(PREFIX)/bin/$(RESULT)" "$(DESTDIR)$(PREFIX)/share/bpop/input_params.h"

