#
# Makefile for BPOP code
# M. Arca Sedda 2025
#


RESULT = BPOP.x
INSTALDIR = /usr/local

VPATH=./build

# Extract cosmological parameters from header file
OMEGA_M := $(shell grep "define OMEGA_M" ./include/input_params.h | awk '{print $$3}')
OMEGA_L := $(shell grep "define OMEGA_L" ./include/input_params.h | awk '{print $$3}')
H0 := $(shell grep "define H0" ./include/input_params.h | awk '{print $$3}')

OUTPUT_DIR = ./include
REDSHIFT_FILE = $(OUTPUT_DIR)/redshift_time.txt

# Generate redshift table before compilation
$(REDSHIFT_FILE): ./include/input_params.h
        python3 ./extra_scripts/generate_redshift_table.py $(OMEGA_M) $(OMEGA_L) $(H0) $@
        chmod 775 $@

$(RESULT) : 
        $(REDSHIFT_FILE)
        $(MAKE) -C ./build

clean:
        rm -f ./build/*.o ./build/*.so ./build/$(RESULT) 
        rm -f *.o *.so $(RESULT) $(REDSHIFT_FILE)
 

install: $(RESULT) $(EXTRAOBJS)
	@[ -d $(INSTALLDIR)/bin ] || mkdir $(INSTALLDIR)/bin
	@[ -d $(INSTALLDIR)/share ] || mkdir $(INSTALLDIR)/share
	@[ -d $(INSTALLDIR)/share/doc ] || mkdir $(INSTALLDIR)/share/doc
	cp ./build/$(RESULT)  $(INSTALLDIR)/bin

uninstall:
	rm -f $(INSTALLDIR)/bin/$(RESULT)
