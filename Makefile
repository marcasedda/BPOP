#
# Makefile for BPOP code
# M. Arca Sedda 2025
#


RESULT = BPOP.x
INSTALDIR = /usr/local

VPATH=./build

$(RESULT) : 
	$(MAKE) -C ./build

clean:
	rm -f ./build/*.o ./build/*.so ./build/$(RESULT) 

install: $(RESULT) $(EXTRAOBJS)
	@[ -d $(INSTALLDIR)/bin ] || mkdir $(INSTALLDIR)/bin
	@[ -d $(INSTALLDIR)/share ] || mkdir $(INSTALLDIR)/share
	@[ -d $(INSTALLDIR)/share/doc ] || mkdir $(INSTALLDIR)/share/doc
	cp ./build/$(RESULT)  $(INSTALLDIR)/bin

uninstall:
	rm -f $(INSTALLDIR)/bin/$(RESULT)
