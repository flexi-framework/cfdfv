#*****************************************************************************!
#
# C F D F V   M A K E F I L E
#
#*****************************************************************************!

include ../Makefile.defs
MAINDIR=../

all: cgns

cgns:
	@# CGNS library
	@ if [ ! -f cgnslib_$(CGNS_VERSION)/build/src/libcgns.a ] ; then \
		if [ ! -d cgnslib_$(CGNS_VERSION) ] ; then \
			tar -xzf cgnslib_$(CGNS_VERSION)-$(CGNS_SUBVERSION).tar.gz 1>/dev/null ; \
		fi ; \
		if [ -d cgnslib_$(CGNS_VERSION)/build ] ; then \
			rm -rf cgnslib_$(CGNS_VERSION)/build ; \
		fi ; \
		mkdir cgnslib_$(CGNS_VERSION)/build && cd cgnslib_$(CGNS_VERSION)/build &&  FC='$(FC)' cmake -DENABLE_FORTRAN=ON -DENABLE_64BIT=ON -DCMAKE_INSTALL_PREFIX=../LINUX64 -DCGNS_BUILD_SHARED=OFF -DENABLE_LEGACY=ON .. && $(MAKE) && $(MAKE) install ; \
	fi ; \

#------------------------------------------------------------------------------
# utility targets
#------------------------------------------------------------------------------
.PHONY: clean veryclean cleanshare

clean:
	@rm -f $(OBJ) *.mod Cfdfv.o *.i
	@rm -f CfdFV.o

veryclean: clean
	@rm -f *~

cleanshare:
	@if [ -d share/cgnslib_$(CGNS_VERSION) ] ; then \
		cd share/cgnslib_$(CGNS_VERSION) && $(MAKE) allclean ; \
	fi
