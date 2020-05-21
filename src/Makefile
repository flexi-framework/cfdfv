### CFDFV Makefile

include ../config.mk
include deplist.mk

AR       = ar
AR_FLAGS = clr

.PHONY: clean veryclean cleanshare

SRC = $(shell ls `cat srcfiles.mk|sed -e 's/EQNSYS/'$(EQNSYS)'/g'`)
OBJ = $(SRC:.f90=.o)

all: builddeps
	@$(MAKE) $(CFDFV_LIB)
	@$(FC) $(FCFLAGS) -c cfdfv.f90 -o cfdfv.o
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ' CFDFV COMPILED, NOW LINKING...'
	@echo '---------------------------------------------------------------'
	@$(FC) $(FLFLAGS) cfdfv.o ../lib/$(CFDFV_LIB) -o ../bin/cfdfv \
		-L. -L../lib -L../share/CGNS-$(CGNS_VERSION)/BUILD/lib -lcgns
	@echo '---------------------------------------------------------------'
	@echo ' CFDFV LINKED!'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

builddeps:
	# Create dependencies
	@echo $(SRC) | sed -e 's/ /\n/g' > srcfiles_pp.mk
	@python3 ../utils/builddeps --source=srcfiles_pp.mk deplist.mk

$(CFDFV_LIB): $(OBJ)
	$(AR) $(AR_FLAGS) ../lib/$(CFDFV_LIB) $(OBJ)

$(OBJ): %.o: %.f90
	$(FC) $(FCFLAGS) $(INCDIRS) -c $< -o $@

clean:
	-@rm -f $(OBJ) *.mod cfdfv.o *.i
	-@rm -f cfdfv.o
	-@rm -f ../bin/cfdfv

veryclean: clean
	-@rm -f *~
	-@rm -f ../lib/libcfdfv.a
