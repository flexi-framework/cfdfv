### CFDFV build configuration file

# select compiler [gnu, intel]
COMPILER = gnu

# select equation system [euler, navierstokes]
EQNSYS = euler

# multithreading flag [on, off]
PARALLEL = off

# debugging flag [on, off]
DEBUG = off

# profiling flag [on, off]
PROF = off

# CGNS library setup
CGNS_VERSION    = 3.4.1
CFDFV_LIB       = libcfdfv.a
INCDIRS         = -I../share/CGNS-$(CGNS_VERSION)/BUILD/include

### Compiler setup

# gnu
ifeq ($(COMPILER), gnu)
  FC      = gfortran-8
  FCFLAGS = -fdefault-real-8 -fdefault-double-8 -fbackslash
  FLFLAGS = -fdefault-real-8 -fdefault-double-8 -fbackslash
  ifeq ($(DEBUG)), on)
    FCFLAGS += -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace -Wunused \
    	       -ffpe-trap=invalid,zero,overflow,underflow,denormal
    FLFLAGS += -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace -Wunused \
    	       -ffpe-trap=invalid,zero,overflow,underflow,denormal
  else
    FCFLAGS += -Ofast -march=native
    FLFLAGS += -Ofast -march=native
  endif
  ifeq ($(PARALLEL)), on)
    FCFLAGS += -fopenmp
    FLFLAGS += -fopenmp
  endif
  ifeq ($(PROF)), on)
    FCFLAGS += -pg
    FLFLAGS += -pg
  endif
endif

# intel
ifeq ($(COMPILER), intel)
  FC      = ifort
  FCFLAGS = -fpp -assume bscc -r8 -i4 -traceback -warn all
  FLFLAGS = -r8 -i4 -traceback -assume bscc
  ifeq ($(DEBUG)), on)
    FCFLAGS += -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,\
	       nooutput_conversion,pointer,bounds,uninit
    FLFLAGS += -g -O0
  else
    FCFLAGS += -O3 -xhost -vec-report0 -no-prec-div
    FLFLAGS += -O3 -xhost -vec-report0 -no-prec-div
  endif
  ifeq ($(PARALLEL)), on)
    FCFLAGS += -openmp
    FLFLAGS += -openmp
  endif
  ifeq ($(PROF)), on)
    FCFLAGS += -p
    FLFLAGS += -p
  endif
endif
