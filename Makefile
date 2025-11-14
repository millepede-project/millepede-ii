# #################################################################
# Makefile for MillePede II (Fortran90) with possible input from C
# 
# Author Volker Blobel, University Hamburg, 2005-2009
# Updates by Gero Flucke, CERN, Claus Kleinwort, DESY
#
# Tested on - 64-bit SL5          with gcc version 4.4.4.
#           - 64 bit Ubuntu 11.10 with gcc version 4.6.1.
#
# See comments about different gcc versions inline to get a
# hint about the necessary adjustments.
# #################################################################
#
# ### Define gcc to be used ###
#
GCC=gcc
GCCVERS=$(shell $(GCC) -dumpversion)
# On SL5 replace default (gcc41) by gcc44 
ifeq ($(filter 4.1.,$(GCCVERS)),4.1.)
  GCC=gcc44
  GCCVERS=$(shell $(GCC) -dumpversion)  
endif
$(info Using gcc version $(GCCVERS))
#
# #################################################################
#
# ### Options ###
#
# All but 'yes' disables support of reading C-binaries:
SUPPORT_READ_C = yes
# If yes (and if SUPPORT_READ_C is yes), uses rfio, i.e. shift library, for IO,
# requires shift library and header to be installed at places defined here:
#
SUPPORT_C_RFIO =
# default path (by '=') to RFIO - overwrite if needed...
RFIO_INCLUDES_DIR = =
RFIO_LIBS_DIR = =
# yes
#
# If yes (and if SUPPORT_READ_C is yes and SUPPORT_C_RFIO is not yes),
# use zlib to read gzipped binary files:
SUPPORT_ZLIB = yes
# default path (by '=') to ZLIB - overwrite if needed
# requires z library and header to be installed at places defined here:
ZLIB_INCLUDES_DIR = =
ZLIB_LIBS_DIR = =
# LAPACK (64 bit) with Intel MKL
SUPPORT_LAPACK64 = 
#yes
LAPACK64 = OPENBLAS
LAPACK64_LIBS_DIR = /usr/lib64
LAPACK64_LIB = openblasp64
#MKLROOT= <path to mkl>
#LAPACK64 = MKL
#LAPACK64_LIBS_DIR = $(MKLROOT)/lib/intel64
#LAPACK64_LIB = mkl_rt
# Intel oneMKL PARDISO
SUPPORT_PARDISO = 
#yes
#
# If yes use multithreading with OpenMP (TM)
SUPPORT_OPENMP = yes
# Instrumentation/profiling with Score-P (score-p.org)
#SCOREPROOT = <path to scorep>
# POMP replaced by 'user regions'
SCOREP = 
#$(SCOREPROOT)/bin/scorep --user
#
# make install copies the binary to $(PREFIX)/bin
PREFIX = .
#
# #################################################################
#
FCOMP = $(SCOREP) $(GCC)
F_FLAGS = -Wall -fautomatic -fno-backslash -O3 -cpp
#
CCOMP = $(SCOREP) $(GCC) 
C_FLAGS = -Wall -Wno-unused-function -O3
F_INCLUDEDIRS = # e.g. -I
C_INCLUDEDIRS = # e.g. -I
#
LOADER = $(SCOREP) $(GCC)
L_FLAGS = -Wall -O3
#.
ifeq ($(filter 4.4., $(GCCVERS)), 4.4.)
# gcc44: 
  C_LIBS = -lgfortran -lgfortranbegin
else  
# gcc45, gcc46:
  C_LIBS = -lgfortran -lm
# math library -lm or -lquadmath may be required  
endif
DEBUG =          # e.g. -g
#
ifeq ($(SUPPORT_OPENMP),yes)
# Multithreading with OpenMP (TM)
  L_FLAGS += -fopenmp
  F_FLAGS += -fopenmp
endif
# for Score-P user regions
ifdef SCOREP
  F_INCLUDEDIRS += -I $(SCOREPROOT)/include
endif
#
ifeq ($(SUPPORT_LAPACK64),yes)
# Using LAPACK64, tested: Intel MKL, OPENBLAS
  C_LIBS  += -L$(LAPACK64_LIBS_DIR) -l$(LAPACK64_LIB)
  F_FLAGS += -DLAPACK64=\"$(LAPACK64)\"
  ifeq ($(LAPACK64),MKL)
    $(info Using MKL for LAPACK64)
    $(info - Provide $(LAPACK64_LIB) in runtime environment)
     ifeq ($(SUPPORT_OPENMP),yes)
      $(info - Set number of MKL threads in OMP_NUM_THREADS)
     else 
      $(info - Set number of MKL threads in MKL_NUM_THREADS)
    endif
    $(info - Set MKL_THREADING_LAYER=GNU)
    $(info )
     ifeq ($(SUPPORT_PARDISO),yes)
      $(info - Using Intel oneMKL PARDISO)
      $(info )
      F_FLAGS += -DPARDISO
      F_INCLUDEDIRS += -I $(MKLROOT)/include
    endif
  endif
  ifeq ($(LAPACK64),OPENBLAS)
    $(info Using OPENBLAS for LAPACK64)
    $(info - Provide $(LAPACK64_LIB) in runtime environment)
     ifeq ($(SUPPORT_OPENMP),yes)
      $(info - Set number of OPENBLAS threads in OMP_NUM_THREADS)
     else 
      $(info - Set number of OPENBLAS threads in OPENBLAS_NUM_THREADS)
    endif
    $(info )
  endif  
endif
#
# objects for this project
#
USER_OBJ_PEDE = mpdef.o mpdalc.o mpmod.o mpmon.o mpbits.o mpqldec.o mptest1.o mptest2.o mille.o mpnum.o mptext.o mphistab.o \
	minresDataModule.o minresModule.o minresqlpDataModule.o minresqlpBlasModule.o minresqlpModule.o \
        randoms.o vertpr.o linesrch.o mppar.o Dbandmatrix.o pede.o
#
# Chose flags/object files for C-binary support:
#
ifeq ($(SUPPORT_READ_C),yes)
  F_FLAGS += -DREAD_C_FILES
  USER_OBJ_PEDE += readc.o
  ifeq ($(SUPPORT_C_RFIO),yes)
    C_FLAGS += -DUSE_SHIFT_RFIO -I$(RFIO_INCLUDES_DIR)
    C_LIBS += -L$(RFIO_LIBS_DIR) -lshift
  else
    ifeq ($(SUPPORT_ZLIB),yes)
      C_FLAGS += -DUSE_ZLIB -I$(ZLIB_INCLUDES_DIR)
      C_LIBS += -L$(ZLIB_LIBS_DIR) -lz
    endif
  endif
endif
#
# Make the executables
EXECUTABLES = pede readMilleBinary
#

all:	$(EXECUTABLES)

pede : 	${USER_OBJ_PEDE} Makefile
	$(LOADER) $(L_FLAGS) \
		-o $@ ${USER_OBJ_PEDE} $(C_LIBS) 
#
readMilleBinary : tools/readMilleBinary.py
	cp tools/readMilleBinary.py readMilleBinary
#
clean:
	rm -f *.o *~ */*.o */*~ *.mod */*.mod
#
clobber: clean 
	rm -f $(EXECUTABLES)

install: $(EXECUTABLES) #clean
	mkdir -p $(PREFIX)/bin
	mv $(EXECUTABLES) $(PREFIX)/bin

# Make the object files - depend on source and include file 
#
%.o : %.f90 Makefile
	${FCOMP} ${F_FLAGS} $(F_INCLUDEDIRS) -c $< -o $@
%.o: %.c Makefile
	$(CCOMP) -c $(C_FLAGS) $(DEFINES) $(C_INCLUDEDIRS) $(DEBUG) -o $@ $<
#
# ##################################################################
# Module dependencies
mpbits.o:  mpdef.o mpdalc.o
mpdalc.o:  mpdef.o
mpmod.o:   mpdef.o
mpnum.o:   mpdef.o
mppar.o:   mpdef.o mpdalc.o
mpqldec.o: mpdef.o mpdalc.o
pede.o:    mpdef.o mpmod.o mpdalc.o mptest1.o mptest2.o mptext.o
# ##################################################################
# END
# ##################################################################
