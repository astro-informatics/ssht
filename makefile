# ======== COMPILER ========

CC      = gcc
OPT		= -Wall -O3 -fopenmp -DSSHT_VERSION=\"1.0b1\" -DSSHT_BUILD=\"`git rev-parse HEAD`\"
#OPT	= -Wall -g -fopenmp -DSSHT_VERSION=\"1.0b1\" -DSSHT_BUILD=\"`git rev-parse HEAD`\"


# ======== LINKS ========

UNAME := $(shell uname)
PROGDIR = ..

ifeq ($(UNAME), Linux)
  MLAB		= /usr/local/MATLAB/R2013a
  MLABINC	= ${MLAB}/extern/include
  MLABLIB	= ${MLAB}/extern/lib

  MEXEXT	= mexa64
  MEX 		= ${MLAB}/bin/mex
  MEXFLAGS	= -cxx
endif
ifeq ($(UNAME), Darwin)
  MLAB		= /Applications/MATLAB_R2013a.app
  MLABINC	= ${MLAB}/extern/include
  MLABLIB	= ${MLAB}/extern/lib

  MEXEXT	= mexmaci64
  MEX 		= ${MLAB}/bin/mex
  MEXFLAGS	= -cxx
endif

SSHTDIR  = $(PROGDIR)/ssht
SSHTLIB  = $(SSHTDIR)/lib/c
SSHTLIBNM= ssht
SSHTSRC  = $(SSHTDIR)/src/c
SSHTBIN  = $(SSHTDIR)/bin/c
SSHTOBJ  = $(SSHTSRC)
SSHTINC  = $(SSHTDIR)/include/c
SSHTDOC  = $(SSHTDIR)/doc/c

ifeq ($(UNAME), Linux)
  FFTWDIR      = $(PROGDIR)/fftw-3.2.2_fPIC
endif
ifeq ($(UNAME), Darwin)
  FFTWDIR      = $(PROGDIR)/fftw
endif

FFTWINC	     = $(FFTWDIR)/include
FFTWLIB      = $(FFTWDIR)/lib
FFTWLIBNM    = fftw3
FFTWOMPLIBNM = fftw3_threads

SSHTSRCMAT	= $(SSHTDIR)/src/matlab
SSHTOBJMAT  	= $(SSHTSRCMAT)
SSHTOBJMEX  	= $(SSHTSRCMAT)


# ======== SOURCE LOCATIONS ========

vpath %.c $(SSHTSRC)
vpath %.h $(SSHTSRC)
vpath %_mex.c $(SSHTSRCMAT)


# ======== FFFLAGS ========

FFLAGS  = -I$(FFTWINC) -I$(SSHTINC)
ifeq ($(UNAME), Linux)
  # Add -fPIC flag (required for mex build).
  # (Note that fftw must also be built with -fPIC.)
  FFLAGS += -fPIC
endif

# ======== LDFLAGS ========

LDFLAGS = -L$(SSHTLIB) -l$(SSHTLIBNM) -L$(FFTWLIB) -l$(FFTWOMPLIBNM) -l$(FFTWLIBNM) -lm

LDFLAGSMEX = -L$(SSHTLIB) -l$(SSHTLIBNM) -L$(FFTWLIB) -l$(FFTWOMPLIBNM) -l$(FFTWLIBNM)


# ======== OBJECT FILES TO MAKE ========

SSHTOBJS = $(SSHTOBJ)/ssht_sampling.o    \
           $(SSHTOBJ)/ssht_dl.o          \
           $(SSHTOBJ)/ssht_core.o        \
           $(SSHTOBJ)/ssht_adjoint.o

SSHTHEADERS = ssht_types.h     \
              ssht_error.h     \
	      ssht_sampling.h  \
	      ssht_dl.h        \
	      ssht_core.h      \
	      ssht_adjoint.h

SSHTOBJSMAT = $(SSHTOBJMAT)/ssht_sampling_mex.o        \
              $(SSHTOBJMAT)/ssht_dl_mex.o              \
              $(SSHTOBJMAT)/ssht_dln_mex.o             \
              $(SSHTOBJMAT)/ssht_forward_mex.o         \
              $(SSHTOBJMAT)/ssht_inverse_mex.o         \
              $(SSHTOBJMAT)/ssht_forward_adjoint_mex.o \
              $(SSHTOBJMAT)/ssht_inverse_adjoint_mex.o

SSHTOBJSMEX = $(SSHTOBJMEX)/ssht_sampling_mex.$(MEXEXT)        \
              $(SSHTOBJMEX)/ssht_dl_mex.$(MEXEXT)              \
              $(SSHTOBJMEX)/ssht_dln_mex.$(MEXEXT)             \
              $(SSHTOBJMEX)/ssht_forward_mex.$(MEXEXT)         \
              $(SSHTOBJMEX)/ssht_inverse_mex.$(MEXEXT)         \
              $(SSHTOBJMEX)/ssht_forward_adjoint_mex.$(MEXEXT) \
              $(SSHTOBJMEX)/ssht_inverse_adjoint_mex.$(MEXEXT)


# ======== MAKE RULES ========

$(SSHTOBJ)/%.o: %.c $(SSHTHEADERS)
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

.PHONY: default
default: lib test about

.PHONY: test
test: $(SSHTBIN)/ssht_test about
$(SSHTBIN)/ssht_test: $(SSHTOBJ)/ssht_test.o $(SSHTLIB)/lib$(SSHTLIBNM).a
	$(CC) $(OPT) $< -o $(SSHTBIN)/ssht_test $(LDFLAGS)

.PHONY: about
about: $(SSHTBIN)/ssht_about
$(SSHTBIN)/ssht_about: $(SSHTOBJ)/ssht_about.o
	$(CC) $(OPT) $< -o $(SSHTBIN)/ssht_about

.PHONY: runtest
runtest: test
	$(SSHTBIN)/ssht_test 64 0 32

.PHONY: all
all: lib test about matlab


# Library

.PHONY: lib
lib: $(SSHTLIB)/lib$(SSHTLIBNM).a
$(SSHTLIB)/lib$(SSHTLIBNM).a: $(SSHTOBJS)
	ar -r $(SSHTLIB)/lib$(SSHTLIBNM).a $(SSHTOBJS)


# Matlab

$(SSHTOBJMAT)/%_mex.o: %_mex.c $(SSHTLIB)/lib$(SSHTLIBNM).a
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@ -I${MLABINC}

$(SSHTOBJMEX)/%_mex.$(MEXEXT): $(SSHTOBJMAT)/%_mex.o $(SSHTLIB)/lib$(SSHTLIBNM).a
	$(MEX) $< -o $@ $(LDFLAGSMEX) $(MEXFLAGS) -L$(MLABLIB)

.PHONY: matlab
matlab: $(SSHTOBJSMEX)


# Documentation

.PHONY: doc
doc:
	doxygen $(SSHTSRC)/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -f $(SSHTDOC)/html/*


# Cleaning up

.PHONY: clean
clean:	tidy
	rm -f $(SSHTOBJ)/*.o
	rm -f $(SSHTLIB)/lib$(SSHTLIBNM).a
	rm -f $(SSHTBIN)/ssht_test
	rm -f $(SSHTBIN)/ssht_about
	rm -f $(SSHTOBJMAT)/*.o
	rm -f $(SSHTOBJMEX)/*.$(MEXEXT)

.PHONY: tidy
tidy:
	rm -f *~

.PHONY: cleanall
cleanall: clean cleandoc
