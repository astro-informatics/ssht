# ======== COMPILER ========

CC      = gcc
OPT	= -Wall -g
#OPT	= -Wall -O3

MLAB	= /Applications/MATLAB_R2010a.app
MLABINC	= ${MLAB}/extern/include
MLABLIB	= ${MLAB}/extern/lib

MEXEXT	= mexmaci64
MEX 	= ${MLAB}/bin/mex
MEXFLAGS= -cxx


# ======== LINKS ========

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
  PROGDIR = /home/jdm57/src
endif
ifeq ($(UNAME), Darwin)
  PROGDIR = /Users/mcewen/src
endif

SSHTDIR  = $(PROGDIR)/ssht
SSHTLIB  = $(SSHTDIR)/lib/c
SSHTLIBNM= ssht
SSHTSRC  = $(SSHTDIR)/src/c
SSHTBIN  = $(SSHTDIR)/bin/c
SSHTOBJ  = $(SSHTSRC)
SSHTINC  = $(SSHTDIR)/include/c
SSHTDOC  = $(SSHTDIR)/doc/c

FFTWDIR      = $(PROGDIR)/fftw
FFTWINC	     = $(FFTWDIR)/include
FFTWLIB      = $(FFTWDIR)/lib
FFTWLIBNM    = fftw3

SSHTSRCMAT	= $(SSHTDIR)/src/matlab
SSHTOBJMAT  	= $(SSHTSRCMAT)
SSHTOBJMEX  	= $(SSHTSRCMAT)


# ======== SOURCE LOCATIONS ========

vpath %.c $(SSHTSRC)
vpath %.h $(SSHTSRC)
vpath %_mex.c $(SSHTSRCMAT)


# ======== FFFLAGS ========

FFLAGS  = -I$(FFTWINC) -I$(SSHTINC)


# ======== LDFLAGS ========

LDFLAGS = -L$(SSHTLIB) -l$(SSHTLIBNM) -L$(FFTWLIB) -l$(FFTWLIBNM) -lm

LDFLAGSMEX = -L$(SSHTLIB) -l$(SSHTLIBNM) $(FFTWLIB)/lib$(FFTWLIBNM).a


# ======== OBJECT FILES TO MAKE ========

SSHTOBJS = $(SSHTOBJ)/ssht_sampling.o    \
           $(SSHTOBJ)/ssht_dl.o          \
           $(SSHTOBJ)/ssht_core.o

SSHTHEADERS = ssht_types.h     \
              ssht_error.h     \
	      ssht_sampling.h  \
	      ssht_dl.h        \
	      ssht_core.h

SSHTOBJSMAT = $(SSHTOBJMAT)/ssht_sampling_mex.o    \
              $(SSHTOBJMAT)/ssht_forward_mex.o     \
              $(SSHTOBJMAT)/ssht_inverse_mex.o

SSHTOBJSMEX = $(SSHTOBJMEX)/ssht_sampling_mex.$(MEXEXT)    \
              $(SSHTOBJMEX)/ssht_forward_mex.$(MEXEXT)     \
              $(SSHTOBJMEX)/ssht_inverse_mex.$(MEXEXT)


# ======== MAKE RULES ========

$(SSHTOBJ)/%.o: %.c $(SSHTHEADERS)
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

.PHONY: test
test: $(SSHTBIN)/ssht_test
$(SSHTBIN)/ssht_test: $(SSHTOBJ)/ssht_test.o $(SSHTLIB)/lib$(SSHTLIBNM).a
	$(CC) $(OPT) $< -o $(SSHTBIN)/ssht_test $(LDFLAGS) 

.PHONY: runtest
runtest: test
	$(SSHTBIN)/ssht_test 64 0

.PHONY: default
default: test

.PHONY: all
all: lib test matlab


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
	rm -f $(SSHTOBJMAT)/*.o
	rm -f $(SSHTOBJMEX)/*.$(MEXEXT)

.PHONY: tidy
tidy:
	rm -f *~ 

.PHONY: cleanall
cleanall: clean cleandoc