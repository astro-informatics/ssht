# ======== COMPILER ========

CC      = gcc
OPT	= -Wall -g
#OPT	= -Wall -O3


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

vpath %.c $(SSHTSRC)
vpath %.h $(SSHTSRC)

# ======== FFFLAGS ========

FFLAGS  = -I$(FFTWINC) -I$(SSHTINC)



# ======== LDFLAGS ========

LDFLAGS = -L$(SSHTLIB) -l$(SSHTLIBNM) -L$(FFTWLIB) -l$(FFTWLIBNM) -lm


# ======== OBJECT FILES TO MAKE ========

SSHTOBJS = $(SSHTOBJ)/ssht_sampling.o    \
          $(SSHTOBJ)/ssht_dl.o          \
          $(SSHTOBJ)/ssht_core.o

SSHTHEADERS = ssht_types.h     \
              ssht_error.h     \
	      ssht_sampling.h  \
	      ssht_dl.h        \
	      ssht_core.h


# ======== MAKE RULES ========

$(SSHTOBJ)/%.o: %.c $(SSHTHEADERS)
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@

.PHONY: test
test: $(SSHTBIN)/ssht_test
$(SSHTBIN)/ssht_test: $(SSHTOBJ)/ssht_test.o $(SSHTLIB)/lib$(SSHTLIBNM).a
	$(CC) $(OPT) $< -o $(SSHTBIN)/ssht_test $(LDFLAGS) 

.PHONY: default
default: all

.PHONY: all
all: lib test


# Library

.PHONY: lib
lib: $(SSHTLIB)/lib$(SSHTLIBNM).a
$(SSHTLIB)/lib$(SSHTLIBNM).a: $(SSHTOBJS)
	ar -r $(SSHTLIB)/lib$(SSHTLIBNM).a $(SSHTOBJS)


# Documentation 

.PHONY: doc
doc:
	doxygen $(SSHTSRC)/doxygen.config
.PHONY: cleandoc
cleandoc:
	rm -rf $(SSHTDOC)/html


# Cleaning up

.PHONY: clean
clean:	tidy
	rm -f $(SSHTOBJ)/*.o
	rm -f $(SSHTLIB)/lib$(SSHTLIBNM).a
	rm -f $(SSHTBIN)/ssht_test

.PHONY: tidy
tidy:
	rm -f *~ 

.PHONY: cleanall
cleanall: clean cleandoc