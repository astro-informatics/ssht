# ======== COMPILER ========

CC      = gcc
OPT	= -Wall -g # -O3 #-g
# OPT 	= -m64 -O3

# ======== LINKS ========

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
  PROGDIR = /home/jdm57/src
endif
ifeq ($(UNAME), Darwin)
  PROGDIR = /Users/mcewen/src
endif

SSHTDIR  = $(PROGDIR)/ssht
SSHTLIB  = $(SSHTDIR)/lib
SSHTLIBNM= ssht
SSHTINC  = $(SSHTDIR)/include
SSHTSRC  = $(SSHTDIR)/src/f90/mod
SSHTPROG = $(SSHTDIR)/src/f90/prog
SSHTBIN  = $(SSHTDIR)/bin
SSHTDOC  = $(SSHTDIR)/doc

FFTWDIR      = $(PROGDIR)/fftw
FFTWINC	     = $(FFTWDIR)/include
FFTWLIB      = $(FFTWDIR)/lib
FFTWLIBNM    = fftw3


# ======== FFFLAGS ========

FFLAGS  = -I$(FFTWINC)


# ======== LDFLAGS ========

LDFLAGS = -L$(FFTWLIB) -l$(FFTWLIBNM) -lm


# ======== OBJECT FILES TO MAKE ========

SSHTOBJ = ssht_sampling.o    \
          ssht_dl.o          \
          ssht_core.o


# ======== MAKE RULES ========

%.o: %.c *.h
	$(CC) $(OPT) $(FFLAGS) -c $< -o $@ 

ssht_test: ssht_test.o $(SSHTOBJ)
	$(CC) $(OPT) ssht_test.o $(SSHTOBJ) -o ssht_test $(LDFLAGS)

default: ssht_test

all: lib ssht_test


# Library

lib: libssht.a
libssht.a: $(SSHTOBJ)
	ar -r libssht.a $(SSHTOBJ)


# Cleaning up

clean:	tidy
	rm -f *.o
	rm -f libssht.a
	rm -f ssht_test

tidy:
	rm -f *~ 