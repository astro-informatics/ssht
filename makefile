# ======== COMPILER ========

#FC      = nagfor
#FC	= /usr/bin/gfortran-4.3
FC      = gfortran
#FC      = g95

ifeq ($(FC),nagfor)
  OPTNAGFOR = -w=x95 -DNAGFOR
endif
ifeq ($(FC),gfortran)
  OPTGFORTRAN = -m64
endif

OPT = $(OPTNAGFOR)  $(OPTGFORTRAN) #-g3 -ggdb #3-ggdb #  -m64 -O3


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

FFTWLIB      = $(PROGDIR)/fftw/lib
#FFTWLIB      = $(PROGDIR)/fftw-3.2.2_m32/lib
FFTWLIBNM    = fftw3

HPIXDIR = $(PROGDIR)/Healpix
HPIXLIB = $(HPIXDIR)/lib
HPIXLIBNM= healpix
HPIXINC = $(HPIXDIR)/include

S2DIR  = $(PROGDIR)/s2
S2LIB  = $(S2DIR)/lib
S2LIBNM= s2
S2INC  = $(S2DIR)/include
S2DOC  = $(S2DIR)/doc

#CFITSIOLIB   = $(PROGDIR)/cfitsio/lib
#CFITSIOLIBNM = cfitsio


# ======== FFFLAGS ========

FFLAGS  = -I$(SSHTINC)
FFLAGSPROG = -I$(HPIXINC) -I$(S2INC)


# ======== LDFLAGS ========

LDFLAGS = -L$(SSHTLIB) -l$(SSHTLIBNM) \
          -L$(FFTWLIB) -l$(FFTWLIBNM)
#         -L$(CFITSIOLIB) -l$(CFITSIOLIBNM) 

LDFLAGSPROG = -L$(S2LIB) -l$(S2LIBNM) \
           -L$(HPIXLIB) -l$(HPIXLIBNM) 


# ======== PPFLAGS ========

ifeq ($(FC),nagfor)
  PPFLAGS = -fpp $(OPT)
else ifeq ($(FC),g95)
  PPFLAGS = -cpp $(OPT)
else ifeq ($(FC),gfortran)
  PPFLAGS = -x f95-cpp-input $(OPT)
endif


# ======== OBJECT FILES TO MAKE ========

SSHTOBJ = $(SSHTINC)/ssht_types_mod.o    \
          $(SSHTINC)/ssht_error_mod.o   \
          $(SSHTINC)/ssht_dl_mod.o      \
          $(SSHTINC)/ssht_sampling_mod.o  \
          $(SSHTINC)/ssht_core_mod.o   


# ======== MAKE RULES ========

default: lib

all:     lib prog test

lib:	 $(SSHTLIB)/lib$(SSHTLIBNM).a

test:    $(SSHTBIN)/ssht_test $(SSHTBIN)/ssht_trapani

runtest: test
	$(SSHTBIN)/ssht_test 64 0

prog:    $(SSHTBIN)/ssht_forward $(SSHTBIN)/ssht_inverse

$(SSHTINC)/%.o: $(SSHTSRC)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 
	mv *.mod $(SSHTINC)

$(SSHTINC)/ssht_test.o:     $(SSHTPROG)/ssht_test.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 

$(SSHTINC)/%.o: $(SSHTPROG)/%.f90
	$(FC) $(FFLAGS) $(FFLAGSPROG) $(PPFLAGS) -I$(S2INC) -c $< -o $@ 


# Library

$(SSHTLIB)/lib$(SSHTLIBNM).a: $(SSHTOBJ)
	ar -r $(SSHTLIB)/lib$(SSHTLIBNM).a $(SSHTOBJ)


# Documentation

docs:
	./f90doc_fpp $(SSHTSRC)/*.f90
	./f90doc_fpp $(SSHTPROG)/*.f90
	./ln_multi $(S2DOC)/s2_*
	./ln_multi $(S2DOC)/index_s2.html
	mv *.html $(SSHTDOC)/.
	./addstyle $(SSHTDOC)/ssht_*

cleandocs:
	rm -f $(SSHTDOC)/ssht_*.html
	rm -f $(SSHTDOC)/gasdev2_dp.html $(SSHTDOC)/ran2_dp.html
	rm -f $(SSHTDOC)/s2_*.html $(SSHTDOC)/index_s2.html

# Cleaning up

clean:	tidy
	rm -f $(SSHTINC)/*.mod
	rm -f $(SSHTINC)/*.o
	rm -f $(SSHTLIB)/lib$(SSHTLIBNM).a
	rm -f $(SSHTBIN)/*

tidy:
	rm -f *.mod
	rm -f $(SSHTSRC)/*~ 
	rm -f $(SSHTPROG)/*~ 


# Module dependencies

$(SSHTINC)/ssht_types_mod.o: $(SSHTSRC)/ssht_types_mod.f90
$(SSHTINC)/ssht_error_mod.o: $(SSHTSRC)/ssht_error_mod.f90  \
                           $(SSHTINC)/ssht_types_mod.o
$(SSHTINC)/ssht_dl_mod.o:    $(SSHTSRC)/ssht_dl_mod.f90     \
                           $(SSHTINC)/ssht_types_mod.o
$(SSHTINC)/ssht_sampling_mod.o:  $(SSHTSRC)/ssht_sampling_mod.f90   \
                           $(SSHTINC)/ssht_types_mod.o    \
                           $(SSHTINC)/ssht_error_mod.o
$(SSHTINC)/ssht_core_mod.o:  $(SSHTSRC)/ssht_core_mod.f90   \
                           $(SSHTINC)/ssht_types_mod.o    \
                           $(SSHTINC)/ssht_error_mod.o    \
                           $(SSHTINC)/ssht_sampling_mod.o    \
                           $(SSHTINC)/ssht_dl_mod.o       


# Program dependencies and compilation

$(SSHTINC)/ssht_test.o:     $(SSHTPROG)/ssht_test.f90 lib
$(SSHTBIN)/ssht_test:       $(SSHTINC)/ssht_test.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_test                          \
	$(SSHTINC)/ssht_test.o $(LDFLAGS) $(PPFLAGS)

$(SSHTINC)/ssht_trapani.o:     $(SSHTPROG)/ssht_trapani.f90 lib
$(SSHTBIN)/ssht_trapani:       $(SSHTINC)/ssht_trapani.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_trapani                          \
	$(SSHTINC)/ssht_trapani.o $(LDFLAGS) $(PPFLAGS)

$(SSHTINC)/ssht_forward.o:     $(SSHTPROG)/ssht_forward.f90 lib
$(SSHTBIN)/ssht_forward:       $(SSHTINC)/ssht_forward.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_forward                       \
	$(SSHTINC)/ssht_forward.o $(LDFLAGSPROG) $(LDFLAGS) $(PPFLAGS)

$(SSHTINC)/ssht_inverse.o:     $(SSHTPROG)/ssht_inverse.f90 lib
$(SSHTBIN)/ssht_inverse:       $(SSHTINC)/ssht_inverse.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_inverse                       \
	$(SSHTINC)/ssht_inverse.o $(LDFLAGSPROG) $(LDFLAGS) $(PPFLAGS)


