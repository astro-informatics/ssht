# ======== COMPILER ========

FC      = gfortran
#FC      = nagfor
#FC      = g95

ifeq ($(FC),nagfor)
  OPTNAGFOR = -w=x95 -DNAGFOR
endif
ifeq ($(FC),gfortran)
  OPTGFORTRAN = -m64
endif

OPT = $(OPTNAGFOR) $(OPTGFORTRAN) -O3 #-g3 -ggdb


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


# ======== FFFLAGS ========

FFLAGS  = -I$(SSHTINC)


# ======== LDFLAGS ========

LDFLAGS = -L$(SSHTLIB) -l$(SSHTLIBNM) \
          -L$(FFTWLIB) -l$(FFTWLIBNM)


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
          $(SSHTINC)/ssht_error_mod.o    \
          $(SSHTINC)/ssht_dl_mod.o       \
          $(SSHTINC)/ssht_sampling_mod.o \
          $(SSHTINC)/ssht_core_mod.o   


# ======== MAKE RULES ========

default: lib

all:     lib prog test

lib:	 $(SSHTLIB)/lib$(SSHTLIBNM).a

test:    $(SSHTBIN)/ssht_test

runtest: test
	$(SSHTBIN)/ssht_test 64 0

prog:    $(SSHTBIN)/ssht_forward $(SSHTBIN)/ssht_inverse

$(SSHTINC)/%.o: $(SSHTSRC)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 
	mv *.mod $(SSHTINC)

$(SSHTINC)/ssht_test.o:     $(SSHTPROG)/ssht_test.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 

$(SSHTINC)/%.o: $(SSHTPROG)/%.f90
	$(FC) $(FFLAGS) $(FFLAGSPROG) $(PPFLAGS) -c $< -o $@ 


# Library

$(SSHTLIB)/lib$(SSHTLIBNM).a: $(SSHTOBJ)
	ar -r $(SSHTLIB)/lib$(SSHTLIBNM).a $(SSHTOBJ)


# Documentation
.PHONY: doc
doc:	
	./f90doc_fpp $(SSHTSRC)/*.f90
	./f90doc_fpp $(SSHTPROG)/*.f90
	mv *.html $(SSHTDOC)/.
	./addstyle $(SSHTDOC)/ssht_*

.PHONY: cleandoc
cleandoc:
	rm -f $(SSHTDOC)/ssht_*.html
	rm -f $(SSHTDOC)/ran2_dp.html

# Cleaning up

.PHONY: clean
clean:	tidy
	rm -f $(SSHTINC)/*.mod
	rm -f $(SSHTINC)/*.o
	rm -f $(SSHTLIB)/lib$(SSHTLIBNM).a
	rm -f $(SSHTBIN)/*

.PHONY: tidy
tidy:
	rm -f *.mod
	rm -f $(SSHTSRC)/*~ 
	rm -f $(SSHTPROG)/*~ 


# Module dependencies

$(SSHTINC)/ssht_types_mod.o: $(SSHTSRC)/ssht_types_mod.f90
$(SSHTINC)/ssht_error_mod.o: $(SSHTSRC)/ssht_error_mod.f90          \
                           $(SSHTINC)/ssht_types_mod.o
$(SSHTINC)/ssht_dl_mod.o:    $(SSHTSRC)/ssht_dl_mod.f90             \
                           $(SSHTINC)/ssht_types_mod.o
$(SSHTINC)/ssht_sampling_mod.o:  $(SSHTSRC)/ssht_sampling_mod.f90   \
                           $(SSHTINC)/ssht_types_mod.o              \
                           $(SSHTINC)/ssht_error_mod.o
$(SSHTINC)/ssht_core_mod.o:  $(SSHTSRC)/ssht_core_mod.f90           \
                           $(SSHTINC)/ssht_types_mod.o              \
                           $(SSHTINC)/ssht_error_mod.o              \
                           $(SSHTINC)/ssht_sampling_mod.o           \
                           $(SSHTINC)/ssht_dl_mod.o       


# Program dependencies and compilation

$(SSHTINC)/ssht_test.o:     $(SSHTPROG)/ssht_test.f90 lib
$(SSHTBIN)/ssht_test:       $(SSHTINC)/ssht_test.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_test                        \
	$(SSHTINC)/ssht_test.o                         \
	$(LDFLAGS) $(PPFLAGS)

$(SSHTINC)/ssht_forward.o:     $(SSHTPROG)/ssht_forward.f90 lib
$(SSHTBIN)/ssht_forward:       $(SSHTINC)/ssht_forward.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_forward                     \
	$(SSHTINC)/ssht_forward.o                      \
	$(LDFLAGS) $(PPFLAGS)

$(SSHTINC)/ssht_inverse.o:     $(SSHTPROG)/ssht_inverse.f90 lib
$(SSHTBIN)/ssht_inverse:       $(SSHTINC)/ssht_inverse.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_inverse                     \
	$(SSHTINC)/ssht_inverse.o                      \
	$(LDFLAGS) $(PPFLAGS)


