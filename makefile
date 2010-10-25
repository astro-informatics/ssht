# ======== COMPILER ========

#FC      = f95
FC      = gfortran
#FC      = g95

ifeq ($(FC),f95)
  OPTF95 = -w=x95
endif

OPT = $(OPTF95) -m64


# ======== LINKS ========

PROGDIR = /Users/mcewen/src

SSHTDIR  = $(PROGDIR)/ssht
SSHTLIB  = $(SSHTDIR)/lib
SSHTLIBNM= ssht
SSHTINC  = $(SSHTDIR)/include
SSHTSRC  = $(SSHTDIR)/src/mod
SSHTPROG = $(SSHTDIR)/src/prog
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

CFITSIOLIB   = $(PROGDIR)/cfitsio/lib
CFITSIOLIBNM = cfitsio


# ======== FFFLAGS ========

FFLAGS  = -I$(SSHTINC)
FFLAGSPROG = -I$(HPIXINC) -I$(S2INC)


# ======== LDFLAGS ========

LDFLAGS = -L$(SSHTLIB) -l$(SSHTLIBNM) \
          -L$(FFTWLIB) -l$(FFTWLIBNM) \
          -L$(CFITSIOLIB) -l$(CFITSIOLIBNM) 

LDFLAGSPROG = -L$(S2LIB) -l$(S2LIBNM) \
           -L$(HPIXLIB) -l$(HPIXLIBNM) 


# ======== PPFLAGS ========

ifeq ($(FC),f95)
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
          $(SSHTINC)/ssht_fileio_mod.o  \
          $(SSHTINC)/ssht_core_mod.o   


# ======== MAKE RULES ========

default: lib

all:     lib prog test

lib:	 $(SSHTLIB)/lib$(SSHTLIBNM).a

test:    $(SSHTBIN)/ssht_test

runtest: test
	$(SSHTBIN)/ssht_test 64

prog:    $(SSHTBIN)/ssht_wav2sky $(SSHTBIN)/ssht_analysis $(SSHTBIN)/ssht_synthesis $(SSHTBIN)/ssht_wavplot $(SSHTBIN)/ssht_mat2fits $(SSHTBIN)/ssht_fits2mat

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
$(SSHTINC)/ssht_fileio_mod.o:  $(SSHTSRC)/ssht_fileio_mod.f90   \
                           $(SSHTINC)/ssht_types_mod.o    \
                           $(SSHTINC)/ssht_error_mod.o    \
                           $(SSHTINC)/ssht_core_mod.o    
$(SSHTINC)/ssht_core_mod.o:  $(SSHTSRC)/ssht_core_mod.f90   \
                           $(SSHTINC)/ssht_types_mod.o    \
                           $(SSHTINC)/ssht_error_mod.o    \
                           $(SSHTINC)/ssht_dl_mod.o       


# Program dependencies and compilation

$(SSHTINC)/ssht_test.o:     $(SSHTPROG)/ssht_test.f90 lib
$(SSHTBIN)/ssht_test:       $(SSHTINC)/ssht_test.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_test                          \
	$(SSHTINC)/ssht_test.o $(LDFLAGS) $(PPFLAGS)

$(SSHTINC)/ssht_wav2sky.o:     $(SSHTPROG)/ssht_wav2sky.f90 lib
$(SSHTBIN)/ssht_wav2sky:       $(SSHTINC)/ssht_wav2sky.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_wav2sky                       \
	$(SSHTINC)/ssht_wav2sky.o $(LDFLAGSPROG) $(LDFLAGS) $(PPFLAGS)

$(SSHTINC)/ssht_analysis.o:     $(SSHTPROG)/ssht_analysis.f90 lib
$(SSHTBIN)/ssht_analysis:       $(SSHTINC)/ssht_analysis.o
	$(FC)                                          \
	 -o $(SSHTBIN)/ssht_analysis                     \
	$(SSHTINC)/ssht_analysis.o $(LDFLAGSPROG) $(LDFLAGS) $(PPFLAGS)               

$(SSHTINC)/ssht_synthesis.o:     $(SSHTPROG)/ssht_synthesis.f90 lib
$(SSHTBIN)/ssht_synthesis:       $(SSHTINC)/ssht_synthesis.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_synthesis                     \
	$(SSHTINC)/ssht_synthesis.o $(LDFLAGSPROG) $(LDFLAGS) $(PPFLAGS)                

$(SSHTINC)/ssht_wavplot.o:     $(SSHTPROG)/ssht_wavplot.f90 lib
$(SSHTBIN)/ssht_wavplot:       $(SSHTINC)/ssht_wavplot.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_wavplot                       \
	$(SSHTINC)/ssht_wavplot.o $(LDFLAGSPROG) $(LDFLAGS) $(PPFLAGS)                       

$(SSHTINC)/ssht_mat2fits.o:     $(SSHTPROG)/ssht_mat2fits.f90 lib
$(SSHTBIN)/ssht_mat2fits:       $(SSHTINC)/ssht_mat2fits.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_mat2fits                       \
	$(SSHTINC)/ssht_mat2fits.o $(LDFLAGSPROG) $(LDFLAGS) $(PPFLAGS)  

$(SSHTINC)/ssht_fits2mat.o:     $(SSHTPROG)/ssht_fits2mat.f90 lib
$(SSHTBIN)/ssht_fits2mat:       $(SSHTINC)/ssht_fits2mat.o
	$(FC)                                          \
	-o $(SSHTBIN)/ssht_fits2mat                       \
	$(SSHTINC)/ssht_fits2mat.o $(LDFLAGSPROG) $(LDFLAGS) $(PPFLAGS)


