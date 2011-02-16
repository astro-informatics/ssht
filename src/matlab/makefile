


MLAB 		= /Applications/MATLAB_R2010a.app
#MLAB 		= /usr/local/MATLAB/R2010b
MLABINC		= ${MLAB}/extern/include
MLABLIB		= ${MLAB}/extern/lib

MXSFX 		= mexmaci64
#MXSFX 		= mexa64
MEX 		= ${MLAB}/bin/mex
MFLAGS 		= -cxx


#FC = /usr/bin/gfortran-4.3 -m64 -fPIC
CC = gcc -Wall #-no-cpp-precomp
C_OPTIMIZE_SWITCH = # -O3 ## For GCC
INC = -I${MLABINC} -I../../include/c -I../../../fftw/include
LIB = -L${MLABLIB} -L../../lib/c -lssht ../../../fftw/lib/libfftw3.a
#LIB = -L${MLABLIB} -L../lib/c -lssht -L../../fftw/lib -lfftw3
CFLAGS = ${C_OPTIMISE_SWITCH} ${INC} 

all: ssht_sampling_mex.${MXSFX} ssht_inverse_mex.${MXSFX} ssht_forward_mex.${MXSFX}

ssht_sampling_mex.${MXSFX}: ssht_sampling_mex.o
	${MEX} ${MFLAGS} ssht_sampling_mex.o ${LIB}

ssht_inverse_mex.${MXSFX}: ssht_inverse_mex.o
	${MEX} ${MFLAGS} ssht_inverse_mex.o ${LIB}

ssht_forward_mex.${MXSFX}: ssht_forward_mex.o
	${MEX} ${MFLAGS} ssht_forward_mex.o ${LIB}


%.o : %.c
	${CC} -c ${CFLAGS} $< ${INC}

clean:
	rm -f ssht_sampling_mex.o ssht_sampling_mex.${MXSFX}
	rm -f ssht_inverse_mex.o ssht_inverse_mex.${MXSFX}
	rm -f ssht_forward_mex.o ssht_forward_mex.${MXSFX}
