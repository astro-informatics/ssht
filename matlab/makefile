


MLAB 		= /Applications/MATLAB_R2010b.app
#MLAB 		= /usr/local/MATLAB/R2010b
MLABINC		= ${MLAB}/extern/include
MLABLIB		= ${MLAB}/extern/lib

MXSFX 		= mexmaci64
#MXSFX 		= mexa64
MEX 		= ${MLAB}/bin/mex
MFLAGS 		= -cxx


#FC = /usr/bin/gfortran-4.3 -m64 -fPIC
CC = gcc #-Wall -pedantic -no-cpp-precomp
C_OPTIMIZE_SWITCH = -O2 ## For GCC
INC = -I${MLABINC} -I../include/c
LIB = -L${MLABLIB} -L../lib/c -lssht
CFLAGS = ${C_OPTIMISE_SWITCH} ${INC} ${LIB}

all: arrayProduct.${MXSFX} ssht_sampling_mex.${MXSFX}

arrayProduct.${MXSFX}: arrayProduct.o
	${MEX} ${MFLAGS} arrayProduct.o -L../lib/c -lssht

ssht_sampling_mex.${MXSFX}: ssht_sampling_mex.o
	${MEX} ${MFLAGS} ssht_sampling_mex.o -L../lib/c -lssht


%.o : %.c
	${CC} -c ${CFLAGS} $<

#timestwo.o: timestwo.F
#	${FC} -c ${CFLAGS} timestwo.F	

clean:
	rm -f arrayProduct.o arrayProduct.${MXSFX}
	rm -f ssht_sampling_mex.o ssht_sampling_mex.${MXSFX}