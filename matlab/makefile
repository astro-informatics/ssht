


MLAB 		= /Applications/MATLAB_R2010a.app
MLAB 		= /usr/local/MATLAB/R2010b
MLABINC		= ${MLAB}/extern/include
MLABLIB		= ${MLAB}/extern/lib

#MXSFX 		= mexmaci64
MXSFX 		= mexa64
MEX 		= ${MLAB}/bin/mex
MFLAGS 		= -fortran


FC = /usr/bin/gfortran-4.3 -m64 -fPIC
#cc -Wall -pedantic -no-cpp-precomp
C_OPTIMIZE_SWITCH = -O2 ## For GCC
INC = -I${MLABINC} -I../include
LIB = -L${MLABLIB} -L../lib -lssht
CFLAGS = ${INC} ${LIB}

all: timestwo.${MXSFX}

timestwo.${MXSFX}: timestwo.o
	${MEX} ${MFLAGS} timestwo.o -L../lib -lssht



timestwo.o: timestwo.F
	${FC} -c ${CFLAGS} timestwo.F	

clean:
	rm -f timestwo.o timestwo.${MXSFX}