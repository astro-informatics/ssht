


MLAB 		= /Applications/MATLAB_R2010a.app
MLABINC		= ${MLAB}/extern/include
MLABLIB		= ${MLAB}/extern/lib

MXSFX 		= mexmaci64
MEX 		= ${MLAB}/bin/mex
MFLAGS 		= -fortran -v


FC = gfortran -m64
#cc -Wall -pedantic -no-cpp-precomp
C_OPTIMIZE_SWITCH = -O2 ## For GCC
INC = -I${MLABINC} 
LIB = -L${MLABLIB} 
CFLAGS = ${INC} ${LIB}

all: timestwo.${MXSFX}

timestwo.${MXSFX}: timestwo.o
	${MEX} ${MFLAGS} timestwo.o



timestwo.o: timestwo.F
	${FC} -c ${CFLAGS} timestwo.F	

clean:
	rm -f timestwo.o timestwo.${MXSFX}