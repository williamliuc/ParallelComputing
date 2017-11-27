ALL = test1
CC=icc
LD=icc
CFLAGS= -O3 -qopenmp -Ofast -ipo
MICFLAGS= -xMIC-AVX512 -no-fma

test1: LUSolve.o LE_SymSprsMat.o LE_FBackwardSym.o
	${LD} ${CFLAGS} LUSolve.o LE_SymSprsMat.o LE_FBackwardSym.o -o test1
LUSolve.o: LUSolve.cpp LE_SymSprsMatDef.h LE_SymSprsMatFunc.h
	${CC} ${CFLAGS} -c LUSolve.cpp
LE_SymSprsMat.o: LE_SymSprsMat.cpp LE_SymSprsMatDef.h LE_SymSprsMatFunc.h
	${CC} ${CFLAGS} -c LE_SymSprsMat.cpp
LE_FBackwardSym.o: LE_FBackwardSym.cpp LE_SymSprsMatDef.h LE_SymSprsMatFunc.h
	${CC} ${CFLAGS} ${MICFLAGS} -c LE_FBackwardSym.cpp

clean:
	rm -f $(ALL) *.o
tidy:
	rm -f $(ALL) *.o *.optrpt  log.txt

