ALL: FKRPROG wvint9 wvint9d

FKRPROG: FKRPROG.f
	gfortran -m64 -o FKRPROG FKRPROG.f -lm
	\mv FKRPROG ../BIN_Linux

wvint9: wvint9.f cwrite.c
	gcc -m64 -c cwrite.c
	gfortran -m64 -o wvint9 wvint9.f cwrite.o -lm
	\mv wvint9 ../BIN_Linux

wvint9d: wvint9d.f cwrite.c
	gcc -m64 -c cwrite.c
	gfortran -m64 -o wvint9d wvint9d.f cwrite.o -lm
	\mv wvint9d ../BIN_Linux

clean:
	\rm *.o
