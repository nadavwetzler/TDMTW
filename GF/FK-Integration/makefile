ALL: FKRPROG wvint9

FKRPROG: FKRPROG.f
	gfortran -m64 -o FKRPROG FKRPROG.f -lm
	\mv FKRPROG ../BIN

wvint9: wvint9.f cwrite.c
	gcc -m64 -c cwrite.c
	gfortran -m64 -o wvint9 wvint9.f cwrite.o -lm
	\mv wvint9 ../BIN

clean:
	\rm *.o
