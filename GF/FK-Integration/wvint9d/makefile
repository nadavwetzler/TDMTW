ALL: wvint9d

wvint9: wvint9.f cwrite.c
	gcc -m64 -c cwrite.c
	gfortran -m64 -o wvint9d wvint9d.f cwrite.o -lm
	\mv wvint9d ../../BIN

clean:
	\rm *.o
