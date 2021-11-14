ALLFLAGS= -m64
CFLAGS	= 
FFLAGS 	= 

ALL	= get_resp sac2helm window bin2sac sac2bin fromHelm mkHelm putmech putmt decon
PROGLIST= get_resp sac2helm window bin2sac sac2bin fromHelm mkHelm putmech putmt 

all:	$(ALL)

get_resp: get_resp.c getenv_utils.c getpar.c
	gcc $(ALLFLAGS) -o get_resp get_resp.c getenv_utils.c getpar.c
	\mv get_resp ../../BIN_Linux

sac2helm: sac2helm.c getpar.c
	gcc $(ALLFLAGS) -o sac2helm sac2helm.c getpar.c
	\mv sac2helm ../../BIN_Linux

sac2bin: sac2bin.c getpar.c
	gcc $(ALLFLAGS) -o sac2bin sac2bin.c getpar.c
	\mv sac2bin ../../BIN_Linux

putmech:  putmech_iso.c readhelm.c getpar.c
	gcc $(ALLFLAGS) -g -o putmech putmech_iso.c readhelm.c getpar.c -lm
	\mv putmech ../../BIN_Linux

putmt:  putmt.c readhelm.c getpar.c
	gcc $(ALLFLAGS) -g -o putmt putmt.c readhelm.c getpar.c -lm
	\mv putmt ../../BIN_Linux

window: window.c getpar.c
	gcc $(ALLFLAGS) -o window window.c getpar.c 
	\mv window ../../BIN_Linux


fromHelm: fromHelm_new.c
	gcc $(ALLFLAGS) -o fromHelm fromHelm_new.c
	\mv fromHelm ../../BIN_Linux

mkHelm: mkHelm2.c getpar.c
	gcc $(ALLFLAGS) -o mkHelm mkHelm2.c getpar.c
	\mv mkHelm ../../BIN_Linux

decon: deconnew.c deconnew_subs2.c getpar.c
	gcc $(ALLFLAGS) -o decon deconnew.c deconnew_subs2.c getpar.c -lm
	\mv decon ../../BIN_Linux

bin2sac: bin2sac.c getpar.c
	gcc $(ALLFLAGS) -o bin2sac bin2sac.c getpar.c
	\mv bin2sac ../../BIN_Linux

##Cleanup routine

clean:
	\rm *.o
