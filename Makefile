include ${ATOM_ROOT}/platform/build/make.inc.${ATOM_PLATFORM}

export EXTRA_LIBS=

all:
	cd dfftpack ; make
	cd Direct ; make

clean:
	cd dfftpack ; make clean 
	cd Direct ; make clean 

