include ${ATOM_ROOT}/platform/build/make.inc.${ATOM_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG} 
else
   FFLAGS=${FOPT}
endif

EXEC=gem_main

OBJECTS = gem_fft_wrapper.o \
          gem_pputil.o \
          gem_com.o \
          gem_equil.o \
          gem_fcnt.o \
          gem_gkps_adi.o \
          gem_outd.o \
          $(EXEC).o

.SUFFIXES : .o .f90

all: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(LMATH) ../dfftpack/libdfftpack.a

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o *.lst $(EXEC)
	cd ${GACODE_ROOT}/modules ; rm -f gem*.mod
