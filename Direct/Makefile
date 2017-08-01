SRCS =	erf.f gem_com.f90 equil.f90 gem.f90 outd.f90 fcnt.f fft_wrapper.f

OBJS =	erf.o gem_com.o equil.o gem.o outd.o fcnt.o fft_wrapper.o

LIBS = $(DFFTPACK) -mkl
PLIB = pputil.o

F90 = ftn
#OPT = -FR -r8 -heap-arrays -O2 -g -traceback -check bounds
OPT = -FR -r8 -O2
LDFLAGS = 

#all : gem

gem: equil.o gem.o outd.o fcnt.o pputil.o gem_com.o erf.o fft_wrapper.o
	$(F90)  -o gem $(OPT) $(OBJS) $(PLIB) $(LIBS) 

pputil.o: pputil.f
	$(F90) -c $(OPT) pputil.f

gem_com.o: gem_com.f90 pputil.o
	$(F90) -c $(OPT) gem_com.f90

equil.o: equil.f90 pputil.o
	$(F90) -c $(OPT) equil.f90

gem.o: gem.f90 fft_wrapper.o pputil.o gem_com.o equil.o
	$(F90) -c $(OPT) gem.f90

outd.o: outd.f90 fft_wrapper.o pputil.o gem_com.o equil.o
	$(F90) -c $(OPT) outd.f90

fcnt.o: fcnt.f
	$(F90) -c $(OPT) fcnt.f

erf.o: erf.f
	$(F90) -c $(OPT) -FI erf.f

fft_wrapper.o: fft_wrapper.f
	$(F90) -c $(OPT) fft_wrapper.f
