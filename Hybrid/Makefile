SRCS =	erf.f gem_com.f equil.f gem.f outd.f fcnt.f fft_wrapper.f

OBJS =	erf.o gem_com.o equil.o gem.o outd.o fcnt.o fft_wrapper.o

LIBS = $(DFFTPACK) -mkl
PLIB = pputil.o

F90 = ftn

#OPT = -FR -r8 -O2 -g -traceback
#OPT = -FR -r8 -heap-arrays -O2 -g -traceback -check bounds
OPT = -FR -r8 -heap-arrays

LDFLAGS = 

#all : gem

gem: equil.o gem.o outd.o fcnt.o pputil.o gem_com.o erf.o fft_wrapper.o
	$(F90)  -o gem $(OPT) $(OBJS) $(PLIB) $(LIBS) 

pputil.o: pputil.f fft_wrapper.o
	$(F90) -c $(OPT) pputil.f

gem_com.o: gem_com.f pputil.o
	$(F90) -c $(OPT) gem_com.f

equil.o: equil.f pputil.o
	$(F90) -c $(OPT)  equil.f

gem.o: gem.f fft_wrapper.o pputil.o gem_com.o equil.o
	$(F90) -c $(OPT) gem.f

outd.o: outd.f fft_wrapper.o pputil.o gem_com.o equil.o
	$(F90) -c $(OPT) outd.f

fcnt.o: fcnt.f
	$(F90) -c $(OPT) fcnt.f

erf.o: erf.f
	$(F90) -c $(OPT) -FI erf.f

fft_wrapper.o: fft_wrapper.f
	$(F90) -c $(OPT) fft_wrapper.f
