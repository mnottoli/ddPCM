#
#   Makefile for ddCOSMO
#
#RunF77 = ifort
#FFLAGS = -O3 -xHost -qopenmp -mkl
RunF77 = gfortran
FFLAGS = -O3 -march=native -fopenmp -llapack -lblas
#RunF77 = pgfortran
#FFLAGS = -O3 -mp

MODS   = ddcosmo.o ddpcm_lib.o
OBJS   = mkrhs.o llgnew.o main.o ddcosmo.o ddpcm_lib.o forces_dd.o efld.o\
	matvec.o cosmo.o jacobi_diis.o
#
all:    $(MODS) $(OBJS)
	$(RunF77) $(FFLAGS) -o main.exe $(OBJS)
#
%.o: %.f
	$(RunF77) $(FFLAGS) -c $*.f
%.o: %.f90
	$(RunF77) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(OBJS) *.exe *.mod
