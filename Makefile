#
#   Makefile for ddCOSMO
#
#RunF77 = ifort
#FFLAGS = -O3 -xHost -qopenmp
#RunF77 = gfortran
#FFLAGS = -O3 -march=native -llapack -lblas -fopenmp
RunF77 = pgfortran
FFLAGS = -O3 -mp -llapack -lblas

MODS   = ddcosmo.o newschwarz.o ddpcm_lib.o
OBJS   = mkrhs.o llgnew.o main.o ddcosmo.o ddpcm_lib.o forces_dd.o efld.o\
	matvec.o cosmo.o jacobi_diis.o newschwarz.o
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
