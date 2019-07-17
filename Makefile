#
#   Makefile for ddCOSMO
#
#RunF77 = ifort
#FFLAGS = -O3 -xHost -qopenmp
RunF77 = gfortran
FFLAGS = -O3 -march=native -llapack -lblas
#RunF77 = pgfortran
#FFLAGS = -O3 -mp

MODS   = bessel.o ddcosmo.o ddlpb.o
OBJS   = mkrhs.o llgnew.o ddlpb.o ddcosmo.o forces_dd.o efld.o\
	matvec.o cosmo.o jacobi_diis.o main.o bessel.o gmres.o
#
all:    $(MODS) $(OBJS)
	$(RunF77) $(FFLAGS) -o main.exe $(OBJS)
%.o: %.f
	$(RunF77) $(FFLAGS) -c $*.f
%.o: %.f90
	$(RunF77) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(OBJS) *.exe *.mod
