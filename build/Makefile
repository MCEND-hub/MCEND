FC = ifort
FFLAGS = -O3 -fopenmp -c -fast  -traceback -i8 -I/opt/intel/oneapi/compiler/2023.2.0/linux/compiler/include/ -I${MKLROOT}/lib/intel64 -pg -g -qmkl=sequential
MKL_LIBS = $(MKLROOT)/lib/intel64
LINK_FLAGS =  -O3 -qopenmp -L${MKLROOT}/lib/intel64 -traceback -Wl,--start-group -lfftw3 -liomp5 -lpthread -lm -ldl -lmkl_intel_lp64 -lmkl_sequential -lmkl_core


BINDIR=../bin
SRCDIR=../src
OBJ = read-input.o generic_data.o akima_interpol.o utils.o combinations.o init.o analysis.o propa.o mcend-2.1.o

all: $(BINDIR)/mcend

$(BINDIR)/mcend: $(OBJ)
	$(FC) -o $(BINDIR)/mcend $(OBJ) $(LINK_FLAGS)

debug: FFLAGS += -O0 -i8 -warn all -check all -traceback -fp-stack-check 
debug: $(BINDIR)/mcend

gfort: FC = gfortran
gfort: FFLAGS = -O3 -c -pg -g
#gfort: FFLAGS = -O1 -c -pg -g
gfort: LINK_FLAGS = -llapack -lblas -lfftw3 -fopenmp -pg -g
gfort: $(BINDIR)/mcend

debug_gfort: FC = gfortran
debug_gfort: FFLAGS = -O0 -c
debug_gfort: LINK_FLAGS = -llapack -lblas -lfftw3 -fopenmp
#debug_gfort: FFLAGS += -O0 -c -g -p -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -Wline-truncation -Wall -fbacktrace -fbounds-check
debug_gfort: FFLAGS += -O0 -c -g -p -ffree-line-length-0 -fcheck=all -Wline-truncation -Wall -fbacktrace -fbounds-check
debug_gfort: $(BINDIR)/mcend

clean:
	rm -f *.o *.mod *.inc *.a

%.o:	$(SRCDIR)/%.f90
	$(FC) $(FFLAGS)  $<



## with hdf5
#CFLAGS	= -O3 -c -g -p -traceback -qopenmp  -I$(HDF5)/include -pg
#LINK_FLAGS = -O3 -traceback -qopenmp -mkl -Wl,-rpath, -lfftw3 $(HDF5)/lib/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a -lz -lrt -lm -pg
#CFLAGS	   = -O3 -c -g -p -traceback -qopenmp
#LINK_FLAGS = -O3 -traceback -qopenmp -mkl -Wl,-rpath, -lfftw3 -pg
# -multiple-processes=2
#CFLAGS =  -O0 -i8 -warn all -check all -traceback -fp-stack-check -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -c -p
#LINK_FLAGS = -O0 -traceback -pg -warn all -check all ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lfftw3 -liomp5 -lpthread -lm -ldl
#LINK_FLAGS = -O3 -traceback -pg ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lfftw3 -liomp5 -lpthread -lm -ldl

## cruntch3/lucas@ubuntu flags
#-multiple-processes=4
#CFLAGS = -O3 -c -fast -qopenmp -xCORE-AVX2 -qopt-matmul  -fpp -traceback -i8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -pg
#CFLAGS = -O3 -c -fast -qopenmp -xCORE-AVX2 -traceback -i8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -pg
#CFLAGS = -O3 -fast -xCORE-AVX2 -traceback -i8 -qopenmp -qopenmp-link=static -threads -qopt-matmul -multiple-processes=4 -fpp -mkl -static-intel -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -pg

#CFLAGS = -O3 -qopenmp -c -fast  -traceback -i8   -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -pg -g -mkl=parallel
#CFLAGS = -O3 -qopenmp -c -fast  -traceback -i8   -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include -pg -g  -mkl=sequential
#MKL_LIBS = $(MKLROOT)/lib/intel64
#LINK_FLAGS = -O3 -qopenmp -traceback ${MKL_LIBS}/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKL_LIBS}/libmkl_intel_ilp64.a ${MKL_LIBS}/libmkl_intel_thread.a ${MKL_LIBS}/libmkl_core.a ${MKL_LIBS}/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lfftw3 -liomp5 -lpthread -lm -ldl
#
#CFLAGS = -O0 -c -warn all -check all -traceback -fp-stack-check -i8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
# LINK_FLAGS = -O3 -traceback -pg ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lfftw3 -liomp5 -lpthread -lm -ldl

#####--GFORTAN--#####
# -O3  -stand f03 -assume realloc_lhs -check all -traceback -warn all -fstack-protector -assume protect_parens -implicitnone
# zora flags, complete debug
#CFLAGS  = -O3 -c  -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan  
# CFLAGS  = -O0 -c -g -p -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
#LINK_FLAGS = -O0 -fdefault-integer-8 -fbacktrace -llapack -lblas -fopenmp -lfftw3 -lgomp -lm -ldl -pg
## for hdf5 compiling
# CFLAGS  = -O3 -c -fbacktrace -fopenmp -I$(HDF5)/include -pg
# LINK_FLAGS = -O3 -fbacktrace -llapack -lblas -fopenmp -lfftw3 -lgomp $(HDF5)/lib/libhdf5_fortran.a $(HDF5)/lib/libhdf5.a -lz -ldl -lm
#
# CFLAGS  = -O0 -c -fopenmp -fdefault-integer-8 -fbacktrace -fstack-arrays -lgomp -pg
# LINK_FLAGS = -O0 -fopenmp -fdefault-integer-8 -fbacktrace -llapack -lblas  -lfftw3 -lgomp -lm -ldl -pg

#CFLAGS  = -O3 -c -fopenmp -fbacktrace
#LINK_FLAGS = -O3 -fopenmp -fbacktrace -llapack -lblas -L/usr/include/ -lfftw3 -fstack-arrays -lgomp -lm -ldl -pg
# LINK_FLAGS = -O3 -fopenmp -fbacktrace -llapack -lblas -lfftw3 -fstack-arrays 
#LINK_FLAGS = -O3 -fopenmp -fbacktrace -llapack -lblas -lfftw3 -fstack-arrays
# -lm -ldl -pg
#LINK_FLAGS = -O3 -fopenmp -fbacktrace /opt/software/ScaLAPACK/2.1.0-gompi-2020a/lib/libscalapack.a -lopenblas -lfftw3 -fstack-arrays
#OBJ = read-input.o generic_data.o akima_interpol.o init.o mcend-2.1.o reinit.o propa.o analysis.o
#OBJ = read-input.o generic_data.o akima_interpol.o utils.o combinations.o init.o analysis.o propa.o mcend-2.1.o

#all:	$(OBJ)
#	$(CC) -p -o mcend $(OBJ) $(LINK_FLAGS)
#        $(CLINK) -p -o mcend $(OBJ) $(LINK_FLAGS)

#ggclean:
#	rm -f *.o *.mod density.* *.Rt *.t *.ij *.chk ssqmat.diag *.log finalpsi* startpsi* finalrho *.out *.txt expec*

#%.o: %.f90
#	$(CC) $(CFLAGS) $<

#.SUFFIXES:.o.c.f90
