
ifneq ($(CUSTOM),yes)
FC = ifort
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf
ifneq ($(NCCLIB),yes)
LIBS += -lnetcdff
endif
INC = -I $(NETCDF_ROOT)/include
FFLAGS = -O3  -xHost -fp-model precise -traceback
ifeq ($(ZEN3),yes)
FFLAGS = -O3 -axCORE-AVX2 -fp-model precise -traceback
endif
PPFLAG90 = -fpp
PPFLAG77 = -fpp
DEBUGFLAG = -check all -debug all -traceback -fpe0
endif

ifeq ($(GFORTRAN),yes)
FC = gfortran
FFLAGS = -O2 -mtune=native -march=native -I $(NETCDF_ROOT)/include
PPFLAG90 = -x f95-cpp-input
PPFLAG77 = -x f77-cpp-input
DEBUGFLAG = -g -Wall -Wextra -fbounds-check -fbacktrace
endif

ifeq ($(CRAY),yes)
FC = ftn
FFLAGS = -h noomp
PPFLAG90 = -eZ
PPFLAG77 = -eZ
DEBUGFLAG =
endif

# Testing - I/O and fpmodel
ifeq ($(TEST),yes)
FFLAGS += $(DEBUGFLAG)
endif

ifeq ($(NCCLIB),yes)
FFLAGS += -Dncclib
endif

OBJ2= findxn.o filt.o sintp16.o \
	one.o amap.o lconset.o \
	maxmin.o fill.o outcdf.o prt_pan.o \
	setxyz_m.o ccinterp.o jimcc_m.o \
	latltoij_m.o xyzinfo_m.o newmpar_m.o indices_m.o \
	parm_m.o precis_m.o ind_m.o jimco_m.o jim_utils.o nfft_m.o \
	latlong_m.o cll_m.o sigdata_m.o smooth.o fillzonal.o fillj.o

one : $(OBJ2)
	$(FC) $(FFLAGS) $(OBJ2) $(LIBS) -o one

clean:
	rm -f *.o core one *.mod

.SUFFIXES:.f90
.f.o:
	$(FC) -c $(FFLAGS) $(INC) $(PPFLAG77) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INC) $(PPFLAG95) $<
%.o : %.mod

one.o sintp16.o setxyz.o: latlong_m.o
one.o outcdf.o : cll_m.o
one.o outcdf.o : sigdata_m.o
one.o : ccinterp.o
utilities.o : utilities.f90
ccinterp.o : ccinterp.f90 setxyz_m.o xyzinfo_m.o latltoij_m.o newmpar_m.o indices_m.o
latltoij_m.o : latltoij_m.f90 xyzinfo_m.o newmpar_m.o
setxyz_m.o : setxyz_m.f90 newmpar_m.o indices_m.o parm_m.o precis_m.o ind_m.o xyzinfo_m.o jimco_m.o jimcc_m.o 
xyzinfo_m.o : xyzinfo_m.f90 precis_m.o
newmpar_m.o : newmpar_m.f90 
precis_m.o : precis_m.f90
indices_m.o : indices_m.f90
parm_m.o : parm_m.f90 precis_m.o 
ind_m.o : ind_m.f90 newmpar_m.o 
jimcc_m.o : jimcc_m.f90 parm_m.o precis_m.o 
jimco_m.o : jimco_m.f90 precis_m.o jim_utils.o nfft_m.o 
jim_utils.o : jim_utils.f90 precis_m.o 
nfft_m.o : nfft_m.f90 precis_m.o 
outcdf.o : xyzinfo_m.o
