SHELL=/bin/sh

CMP = ifort
FFLAGS = -O -fpp 
INC = -I $(NETCDF_ROOT)/include
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf -lnetcdff

OBJ2= findxn.o filt.o sintp16.o \
	one.o amap.o lconset.o \
	maxmin.o fill.o outcdf.o prt_pan.o \
	setxyz_m.o ccinterp.o jimcc_m.o \
	latltoij_m.o xyzinfo_m.o newmpar_m.o indices_m.o netcdf_m.o\
	parm_m.o precis_m.o ind_m.o jimco_m.o jim_utils.o nfft_m.o \
	latlong_m.o cll_m.o sigdata_m.o smooth.o fillzonal.o fillj.o

one : $(OBJ2)
	$(CMP) $(FFLAGS) $(OBJ2) $(LIBS) -o one
#	$(CMP) $(FFLAGS) $(OBJ2) $(LIBS) -o one.C$(RES)L$(LEV)

clean:
	rm -f *.o core one *.mod

.SUFFIXES:.f90
.f.o:
	$(CMP) -c $(FFLAGS) $(INC) $<
.f90.o:
	$(CMP) -c $(FFLAGS) $(INC) $<
%.o : %.mod

one.o sintp16.o setxyz.o: latlong_m.o netcdf_m.o
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
