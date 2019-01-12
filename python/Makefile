# Makefile to generate proposed pyf files from FORTRAN modules
# ------------------------------------------------------------
# J.-M. Brankart, January 2019
# 
include ../macro/make.cal1_ifort
#
OBJ = sesam.so
#
all: $(OBJ)
	@echo "Building shared libraries: $(OBJ)"
	-mv $(OBJ) modules
#
clean:
	-rm -f $(OBJ)
	-rm -fr build
	-mkdir build
#
sesam.so: signatures/sesam.pyf ../src/sesam.F90
	-ln -s ../build/*.o . ; rm sesam.o
	$(F2PY) $(F2PYFLAGS) $(LDFLAGS) --build-dir build --f90flags="$(CPPFLAGS) $(FFLAGS)" -c signatures/sesam.pyf ../src/sesam.F90 *.o $(NCDFLIB)/libnetcdff.so $(NCDFLIB)/libnetcdf.so
	-rm *.o

