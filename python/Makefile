# Makefile to generate proposed pyf files from FORTRAN modules
# ------------------------------------------------------------
# J.-M. Brankart, January 2019
# 
include ../macro/make.ada_ifort
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
sesam.so: signatures/sesam.pyf ../src/sesam.F
	-ln -s ../build/*.o . ; rm sesam.o
	$(F2PY) $(F2PYFLAGS) $(LDFLAGS) --build-dir build --f90flags="$(CPPFLAGS) $(FFLAGS)" -c signatures/sesam.pyf ../src/sesam.F *.o $(NCDF)/libnetcdff.so $(NCDF)/libnetcdf.so
	-rm *.o

