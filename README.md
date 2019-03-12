# SeSAM
System of Sequential Assimilation Modules

SeSAM is a collection of tools to perform the various basic operations
that are required in sequential data assimilation systems.

SeSAM is distributed under the terms of the CeCILL free software license agreement.
See LICENSE.txt for more information.

SeSAM is written in Fortran.
Interfaces to python are in development.

### Installation of SeSAM

You need a FORTRAN-90 compiler and the NetCDF library (with f90 support) installed.

You also need installing the EnsDAM library, from https://github.com/brankart/ensdam

To compile SeSAM :

- create a 'make.macro' file corresponding to your compiler in the '../macro' directory.
  This is the Makefile configurable part, which specifies
  your compiler options and where to find the NetCDF and EnsDAM libraries.

```bash
cd build
ln -sf ../macro/make.(MY_MACHINE) Makefile.macro
```

- compile with:

```bash
cd build
make
```

### SeSAM documentation

SeSAM documentation is available from:
http://pp.ige-grenoble.fr/pageperso/brankarj/SESAM/
