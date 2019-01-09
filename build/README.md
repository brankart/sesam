## SeSAM installation directory

To compile SeSAM :

- create a 'make.macro' file corresponding to your compiler in the '../macro' directory.
  This is the Makefile configurable part, which specifies
  your compiler options and where to find the NetCDF and EnsDAM libraries.

- edit the file 'Makefile.template' to include this 'make.macro' file (first line below the title)

- compile with:

```bash
make
```

To update the Makefile (if the source are modified) :

```bash
platform="linux_ifort"    # define the macro file to use
target="$HOME/bin/sesam"  # define the name of the executable
./mkmf -t macro/make.$platform -p $target ../src/*.[Ffh]
```