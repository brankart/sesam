## SeSAM installation directory

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
make
```

To update the Makefile (if the source are modified) :

```bash
target="$HOME/bin/sesam"  # define the name of the executable
./mkmf -t Makefile.macro -p $target ../src/*.[Ffh]*
```
