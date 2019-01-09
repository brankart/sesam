## Python interface to EnsDAM

### Warning: This interface is still experimental.

Only the main SeSAM routine has been interfaced yet,
with SeSAM commandline options as argument.

### To generate the Python interface (with f2py) :

- edit the 'make.macro' file corresponding to your compiler in the 'macro' directory.
This is the Makefile configurable part, which specifies options to pass to f2py.

- edit the Makefile to include this 'make.macro' file (first line below the title)

- compile with:

```bash
cd python
make
```

- if everything goes well, the SeSAM shared librares should
have been created in the 'modules' directory,
ready to be imported in the python interface module (ensdam.py).

### To import EnsDAM in python:

```python
import sesam
```

