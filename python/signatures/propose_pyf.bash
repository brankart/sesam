#!/bin/bash
#

SRCDIR="../../src"

ffile="${SRCDIR}/sesam.F90" ; modname='sesam'
f2py --overwrite-signature -m ${modname} -h ${modname}.pyf-proposed $ffile

