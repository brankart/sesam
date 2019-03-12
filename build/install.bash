#!/bin/bash
#

target="$HOME/bin/sesam"

./mkmf -t Makefile.macro -p $target ../src/*.[Ffh]90

make

#make install

