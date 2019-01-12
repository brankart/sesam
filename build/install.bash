#!/bin/bash
#

platform="cal1_ifort"
target="$HOME/bin/sesam"

./mkmf -t ../macro/make.$platform -p $target ../src/*.[Ffh]90

make

#make install

