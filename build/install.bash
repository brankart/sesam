#!/bin/bash
#

platform="cal1_ifort"
target="$HOME/bin/sesam"

./mkmf -t ../macro/make.$platform -p $target ../src/*.[Ffh]

make

#make install

