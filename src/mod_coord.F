C Copyright: CNRS - Université de Grenoble
C
C Contributors : Jean-Michel Brankart, Charles-Emmanuel Testut, Laurent Parent,
C                Emmanuel Cosme, Claire Lauvernet, Frédéric Castruccio
C
C Jean-Michel.Brankart@hmg.inpg.fr
C
C This software is governed by the CeCILL license under French law and
C abiding by the rules of distribution of free software.  You can  use,
C modify and/ or redistribute the software under the terms of the CeCILL
C license as circulated by CEA, CNRS and INRIA at the following URL
C "http://www.cecill.info".
C
CCC---------------------------------------------------------------------
CCC
CCC                        MODULE MOD_COORD
CCC
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Module with SESAM grid arrays
CCC  -------
CCC
CCC  original     : 99/05 (C.E. Testut)
CCC  original     : 03/02 (J.M. Brankart)
CCC---------------------------------------------------------------------
#include "config.main.h"
C
      MODULE mod_coord
C
      use mod_main
C
      IMPLICIT NONE
C
C -1.- Arrays with SESAM grids
C ----------------------------
C longi  : X coordinate array (special case: x(i,j,k)=x(i))
C latj   : Y coordinate array (special case: y(i,j,k)=y(j))
C levk   : Z coordinate array (special case: z(i,j,k)=z(k))
C gridij : (X,Y) coordinate array (special case:
C          x(i,j,k)=x(i,j),  y(i,j,k)=y(i,j))
      BIGREAL, dimension(:), allocatable, save :: longi,latj,levk,time
      TYPE (type_gridij), dimension(:,:), allocatable, save :: gridij
C
C -2.- Arrays with SESAM grid connections
C ---------------------------------------
C jpcon   : number of pairs of connection lines
C           (2 lines of a pair will be connected)
C dir1con : direction of connection lines 
C dir2con : orthogonal direction
C signcon : side of the line to connect
C idxcon  : position of the line to connect
C mincon  : beginning of the segment of the line to connect
C maxcon  : end of the segment of the line to connect
      INTEGER, save :: jpcon
      INTEGER, dimension(:,:,:,:,:), allocatable, save :: pty
      INTEGER, dimension(:,:,:), allocatable, save :: ptycon
      INTEGER, dimension(:,:), allocatable, save :: dir1con,dir2con
      INTEGER, dimension(:,:), allocatable, save :: signcon,idxcon
      INTEGER, dimension(:,:), allocatable, save :: mincon,maxcon
C ptbubcon   : pointers indicating to which point of the grid
C              corresponds each element of a local data section
C              (including the existence of grid connection)
C ptbubmask  : mask indicating which points of a local data section
C              falls outside of the grid
C ptyconmask : mask array used to check which connection line is
C              the target and which is the source
      INTEGER, dimension(:,:), allocatable, save :: ptbubcon
      LOGICAL, dimension(:), allocatable, save :: ptbubmask
      LOGICAL, dimension(:), allocatable, save :: ptyconmask
C
      END MODULE mod_coord
CC----------------------------------------------------------------------
