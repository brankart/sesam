! Copyright: CNRS - Université de Grenoble
!
! Contributors : Jean-Michel Brankart, Charles-Emmanuel Testut, Laurent Parent,
!                Emmanuel Cosme, Claire Lauvernet, Frédéric Castruccio
!
! Jean-Michel.Brankart@hmg.inpg.fr
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software.  You can  use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
!---------------------------------------------------------------------
!
!                        MODULE MOD_COORD
!
!---------------------------------------------------------------------
!
!  Purpose : Module with SESAM grid arrays
!  -------
!
!  original     : 99/05 (C.E. Testut)
!  original     : 03/02 (J.M. Brankart)
!---------------------------------------------------------------------
#include "config.main.h90"
!
      MODULE mod_coord
!
      use mod_main
!
      IMPLICIT NONE
!
! -1.- Arrays with SESAM grids
! ----------------------------
! longi  : X coordinate array (special case: x(i,j,k)=x(i))
! latj   : Y coordinate array (special case: y(i,j,k)=y(j))
! levk   : Z coordinate array (special case: z(i,j,k)=z(k))
! gridij : (X,Y) coordinate array (special case:
!          x(i,j,k)=x(i,j),  y(i,j,k)=y(i,j))
      BIGREAL, dimension(:), allocatable, save :: longi,latj,levk,time
      TYPE (type_gridij), dimension(:,:), allocatable, save :: gridij
!
! -2.- Arrays with SESAM grid connections
! ---------------------------------------
! jpcon   : number of pairs of connection lines
!           (2 lines of a pair will be connected)
! dir1con : direction of connection lines 
! dir2con : orthogonal direction
! signcon : side of the line to connect
! idxcon  : position of the line to connect
! mincon  : beginning of the segment of the line to connect
! maxcon  : end of the segment of the line to connect
      INTEGER, save :: jpcon
      INTEGER, dimension(:,:,:,:,:), allocatable, save :: pty
      INTEGER, dimension(:,:,:), allocatable, save :: ptycon
      INTEGER, dimension(:,:), allocatable, save :: dir1con,dir2con
      INTEGER, dimension(:,:), allocatable, save :: signcon,idxcon
      INTEGER, dimension(:,:), allocatable, save :: mincon,maxcon
! ptbubcon   : pointers indicating to which point of the grid
!              corresponds each element of a local data section
!              (including the existence of grid connection)
! ptbubmask  : mask indicating which points of a local data section
!              falls outside of the grid
! ptyconmask : mask array used to check which connection line is
!              the target and which is the source
      INTEGER, dimension(:,:), allocatable, save :: ptbubcon
      LOGICAL, dimension(:), allocatable, save :: ptbubmask
      LOGICAL, dimension(:), allocatable, save :: ptyconmask
!
      END MODULE mod_coord
!----------------------------------------------------------------------
