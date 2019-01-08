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
CCC                        MODULE MOD_SPACEXYO
CCC
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Module with sesam object common arrays
CCC  -------   in Vx, Vy, Vo, Vz spaces, in reduced space
CCC            and in Io space
CCC
CCC  original     : 99/05 (C.E. Testut)
CCC  modification : 03/02 (J.M. Brankart)
CCC  modification : 07/09 (F. Castruccio)
CCC---------------------------------------------------------------------
#include "config.main.h"
C
      MODULE mod_spacexyo
C
      use mod_main
C
      IMPLICIT NONE
C
C -1.- Vx space (size = jpxend)
C -----------------------------
C spvalvar : default special value in (structured) var files
C jpxend   : total size of Vx vectors
C jpx      : size of current block (segment) of Vx vectors
C jpnxend  : number of blocks (segments) in Vx vectors
      BIGREAL, parameter :: spvalvar = 999999.0_kr
      INTEGER, save :: jpx,jpxend,jpnxend
C
C arraynx_jindxbeg   : indices pointing to the beginning of Vx segments in Vx vector
C arraynx_jpindxend  : sizes of Vx segments in Vx vector
C arraynx_jindvarbeg : index of first variable in Vx segments
C arraynx_jindkbeg   : index of first vertical level in Vx segments
C arraynx_jindtbeg   : index of first time index in Vx segments
      INTEGER, dimension(:), allocatable, save :: arraynx_jindxbeg,
     $     arraynx_jpindxend,arraynx_jindvarbeg,
     $     arraynx_jindkbeg,arraynx_jindtbeg
C
C -2.- Vy space (size = jpyend)
C -----------------------------
C spvaldta : default special value in (structured) dta files
C jpyend   : total size of Vy vectors
C jpy      : current size of Vy vectors
      BIGREAL, parameter :: spvaldta = 999999.0_kr
      INTEGER, save :: jpy,jpyend
C
C -3.- Vo space (size = jpoend)
C -----------------------------
C spvalobs : default special value in obs files
C jpoend   : total size of Vo vectors
C jpo      : current size of Vo vectors
C jpitpend : total number of interpolation points in observation operator
C jptp     : current number of interpolation points in observation operator
      BIGREAL, parameter :: spvalobs = 999999.0_kr
      INTEGER, save :: jpo,jpoend,jpitp,jpitpend
C gridijobs  : observation locations (2D case)
C gridijkobs : observation locations (3D case)
C poscoefobs : observation operator (interpolation points and coefficients)
      TYPE (type_gridij), dimension(:), allocatable :: gridijobs
      TYPE (type_gridijk), dimension(:), allocatable :: gridijkobs
      TYPE (type_poscoef), dimension(:,:), allocatable :: poscoefobs
C
C -4.- Vz space (size = jpz)
C --------------------------
C jpz      : total number of subsystems (=number of local data sections)
      INTEGER, save :: jpz
C pt1bubidx,             : available arrays (size = jpz x dtaend)
C pt2bubidx, pt3bubidx,  : to store index of local data section (influence
C pt4bubidx              : bubble) corresponding to each subsystem
      INTEGER, dimension (:,:), allocatable, save :: pt1bubidx
      INTEGER, dimension (:,:), allocatable, save :: pt2bubidx
      INTEGER, dimension (:,:), allocatable, save :: pt3bubidx
      INTEGER, dimension (:,:), allocatable, save :: pt4bubidx
C pt1dtalon   : available arrays (size = jpz x dtaend)
C pt1dtalat   : to store longitude, latitude, depth or time index
C pt1dtadepth : in global structured 4D variable array
C pt1dtatime  : where to paste local data section
      INTEGER, dimension (:,:), allocatable, save :: pt1dtalon, pt1dtalat
      INTEGER, dimension (:,:), allocatable, save :: pt1dtadepth, pt1dtatime
C pt1bublon   : available arrays (size = jpz x dtaend)
C pt1bublat   : to store longitude, latitude, depth or time index
C pt1bubdepth : in local data section corresponding 
C pt1bubtime  : to global index (stored in previous arrays)
      INTEGER, dimension (:,:), allocatable, save :: pt1bublon, pt1bublat
      INTEGER, dimension (:,:), allocatable, save :: pt1bubdepth, pt1bubtime
C bubblk1 : available arrays (size = zon_jpi x zon_jpj x zon_jpk
C bubblk2 :     x zon_jpt x dtaend x jpnbub) to store series of
C bubblk3 :     local data section arrays (influence bubbles,
C bubblk4 :     local error modes, ...)
      BIGREAL, dimension (:,:,:,:,:,:), allocatable, save :: bubblk1
      BIGREAL, dimension (:,:,:,:,:,:), allocatable, save :: bubblk2
      BIGREAL, dimension (:,:,:,:,:,:), allocatable, save :: bubblk3
      BIGREAL, dimension (:,:,:,:,:,:), allocatable, save :: bubblk4
C ptbub     : available array (size = zon_jpi x zon_jpj x zon_jpk
C           :     x zon_jpt x dtaend) storing pointers positionning all
C           :     variables of current local data section in Vy vector
C vectptbub : available array (size = jpnbub) to store indices
C           :     of local data sections (bubidx) corresponding to 
C           :     all subsystems currently being resolved
      INTEGER, dimension (:,:,:,:,:), allocatable, save :: ptbub
      INTEGER, dimension (:), allocatable, save :: vectptbub
C
C -5.- Vr space (size = jprend) and Vz space (size = jpz)
C -------------------------------------------------------
C jprend   : total size of Vr vectors (reduced space)
C jpperc   : number of percentiles for anamorphosis
      INTEGER, save :: jprend, jpperc
C
C -6.- Vm space (size = jpmend)
C -----------------------------
C jpmend   : total size of Vm vectors (number of constraints)
      INTEGER, save :: jpmend, jpsmplend
C
C -7.- Io space (size = jpdbsend)
C -------------------------------
C jpdbsend   : total size of Io vectors (observation data base)
C jpdb       : current size of Io vectors
      INTEGER, save :: jpdbsend,jpdbs
C
      END MODULE mod_spacexyo
CC----------------------------------------------------------------------
