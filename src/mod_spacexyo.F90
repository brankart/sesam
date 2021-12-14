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
!                        MODULE MOD_SPACEXYO
!
!---------------------------------------------------------------------
!
!  Purpose : Module with sesam object common arrays
!  -------   in Vx, Vy, Vo, Vz spaces, in reduced space
!            and in Io space
!
!  original     : 99/05 (C.E. Testut)
!  modification : 03/02 (J.M. Brankart)
!  modification : 07/09 (F. Castruccio)
!---------------------------------------------------------------------
#include "config.main.h90"
!
      MODULE mod_spacexyo
!
      use mod_main
!
      IMPLICIT NONE
!
! -1.- Vx space (size = jpxend)
! -----------------------------
! spvalvar : default special value in (structured) var files
! jpxend   : total size of Vx vectors
! jpx      : size of current block (segment) of Vx vectors
! jpnxend  : number of blocks (segments) in Vx vectors
      BIGREAL, parameter :: spvalvar = 999999.0_kr
      INTEGER, save :: jpx,jpxend,jpnxend
!
! arraynx_jindxbeg   : indices pointing to the beginning of Vx segments in Vx vector
! arraynx_jpindxend  : sizes of Vx segments in Vx vector
      INTEGER, dimension(:), allocatable, save :: arraynx_jindxbeg, &
     &     arraynx_jpindxend
!
! -2.- Vy space (size = jpyend)
! -----------------------------
! spvaldta : default special value in (structured) dta files
! jpyend   : total size of Vy vectors
! jpy      : current size of Vy vectors
      BIGREAL, parameter :: spvaldta = 999999.0_kr
      INTEGER, save :: jpy,jpyend
!
! -3.- Vo space (size = jpoend)
! -----------------------------
! spvalobs : default special value in obs files
! jpoend   : total size of Vo vectors
! jpo      : current size of Vo vectors
! jpitpend : total number of interpolation points in observation operator
! jptp     : current number of interpolation points in observation operator
      BIGREAL, parameter :: spvalobs = 999999.0_kr
      INTEGER, save :: jpo,jpoend,jpnoend,jpitp,jpitpend
! gridijobs  : observation locations (2D case)
! gridijkobs : observation locations (3D case)
! poscoefobs : observation operator (interpolation points and coefficients)
      TYPE (type_gridij), dimension(:), allocatable :: gridijobs
      TYPE (type_gridijk), dimension(:), allocatable :: gridijkobs
      TYPE (type_poscoef), dimension(:,:), allocatable :: poscoefobs
!
! vo_idxbeg  : indices pointing to the beginning of Vo segments in Vo vector
! vo_idxend  : sizes of Vo segments in Vo vector
      INTEGER, dimension(:), allocatable, save :: vo_idxbeg, vo_idxend
!
! -4.- Vz space (size = jpz)
! --------------------------
! jpz      : total number of subsystems (=number of local data sections)
      INTEGER, save :: jpz
! pt1bubidx,             : available arrays (size = jpz x dtaend)
! pt2bubidx, pt3bubidx,  : to store index of local data section (influence
! pt4bubidx              : bubble) corresponding to each subsystem
      INTEGER, dimension (:,:), allocatable, save :: pt1bubidx
      INTEGER, dimension (:,:), allocatable, save :: pt2bubidx
      INTEGER, dimension (:,:), allocatable, save :: pt3bubidx
      INTEGER, dimension (:,:), allocatable, save :: pt4bubidx
! pt1dtalon   : available arrays (size = jpz x dtaend)
! pt1dtalat   : to store longitude, latitude, depth or time index
! pt1dtadepth : in global structured 4D variable array
! pt1dtatime  : where to paste local data section
      INTEGER, dimension (:,:), allocatable, save :: pt1dtalon, pt1dtalat
      INTEGER, dimension (:,:), allocatable, save :: pt1dtadepth, pt1dtatime
! pt1bublon   : available arrays (size = jpz x dtaend)
! pt1bublat   : to store longitude, latitude, depth or time index
! pt1bubdepth : in local data section corresponding 
! pt1bubtime  : to global index (stored in previous arrays)
      INTEGER, dimension (:,:), allocatable, save :: pt1bublon, pt1bublat
      INTEGER, dimension (:,:), allocatable, save :: pt1bubdepth, pt1bubtime
! bubblk1 : available arrays (size = zon_jpi x zon_jpj x zon_jpk
! bubblk2 :     x zon_jpt x dtaend x jpnbub) to store series of
! bubblk3 :     local data section arrays (influence bubbles,
! bubblk4 :     local error modes, ...)
      BIGREAL, dimension (:,:,:,:,:,:), allocatable, save :: bubblk1
      BIGREAL, dimension (:,:,:,:,:,:), allocatable, save :: bubblk2
      BIGREAL, dimension (:,:,:,:,:,:), allocatable, save :: bubblk3
      BIGREAL, dimension (:,:,:,:,:,:), allocatable, save :: bubblk4
! ptbub     : available array (size = zon_jpi x zon_jpj x zon_jpk
!           :     x zon_jpt x dtaend) storing pointers positionning all
!           :     variables of current local data section in Vy vector
! vectptbub : available array (size = jpnbub) to store indices
!           :     of local data sections (bubidx) corresponding to 
!           :     all subsystems currently being resolved
      INTEGER, dimension (:,:,:,:,:), allocatable, save :: ptbub
      INTEGER, dimension (:), allocatable, save :: vectptbub
!
! -5.- Vr space (size = jprend) and Vz space (size = jpz)
! -------------------------------------------------------
! jprend   : total size of Vr vectors (reduced space)
! jpperc   : number of percentiles for anamorphosis
      INTEGER, save :: jprend, jpperc
!
! -6.- Vm space (size = jpmend)
! -----------------------------
! jpmend   : total size of Vm vectors (number of constraints)
      INTEGER, save :: jpmend, jpsmplend
!
! -7.- Io space (size = jpdbsend)
! -------------------------------
! jpdbsend   : total size of Io vectors (observation data base)
! jpdb       : current size of Io vectors
      INTEGER, save :: jpdbsend,jpdbs
!
      END MODULE mod_spacexyo
!----------------------------------------------------------------------
