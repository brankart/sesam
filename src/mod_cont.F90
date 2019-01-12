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
!                        MODULE MOD_CONT
!
!---------------------------------------------------------------------
!
!  Purpose : Module with SESAM contour and filter definition arrays
!  -------
!
!  original     : 99/05 (J.M. Brankart)
!  modification : 03/02 (J.M. Brankart)
!---------------------------------------------------------------------
#include "config.main.h90"
!
      MODULE mod_cont
!
      use mod_main
!
      IMPLICIT NONE
!
! -1.- Arrays with contour definition
! -----------------------------------
! jpc        : Number of horizontal contours (defining vertical cylinders)
! jppend     : Maximal number of vertices defining a polygonal contours
! jplayend   : Maximum number of horizontal slices in cylinders defined 
!              by the horizontal polygonal contours
! jpp        : Number of vertices in every polygonal contours
! jplay      : Number of horizontal slices in every cylinders
! contij     : Horizontal position of contour vertices
! contidx    : Index of element of the partition defined by the contours and slices
! contlevmin : Minimum vertical level for each slice of every cylinders
! contlevmax : Maximum vertical level for each slice of every cylinders
      INTEGER, save :: jpc, jppend, jplayend
      INTEGER, dimension(:), allocatable, save :: jpp, jplay
      TYPE (type_gridij), dimension(:,:), allocatable :: contij
      INTEGER, dimension(:,:), allocatable :: contidx
      BIGREAL, dimension(:,:), allocatable :: contlevmin
      BIGREAL, dimension(:,:), allocatable :: contlevmax
! Piecewise function definition on parallepipedic elements
! f(x,y,z) = f0 * fx(x) * fy(y) * fz(z)
! contf0 : f0 value for each slice of every cylinders
! contfx : fx function type for each slice of every cylinders
! contfy : fy function type for each slice of every cylinders
! contfz : fz function type for each slice of every cylinders
      BIGREAL, dimension(:,:), allocatable :: contf0
      INTEGER, dimension(:,:), allocatable :: contfx
      INTEGER, dimension(:,:), allocatable :: contfy
      INTEGER, dimension(:,:), allocatable :: contfz
!
! -2.- Arrays with filter definition
! ----------------------------------
! filt_typ  : type of low-pass filter to apply
! filt_ord  : order of spline filter
! filt_coef : coeffcients for spline filter
! filt_length : horizontal cutoff length scale
! filt_threshold : minimal value in filter matrix
! filt_bub : filter matrix (as a local data section)
      CHARACTER(len=word80), save :: filt_typ
      INTEGER, save :: filt_ord
      BIGREAL, save :: filt_length, filt_threshold
      BIGREAL, dimension (:), allocatable, save :: filt_coef
      BIGREAL, dimension (:,:,:,:,:), allocatable :: filt_bub
!
      END MODULE mod_cont
!----------------------------------------------------------------------
