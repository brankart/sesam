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
CCC                        MODULE MOD_CONT
CCC
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Module with SESAM contour and filter definition arrays
CCC  -------
CCC
CCC  original     : 99/05 (J.M. Brankart)
CCC  modification : 03/02 (J.M. Brankart)
CCC---------------------------------------------------------------------
#include "config.main.h"
C
      MODULE mod_cont
C
      use mod_main
C
      IMPLICIT NONE
C
C -1.- Arrays with contour definition
C -----------------------------------
C jpc        : Number of horizontal contours (defining vertical cylinders)
C jppend     : Maximal number of vertices defining a polygonal contours
C jplayend   : Maximum number of horizontal slices in cylinders defined 
C              by the horizontal polygonal contours
C jpp        : Number of vertices in every polygonal contours
C jplay      : Number of horizontal slices in every cylinders
C contij     : Horizontal position of contour vertices
C contidx    : Index of element of the partition defined by the contours and slices
C contlevmin : Minimum vertical level for each slice of every cylinders
C contlevmax : Maximum vertical level for each slice of every cylinders
      INTEGER, save :: jpc, jppend, jplayend
      INTEGER, dimension(:), allocatable, save :: jpp, jplay
      TYPE (type_gridij), dimension(:,:), allocatable :: contij
      INTEGER, dimension(:,:), allocatable :: contidx
      BIGREAL, dimension(:,:), allocatable :: contlevmin
      BIGREAL, dimension(:,:), allocatable :: contlevmax
C Piecewise function definition on parallepipedic elements
C f(x,y,z) = f0 * fx(x) * fy(y) * fz(z)
C contf0 : f0 value for each slice of every cylinders
C contfx : fx function type for each slice of every cylinders
C contfy : fy function type for each slice of every cylinders
C contfz : fz function type for each slice of every cylinders
      BIGREAL, dimension(:,:), allocatable :: contf0
      INTEGER, dimension(:,:), allocatable :: contfx
      INTEGER, dimension(:,:), allocatable :: contfy
      INTEGER, dimension(:,:), allocatable :: contfz
C
C -2.- Arrays with filter definition
C ----------------------------------
C filt_typ  : type of low-pass filter to apply
C filt_ord  : order of spline filter
C filt_coef : coeffcients for spline filter
C filt_length : horizontal cutoff length scale
C filt_threshold : minimal value in filter matrix
C filt_bub : filter matrix (as a local data section)
      CHARACTER(len=word80), save :: filt_typ
      INTEGER, save :: filt_ord
      BIGREAL, save :: filt_length, filt_threshold
      BIGREAL, dimension (:), allocatable, save :: filt_coef
      BIGREAL, dimension (:,:,:,:,:), allocatable :: filt_bub
C
      END MODULE mod_cont
CC----------------------------------------------------------------------
