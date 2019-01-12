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
!                        MODULE MOD_MASK
!
!---------------------------------------------------------------------
!
!  Purpose : module with SESAM mask arrays
!  -------
!
!  original    : 99/05 (C.E. Testut)
!---------------------------------------------------------------------
#include "config.main.h90"
! --- serie coord
      MODULE mod_mask
!
      IMPLICIT NONE
!
! mask        : SESAM mask array
! tabindxtoy  : table of pointers from Vy vector to Vx vector
!               (telling which element of Vx vector corresponds to
!               each element of Vy vector)
!
      SMALLINTGR2, dimension (:,:,:,:), allocatable, save :: mask
      INTEGER, dimension (:), allocatable, save :: tabindxtoy
!
      END MODULE mod_mask
!----------------------------------------------------------------------
