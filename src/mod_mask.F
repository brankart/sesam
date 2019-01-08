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
CCC                        MODULE MOD_MASK
CCC
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : module with SESAM mask arrays
CCC  -------
CCC
CCC  original    : 99/05 (C.E. Testut)
CCC---------------------------------------------------------------------
#include "config.main.h"
C --- serie coord
      MODULE mod_mask
C
      IMPLICIT NONE
C
C mask        : SESAM mask array
C tabindxtoy  : table of pointers from Vy vector to Vx vector
C               (telling which element of Vx vector corresponds to
C               each element of Vy vector)
C
      SMALLINTGR2, dimension (:,:,:,:), allocatable, save :: mask
      INTEGER, dimension (:), allocatable, save :: tabindxtoy
C
      END MODULE mod_mask
CC----------------------------------------------------------------------
