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
CC---------------------------------------------------------------------
CC
CC                        MODULE MOD_MAIN
CC
CC---------------------------------------------------------------------
CC
CC  Purpose : Main module with basic definitions needed everywhere
CC  -------   in SESAM. This module is loaded in most SESAM routines.
CC
CC  original     : 99/05 (C.E. Testut)
CC  modification : 03-02 (J.M. Brankart)
CC---------------------------------------------------------------------
#include "config.main.h"
C
      MODULE mod_main
C
      IMPLICIT NONE
C
#include "param.control.h"
C
C Variables for SESAM parallelisation
C
      INTEGER, save :: jpproc=1
      INTEGER, save :: jproc=0
C
#if defined MPI
      INTEGER, save :: mpi_code
      include "mpif.h"
#endif
C
C Variable type definitions
C
      TYPE type_swiarg
          CHARACTER(len=swilg) :: swi
          CHARACTER(len=bgword) :: arg
      END TYPE type_swiarg
C
      TYPE type_gridij
          BIGREAL :: longi
          BIGREAL :: latj
      END TYPE type_gridij
C
      TYPE type_gridijk
          BIGREAL :: longi
          BIGREAL :: latj
          BIGREAL :: levk
      END TYPE type_gridijk
C
      TYPE type_grid4d
          BIGREAL :: lon
          BIGREAL :: lat
          BIGREAL :: dep
          BIGREAL :: tim
      END TYPE type_grid4d
C
      TYPE type_poscoef
          INTEGER :: pos
          BIGREAL :: coef
      END TYPE type_poscoef
C
C Interfaces for basic utility routines
C
      INTERFACE
C
         CHARACTER(len=2048) FUNCTION readnextline(knumfil,
     $        ktextexclusion) 
         INTEGER, intent(in) :: knumfil
         CHARACTER(len=*), intent(in) :: ktextexclusion
         END FUNCTION readnextline
C
         INTEGER FUNCTION lenv(ctext) 
         CHARACTER(len=*), intent(in) :: ctext
         END FUNCTION lenv
C
         INTEGER FUNCTION indext(ctext,ctab,imax)
         INTEGER, intent(in) :: imax
         CHARACTER(len=*), intent(in) :: ctext
         CHARACTER(len=*), dimension(:), intent(in) :: ctab
         END FUNCTION indext
C
         INTEGER FUNCTION posit(cwrd,ctext)
         CHARACTER(len=*), intent(in) :: cwrd, ctext
         END FUNCTION posit
C
         INTEGER FUNCTION mkint(ctext) 
         CHARACTER(len=*), intent(in) :: ctext
         END FUNCTION mkint
C
         SUBROUTINE PRINTERROR2 (sesamhelp,sesamiost,errortype,
     $     namefile,nameroutine,comment)
         INTEGER, intent(in) :: sesamhelp,sesamiost
         CHARACTER(len=*), intent(in) :: namefile,nameroutine
         INTEGER, intent(in) :: errortype
         CHARACTER(len=*), optional, intent(in) :: comment
         END SUBROUTINE PRINTERROR2
C
         SUBROUTINE PRINTCOMMENT (sesamiost,ctext1,ctext2,ctext3)
         INTEGER, intent(in) :: sesamiost
         CHARACTER(len=*), intent(out) :: ctext1,ctext2,ctext3
         END SUBROUTINE PRINTCOMMENT
C
         SUBROUTINE SHELLORDER (ctext)
         CHARACTER(len=*), intent(in) :: ctext
         END SUBROUTINE SHELLORDER
C
      END INTERFACE
C
      END MODULE mod_main
CC----------------------------------------------------------------------
