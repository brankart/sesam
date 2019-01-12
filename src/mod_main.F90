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
!                        MODULE MOD_MAIN
!
!---------------------------------------------------------------------
!
!  Purpose : Main module with basic definitions needed everywhere
!  -------   in SESAM. This module is loaded in most SESAM routines.
!
!  original     : 99/05 (C.E. Testut)
!  modification : 03-02 (J.M. Brankart)
!---------------------------------------------------------------------
#include "config.main.h90"
!
      MODULE mod_main
!
      IMPLICIT NONE
!
#include "param.control.h90"
!
! Variables for SESAM parallelisation
!
      INTEGER, save :: jpproc=1
      INTEGER, save :: jproc=0
!
#if defined MPI
      INTEGER, save :: mpi_code
      include "mpif.h90"
#endif
!
! Variable type definitions
!
      TYPE type_swiarg
          CHARACTER(len=swilg) :: swi
          CHARACTER(len=bgword) :: arg
      END TYPE type_swiarg
!
      TYPE type_gridij
          BIGREAL :: longi
          BIGREAL :: latj
      END TYPE type_gridij
!
      TYPE type_gridijk
          BIGREAL :: longi
          BIGREAL :: latj
          BIGREAL :: levk
      END TYPE type_gridijk
!
      TYPE type_grid4d
          BIGREAL :: lon
          BIGREAL :: lat
          BIGREAL :: dep
          BIGREAL :: tim
      END TYPE type_grid4d
!
      TYPE type_poscoef
          INTEGER :: pos
          BIGREAL :: coef
      END TYPE type_poscoef
!
! Interfaces for basic utility routines
!
      INTERFACE
!
         CHARACTER(len=2048) FUNCTION readnextline(knumfil, &
     &        ktextexclusion) 
         INTEGER, intent(in) :: knumfil
         CHARACTER(len=*), intent(in) :: ktextexclusion
         END FUNCTION readnextline
!
         INTEGER FUNCTION lenv(ctext) 
         CHARACTER(len=*), intent(in) :: ctext
         END FUNCTION lenv
!
         INTEGER FUNCTION indext(ctext,ctab,imax)
         INTEGER, intent(in) :: imax
         CHARACTER(len=*), intent(in) :: ctext
         CHARACTER(len=*), dimension(:), intent(in) :: ctab
         END FUNCTION indext
!
         INTEGER FUNCTION posit(cwrd,ctext)
         CHARACTER(len=*), intent(in) :: cwrd, ctext
         END FUNCTION posit
!
         INTEGER FUNCTION mkint(ctext) 
         CHARACTER(len=*), intent(in) :: ctext
         END FUNCTION mkint
!
         SUBROUTINE PRINTERROR2 (sesamhelp,sesamiost,errortype, &
     &     namefile,nameroutine,comment)
         INTEGER, intent(in) :: sesamhelp,sesamiost
         CHARACTER(len=*), intent(in) :: namefile,nameroutine
         INTEGER, intent(in) :: errortype
         CHARACTER(len=*), optional, intent(in) :: comment
         END SUBROUTINE PRINTERROR2
!
         SUBROUTINE PRINTCOMMENT (sesamiost,ctext1,ctext2,ctext3)
         INTEGER, intent(in) :: sesamiost
         CHARACTER(len=*), intent(out) :: ctext1,ctext2,ctext3
         END SUBROUTINE PRINTCOMMENT
!
         SUBROUTINE SHELLORDER (ctext)
         CHARACTER(len=*), intent(in) :: ctext
         END SUBROUTINE SHELLORDER
!
      END INTERFACE
!
      END MODULE mod_main
!----------------------------------------------------------------------
