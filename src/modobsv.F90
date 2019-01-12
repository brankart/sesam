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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ---                                                           ---
! ---                    MODOBSV.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-12 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                       ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  modobsv
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modobsv
!---------------------------------------------------------------------
!
!  Purpose : Observation management mode
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mkdbstoobs
      use mknullobs
      use mkobstodta
      use mkdtatoobs
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionobsv
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: OBSV  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actionobsv=naction
!
! -2.- Perform required action
! ----------------------------
!
      SELECT CASE (actionobsv)
      CASE (1)
! Action: -indbs *.dbs -outobs *.obs
         CALL dbstoobs (argindbs,argoutobs,argaffectobs)
      CASE (2)
! Action: -nullobs toto -outobs *.obs
         CALL nullobs (argnullobs,argoutobs)
      CASE (3)
! Action: -inobs *.obs -outdta *.dta
         CALL obstodta (arginobs,argoutdta)
      CASE (4)
! Action: -indta *.dta -outobs *.obs
         CALL mkdtatoobswithoutobs (argindta,argoutobs)
      CASE (5)
! Action: -inobs *.obs -indta *.dta -outobs *.obs
         CALL mkdtatoobsbyobs (arginobs,argindta,argoutobs)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module OBSV    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modobsv','modobsv')
!
      END

