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
! ---                    MODFILT.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-05 (J.M. Brankart)                      ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- original     : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  modfilt
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modfilt
!---------------------------------------------------------------------
!
!  Purpose : Low-pass filter
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use algofilt
      use mkfiltzon
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionfilt, jorder, jarg, flagyo
      CHARACTER(len=bgword) :: inyo, configo
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: FILT  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actionfilt=naction
!
! -A- flagyo
! = 2 => operate on dta files
! = 3 => operate on obs files
      flagyo=0
      SELECT CASE (actionfilt)
      CASE (1,2)
         flagyo=2
      CASE (3)
         flagyo=2
      CASE (4)
         flagyo=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -2.- Print description of operations to standard output
! -------------------------------------------------------
!
! --- inyo,configo
      SELECT CASE (flagyo)
      CASE (2)
         inyo=argindta
         configo=' '
      CASE (3)
         inyo=arginobs
         configo=argconfigobs
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -3.- Perform required action
! ----------------------------
!
      SELECT CASE (actionfilt)
      CASE (1)
! Action: (1) -incfg *.cfg -outzon *.zon
         CALL filtzon (argincfg,argoutzon)
      CASE (2)
! Action: (2) -incfg *.cfg -indta *.dta -outzon *.zon
         CALL filtzon (argincfg,argoutzon,kargindta=argindta)
      CASE (3,4)
! Action: (3) -indta *.dta -inzon *.zon -outdta *.dta
! Action: (4) -inobs *.obs -inzon *.zon -outdta *.dta -configobs *.obs
         CALL calcfilt (argoutdta,flagyo,inyo,arginzon,configo)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'modfilt','modfilt')
!
      END
