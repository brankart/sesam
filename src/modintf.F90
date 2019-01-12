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
! ---                    MODINTF.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  modintf
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modintf 
!---------------------------------------------------------------------
!
!  Purpose : Interface an object from a file format to another
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use algointf
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: flaginxyoz,flagxyoz,flagbas,actionintf
      CHARACTER(len=bgword) :: inxyoz,outxyoz,configoz
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: INTF  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
!----------------------------------------------------------------
!
      actionintf=naction
!
! -A- flaginxyoz
! = 1 => input var files
! = 2 => input dta files
! = 3 => input obs files
! = 4 => input zon files
      SELECT CASE (actionintf)
      CASE (1,2,3,4,10,11,12,13)
         flaginxyoz=1
      CASE (5,6,7,14,15,16)
         flaginxyoz=2
      CASE (8,17)
         flaginxyoz=3
      CASE (9,18)
         flaginxyoz=4
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -B- flagxyoz
! = 1 => output var files
! = 2 => output dta files
! = 3 => output obs files
! = 4 => output zon files
      SELECT CASE (actionintf)
      CASE (1,10)
         flagxyoz=1
      CASE (2,5,11,14)
         flagxyoz=2
      CASE (3,6,8,12,15,17)
         flagxyoz=3
      CASE (4,7,9,13,16,18)
         flagxyoz=4
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -C- flagbas
! = 0 => interface vector object
! = 1 => interface covariance matrix object
      SELECT CASE (actionintf)
      CASE (1,2,3,4,5,6,7,8,9)
         flagbas=0
      CASE (10,11,12,13,14,15,16,17,18)
         flagbas=1
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -2.- Define input and output file or directory names
!-----------------------------------------------------
!
      IF (flagbas.EQ.0) THEN
         SELECT CASE (flaginxyoz)
         CASE (1)
            inxyoz=arginvar
         CASE (2)
            inxyoz=argindta
         CASE (3)
            inxyoz=arginobs
         CASE (4)
            inxyoz=arginzon
         CASE DEFAULT
            GOTO 1000
         END SELECT
         SELECT CASE (flagxyoz)
         CASE (1)
            outxyoz=argoutvar
         CASE (2)
            outxyoz=argoutdta
         CASE (3)
            outxyoz=argoutobs
         CASE (4)
            outxyoz=argoutzon
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ELSE
         SELECT CASE (flaginxyoz)
         CASE (1)
            inxyoz=arginxbas
         CASE (2)
            inxyoz=arginybas
         CASE (3)
            inxyoz=arginobas
         CASE (4)
            inxyoz=arginzbas
         CASE DEFAULT
            GOTO 1000
         END SELECT
         SELECT CASE (flagxyoz)
         CASE (1)
            outxyoz=argoutxbas
         CASE (2)
            outxyoz=argoutybas
         CASE (3)
            outxyoz=argoutobas
         CASE (4)
            outxyoz=argoutzbas
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ENDIF
!
      IF (inxyoz.EQ.outxyoz) GOTO 110
!
      SELECT CASE (flagxyoz)
      CASE (1,2)
         configoz=''
      CASE (3)
         configoz=argconfigobs
      CASE (4)
         configoz=argconfigzon
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -3.- Perform required action
! ----------------------------
!
      IF  (flagbas.EQ.0) THEN
         SELECT CASE (flagxyoz)
         CASE (1,2)
            CALL algointfxyoz (inxyoz,outxyoz,flagxyoz)
         CASE (3,4)
            CALL algointfxyoz (inxyoz,outxyoz, &
     &                   flagxyoz,kconfigoz=configoz)
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ELSE
         SELECT CASE (flagxyoz)
         CASE (1,2)
            CALL algointfxyozbas (inxyoz,outxyoz,flagxyoz)
         CASE (3,4)
            CALL algointfxyozbas (inxyoz,outxyoz, &
     &                   flagxyoz,kconfigoz=configoz)
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ENDIF
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module INTF    &'          
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modintf','modintf')
!
 110  WRITE (texterror,*) 'Identical input and output names'
      CALL printerror2(0,110,3,'modintf','modintf',comment=texterror)
!
      END
