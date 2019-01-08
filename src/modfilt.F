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
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C ---                                                           ---
C ---                    MODFILT.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-05 (J.M. Brankart)                      ---
C --- modification : 01-06 (C.E. Testut)                        ---
C --- original     : 03-02 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  modfilt
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modfilt
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Low-pass filter
CCC  -------
CCC  Method : Identify action and perform required operations
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use algofilt
      use mkfiltzon
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: actionfilt, jorder, jarg, flagyo
      CHARACTER(len=bgword) :: inyo, configo
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: FILT  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C ---------------------------------------------------------------
C
      actionfilt=naction
C
C -A- flagyo
C = 2 => operate on dta files
C = 3 => operate on obs files
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
C
C -2.- Print description of operations to standard output
C -------------------------------------------------------
C
C --- inyo,configo
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
C
C -3.- Perform required action
C ----------------------------
C
      SELECT CASE (actionfilt)
      CASE (1)
C Action: (1) -incfg *.cfg -outzon *.zon
         CALL filtzon (argincfg,argoutzon)
      CASE (2)
C Action: (2) -incfg *.cfg -indta *.dta -outzon *.zon
         CALL filtzon (argincfg,argoutzon,kargindta=argindta)
      CASE (3,4)
C Action: (3) -indta *.dta -inzon *.zon -outdta *.dta
C Action: (4) -inobs *.obs -inzon *.zon -outdta *.dta -configobs *.obs
         CALL calcfilt (argoutdta,flagyo,inyo,arginzon,configo)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'modfilt','modfilt')
C
      END
