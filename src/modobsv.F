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
C ---                    MODOBSV.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-12 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                       ---
C --- modification : 03-02 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  modobsv
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modobsv
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Observation management mode
CCC  -------
CCC  Method : Identify action and perform required operations
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mkdbstoobs
      use mknullobs
      use mkobstodta
      use mkdtatoobs
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: actionobsv
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: OBSV  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C ---------------------------------------------------------------
C
      actionobsv=naction
C
C -2.- Perform required action
C ----------------------------
C
      SELECT CASE (actionobsv)
      CASE (1)
C Action: -indbs *.dbs -outobs *.obs
         CALL dbstoobs (argindbs,argoutobs,argaffectobs)
      CASE (2)
C Action: -nullobs toto -outobs *.obs
         CALL nullobs (argnullobs,argoutobs)
      CASE (3)
C Action: -inobs *.obs -outdta *.dta
         CALL obstodta (arginobs,argoutdta)
      CASE (4)
C Action: -indta *.dta -outobs *.obs
         CALL mkdtatoobswithoutobs (argindta,argoutobs)
      CASE (5)
C Action: -inobs *.obs -indta *.dta -outobs *.obs
         CALL mkdtatoobsbyobs (arginobs,argindta,argoutobs)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module OBSV    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'modobsv','modobsv')
C
      END

