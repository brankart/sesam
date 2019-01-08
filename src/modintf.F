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
C ---                    MODINTF.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                        ---
C --- modification : 03-02 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  modintf
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modintf 
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Interface an object from a file format to another
CCC  -------
CCC  Method : Identify action and perform required operations
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use algointf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: flaginxyoz,flagxyoz,flagbas,actionintf
      CHARACTER(len=bgword) :: inxyoz,outxyoz,configoz
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: INTF  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C----------------------------------------------------------------
C
      actionintf=naction
C
C -A- flaginxyoz
C = 1 => input var files
C = 2 => input dta files
C = 3 => input obs files
C = 4 => input zon files
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
C
C -B- flagxyoz
C = 1 => output var files
C = 2 => output dta files
C = 3 => output obs files
C = 4 => output zon files
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
C
C -C- flagbas
C = 0 => interface vector object
C = 1 => interface covariance matrix object
      SELECT CASE (actionintf)
      CASE (1,2,3,4,5,6,7,8,9)
         flagbas=0
      CASE (10,11,12,13,14,15,16,17,18)
         flagbas=1
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -2.- Define input and output file or directory names
C-----------------------------------------------------
C
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
C
      IF (inxyoz.EQ.outxyoz) GOTO 110
C
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
C
C -3.- Perform required action
C ----------------------------
C
      IF  (flagbas.EQ.0) THEN
         SELECT CASE (flagxyoz)
         CASE (1,2)
            CALL algointfxyoz (inxyoz,outxyoz,flagxyoz)
         CASE (3,4)
            CALL algointfxyoz (inxyoz,outxyoz,
     $                   flagxyoz,kconfigoz=configoz)
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ELSE
         SELECT CASE (flagxyoz)
         CASE (1,2)
            CALL algointfxyozbas (inxyoz,outxyoz,flagxyoz)
         CASE (3,4)
            CALL algointfxyozbas (inxyoz,outxyoz,
     $                   flagxyoz,kconfigoz=configoz)
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ENDIF
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module INTF    &'          
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'modintf','modintf')
C
 110  WRITE (texterror,*) 'Identical input and output names'
      CALL printerror2(0,110,3,'modintf','modintf',comment=texterror)
C
      END
