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
C ---                    MODOERR.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 00-01 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                        ---
C --- modification : 03-02 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  modoerr
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modoerr
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Observation error management
CCC  -------
CCC  Method : Identify action and perform required operations
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use algooerr
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: actionoerr
      CHARACTER(len=bgword) :: textoper
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: OERR  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C ---------------------------------------------------------------
C
      actionoerr=naction
C
      SELECT CASE (actionoerr)
      CASE (1,2,5,6,7,8,9)
C --- Nothing
      CASE (3,4)
         textoper=argtypeoper
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -3.- Perform required action
C ----------------------------
C
      SELECT CASE (actionoerr)
      CASE (1)
C Action: -outerrdta *.dta -incfg *.cfg
         CALL mkerrtodta (argouterrdta,argincfg)
      CASE (2)
C Action: -inerrdta *.dta -outerrdta *.dta -incfg *.cfg
         CALL mkerrtodta (argouterrdta,argincfg,kfninerrdta=arginerrdta)
      CASE (3)
C Action: -instddta <file_xy> -indta <file_xy> -outstddta <file_y>
C         -typeoper <operation>
         CALL mkoperstddta (arginstddta,argindta,argoutstddta,textoper)
      CASE (4)
C Action: -instddta <file_xy> -inobs <file_xyo> -configobs <file_o>
C         -outstddta <file_y> -typeoper <operation> 
         CALL mkoperstddta (arginstddta,arginobs,argoutstddta,
     $                      textoper,kconfigo=argconfigobs)
      CASE (5)
C Action: -inobas *.obs.bas -outerrdta *.dta -outdtaref *.dta
         CALL mkerrobastodta (arginobas,argouterrdta)
      CASE (6)
C Action: -inobas *.obs.bas -inybas *.dta.bas -outerrdta *.dta
C         -outstddta *.dta
         CALL mkerrobastodta (arginobas,argouterrdta,
     $             kfninybas=arginybas,kfnoutstddta=argoutstddta)
      CASE (7)
C Action: -inobas *.obs.bas -inybas *.dta.bas -inybasref *.dta.bas
C         -outerrdta *.dta -outstddta *.dta
         CALL mkerrobastodta (arginobas,argouterrdta,kfninybas=arginybas,
     $        kfnoutstddta=argoutstddta,kfninybasref=arginybasref)
      CASE (8)
C Action: -inobas *.obs.bas -inybas *.dta.bas -outbiasdta *.dta
         CALL mkbiastodta (arginobas,arginybas,argoutbiasdta)
      CASE (9)
C Action: -inobas *.obs.bas -outdtaref *.dta -outstddta *.dta
         CALL mkerrobastodta (arginobas,argoutdtaref,
     $        kfnoutstddta=argoutstddta)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module OERR    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'modoerr','modoerr')
C
      END

