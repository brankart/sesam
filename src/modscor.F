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
C ---                    MODSCOR.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2015-03 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE modscor
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modscor
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : probabilistic scores (CRPS,...)
CCC  -------
CCC  Method : Identify action and perform required operations
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpmend
      use algoscor
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: actionscor, flagxyo
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: SCOR  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C ---------------------------------------------------------------
C
      actionscor=naction
C
C -2.- Perform required action
C ----------------------------
C
      SELECT CASE (actionscor)
      CASE (1)
         flagxyo=1
      CASE (2)
         flagxyo=2
      CASE (3)
         flagxyo=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      SELECT CASE (actionscor)
      CASE (1)
C Action: -inxbas <file_xbas> -invar <file_x> -typeoper <operation>
         CALL calcscor(arginxbas,arginvar,flagxyo,argtypeoper,
     $                 arginpartvar,argconfigobs)
      CASE (2)
C Action: -inxbas <file_xybas> -indta <file_xy> -typeoper <operation>
         CALL calcscor(arginybas,argindta,flagxyo,argtypeoper,
     $                 arginpartvar,argconfigobs)
      CASE (3)
C Action: -inxbas <file_xyobas> -inobs <file_xyo> -configobs <file_o> -typeoper <operation>
         CALL calcscor(arginobas,arginobs,flagxyo,argtypeoper,
     $                 arginpartvar,argconfigobs)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module SCOR    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'modscor','modscor')
C
      END
