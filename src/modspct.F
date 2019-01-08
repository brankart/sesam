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
C ---                    MODSPCT.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2016-06 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE modspct
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modspct
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : spectral transformation
CCC  -------
CCC  Method : projection on spherical harmonics
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpmend
      use algospct
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: actionspct, flagxyo
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: SPCT  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C ---------------------------------------------------------------
C
      actionspct=naction
C
C -2.- Perform required action
C ----------------------------
C
      SELECT CASE (actionspct)
      CASE (1,4)
         flagxyo=1
      CASE (2,5)
         flagxyo=2
      CASE (3,6)
         flagxyo=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      SELECT CASE (actionspct)
      CASE (1)
C Action: -invar <file_x> -outvar <file_x> -typeoper <operation>
         CALL algospctvct(arginvar,argoutvar,
     $                    flagxyo,argtypeoper,argconfigobs)
      CASE (2)
C Action: -indta <file_xy> -outdta <file_y> -typeoper <operation>
         CALL algospctvct(argindta,argoutdta,
     $                    flagxyo,argtypeoper,argconfigobs)
      CASE (3)
C Action: -inobs <file_xyo> -configobs <file_o> -outdta <file_y> -typeoper <operation>
         CALL algospctvct(arginobs,argoutdta,
     $                    flagxyo,argtypeoper,argconfigobs)
      CASE (4)
C Action: -inxbas <file_xbas> -outxbas <file_xbas> -typeoper <operation>
         CALL algospctbas(arginxbas,argoutxbas,
     $                    flagxyo,argtypeoper,argconfigobs)
      CASE (5)
C Action: -inybas <file_xybas> -outybas <file_ybas> -typeoper <operation>
         CALL algospctbas(arginybas,argoutybas,
     $                    flagxyo,argtypeoper,argconfigobs)
      CASE (6)
C Action: -inobas <file_xyobas> -configobs <file_o> -outybas <file_ybas> -typeoper <operation>
         CALL algospctbas(arginybas,argoutybas,
     $                    flagxyo,argtypeoper,argconfigobs)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module SPCT    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'modspct','modspct')
C
      END
