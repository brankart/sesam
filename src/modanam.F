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
C ---                    MODANAM.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2008-03 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE modanam
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modanam
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : anamorphosis operations
CCC  -------
CCC  Method : Identify action and perform required operations
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpmend
      use algoanam
      use algoperc
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: actionanam, flagxyo
      LOGICAL :: flaganam
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: ANAM  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C ---------------------------------------------------------------
C
      actionanam=naction
C
C -2.- Perform required action
C ----------------------------
C
      SELECT CASE (actionanam)
      CASE (1,4,7)
         flagxyo=1
      CASE (2,5,8,10)
         flagxyo=2
      CASE (3,6,9,11)
         flagxyo=3
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      IF (largtypeoper) THEN
        flaganam=(argtypeoper(1:1).EQ.'+')
      ENDIF
C
      SELECT CASE (actionanam)
      CASE (1)
C Action: -inxbas <file_xbas> -outxbasref <file_xbas>
         CALL calcperc(arginxbas,argoutxbasref,flagxyo,argconfigobs)
      CASE (2)
C Action: -inxbas <file_xybas> -outxbasref <file_ybas>
         CALL calcperc(arginybas,argoutybasref,flagxyo,argconfigobs)
      CASE (3)
C Action: -inxbas <file_xyobas> -outxbasref <file_obas> -configobs <file_o>
         CALL calcperc(arginobas,argoutobasref,flagxyo,argconfigobs)
      CASE (4)
C Action: -invar <file_x> -inxbasref <file_xbas> -outvar <file_x>
C         -typeoper <operation>
         CALL algoanamvct(arginvar,argoutvar,arginxbasref,
     $                    flagxyo,flaganam,argconfigobs)
      CASE (5)
C Action: -indta <file_xy> -inybasref <file_xybas> -outdta <file_y>
C         -typeoper <operation>
         CALL algoanamvct(argindta,argoutdta,arginybasref,
     $                    flagxyo,flaganam,argconfigobs)
      CASE (6)
C Action: -inobs <file_x> -inobasref <file_xbas> -outobs <file_o>
C          -typeoper <operation> -configobs <file_o>
         CALL algoanamvct(arginobs,argoutobs,arginobasref,
     $                    flagxyo,flaganam,argconfigobs)
      CASE (7)
C Action: -inxbas <file_xbas> -inxbasref <file_xbas> -outxbas <file_xbas>
C         -typeoper <operation>
         CALL algoanambas(arginxbas,argoutxbas,arginxbasref,
     $                    flagxyo,flaganam,argconfigobs)
      CASE (8)
C Action: -inybas <file_xybas> -inybasref <file_xybas> -outybas <file_ybas>
C         -typeoper <operation>
         CALL algoanambas(arginybas,argoutybas,arginybasref,
     $                    flagxyo,flaganam,argconfigobs)
      CASE (9)
C Action: -inobas <file_xyobas> -inobasref <file_xyobas> -outobas <file_obas>
C         -typeoper <operation> -configobs <file_o>
         CALL algoanambas(arginobas,argoutobas,arginobasref,
     $                    flagxyo,flaganam,argconfigobs)
      CASE (10)
C Action: -inobs <file_xyo> -inobas <file_xyobas> -outobas <file_obas>
         CALL algoanamobs(argindta,arginybas,argoutybas,
     $                    argconfigobs,flagxyo)
      CASE (11)
C Action: -inobs <file_xyo> -inobas <file_xyobas> -outobas <file_obas>
C         -configobs <file_o>
         CALL algoanamobs(arginobs,arginobas,argoutobas,
     $                    argconfigobs,flagxyo)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module ANAM    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'modanam','modanam')
C
      END
