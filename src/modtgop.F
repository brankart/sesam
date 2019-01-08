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
C ---                    MODTGOP.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2006-09 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE modtgop
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modtgop
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : truncated Gaussian operations
CCC  -------
CCC  Method : Identify action and perform required operations
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpmend
      use algotgest
      use mktgsmpl
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: actiontgop
      CHARACTER(len=bgword) :: fnamein, fnameout
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: TGOP  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C ---------------------------------------------------------------
C
      actiontgop=naction
C
C -2.- Perform required action
C ----------------------------
C
      SELECT CASE (actiontgop)
      CASE (1)
C Action: -outvar <file_x> -invar <file_x> -inxbas <file_xbas> -incstr <file_xbas>
         CALL calctgest(arginxbas,argincstr,arginvar,argoutvar)
      CASE (2)
C Action: -outvar <file_x> -invar <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -incstr <file_xbas>
         CALL calctgest(arginxbas,argincstr,arginvar,
     $                  argoutvar,argoutxbas)
      CASE (3)
C Action: -invar <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -incstr <file_xbas>
         IF (jpmend.EQ.0) THEN
           CALL gsmpl(arginxbas,arginvar,argoutxbas)
         ELSE
           CALL tgsmpl(arginxbas,argincstr,arginvar,argoutxbas)
         ENDIF
      CASE (4)
C Action: -outvar <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -incstr <file_xbas>
         STOP 'This option is no more available'
C        CALL algotgpar(arginxbas,argincstr,argoutvar,argoutxbas)
      CASE (5)
C Action: -outvar <file_x> -inxbas <file_xbas> -inxbasref <file_xbas>
         STOP 'This option is no more available'
C        CALL algoKStest(arginxbas,arginxbasref,argoutvar)
      CASE (6)
C Action: -outvar <file_x> -invar <file_x> -inxbas <file_xbas> -incstr <file_xbas> -inpartvar <file_x>
         CALL calctgest(arginxbas,argincstr,arginvar,
     $                  argoutvar,karginpartvar=arginpartvar)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module TGOP    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'modtgop','modtgop')
C
      END
