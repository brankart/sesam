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
! ---                    MODANAM.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2008-03 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE modanam
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modanam
!---------------------------------------------------------------------
!
!  Purpose : anamorphosis operations
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpmend
      use algoanam
      use algoperc
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionanam, flagxyo
      LOGICAL :: flaganam
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: ANAM  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actionanam=naction
!
! -2.- Perform required action
! ----------------------------
!
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
!
      IF (largtypeoper) THEN
        flaganam=(argtypeoper(1:1).EQ.'+')
      ENDIF
!
      SELECT CASE (actionanam)
      CASE (1)
! Action: -inxbas <file_xbas> -outxbasref <file_xbas>
         CALL calcperc(arginxbas,argoutxbasref,flagxyo,argconfigobs)
      CASE (2)
! Action: -inxbas <file_xybas> -outxbasref <file_ybas>
         CALL calcperc(arginybas,argoutybasref,flagxyo,argconfigobs)
      CASE (3)
! Action: -inxbas <file_xyobas> -outxbasref <file_obas> -configobs <file_o>
         CALL calcperc(arginobas,argoutobasref,flagxyo,argconfigobs)
      CASE (4)
! Action: -invar <file_x> -inxbasref <file_xbas> -outvar <file_x>
!         -typeoper <operation>
         CALL algoanamvct(arginvar,argoutvar,arginxbasref, &
     &                    flagxyo,flaganam,argconfigobs)
      CASE (5)
! Action: -indta <file_xy> -inybasref <file_xybas> -outdta <file_y>
!         -typeoper <operation>
         CALL algoanamvct(argindta,argoutdta,arginybasref, &
     &                    flagxyo,flaganam,argconfigobs)
      CASE (6)
! Action: -inobs <file_x> -inobasref <file_xbas> -outobs <file_o>
!          -typeoper <operation> -configobs <file_o>
         CALL algoanamvct(arginobs,argoutobs,arginobasref, &
     &                    flagxyo,flaganam,argconfigobs)
      CASE (7)
! Action: -inxbas <file_xbas> -inxbasref <file_xbas> -outxbas <file_xbas>
!         -typeoper <operation>
         CALL algoanambas(arginxbas,argoutxbas,arginxbasref, &
     &                    flagxyo,flaganam,argconfigobs)
      CASE (8)
! Action: -inybas <file_xybas> -inybasref <file_xybas> -outybas <file_ybas>
!         -typeoper <operation>
         CALL algoanambas(arginybas,argoutybas,arginybasref, &
     &                    flagxyo,flaganam,argconfigobs)
      CASE (9)
! Action: -inobas <file_xyobas> -inobasref <file_xyobas> -outobas <file_obas>
!         -typeoper <operation> -configobs <file_o>
         CALL algoanambas(arginobas,argoutobas,arginobasref, &
     &                    flagxyo,flaganam,argconfigobs)
      CASE (10)
! Action: -inobs <file_xyo> -inobas <file_xyobas> -outobas <file_obas>
         CALL algoanamobs(argindta,arginybas,argoutybas, &
     &                    argconfigobs,flagxyo)
      CASE (11)
! Action: -inobs <file_xyo> -inobas <file_xyobas> -outobas <file_obas>
!         -configobs <file_o>
         CALL algoanamobs(arginobs,arginobas,argoutobas, &
     &                    argconfigobs,flagxyo)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module ANAM    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modanam','modanam')
!
      END
