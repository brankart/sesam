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
! ---                    MODSPCT.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2016-06 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE modspct
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modspct
!---------------------------------------------------------------------
!
!  Purpose : spectral transformation
!  -------
!  Method : projection on spherical harmonics
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpmend
      use algospct
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionspct, flagxyo
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: SPCT  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actionspct=naction
!
! -2.- Perform required action
! ----------------------------
!
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
!
      SELECT CASE (actionspct)
      CASE (1)
! Action: -invar <file_x> -outvar <file_x> -typeoper <operation>
         CALL algospctvct(arginvar,argoutvar, &
     &                    flagxyo,argtypeoper,argconfigobs)
      CASE (2)
! Action: -indta <file_xy> -outdta <file_y> -typeoper <operation>
         CALL algospctvct(argindta,argoutdta, &
     &                    flagxyo,argtypeoper,argconfigobs)
      CASE (3)
! Action: -inobs <file_xyo> -configobs <file_o> -outdta <file_y> -typeoper <operation>
         CALL algospctvct(arginobs,argoutdta, &
     &                    flagxyo,argtypeoper,argconfigobs)
      CASE (4)
! Action: -inxbas <file_xbas> -outxbas <file_xbas> -typeoper <operation>
         CALL algospctbas(arginxbas,argoutxbas, &
     &                    flagxyo,argtypeoper,argconfigobs)
      CASE (5)
! Action: -inybas <file_xybas> -outybas <file_ybas> -typeoper <operation>
         CALL algospctbas(arginybas,argoutybas, &
     &                    flagxyo,argtypeoper,argconfigobs)
      CASE (6)
! Action: -inobas <file_xyobas> -configobs <file_o> -outybas <file_ybas> -typeoper <operation>
         CALL algospctbas(arginybas,argoutybas, &
     &                    flagxyo,argtypeoper,argconfigobs)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module SPCT    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modspct','modspct')
!
      END
