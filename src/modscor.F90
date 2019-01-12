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
! ---                    MODSCOR.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2015-03 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE modscor
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modscor
!---------------------------------------------------------------------
!
!  Purpose : probabilistic scores (CRPS,...)
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpmend
      use algoscor
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionscor, flagxyo
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: SCOR  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actionscor=naction
!
! -2.- Perform required action
! ----------------------------
!
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
!
      SELECT CASE (actionscor)
      CASE (1)
! Action: -inxbas <file_xbas> -invar <file_x> -typeoper <operation>
         CALL calcscor(arginxbas,arginvar,flagxyo,argtypeoper, &
     &                 arginpartvar,argconfigobs)
      CASE (2)
! Action: -inxbas <file_xybas> -indta <file_xy> -typeoper <operation>
         CALL calcscor(arginybas,argindta,flagxyo,argtypeoper, &
     &                 arginpartvar,argconfigobs)
      CASE (3)
! Action: -inxbas <file_xyobas> -inobs <file_xyo> -configobs <file_o> -typeoper <operation>
         CALL calcscor(arginobas,arginobs,flagxyo,argtypeoper, &
     &                 arginpartvar,argconfigobs)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module SCOR    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modscor','modscor')
!
      END
