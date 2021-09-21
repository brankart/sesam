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
! ---                    MODRANK.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2021-09 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE modrank
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modrank
!---------------------------------------------------------------------
!
!  Purpose : computation of rank of data within an ensemble
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpmend
      use algorank
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionrank, flagxyo
      LOGICAL :: flagrank
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: RANK  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actionrank=naction
!
! -2.- Perform required action
! ----------------------------
!
      SELECT CASE (actionrank)
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
      SELECT CASE (actionrank)
      CASE (1)
! Action: -invar <file_x> -inxbasref <file_xbas> -outvar <file_x>
         CALL calcrank(arginvar,argoutvar,arginxbasref, &
     &                 flagxyo,argconfigobs)
      CASE (2)
! Action: -indta <file_xy> -inybasref <file_xybas> -outdta <file_y>
         CALL calcrank(argindta,argoutdta,arginybasref, &
     &                 flagxyo,argconfigobs)
      CASE (3)
! Action: -inobs <file_x> -inobasref <file_xbas> -outobs <file_o>
!         -configobs <file_o>
         CALL calcrank(arginobs,argoutobs,arginobasref, &
     &                 flagxyo,argconfigobs)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module RANK    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modrank','modrank')
!
      END
