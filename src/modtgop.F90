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
! ---                    MODTGOP.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2006-09 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE modtgop
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modtgop
!---------------------------------------------------------------------
!
!  Purpose : truncated Gaussian operations
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpmend
      use algotgest
      use mktgsmpl
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actiontgop
      CHARACTER(len=bgword) :: fnamein, fnameout
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: TGOP  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actiontgop=naction
!
! -2.- Perform required action
! ----------------------------
!
      SELECT CASE (actiontgop)
      CASE (1)
! Action: -outvar <file_x> -invar <file_x> -inxbas <file_xbas> -incstr <file_xbas>
         CALL calctgest(arginxbas,argincstr,arginvar,argoutvar)
      CASE (2)
! Action: -outvar <file_x> -invar <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -incstr <file_xbas>
         CALL calctgest(arginxbas,argincstr,arginvar, &
     &                  argoutvar,argoutxbas)
      CASE (3)
! Action: -invar <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -incstr <file_xbas>
         IF (jpmend.EQ.0) THEN
           CALL gsmpl(arginxbas,arginvar,argoutxbas)
         ELSE
           CALL tgsmpl(arginxbas,argincstr,arginvar,argoutxbas)
         ENDIF
      CASE (4)
! Action: -outvar <file_x> -inxbas <file_xbas> -outxbas <file_xbas> -incstr <file_xbas>
         STOP 'This option is no more available'
!        CALL algotgpar(arginxbas,argincstr,argoutvar,argoutxbas)
      CASE (5)
! Action: -outvar <file_x> -inxbas <file_xbas> -inxbasref <file_xbas>
         STOP 'This option is no more available'
!        CALL algoKStest(arginxbas,arginxbasref,argoutvar)
      CASE (6)
! Action: -outvar <file_x> -invar <file_x> -inxbas <file_xbas> -incstr <file_xbas> -inpartvar <file_x>
         CALL calctgest(arginxbas,argincstr,arginvar, &
     &                  argoutvar,karginpartvar=arginpartvar)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module TGOP    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modtgop','modtgop')
!
      END
