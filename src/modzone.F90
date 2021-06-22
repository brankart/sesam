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
! ---                    MODZONE.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-11 (J.M. Brankart)                      ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE modzone
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modzone
!---------------------------------------------------------------------
!
!  Purpose : Local data section management
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpz
      use mkzontodta
      use algozone
      use algolocal
      use mkcorrzbas
      use mkpartvar
      use mkreducevar
      use mkcnttozon
      use mkpttozon
      use mkrztovar
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionzone, ios, jz
      CHARACTER(len=bgword) :: fnamein, fnameout
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: ZONE  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actionzone=naction
!
! -2.- Get index of local data section to extract (for action 1)
! --------------------------------------------------------------
!
      SELECT CASE (actionzone)
      CASE (1)
         READ(argzonindex,'(I8)',IOSTAT=ios) jz
         IF (ios.NE.0) GOTO 113
         IF ((jz.LE.0).OR.(jz.GT.jpz)) GOTO 113
      CASE (2,3,4,5,6,7,8)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -3.- Perform required action
! ----------------------------
!
      SELECT CASE (actionzone)
      CASE (1)
! Action: -inzon *.zon -outdta *.dta
         CALL zontodta(arginzon,argoutdta,jz)
      CASE (2)
! Action: -outpartvar *.var -outzon *.zon
         IF (traditional) THEN
           CALL calczone(argoutpartvar,argoutzon,argincfg)
         ELSE
           CALL calclocal(argoutpartvar,argoutzon,argincfg)
         ENDIF
      CASE (3)
! Action: -inzon *.zon -incfg *.cfg -outzbas *.zon.bas
         CALL corrzbas(arginzon,argincfg,argoutzbas)
      CASE (4)
! Action: -incfg *.cfg -outpartvar *.var
         CALL partvar(argincfg,argoutpartvar)
      CASE (5)
! Action: -incfg *.cfg -outvar *.var
         CALL reducevar(argincfg,argoutvar)
      CASE (6)
! Action: -incfg *.cfg -outzon *.zon
         CALL cnttozon(argincfg,argoutzon)
      CASE (7)
! Action: -inzon *.zon -outptzon *.zon
         CALL pttozon(arginzon,argoutptzon)
      CASE (8)
! Action: -inrz *.crz -inpartvar *.cdf -outvar *.cdf
         CALL rztovar(arginrz,arginpartvar,argoutvar,argincfg)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module ZONE    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modzone','modzone')
!
 113  WRITE (texterror,*) 'Invalid local data section index'
      CALL printerror2(0,113,3,'modzone','modzone',comment=texterror)
!
      END
