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
C ---                    MODZONE.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-11 (J.M. Brankart)                      ---
C --- modification : 01-06 (C.E. Testut)                        ---
C --- modification : 03-02 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE modzone
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modzone
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Local data section management
CCC  -------
CCC  Method : Identify action and perform required operations
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpz
      use mkzontodta
      use algozone
      use mkcorrzbas
      use mkpartvar
      use mkreducevar
      use mkcnttozon
      use mkpttozon
      use mkrztovar
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: actionzone, ios, jz
      CHARACTER(len=bgword) :: fnamein, fnameout
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: ZONE  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Define operations to perform as a function of action index
C ---------------------------------------------------------------
C
      actionzone=naction
C
C -2.- Get index of local data section to extract (for action 1)
C --------------------------------------------------------------
C
      SELECT CASE (actionzone)
      CASE (1)
         READ(argzonindex,'(I8)',IOSTAT=ios) jz
         IF (ios.NE.0) GOTO 113
         IF ((jz.LE.0).OR.(jz.GT.jpz)) GOTO 113
      CASE (2,3,4,5,6,7,8)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -3.- Perform required action
C ----------------------------
C
      SELECT CASE (actionzone)
      CASE (1)
C Action: -inzon *.zon -outdta *.dta
         CALL zontodta(arginzon,argoutdta,jz)
      CASE (2)
C Action: -outpartvar *.var -outzon *.zon
         CALL calczone(argoutpartvar,argoutzon,argincfg)
      CASE (3)
C Action: -inzon *.zon -incfg *.cfg -outzbas *.zon.bas
         CALL corrzbas(arginzon,argincfg,argoutzbas)
      CASE (4)
C Action: -incfg *.cfg -outpartvar *.var
         CALL partvar(argincfg,argoutpartvar)
      CASE (5)
C Action: -incfg *.cfg -outvar *.var
         CALL reducevar(argincfg,argoutvar)
      CASE (6)
C Action: -incfg *.cfg -outzon *.zon
         CALL cnttozon(argincfg,argoutzon)
      CASE (7)
C Action: -inzon *.zon -outptzon *.zon
         CALL pttozon(arginzon,argoutptzon)
      CASE (8)
C Action: -inrz *.crz -inpartvar *.cdf -outvar *.cdf
         CALL rztovar(arginrz,arginpartvar,argoutvar,argincfg)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module ZONE    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'modzone','modzone')
C
 113  WRITE (texterror,*) 'Invalid local data section index'
      CALL printerror2(0,113,3,'modzone','modzone',comment=texterror)
C
      END
