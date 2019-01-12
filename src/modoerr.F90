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
! ---                    MODOERR.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 00-01 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  modoerr
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modoerr
!---------------------------------------------------------------------
!
!  Purpose : Observation error management
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use algooerr
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionoerr
      CHARACTER(len=bgword) :: textoper
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: OERR  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actionoerr=naction
!
      SELECT CASE (actionoerr)
      CASE (1,2,5,6,7,8,9)
! --- Nothing
      CASE (3,4)
         textoper=argtypeoper
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -3.- Perform required action
! ----------------------------
!
      SELECT CASE (actionoerr)
      CASE (1)
! Action: -outerrdta *.dta -incfg *.cfg
         CALL mkerrtodta (argouterrdta,argincfg)
      CASE (2)
! Action: -inerrdta *.dta -outerrdta *.dta -incfg *.cfg
         CALL mkerrtodta (argouterrdta,argincfg,kfninerrdta=arginerrdta)
      CASE (3)
! Action: -instddta <file_xy> -indta <file_xy> -outstddta <file_y>
!         -typeoper <operation>
         CALL mkoperstddta (arginstddta,argindta,argoutstddta,textoper)
      CASE (4)
! Action: -instddta <file_xy> -inobs <file_xyo> -configobs <file_o>
!         -outstddta <file_y> -typeoper <operation> 
         CALL mkoperstddta (arginstddta,arginobs,argoutstddta, &
     &                      textoper,kconfigo=argconfigobs)
      CASE (5)
! Action: -inobas *.obs.bas -outerrdta *.dta -outdtaref *.dta
         CALL mkerrobastodta (arginobas,argouterrdta)
      CASE (6)
! Action: -inobas *.obs.bas -inybas *.dta.bas -outerrdta *.dta
!         -outstddta *.dta
         CALL mkerrobastodta (arginobas,argouterrdta, &
     &             kfninybas=arginybas,kfnoutstddta=argoutstddta)
      CASE (7)
! Action: -inobas *.obs.bas -inybas *.dta.bas -inybasref *.dta.bas
!         -outerrdta *.dta -outstddta *.dta
         CALL mkerrobastodta (arginobas,argouterrdta,kfninybas=arginybas, &
     &        kfnoutstddta=argoutstddta,kfninybasref=arginybasref)
      CASE (8)
! Action: -inobas *.obs.bas -inybas *.dta.bas -outbiasdta *.dta
         CALL mkbiastodta (arginobas,arginybas,argoutbiasdta)
      CASE (9)
! Action: -inobas *.obs.bas -outdtaref *.dta -outstddta *.dta
         CALL mkerrobastodta (arginobas,argoutdtaref, &
     &        kfnoutstddta=argoutstddta)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module OERR    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modoerr','modoerr')
!
      END

