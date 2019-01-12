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
! ---                    MODOPER.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 00-02 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  modoper
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modoper
!---------------------------------------------------------------------
!
!  Purpose : Arithmetic operations on vector objetcs
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use algooper
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionoper,flaganlxyoz,flagnbxyoz
      CHARACTER(len=bgword) :: outxyoz,inxyoz,inrefxyoz,incfg, &
     &     textoper,configoz, &
     &     printoutxyoz,printincfg,printinxyoz, &
     &     printinrefxyoz,printconfigoz
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: OPER  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
!----------------------------------------------------------------
!
      actionoper=naction
!
! -A- flaganlxyoz
! = 1 => operate on var files
! = 2 => operate on dta files
! = 3 => operate on obs files
! = 4 => operate on zon files
      SELECT CASE (actionoper)
      CASE (1,5,9)
         flaganlxyoz=1
      CASE (2,6,10)
         flaganlxyoz=2
      CASE (3,7,11)
         flaganlxyoz=3
      CASE (4,8,12)
         flaganlxyoz=4
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -B- flagnbxyoz
! = 0 => operate with 0 input vector
! = 1 => operate with 1 input vector
! = 2 => operate with 2 input vectors
      SELECT CASE (actionoper)
      CASE (1,2,3,4)
         flagnbxyoz=0
      CASE (5,6,7,8)
         flagnbxyoz=1
      CASE (9,10,11,12)
         flagnbxyoz=2
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -C- operation to perform
      textoper=argtypeoper
!
! -2.- Print description of operations to standard output
! -------------------------------------------------------
! --- outxyoz
      configoz=''
      SELECT CASE (flaganlxyoz)
      CASE (1)
         outxyoz=argoutvar
         printoutxyoz='Xout'
         printconfigoz=''
      CASE (2)
         outxyoz=argoutdta
         printoutxyoz='Yout'
         printconfigoz='Hyx'
      CASE (3)
         outxyoz=argoutobs
         printoutxyoz='Oout'
         configoz=argconfigobs
         printconfigoz='Hox'
      CASE (4)
         outxyoz=argoutzon
         printoutxyoz='Zout'
         configoz=argconfigzon
         printconfigoz='Hz'
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (flagnbxyoz.EQ.0) THEN
         incfg=argincfg
         printincfg='operation'
      ENDIF
!
! --- inxyoz
      IF (flagnbxyoz.GE.1) THEN
         SELECT CASE (flaganlxyoz)
         CASE (1)
            inxyoz=arginvar
            printinxyoz='Xin'
         CASE (2)
            inxyoz=argindta
            printinxyoz='Yin'
         CASE (3)
            inxyoz=arginobs
            printinxyoz='Oin'
         CASE (4)
            inxyoz=arginzon
            printinxyoz='Zin'
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ENDIF
!
! --- inrefxyoz
      IF (flagnbxyoz.EQ.2) THEN
         SELECT CASE (flaganlxyoz)
         CASE (1)
            inrefxyoz=arginvarref
            printinrefxyoz='Xinref'
         CASE (2)
            inrefxyoz=argindtaref
            printinrefxyoz='Yinref'
         CASE (3)
            inrefxyoz=arginobsref
            printinrefxyoz='Oinref'
         CASE (4)
            inrefxyoz=arginzonref
            printinrefxyoz='Zinref'
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ENDIF
!
! -3.- Perform required action
! ----------------------------
!
      SELECT CASE (flaganlxyoz)
      CASE (1,2,3)
         CALL algoopervct (flaganlxyoz,outxyoz,flagnbxyoz,incfg,inxyoz, &
     &        inrefxyoz,textoper,configoz)
      CASE (4)
         CALL algooperzon (outxyoz,flagnbxyoz,incfg,inxyoz, &
     &        inrefxyoz,textoper,configoz)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module OPER    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modoper','modoper')
!
      END
