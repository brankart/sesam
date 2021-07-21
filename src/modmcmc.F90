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
! ---                    MODMCMC.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2021-07 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE modmcmc
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE modmcmc
!---------------------------------------------------------------------
!
!  Purpose : MCMC sampler
!  -------
!  Method : Identify action and perform required operations
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpmend
      use algomcmc
      IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: actionmcmc, flagxyo
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  Running SESAM module: MCMC  &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Define operations to perform as a function of action index
! ---------------------------------------------------------------
!
      actionmcmc=naction
!
! -2.- Perform required action
! ----------------------------
!
      SELECT CASE (actionmcmc)
      CASE (1)
! Action: -inxbas <file_xbas> -outxbas <file_xbas> -iterate <iteration number>
         flagxyo=1
         CALL calcmcmc(arginxbas,argoutxbas,flagxyo)
      CASE (2)
! Action: -inybas <file_xybas> -outybas <file_ybas> -iterate <iteration number>
         flagxyo=2
         CALL calcmcmc(arginybas,argoutybas,flagxyo)
      CASE (3)
! Action: -inobas <file_xyobas> -outobas <file_obas> -iterate <iteration number> -configobs <file_o>
         flagxyo=3
         CALL calcmcmc(arginobas,argoutobas,flagxyo, &
     &                 kconfigo=argconfigobs)
      CASE (4)
! Action: -inxbas <file_xbas> -outxbas <file_xbas> -iterate <iteration number> -inobs <file_o> -configobs <file_o>
         flagxyo=1
         CALL calcmcmc(arginxbas,argoutxbas,flagxyo, &
     &                 kconfigo=argconfigobs,kinobs=arginobs)
      CASE (5)
! Action: -inybas <file_xybas> -outybas <file_ybas> -iterate <iteration number> -inobs <file_o> -configobs <file_o>
         flagxyo=2
         CALL calcmcmc(arginybas,argoutybas,flagxyo, &
     &                 kconfigo=argconfigobs,kinobs=arginobs)
      CASE (6)
! Action: -inobas <file_xyobas> -outobas <file_obas> -iterate <iteration number> -inobs <file_o> -configobs <file_o>
         flagxyo=3
         CALL calcmcmc(arginobas,argoutobas,flagxyo, &
     &                 kconfigo=argconfigobs,kinobs=arginobs)
      CASE (7)
! Action: -inxbas <file_xbas> -outxbas <file_xbas> -iterate <iteration number> -inobs <file_o> -configobs <file_o> -inobas <dir_xyobas> -outobas <dir_obas>
         flagxyo=1
         CALL calcmcmc(arginxbas,argoutxbas,flagxyo, &
     &                 kconfigo=argconfigobs,kinobs=arginobs, &
     &                 kinobas=arginobas,koutobas=argoutobas)
      CASE (8)
! Action: -inybas <file_xybas> -outybas <file_ybas> -iterate <iteration number> -inobs <file_o> -configobs <file_o> -inobas <dir_xyobas> -outobas <dir_obas>
         flagxyo=2
         CALL calcmcmc(arginybas,argoutybas,flagxyo, &
     &                 kconfigo=argconfigobs,kinobs=arginobs, &
     &                 kinobas=arginobas,koutobas=argoutobas)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '&'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&  End of SESAM module MCMC    &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'modmcmc','modmcmc')
!
      END
