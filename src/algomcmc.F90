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
! ---                   ALGOMCMC.F90                            ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2021-07 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE calcmcmc
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algomcmc
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcmcmc

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcmcmc (kinbas,kinbas2,koutbas,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : apply MCMC sampler
!
!  Method : Accumulate data region by region
!  ------   Compute requested score (CRPS, RCRV)
!
!  Input : kinbas     : Input ensemble
!  -----   kinbas2    : Second input ensemble (to compute Schur products)
!          kflagxyo   : Vector type (1=Vx,2=Vy,3=Vo)
!          koutbas    : Output ensemble
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : &
     &     jpoend,jpitpend,jpx,jpxend,jpyend,jprend,jpsmplend, &
     &     poscoefobs,arraynx_jpindxend
      use hioxyo
      use hiobas
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinbas,kinbas2,koutbas
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:,:), allocatable, save :: inens
      BIGREAL, dimension(:,:), allocatable, save :: inens2
      BIGREAL, dimension(:,:), allocatable, save :: upens
!
      INTEGER :: allocok,jpssize,jpitpsize,jprsize,jpsmpl
      INTEGER :: jnxyo,js,jr,jsend
      LOGICAL :: lectinfo
      INTEGER :: jrbasdeb,jrbasfin
!----------------------------------------------------------------------
!
      jprsize=jprend
      jpsmpl=jpsmplend
      jpitpsize=1
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modscor/calcmcmc :'
         WRITE(numout,*) '         MCMC sampler algorithm'
      ENDIF
!
      SELECT CASE (kflagxyo)
      CASE (1)
         jpssize=jpx
      CASE (2)
         jpssize=jpyend
      CASE (3)
         jpssize=jpoend
         jpitpsize=jpitpend
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Allocate Cxyo array
      allocate ( inens(1:jpssize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      inens(:,:) = FREAL(0.0)
!
      allocate ( inens2(1:jpssize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      inens2(:,:) = FREAL(0.0)
!
      allocate ( upens(1:jpssize,1:jpsmpl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      inens2(:,:) = FREAL(0.0)
!
! Define segment of vector to read by current processor
      jnxyo=1+jproc
      IF (kflagxyo.EQ.1) THEN
        jsend=arraynx_jpindxend(jnxyo)
      ELSE
        jsend=jpssize
      ENDIF
!
! -1.- Read input ensemble
! ------------------------
!
      IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
        WRITE(numout,*) '    ==> READING input ensemble'
      ENDIF
!
      jrbasdeb=1
      jrbasfin=jprsize
      lectinfo=.FALSE.
      SELECT CASE (kflagxyo)
      CASE (1,2)
        CALL readbas(kinbas,inens(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo)
      CASE (3)
!       CALL readbas(kinbas,inens(:,:),jnxyo,jrbasdeb,jrbasfin, &
!    &           lectinfo,kflagxyo,poscoefobs(:,:))
      CASE DEFAULT
        GOTO 1000
      END SELECT
!     
! -2.- Read second input ensmeble
! -------------------------------
!
      IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
        WRITE(numout,*) '    ==> READING second input ensemble'
      ENDIF
!
      jrbasdeb=1
      jrbasfin=jprsize
      lectinfo=.FALSE.
      SELECT CASE (kflagxyo)
      CASE (1,2)
        CALL readbas(kinbas2,inens2(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo)
      CASE (3)
!       CALL readbas(kinbas2,inens2(:,:),jnxyo,jrbasdeb,jrbasfin, &
!    &           lectinfo,kflagxyo,poscoefobs(:,:))
      CASE DEFAULT
        GOTO 1000
      END SELECT
!
! -4.- Write output ensemble
! --------------------------
!
      IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
        WRITE(numout,*) '    ==> WRITING output ensemble'
      ENDIF
!
      jrbasdeb=1
      jrbasfin=jpsmpl
      SELECT CASE (kflagxyo)
      CASE (1)
        CALL writebas(koutbas,upens(:,:), &
     &             jnxyo,jrbasdeb,jrbasfin)
      CASE (2)
        CALL writeyobas(koutbas,upens(:,:), &
     &             jrbasdeb,jrbasfin,kflagxyo)
      CASE (3)
!       CALL writeyobas(koutbas,upens(:,:), &
!    &             jrbasdeb,jrbasfin,kflagxyo,vectorms(:), &
!    &             gridijkobs(:),poscoefobs(:,:))
      CASE DEFAULT
        GOTO 1000
      END SELECT
!
! --- deallocation
      IF (allocated(inens)) deallocate(inens)
      IF (allocated(inens2)) deallocate(inens2)
      IF (allocated(upens)) deallocate(upens)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algo','algo')
 1001 CALL printerror2(0,1001,3,'algo','algo')
!
!
 101  WRITE (texterror,*) 'Bad ...'
      CALL printerror2(0,101,3,'algomcmc','algomcmc', &
     &     comment=texterror)
 102  WRITE (texterror,*) 'Bad ...'
      CALL printerror2(0,102,3,'algomcmc','algomcmc', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algomcmc
