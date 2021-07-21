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
      use ensdam_mcmc_update
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcmcmc

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcmcmc (kinbas,koutbas,kflagxyo, &
     &                     kconfigo,kinobs,kinobas,koutobas)
!---------------------------------------------------------------------
!
!  Purpose : apply MCMC sampler
!
!  Method : Accumulate data region by region
!  ------   Compute requested score (CRPS, RCRV)
!
!  Input : kinbas     : Input ensemble
!  -----   koutbas    : Output ensemble
!          kflagxyo   : Vector type (1=Vx,2=Vy,3=Vo)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilvalid
      use mod_spacexyo , only : &
     &     jpoend,jpitpend,jpx,jpxend,jpyend,jprend,jpsmplend, &
     &     poscoefobs,gridijkobs,arraynx_jpindxend
      use hioxyo
      use hiobas
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinbas,koutbas
      INTEGER, intent(in) :: kflagxyo
      CHARACTER(len=*), intent(in), optional :: kconfigo
      CHARACTER(len=*), intent(in), optional :: kinobs
      CHARACTER(len=*), intent(in), optional :: kinobas
      CHARACTER(len=*), intent(in), optional :: koutobas
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:,:,:), allocatable, save :: inens
      BIGREAL, dimension(:,:), allocatable, save :: upens
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      INTEGER :: allocok,jpssize,jpitpsize,jprsize,jpsmpl
      INTEGER :: jnxyo,js,jr,jsend,jscl,flagcfg
      LOGICAL :: lectinfo
      INTEGER :: jrbasdeb,jrbasfin
      CHARACTER(len=bgword) :: dirname
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
         IF (.NOT.PRESENT(kconfigo)) GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT

      IF (kflagxyo.EQ.3) THEN
! Operation performed in observation space -> read observation features
! Read poscoefobs, vectorms and gridijkobs arrays
!
        allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
!
        allocate ( gridijkobs(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
!
        allocate ( vectorms(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectorms(:) = FREAL(0.0)
!
        flagcfg=1
        CALL readcfgobs (kconfigo,flagcfg, &
     &        kvectorms=vectorms(:))
        flagcfg=2
        CALL readcfgobs (kconfigo,flagcfg, &
     &        kgridijkobs=gridijkobs(:))
        flagcfg=3
        CALL readcfgobs (kconfigo,flagcfg, &
     &        kposcoefobs=poscoefobs(:,:))

      ENDIF
!
! Allocate Cxyo array
      allocate ( inens(1:jpssize,1:jprsize,1:jpscl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      inens(:,:,:) = FREAL(0.0)
!
      allocate ( upens(1:jpssize,1:jpsmpl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      upens(:,:) = FREAL(0.0)
!
! Define segment of vector to read by current processor
      jnxyo=1+jproc
      IF (kflagxyo.EQ.1) THEN
        jsend=arraynx_jpindxend(jnxyo)
      ELSE
        jsend=jpssize
      ENDIF
!
! -1.- Read multiple scale input ensemble
! ---------------------------------------
!
      IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
        WRITE(numout,*) '    ==> READING multiple scale input ensemble'
      ENDIF
!
      jrbasdeb=1
      jrbasfin=jprsize
      lectinfo=.FALSE.
!
! Loop on scales
      DO jscl =1,jpscl

        ! Set ensemble directory name for this scale
        CALL fildirnam(dirname,kinbas,jscl)

        ! Read ensemble in this directory
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(dirname,inens(:,:,jscl),jnxyo,jrbasdeb,jrbasfin, &
     &                 lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(dirname,inens(:,:,jscl),jnxyo,jrbasdeb,jrbasfin, &
     &                 lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT

      ENDDO
!
! -2.- Read initial condition for MCMC sampler (if needed)
! --------------------------------------------------------
!
      IF (.NOT.mcmc_zero_start) THEN

        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING initial condition', &
     &                                  ' from the output directory'
        ENDIF

        jrbasdeb=1
        jrbasfin=jpsmpl

        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(koutbas,upens(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &                 lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(koutbas,upens(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &                 lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT

      ENDIF
!     
! -3.- MCMC iteration
! -------------------
!
      CALL mcmc_iteration( maxiter, upens, inens, scl_mult(1:jpscl), cost_jo )
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
        CALL writeyobas(koutbas,upens(:,:), &
     &             jrbasdeb,jrbasfin,kflagxyo,vectorms(:), &
     &             gridijkobs(:),poscoefobs(:,:))
      CASE DEFAULT
        GOTO 1000
      END SELECT
!
! --- deallocation
      IF (allocated(inens)) deallocate(inens)
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
      FUNCTION cost_jo(state)
!---------------------------------------------------------------------
!
!  Purpose : observation cost function
!
!  Method : Callnack routine to provide to the MCMC sampler
!  ------   Jo = -log p(yo|hx), computed using state vector as argument
!
!  Input : state   : state vector
!  -----
!---------------------------------------------------------------------
! modules
! =======
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      REAL(kind=8), dimension(:), intent(in) :: state
      REAL(kind=8) :: cost_jo
!----------------------------------------------------------------------
! local declarations
! ==================

      cost_jo = 0.

      END FUNCTION cost_jo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algomcmc
