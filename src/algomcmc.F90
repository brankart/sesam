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
      use mod_mask
      use mod_spacexyo , only : &
     &     jpoend,jpitpend,jpx,jpxend,jpyend,jprend,jpsmplend, &
     &     jpperc,poscoefobs,gridijkobs,arraynx_jpindxend
      use utilmkh
      use utilroa
      use utilfiles
      use ensdam_mcmc_update
      use ensdam_anatra
      use ensdam_obserror
      use ensdam_score_optimality
      use ensdam_storng
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcmcmc

      LOGICAL, save :: obs_constraint, obs_ensemble, obs_anam
      INTEGER, save :: flagxyo
      BIGREAL, dimension(:), allocatable, save :: obs
      BIGREAL, dimension(:), allocatable, save :: oestd
      BIGREAL, dimension(:), allocatable, save :: obseq
      BIGREAL, dimension(:,:), allocatable, save :: quantiles_ens
      BIGREAL, dimension(:), allocatable, save :: quantiles_ref

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
!          kconfigo   : Observation operator
!          kinobs     : Observation vector
!          kinobas    : Input ensemble in observation space
!          koutobas   : Output ensemble in observation space
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilvalid
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
      BIGREAL, dimension(:,:,:), allocatable, save :: inensobs
      BIGREAL, dimension(:,:), allocatable, save :: upens
      BIGREAL, dimension(:,:), allocatable, save :: upensobs
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      INTEGER :: allocok,jpssize,jpitpsize,jprsize,jpsmpl,jposize
      INTEGER :: jnxyo,jnobs,js,jr,jsend,jscl,flagcfg,flago
      LOGICAL :: lectinfo,filexists
      INTEGER :: jrbasdeb,jrbasfin,numidx
      CHARACTER(len=bgword) :: dirname, fname
!----------------------------------------------------------------------
!
      CALL kiss_load()
!
      jprsize=jprend
      jpsmpl=jpsmplend
      jpitpsize=jpitpend
      jposize=jpoend
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modscor/calcmcmc :'
         WRITE(numout,*) '         MCMC sampler algorithm'
      ENDIF
!
      flagxyo = kflagxyo
      obs_constraint = PRESENT(kinobs)
      obs_ensemble = PRESENT(kinobas)
      obs_anam = larganamorphosis
!
      SELECT CASE (kflagxyo)
      CASE (1)
         jpssize=jpx
      CASE (2)
         jpssize=jpyend
      CASE (3)
         jpssize=jpoend
         IF (.NOT.PRESENT(kconfigo)) GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT

!
! -1.- Read observation related information
! -----------------------------------------
!
      IF ((kflagxyo.EQ.3).OR.obs_constraint) THEN
! Operation performed in observation space -> read observation features
! Read poscoefobs, vectorms and gridijkobs arrays
!
        allocate ( poscoefobs(1:jposize,1:jpitpsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
!
        allocate ( gridijkobs(1:jposize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
!
        allocate ( vectorms(1:jposize), stat=allocok )
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
      IF (obs_constraint) THEN
! Read observations
        allocate ( obs(1:jposize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        obs(:)=FREAL(0.0)
!
        CALL readvalobs(kinobs,obs(:))

! Read observation error standard deviation
        obserror_type = obserror_type_sesam

        allocate ( oestd(1:jposize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        oestd(:)=FREAL(0.0)
!
        flago=3 ; jnobs=1 ; lectinfo=.FALSE.
        IF (largoestd) THEN
          IF ((validextvar(argoestd)).OR.(validextdta(argoestd)) &
     &        .OR.(validextobs(argoestd))) THEN
            CALL readxyo(argoestd,oestd(:), &
     &           jnobs,lectinfo,flago,poscoefobs(:,:))
          ELSE
            GOTO 1000
          ENDIF
        ELSE
          CALL mkyorms (oestd(:),flago)
        ENDIF

! Allocate arrays for observation equivalent to state vector
        allocate ( obseq(1:jposize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        obseq(:)=FREAL(0.0)

! Read anamorphosis transformation of observation
        IF (obs_anam) THEN

          allocate ( quantiles_ens(1:jposize,1:jpperc), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
          quantiles_ens(:,:)=FREAL(0.0)

          allocate ( quantiles_ref(1:jpperc), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
          quantiles_ref(:)=FREAL(0.0)

          jrbasdeb=1
          jrbasfin=jpperc
          flago=3 ; jnobs=1 ; lectinfo=.FALSE.
          CALL readbas(arganamorphosis,quantiles_ens(:,:),jnobs, &
     &                 jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))

          CALL readscalbas(arganamorphosis,'percref',quantiles_ref)

        ENDIF

      ENDIF
!
! -2.- Read multiple scale input ensemble
! ---------------------------------------
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
      IF (obs_ensemble) THEN

        allocate ( inensobs(1:jposize,1:jprsize,1:jpscl), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        inensobs(:,:,:) = FREAL(0.0)
!
        allocate ( upensobs(1:jposize,1:jpsmpl), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        upensobs(:,:) = FREAL(0.0)

      ENDIF
! Define segment of vector to read by current processor
      jnxyo=1+jproc
      IF (kflagxyo.EQ.1) THEN
        jsend=arraynx_jpindxend(jnxyo)
      ELSE
        jsend=jpssize
      ENDIF
!
      IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
        WRITE(numout,*) '    ==> READING multiple scale input ensemble'
      ENDIF
!
      jrbasdeb=1
      jrbasfin=jprsize
      lectinfo=.FALSE.
!
! Loop on scales to read multiscale ensemble
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

        IF (obs_ensemble) THEN
          flago=3 ; jnobs=1
          ! Read the same ensemble in observation space
          CALL fildirnam(dirname,kinobas,jscl)
          CALL readbas(dirname,inensobs(:,:,jscl),jnobs, &
     &                 jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))
        ENDIF

      ENDDO
!
! -3.- Read initial condition for MCMC sampler (if needed)
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

!       Read restart index in MCMC chain
        WRITE(fname,'("./",A,"/",A)') koutbas(1:lenv(koutbas)), &
     &                                'mcmc_index.txt'
        INQUIRE (FILE=fname,EXIST=filexists)
        IF (filexists) THEN
          CALL openfile(numidx,fname)
          READ(numidx,*) mcmc_index
          CLOSE(UNIT=numidx)
        ELSE
          mcmc_index=1
        ENDIF

      ENDIF
!     
! -4.- MCMC iteration
! -------------------
!
      IF (obs_ensemble) THEN
        CALL mcmc_iteration( maxiter, upensobs, inensobs, &
     &                       scl_mult(1:jpscl), cost_jo, &
     &                       upxens=upens, xens=inens, &
     &                       my_test=convergence_test )
      ELSE
        CALL mcmc_iteration( maxiter, upens, inens, &
     &                       scl_mult(1:jpscl), cost_jo, &
     &                       my_test=convergence_test )
      ENDIF
!     
! -5.- Write output ensemble
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

!     Write restart index in MCMC chain
      WRITE(fname,'("./",A,"/",A)') koutbas(1:lenv(koutbas)), &
     &                              'mcmc_index.txt'
      CALL openfile(numidx,fname)
      WRITE(numidx,*) mcmc_index
      CLOSE(UNIT=numidx)
!
! --- deallocation
      IF (allocated(inens)) deallocate(inens)
      IF (allocated(upens)) deallocate(upens)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(obs)) deallocate(obs)
!
      CALL kiss_save()
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algomcmc','algomcmc')
 1001 CALL printerror2(0,1001,3,'algomcmc','algomcmc')
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
!  Input : state   : state vector (x or hx)
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

      cost_jo = 0.

      IF (obs_constraint) THEN
! Compute observation equivalent from input state
        IF (obs_ensemble) THEN
          obseq = state
        ELSE
          SELECT CASE (flagxyo)
          CASE (1)
            CALL mkhytoo(state(tabindxtoy(:)),obseq,poscoefobs)
          CASE (2)
            CALL mkhytoo(state,obseq,poscoefobs)
          CASE (3)
            obseq = state
          END SELECT
        ENDIF
! Perform backward anamorphosis (if requested)
        IF (obs_anam) THEN
          call ana_backward( obseq, quantiles_ens, quantiles_ref )
        ENDIF
! Evaluate observation cost function
        cost_jo = obserror_logpdf( obs, obseq, oestd )
      ENDIF

      END FUNCTION cost_jo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION convergence_test(upens,upxens)
!---------------------------------------------------------------------
!
!  Purpose : check convergence of the MCMC chain
!
!  Method : callback routine providing the current ensemble
!  ------
!
!  Input : ensemble as produced by the MCMC chain
!  -----
!---------------------------------------------------------------------
! modules
! =======
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      REAL(KIND=8), DIMENSION(:,:), intent(in) :: upens
      REAL(KIND=8), DIMENSION(:,:), intent(in), optional :: upxens
      LOGICAL :: convergence_test
!----------------------------------------------------------------------
      REAL(KIND=8) :: score
      REAL(KIND=8), dimension(:,:), allocatable :: ensobseq
      INTEGER :: jposize,jprsize,jr,allocok
!
      convergence_test = .FALSE.

! Allocate observation equivalent of input ensemble
      jposize = SIZE(obs,1) ; jprsize = SIZE(upens,2)
      allocate ( ensobseq(1:jposize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001

! Compute observation equivalent of input ensemble
      IF (obs_ensemble) THEN
        ensobseq = upens
      ELSE
        SELECT CASE (flagxyo)
        CASE (1)
          DO jr=1,jprsize
            CALL mkhytoo(upens(tabindxtoy(:),jr),ensobseq(:,jr),poscoefobs)
          ENDDO
        CASE (2)
          DO jr=1,jprsize
            CALL mkhytoo(upens(:,jr),ensobseq(:,jr),poscoefobs)
          ENDDO
        CASE (3)
          ensobseq = upens
        END SELECT
      ENDIF

! Perform backward anamorphosis (if requested)
      IF (obs_anam) THEN
        call ana_backward( ensobseq, quantiles_ens, quantiles_ref )
      ENDIF

! Compute optimality score
      !CALL optimality_score(score, ensobseq, obs, cdf_obs)
      CALL optimality_score(score, ensobseq, obs, oestd)
      IF (jproc.eq.0)  PRINT *, 'OPTIMALITY:',mcmc_index,score

! Deallocate array
      IF (allocated(ensobseq)) deallocate(ensobseq)

      RETURN
!
 1001 CALL printerror2(0,1001,3,'algomcmc','algomcmc')

      END FUNCTION convergence_test
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION cdf_obs(o,y,obs_idx)
!---------------------------------------------------------------------
!
!  Purpose : callback routine to compute observation error cdf
!
!  Method : compute the rank of observation yo in p(yo|Hx), given y=Hx
!
!---------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
      REAL(KIND=8), intent(in) :: o
      REAL(KIND=8), intent(in) :: y
      INTEGER, intent(in) :: obs_idx
      REAL(KIND=8) :: cdf_obs
!----------------------------------------------------------------------

      cdf_obs = obserror_cdf( o, y, oestd(obs_idx) )

      END FUNCTION cdf_obs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algomcmc
