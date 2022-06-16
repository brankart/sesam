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
      use mod_coord
      use mod_spacexyo , only : &
     &     jpo,jpoend,jpitpend,jpx,jpxend,jpyend,jprend,jpsmplend, &
     &     jpperc,poscoefobs,gridijkobs, &
     &     arraynx_jindxbeg,arraynx_jpindxend, &
     &     vo_idxbeg,vo_idxend
      use mod_mpitime

      use hiogrd
      use utilmkh
      use utilroa
      use utilfiles
      use utilvct
      use utilconstraint
      use ensdam_mcmc_update
      use ensdam_anatra
      use ensdam_obserror
      use ensdam_score_optimality
      use ensdam_storng
      use flowsampler_adv_constraint
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcmcmc

      BIGREAL, PUBLIC, save :: oestd_inflation=1. ! Observation error inflation factor

      LOGICAL, save :: obs_constraint=.FALSE.  ! Observation constraint
      LOGICAL, save :: obs_ensemble=.FALSE.    ! Ensemble in obs space
      LOGICAL, save :: obs_anam=.FALSE.        ! Anamorphosis in obs operator
      LOGICAL, save :: rebuild_state=.FALSE.   ! Rebuild in cost function
      LOGICAL, save :: all_ensemble=.FALSE.    ! Ensemble concatenating state and obs
      LOGICAL, save :: center_reduce=.FALSE.   ! Ensemble concatenating state and obs

      INTEGER, save :: flagxyo, jnxyo
      BIGREAL, dimension(:), allocatable, save :: obs
      BIGREAL, dimension(:), allocatable, save :: oestd
      BIGREAL, dimension(:), allocatable, save :: obseq
      BIGREAL, dimension(:,:), allocatable, save :: quantiles_ens
      BIGREAL, dimension(:), allocatable, save :: quantiles_ref
      BIGREAL, dimension(:), allocatable, save :: ens_mean
      BIGREAL, dimension(:), allocatable, save :: ens_std

      ! Storage for rebulding full state vector in cost function
      BIGREAL, dimension(:), allocatable, save :: fullstate
      BIGREAL, dimension(:), allocatable, save :: tmpstate

      ! Vector sizes
      INTEGER, save :: jpssize ! Size of state vector
      INTEGER, save :: jposize ! Size of observation vector
      INTEGER, save :: jpasize ! Size of "all" vector

      INTEGER, save :: njo=0 ! totel number of evaluations of the cost function

      ! Initial and target value of cost function
      BIGREAL, save :: cost_jini, cost_jopt, cost_jtest
      LOGICAL, save :: diag_in_j = .FALSE.
      ! Weights for MCMC schedule
      BIGREAL, save :: alpha, beta
      BIGREAL, save :: mcmc_schedule_factor=0.
      ! Number of last iterations over which weights are "averaged"
      INTEGER, save :: mcmc_schedule_iter=5
      ! Debugging option
      BIGREAL, save :: dyn_constraint_fac = 1.0
      LOGICAL, save :: debug_constraint = .FALSE.

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
      INTEGER :: allocok,jpitpsize,jprsize,jpsmpl,nbr,jidx
      INTEGER :: js,jr,jscl,jsmpl,flagcfg,flago,jpisize,jpjsize,jptsize
      LOGICAL :: lectinfo,mcmc_restart
      INTEGER :: jrbasdeb,jrbasfin,numidx
      CHARACTER(len=bgword) :: dirname, fname
!----------------------------------------------------------------------
!
      CALL kiss_load()
!
! Define segment of vector to read by current processor
      jnxyo=1+jproc
!
      jprsize=jprend
      jpsmpl=jpsmplend
      jpitpsize=jpitpend
      jposize=jpoend
      IF (lsplitobs) jposize=vo_idxend(jnxyo)
      IF ((lsplitobs).AND.(.NOT.lsplitstate)) GOTO 1000
!
      SELECT CASE (kflagxyo)
      CASE (1)
         jpssize=jpxend
         IF (lsplitstate) jpssize=jpx
      CASE (2)
         jpssize=jpyend
         IF (lsplitobs) GOTO 1004
      CASE (3)
         jpssize=jpoend
         IF (lsplitobs) jpssize=vo_idxend(jnxyo)
         IF (.NOT.PRESENT(kconfigo)) GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
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
      rebuild_state = dyn_constraint .OR. (.NOT.obs_ensemble)
      rebuild_state = rebuild_state .AND. lsplitstate
      all_ensemble = rebuild_state .AND. obs_ensemble
      center_reduce = largbias

      flago=3 ; lectinfo=.FALSE.
      IF (center_reduce.AND.(.NOT.largreducevar)) GOTO 101

      jpasize = jpssize
      IF (all_ensemble) jpasize = jpssize + jposize
!
! -0.- Read/initialize information related to the dynamical constraint
! --------------------------------------------------------------------
!
      IF (dyn_constraint) THEN
        ! Read horizontal grid
        IF (MAXVAL(varngrd(1:varend)).GT.1) GOTO 1000

        jpisize=MAXVAL(var_jpi(1:varend))
        jpjsize=MAXVAL(var_jpj(1:varend))

        allocate (longi(1:jpisize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        longi(:) = FREAL(0.0)

        allocate (latj(1:jpjsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        latj(:) = FREAL(0.0)

        CALL readgrd(1,1)

        ! Define parameters in flowsampler
        physical_units=.TRUE.
        ellipsoid_correction=.FALSE.
        spherical_delta=.FALSE.
        normalize_residual=.FALSE.
        qg_model=.TRUE.
        pv_model=.FALSE.
        adt_ref=-1.
        rossby_radius=3.d4
        dissip_rate=-0.3d-11
        dissip_rate=-1.d-21
        dissip_rate=0.

        ! Dynamical constraint uncertainty
        dyn_constraint_fac = 0.5_kr / ( dyn_constraint_std * dyn_constraint_std )
        IF (.NOT.normalize_residual) THEN
          dyn_constraint_fac = dyn_constraint_fac * ( secs_in_day * secs_in_day )
          dyn_constraint_fac = dyn_constraint_fac * ( secs_in_day * secs_in_day )
        ENDIF

        ! Initialize grid in flowsampler
        nlon = jpisize
        nlat = jpjsize
        lonmin = longi(1)
        lonmax = longi(jpisize)
        latmin = latj(1)
        latmax = latj(jpjsize)
        CALL defgrid()

        ! Check validity of configuration with constraint
        IF (varend.NE.1) GOTO 101
        IF (var_nam(1).NE.'ADT  ') GOTO 101
        IF (MAXVAL(varngrd(1:varend)).GT.1) GOTO 101
        IF (obs_constraint.AND.(.NOT.all_ensemble)) GOTO 101
        IF (.NOT.center_reduce) GOTO 101

        jptsize = MAXVAL(var_jpt(1:varend))
        jidx = MOD(jproc,jpproc/jptsize) + 1  ! index in current timestep
        IF (MOD(jpxend,jptsize).NE.0) GOTO 101
        IF (MOD(jpproc,jptsize).NE.0) GOTO 102
        IF (MOD(jpisize*jpjsize,jpproc/jptsize).NE.0) GOTO 102
        IF (jpxend.NE.jpisize*jpjsize*jptsize) THEN
          print *, 'Land mask with dynamical constraint: not yet coded'
          GOTO 101
        ENDIF

        ! Initialize constraint arrays
        CALL constraint_init(jptsize,dyn_constraint_dt)
      ENDIF

      IF (rebuild_state) THEN
        IF (dyn_constraint) THEN
          ! Allocate state vector for one timestep for temporary work
          allocate ( fullstate(1:jpxend/jptsize), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
          ! Allocate state vector block for temporary work
          allocate ( tmpstate(1:jpssize), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
        ELSE
          ! Allocate full state vector for temporary work
          allocate ( fullstate(1:jpxend), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
        ENDIF
      ENDIF

      IF (center_reduce) THEN
        ! Allocate ensemble mean and std
        allocate ( ens_mean(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        allocate ( ens_std(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001

        ! Read ensemble mean and std
        CALL readxyo (argbias,ens_mean(:),jnxyo,lectinfo,flagxyo)
        CALL readxyo (argreducevar,ens_std(:),jnxyo,lectinfo,flagxyo)

      ENDIF

      IF (dyn_constraint.AND.debug_constraint) THEN
        CALL readxyo ('debug_in.cpak',tmpstate(:),jnxyo,lectinfo,flagxyo)

        CALL rebuild_fullstate(tmpstate(1:jpssize))
        cost_jtest=eval_constraint(fullstate)
        cost_jtest=cost_jtest*dyn_constraint_fac
        print *, jproc,cost_jtest

        !For debug
        !CALL mk8vct(fullstate(:),zeta1(:,:),1,1,1,nbr,1)

        nbr = arraynx_jpindxend(1)
        tmpstate(1:nbr) = fullstate((jidx-1)*nbr+1:jidx*nbr)
        CALL writevar('debug_out.cpak',tmpstate(:),jnxyo)

        STOP 'Stopping after debug'
      ENDIF
!
! -1.- Read information related to observations
! ---------------------------------------------
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
        IF (lsplitobs) THEN
          CALL readpartcfgobs (kconfigo,jnxyo, &
     &                         kvectorms=vectorms(:))
          CALL readpartcfgobs (kconfigo,jnxyo, &
     &                         kgridijkobs=gridijkobs(:))
          CALL readpartcfgobs (kconfigo,jnxyo, &
     &                         kposcoefobs=poscoefobs(:,:))
        ELSE
          CALL readcfgobs (kconfigo,flagcfg, &
     &                     kvectorms=vectorms(:))
          CALL readcfgobs (kconfigo,flagcfg, &
     &                     kgridijkobs=gridijkobs(:))
          CALL readcfgobs (kconfigo,flagcfg, &
     &                     kposcoefobs=poscoefobs(:,:))
        ENDIF

      ENDIF
!
      IF (obs_constraint) THEN
! Read observations
        allocate ( obs(1:jposize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        obs(:)=FREAL(0.0)
!
        CALL readxyo(kinobs,obs(:), &
     &           jnxyo,lectinfo,flago,poscoefobs(:,:))

! Read observation error standard deviation
        obserror_type = obserror_type_sesam

        allocate ( oestd(1:jposize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        oestd(:)=FREAL(0.0)
!
        IF (largoestd) THEN
          IF ((validextvar(argoestd)).OR.(validextdta(argoestd)) &
     &        .OR.(validextobs(argoestd))) THEN
            CALL readxyo(argoestd,oestd(:), &
     &           jnxyo,lectinfo,flago,poscoefobs(:,:))
          ELSE
            GOTO 1000
          ENDIF
        ELSE
          IF (lsplitobs) THEN
            CALL mkpartorms (oestd(:),jnxyo)
          ELSE
            CALL mkyorms (oestd(:),flago)
          ENDIF
        ENDIF

! Observation error standard deviation inflation
        oestd(:)=oestd(:)*oestd_inflation

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
          CALL readbas(arganamorphosis,quantiles_ens(:,:),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))

          CALL readscalbas(arganamorphosis,'percref',quantiles_ref)

        ENDIF

      ENDIF
!
! -2.- Read multiple scale input ensemble
! ---------------------------------------
!
! Allocate Cxyo array
      allocate ( inens(1:jpasize,1:jprsize,1:jpscl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      inens(:,:,:) = FREAL(0.0)
!
      allocate ( upens(1:jpasize,1:jpsmpl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      upens(:,:) = FREAL(0.0)
!
      IF ((obs_ensemble).AND.(.NOT.all_ensemble)) THEN

        allocate ( inensobs(1:jposize,1:jprsize,1:jpscl), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        inensobs(:,:,:) = FREAL(0.0)
!
        allocate ( upensobs(1:jposize,1:jpsmpl), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        upensobs(:,:) = FREAL(0.0)

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
          CALL readbas(dirname,inens(1:jpssize,:,jscl),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(dirname,inens(1:jpssize,:,jscl),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo, &
     &                 kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT

        IF (obs_ensemble) THEN
          ! Read the same ensemble in observation space
          CALL fildirnam(dirname,kinobas,jscl)
          IF (all_ensemble) THEN
            CALL readbas(dirname,inens(jpssize+1:jpasize,:,jscl),jnxyo,&
     &                 jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))
          ELSE
            CALL readbas(dirname,inensobs(:,:,jscl),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))
          ENDIF
        ENDIF

      ENDDO
!
! -3.- Read initial condition for MCMC sampler (if needed)
! --------------------------------------------------------
! Check MCMC restart file
      WRITE(fname,'("./",A,"/",A)') koutbas(1:lenv(koutbas)), &
     &                                'mcmc_restart.txt'
      INQUIRE (FILE=fname,EXIST=mcmc_restart)

      IF (mcmc_restart) THEN

!       Read MCMC restart file
        numidx=10
        CALL openfile(numidx,fname)
        READ(numidx,*) mcmc_index
        READ(numidx,*) cost_jini, cost_jopt
        READ(numidx,*) mcmc_schedule, alpha, beta
        CLOSE(UNIT=numidx)

!       Read MCMC initial condition

        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING initial condition', &
     &                                  ' from the output directory'
        ENDIF

        jrbasdeb=1
        jrbasfin=jpsmpl

        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(koutbas,upens(1:jpssize,:),jnxyo,jrbasdeb,jrbasfin, &
     &                 lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(koutbas,upens(1:jpssize,:),jnxyo,jrbasdeb,jrbasfin, &
     &                 lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT

        IF (obs_ensemble) THEN
          ! Read the same ensemble in observation space
          IF (all_ensemble) THEN
            CALL readbas(koutobas,upens(jpssize+1:jpasize,:),jnxyo, &
     &                   jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))
          ELSE
            CALL readbas(koutobas,upensobs(:,:),jnxyo, &
     &                   jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))
          ENDIF
        ENDIF

      ELSE

          mcmc_index=1

      ENDIF

! Initializations at first iteration
      IF (mcmc_index.EQ.1) THEN

!       Initialize cost_jini
        diag_in_j=.TRUE.
        cost_jini = 0.
        DO jsmpl=1,jpsmpl
          IF (obs_ensemble.AND.(.NOT.all_ensemble)) THEN
            cost_jini = cost_jini + cost_jo( upensobs(:,jsmpl) )
          ELSE
            cost_jini = cost_jini + cost_jo( upens(:,jsmpl) )
          ENDIF
        ENDDO
        cost_jini = cost_jini / jpsmpl
        diag_in_j=.FALSE.
        IF (jproc.eq.0) PRINT *, 'Initial J:',cost_jini

!       Initialize cost_jopt
        cost_jopt = FREAL(jpoend) / 2._kr
        IF (dyn_constraint) THEN
          cost_jtest = FREAL((jpisize-2)*(jpjsize-2)*(jptsize-1)) / 2._kr
          IF (jproc.eq.0) PRINT *, 'Target Jobs:',cost_jopt
          IF (jproc.eq.0) PRINT *, 'Target Jdyn:',cost_jtest
          cost_jopt = cost_jopt + cost_jtest
        ENDIF
        IF (jproc.eq.0) PRINT *, 'Target J:',cost_jopt

!       Initialize mcmc_schedule
        mcmc_schedule = 0.0_kr
        beta = 1._kr / FREAL(mcmc_schedule_iter*jpsmpl)
        alpha = 1._kr - beta

      ENDIF
!     
! -4.- MCMC iteration
! -------------------
#if defined MPI
      call MPI_TIMER(1)
      call MPI_TIMER(0)
#endif
      IF ((obs_ensemble).AND.(.NOT.all_ensemble)) THEN
        CALL mcmc_iteration( maxiter, upensobs, inensobs, &
     &                       scl_mult(1:jpscl), cost_jo, &
     &                       upxens=upens, xens=inens, &
     &                       my_test=convergence_test )
      ELSE
        CALL mcmc_iteration( maxiter, upens, inens, &
     &                       scl_mult(1:jpscl), cost_jo, &
     &                       my_test=convergence_test )
      ENDIF
#if defined MPI
      call MPI_TIMER(1)
      call MPI_TIMER(0)
#endif
      IF (jproc.eq.0) PRINT *, 'Evaluations of cost function:', njo
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
        CALL writebas(koutbas,upens(1:jpssize,:), &
     &             jnxyo,jrbasdeb,jrbasfin)
      CASE (2)
        CALL writeyobas(koutbas,upens(1:jpssize,:), &
     &             jrbasdeb,jrbasfin,kflagxyo)
      CASE (3)
        IF (lsplitobs) THEN
          CALL writepartobas(koutbas,upens(1:jpssize,:),jnxyo, &
     &             jrbasdeb,jrbasfin,vectorms(:), &
     &             gridijkobs(:),poscoefobs(:,:))
        ELSE
          CALL writeyobas(koutbas,upens(1:jpssize,:), &
     &             jrbasdeb,jrbasfin,kflagxyo,vectorms(:), &
     &             gridijkobs(:),poscoefobs(:,:))
        ENDIF
      CASE DEFAULT
        GOTO 1000
      END SELECT

      IF ((obs_ensemble).AND.(.NOT.mcmc_restart)) THEN
        ! write the same ensemble in observation space
        IF (all_ensemble) THEN
          IF (lsplitobs) THEN
            CALL writepartobas(koutobas,upens(jpssize+1:jpasize,:),jnxyo, &
     &               jrbasdeb,jrbasfin,vectorms(:), &
     &               gridijkobs(:),poscoefobs(:,:))
          ELSE
            CALL writeyobas(koutobas,upens(jpssize+1:jpasize,:), &
     &               jrbasdeb,jrbasfin,flago,vectorms(:), &
     &               gridijkobs(:),poscoefobs(:,:))
          ENDIF
        ELSE
          IF (lsplitobs) THEN
            CALL writepartobas(koutobas,upensobs(:,:),jnxyo, &
     &               jrbasdeb,jrbasfin,vectorms(:), &
     &               gridijkobs(:),poscoefobs(:,:))
          ELSE
            CALL writeyobas(koutobas,upensobs(:,:), &
     &               jrbasdeb,jrbasfin,flago,vectorms(:), &
     &               gridijkobs(:),poscoefobs(:,:))
          ENDIF
        ENDIF
      ENDIF

!     Write restart index in MCMC chain
      WRITE(fname,'("./",A,"/",A)') koutbas(1:lenv(koutbas)), &
     &                              'mcmc_restart.txt'
      numidx=10
      CALL openfile(numidx,fname,kstatus=clunk)
      WRITE(numidx,*) mcmc_index
      WRITE(numidx,*) cost_jini, cost_jopt
      WRITE(numidx,*) mcmc_schedule, alpha, beta
      CLOSE(UNIT=numidx)
!
! --- deallocation
      IF (allocated(inens)) deallocate(inens)
      IF (allocated(upens)) deallocate(upens)
      IF (allocated(inensobs)) deallocate(inensobs)
      IF (allocated(upensobs)) deallocate(upensobs)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(obs)) deallocate(obs)
      IF (allocated(fullstate)) deallocate(fullstate)
      IF (allocated(ens_mean)) deallocate(ens_mean)
      IF (allocated(ens_std)) deallocate(ens_std)
      IF (allocated(longi)) deallocate(longi)
      IF (allocated(latj)) deallocate(latj)
!
      CALL kiss_save()
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algomcmc','algomcmc')
 1001 CALL printerror2(0,1001,3,'algomcmc','algomcmc')
 1004 CALL printerror2(0,1004,3,'algomcmc','algomcmc')
!
!
 101  WRITE (texterror,*) 'Bad dynamical constraint configuration'
      CALL printerror2(0,101,3,'algomcmc','algomcmc', &
     &     comment=texterror)
 102  WRITE (texterror,*) 'Incompatible number of processors'
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
      REAL(kind=8) :: cost_jobs, cost_jdyn, cost_jdis, cost_ratio

      cost_jobs = 0. ; cost_jdyn = 0.

! Rebuild full state vector if needed
      IF (rebuild_state) THEN
        CALL rebuild_fullstate(state(1:jpssize))
      ENDIF

! Apply dynamical constraint
      IF (dyn_constraint) THEN
        IF (dissip_rate.NE.0.) THEN
          cost_jdyn=eval_constraint(fullstate,cost_jdis)
          cost_jdis=cost_jdis*dyn_constraint_fac
        ELSE
          cost_jdyn=eval_constraint(fullstate)
        ENDIF
        cost_jdyn=cost_jdyn*dyn_constraint_fac
      ENDIF

      IF (obs_constraint) THEN
! Compute observation equivalent from input state
        IF (obs_ensemble) THEN
          IF (all_ensemble) THEN
            obseq = state(jpssize+1:jpasize)
          ELSE
            obseq = state
          ENDIF
        ELSE
          SELECT CASE (flagxyo)
          CASE (1)
            IF (rebuild_state) THEN
              CALL mkhytoo(fullstate(tabindxtoy(:)),obseq,poscoefobs)
            ELSE
              CALL mkhytoo(state(tabindxtoy(:)),obseq,poscoefobs)
            ENDIF
          CASE (2)
            CALL mkhytoo(state,obseq,poscoefobs)
          CASE (3)
            obseq = state
          END SELECT
        ENDIF
! Perform backward anamorphosis (if requested)
        IF (obs_anam) THEN
          CALL ana_backward( obseq, quantiles_ens, quantiles_ref )
        ENDIF
! Evaluate observation cost function
        cost_jobs = obserror_logpdf( obs, obseq, oestd )
      ENDIF

#if defined MPI
      CALL mpi_allreduce (mpi_in_place, cost_jobs, 1,  &
     &     mpi_double_precision,mpi_sum,mpi_comm_world,mpi_code)
      CALL mpi_allreduce (mpi_in_place, cost_jdyn, 1,  &
     &     mpi_double_precision,mpi_sum,mpi_comm_world,mpi_code)
      IF (dissip_rate.NE.0.) THEN
        CALL mpi_allreduce (mpi_in_place, cost_jdis, 1,  &
     &       mpi_double_precision,mpi_sum,mpi_comm_world,mpi_code)
      ENDIF
#endif

! Modify MCMC schedule as a function of J
      IF (mcmc_schedule_factor.GT.0.) THEN
        cost_ratio = mcmc_schedule_factor &
     &    * ( cost_jo-cost_jopt ) / ( cost_jini-cost_jopt )
        mcmc_schedule = alpha * mcmc_schedule + beta * cost_ratio
      ENDIF

      cost_jo = cost_jobs + cost_jdyn
      IF (dissip_rate.NE.0.) cost_jo = cost_jo + cost_jdis

      IF (diag_in_j) THEN
        IF (dissip_rate.NE.0.) THEN
          IF (jproc==0) PRINT *, 'Jobs, Jdyn, Jdis:',cost_jobs,cost_jdyn,cost_jdis
        ELSE
          IF (jproc==0) PRINT *, 'Jobs, Jdyn:',cost_jobs,cost_jdyn
        ENDIF
      ENDIF

      njo = njo + 1

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
      INTEGER :: jposize,jprsize,jr,jsmpl,jpsmpl,allocok
      REAL(KIND=8) :: cost_diag
!
      convergence_test = .FALSE.

! Check cost function
      diag_in_j=.TRUE.
      cost_diag = 0.
      jpsmpl=SIZE(upens,2)
      DO jsmpl=1,jpsmpl
        cost_diag = cost_diag + cost_jo( upens(:,jsmpl) )
      ENDDO
      cost_diag = cost_diag / jpsmpl
      diag_in_j=.FALSE.
      IF (jproc.eq.0) PRINT *, 'Average cost function:',cost_diag

! Allocate observation equivalent of input ensemble
      jposize = SIZE(obs,1) ; jprsize = SIZE(upens,2)
      allocate ( ensobseq(1:jposize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001

! Compute observation equivalent of input ensemble
      IF (obs_ensemble) THEN
        IF (all_ensemble) THEN
          ensobseq(1:jposize,:) = upens(jpssize+1:jpasize,:)
        ELSE
          ensobseq(:,:) = upens(:,:)
        ENDIF
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
      oestd(:)=oestd(:)/oestd_inflation
      CALL optimality_score(score, ensobseq, obs, oestd)
      oestd(:)=oestd(:)*oestd_inflation
      IF (jproc.eq.0)  PRINT *, 'OPTIMALITY:',mcmc_index,score
      IF (mcmc_schedule_factor.GT.0.) THEN
        IF (jproc.eq.0)  PRINT *, 'SCHEDULE:',1./FREAL(mcmc_index),mcmc_schedule
      ENDIF

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
      SUBROUTINE rebuild_fullstate(state)
!---------------------------------------------------------------------
!
!  Purpose : rebuild full state vector from blocks
!
!  Method : Callnack routine to provide to the MCMC sampler
!  ------
!
!  Input : state     : block of state vector
!  -----
!---------------------------------------------------------------------
! modules
! =======
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      REAL(kind=8), dimension(:), intent(in) :: state
!----------------------------------------------------------------------
      INTEGER :: siz

#if defined MPI
      IF (dyn_constraint) THEN
        IF (center_reduce) THEN
          tmpstate(:) = state(:) * ens_std(:) + ens_mean(:)
        ENDIF
        siz = arraynx_jpindxend(1)
        CALL MPI_ALLGATHER(tmpstate,siz,MPI_DOUBLE_PRECISION, &
     &                    fullstate,siz,MPI_DOUBLE_PRECISION, &
     &                    mpi_comm_timestep,mpi_code)
      ELSE
        fullstate(arraynx_jindxbeg(jnxyo):    &
     &            arraynx_jindxbeg(jnxyo)+    &
     &            arraynx_jpindxend(jnxyo)-1) &
     &    = state(1:arraynx_jpindxend(jnxyo))
        CALL mpi_allreduce(MPI_IN_PLACE,fullstate,jpxend, &
     &       mpi_double_precision,mpi_sum,mpi_comm_world,mpi_code)
      ENDIF
#endif

      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algomcmc
