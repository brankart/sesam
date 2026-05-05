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
! ---                   ALGOMCMC_ADAP.F90                       ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2024-10 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE calcmcmc_adap
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algomcmc_adap
#ifdef OPENACC
      use openacc
#endif
      use mod_main
      use mod_mask
      use mod_coord
      use mod_cfgxyo
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
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcmcmc_adap

      BIGREAL, PUBLIC, save :: mcmc_adap_ratio=0.  ! Observation error inflation factor
      LOGICAL, PUBLIC, save :: mcmc_prior=.FALSE. ! Sample prior distribution (no observation constraint)

      INTEGER, save :: flagxyo, jnxyo
      BIGREAL, dimension(:), allocatable, save :: obs
      BIGREAL, dimension(:), allocatable, save :: oestd
      BIGREAL, dimension(:,:), allocatable, save :: quantiles_ens
      BIGREAL, dimension(:), allocatable, save :: quantiles_ref

      ! Vector sizes
      INTEGER, save :: jpssize  ! Size of state vector
      INTEGER, save :: jposize  ! Size of observation vector
      INTEGER, save :: jpsasize ! Size of augmented state vector (with adapative components)
      INTEGER, save :: jpoasize ! Size of augmented state vector (with adapative components)

      INTEGER, save :: njo=0 ! totel number of evaluations of the cost function

      ! Initial and target value of cost function
      BIGREAL, save :: cost_jini, cost_jopt
      LOGICAL, save :: diag_in_j = .FALSE.

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcmcmc_adap (kinbas,koutbas,kflagxyo, &
     &                          kconfigo,kinobs,kinobas,koutobas)
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
      CHARACTER(len=1) :: ens_type ! type of ensemble to read and write
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
      CASE DEFAULT
         GOTO 1000
      END SELECT

! Define size of augmented state vector and augmented observation vector
      IF (mcmc_adap_type<0) GOTO 101
      IF (mcmc_adap_type>3) GOTO 101
      jpsasize = jpssize * ( mcmc_adap_type + 1 )
      jpoasize = jposize * ( mcmc_adap_type + 1 )
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modscor/calcmcmc_adap :'
         WRITE(numout,*) '         adaptive MCMC sampler algorithm'
         WRITE(numout,*) '         -> adaptive type:', mcmc_adap_type
      ENDIF
!
      flagxyo = kflagxyo
      flago=3 ; lectinfo=.FALSE.
!
! -1.- Read information related to observations
! ---------------------------------------------
!
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
     &                       kvectorms=vectorms(:))
        CALL readpartcfgobs (kconfigo,jnxyo, &
     &                       kgridijkobs=gridijkobs(:))
        CALL readpartcfgobs (kconfigo,jnxyo, &
     &                       kposcoefobs=poscoefobs(:,:))
      ELSE
        CALL readcfgobs (kconfigo,flagcfg, &
     &                   kvectorms=vectorms(:))
        CALL readcfgobs (kconfigo,flagcfg, &
     &                   kgridijkobs=gridijkobs(:))
        CALL readcfgobs (kconfigo,flagcfg, &
     &                   kposcoefobs=poscoefobs(:,:))
      ENDIF

!
! Read observations
      allocate ( obs(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      obs(:)=FREAL(0.0)
!
      CALL readxyo(kinobs,obs(:), &
     &             jnxyo,lectinfo,flago,poscoefobs(:,:))

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
     &         jnxyo,lectinfo,flago,poscoefobs(:,:))
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

!
! -2.- Read multiple scale input ensemble
! ---------------------------------------
!
! Allocate Cxyo array
      allocate ( inens(1:jpsasize,1:jprsize,1:jpscl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      inens(:,:,:) = FREAL(0.0)
!
      allocate ( upens(1:jpsasize,1:jpsmpl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      upens(:,:) = FREAL(0.0)
!
      allocate ( inensobs(1:jpoasize,1:jprsize,1:jpscl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      inensobs(:,:,:) = FREAL(0.0)
!
      allocate ( upensobs(1:jpoasize,1:jpsmpl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      upensobs(:,:) = FREAL(0.0)
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
        CALL readbas(dirname,inens(1:jpssize,:,jscl),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo,kflagxyo)

        ! Read the same ensemble in observation space
        CALL fildirnam(dirname,kinobas,jscl)
        CALL readbas(dirname,inensobs(1:jposize,:,jscl),jnxyo, &
     &               jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))

        ! Read adaptive components of input ensemble in this directory
        IF (mcmc_adap_type>0) THEN
          CALL fildirnam(dirname,kinbas,jscl,ktype='M')
          CALL readbas(dirname,inens(  jpssize+1:2*jpssize ,:,jscl),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo,kflagxyo)
          CALL fildirnam(dirname,kinobas,jscl,ktype='M')
          CALL readbas(dirname,inensobs( jposize+1:2*jposize ,:,jscl),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))
        ENDIF
        IF (mcmc_adap_type>1) THEN
          CALL fildirnam(dirname,kinbas,jscl,ktype='A')
          CALL readbas(dirname,inens(2*jpssize+1:3*jpssize ,:,jscl),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo,kflagxyo)
          CALL fildirnam(dirname,kinobas,jscl,ktype='A')
          CALL readbas(dirname,inensobs( 2*jposize+1:3*jposize ,:,jscl),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))
        ENDIF
        IF (mcmc_adap_type>2) THEN
          CALL fildirnam(dirname,kinbas,jscl,ktype='B')
          CALL readbas(dirname,inens(3*jpssize+1:4*jpssize ,:,jscl),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo,kflagxyo)
          CALL fildirnam(dirname,kinobas,jscl,ktype='B')
          CALL readbas(dirname,inensobs( 3*jposize+1:4*jposize ,:,jscl),jnxyo, &
     &                 jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))
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
        CLOSE(UNIT=numidx)

!       Read MCMC initial condition

        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING initial condition', &
     &                                  ' from the output directory'
        ENDIF

        jrbasdeb=1
        jrbasfin=jpsmpl

        CALL readbas(koutbas,upens(1:jpssize,:),jnxyo,jrbasdeb,jrbasfin, &
     &               lectinfo,kflagxyo)

        ! Read the same ensemble in observation space
        CALL readbas(koutobas,upensobs(1:jposize,:),jnxyo, &
     &               jrbasdeb,jrbasfin,lectinfo,flago,poscoefobs(:,:))

        IF (mcmc_adap_type>0) GOTO 1000

      ELSE

        mcmc_index=1

      ENDIF

      !$acc data copyin(upensobs,obs,oestd)

! Initializations at first iteration
      IF (mcmc_index.EQ.1) THEN

!       Initialize cost_jini
        diag_in_j=.TRUE.
        cost_jini = 0.
        DO jsmpl=1,jpsmpl
          cost_jini = cost_jini + cost_jo( upensobs(:,jsmpl) )
        ENDDO
        cost_jini = cost_jini / jpsmpl
        diag_in_j=.FALSE.
        IF (jproc.eq.0) PRINT *, 'Initial J:',cost_jini

!       Initialize cost_jopt
        IF  ( (obserror_type_sesam.EQ.'gaussian') .OR. &
       &      (obserror_type_sesam.EQ.'lognormal')  ) THEN
          cost_jopt = FREAL(jpoend) / 2._kr
          IF (jproc.eq.0) PRINT *, 'Target J:',cost_jopt
        ELSE
          cost_jopt = 0.
        ENDIF

      ENDIF
!     
! -4.- MCMC iteration
! -------------------
#if defined MPI
      call MPI_TIMER(1)
      call MPI_TIMER(0)
#endif
      CALL mcmc_iteration( maxiter, upensobs, inensobs, &
     &                     scl_mult(1:jpscl), cost_jo, &
     &                     upxens=upens, xens=inens, &
     &                     my_test=convergence_test )

      !$acc end data

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
      ! Set directory name for output ensemble
      jscl=1
      CALL fildirnam(dirname,koutbas,jscl)

      jrbasdeb=1
      jrbasfin=jpsmpl
      CALL writebas(dirname,upens(1:jpssize,:), &
     &              jnxyo,jrbasdeb,jrbasfin)

      ! Write adaptive components of input ensemble in this directory
      IF (mcmc_adap_type>0) THEN
        CALL fildirnam(dirname,koutbas,jscl,ktype='M')
        CALL writebas(dirname,upens( jpssize+1:2*jpssize ,:), &
     &                jnxyo,jrbasdeb,jrbasfin)
      ENDIF
      IF (mcmc_adap_type>1) THEN
        CALL fildirnam(dirname,koutbas,jscl,ktype='A')
        CALL writebas(dirname,upens( 2*jpssize+1:3*jpssize ,:), &
     &                jnxyo,jrbasdeb,jrbasfin)
      ENDIF
      IF (mcmc_adap_type>2) THEN
        CALL fildirnam(dirname,koutbas,jscl,ktype='B')
        CALL writebas(dirname,upens( 3*jpssize+1:4*jpssize ,:), &
     &                jnxyo,jrbasdeb,jrbasfin)
      ENDIF

!     IF (.NOT.mcmc_restart) THEN
!       ! write the same ensemble in observation space
!       CALL fildirnam(dirname,koutobas,jscl)
!       IF (lsplitobs) THEN
!         CALL writepartobas(dirname,upensobs(1:jposize,:),jnxyo, &
!    &               jrbasdeb,jrbasfin,vectorms(:), &
!    &               gridijkobs(:),poscoefobs(:,:))
!       ELSE
!         CALL writeyobas(dirname,upensobs(1:jposize,:), &
!    &               jrbasdeb,jrbasfin,flago,vectorms(:), &
!    &               gridijkobs(:),poscoefobs(:,:))
!       ENDIF
!     ENDIF

      print *, 'End of ensemble writing, proc=:',jproc

!     IF (jproc.eq.0) THEN
!       Write restart index in MCMC chain
!       jscl=1
!       CALL fildirnam(dirname,koutbas,jscl)
!       WRITE(fname,'("./",A,"/",A)') dirname(1:lenv(dirname)), &
!    &                                'mcmc_restart.txt'
!       numidx=10
!       CALL openfile(numidx,fname,kstatus=clunk)
!       WRITE(numidx,*) mcmc_index
!       WRITE(numidx,*) cost_jini, cost_jopt
!       CLOSE(UNIT=numidx)
!     ENDIF
!
!     print *, 'End of restart writing, proc=:',jproc

! --- deallocation
      IF (allocated(inens)) deallocate(inens)
      IF (allocated(upens)) deallocate(upens)
      IF (allocated(inensobs)) deallocate(inensobs)
      IF (allocated(upensobs)) deallocate(upensobs)

      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(obs)) deallocate(obs)

      CALL kiss_save()

      print *, 'End of MCMC routine, proc=:',jproc
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algomcmc_adap','calcmcmc_adap')
 1001 CALL printerror2(0,1001,3,'algomcmc_adap','calcmcmc_adap')
 1004 CALL printerror2(0,1004,3,'algomcmc_adap','calcmcmc_adap')
!
!
 101  WRITE (texterror,*) 'Bad type of adaptive algorithm'
      CALL printerror2(0,101,3,'algomcmc_adap','calcmcmc_adap', &
     &     comment=texterror)
 102  WRITE (texterror,*) 'Incompatible number of processors'
      CALL printerror2(0,102,3,'algomcmc_adap','calcmcmc_adap', &
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
      REAL(kind=8) :: cost_jobs, cost_jdis, cost_ratio
      INTEGER :: jo, jpo
      REAL(kind=8) :: cost_one, cost_alpha

      IF (mcmc_prior) THEN
        cost_jo = 0.
        RETURN
      ENDIF

      jpo = jposize
      cost_jobs  = 0.
      cost_alpha = 0.

      !IF (jproc==0) print *, 'call to jo:',mcmc_index,njo

! Evaluate observation cost function
      cost_jobs = 0.
      IF (mcmc_adap_type==0) THEN
        !$acc data present(state,obs,oestd)
        !$acc parallel loop private(cost_one) reduction(+:cost_jobs)
        DO jo = 1,jpo
          cost_one  = ( state(jo) - obs(jo) ) / oestd(jo)
          cost_jobs = cost_jobs + 0.5 * cost_one * cost_one
        ENDDO
        !$acc end parallel loop
        !$acc end data
      ELSEIF (mcmc_adap_type==1) THEN
        !$acc data present(state,obs,oestd)
        !$acc parallel loop private(cost_one) reduction(+:cost_jobs)
        DO jo = 1,jpo
          cost_one  = ( state(jo) - state(jpo+jo) - obs(jo) ) / oestd(jo)
          cost_jobs = cost_jobs + 0.5 * cost_one * cost_one
        ENDDO
        !$acc end parallel loop
        !$acc end data
      ELSEIF (mcmc_adap_type==2) THEN
        cost_alpha = 0.
        !$acc data present(state,obs,oestd)
        !$acc parallel loop private(cost_one) reduction(+:cost_jobs) reduction(+:cost_alpha)
        DO jo = 1,jpo
          cost_one   = ( state(jo) - state(jpo+jo) - obs(jo) ) / oestd(jo)
          cost_jobs  = cost_jobs + 0.5 * cost_one * cost_one
          !cost_alpha = cost_alpha + 0.5 * state(jo) * state(jo) * ( EXP(-2.0*state(2*jpo+jo)) - 1.0 ) + state(2*jpo+jo)
        ENDDO
        !$acc end parallel loop
        !$acc end data
      ELSEIF (mcmc_adap_type==3) THEN
        cost_alpha = 0.
        !$acc data present(state,obs,oestd)
        !$acc parallel loop private(cost_one) reduction(+:cost_jobs) reduction(+:cost_alpha)
        DO jo = 1,jpo
          cost_one  = ( state(jo) - state(jpo+jo) - obs(jo) ) / ( oestd(jo) * EXP(state(3*jpo+jo)) )
          cost_jobs = cost_jobs + 0.5 * cost_one * cost_one
          cost_alpha = cost_alpha + state(3*jpo+jo)
          !cost_alpha = cost_alpha + 0.5 * state(jo) * state(jo) * ( EXP(-2.0*state(2*jpo+jo)) - 1.0 ) + state(2*jpo+jo)
        ENDDO
        !$acc end parallel loop
        !$acc end data
      ENDIF

      !cost_jobs = 0.

#if defined MPI
      CALL mpi_allreduce (mpi_in_place, cost_jobs, 1,  &
     &     mpi_double_precision,mpi_sum,mpi_comm_world,mpi_code)
      CALL mpi_allreduce (mpi_in_place, cost_alpha, 1,  &
     &     mpi_double_precision,mpi_sum,mpi_comm_world,mpi_code)
#endif

      IF (diag_in_j) THEN
        IF (jproc==0) PRINT *, 'Jobs:',cost_jobs, cost_alpha
      ENDIF

      !cost_jo = cost_jobs + mcmc_adap_ratio * cost_alpha
      cost_jo = cost_jobs + cost_alpha

      !IF (jproc==0) print *, 'end call to jo:',mcmc_index,njo

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
      !jposize = SIZE(obs,1) ; jprsize = SIZE(upens,2)
      !allocate ( ensobseq(1:jposize,1:jprsize), stat=allocok )
      !IF (allocok.NE.0) GOTO 1001

! Compute observation equivalent of input ensemble
      !ensobseq(:,:) = upens(:,:)

! Compute optimality score
      !IF (obserror_type_sesam.EQ.'gaussian') THEN
      !  CALL optimality_score(score, ensobseq, obs, oestd)
      !ELSE
      !  CALL optimality_score(score, ensobseq, obs, cdf_obs)
      !ENDIF

      !IF (jproc.eq.0)  PRINT *, 'OPTIMALITY:',mcmc_index,score

! Deallocate array
      !IF (allocated(ensobseq)) deallocate(ensobseq)

      RETURN
!
 1001 CALL printerror2(0,1001,3,'algomcmc_adap','algomcmc_adap')

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

      call obserror_cdf( o, y, oestd(obs_idx), cdf_obs )

      END FUNCTION cdf_obs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algomcmc_adap
