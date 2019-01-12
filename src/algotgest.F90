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
! ---                   ALGOTGEST.F90                             ---
! ---                                                           ---
! --- original     : 2006-09 (C. Lauvernet)                     ---
! --- modification : 2007-09 (F. Castruccio)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE calctgest
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algotgest
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calctgest

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calctgest(karginxbas,kargincstr,karginvar, &
     &                     kargoutvar,kargoutxbas,karginpartvar)
!---------------------------------------------------------------------
!
!  Purpose : Compute maximum likelihood or variance minimizing estimator
!  -------   of a TG pdf from the TG parameters
!  Method :
!  ------
!  Input :
!  -----
!  Output :
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpx, jprend, jpmend
      use hioxyo
      use hiobas
      use ensdam_storng
      use ensdam_stotge
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) ::  &
     &     karginxbas,kargincstr,karginvar,kargoutvar
      CHARACTER(len=*), intent(in), optional :: kargoutxbas
      CHARACTER(len=*), intent(in), optional :: karginpartvar
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:,:), allocatable, save :: basexr
      BIGREAL, dimension(:,:), allocatable, save :: basexm
      BIGREAL, dimension(:), allocatable, save :: vectxi
      BIGREAL, dimension(:), allocatable, save :: vectxo
      BIGREAL, dimension(:), allocatable, save :: vectxpart
      BIGREAL, dimension(:), allocatable, save :: vectxd
!
      BIGREAL, dimension(:), allocatable, save :: vecbm, vecdm
      BIGREAL, dimension(:,:), allocatable :: matbzm, matArm
!
      BIGREAL, dimension(:), allocatable :: coef_r
      BIGREAL, dimension(:), allocatable :: coef_r1
!
      INTEGER, parameter :: jpsmpl1d=1000
      INTEGER, parameter :: jpiter=60, jpiter1=1000, jpiter2=1000
      BIGREAL, parameter :: eps=0.05_kr, beta=1.5_kr
!
      CHARACTER(len=50) :: filecst
      INTEGER :: flagxyo, jnxyo, allocok
      INTEGER :: jpxsize, jpsmpl
      INTEGER :: jprsize,jrbasdeb,jrbasfin
      INTEGER :: jpmsize,jmbasdeb,jmbasfin
      INTEGER :: jpzsize,jpusize
      INTEGER :: jr, jx, jm, jiter, jz, ju, jsmpl
      LOGICAL :: maxlikelihood,lectinfo,localtg,existence,homogenous
      BIGREAL :: a, b
      BIGREAL :: fac, err, itercount
!
      BIGREAL, dimension(:,:), allocatable :: tgvsmpl
      BIGREAL, dimension(:), allocatable :: vectui,vectuo
      BIGREAL, dimension(:), allocatable :: std_r
      BIGREAL, dimension(:,:), allocatable :: baseur,baseum
!
      LOGICAL, dimension(:), allocatable :: maskxpart
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modtgop/algotgest :'
         WRITE(numout,*) '         compute TG estimator'
      ENDIF
!
! -0.1- Define some parameters :
! ------------------------------
      localtg = present(karginpartvar)
      IF (localtg) THEN
        maxlikelihood = .FALSE.
      ELSE
        maxlikelihood = .NOT. present(kargoutxbas)
      ENDIF
!
      filecst=kargincstr(1:len_trim(kargincstr)) // '/vectb'
      INQUIRE(FILE=filecst(1:len_trim(filecst)),EXIST=existence)
      IF (existence) THEN
        homogenous = .TRUE.
      ELSE
        filecst=kargincstr(1:len_trim(kargincstr)) // '/matzmb'
        INQUIRE(FILE=filecst(1:len_trim(filecst)),EXIST=existence) 
        IF (existence) THEN
          homogenous = .FALSE.
        ELSE
! Return an error if be vector (file vectb or matzmb) is not present
          GOTO 114
        ENDIF
      ENDIF
!
      jpxsize=jpx
      jprsize=jprend
      jpmsize=jpmend
      jpsmpl=jpsmpl1d*(jprsize-1)+1
!      
      jrbasdeb = 2
      jrbasfin = jprsize
      jmbasdeb = 1
      jmbasfin = jpmsize
!
      lectinfo=.FALSE.
      flagxyo=1
      jnxyo=1
!
      CALL kiss_load()
!
! -0.2- Allocate arrays :
! -----------------------
!
! --- allocation vectxi
      allocate ( vectxi(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxi(:) = FREAL(0.0)
! --- allocation vectxo
      allocate ( vectxo(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxo(:) = FREAL(0.0)
! --- allocation basexr
      allocate ( basexr(1:jpxsize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basexr(:,:) = FREAL(0.0)
! --- allocation basexm
      allocate ( basexm(1:jpxsize,1:jpmsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basexm(:,:) = FREAL(0.0)
! --- allocation vecbm
      allocate ( vecbm(1:jpmsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecbm(:) = FREAL(0.0)
! --- allocation vecdm
      allocate ( vecdm(1:jpmsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecdm(:) = FREAL(0.0)
! --- allocation matArm
      allocate ( matArm(1:jprsize-1,1:jpmsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      matArm(:,:) = FREAL(0.0)
! --- allocation coef_r
      allocate ( coef_r(1:jprsize-1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      coef_r(:) = FREAL(0.0)
! --- allocation coef_r1
      allocate ( coef_r1(1:jprsize-1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      coef_r1(:) = FREAL(0.0)
! --- allocation std_r
      allocate ( std_r(1:jprsize-1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      std_r(:) = FREAL(0.0)
! --- allocation tgvsmpl
      allocate ( tgvsmpl(1:jpsmpl,1:jprsize-1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tgvsmpl(:,:) = FREAL(0.0)
!
!
      IF (localtg) THEN
!
! --- allocation vectxpart
        allocate ( vectxpart(1:jpxsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectxpart(:) = FREAL(0.0)
!
! --- Read the partition in subdomains
        CALL readxyo (karginpartvar,vectxpart(:),jnxyo, &
     &            lectinfo,flagxyo)

        jpzsize=MAXVAL(vectxpart)
!
      ELSE
!
        jpzsize=1
!
      ENDIF
!
! --- allocation matbzm
      allocate ( matbzm(1:jpzsize,1:jpmsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      matbzm(:,:) = FREAL(0.0)
!
! -1.- Read input data
! --------------------
!
      CALL readbas(karginxbas,basexr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &             lectinfo,flagxyo)
!
      CALL readbas(kargincstr,basexm(:,:),jnxyo,jmbasdeb,jmbasfin, &
     &             lectinfo,flagxyo)
!
      IF (homogenous) THEN
        CALL readvectb(kargincstr,matbzm(1,:))
        IF (localtg) THEN
          DO jz=2,jpzsize
            matbzm(jz,:)=matbzm(1,:)
          ENDDO
        ENDIF
      ELSE
        CALL readmatzmb(kargincstr,matbzm(:,:))
      ENDIF
!
      CALL readvar(karginvar,vectxi(:),jnxyo,lectinfo,flagxyo)
!
! -2.- Loop on subdomains for computing local estimators
!      Compute the input matrices for each subproblem
! ------------------------------------------------------
!
      DO jz = 1+jproc,jpzsize,jpproc
!
        IF (localtg) THEN
!
! --- allocation vectxd
          allocate ( vectxd(1:jpxsize), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
          vectxd(:) = FREAL(0.0)
!
          WHERE (vectxpart == jz)
            vectxd = 1.
          ELSEWHERE
            vectxd = 0.
          ENDWHERE
!
          jpusize = SUM(vectxd)
!
        ELSE
!
          jpusize = jpxsize
!
        ENDIF
!
!  --- allocation vectui
        allocate ( vectui(1:jpusize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectui(:) = FREAL(0.0)
!
!  --- allocation baseur
        allocate ( baseur(1:jpusize,1:jprsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        baseur(:,:) = FREAL(0.0)
!
!  --- allocation baseum
        allocate ( baseum(1:jpusize,1:jpmsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        baseum(:,:) = FREAL(0.0)
!
!  --- allocation vectuo
        allocate ( vectuo(1:jpusize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectuo(:) = FREAL(0.0)
!
        IF (localtg) THEN
!
!  --- allocation maskxpart
          allocate ( maskxpart(1:jpxsize), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
          maskxpart(:) = .FALSE.
!
          maskxpart = vectxd == 1
          vectui = PACK(vectxi,maskxpart)
          DO jr = 1, jprsize
            baseur(:,jr) = PACK(basexr(:,jr),maskxpart)
          ENDDO
          DO jm = 1, jpmsize
            baseum(:,jm) = PACK(basexm(:,jm),maskxpart)
          ENDDO
!
        ELSE
!
          vectui(:) = vectxi(:)
          baseur(:,:) = basexr(:,:)
          baseum(:,:) = basexm(:,:)
!
        ENDIF
!
! -3.- Transform problem into reduced space
! -----------------------------------------
!
        DO jm=1,jpmsize
! Compute A and b matrix in reduced space
          DO jr=1,jprsize-1
            matArm(jr,jm)=DOT_PRODUCT(baseum(:,jm),baseur(:,jr+1))
          ENDDO
          vecbm(jm)=matbzm(jz,jm)-DOT_PRODUCT(baseum(:,jm),vectui(:))
! Normalize A and b matrix in reduced space
          vecdm(jm)=SQRT(DOT_PRODUCT(matArm(:,jm),matArm(:,jm)))
          IF (vecdm(jm).GT.0.0_kr) THEN
            matArm(:,jm)=matArm(:,jm)/vecdm(jm)
            vecbm(jm)=vecbm(jm)/vecdm(jm)
          ENDIF
        ENDDO
!
! -4.- Compute TG estimator
! -------------------------
!
        IF (maxlikelihood) THEN
          IF (ANY(vecbm(:).LT.0.0_kr)) THEN
!
! -4.1- Compute maximum likelihood TG estimator
! ---------------------------------------------
!
! -4.1.1- Find the maximum of the TG pdf
! by "the method of projections onto convex sets"
! see Numerical Recipes - page 805
!     Stephen Boyd and Jon Dattoro - last equation on page 5)
            fac = 0.5_kr ; itercount=0
            coef_r1(:) = 1.0_kr
            DO jiter=1,jpiter2
! Try to decrease the current vector norm
              coef_r(:)=(1.0_kr-fac)*coef_r1(:)
! Project successively onto unverified constraints
              DO jm=1,jpmsize
                vecdm(jm)=DOT_PRODUCT(matArm(:,jm),coef_r(:))-vecbm(jm)
                IF (vecdm(jm).GT.0.0_kr) THEN
                  coef_r(:)=coef_r(:)-beta*vecdm(jm)*matArm(:,jm)
                ENDIF
              ENDDO
! Update decrease factor
              err = MAXVAL(ABS(coef_r(:)/coef_r1(:)-1.0_kr))
              fac = MAX( err , 0.1_kr*fac )
              fac = MAX( fac , 0.1_kr*eps )
! Exit the loop as soon as required accuracy is reached
! (for 4 successive iterations)
              IF (err.LT.eps) THEN
                itercount=itercount+1
                IF (itercount.GE.4) EXIT
              ELSE
                itercount=0
              ENDIF
! Exit the loop with a warning if the maximum number of iterations is reached
              IF (jiter.EQ.jpiter2) THEN
                print *, 'Warning: no convergence in the max estimate'
              ENDIF
! Save the current vector for next iteration
              coef_r1(:)=coef_r(:)
            ENDDO
!
          ENDIF
        ELSE
!
! -4.2- Compute minimum error variance TG estimator
! -------------------------------------------------
!
! -4.2.1- Find one initial vector verifying all inequalities
! (by "the method of projections onto convex sets",
! see Numerical Recipes - equation 18.5.26
!     Stephen Boyd and Jon Dattoro - last equation on page 5)
!
          coef_r(:) = 0.0_kr
          vecdm(:) = -vecbm(:)
          tgvsmpl(:,:) = 0.0_kr
          DO jiter=1,jpiter1
! Project onto the first unverified constraints
            DO jm=1,jpmsize
              IF (vecdm(jm).GT.0.0_kr) THEN
                coef_r(:)=coef_r(:)-beta*vecdm(jm)*matArm(:,jm) ; EXIT
              ENDIF
            ENDDO
! Update constraint misfits (Ax-b)
            DO jm=1,jpmsize
              vecdm(jm)=DOT_PRODUCT(matArm(:,jm),coef_r(:))-vecbm(jm)
            ENDDO
! Exit the loop as soon as all inequalities are verified
            IF (ALL(vecdm(:).LE.0.0_kr)) EXIT
! Return an error if the maximum number of iterations is reached
            IF (jiter.EQ.jpiter1) GOTO 113
          ENDDO
!
! -4.2.2- Compute the mean of the TG pdf with accuracy = eps
!
! Initialize first iteration with
! - a first estimate of the mean 'coef_r1(:)'
! - a valid initial state 'tgvsmpl(1,:)' for the sampling of the TG pdf
          itercount=0
          coef_r1(:) = 0.0_kr
          std_r(:) = 0.0_kr
          tgvsmpl(1,:) = coef_r(:)
! - improve randomness of the initial state of the Gibbs sampler
!   by computing a first sample of unused values
          CALL ranv_tg(tgvsmpl,matArm,vecbm)
          tgvsmpl(1,:) = tgvsmpl(jpsmpl,:)
          DO jiter=1,jpiter
! Compute new sample of the TG pdf
            CALL ranv_tg(tgvsmpl,matArm,vecbm)
! Compute correction to the mean estimate due to this new sample
            DO jr=1,jprsize-1
              coef_r(jr)=SUM(tgvsmpl(2:jpsmpl,jr))/FREAL(jpsmpl-1)
              coef_r(jr)=(coef_r(jr)-coef_r1(jr))/FREAL(jiter)
            ENDDO
! Compute new mean
            coef_r(:)=coef_r1(:)+coef_r(:)
! Compute the square of the difference to mean
            DO jr=1,jprsize-1
              DO jsmpl=2,jpsmpl
                std_r(jr) = std_r(jr) + ( (tgvsmpl(jsmpl,jr)-coef_r(jr)) &
     &                       * (tgvsmpl(jsmpl,jr)-coef_r(jr)))
              ENDDO
            ENDDO
! Exit the loop if last correction is smaller than tolerance
! (for 4 successive iterations)
!           IF (MAXVAL(ABS(coef_r(:)/coef_r1(:))).LT.eps) THEN
!           IF (1.0_kr/SQRT(FREAL(jiter*jpsmpl)).LT.eps) THEN
!           IF (MAXVAL(ABS(coef_r(:))).LT.eps) THEN
            IF (MAXVAL(SQRT(std_r(:))/FREAL(jiter*(jpsmpl-1))).LT.eps) THEN
              itercount=itercount+1
              IF (itercount.GE.4) EXIT
            ELSE
              itercount=0
            ENDIF
! Exit the loop with a warning if the maximum number of iterations is reached
            IF (jiter.EQ.jpiter) THEN
              print *, 'Warning: no convergence in the mean estimate'
            ENDIF
! Initialize next iteration
            coef_r1(:) = coef_r(:)
            tgvsmpl(1,:) = tgvsmpl(jpsmpl,:)
          ENDDO
!
        ENDIF
!
! -5.- Transform problem back into original space
! -----------------------------------------------
!
        DO ju=1,jpusize
          vectuo(ju)=vectui(ju)+DOT_PRODUCT(baseur(ju,2:jprsize), &
     &                                      coef_r(1:jprsize-1))
        ENDDO

        IF (localtg) THEN
          vectxo = UNPACK(vectuo,maskxpart,vectxo)
        ELSE
          vectxo = vectuo
        ENDIF
!
        IF (allocated(vectui)) deallocate (vectui)
        IF (allocated(vectuo)) deallocate (vectuo)
        IF (allocated(baseur)) deallocate (baseur)
        IF (allocated(baseum)) deallocate (baseum)
        IF (allocated(maskxpart)) deallocate (maskxpart)
        IF (allocated(vectxd)) deallocate (vectxd)
!
! --- End loop on subdomains
!
      ENDDO
!
! -6.- Write output data
! ----------------------
#if defined MPI
      CALL mpi_barrier(mpi_comm_world,mpi_code)
      CALL mpi_reduce(vectxo,vectxo,jpxsize,mpi_double_precision, &
     &                mpi_sum,0,mpi_comm_world,mpi_code)
#endif
!
      IF (jproc.EQ.0) THEN
        CALL writevar(kargoutvar,vectxo(:),jnxyo)
      ENDIF
!
      CALL kiss_save()
!
      IF (allocated(basexr)) deallocate (basexr)
      IF (allocated(vectxi)) deallocate (vectxi)
      IF (allocated(vectxo)) deallocate (vectxo)
      IF (allocated(matbzm)) deallocate (matbzm)
      IF (allocated(vecbm)) deallocate (vecbm)
      IF (allocated(matArm)) deallocate (matArm)
      IF (allocated(coef_r)) deallocate (coef_r)
      IF (allocated(basexm)) deallocate (basexm)
      IF (allocated(vecdm)) deallocate (vecdm)
      IF (allocated(coef_r1)) deallocate (coef_r1)
      IF (allocated(tgvsmpl)) deallocate (tgvsmpl)
      IF (allocated(std_r)) deallocate (std_r)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algotgest','calctgest')
 1001 CALL printerror2(0,1001,3,'algotgest','calctgest')
!
!
 112  WRITE (texterror,*) 'This action is not yet available'
      CALL printerror2(0,112,3,'algotgest','calctgest',comment=texterror)
 113  WRITE (texterror,*) 'Valid initial vector not found'
      CALL printerror2(0,113,3,'algotgest','calctgest',comment=texterror)
 114  WRITE (texterror,*) 'Constraint b vector (vectb or matzmb) not found'
      CALL printerror2(0,114,3,'algotgest','algotgest',comment=texterror)
!
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algotgest
