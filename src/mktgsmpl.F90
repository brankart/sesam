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
! ---                   MKTGSMPL.F90                              ---
! ---                                                           ---
! --- original     : 2006-09 (J.M. Brankart)                    ---
! --- modification : 2008-01 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE tgsmpl : Compute truncated Gaussian sample
! --- SUBROUTINE gsmpl : Compute truncated Gaussian sample
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mktgsmpl
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC tgsmpl, gsmpl

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE tgsmpl(karginxbas,kargincstr,karginvar,kargoutxbas)
!---------------------------------------------------------------------
!
!  Purpose : Compute sample of a TG pdf from the TG parameters
!  -------
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
      use mod_spacexyo , only : jpx, jprend, jpmend, jpsmplend
      use hioxyo
      use hiobas
      use ensdam_storng
      use ensdam_stotge
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) ::  &
     &     karginxbas,kargincstr,karginvar,kargoutxbas
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:,:), allocatable, save :: basexr
      BIGREAL, dimension(:,:), allocatable, save :: basexm
      BIGREAL, dimension(:), allocatable, save :: vectxi
!
      BIGREAL, dimension(:), allocatable, save :: vecbm, vecdm
      BIGREAL, dimension(:,:), allocatable :: matArm
!
      BIGREAL, dimension(:), allocatable :: coef_r
!
      INTEGER, parameter :: jpsmpl1=100
      INTEGER, parameter :: jpiter=10000
      BIGREAL, parameter :: beta=1.5_kr
!
      INTEGER :: flagxyo, jnxyo, allocok
      INTEGER :: jpxsize, jpsmpl
      INTEGER :: jprsize,jrbasdeb,jrbasfin
      INTEGER :: jpmsize,jmbasdeb,jmbasfin
      INTEGER :: jr, jx, jm, jiter, jsmpl
      LOGICAL :: lectinfo
      BIGREAL :: a, b
!
      BIGREAL, dimension(:,:), allocatable :: tgvsmpl
      BIGREAL, dimension(:,:), allocatable :: sample
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modtgop/mktgsmpl :'
         WRITE(numout,*) '         compute TG sample'
      ENDIF
!
! -0.1- Define some parameters :
! ------------------------------
!
      jpxsize=jpx
      jprsize=jprend
      jpmsize=jpmend
      jpsmpl=jpsmplend
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
! -0.2- Allocate arrays :
! -----------------------
!
! --- allocation vectxi
      allocate ( vectxi(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxi(:) = FREAL(0.0)
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
      CALL readvectb(kargincstr,vecbm(:))
!
      CALL readvar(karginvar,vectxi(:),jnxyo,lectinfo,flagxyo)
!
! -2.- Transform problem into reduced space
! -----------------------------------------
!
      DO jm=1,jpmsize
! Compute A and b matrix in reduced space
        DO jr=1,jprsize-1
          matArm(jr,jm)=DOT_PRODUCT(basexm(:,jm),basexr(:,jr+1))
        ENDDO
        vecbm(jm)=vecbm(jm)-DOT_PRODUCT(basexm(:,jm),vectxi(:))
! Normalize A and b matrix in reduced space
        vecdm(jm)=SQRT(DOT_PRODUCT(matArm(:,jm),matArm(:,jm)))
        IF (vecdm(jm).GT.0.0_kr) THEN
          matArm(:,jm)=matArm(:,jm)/vecdm(jm)
          vecbm(jm)=vecbm(jm)/vecdm(jm)
        ENDIF
      ENDDO
      IF (allocated(basexm)) deallocate (basexm)
      IF (allocated(vecdm)) deallocate (vecdm)
!
! -3.- Compute TG sample
! ----------------------
!
! --- allocation vecdm
        allocate ( vecdm(1:jpmsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
!
        CALL kiss_load()
!
! -3.1- Find one initial vector verifying all inequalities
!       (by "the method of projections onto convex sets",
!       see Numerical Recipes - equation 18.5.26
!       Stephen Boyd and Jon Dattoro - last equation on page 5)
!
        coef_r(:) = 0.0_kr
        vecdm(:) = -vecbm(:)
        DO jiter=1,jpiter
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
          IF (jiter.EQ.jpiter) GOTO 113
        ENDDO
!
        IF (allocated(vecdm)) deallocate (vecdm)
!
! -3.2- Compute sample of required size of the TG pdf
!
! --- allocation tgvsmpl
        allocate ( tgvsmpl(1:jpsmpl1,1:jprsize-1), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
! --- allocation sample
        allocate ( sample(1:jpxsize,1:jpsmpl), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
!
! -3.2.1- Compute sample in reduced space
!
        tgvsmpl(1,:) = coef_r(:)
        DO jsmpl=1,jpsmpl
          CALL ranv_tg(tgvsmpl,matArm,vecbm)
          tgvsmpl(1,:) = tgvsmpl(jpsmpl1,:)
          coef_r(:) = tgvsmpl(jpsmpl1,:)
!
! -3.2.2- Transform sample back into original space
!
          DO jx=1,jpxsize
            sample(jx,jsmpl)=vectxi(jx)+ &
     &         DOT_PRODUCT(basexr(jx,2:jprsize),coef_r(1:jprsize-1))
          ENDDO
!
        ENDDO
!
        CALL kiss_save()
!
! -3.3- Write sample in output file
!
        CALL writebas(kargoutxbas,sample(1:jpxsize,1:jpsmpl), &
     &                jnxyo,1,jpsmpl)
!
! --- deallocation tgvsmpl
      IF (allocated(tgvsmpl)) deallocate (tgvsmpl)
      IF (allocated(basexr)) deallocate (basexr)
      IF (allocated(vectxi)) deallocate (vectxi)
      IF (allocated(sample)) deallocate (sample)
      IF (allocated(vecbm)) deallocate (vecbm)
      IF (allocated(matArm)) deallocate (matArm)
      IF (allocated(coef_r)) deallocate (coef_r)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mktgsmpl','tgsmpl')
 1001 CALL printerror2(0,1001,3,'mktgsmpl','tgsmpl')
!
!
 113  WRITE (texterror,*) 'Valid initial vector not found'
      CALL printerror2(0,113,3,'mktgsmpl','tgsmpl',comment=texterror)
!
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE gsmpl(karginxbas,karginvar,kargoutxbas)
!---------------------------------------------------------------------
!
!  Purpose : Compute sample of a Gaussian pdf
!  -------
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
      use mod_spacexyo , only : jpx, jprend, jpsmplend
      use hioxyo
      use hiobas
      use ensdam_storng
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: karginxbas,karginvar,kargoutxbas
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vectxi
      BIGREAL, dimension(:,:), allocatable, save :: basexr
      BIGREAL, dimension(:,:), allocatable, save :: basexm
!
      INTEGER :: flagxyo, jnxyo, allocok
      INTEGER :: jpxsize, jpsmpl
      INTEGER :: jprsize,jr0,jr1,valbase
      INTEGER :: jr, jx, jsmpl
      LOGICAL :: lectinfo,lmodprint
      BIGREAL :: norm
      REAL(KIND=8) :: gran
!
      BIGREAL, dimension(:,:), allocatable :: sample
      BIGREAL, dimension(:,:), allocatable :: coef_ran
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modtgop/mkgsmpl :'
         WRITE(numout,*) '         compute Gaussian sample'
      ENDIF
!
! -0.1- Define some parameters :
! ------------------------------
!
      jpxsize=jpx
      jprsize=jprend
      jpsmpl=jpsmplend
!
      lectinfo=.FALSE.
      flagxyo=1
      jnxyo=1
!
! Read header information from reduced order cov. matrix directory
      CALL readinfobas(karginxbas,valbase)
      IF (valbase.LT.1) THEN
        PRINT *, 'Warning: input directory is an ensemble'
        jr0 = 1 ; jr1 = jprsize
      ELSE
        PRINT *, 'Warning: input directory is a covariance square root'
        jr0 = 2 ; jr1 = jprsize
      ENDIF
!
! -0.2- Allocate arrays :
! -----------------------
!
! --- allocation vectxi
      allocate ( vectxi(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxi(:) = FREAL(0.0)
! --- allocation basexr
      allocate ( basexr(1:jpxsize,jr0:jr1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basexr(:,:) = FREAL(0.0)
! --- allocation coef_ran
      allocate ( coef_ran(1:jpsmpl,jr0:jr1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      coef_ran(:,:) = FREAL(0.0)
! --- allocation sample
      allocate ( sample(1:jpxsize,1:jpsmpl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
! -1.- Generate Gaussian random numbers
! ---------------------------------------------------------
!
      CALL kiss_load()
      DO jsmpl=1,jpsmpl
        DO jr=jr0,jr1
          CALL kiss_gaussian(gran)
          coef_ran(jsmpl,jr) = gran
        ENDDO
        IF (argdisable(1:1).EQ.'T') THEN
          norm=DOT_PRODUCT(coef_ran(jsmpl,jr0:jr1),coef_ran(jsmpl,jr0:jr1))
          norm=norm/FREAL(jr1-jr0+1)
          coef_ran(jsmpl,jr0:jr1)=coef_ran(jsmpl,jr0:jr1)/SQRT(norm)
        ENDIF
        IF (valbase.LT.1) THEN
          coef_ran(jsmpl,jr0:jr1)=coef_ran(jsmpl,jr0:jr1)/SQRT(FREAL(jr1-jr0))
        ENDIF
      ENDDO
      CALL kiss_save()
!
      DO jnxyo=1,limjpnxyo(flagxyo)
        lmodprint=(MOD(jnxyo-1,(limjpnxyo(flagxyo)/10+1)).EQ.0)
        IF (lmodprint) print *,'Memory part number : ', &
     &           jnxyo,'/',limjpnxyo(flagxyo)
!
! -2.- Read input square root covariance matrix or ensemble
! ---------------------------------------------------------
!
        CALL readbas(karginxbas,basexr(:,:),jnxyo,jr0,jr1, &
     &               lectinfo,flagxyo)
!
! -3.- Read mean or compute ensemble mean
! ---------------------------------------
!
        IF (valbase.LT.1) THEN
          DO jx=1,jpxsize
            vectxi(jx) = SUM(basexr(jx,jr0:jr1))/FREAL(jr1-jr0+1)
            basexr(jx,jr0:jr1) = basexr(jx,jr0:jr1) - vectxi(jx)
          ENDDO
        ELSE
          CALL readvar(karginvar,vectxi(:),jnxyo,lectinfo,flagxyo)
        ENDIF
!
! -4.- Compute sample with required mean and square root covariance
! -----------------------------------------------------------------
!
        DO jsmpl=1,jpsmpl
        DO jx=1,jpxsize
          sample(jx,jsmpl)=vectxi(jx)+ &
     &       DOT_PRODUCT(basexr(jx,jr0:jr1),coef_ran(jsmpl,jr0:jr1))
        ENDDO
        ENDDO
!
! -5.- Write sample in output file
! --------------------------------
!
        CALL writebas(kargoutxbas,sample(1:jpxsize,1:jpsmpl), &
     &                jnxyo,1,jpsmpl)
!
      ENDDO
!
! --- deallocation
      IF (allocated(basexr)) deallocate (basexr)
      IF (allocated(vectxi)) deallocate (vectxi)
      IF (allocated(sample)) deallocate (sample)
      IF (allocated(coef_ran)) deallocate (coef_ran)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mktgsmpl','gsmpl')
 1001 CALL printerror2(0,1001,3,'mktgsmpl','gsmpl')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mktgsmpl
