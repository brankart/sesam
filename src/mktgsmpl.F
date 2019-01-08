C Copyright: CNRS - Université de Grenoble
C
C Contributors : Jean-Michel Brankart, Charles-Emmanuel Testut, Laurent Parent,
C                Emmanuel Cosme, Claire Lauvernet, Frédéric Castruccio
C
C Jean-Michel.Brankart@hmg.inpg.fr
C
C This software is governed by the CeCILL license under French law and
C abiding by the rules of distribution of free software.  You can  use,
C modify and/ or redistribute the software under the terms of the CeCILL
C license as circulated by CEA, CNRS and INRIA at the following URL
C "http://www.cecill.info".
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C ---                                                           ---
C ---                   MKTGSMPL.F                              ---
C ---                                                           ---
C --- original     : 2006-09 (J.M. Brankart)                    ---
C --- modification : 2008-01 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE tgsmpl : Compute truncated Gaussian sample
C --- SUBROUTINE gsmpl : Compute truncated Gaussian sample
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mktgsmpl
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC tgsmpl, gsmpl

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE tgsmpl(karginxbas,kargincstr,karginvar,kargoutxbas)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute sample of a TG pdf from the TG parameters
CCC  -------
CCC  Method :
CCC  ------
CCC  Input :
CCC  -----
CCC  Output :
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpx, jprend, jpmend, jpsmplend
      use hioxyo
      use hiobas
      use ensdam_storng
      use ensdam_stotge
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: 
     $     karginxbas,kargincstr,karginvar,kargoutxbas
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:,:), allocatable, save :: basexr
      BIGREAL, dimension(:,:), allocatable, save :: basexm
      BIGREAL, dimension(:), allocatable, save :: vectxi
C
      BIGREAL, dimension(:), allocatable, save :: vecbm, vecdm
      BIGREAL, dimension(:,:), allocatable :: matArm
C
      BIGREAL, dimension(:), allocatable :: coef_r
C
      INTEGER, parameter :: jpsmpl1=100
      INTEGER, parameter :: jpiter=10000
      BIGREAL, parameter :: beta=1.5_kr
C
      INTEGER :: flagxyo, jnxyo, allocok
      INTEGER :: jpxsize, jpsmpl
      INTEGER :: jprsize,jrbasdeb,jrbasfin
      INTEGER :: jpmsize,jmbasdeb,jmbasfin
      INTEGER :: jr, jx, jm, jiter, jsmpl
      LOGICAL :: lectinfo
      BIGREAL :: a, b
C
      BIGREAL, dimension(:,:), allocatable :: tgvsmpl
      BIGREAL, dimension(:,:), allocatable :: sample
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modtgop/mktgsmpl :'
         WRITE(numout,*) '         compute TG sample'
      ENDIF
C
C -0.1- Define some parameters :
C ------------------------------
C
      jpxsize=jpx
      jprsize=jprend
      jpmsize=jpmend
      jpsmpl=jpsmplend
C
      jrbasdeb = 2
      jrbasfin = jprsize
      jmbasdeb = 1
      jmbasfin = jpmsize
C
      lectinfo=.FALSE.
      flagxyo=1
      jnxyo=1
C
C -0.2- Allocate arrays :
C -----------------------
C
C --- allocation vectxi
      allocate ( vectxi(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxi(:) = FREAL(0.0)
C --- allocation basexr
      allocate ( basexr(1:jpxsize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basexr(:,:) = FREAL(0.0)
C --- allocation basexm
      allocate ( basexm(1:jpxsize,1:jpmsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basexm(:,:) = FREAL(0.0)
C --- allocation vecbm
      allocate ( vecbm(1:jpmsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecbm(:) = FREAL(0.0)
C --- allocation vecdm
      allocate ( vecdm(1:jpmsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecdm(:) = FREAL(0.0)
C --- allocation matArm
      allocate ( matArm(1:jprsize-1,1:jpmsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      matArm(:,:) = FREAL(0.0)
C --- allocation coef_r
      allocate ( coef_r(1:jprsize-1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      coef_r(:) = FREAL(0.0)
C
C -1.- Read input data
C --------------------
C
      CALL readbas(karginxbas,basexr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $             lectinfo,flagxyo)
C
      CALL readbas(kargincstr,basexm(:,:),jnxyo,jmbasdeb,jmbasfin,
     $             lectinfo,flagxyo)
C
      CALL readvectb(kargincstr,vecbm(:))
C
      CALL readvar(karginvar,vectxi(:),jnxyo,lectinfo,flagxyo)
C
C -2.- Transform problem into reduced space
C -----------------------------------------
C
      DO jm=1,jpmsize
C Compute A and b matrix in reduced space
        DO jr=1,jprsize-1
          matArm(jr,jm)=DOT_PRODUCT(basexm(:,jm),basexr(:,jr+1))
        ENDDO
        vecbm(jm)=vecbm(jm)-DOT_PRODUCT(basexm(:,jm),vectxi(:))
C Normalize A and b matrix in reduced space
        vecdm(jm)=SQRT(DOT_PRODUCT(matArm(:,jm),matArm(:,jm)))
        IF (vecdm(jm).GT.0.0_kr) THEN
          matArm(:,jm)=matArm(:,jm)/vecdm(jm)
          vecbm(jm)=vecbm(jm)/vecdm(jm)
        ENDIF
      ENDDO
      IF (allocated(basexm)) deallocate (basexm)
      IF (allocated(vecdm)) deallocate (vecdm)
C
C -3.- Compute TG sample
C ----------------------
C
C --- allocation vecdm
        allocate ( vecdm(1:jpmsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
C
        CALL kiss_load()
C
C -3.1- Find one initial vector verifying all inequalities
C       (by "the method of projections onto convex sets",
C       see Numerical Recipes - equation 18.5.26
C       Stephen Boyd and Jon Dattoro - last equation on page 5)
C
        coef_r(:) = 0.0_kr
        vecdm(:) = -vecbm(:)
        DO jiter=1,jpiter
C Project onto the first unverified constraints
          DO jm=1,jpmsize
            IF (vecdm(jm).GT.0.0_kr) THEN
              coef_r(:)=coef_r(:)-beta*vecdm(jm)*matArm(:,jm) ; EXIT
            ENDIF
          ENDDO
C Update constraint misfits (Ax-b)
          DO jm=1,jpmsize
            vecdm(jm)=DOT_PRODUCT(matArm(:,jm),coef_r(:))-vecbm(jm)
          ENDDO
C Exit the loop as soon as all inequalities are verified
          IF (ALL(vecdm(:).LE.0.0_kr)) EXIT
C Return an error if the maximum number of iterations is reached
          IF (jiter.EQ.jpiter) GOTO 113
        ENDDO
C
        IF (allocated(vecdm)) deallocate (vecdm)
C
C -3.2- Compute sample of required size of the TG pdf
C
C --- allocation tgvsmpl
        allocate ( tgvsmpl(1:jpsmpl1,1:jprsize-1), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
C --- allocation sample
        allocate ( sample(1:jpxsize,1:jpsmpl), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
C
C -3.2.1- Compute sample in reduced space
C
        tgvsmpl(1,:) = coef_r(:)
        DO jsmpl=1,jpsmpl
          CALL ranv_tg(tgvsmpl,matArm,vecbm)
          tgvsmpl(1,:) = tgvsmpl(jpsmpl1,:)
          coef_r(:) = tgvsmpl(jpsmpl1,:)
C
C -3.2.2- Transform sample back into original space
C
          DO jx=1,jpxsize
            sample(jx,jsmpl)=vectxi(jx)+
     $         DOT_PRODUCT(basexr(jx,2:jprsize),coef_r(1:jprsize-1))
          ENDDO
C
        ENDDO
C
        CALL kiss_save()
C
C -3.3- Write sample in output file
C
        CALL writebas(kargoutxbas,sample(1:jpxsize,1:jpsmpl),
     $                jnxyo,1,jpsmpl)
C
C --- deallocation tgvsmpl
      IF (allocated(tgvsmpl)) deallocate (tgvsmpl)
      IF (allocated(basexr)) deallocate (basexr)
      IF (allocated(vectxi)) deallocate (vectxi)
      IF (allocated(sample)) deallocate (sample)
      IF (allocated(vecbm)) deallocate (vecbm)
      IF (allocated(matArm)) deallocate (matArm)
      IF (allocated(coef_r)) deallocate (coef_r)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mktgsmpl','tgsmpl')
 1001 CALL printerror2(0,1001,3,'mktgsmpl','tgsmpl')
C
C
 113  WRITE (texterror,*) 'Valid initial vector not found'
      CALL printerror2(0,113,3,'mktgsmpl','tgsmpl',comment=texterror)
C
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE gsmpl(karginxbas,karginvar,kargoutxbas)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute sample of a Gaussian pdf
CCC  -------
CCC  Method :
CCC  ------
CCC  Input :
CCC  -----
CCC  Output :
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpx, jprend, jpsmplend
      use hioxyo
      use hiobas
      use ensdam_storng
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: karginxbas,karginvar,kargoutxbas
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vectxi
      BIGREAL, dimension(:,:), allocatable, save :: basexr
      BIGREAL, dimension(:,:), allocatable, save :: basexm
C
      INTEGER :: flagxyo, jnxyo, allocok
      INTEGER :: jpxsize, jpsmpl
      INTEGER :: jprsize,jr0,jr1,valbase
      INTEGER :: jr, jx, jsmpl
      LOGICAL :: lectinfo,lmodprint
      BIGREAL :: norm
      REAL(KIND=8) :: gran
C
      BIGREAL, dimension(:,:), allocatable :: sample
      BIGREAL, dimension(:,:), allocatable :: coef_ran
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modtgop/mkgsmpl :'
         WRITE(numout,*) '         compute Gaussian sample'
      ENDIF
C
C -0.1- Define some parameters :
C ------------------------------
C
      jpxsize=jpx
      jprsize=jprend
      jpsmpl=jpsmplend
C
      lectinfo=.FALSE.
      flagxyo=1
      jnxyo=1
C
C Read header information from reduced order cov. matrix directory
      CALL readinfobas(karginxbas,valbase)
      IF (valbase.LT.1) THEN
        PRINT *, 'Warning: input directory is an ensemble'
        jr0 = 1 ; jr1 = jprsize
      ELSE
        PRINT *, 'Warning: input directory is a covariance square root'
        jr0 = 2 ; jr1 = jprsize
      ENDIF
C
C -0.2- Allocate arrays :
C -----------------------
C
C --- allocation vectxi
      allocate ( vectxi(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxi(:) = FREAL(0.0)
C --- allocation basexr
      allocate ( basexr(1:jpxsize,jr0:jr1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basexr(:,:) = FREAL(0.0)
C --- allocation coef_ran
      allocate ( coef_ran(1:jpsmpl,jr0:jr1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      coef_ran(:,:) = FREAL(0.0)
C --- allocation sample
      allocate ( sample(1:jpxsize,1:jpsmpl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
C -1.- Generate Gaussian random numbers
C ---------------------------------------------------------
C
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
C
      DO jnxyo=1,limjpnxyo(flagxyo)
        lmodprint=(MOD(jnxyo-1,(limjpnxyo(flagxyo)/10+1)).EQ.0)
        IF (lmodprint) print *,'Memory part number : ',
     $           jnxyo,'/',limjpnxyo(flagxyo)
C
C -2.- Read input square root covariance matrix or ensemble
C ---------------------------------------------------------
C
        CALL readbas(karginxbas,basexr(:,:),jnxyo,jr0,jr1,
     $               lectinfo,flagxyo)
C
C -3.- Read mean or compute ensemble mean
C ---------------------------------------
C
        IF (valbase.LT.1) THEN
          DO jx=1,jpxsize
            vectxi(jx) = SUM(basexr(jx,jr0:jr1))/FREAL(jr1-jr0+1)
            basexr(jx,jr0:jr1) = basexr(jx,jr0:jr1) - vectxi(jx)
          ENDDO
        ELSE
          CALL readvar(karginvar,vectxi(:),jnxyo,lectinfo,flagxyo)
        ENDIF
C
C -4.- Compute sample with required mean and square root covariance
C -----------------------------------------------------------------
C
        DO jsmpl=1,jpsmpl
        DO jx=1,jpxsize
          sample(jx,jsmpl)=vectxi(jx)+
     $       DOT_PRODUCT(basexr(jx,jr0:jr1),coef_ran(jsmpl,jr0:jr1))
        ENDDO
        ENDDO
C
C -5.- Write sample in output file
C --------------------------------
C
        CALL writebas(kargoutxbas,sample(1:jpxsize,1:jpsmpl),
     $                jnxyo,1,jpsmpl)
C
      ENDDO
C
C --- deallocation
      IF (allocated(basexr)) deallocate (basexr)
      IF (allocated(vectxi)) deallocate (vectxi)
      IF (allocated(sample)) deallocate (sample)
      IF (allocated(coef_ran)) deallocate (coef_ran)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mktgsmpl','gsmpl')
 1001 CALL printerror2(0,1001,3,'mktgsmpl','gsmpl')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mktgsmpl
