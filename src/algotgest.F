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
C ---                   ALGOTGEST.F                             ---
C ---                                                           ---
C --- original     : 2006-09 (C. Lauvernet)                     ---
C --- modification : 2007-09 (F. Castruccio)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE calctgest
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algotgest
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calctgest

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calctgest(karginxbas,kargincstr,karginvar,
     $                     kargoutvar,kargoutxbas,karginpartvar)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute maximum likelihood or variance minimizing estimator
CCC  -------   of a TG pdf from the TG parameters
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
      use mod_spacexyo , only : jpx, jprend, jpmend
      use hioxyo
      use hiobas
      use ensdam_storng
      use ensdam_stotge
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: 
     $     karginxbas,kargincstr,karginvar,kargoutvar
      CHARACTER(len=*), intent(in), optional :: kargoutxbas
      CHARACTER(len=*), intent(in), optional :: karginpartvar
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:,:), allocatable, save :: basexr
      BIGREAL, dimension(:,:), allocatable, save :: basexm
      BIGREAL, dimension(:), allocatable, save :: vectxi
      BIGREAL, dimension(:), allocatable, save :: vectxo
      BIGREAL, dimension(:), allocatable, save :: vectxpart
      BIGREAL, dimension(:), allocatable, save :: vectxd
C
      BIGREAL, dimension(:), allocatable, save :: vecbm, vecdm
      BIGREAL, dimension(:,:), allocatable :: matbzm, matArm
C
      BIGREAL, dimension(:), allocatable :: coef_r
      BIGREAL, dimension(:), allocatable :: coef_r1
C
      INTEGER, parameter :: jpsmpl1d=1000
      INTEGER, parameter :: jpiter=60, jpiter1=1000, jpiter2=1000
      BIGREAL, parameter :: eps=0.05_kr, beta=1.5_kr
C
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
C
      BIGREAL, dimension(:,:), allocatable :: tgvsmpl
      BIGREAL, dimension(:), allocatable :: vectui,vectuo
      BIGREAL, dimension(:), allocatable :: std_r
      BIGREAL, dimension(:,:), allocatable :: baseur,baseum
C
      LOGICAL, dimension(:), allocatable :: maskxpart
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modtgop/algotgest :'
         WRITE(numout,*) '         compute TG estimator'
      ENDIF
C
C -0.1- Define some parameters :
C ------------------------------
      localtg = present(karginpartvar)
      IF (localtg) THEN
        maxlikelihood = .FALSE.
      ELSE
        maxlikelihood = .NOT. present(kargoutxbas)
      ENDIF
C
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
C Return an error if be vector (file vectb or matzmb) is not present
          GOTO 114
        ENDIF
      ENDIF
C
      jpxsize=jpx
      jprsize=jprend
      jpmsize=jpmend
      jpsmpl=jpsmpl1d*(jprsize-1)+1
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
      CALL kiss_load()
C
C -0.2- Allocate arrays :
C -----------------------
C
C --- allocation vectxi
      allocate ( vectxi(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxi(:) = FREAL(0.0)
C --- allocation vectxo
      allocate ( vectxo(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxo(:) = FREAL(0.0)
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
C --- allocation coef_r1
      allocate ( coef_r1(1:jprsize-1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      coef_r1(:) = FREAL(0.0)
C --- allocation std_r
      allocate ( std_r(1:jprsize-1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      std_r(:) = FREAL(0.0)
C --- allocation tgvsmpl
      allocate ( tgvsmpl(1:jpsmpl,1:jprsize-1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tgvsmpl(:,:) = FREAL(0.0)
C
C
      IF (localtg) THEN
C
C --- allocation vectxpart
        allocate ( vectxpart(1:jpxsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectxpart(:) = FREAL(0.0)
C
C --- Read the partition in subdomains
        CALL readxyo (karginpartvar,vectxpart(:),jnxyo,
     $            lectinfo,flagxyo)

        jpzsize=MAXVAL(vectxpart)
C
      ELSE
C
        jpzsize=1
C
      ENDIF
C
C --- allocation matbzm
      allocate ( matbzm(1:jpzsize,1:jpmsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      matbzm(:,:) = FREAL(0.0)
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
C
      CALL readvar(karginvar,vectxi(:),jnxyo,lectinfo,flagxyo)
C
C -2.- Loop on subdomains for computing local estimators
C      Compute the input matrices for each subproblem
C ------------------------------------------------------
C
      DO jz = 1+jproc,jpzsize,jpproc
C
        IF (localtg) THEN
C
C --- allocation vectxd
          allocate ( vectxd(1:jpxsize), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
          vectxd(:) = FREAL(0.0)
C
          WHERE (vectxpart == jz)
            vectxd = 1.
          ELSEWHERE
            vectxd = 0.
          ENDWHERE
C
          jpusize = SUM(vectxd)
C
        ELSE
C
          jpusize = jpxsize
C
        ENDIF
C
C  --- allocation vectui
        allocate ( vectui(1:jpusize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectui(:) = FREAL(0.0)
C
C  --- allocation baseur
        allocate ( baseur(1:jpusize,1:jprsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        baseur(:,:) = FREAL(0.0)
C
C  --- allocation baseum
        allocate ( baseum(1:jpusize,1:jpmsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        baseum(:,:) = FREAL(0.0)
C
C  --- allocation vectuo
        allocate ( vectuo(1:jpusize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectuo(:) = FREAL(0.0)
C
        IF (localtg) THEN
C
C  --- allocation maskxpart
          allocate ( maskxpart(1:jpxsize), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
          maskxpart(:) = .FALSE.
C
          maskxpart = vectxd == 1
          vectui = PACK(vectxi,maskxpart)
          DO jr = 1, jprsize
            baseur(:,jr) = PACK(basexr(:,jr),maskxpart)
          ENDDO
          DO jm = 1, jpmsize
            baseum(:,jm) = PACK(basexm(:,jm),maskxpart)
          ENDDO
C
        ELSE
C
          vectui(:) = vectxi(:)
          baseur(:,:) = basexr(:,:)
          baseum(:,:) = basexm(:,:)
C
        ENDIF
C
C -3.- Transform problem into reduced space
C -----------------------------------------
C
        DO jm=1,jpmsize
C Compute A and b matrix in reduced space
          DO jr=1,jprsize-1
            matArm(jr,jm)=DOT_PRODUCT(baseum(:,jm),baseur(:,jr+1))
          ENDDO
          vecbm(jm)=matbzm(jz,jm)-DOT_PRODUCT(baseum(:,jm),vectui(:))
C Normalize A and b matrix in reduced space
          vecdm(jm)=SQRT(DOT_PRODUCT(matArm(:,jm),matArm(:,jm)))
          IF (vecdm(jm).GT.0.0_kr) THEN
            matArm(:,jm)=matArm(:,jm)/vecdm(jm)
            vecbm(jm)=vecbm(jm)/vecdm(jm)
          ENDIF
        ENDDO
C
C -4.- Compute TG estimator
C -------------------------
C
        IF (maxlikelihood) THEN
          IF (ANY(vecbm(:).LT.0.0_kr)) THEN
C
C -4.1- Compute maximum likelihood TG estimator
C ---------------------------------------------
C
C -4.1.1- Find the maximum of the TG pdf
C by "the method of projections onto convex sets"
C see Numerical Recipes - page 805
C     Stephen Boyd and Jon Dattoro - last equation on page 5)
            fac = 0.5_kr ; itercount=0
            coef_r1(:) = 1.0_kr
            DO jiter=1,jpiter2
C Try to decrease the current vector norm
              coef_r(:)=(1.0_kr-fac)*coef_r1(:)
C Project successively onto unverified constraints
              DO jm=1,jpmsize
                vecdm(jm)=DOT_PRODUCT(matArm(:,jm),coef_r(:))-vecbm(jm)
                IF (vecdm(jm).GT.0.0_kr) THEN
                  coef_r(:)=coef_r(:)-beta*vecdm(jm)*matArm(:,jm)
                ENDIF
              ENDDO
C Update decrease factor
              err = MAXVAL(ABS(coef_r(:)/coef_r1(:)-1.0_kr))
              fac = MAX( err , 0.1_kr*fac )
              fac = MAX( fac , 0.1_kr*eps )
C Exit the loop as soon as required accuracy is reached
C (for 4 successive iterations)
              IF (err.LT.eps) THEN
                itercount=itercount+1
                IF (itercount.GE.4) EXIT
              ELSE
                itercount=0
              ENDIF
C Exit the loop with a warning if the maximum number of iterations is reached
              IF (jiter.EQ.jpiter2) THEN
                print *, 'Warning: no convergence in the max estimate'
              ENDIF
C Save the current vector for next iteration
              coef_r1(:)=coef_r(:)
            ENDDO
C
          ENDIF
        ELSE
C
C -4.2- Compute minimum error variance TG estimator
C -------------------------------------------------
C
C -4.2.1- Find one initial vector verifying all inequalities
C (by "the method of projections onto convex sets",
C see Numerical Recipes - equation 18.5.26
C     Stephen Boyd and Jon Dattoro - last equation on page 5)
C
          coef_r(:) = 0.0_kr
          vecdm(:) = -vecbm(:)
          tgvsmpl(:,:) = 0.0_kr
          DO jiter=1,jpiter1
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
            IF (jiter.EQ.jpiter1) GOTO 113
          ENDDO
C
C -4.2.2- Compute the mean of the TG pdf with accuracy = eps
C
C Initialize first iteration with
C - a first estimate of the mean 'coef_r1(:)'
C - a valid initial state 'tgvsmpl(1,:)' for the sampling of the TG pdf
          itercount=0
          coef_r1(:) = 0.0_kr
          std_r(:) = 0.0_kr
          tgvsmpl(1,:) = coef_r(:)
C - improve randomness of the initial state of the Gibbs sampler
C   by computing a first sample of unused values
          CALL ranv_tg(tgvsmpl,matArm,vecbm)
          tgvsmpl(1,:) = tgvsmpl(jpsmpl,:)
          DO jiter=1,jpiter
C Compute new sample of the TG pdf
            CALL ranv_tg(tgvsmpl,matArm,vecbm)
C Compute correction to the mean estimate due to this new sample
            DO jr=1,jprsize-1
              coef_r(jr)=SUM(tgvsmpl(2:jpsmpl,jr))/FREAL(jpsmpl-1)
              coef_r(jr)=(coef_r(jr)-coef_r1(jr))/FREAL(jiter)
            ENDDO
C Compute new mean
            coef_r(:)=coef_r1(:)+coef_r(:)
C Compute the square of the difference to mean
            DO jr=1,jprsize-1
              DO jsmpl=2,jpsmpl
                std_r(jr) = std_r(jr) + ( (tgvsmpl(jsmpl,jr)-coef_r(jr))
     $                       * (tgvsmpl(jsmpl,jr)-coef_r(jr)))
              ENDDO
            ENDDO
C Exit the loop if last correction is smaller than tolerance
C (for 4 successive iterations)
C           IF (MAXVAL(ABS(coef_r(:)/coef_r1(:))).LT.eps) THEN
C           IF (1.0_kr/SQRT(FREAL(jiter*jpsmpl)).LT.eps) THEN
C           IF (MAXVAL(ABS(coef_r(:))).LT.eps) THEN
            IF (MAXVAL(SQRT(std_r(:))/FREAL(jiter*(jpsmpl-1))).LT.eps) THEN
              itercount=itercount+1
              IF (itercount.GE.4) EXIT
            ELSE
              itercount=0
            ENDIF
C Exit the loop with a warning if the maximum number of iterations is reached
            IF (jiter.EQ.jpiter) THEN
              print *, 'Warning: no convergence in the mean estimate'
            ENDIF
C Initialize next iteration
            coef_r1(:) = coef_r(:)
            tgvsmpl(1,:) = tgvsmpl(jpsmpl,:)
          ENDDO
C
        ENDIF
C
C -5.- Transform problem back into original space
C -----------------------------------------------
C
        DO ju=1,jpusize
          vectuo(ju)=vectui(ju)+DOT_PRODUCT(baseur(ju,2:jprsize),
     $                                      coef_r(1:jprsize-1))
        ENDDO

        IF (localtg) THEN
          vectxo = UNPACK(vectuo,maskxpart,vectxo)
        ELSE
          vectxo = vectuo
        ENDIF
C
        IF (allocated(vectui)) deallocate (vectui)
        IF (allocated(vectuo)) deallocate (vectuo)
        IF (allocated(baseur)) deallocate (baseur)
        IF (allocated(baseum)) deallocate (baseum)
        IF (allocated(maskxpart)) deallocate (maskxpart)
        IF (allocated(vectxd)) deallocate (vectxd)
C
C --- End loop on subdomains
C
      ENDDO
C
C -6.- Write output data
C ----------------------
#if defined MPI
      CALL mpi_barrier(mpi_comm_world,mpi_code)
      CALL mpi_reduce(vectxo,vectxo,jpxsize,mpi_double_precision,
     $                mpi_sum,0,mpi_comm_world,mpi_code)
#endif
C
      IF (jproc.EQ.0) THEN
        CALL writevar(kargoutvar,vectxo(:),jnxyo)
      ENDIF
C
      CALL kiss_save()
C
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
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algotgest','calctgest')
 1001 CALL printerror2(0,1001,3,'algotgest','calctgest')
C
C
 112  WRITE (texterror,*) 'This action is not yet available'
      CALL printerror2(0,112,3,'algotgest','calctgest',comment=texterror)
 113  WRITE (texterror,*) 'Valid initial vector not found'
      CALL printerror2(0,113,3,'algotgest','calctgest',comment=texterror)
 114  WRITE (texterror,*) 'Constraint b vector (vectb or matzmb) not found'
      CALL printerror2(0,114,3,'algotgest','algotgest',comment=texterror)
C
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algotgest
