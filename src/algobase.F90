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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ---                                                           ---
! ---                    ALGOBASE.F90                             ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 00-01 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  calceofs
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algobase
      use mod_main
      use algortho
      use mkmeanect
      use mkconnect
      IMPLICIT NONE
      PRIVATE

      PUBLIC calceofs

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calceofs(kflaganlxyo,kinbasxyo,koutbasxyo, &
     &     kflagcovxyo, &
     &     kflagincovxyoz,kincovxyoz, &
     &     kflagbicovoz,koutcovoz, &
     &     kflaggloloc,kinpartxyo,kinzon,kconfigo, &
     &     ktextdisable,kinbasrefxyoz)
!---------------------------------------------------------------------
!
!  Purpose : Implement a global & local version of the algorithm
!            which computes EOFs 
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
      use mod_spacexyo , only : jpoend,jpitpend,jprend,jpx,jpxend, &
     &     jpyend,jpz,poscoefobs,gridijkobs, &
     &     pt1bubidx, pt2bubidx, pt3bubidx, pt4bubidx,  &
     &     pt1dtalon, pt1dtalat, pt1dtadepth, pt1dtatime, &
     &     pt1bublon, pt1bublat, pt1bubdepth, pt1bubtime, &
     &     vectptbub, ptbub, &
     &     bubblk1, bubblk2, bubblk3, bubblk4,spvaldta
      use hioxyo
      use hiobas
      use hiozon
      use utilroa
      use utilmkh
      use utilmkto
      use utilclc
      use utilprint
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(in) :: kflaganlxyo,kflagcovxyo, &
     &     kflagincovxyoz,kflagbicovoz,kflaggloloc
      CHARACTER(len=*), intent(in) :: kinbasxyo,koutbasxyo, &
     &     kincovxyoz,koutcovoz,kinpartxyo,kinzon,kconfigo, &
     &     ktextdisable
      CHARACTER(len=*), optional, intent(in) :: kinbasrefxyoz
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vectsweight
      BIGREAL, dimension(:), allocatable, save :: vectsweightz
      BIGREAL, dimension(:), allocatable, save :: vectsmean
      BIGREAL, dimension(:), allocatable, save :: vectsect
      BIGREAL, dimension(:,:), allocatable, save :: basesr
      BIGREAL, dimension(:,:), allocatable, save :: mat_yr
      BIGREAL, dimension(:), allocatable, save :: vecty
      BIGREAL, dimension(:), allocatable, save :: vectspart
      BIGREAL, dimension(:), allocatable, save :: vectybub
      BIGREAL, dimension(:), allocatable, save :: vectsstd
      BIGREAL, dimension(:), allocatable, save :: vectsmode
      BIGREAL, dimension(:,:), allocatable, save :: basesd
      BIGREAL, dimension(:,:), allocatable, save :: basexr
      BIGREAL, dimension(:), allocatable, save :: vectxmean
      BIGREAL, dimension(:), allocatable, save :: vectxect
      BIGREAL, dimension(:), allocatable, save :: vectxpart
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      BIGREAL, dimension(:,:), allocatable :: lambda
      BIGREAL, dimension(:,:), allocatable :: amplitude
      BIGREAL, dimension(:,:), allocatable :: covrd
      BIGREAL, dimension(:,:), allocatable, save :: covrr
!
      BIGREAL, dimension(:,:,:), allocatable :: matUrrz
      BIGREAL, dimension(:,:,:), allocatable :: covrrz
      BIGREAL, dimension(:,:,:), allocatable :: covrdz
!
      BIGREAL, dimension(:,:), allocatable :: matBqr, matBqd
!
      INTEGER :: jprbasrefin,valbasref,jdbasdeb,jdbasfin,jd
!
      INTEGER :: allocok,jpxsize,jpssize,jpysize
      INTEGER :: jprsize,jpitpsize,jpzsize,jpdsize,jpqsize
      INTEGER :: jprbasin,jprbasout,jpnbubsize
      INTEGER :: jnxyo,jr,jz,ios,valbase,nvpnull,nvectindep
      INTEGER :: jrbasdeb,jrbasdeb1,jrbasfin,jrbas,jrbas1,jrbas2,jext
      INTEGER :: serie,numjr,flagcfg,flagxyo
      CHARACTER(len=bgword) :: fnamein,fnameout,vctnamin,vctnamout
      CHARACTER(len=bgword), dimension (:), allocatable :: tabfnout, &
     &     tabfnin
      CHARACTER(len=bgword) :: titre,variable,dirnambas,text
      LOGICAL :: lmodprintjnxyo,lmodprintjnz,zlectinfo,lectinfo
      LOGICAL :: lmoyectold,all
      LOGICAL :: lPaout,lWoestd,lweighto,lweightoe,lweight,lweightz
      INTEGER :: jbub,jdta,jprbasout1,jprbasin1,ji,jj,jk,jt,js
      INTEGER :: zon_jpi1,zon_jpj1,zon_jpk1,zon_jpt1,jpbub1,jpz1
      INTEGER :: jzfin,jzend,jbub1,jbubdeb,jbubfin,jz1
      INTEGER :: jpbubsize1,jbubend,jpu
      INTEGER :: jnz, jz2, jsbasfin, jq, jqbas1, jqmat1, jqbas2, jdbas2
      LOGICAL :: allmemjpzr
      INTEGER :: jqbasdeb,jqbasfin,jqmatdeb,jqmatfin
      BIGREAL :: coeflimite
      INTEGER, dimension (:,:), allocatable :: ptbubidx
      INTEGER, dimension (:,:), allocatable :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), allocatable :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), allocatable :: ptbublon, ptbublat
      INTEGER, dimension (:,:), allocatable :: ptbubdepth, ptbubtime
      INTEGER, dimension (:,:), allocatable :: ptbubidx1
      INTEGER, dimension (:,:), allocatable :: ptdtalon1, ptdtalat1
      INTEGER, dimension (:,:), allocatable :: ptdtadepth1, ptdtatime1
      INTEGER, dimension (:,:), allocatable :: ptbublon1, ptbublat1
      INTEGER, dimension (:,:), allocatable :: ptbubdepth1, ptbubtime1
!
      INTEGER, dimension(:), allocatable :: vectynbpoint
      BIGREAL, dimension (:,:,:,:,:), allocatable :: bubmode,bubstd
      BIGREAL, dimension (:,:,:,:,:,:), allocatable :: bubr
!----------------------------------------------------------------------
      INTEGER, dimension(:), allocatable :: ptu
      BIGREAL, dimension(:,:), allocatable :: baseur,baseud,baseuq
      BIGREAL, dimension(:), allocatable :: vectuweight
      BIGREAL, dimension(:), allocatable :: vectumean
      BIGREAL, dimension(:), allocatable :: vectuect
      TYPE (type_poscoef), dimension(:,:), allocatable :: poscoefuobs
!----------------------------------------------------------------------
      INTEGER :: fixjpbub
!----------------------------------------------------------------------
      lPaout=(ktextdisable(1:1).EQ.'T')
      lWoestd=(ktextdisable(2:2).EQ.'F')
!
      lweightz = (kflaggloloc.NE.0)
      lweighto = ( (largoestd) .OR. &
     &     ((lWoestd).AND.(kflagcovxyo.GE.2)) )
      lweightoe=(lweightz.AND.largoecorrel.AND.lweighto)
      lweight  = ( (lweighto) .OR. (largweight)  &
     &     .OR. (lweightoe) )
!
      jprsize=jprend
      jpysize=jpyend
      SELECT CASE (kflaganlxyo)
      CASE(1)
         jpxsize=jpx
      CASE(2)
         jpxsize=jpyend
      CASE(3)
         jpxsize=jpoend
      CASE DEFAULT
         GOTO 1000
      END SELECT
      SELECT CASE (kflagcovxyo)
      CASE(1)
         jpssize=jpx
         jpitpsize=1
      CASE(2)
         jpssize=jpyend
         jpitpsize=1
      CASE(3)
         jpssize=jpoend
         jpitpsize=jpitpend
      CASE DEFAULT
         GOTO 1000
      END SELECT
      SELECT CASE (kflaggloloc)
      CASE(0)
         jpzsize=1
      CASE (2,3,4)
         jpzsize=jpz
      CASE DEFAULT
         GOTO 1000
      END SELECT
! --- allocation covrr
      allocate ( covrr(1:jprsize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      covrr(:,:) = FREAL(0.0)
! --- allocation covrrz
      allocate ( covrrz(1:jprsize,1:jprsize,0:jpzsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      covrrz(:,:,0:) = FREAL(0.0)
! --- allocation matUrrz
      allocate ( matUrrz(1:jprsize,1:jprsize,0:jpzsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      matUrrz(:,:,0:) = FREAL(0.0)
! --- allocation lambda
      allocate ( lambda(1:jprsize,0:jpzsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lambda(:,0:) = FREAL(0.0)
! --- allocation amplitude
      allocate ( amplitude(1:jprsize,0:jpzsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      amplitude(:,0:) = FREAL(0.0)
! --- allocation basesr
      allocate ( basesr(1:jpssize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basesr(:,:) = FREAL(0.0)
! --- allocation vectsmean
      allocate ( vectsmean(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectsmean(:) = FREAL(0.0)
! --- allocation vectsect
      allocate ( vectsect(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectsect(:) = FREAL(0.0)
      IF ((lweight).OR.(lweightz)) THEN
! --- allocation vectsweight
         allocate ( vectsweight(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectsweight(:) = FREAL(0.0)
! --- allocation vectsweightz
         allocate ( vectsweightz(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectsweightz(:) = FREAL(0.0)
      ENDIF
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine algobase &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
! -0.1- reading or define config.obs :
! ------------------------------------
!
      IF (kflagcovxyo.EQ.3) THEN
!
! --- allocation poscoefobs
         allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
         flagcfg=3
         CALL readcfgobs (kconfigo,flagcfg, &
     &        kposcoefobs=poscoefobs(:,:))
!
      ENDIF
!
      IF ((lPaout).AND.((kflagcovxyo.EQ.3).AND. &
     &     ((kflagbicovoz.EQ.3).OR.(kflaganlxyo.EQ.3)))) THEN
!
! --- allocation gridijkobs
         allocate ( gridijkobs(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
         flagcfg=2
         CALL readcfgobs (kconfigo,flagcfg, &
     &        kgridijkobs=gridijkobs(:))
! --- allocation vectorms
         allocate ( vectorms(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectorms(:)=FREAL(0.0)
         flagcfg=1
         CALL readcfgobs (kconfigo,flagcfg, &
     &        kvectorms=vectorms(:))
!
      ENDIF
      IF (present(kinbasrefxyoz)) THEN
         serie=0
         numjr=1
         CALL fildirbas (fnamein,kinbasrefxyoz,jprbasrefin,numjr,serie)
         CALL readinfobas(kinbasrefxyoz,valbasref)
         SELECT CASE (valbasref)
         CASE (-1)
            jpdsize=jprbasrefin
            jdbasdeb=1
            jdbasfin=jpdsize
         CASE (0)
            jpdsize=jprbasrefin
            jdbasdeb=1
            jdbasfin=jpdsize
         CASE (1:)
            jpdsize=jprbasrefin
            jdbasdeb=2
            jdbasfin=jprbasrefin
         END SELECT
! --- allocation basesd
         allocate ( basesd(1:jpssize,1:jpdsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         basesd(:,:) = FREAL(0.0)
! --- allocation covrd
         allocate ( covrd(1:jprsize,1:jpdsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         covrd(:,:) = FREAL(0.0)
! --- allocation covrdz
         allocate ( covrdz(1:jprsize,1:jpdsize,0:jpzsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         covrdz(:,:,0:) = FREAL(0.0)
      ENDIF
!
      serie=0
      numjr=1
      CALL fildirbas (fnamein,kinbasxyo,jprbasin,numjr,serie)
      CALL fildirbas (fnameout,koutbasxyo,jprbasout,numjr,serie)
      IF (jprbasin.GT.jprend) GOTO 1000
      IF (jprbasout.GT.jprend) GOTO 1000
      CALL readinfobas(kinbasxyo,valbase)
!
!      IF (valbase.GE.1) THEN
!         jrbasdeb = 2
!         jrbasfin = jprbasin
!         all=.TRUE.
!         titre='INITIALES (ti-1)'
!         IF (nprint.GE.1) CALL printrapport 
!     $        (matUrrz(:,jrbasdeb:jrbasfin,1),
!     $        lambda(jrbasdeb:jrbasfin,1),titre,jrbasdeb,all)
!      ENDIF
!
! -0.3- define parameter :
! ------------------------
!
      SELECT CASE (valbase)
      CASE(1:3)
! --- the new orthogonalisations
         jrbasdeb = 2
         jrbasfin = jprsize
         nvpnull = 0
         nvectindep = 1
      CASE(-1:0)
! --- the orthogonalisations
         jrbasdeb = 1
         jrbasfin = jprsize
         nvpnull = 1
         nvectindep = jprsize - nvpnull
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      IF (lweightoe) THEN
         serie=0
         numjr=1
         CALL fildirbas (fnamein,argoecorrel,jpqsize,numjr,serie)
         jqbasdeb = 2
         jqbasfin = jpqsize
! --- allocation matBqr
         allocate ( matBqr(1:jpqsize,1:jprsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         matBqr(:,:) = FREAL(0.0)
      ENDIF
!
      lectinfo = .FALSE.
      covrr(:,:) = FREAL(0.0)
      IF (present(kinbasrefxyoz)) covrd(:,:) = FREAL(0.0)
      DO jnxyo=1,limjpnxyo(kflagcovxyo)
         lmodprintjnxyo=(MOD(jnxyo-1,(limjpnxyo(kflagcovxyo)/3+1)).EQ.0)
         IF (lmodprintjnxyo) &
     &        print *,'Memory part number : ', &
     &           jnxyo,'/',limjpnxyo(kflagcovxyo)
         IF (lmodprintjnxyo) print *, &
     &        '--> Loading of objects ...'
!
! -1.1- Read weight vector :
! ---------------------------
!
         IF (largweight) THEN
!
            IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE : ', &
     &              'loading weight vector (W)'
            ENDIF
            IF (kflagcovxyo.EQ.3) THEN
               IF ((validextvar(argweight)) &
     &              .OR.(validextdta(argweight)) &
     &              .OR.(validextobs(argweight))) THEN
                  CALL readxyo(argweight,vectsweight(:), &
     &                 jnxyo,lectinfo,kflagcovxyo,poscoefobs(:,:))   
               ELSE 
                  GOTO 1000
               ENDIF   
            ELSE
               IF ((validextvar(argweight)) &
     &              .OR.(validextdta(argweight))) THEN
                  CALL readxyo(argweight,vectsweight(:), &
     &                 jnxyo,lectinfo,kflagcovxyo)   
               ELSE 
                  GOTO 1000
               ENDIF   
            ENDIF
            vectsweight(:) = ABS(vectsweight(:))
            IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
               titre='A few elements of sqrt(W)'
               variable='sqrt(W(*))'
               CALL printtab_r (vectsweight(1:MIN0(jprsize,5)), &
     &              titre,variable,1) 
            ENDIF
!
         ENDIF
!
! -1.2- Read (or fill) and inverse observation error std vector :
! ---------------------------------------------------------------
!
         IF (lweighto) THEN
!
            IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE : ', &
     &              'filling observation error std vector'
            ENDIF
            IF (largoestd) THEN
               IF (kflagcovxyo.EQ.3) THEN
                  IF ((validextvar(argoestd)) &
     &                 .OR.(validextdta(argoestd)) &
     &                 .OR.(validextobs(argoestd))) THEN
                     CALL readxyo(argoestd,vectsweightz(:), &
     &                    jnxyo,lectinfo,kflagcovxyo,poscoefobs(:,:))   
                  ELSE 
                     GOTO 1000
                  ENDIF   
               ELSE
                  IF ((validextvar(argoestd)) &
     &                 .OR.(validextdta(argoestd))) THEN
                     CALL readxyo(argoestd,vectsweightz(:), &
     &                    jnxyo,lectinfo,kflagcovxyo)   
                  ELSE 
                     GOTO 1000
                  ENDIF   
               ENDIF
            ELSE
               CALL mkyorms (vectsweightz(:),kflagcovxyo)
            ENDIF
            IF (ANY(vectsweightz(:).EQ.0)) GOTO 1000
            vectsweightz(:) = FREAL(1.0)/(ABS(vectsweightz(:)))
            IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
               titre='A few elements of sqrt(Woestd)'
               variable='sqrt(Woestd(*))'
               CALL printtab_r (vectsweightz(1:MIN0(jprsize,5)), &
     &              titre,variable,1) 
            ENDIF
!
            IF (largweight) THEN
               vectsweight(:) = vectsweight(:)*vectsweightz(:)
            ELSE
               vectsweight(:) = vectsweightz(:)
            ENDIF
            IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
               titre='A few elements of final sqrt(W)'
               variable='final sqrt(W(*))'
               CALL printtab_r (vectsweight(1:MIN0(jprsize,5)), &
     &              titre,variable,1) 
            ENDIF
!
!
         ENDIF
!
! -1.3- Compute the weight correction to the observation correlation matrix
! -------------------------------------------------------------------------
!
         IF (lweightoe) THEN
            IF (kflagcovxyo.EQ.3) THEN
!
! --- allocation vecty
               allocate ( vecty(0:jpysize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               vecty(:) = FREAL(0.0)
! --- allocation vectynbpoint
               allocate ( vectynbpoint(0:jpysize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               vectynbpoint(:) = 0
!
               coeflimite=FREAL(0.75)
               CALL mkhotoy(vectsweight(:),vecty(:), &
     &              poscoefobs(:,:),spvaldta,kcoeflimite=coeflimite, &
     &              kvectynbpoint=vectynbpoint(:))
               vecty(:) = FREAL(vectynbpoint(:))
               CALL mkhytoo(vecty(1:jpysize),vectsweightz(1:jpssize), &
     &              poscoefobs(:,:))
!     
               vectsweight(:) = vectsweight(:)  &
     &              / MAX(FREAL(1.0),SQRT(vectsweightz(:)))
!     
               vectsweightz(:) = FREAL(0.0)
               IF (allocated(vecty)) deallocate(vecty)
               IF (allocated(vectynbpoint)) deallocate(vectynbpoint)
!
            ENDIF
         ENDIF
!
! -3.1- Read necessary part of inbasref for covariance computation :
! ------------------------------------------------------------------
!
         IF (present(kinbasrefxyoz).AND.(kflagincovxyoz.NE.4)) THEN
!
            IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE : ', &
     &              'loading part of D for covariance computation'
            ENDIF
            IF (kflagcovxyo.EQ.3) THEN
! --- load base.o
               CALL readbas(kinbasrefxyoz,basesd(:,:),jnxyo, &
     &              jdbasdeb,jdbasfin,lectinfo, &
     &              kflagcovxyo,poscoefobs(:,:))
            ELSE
! --- load base.xy
               CALL readbas(kinbasrefxyoz,basesd(:,:),jnxyo, &
     &              jdbasdeb,jdbasfin,lectinfo,kflagcovxyo)
            ENDIF
            IF (valbasref.EQ.0) THEN
               CALL meanect(vectsmean(:), &
     &           vectsect(:),basesd(:,jdbasdeb:jdbasfin))
               DO jd=jdbasdeb,jdbasfin
                  basesd(:,jd) = basesd(:,jd) - vectsmean(:)
               ENDDO
            ENDIF
!
         ENDIF
!
! -3.2- Read necessary part of A for covariance computation :
! ------------------------------------------------------
!
         IF (kflagincovxyoz.NE.4) THEN
!
            IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE : ', &
     &              'loading part of A for covariance computation'
            ENDIF
            IF (kflagcovxyo.EQ.3) THEN
! --- load base.o
               CALL readbas(kincovxyoz,basesr(:,:),jnxyo, &
     &              jrbasdeb,jrbasfin,lectinfo, &
     &              kflagcovxyo,poscoefobs(:,:))
            ELSE
! --- load base.xy
               CALL readbas(kincovxyoz,basesr(:,:),jnxyo, &
     &              jrbasdeb,jrbasfin,lectinfo,kflagcovxyo)
            ENDIF
            IF (valbase.LE.0) CALL meanect(vectsmean(:), &
     &           vectsect(:),basesr(:,jrbasdeb:jrbasfin))
            IF (valbase.EQ.0) THEN
               DO jr=jrbasdeb,jrbasfin
                  basesr(:,jr) = basesr(:,jr) - vectsmean(:)
               ENDDO
            ENDIF
!
         ENDIF
!
! -4.- Prepare arrays for the management of the influence bubbles
! ----------------------------------------------------------------
!
         fixjpbub = jpfixjpz*dtaend
         limjpnz(:) = 1
         IF (kflaggloloc.NE.0) THEN
            CALL evalhdrzon(kinzon,zon_jpi,zon_jpj, &
     &           zon_jpk,zon_jpt,jpbubend(1),jpz)
            IF (jpbubend(1).GT.fixjpbub) THEN
               jpnbub(1) = fixjpbub
               limjpnz(1) = 1 + jpzsize * dtaend / fixjpbub
            ELSE
               jpnbub(1) = jpbubend(1)
               limjpnz(1) = 1
            ENDIF
!
            IF (kflagincovxyoz.EQ.4) THEN
               CALL readhdrzbas(kincovxyoz, &
     &              zon_jpi,zon_jpj,zon_jpk,zon_jpt, &
     &              jpbubend(2),jpz,jrbasdeb,jrbasfin)
               IF (jpbubend(2).GT.fixjpbub) THEN
                  jpnbub(2) = fixjpbub
                  limjpnz(2) = 1 + jpzsize * dtaend / fixjpbub
               ELSE
                  jpnbub(2) = jpbubend(2)
                  limjpnz(2) = 1
               ENDIF
            ENDIF
!
            IF (lweightoe) THEN
               CALL readhdrzbas(argoecorrel, &
     &              zon_jpi,zon_jpj,zon_jpk,zon_jpt, &
     &              jpbubend(3),jpz,jqbasdeb,jqbasfin)
               IF (jpbubend(3).GT.fixjpbub) THEN
                  jpnbub(3) = fixjpbub
                  limjpnz(3) = 1 + jpzsize * dtaend / fixjpbub
               ELSE
                  jpnbub(3) = jpbubend(3)
                  limjpnz(3) = 1
               ENDIF
            ENDIF
            IF (present(kinbasrefxyoz).AND.(kflagincovxyoz.NE.4)) THEN
               CALL readhdrzbas(kinbasrefxyoz, &
     &              zon_jpi,zon_jpj,zon_jpk,zon_jpt, &
     &              jpbubend(4),jpz,jdbasdeb,jdbasfin)
               IF (jpbubend(4).GT.fixjpbub) THEN
                  jpnbub(4) = fixjpbub
                  limjpnz(4) = 1 + jpzsize * dtaend / fixjpbub
               ELSE
                  jpnbub(4) = jpbubend(4)
                  limjpnz(4) = 1
               ENDIF
            ENDIF
!
! --- allocation vectptbub
            allocate ( vectptbub(1:MAXVAL(jpnbub)), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vectptbub(:) = 0
!
! --- allocation vectybub
            allocate ( vectybub(0:jpysize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vectybub(0:) = FREAL(0.0)
!
! --- allocation bubblk1
            allocate ( bubblk1(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &           1:jpnbub(1),1:1), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            bubblk1(:,:,:,:,:,:) = FREAL(0.0)
! --- allocation bubblk2
            IF (kflaggloloc.EQ.4) THEN
               allocate ( bubblk2(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &              1:jpnbub(2),1:jprsize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               bubblk2(:,:,:,:,:,:) = FREAL(0.0)
            ENDIF
! --- allocation bubblk3
            IF (lweightoe) THEN
               allocate ( bubblk3(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &              1:jpnbub(3),1:jpqsize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               bubblk3(:,:,:,:,:,:) = FREAL(0.0)
            ENDIF
! --- allocation bubblk4
            IF (present(kinbasrefxyoz).AND.(kflagincovxyoz.NE.4)) THEN
               allocate ( bubblk4(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &              1:jpnbub(4),1:jpdsize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               bubblk4(:,:,:,:,:,:) = FREAL(0.0)
            ENDIF
!
! --- allocation ptbub
            allocate ( ptbub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &           1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptbub(:,:,:,:,:) = 0
!
! --- allocation zone pointers
            allocate ( pt1dtalon(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt1dtalon(:,:) = 0
            allocate ( pt1dtalat(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt1dtalat(:,:) = 0
            allocate ( pt1dtadepth(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt1dtadepth(:,:) = 0
            allocate ( pt1dtatime(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt1dtatime(:,:) = 0
!
            allocate ( pt1bublon(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt1bublon(:,:) = 0
            allocate ( pt1bublat(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt1bublat(:,:) = 0
            allocate ( pt1bubdepth(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt1bubdepth(:,:) = 0
            allocate ( pt1bubtime(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt1bubtime(:,:) = 0
!     
            allocate ( pt1bubidx(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt1bubidx(:,:) = 0
!
            IF (kflagincovxyoz.EQ.4) THEN
               allocate ( pt2bubidx(1:jpzsize,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               pt2bubidx(:,:) = 0
            ENDIF
!
            IF (lweightoe) THEN
               allocate ( pt3bubidx(1:jpzsize,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               pt3bubidx(:,:) = 0
            ENDIF
!
            IF (present(kinbasrefxyoz).AND.(kflagincovxyoz.NE.4)) THEN
               allocate ( pt4bubidx(1:jpzsize,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               pt4bubidx(:,:) = 0
            ENDIF
!     
            IF (kflagincovxyoz.EQ.4) THEN
               CALL readptzbas (kincovxyoz,pt2bubidx, &
     &              pt1dtalon,pt1dtalat,pt1dtadepth,pt1dtatime, &
     &              pt1bublon,pt1bublat,pt1bubdepth, &
     &              pt1bubtime,jrbasdeb,jrbasfin)
            ENDIF
!     
            IF (lweightoe) THEN
               CALL readptzbas (argoecorrel,pt3bubidx, &
     &              pt1dtalon,pt1dtalat,pt1dtadepth,pt1dtatime, &
     &              pt1bublon,pt1bublat,pt1bubdepth, &
     &              pt1bubtime,jqbasdeb,jqbasfin)
            ENDIF
!
            IF (present(kinbasrefxyoz).AND.(kflagincovxyoz.NE.4)) THEN
               CALL readptzbas (kinbasrefxyoz,pt4bubidx, &
     &              pt1dtalon,pt1dtalat,pt1dtadepth,pt1dtatime, &
     &              pt1bublon,pt1bublat,pt1bubdepth, &
     &              pt1bubtime,jdbasdeb,jdbasfin)
            ENDIF
!
            CALL readptzon (kinzon,pt1bubidx,pt1dtalon,pt1dtalat, &
     &           pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth, &
     &           pt1bubtime)
!
         ENDIF
!
! -5.- Loop on the blocks of subdomains to compute local BASE covariance
! -------------------------------------------------------------------------
!
         IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : beginning loop', &
     &           'over subdomains for covariance computation'
            WRITE(numout,*) '              n = ',jpzsize
         ENDIF
         IF (lmodprintjnxyo) print *, &
     &        '--> Computation of (cov) ...'
         covrrz(:,:,0:) = FREAL(0.0)
         IF (present(kinbasrefxyoz)) covrdz(:,:,0:) = FREAL(0.0)
         DO jnz = 1,MAXVAL(limjpnz)
!
         jz1 = 1
         jz2 = jpzsize
!
         IF (kflaggloloc.NE.0) THEN
!
            IF (ANY(jpnbub(:).NE.jpbubend(:))) THEN
               jz1 = (jnz - 1) * fixjpbub / dtaend + 1
               jz2 = MIN(jnz*fixjpbub/dtaend,jpzsize)
               jpnbubsize = ( jz2 - jz1 + 1 ) * dtaend
            ENDIF
!
            zlectinfo = .FALSE.
            IF (jnz.EQ.1) zlectinfo = .TRUE.
!
            IF (jpnbub(1).EQ.jpbubend(1)) THEN
               IF (jnz.EQ.1) THEN
                  vectptbub(1:jpnbub(1)) =  &
     &                 (/ (jbub, jbub = 1,jpnbub(1)) /)
                  CALL readnbubzon(kinzon,vectptbub(1:jpnbub(1)), &
     &                 bubblk1(:,:,:,:,1:jpnbub(1),1),zlectinfo)
               ENDIF
            ELSE
               vectptbub(1:jpnbubsize) = (/ ((pt1bubidx(jz,jdta), &
     &              jdta=1,dtaend),jz=jz1,jz2) /)
               CALL readnbubzon(kinzon,vectptbub(1:jpnbubsize), &
     &              bubblk1(:,:,:,:,1:jpnbubsize,1),zlectinfo)
            ENDIF
!
            IF (kflagincovxyoz.EQ.4) THEN
               IF (jpnbub(2).EQ.jpbubend(2)) THEN
                  IF (jnz.EQ.1) THEN
                     vectptbub(1:jpnbub(2)) =  &
     &                    (/ (jbub, jbub = 1,jpnbub(2)) /)
                     CALL readnbubzbas(kincovxyoz, &
     &                    bubblk2(:,:,:,:,1:jpnbub(2),:), &
     &                    vectptbub(1:jpnbub(2)), &
     &                    jrbasdeb,jrbasfin,zlectinfo)
                  ENDIF
               ELSE
                  vectptbub(1:jpnbubsize) = (/ ((pt2bubidx(jz,jdta), &
     &                 jdta=1,dtaend),jz=jz1,jz2) /)
                  CALL readnbubzbas(kincovxyoz, &
     &                 bubblk2(:,:,:,:,1:jpnbubsize,:), &
     &                 vectptbub(1:jpnbubsize), &
     &                 jrbasdeb,jrbasfin,zlectinfo)
               ENDIF
            ENDIF
!
            IF (lweightoe) THEN
               IF (jpnbub(3).EQ.jpbubend(3)) THEN
                  IF (jnz.EQ.1) THEN
                     vectptbub(1:jpnbub(3)) =  &
     &                    (/ (jbub, jbub = 1,jpnbub(3)) /)
                     CALL readnbubzbas(argoecorrel, &
     &                    bubblk3(:,:,:,:,1:jpnbub(3),:), &
     &                    vectptbub(1:jpnbub(3)), &
     &                    jqbasdeb,jqbasfin,zlectinfo)
                  ENDIF
               ELSE
                  vectptbub(1:jpnbubsize) = (/ ((pt3bubidx(jz,jdta), &
     &                 jdta=1,dtaend),jz=jz1,jz2) /)
                  CALL readnbubzbas(argoecorrel, &
     &                 bubblk3(:,:,:,:,1:jpnbubsize,:), &
     &                 vectptbub(1:jpnbubsize), &
     &                 jqbasdeb,jqbasfin,zlectinfo)
               ENDIF
            ENDIF
!
            IF (present(kinbasrefxyoz).AND.(kflagincovxyoz.NE.4)) THEN
               IF (jpnbub(4).EQ.jpbubend(4)) THEN
                  IF (jnz.EQ.1) THEN
                     vectptbub(1:jpnbub(4)) =  &
     &                    (/ (jbub, jbub = 1,jpnbub(4)) /)
                     CALL readnbubzbas(kinbasrefxyoz, &
     &                    bubblk4(:,:,:,:,1:jpnbub(4),:), &
     &                    vectptbub(1:jpnbub(4)), &
     &                    jdbasdeb,jdbasfin,zlectinfo)
                  ENDIF
               ELSE
                  vectptbub(1:jpnbubsize) = (/ ((pt4bubidx(jz,jdta), &
     &                 jdta=1,dtaend),jz=jz1,jz2) /)
                  CALL readnbubzbas(kinbasrefxyoz, &
     &                 bubblk4(:,:,:,:,1:jpnbubsize,:), &
     &                 vectptbub(1:jpnbubsize), &
     &                 jdbasdeb,jdbasfin,zlectinfo)
               ENDIF
            ENDIF
         ENDIF
!
! -5.1- Loop on the subdomains
! ----------------------------
! define the obs influence zone for current subdomain :
! -----------------------------------------------------
!
         DO jz=jz1,jz2
!
! -5.1- Define or read obs influence zone for current subdomain :
! ---------------------------------------------------------------
!
            IF (kflaggloloc.NE.0) THEN
!
               vectybub(0:jpysize) = FREAL(0.0)
               DO jdta=1,dtaend
                  IF (jpnbub(1).EQ.jpbubend(1)) THEN
                     jbub = pt1bubidx(jz,jdta)
                  ELSE
                     jbub = (jz-jz1) * dtaend + jdta
                  ENDIF
!
                  CALL mkptbub (ptbub(:,:,:,:,jdta),jdta, &
     &                     pt1dtalon(jz,jdta),pt1dtalat(jz,jdta), &
     &                     pt1dtadepth(jz,jdta),pt1dtatime(jz,jdta), &
     &                     pt1bublon(jz,jdta),pt1bublat(jz,jdta), &
     &                     pt1bubdepth(jz,jdta),pt1bubtime(jz,jdta))
!
                  DO jt=1,zon_jpt
                  DO jk=1,zon_jpk
                  DO jj=1,zon_jpj
                     vectybub(ptbub(:,jj,jk,jt,jdta)) =  &
     &                    bubblk1(:,jj,jk,jt,jbub,1)
                  ENDDO
                  ENDDO
                  ENDDO
               ENDDO
!
               IF (largconnect) CALL mkconnecty(vectybub(0:jpysize))
!
               IF (kflagcovxyo.EQ.3) THEN
                  CALL mkhytoo(vectybub(1:jpysize), &
     &                 vectsweightz(1:jpssize),poscoefobs(:,:))
               ELSE
                  vectsweightz(:)=vectybub(1:jpysize)
               ENDIF
!                 
               jpu = COUNT(vectsweightz(:) .NE. FREAL(0.0))
! --- allocation ptu
               allocate ( ptu(1:jpu), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               ptu(:) = PACK( (/ (js, js=1,jpssize) /) , &
     &              vectsweightz(:) .NE. FREAL(0.0) )
! --- allocation baseur
               allocate ( baseur(1:jpu,1:jprsize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               baseur(:,:) = FREAL(0.0)
! --- allocation baseur
               IF (present(kinbasrefxyoz)) THEN
                  allocate ( baseud(1:jpu,1:jpdsize), stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  baseud(:,:) = FREAL(0.0)
               ENDIF
! --- allocation vectuweight
               allocate ( vectuweight(1:jpu), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               vectuweight(:) = FREAL(0.0)
! --- allocation vectumean
               allocate ( vectumean(1:jpu), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               vectumean(:) = FREAL(0.0)
! --- allocation vectuect
               allocate ( vectuect(1:jpu), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               vectuect(:) = FREAL(0.0)
! --- allocation poscoefuobs
               IF ((lweightoe).OR.(kflagincovxyoz.EQ.4)) THEN
                  allocate ( poscoefuobs(1:jpu,1:jpitpsize),  &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  poscoefuobs(1:jpu,1:jpitpsize) =  &
     &                 poscoefobs(ptu(1:jpu),1:jpitpsize)
               ENDIF
! --- allocation baseuq
               IF (lweightoe) THEN
                  allocate ( baseuq(1:jpu,1:jpqsize), stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  baseuq(:,:) = FREAL(0.0)
               ENDIF
            ENDIF
!
            lmodprintjnz= &
     &           (lmodprintjnxyo.AND.(MOD(jz-1,(jpzsize/5+1)).EQ.0))
            IF (lmodprintjnz) print *,  &
     &           'Zone number : ',jz,'/',jpzsize,'( Size =',jpu,')'
!
! -5.2.1- Read necessary part of inbasref for covariance computation (case kflagincovxyoz.EQ.4) :
! -----------------------------------------------------------------------------------------------
!
            IF (present(kinbasrefxyoz)) THEN
            SELECT CASE (kflagincovxyoz)
            CASE (1)
! --- nothing
            CASE (2,3)
               baseud(1:jpu,1:jpdsize)=basesd(ptu(1:jpu),1:jpdsize)
            CASE (4)
!
               IF ((nprint.GE.2).AND.(jnxyo.EQ.1).AND.(jz.EQ.1)) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOBASE : ', &
     &                 'loading part of inbasref for covariance ', &
     &                 'computation (case kflagincovxyoz.EQ.4)'
               ENDIF
! --- load base.yo
               DO jd = jdbasdeb,jsbasfin
!
                  vectybub(0:jpysize) = FREAL(0.0)
                  DO jdta=1,dtaend
                     IF (jpnbub(4).EQ.jpbubend(4)) THEN
                        jbub = pt4bubidx(jz,jdta)
                     ELSE
                        jbub = (jz-jz1) * dtaend + jdta
                     ENDIF
!
                     DO jt=1,zon_jpt
                     DO jk=1,zon_jpk
                     DO jj=1,zon_jpj
                        vectybub(ptbub(:,jj,jk,jt,jdta)) =  &
     &                       bubblk4(:,jj,jk,jt,jbub,jd)
                     ENDDO
                     ENDDO
                     ENDDO
                  ENDDO
!
                  IF (largconnect) CALL mkconnecty(vectybub(0:jpysize))
!
                  CALL mkhytou(vectybub(1:jpysize),baseud(:,jd), &
     &                 poscoefuobs(:,:))
!
               ENDDO
!                  flagxyo = 3
!                  CALL readzbas(kinbasrefxyoz,baseud(:,:),
!     $                 jz,jdbasdeb,jdbasfin,
!     $                 lectinfo,flagxyo,poscoefuobs)
!
               IF (valbasref.EQ.0) THEN
                  CALL meanect(vectumean(:), &
     &                 vectuect(:),baseud(:,jdbasdeb:jdbasfin))
                  DO jd=jdbasdeb,jdbasfin
                     baseud(:,jd) = baseud(:,jd) - vectumean(:)
                  ENDDO
               ENDIF
            CASE DEFAULT
               GOTO 1000
            END SELECT
!
            ENDIF
!
! -5.2.2- Read necessary part of A for covariance computation (case kflagincovxyoz.EQ.4) :
! ------------------------------------------------------------------------------
!
            SELECT CASE (kflagincovxyoz)
            CASE (1)
! --- nothing
            CASE (2,3)
               baseur(1:jpu,1:jprsize)=basesr(ptu(1:jpu),1:jprsize)
            CASE (4)
!
               IF ((nprint.GE.2).AND.(jnxyo.EQ.1).AND.(jz.EQ.1)) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOBASE : ', &
     &                 'loading part of A for covariance computation', &
     &                 ' (case kflagincovxyoz.EQ.4)'
               ENDIF
! --- load base.yo
               DO jr = jrbasdeb,jrbasfin
!
                  vectybub(0:jpysize) = FREAL(0.0)
                  DO jdta=1,dtaend
                     IF (jpnbub(2).EQ.jpbubend(2)) THEN
                        jbub = pt2bubidx(jz,jdta)
                     ELSE
                        jbub = (jz-jz1) * dtaend + jdta
                     ENDIF
!
                     DO jt=1,zon_jpt
                     DO jk=1,zon_jpk
                     DO jj=1,zon_jpj
                        vectybub(ptbub(:,jj,jk,jt,jdta)) =  &
     &                       bubblk2(:,jj,jk,jt,jbub,jr)
                     ENDDO
                     ENDDO
                     ENDDO
                  ENDDO
!
                  IF (largconnect) CALL mkconnecty(vectybub(0:jpysize))
!
                  CALL mkhytou(vectybub(1:jpysize),baseur(:,jr), &
     &                 poscoefuobs(:,:))
!
               ENDDO
!               flagxyo = 3
!               CALL readzbas(kincovxyoz,baseur(:,:),
!     $              jz,jrbasdeb,jrbasfin,
!     $              lectinfo,flagxyo,poscoefuobs)
!
               IF (valbase.EQ.0) CALL meanect(vectumean(:), &
     &              vectuect(:),baseur(:,jrbasdeb:jrbasfin))
               IF (valbase.EQ.0) THEN
                  DO jr=jrbasdeb,jrbasfin
                     baseur(:,jr) = baseur(:,jr) - vectumean(:)
                  ENDDO
               ENDIF
            CASE DEFAULT
               GOTO 1000
            END SELECT
!
! -4.5- Load the inverse observation error correlation matrix 
! -----------------------------------------------------------
!
            IF (lweightoe) THEN
               IF (jz.EQ.1) THEN
                  IF (kflaggloloc.EQ.0) GOTO 1000
                  IF (nprint.GE.1) THEN
                     WRITE(numout,*)
                     WRITE(numout,*) 'ALGOBASE : ', &
     &                    'loading observation error correlation matrix'
                  ENDIF
               ENDIF
! --- load base.yo
               DO jq = jqbasdeb,jqbasfin
                  vectybub(0:jpysize) = FREAL(0.0)
!
                  DO jdta=1,dtaend
                     IF (jpnbub(3).EQ.jpbubend(3)) THEN
                        jbub = pt3bubidx(jz,jdta)
                     ELSE
                        jbub = (jz-jz1) * dtaend + jdta
                     ENDIF
!
                     DO jt=1,zon_jpt
                     DO jk=1,zon_jpk
                     DO jj=1,zon_jpj
                        vectybub(ptbub(:,jj,jk,jt,jdta)) =  &
     &                       bubblk3(:,jj,jk,jt,jbub,jq)
                     ENDDO
                     ENDDO
                     ENDDO
                  ENDDO
!
                  IF (largconnect) CALL mkconnecty(vectybub(0:jpysize))
!
                  CALL mkhytou(vectybub(1:jpysize),baseuq(:,jq), &
     &                 poscoefuobs(:,:))
!
               ENDDO
            ENDIF
!
! -5.3- Weight using the influence zone shape function
! -----------------------------------------------------------
!
            IF (kflaggloloc.EQ.0) THEN
               IF (lweight) THEN
                  vectsweightz(:) = vectsweight(:)
               ENDIF
            ELSE
               IF (lweight) THEN
                  vectuweight(1:jpu) = vectsweightz(ptu(1:jpu))  &
     &                 * vectsweight(ptu(1:jpu))
               ELSE
                  vectuweight(1:jpu) = vectsweightz(ptu(1:jpu))
               ENDIF
!
            ENDIF
!
! -5.4- Covariance computing :
! -------------------------
!
            IF ((lweightz).OR.(lweight)) THEN
               IF (kflaggloloc.EQ.0) THEN
                  DO jrbas1=jrbasdeb,jrbasfin
                  DO jrbas2=jrbas1,jrbasfin
                     covrrz(jrbas2,jrbas1,jz)= &
     &                    DOT_PRODUCT(basesr(:,jrbas1)*vectsweightz(:), &
     &                    basesr(:,jrbas2)*vectsweightz(:))
                  ENDDO
                  ENDDO
               ELSE
                  IF (lweightoe) THEN
                     DO jrbas2=jrbasdeb,jrbasfin
                     DO jqbas1=jqbasdeb,jqbasfin
                        jqmat1 = jqbas1 - 1
                        matBqr(jqmat1,jrbas2) =  &
     &                       DOT_PRODUCT(baseuq(:,jqbas1), &
     &                       baseur(:,jrbas2)*vectuweight(:))
                     ENDDO
                     ENDDO
                     DO jrbas1=jrbasdeb,jrbasfin
                     DO jrbas2=jrbas1,jrbasfin
                        covrrz(jrbas2,jrbas1,jz) =  &
     &                       DOT_PRODUCT( &
     &                       matBqr(jqmatdeb:jqmatfin,jrbas2), &
     &                       matBqr(jqmatdeb:jqmatfin,jrbas1))
                     ENDDO
                     ENDDO
                  ELSE
                     DO jrbas1=jrbasdeb,jrbasfin
                     DO jrbas2=jrbas1,jrbasfin
                        covrrz(jrbas2,jrbas1,jz)= &
     &                       DOT_PRODUCT(baseur(:,jrbas1)*vectuweight(:), &
     &                       baseur(:,jrbas2)*vectuweight(:))
                     ENDDO
                     ENDDO
                  ENDIF
               ENDIF
            ELSE
               IF (kflaggloloc.EQ.0) THEN
                  DO jrbas1=jrbasdeb,jrbasfin
                  DO jrbas2=jrbas1,jrbasfin
                     covrrz(jrbas2,jrbas1,jz)= &
     &                    DOT_PRODUCT(basesr(:,jrbas1), &
     &                    basesr(:,jrbas2))
                  ENDDO
                  ENDDO
               ELSE
                  DO jrbas1=jrbasdeb,jrbasfin
                  DO jrbas2=jrbas1,jrbasfin
                     covrrz(jrbas2,jrbas1,jz)= &
     &                    DOT_PRODUCT(baseur(:,jrbas1), &
     &                    baseur(:,jrbas2))
                  ENDDO
                  ENDDO
               ENDIF
            ENDIF
            DO jrbas=jrbasdeb,jrbasfin
!               covrrz(jrbas:jprsize,jrbas,jz) =
!     $              covrrz(jrbas,jrbas:jrbasfin,jz)
               covrrz(jrbas,jrbas:jrbasfin,jz) = &
     &              covrrz(jrbas:jrbasfin,jrbas,jz)
            ENDDO
!
!
! -5.5- inbasref-covariance computing :
! ----------------------------------
!
            IF (present(kinbasrefxyoz)) THEN
               IF ((lweightz).OR.(lweight)) THEN
                  covrdz(:,:,jz) = FREAL(0.0)
                  IF (kflaggloloc.EQ.0) THEN
                     DO jr=jrbasdeb,jrbasfin
                     DO jd=jdbasdeb,jdbasfin
                        covrdz(jr,jd,jz) = &
     &                       DOT_PRODUCT(basesr(:,jr)*vectsweightz(:), &
     &                       basesd(:,jd)*vectsweightz(:))
                     ENDDO
                     ENDDO
                  ELSE
                     IF (lweightoe) THEN
                        DO jdbas2=jdbasdeb,jdbasfin
                        DO jqbas1=jqbasdeb,jqbasfin
                           jqmat1 = jqbas1 - 1
                           matBqd(jqmat1,jdbas2) =  &
     &                          DOT_PRODUCT(baseuq(:,jqbas1), &
     &                          baseud(:,jdbas2)*vectuweight(:))
                        ENDDO
                        ENDDO
                        DO jd=jdbasdeb,jdbasfin
                        DO jr=jrbasdeb,jrbasfin
                           covrdz(jr,jd,jz) =  &
     &                         DOT_PRODUCT(matBqr(jqmatdeb:jqmatfin,jr), &
     &                          matBqd(jqmatdeb:jqmatfin,jd))
                        ENDDO
                        ENDDO
                     ELSE
                        DO jr=jrbasdeb,jrbasfin
                        DO jd=jdbasdeb,jdbasfin
                           covrdz(jr,jd,jz) = &
     &                          DOT_PRODUCT(baseur(:,jr)*vectuweight(:), &
     &                          baseud(:,jd)*vectuweight(:))
                        ENDDO
                        ENDDO
                     ENDIF               
                  ENDIF           
               ELSE
                  covrdz(:,:,jz) = FREAL(0.0)
                  IF (kflaggloloc.EQ.0) THEN
                     DO jr=jrbasdeb,jrbasfin
                     DO jd=jdbasdeb,jdbasfin
                        covrdz(jr,jd,jz) = &
     &                       DOT_PRODUCT(basesr(:,jr), &
     &                       basesd(:,jd))
                     ENDDO
                     ENDDO
                  ELSE
                     DO jr=jrbasdeb,jrbasfin
                     DO jd=jdbasdeb,jdbasfin
                        covrdz(jr,jd,jz) = &
     &                       DOT_PRODUCT(baseur(:,jr), &
     &                       baseud(:,jd))
                     ENDDO
                     ENDDO
                  ENDIF               
               ENDIF               
            ENDIF
!
! --- desallocate some arrays
            IF (allocated(ptu))  deallocate (ptu)
            IF (allocated(baseur))  deallocate (baseur)
            IF (allocated(vectuweight))  deallocate (vectuweight)
            IF (allocated(vectuect))  deallocate (vectuect)
            IF (allocated(vectumean))  deallocate (vectumean)
            IF (allocated(poscoefuobs))  deallocate (poscoefuobs)
!
         ENDDO
         ENDDO
!
! -6.- next step of covariance computing :
! ----------------------------------------
!
         IF ((kflagcovxyo.EQ.1).AND.(limjpnxyo(kflagcovxyo).NE.1)) THEN
            covrr(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin) =  &
     &           covrrz(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin,1) +  &
     &           covrr(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin)
         ENDIF
         IF (present(kinbasrefxyoz)) THEN
            IF ((kflagcovxyo.EQ.1) &
     &           .AND.(limjpnxyo(kflagcovxyo).NE.1)) THEN
               covrd(jrbasdeb:jrbasfin,jdbasdeb:jdbasfin) =  &
     &              covrdz(jrbasdeb:jrbasfin,jdbasdeb:jdbasfin,1) +  &
     &              covrd(jrbasdeb:jrbasfin,jdbasdeb:jdbasfin)
            ENDIF
         ENDIF
!
      ENDDO
!
      IF (limjpnxyo(kflagcovxyo).NE.1) THEN
         titre = 'Tableau covariance'
         variable = 'covrr(:,:)'
         CALL printtab_rr (covrr(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin), &
     &        titre,variable,jrbasdeb,jrbasfin)
      ENDIF
!
! -7.- End of covariance computing :
! ----------------------------------
!
      IF ((kflagcovxyo.EQ.1).AND.(limjpnxyo(kflaganlxyo).NE.1)) THEN
         covrrz(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin,1) =  &
     &        covrr(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin)
         IF (present(kinbasrefxyoz)) THEN
            covrdz(jrbasdeb:jrbasfin,jdbasdeb:jdbasfin,1) =  &
     &           covrd(jrbasdeb:jrbasfin,jdbasdeb:jdbasfin)
         ENDIF
      ENDIF
      covrrz(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin,1:) =  &
     &     covrrz(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin,1:)  &
     &     / nvectindep
!
      titre = 'Tableau covariance'
      variable = 'covrrz(:,:,1)'
      CALL printtab_rr (covrrz(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin,1), &
     &     titre,variable,jrbasdeb,jrbasfin)
!
! -8.- Loop on subdomains to compute local BASE rs
! ----------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ALGOBASE : ', &
     &        'beginning loop over subdomains for rs computation'
         WRITE(numout,*) '              n = ',jpzsize
      ENDIF
      IF (nprint.GE.1) print *,'--> Computation of (U,lambda) ...'
      DO jz=1,jpzsize
         lmodprintjnz=(MOD(jz-1,(jpzsize/2+1)).EQ.0)
         IF ((nprint.GE.1).AND.((lmodprintjnz).OR.(kflaggloloc.EQ.0))) print *,  &
     &        'Zone number : ',jz,'/',jpzsize
!
         CALL mkrs (covrrz(jrbasdeb:,jrbasdeb:,jz), &
     &        matUrrz(jrbasdeb:,jrbasdeb:,jz), &
     &        lambda(jrbasdeb:,jz),nvpnull)
!
      ENDDO
!
      IF (nprint.GE.1) print *,'--> Print of (U,lambda) ...'
! --- Print control of eigen values - amplitudes - eigen vectors :
      DO jz=1,jpzsize
         lmodprintjnz=(MOD(jz-1,(jpzsize/10+1)).EQ.0)
         IF ((nprint.GE.1).AND.((lmodprintjnz).OR.(kflaggloloc.EQ.0))) THEN
            all=.TRUE.
            titre = &
     &      'ALGOBASE : eigen values,eigen vectors UlambdaU^=(HA)^W(HA)'
            CALL printrapport (matUrrz(:,jrbasdeb:jrbasfin,jz), &
     &           lambda(jrbasdeb:jrbasfin,jz),titre,jrbasdeb,all)
            IF (kflaggloloc.EQ.0) THEN
              CALL printamplitude (matUrrz(:,jrbasdeb:jrbasfin,jz), &
     &                             koutbasxyo)
            ENDIF
         ENDIF
      ENDDO
!
! --- U = U . (1/(jpr-1))
      matUrrz(:,:,1:) = matUrrz(:,:,1:) / SQRT(FREAL(nvectindep))
!
      DO jz=1,jpzsize
         lmodprintjnz=(MOD(jz-1,(jpzsize/10+1)).EQ.0)
         IF ((nprint.GE.1).AND.((lmodprintjnz).OR.(kflaggloloc.EQ.0))) THEN
            titre='ALGOBASE: U eigen basis'
            variable='kmatUrr/SQRT(nvectindep)'
            CALL printtab_rr ( &
     &           matUrrz(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin,jz), &
     &           titre,variable,jrbasdeb,jrbasdeb)
         ENDIF
      ENDDO
!
! -8.1- Loop on subdomains to compute local inbasref cov
! ----------------------------------------------------------
!
      IF (present(kinbasrefxyoz)) THEN
         DO jd=jdbasdeb,jdbasfin
         DO jz=1,jpzsize
            lmodprintjnz=(MOD(jz-1,(jpzsize/5+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprintjnz).OR.(kflaggloloc.EQ.0))) THEN
               titre='ALGOBASE: elements of covrdz=(HA)W(HB)'
               WRITE (variable,'(A,I3,A)') 'covrdz(:,',jd,',jz)'
               CALL printtab_r (covrdz(:,jd,jz), &
     &              titre,variable,1) 
            ENDIF
         ENDDO
         ENDDO
         IF (nprint.GE.2) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : beginning loop', &
     &           ' over subdomains for end regression computation'
            WRITE(numout,*) '              n = ',jpzsize
         ENDIF
         IF (nprint.GE.2) &
     &      print *,'--> Computation of U.lambda~.U^.covrdz ...'
         DO jz=1,jpzsize
            lmodprintjnz=(MOD(jz-1,(jpzsize/2+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprintjnz).OR.(kflaggloloc.EQ.0))) print *,  &
     &           'Zone number : ',jz,'/',jpzsize
! --- Urr^ covrd
            covrd(:,:) = FREAL(0.0)
            DO jd=jdbasdeb,jdbasfin
            DO jr=2,jrbasfin
               covrd(jr,jd) =  &
     &              DOT_PRODUCT(matUrrz(:,jr,jz), &
     &              covrdz(:,jd,jz))
            ENDDO
            ENDDO
            covrdz(:,:,jz) = covrd(:,:)
         ENDDO
         DO jd=jdbasdeb,jdbasfin
         DO jz=1,jpzsize
            lmodprintjnz=(MOD(jz-1,(jpzsize/5+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprintjnz).OR.(kflaggloloc.EQ.0))) THEN
               titre='ALGOBASE: elements of covrdz=U^.(HA)W(HB)'
               WRITE (variable,'(A,I3,A)') 'covrdz(:,',jd,',jz)'
               CALL printtab_r (covrdz(:,jd,jz), &
     &              titre,variable,1) 
            ENDIF
         ENDDO
         ENDDO
!
         DO jz=1,jpzsize
            lmodprintjnz=(MOD(jz-1,(jpzsize/2+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprintjnz).OR.(kflaggloloc.EQ.0)))  &
     &           print *,'Zone number : ',jz,'/',jpzsize
! --- lambda~ Urr^ covrd
            DO jd=jdbasdeb,jdbasfin
               covrdz(2:jrbasfin,jd,jz) = covrdz(2:jrbasfin,jd,jz) / &
     &              lambda (2:jrbasfin,jz)
            ENDDO
         ENDDO
         DO jd=jdbasdeb,jdbasfin
         DO jz=1,jpzsize
            lmodprintjnz=(MOD(jz-1,(jpzsize/5+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprintjnz).OR.(kflaggloloc.EQ.0))) THEN
               titre='ALGOBASE: elements of covrdz=lambda~ U^.(HA)W(HB)'
               WRITE (variable,'(A,I3,A)') 'covrdz(:,',jd,',jz)'
               CALL printtab_r (covrdz(:,jd,jz), &
     &              titre,variable,1) 
            ENDIF
         ENDDO
         ENDDO
!
         DO jz=1,jpzsize
            lmodprintjnz=(MOD(jz-1,(jpzsize/2+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprintjnz).OR.(kflaggloloc.EQ.0))) print *,  &
     &           'Zone number : ',jz,'/',jpzsize
! --- U lambda~ Urr^ covrd
            covrd(:,:) = FREAL(0.0)
            DO jd=jdbasdeb,jdbasfin
            DO jr=jrbasdeb,jrbasfin
               covrd(jr,jd) =  &
     &              DOT_PRODUCT(matUrrz(jr,2:jrbasfin,jz), &
     &              covrdz(2:jrbasfin,jd,jz))
            ENDDO
            ENDDO
            covrdz(:,:,jz) = covrd(:,:)
         ENDDO
         DO jd=jdbasdeb,jdbasfin
         DO jz=1,jpzsize
            lmodprintjnz=(MOD(jz-1,(jpzsize/5+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprintjnz).OR.(kflaggloloc.EQ.0))) THEN
               titre='ALGOBASE: elements of covrdz=UlambdaU^.(HA)W(HB)'
               WRITE (variable,'(A,I3,A)') 'covrdz(:,',jd,',jz)'
               CALL printtab_r (covrdz(:,jd,jz), &
     &              titre,variable,1) 
            ENDIF
         ENDDO
         ENDDO
!
      ENDIF
!
! -10.- Allocation/desallocation
! ------------------------------
!
! --- desallocate some arrays
      IF (allocated(covrd))  deallocate (covrd)
      IF (allocated(covrrz))  deallocate (covrrz)
      IF (allocated(covrr))  deallocate (covrr)
      IF (allocated(lambda))  deallocate (lambda)
      IF (allocated(amplitude))  deallocate (amplitude)
      IF (allocated(vectsweightz))  deallocate (vectsweightz)
      IF (allocated(vectsweight))  deallocate (vectsweight)
      IF (allocated(basesd))  deallocate (basesd)
!
      IF (allocated(ptbub)) deallocate (ptbub)
      IF (allocated(vectybub)) deallocate (vectybub)
      IF (allocated(bubblk1)) deallocate (bubblk1)
      IF (allocated(bubblk2)) deallocate (bubblk2)
      IF (allocated(bubblk3)) deallocate (bubblk3)
      IF (allocated(matBqr)) deallocate (matBqr)
      IF (allocated(pt3bubidx)) deallocate (pt3bubidx)
         IF (allocated(vectptbub)) deallocate (vectptbub)
!
      IF (kflagbicovoz.NE.4) THEN
         IF (allocated(pt1bubidx)) deallocate (pt1bubidx)
         IF (allocated(pt2bubidx)) deallocate (pt2bubidx)
         IF (allocated(pt1bublon)) deallocate (pt1bublon)
         IF (allocated(pt1bublat)) deallocate (pt1bublat)
         IF (allocated(pt1bubdepth)) deallocate (pt1bubdepth)
         IF (allocated(pt1bubtime)) deallocate (pt1bubtime)
         IF (allocated(pt1dtalon)) deallocate (pt1dtalon)
         IF (allocated(pt1dtalat)) deallocate (pt1dtalat)
         IF (allocated(pt1dtadepth)) deallocate (pt1dtadepth)
         IF (allocated(pt1dtatime)) deallocate (pt1dtatime)
      ENDIF
!
      IF (.NOT.(lPaout)) RETURN
!
      IF (kflagbicovoz.EQ.4) THEN
         IF (allocated(basesr))  deallocate (basesr)
         IF (allocated(vectsmean))  deallocate (vectsmean)
         IF (allocated(vectsect))  deallocate (vectsect)
! --- allocation basesr
         allocate ( basesr(0:jpysize,1:jprsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         basesr(:,:) = FREAL(0.0)
! --- allocation vectsmean
         allocate ( vectsmean(0:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectsmean(:) = FREAL(0.0)
! --- allocation vectsect
         allocate ( vectsect(0:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectsect(:) = FREAL(0.0)
      ENDIF
      IF (kflagbicovoz.EQ.4) THEN
! --- allocation vectsmode
         allocate ( vectsmode(0:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectsmode(:) = FREAL(0.0)
! --- allocation vectsstd
         allocate ( vectsstd(0:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectsstd(:) = FREAL(0.0)
      ENDIF
!
! --- allocate an additional state vector
! --- allocation mat_yr
      allocate ( mat_yr(1:MIN(jpysize,jpxsize),1:jprsize),stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      mat_yr(:,:) = FREAL(0.0)
      IF (kflagbicovoz.NE.0) THEN
! --- allocation vectspart
         allocate ( vectspart(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectspart(:) = FREAL(0.0)
      ENDIF
!
! -11.- bi-cov treatment
! -----------------------
!
      IF ((nprint.GE.1).AND.(kflagbicovoz.NE.0)) print *, &
     &     '--> Computation of (bi-bas) ...'
      jnxyo=1
      lectinfo=.FALSE.
!
      SELECT CASE (kflagbicovoz)
      CASE (0)
! --- Nothing
      CASE (3)
!
! -11.1- Compute and write the bi-initial error covariance matrix (HS) in case outobas
! -------------------------------------------------------------------------------------
         IF (nprint.GE.2) THEN
!               print *,'-11.1- Compute and write the bi-initial error',
!     $              ' covariance matrix (HS) in case outobas'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : ', &
     &           'compute and write the bi-initial', &
     &           ' error covariance matrix (HS) in case outobas'
         ENDIF
!     
         IF (kflaggloloc.EQ.3) THEN
!
! -11.1.1- Read the partition in subdomains
! -------------------------------------------
            IF (nprint.GE.2) THEN
!     print *,'-11.1.1- Read the partition in subdomains for bicov'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE: ', &
     &              'reading the partition in subdomains for bicov'
            ENDIF
! --- allocation vecty
            allocate ( vecty(1:jpysize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vecty(:) = FREAL(0.0)
            flagxyo=2
            lmoyectold=lmoyect
            lmoyect=.FALSE.
            CALL readxyo (kinpartxyo,vecty(:),jnxyo, &
     &           lectinfo,flagxyo)
            lmoyect=lmoyectold
            vectspart(:)=vecty(poscoefobs(:,1)%pos)
            IF (allocated(vecty))  deallocate (vecty)
         ENDIF
!
! -11.1.2- Compute the bi-initial error covariance matrix (HS)
! -----------------------------------------------------------
         IF (nprint.GE.2) THEN
!               print *,'-11.2.1- Compute the bi-initial error',
!     $              ' covariance matrix (HS)'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : ', &
     &           'compute the bi-initial', &
     &           ' error covariance matrix (S)'
         ENDIF
         numjr=0
         serie=1
         CALL fildirbas (vctnamout,koutcovoz,jprbasout,numjr,serie)
         IF (kflaggloloc.EQ.0) THEN
            CALL prodmat_zr_rr_vias (basesr(1:,jrbasdeb:jrbasfin), &
     &           mat_yr(:,jrbasdeb:jrbasfin), &
     &           matUrrz(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin,1), &
     &           (2-(1-nvpnull)),(jprbasout-(1-nvpnull)))
         ELSE
            CALL prodmat_wr_rrz_vias(basesr(1:,jrbasdeb:jrbasfin), &
     &           mat_yr(:,jrbasdeb:jrbasfin), &
     &           matUrrz(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin,0:), &
     &           vectspart(:), &
     &           (2-(1-nvpnull)),(jprbasout-(1-nvpnull)))  
         ENDIF
!
! -11.1.3- Transfer or compute mean file into output bibase directory :
! ---------------------------------------------------------------------
         IF (nprint.GE.2) THEN
!            print *,'-11.1.3- Transfer or compute mean file into',
!     $           ' output bibase directory'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : ', &
     &           'Transfer mean file into output bibase directory'
         ENDIF
         IF (valbase.GE.1) THEN
            numjr=1
            serie=1
            CALL fildirbas (vctnamin,kincovxyoz,jprbasin,numjr,serie)
            WRITE(fnamein,'("./",A,"/",A)')  &
     &           kincovxyoz(1:lenv(kincovxyoz)), &
     &           vctnamin(1:lenv(vctnamin))
            CALL readxyo(fnamein,basesr(1:,1), &
     &           jnxyo,lectinfo,kflagbicovoz,poscoefobs(:,:))
         ELSE
            basesr(1:,1) = vectsmean(1:)
         ENDIF
!     
! -11.1.4- Writing output base directory 
! ---------------------------------------
         IF (nprint.GE.2) THEN
!            print *,'-11.1.4- Writing',
!     $           ' output bibase directory'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : ', &
     &           'writing the bicovariance (HS)'
         ENDIF
         numjr=1
         serie=1
         CALL fildirbas (vctnamin,koutcovoz,jprbasout,numjr,serie)
         jrbasdeb1=1
         CALL writeyobas(koutcovoz,basesr(1:,:), &
     &        jrbasdeb1,jprbasout,kflagbicovoz,vectorms(:), &
     &        gridijkobs(:),poscoefobs(:,:))
         CALL writeinfobas(koutcovoz,kflagbicovoz)
!
! -11.1.5- Compute and write diagonal Pa into output base directory    
! ------------------------------------------------------------------
         IF (nprint.GE.2) THEN
!            print *,'-11.1.5- Compute and write',
!     $           ' diagonal Pa into output bibase directory'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : Compute and write diagonal', &
     &           ' HPH into output bibase directory'
         ENDIF
         IF (valbase.EQ.0) THEN
! --- write ect init of the population => *SQRT(jpr/(jpr-1))
            vectsect(:) = vectsect(:) &
     &           * SQRT(FREAL(jprsize)/FREAL(jprsize-1))
            numjr=-1
            serie=1
            CALL fildirbas (vctnamout,koutcovoz,jprbasout,numjr,serie)
            jext=indext(vctnamout,extobstab,nbextobs)
            IF (extobsunit(jext)) THEN
               WRITE (text,'(A)') extobstab(jext)
            ELSE
               WRITE (text,'("#",A)') extobstab(jext)
            ENDIF
            WRITE(fnameout,'("./",A,"/vctectinit",A)')  &
     &           koutcovoz(1:lenv(koutcovoz)), &
     &           text(1:lenv(text))
            CALL writeobs(fnameout,vectsect(1:),vectorms(:), &
     &           gridijkobs(:),poscoefobs(:,:))
         ENDIF
! --- compute ect end
         basesr(:,1)=FREAL(0.0)
         DO jr=(jrbasdeb+nvpnull),jprbasout
            basesr(1:,1)= basesr(1:,jr)*basesr(1:,jr)+basesr(1:,1)
         ENDDO
         basesr(1:,1)=SQRT(basesr(1:,1))
! --- write ect end
         numjr=0
         serie=1
         CALL fildirbas (vctnamout,koutcovoz,jprbasout,numjr,serie)
         WRITE(fnameout,'("./",A,"/",A)')  &
     &        koutcovoz(1:lenv(koutcovoz)), &
     &        vctnamout(1:lenv(vctnamout))
         CALL writeobs(fnameout,basesr(1:,1),vectorms(:), &
     &        gridijkobs(:),poscoefobs(:,:))
!
! -11.1.5- Algortho :
! --------------------
! --- orthogonalite computation
      IF ((kflaggloloc.EQ.0).AND.(.NOT.(present(kinbasrefxyoz)))) THEN
         print *,'Algortho :  orthogonalite computation...' 
         IF (largweight) THEN
            CALL chkortho(koutcovoz,kflagbicovoz,basesr(1:,:), &
     &           2,jprbasout,largweight, &
     &           kargxyoweight=argweight,kvectsweight=basesr(1:,1))
         ELSE
            CALL chkortho(koutcovoz,kflagbicovoz,basesr(1:,:), &
     &           2,jprbasout,largweight)
         ENDIF
      ENDIF
!
      CASE (4)
!
! -11.2- Compute and write the local (bi-)initial error covariance matrix (HS) in case outzbas
! --------------------------------------------------------------------------------------------
         IF (nprint.GE.2) THEN
!               print *,'-11.2- Compute and write the local (bi-)initial error',
!     $              ' covariance matrix (HS) in case outzbas'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : ', &
     &           'compute and write the local (bi-)initial', &
     &           ' error covariance matrix (HS) in case outzbas'
         ENDIF
!
! -11.2.0- Read necessary part of A for covariance computation :
! ---------------------------------------------------------------
!
         IF (kflagincovxyoz.NE.4) THEN
!
            IF (nprint.GE.2) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE : ', &
     &              'loading part of A for covariance computation'
            ENDIF
            IF (kflagcovxyo.LE.2) THEN
! --- load base.xy
               flagxyo=2
               CALL readbas(kincovxyoz,basesr(1:,:),jnxyo, &
     &              jrbasdeb,jrbasfin,lectinfo,flagxyo)
            ELSE
               GOTO 1000
            ENDIF
            SELECT CASE (valbase)
            CASE (:0)
               CALL meanect(vectsmean(1:), &
     &              vectsect(1:),basesr(1:,jrbasdeb:jrbasfin))
! --- ect init of the population => *SQRT(jpr/(jpr-1))
               vectsect(1:) = vectsect(1:) &
     &              * SQRT(FREAL(jprsize)/FREAL(jprsize-1))
            CASE (1:)
               numjr=1
               serie=1
               CALL fildirbas (vctnamin,kinbasxyo,jprbasin,numjr,serie)
               WRITE(fnamein,'("./",A,"/",A)')  &
     &              kinbasxyo(1:lenv(kinbasxyo)), &
     &              vctnamin(1:lenv(vctnamin))
               CALL readxyo(fnamein,vectsmean(1:), &
     &              jnxyo,lectinfo,flagxyo)
               numjr=0
               serie=1
               CALL fildirbas (vctnamin,kinbasxyo,jprbasin,numjr,serie)
               WRITE(fnamein,'("./",A,"/",A)')  &
     &              kinbasxyo(1:lenv(kinbasxyo)), &
     &              vctnamin(1:lenv(vctnamin))
               CALL readxyo(fnamein,vectsect(1:), &
     &              jnxyo,lectinfo,flagxyo)
            CASE DEFAULT 
               GOTO 1000
            END SELECT
            IF (valbase.EQ.0) THEN
               DO jr=jrbasdeb,jrbasfin
                  basesr(1:,jr) = basesr(1:,jr) - vectsmean(1:)
               ENDDO
            ENDIF
         ENDIF
!
! -11.2.0- Evaluation of the filename :
! -------------------------------------
!
         IF (kflagincovxyoz.EQ.4) THEN
!
            CALL fildirbas (vctnamin,kincovxyoz,jprbasin,numjr,serie)
            allocate ( tabfnin(0:jprbasin), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            tabfnin(:) = ''
            DO numjr=0,jprbasin
               serie=1
               CALL fildirbas (vctnamin,kincovxyoz,jprbasin1,numjr,serie)
               WRITE(tabfnin(numjr),'("./",A,"/",A)')  &
     &              kincovxyoz(1:lenv(kincovxyoz)), &
     &              vctnamin(1:lenv(vctnamin))
            ENDDO
         ENDIF
!
         CALL fildirbas (vctnamout,koutcovoz,jprbasout,numjr,serie)
         allocate ( tabfnout(-1:jprbasout), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         tabfnout(:) = ''
         DO numjr=0,jprbasout
            serie=1
            CALL fildirbas (vctnamout,koutcovoz,jprbasout1,numjr,serie)
            WRITE(tabfnout(numjr),'("./",A,"/",A)')  &
     &           koutcovoz(1:lenv(koutcovoz)), &
     &           vctnamout(1:lenv(vctnamout))
         ENDDO
         numjr=-1
         serie=1
         CALL fildirbas (vctnamout,koutcovoz,jprbasout,numjr,serie)
         jext=indext(vctnamout,extzontab,nbextzon)
         IF (extzonunit(jext)) THEN
            WRITE (text,'(A)') extzontab(jext)
         ELSE
            WRITE (text,'("#",A)') extzontab(jext)
         ENDIF
         WRITE(tabfnout(-1),'("./",A,"/vctectinit",A)')  &
     &        koutcovoz(1:lenv(koutcovoz)), &
     &        text(1:lenv(text))
!
! -11.2.1- Write the info of the local (bi-)initial error covariance matrix (HzxA)
! --------------------------------------------------------------------------------
!
         CALL evalhdrzon (kinzon,zon_jpi,zon_jpj,zon_jpk, &
     &        zon_jpt,jpbub,jpz)
!         allmemjpzweight=(fixjpbub.GE.jpbub)
!
         IF (kflagincovxyoz.EQ.4) THEN
!
            CALL evalhdrzon (tabfnin(1),zon_jpi1,zon_jpj1,zon_jpk1, &
     &        zon_jpt1,jpbub1,jpz1)
            IF (zon_jpi1.NE.zon_jpi) GOTO 102
            IF (zon_jpj1.NE.zon_jpj) GOTO 102
            IF (zon_jpk1.NE.zon_jpk) GOTO 102
            IF (zon_jpt1.NE.zon_jpt) GOTO 102
            IF (jpz1.NE.jpz) GOTO 102
         ELSE
            jpbub1=jpz*dtaend
         ENDIF
         allmemjpzr=(fixjpbub.GE.jpbub1)
!
!         IF (allmemjpzweight) THEN
!            jpbubsize=jpbub
!         ELSE
!            jpbubsize=((fixjpbub/dtaend)+1)*dtaend
!         ENDIF
         
         IF (allmemjpzr) THEN
            jpbubsize1=jpbub1
         ELSE
            jpbubsize1=(fixjpbub/dtaend)*dtaend
         ENDIF
! --- allocation zone pointers
         allocate ( ptbubidx(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptbubidx(:,:) = 0
         allocate ( ptdtalon(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptdtalon(:,:) = 0
         allocate ( ptdtalat(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptdtalat(:,:) = 0
         allocate ( ptdtadepth(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptdtadepth(:,:) = 0
         allocate ( ptdtatime(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptdtatime(:,:) = 0
         allocate ( ptbublon(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptbublon(:,:) = 0
         allocate ( ptbublat(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptbublat(:,:) = 0
         allocate ( ptbubdepth(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptbubdepth(:,:) = 0
         allocate ( ptbubtime(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptbubtime(:,:) = 0
         allocate ( vectptbub(1:jpbubsize1), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectptbub(:) = 0
!
         CALL readptzon (kinzon,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime)
!
! --- allocation zone pointers
         allocate ( ptbubidx1(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptbubidx1(:,:) = 0
!
         IF (kflagincovxyoz.NE.4) THEN
            jpbub1 = jpzsize*dtaend
            jbub=0
            DO jz = 1,jpzsize
            DO jdta = 1,dtaend
               jbub=jbub+1
               ptbubidx1(jz,jdta)=jbub
            ENDDO
            ENDDO  
         ELSE
!
            allocate ( ptdtalon1(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptdtalon1(:,:) = 0
            allocate ( ptdtalat1(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptdtalat1(:,:) = 0
            allocate ( ptdtadepth1(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptdtadepth1(:,:) = 0
            allocate ( ptdtatime1(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptdtatime1(:,:) = 0
!
            allocate ( ptbublon1(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptbublon1(:,:) = 0
            allocate ( ptbublat1(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptbublat1(:,:) = 0
            allocate ( ptbubdepth1(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptbubdepth1(:,:) = 0
            allocate ( ptbubtime1(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptbubtime1(:,:) = 0
            CALL readptzon (tabfnin(1),ptbubidx1,ptdtalon1,ptdtalat1, &
     &           ptdtadepth1,ptdtatime1,ptbublon1,ptbublat1,ptbubdepth1, &
     &           ptbubtime1)
!
            IF (ANY(ptdtalon1(:,:).NE.ptdtalon(:,:))) GOTO 102
            IF (ANY(ptdtalat1(:,:).NE.ptdtalat(:,:))) GOTO 102
            IF (ANY(ptdtadepth1(:,:).NE.ptdtadepth(:,:))) GOTO 102
            IF (ANY(ptdtatime1(:,:).NE.ptdtatime(:,:))) GOTO 102
            IF (ANY(ptbublon1(:,:).NE.ptbublon(:,:))) GOTO 102
            IF (ANY(ptbublat1(:,:).NE.ptbublat(:,:))) GOTO 102
            IF (ANY(ptbubdepth1(:,:).NE.ptbubdepth(:,:))) GOTO 102
            IF (ANY(ptbubtime1(:,:).NE.ptbubtime(:,:))) GOTO 102
!
            IF (allocated(ptdtalon1)) deallocate(ptdtalon1)
            IF (allocated(ptdtalat1)) deallocate(ptdtalat1)
            IF (allocated(ptdtadepth1)) deallocate(ptdtadepth1)
            IF (allocated(ptdtatime1)) deallocate(ptdtatime1)
            IF (allocated(ptbublon1)) deallocate(ptbublon1)
            IF (allocated(ptbublat1)) deallocate(ptbublat1)
            IF (allocated(ptbubdepth1)) deallocate(ptbubdepth1)
            IF (allocated(ptbubtime1)) deallocate(ptbubtime1)
!
            jpbub1 = jpzsize*dtaend
            jbub=0
            DO jz = 1,jpzsize
            DO jdta = 1,dtaend
               jbub=jbub+1
               IF (ptbubidx1(jz,jdta).NE.jbub) GOTO 102
            ENDDO
            ENDDO  
!            
         ENDIF
!
         jrbasdeb1=0
         CALL writehdrzbas(koutcovoz,zon_jpi,zon_jpj,zon_jpk, &
     &        zon_jpt,jpbub1,jpzsize,jrbasdeb1,jprbasout)
         CALL writeptzbas(koutcovoz, &
     &        ptbubidx1,ptdtalon,ptdtalat, &
     &        ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &        ptbubtime,jrbasdeb1,jprbasout)
! --- write header std init 
         CALL writehdrzon (tabfnout(-1),zon_jpi,zon_jpj,zon_jpk, &
     &        zon_jpt,jpbub1,jpzsize)
         CALL writeptzon (tabfnout(-1),ptbubidx1,ptdtalon,ptdtalat, &
     &        ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &        ptbubtime)
!
! --- allocation bubr
         allocate ( bubr(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &        1:jpbubsize1,1:jprsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         bubr(:,:,:,:,:,:) = FREAL(0.0)
! --- allocation bubmode
         allocate ( bubmode(1:zon_jpi,1:zon_jpj,1:zon_jpk, &
     &        1:zon_jpt,1:jpbubsize1), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         bubmode(:,:,:,:,:) = FREAL(0.0)
! --- allocation bubstd
         allocate ( bubstd(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &        1:jpbubsize1), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         bubstd(:,:,:,:,:) = FREAL(0.0)
! --- allocation ptbub
         allocate ( ptbub(1:zon_jpi,1:zon_jpj,1:zon_jpk, &
     &        1:zon_jpt,1:jpbubsize1), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptbub(:,:,:,:,:) = 0
!
         DO jz=1,jpzsize,(jpbubsize1/dtaend)
            lmodprintjnz=(MOD((jz-1)/(jpbubsize1/dtaend), &
     &           (jpzsize/5/(jpbubsize1/dtaend)+1)).EQ.0)
            IF (lmodprintjnz) print *,  &
     &           'Zone number : ',jz,'/',jpzsize
!
            jbub=1+(jz-1)*dtaend
!
            jbubfin=MIN(jpbubsize1,jpbub1-(jbub-1))+jbub-1
            jbubend=jbubfin-(jbub-1)
            jzfin=jz-1+MIN((jpbubsize1/dtaend),jpzsize-(jz-1))
            jzend=jzfin-jz+1
!
! -11.2.2- Computing the ptbub
! -----------------------------
            
            DO jz1=jz,jzfin
            DO jdta=1,dtaend
               jbub1=ptbubidx1(jz1,jdta)-(jbub-1)
               CALL mkptbub (ptbub(:,:,:,:,jbub1),jdta, &
     &              ptdtalon(jz1,jdta),ptdtalat(jz1,jdta), &
     &              ptdtadepth(jz1,jdta),ptdtatime(jz1,jdta), &
     &              ptbublon(jz1,jdta),ptbublat(jz1,jdta), &
     &              ptbubdepth(jz1,jdta),ptbubtime(jz1,jdta))
            ENDDO
            ENDDO
!
! -11.2.3- Loading the bi-initial error covariance matrix (HzxA)
! -----------------------------------------------------------
            IF (kflagincovxyoz.EQ.4) THEN
!
               IF ((nprint.GE.2).AND.(jz.EQ.1)) THEN
!                  print *,'-11.2.1- Loading the bi-initial error ',
!     $                 'covariance matrix'
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOBASE : ', &
     &                 'loading part of A for covariance computation', &
     &                 ' (case kflagincovxyoz.EQ.4)'
               ENDIF
!
               vectptbub(1:jbubend) =  &
     &              (/ ((ptbubidx1(jz1,jdta), &
     &              jdta=1,dtaend),jz1=jz,jzfin) /)
!
               bubr(:,:,:,:,:,:)=FREAL(0.0)
               DO numjr=jrbasdeb,jrbasfin
                  CALL readnbubzon(tabfnin(numjr),vectptbub(1:jbubend), &
     &                 bubr(:,:,:,:,1:jbubend,numjr),lectinfo)
               ENDDO
               bubmode(:,:,:,:,:)=FREAL(0.0)
               bubstd(:,:,:,:,:)=FREAL(0.0)
               SELECT CASE (valbase)
               CASE (:0)
                  CALL meanstdnbub (bubmode(:,:,:,:,1:jbubend), &
     &                 bubstd(:,:,:,:,1:jbubend), &
     &                 bubr(:,:,:,:,1:jbubend,jrbasdeb:jrbasfin))
! --- ect init of the population => *SQRT(jpr/(jpr-1))
                  bubstd(:,:,:,:,1:jbubend) = bubstd(:,:,:,:,1:jbubend) &
     &                 * SQRT(FREAL(jprsize)/FREAL(jprsize-1))
               CASE (1:)
                  CALL readnbubzon(tabfnin(1),vectptbub(1:jbubend), &
     &                 bubmode(:,:,:,:,1:jbubend),lectinfo)
                  CALL readnbubzon(tabfnin(0),vectptbub(1:jbubend), &
     &                 bubstd(:,:,:,:,1:jbubend),lectinfo)
               CASE DEFAULT 
                  GOTO 1000
               END SELECT
               IF (valbase.EQ.0) THEN
                  DO jr=jrbasdeb,jrbasfin
                     bubr(:,:,:,:,1:jbubend,jr) =  &
     &                    bubr(:,:,:,:,1:jbubend,jr) -  &
     &                    bubmode(:,:,:,:,1:jbubend)
                  ENDDO
               ENDIF
            IF (.NOT.present(kinbasrefxyoz)) THEN
!
! -11.2.4- Writing mean file into output bibase directory :
! ---------------------------------------------------------------------
               IF ((nprint.GE.2).AND.(jz.EQ.1)) THEN
                  print *,'-11.2.4- Writing', &
     &                 ' output bibase directory'
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOBASE : ', &
     &                 'writing the bicovariance (HS)'
               ENDIF
               CALL writenbubzon(tabfnout(1), &
     &              bubmode(:,:,:,:,1:jbubend),jbub)
!
! -11.2.5- Write std init into output base directory    
! ---------------------------------------------------
               IF ((nprint.GE.2).AND.(jz.EQ.1)) THEN
!                  print *,'-11.2.5- Write',
!     $           ' std init into output bibase directory'
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOBASE : Write std init', &
     &                 ' into output bibase directory'
               ENDIF
               IF (valbase.EQ.0) THEN
! --- write std init 
                  CALL writenbubzon(tabfnout(-1), &
     &                 bubstd(:,:,:,:,1:jbubend),jbub)
               ENDIF
            ENDIF
!
            ELSE
!
! -11.2.6- Computing the local (bi-)initial error covariance matrix (HzxA)
! ----------------------------------------------------------------
               IF ((nprint.GE.2).AND.(jz.EQ.1)) THEN
!                  print *,'-11.2.4- Computing the local (bi-)initial',
!     $                 ' error covariance matrix (HzxA)'
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOBASE : Computing ', &
     &                 'the bi-initial error covariance matrix (HzxA)'
               ENDIF

               bubr(:,:,:,:,:,:) = FREAL(0.0)
!
               DO jz1=1,jzend
               DO jdta=1,dtaend
                  jbub1=jdta+(jz1-1)*dtaend
!CDIR NODEP
               DO jt=1,zon_jpt
!CDIR NODEP
               DO jk=1,zon_jpk
!CDIR NODEP
               DO jj=1,zon_jpj
!CDIR NODEP
               DO jr=jrbasdeb,jrbasfin
!CDIR NODEP
                  bubr(:,jj,jk,jt,jbub1,jr)  &
     &                 =basesr(ptbub(:,jj,jk,jt,jbub1),jr)
               ENDDO
               ENDDO
               ENDDO
               ENDDO
               ENDDO
               ENDDO
               IF (jbubend.NE.jbub1) GOTO 1000
!
! -11.2.7- Writing mean file into output bibase directory :
! ---------------------------------------------------------------------
            IF (.NOT.present(kinbasrefxyoz)) THEN
               IF ((nprint.GE.2).AND.(jz.EQ.1)) THEN
!                  print *,'-11.2.4- Writing',
!     $                 ' output bibase directory'
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOBASE : ', &
     &                 'writing the bicovariance (HS)'
               ENDIF
               bubmode(:,:,:,:,:) = FREAL (0.0)
!
               DO jz1=1,jzend
               DO jdta=1,dtaend
                  jbub1=jdta+(jz1-1)*dtaend
!CDIR NODEP
               DO jt=1,zon_jpt
!CDIR NODEP
               DO jk=1,zon_jpk
!CDIR NODEP
               DO jj=1,zon_jpj
!CDIR NODEP
                  bubmode(:,jj,jk,jt,jbub1)  &
     &                 = vectsmean(ptbub(:,jj,jk,jt,jbub1))
               ENDDO
               ENDDO
               ENDDO
               ENDDO
               ENDDO
               IF (jbubend.NE.jbub1) GOTO 1000
               CALL writenbubzon(tabfnout(1), &
     &              bubmode(:,:,:,:,1:jbubend),jbub)
!
! -11.2.8- Write std init into output base directory    
! ---------------------------------------------------
               IF ((nprint.GE.2).AND.(jz.EQ.1)) THEN
!                  print *,'-11.2.5- Write',
!     $                 ' std init into output bibase directory'
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOBASE : Write std init', &
     &                 ' into output bibase directory'
               ENDIF
               IF (valbase.EQ.0) THEN
! --- write std init
                  bubmode(:,:,:,:,:) = FREAL (0.0)
!
                  DO jz1=1,jzend
                  DO jdta=1,dtaend
                     jbub1=jdta+(jz1-1)*dtaend
!CDIR NODEP
                  DO jt=1,zon_jpt
!CDIR NODEP
                  DO jk=1,zon_jpk
!CDIR NODEP
                  DO jj=1,zon_jpj
!CDIR NODEP
                     bubmode(:,jj,jk,jt,jbub1)  &
     &                    = vectsect(ptbub(:,jj,jk,jt,jbub1))
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
                  IF (jbubend.NE.jbub1) GOTO 1000
                  CALL writenbubzon(tabfnout(-1), &
     &                 bubmode(:,:,:,:,1:jbubend),jbub)
               ENDIF
            ENDIF
!
            ENDIF
!
! -11.2.9- Compute and write the bi-initial error covariance matrix (HzxS)
! ------------------------------------------------------------------------
            IF ((nprint.GE.2).AND.(jz.EQ.1)) THEN
!               print *,'-11.2.2- Compute the bi-initial error',
!     $              ' covariance matrix (HzxS)'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE : ', &
     &              'compute the bi-initial', &
     &              ' error covariance matrix (HzxS)'
            ENDIF
            IF (present(kinbasrefxyoz)) THEN
               DO jd=jdbasdeb,jdbasfin
                  bubmode(:,:,:,:,1:jbubend) = FREAL (0.0)
                  DO jz1=1,jzend
                  DO jdta=1,dtaend
                     jbub1=jdta+(jz1-1)*dtaend
                  DO jrbas=jrbasdeb,jrbasfin
                     bubmode(:,:,:,:,jbub1) = bubmode(:,:,:,:,jbub1) + &
     &                    bubr(:,:,:,:,jbub1,jrbas)*covrdz(jrbas,jd,jz1+jz-1)
                  ENDDO
                  ENDDO
                  ENDDO
                  CALL writenbubzon(tabfnout(jd), &
     &                 bubmode(:,:,:,:,1:jbubend),jbub)
               ENDDO
            ELSE
               bubstd(:,:,:,:,1:jbubend) = FREAL (0.0)
               DO numjr=2,jprbasout
                  jr=numjr-(1-nvpnull)
                  bubmode(:,:,:,:,:) = FREAL (0.0)
                  DO jz1=1,jzend
                  DO jdta=1,dtaend
                     jbub1=jdta+(jz1-1)*dtaend
                  DO jrbas=jrbasdeb,jrbasfin
                     bubmode(:,:,:,:,jbub1) = bubmode(:,:,:,:,jbub1) + &
     &                    bubr(:,:,:,:,jbub1,jrbas)* &
     &                    matUrrz(jrbas,jr,jz1+jz-1)
                  ENDDO
                  ENDDO
                  ENDDO
                  bubstd(:,:,:,:,1:jbubend) = bubstd(:,:,:,:,1:jbubend) +  &
     &                 bubmode(:,:,:,:,1:jbubend)* &
     &                 bubmode(:,:,:,:,1:jbubend)
                  CALL writenbubzon(tabfnout(numjr), &
     &                 bubmode(:,:,:,:,1:jbubend),jbub)
               ENDDO
            ENDIF
!
! -11.2.6- Compute and write diagonal Pa into output base directory    
! ------------------------------------------------------------------
            IF (.NOT.present(kinbasrefxyoz)) THEN
            IF ((nprint.GE.2).AND.(jz.EQ.1)) THEN
!            print *,'-11.2.5- Compute and write',
!     $           ' diagonal Pa into output bibase directory'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE : Compute and write diagonal', &
     &              ' HPH into output bibase directory'
            ENDIF
! --- write ect end
            bubstd(:,:,:,:,1:jbubend) = SQRT(bubstd(:,:,:,:,1:jbubend))
            CALL writenbubzon(tabfnout(0),bubstd(:,:,:,:,1:jbubend),jbub)
            ENDIF
!
         ENDDO
!     
         CALL writeinfobas(koutcovoz,kflagbicovoz)
!
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -12.- Allocation/desallocation
! ------------------------------
!
! --- desallocate some arrays
      IF (allocated(vectspart))  deallocate (vectspart)
      IF (allocated(basesr))  deallocate (basesr)
      IF (allocated(vectsmean))  deallocate (vectsmean)
      IF (allocated(vectsect))  deallocate (vectsect)
      IF (allocated(vectsmode))  deallocate (vectsmode)
      IF (allocated(vectsstd))  deallocate (vectsstd)
      IF (allocated(bubmode)) deallocate(bubmode)
      IF (allocated(bubstd)) deallocate(bubstd)
      IF (allocated(bubr)) deallocate(bubr)
      IF (allocated(ptbub)) deallocate(ptbub)
      IF (allocated(ptbubidx)) deallocate(ptbubidx)
      IF (allocated(ptbubidx1)) deallocate(ptbubidx1)
      IF (allocated(ptdtalon)) deallocate(ptdtalon)
      IF (allocated(ptdtalat)) deallocate(ptdtalat)
      IF (allocated(ptdtadepth)) deallocate(ptdtadepth)
      IF (allocated(ptdtatime)) deallocate(ptdtatime)
      IF (allocated(ptbublon)) deallocate(ptbublon)
      IF (allocated(ptbublat)) deallocate(ptbublat)
      IF (allocated(ptbubdepth)) deallocate(ptbubdepth)
      IF (allocated(ptbubtime)) deallocate(ptbubtime)
!
! --- allocation vectxpart
      allocate ( vectxpart(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxpart(:) = FREAL(0.0)
! --- allocation basexr
      allocate ( basexr(1:jpxsize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basexr(:,:) = FREAL(0.0)
! --- allocation vectxmean
      allocate ( vectxmean(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxmean(:) = FREAL(0.0)
! --- allocation vectxect
      allocate ( vectxect(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxect(:) = FREAL(0.0)
!
! -13.- Treatment of the initial error covariance matrix (S)
! ----------------------------------------------------------
!
      IF (nprint.GE.1) print *, &
     &     '--> Computation of (P) ...'
      lectinfo = .FALSE.
      DO jnxyo=1,limjpnxyo(kflaganlxyo)
         lmodprintjnxyo=(MOD(jnxyo-1,(limjpnxyo(kflaganlxyo)/3+1)).EQ.0)
         IF ((lmodprintjnxyo).AND.(nprint.GE.1)) &
     &        print *,'Memory part number : ', &
     &           jnxyo,'/',limjpnxyo(kflaganlxyo)
!
! -13.1- Read the partition in subdomains
! --------------------------------------
         IF (kflaggloloc.NE.0) THEN
            IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
!               print *,'-13.1- Read the partition in subdomains'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE: ', &
     &              'reading the partition in subdomains'
            ENDIF
            IF (kflaganlxyo.NE.3) THEN
               lmoyectold=lmoyect
               lmoyect=.FALSE.
               CALL readxyo (kinpartxyo,vectxpart(:),jnxyo, &
     &              lectinfo,kflaganlxyo)
               lmoyect=lmoyectold
            ELSE
! --- allocation vecty
               allocate ( vecty(1:jpysize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               vecty(:) = FREAL(0.0)
               flagxyo=2
               lmoyectold=lmoyect
               lmoyect=.FALSE.
               CALL readxyo (kinpartxyo,vecty(:),jnxyo, &
     &              lectinfo,flagxyo,poscoefobs(:,:))
               lmoyect=lmoyectold
               vectxpart(:)=vecty(poscoefobs(:,1)%pos)
               IF (allocated(vecty))  deallocate (vecty)
            ENDIF
         ENDIF
!
! -13.2- Read the input basis A
! -----------------------------
         IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
!            print *,'-13.2- Read the input basis A'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : ', &
     &           'reading the input basis A'
         ENDIF
         IF (kflagcovxyo.EQ.3) THEN
! --- load base.o
            CALL readbas (kinbasxyo,basexr(:,:), &
     &           jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflaganlxyo,poscoefobs(:,:))
         ELSE
! --- load base.xy
            CALL readbas (kinbasxyo,basexr(:,:), &
     &           jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflaganlxyo)
         ENDIF
         IF (valbase.LE.0) CALL meanect(vectxmean(:),vectxect(:), &
     &           basexr(:,jrbasdeb:jrbasfin))
         IF (valbase.EQ.0) THEN
            DO jr=jrbasdeb,jrbasfin
               basexr(:,jr) = basexr(:,jr) - vectxmean(:)
            ENDDO
         ENDIF
!
! -13.3- Compute the initial error covariance matrix (S)
! ----------------------------------------------------
         IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
!            print *,'-13.- Compute the ',
!     $           'error covariance matrix (S)'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : ', &
     &           'computing the initial error covariance matrix (S)'
         ENDIF
         numjr=0
         serie=1
         CALL fildirbas (vctnamout,koutbasxyo,jprbasout,numjr,serie)
         IF (present(kinbasrefxyoz)) THEN
            IF (kflaggloloc.EQ.0) THEN
               CALL prodmat_zr_rd_vias (basexr(:,jrbasdeb:jrbasfin), &
     &              mat_yr(:,jrbasdeb:jrbasfin), &
     &              covrdz(jrbasdeb:jrbasfin,jdbasdeb:jdbasfin,1), &
     &              jdbasdeb,jdbasfin)
            ELSE
               CALL prodmat_wr_rdz_vias(basexr(:,jrbasdeb:jrbasfin), &
     &              mat_yr(:,jrbasdeb:jrbasfin), &
     &              covrdz(jrbasdeb:jrbasfin,jdbasdeb:jdbasfin,0:), &
     &              vectxpart(:), &
     &              jdbasdeb,jdbasfin)  
            ENDIF
         ELSE
            IF (kflaggloloc.EQ.0) THEN
               CALL prodmat_zr_rr_vias (basexr(:,jrbasdeb:jrbasfin), &
     &              mat_yr(:,jrbasdeb:jrbasfin), &
     &              matUrrz(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin,1), &
     &              (2-(1-nvpnull)),(jprbasout-(1-nvpnull)))
!
            ELSE
               CALL prodmat_wr_rrz_vias(basexr(:,jrbasdeb:jrbasfin), &
     &              mat_yr(:,jrbasdeb:jrbasfin), &
     &              matUrrz(jrbasdeb:jrbasfin,jrbasdeb:jrbasfin,0:), &
     &              vectxpart(:), &
     &              (2-(1-nvpnull)),(jprbasout-(1-nvpnull)))  
            ENDIF
         ENDIF
!
! -13.4- Write or transfer mean file into output base directory 
! --------------------------------------------------------------
         IF (.NOT.(present(kinbasrefxyoz))) THEN
            IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
!               print *,'-13.4- Write or transfer mean'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE : ', &
     &              'Write or transfer mean file ', &
     &              'into output base directory'
            ENDIF
            numjr=1
            serie=1
            CALL fildirbas (vctnamin,kinbasxyo,jprbasin,numjr,serie)
            WRITE(fnamein,'("./",A,"/",A)') kinbasxyo(1:lenv(kinbasxyo)), &
     &           vctnamin(1:lenv(vctnamin))
            IF (valbase.GE.1) THEN
               IF (kflaganlxyo.EQ.3) THEN
                  CALL readxyo(fnamein,basexr(:,1), &
     &                 jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
               ELSE
                  CALL readxyo(fnamein,basexr(:,1), &
     &                 jnxyo,lectinfo,kflaganlxyo)
               ENDIF
            ELSE
               basexr(:,1) = vectxmean(:)
            ENDIF
         ENDIF
!
! -13.5- Writing output base directory 
! -------------------------------------
         IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
!            print *,'-13.5- Writing output base directory'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOBASE : writing the initial', &
     &           ' error covariance matrix (S) (S)'
         ENDIF
         jrbasdeb1=1
         IF (present(kinbasrefxyoz)) THEN
            jrbasdeb1=jdbasdeb
            jprbasout=jdbasfin
         ENDIF
         SELECT CASE (kflaganlxyo)
         CASE (1)
            CALL writebas(koutbasxyo,basexr(:,:), &
     &           jnxyo,jrbasdeb1,jprbasout)
         CASE (2)
            CALL writeyobas(koutbasxyo,basexr(:,:), &
     &           jrbasdeb1,jprbasout,kflaganlxyo)
         CASE (3)
            CALL writeyobas(koutbasxyo,basexr(:,:), &
     &           jrbasdeb1,jprbasout,kflaganlxyo,vectorms(:), &
     &           gridijkobs(:),poscoefobs(:,:))
         CASE DEFAULT
            GOTO 1000
         END SELECT
! --- affectation of outbas
         IF (jnxyo.EQ.1) CALL writeinfobas(koutbasxyo,kflaganlxyo)
!
! -13.6- Compute and write diagonal Pa into output base directory
! -----------------------------------------------------------------
         IF (.NOT.(present(kinbasrefxyoz))) THEN
            IF (nprint.GE.2) THEN
!               print *,'-13.6- Compute and write diagonal'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOBASE : ', &
     &              'Compute and write diagonal ', &
     &              'Pa into output base directory'
            ENDIF
            IF (valbase.EQ.0) THEN
!     --- write ect init of the population => *SQRT(jpr/(jpr-1))
               vectxect(:) = vectxect(:) &
     &              * SQRT(FREAL(jprsize)/FREAL(jprsize-1))
               numjr=-1
               serie=1
               CALL fildirbas (vctnamout,koutbasxyo, &
     &              jprbasout,numjr,serie)
               SELECT CASE (kflaganlxyo)
               CASE (1)
                  jext=indext(vctnamout,extvartab,nbextvar)
                  IF (extvarunit(jext)) THEN
                     WRITE (text,'(A)') extvartab(jext)
                  ELSE
                     WRITE (text,'("#",A)') extvartab(jext)
                  ENDIF
               CASE (2)
                  jext=indext(vctnamout,extdtatab,nbextdta)
                  IF (extdtaunit(jext)) THEN
                     WRITE (text,'(A)') extdtatab(jext)
                  ELSE
                     WRITE (text,'("#",A)') extdtatab(jext)
                  ENDIF
               CASE (3)
                  jext=indext(vctnamout,extobstab,nbextobs)
                  IF (extobsunit(jext)) THEN
                     WRITE (text,'(A)') extobstab(jext)
                  ELSE
                     WRITE (text,'("#",A)') extobstab(jext)
                  ENDIF
               CASE DEFAULT
                  GOTO 1000
               END SELECT
               WRITE(fnameout,'("./",A,"/vctectinit",A)')  &
     &              koutbasxyo(1:lenv(koutbasxyo)), &
     &              text(1:lenv(text))
               SELECT CASE (kflaganlxyo)
               CASE (1)
                  CALL writevar (fnameout,vectxect(:),jnxyo)
               CASE (2)
                  CALL writedta (fnameout,vectxect(:))
               CASE (3)
                  CALL writeobs(fnameout,vectxect(:),vectorms(:), &
     &                 gridijkobs(:),poscoefobs(:,:))
               CASE DEFAULT
                  GOTO 1000
               END SELECT
            ENDIF
!     --- compute ect end
            basexr(:,1)=FREAL(0.0)
            DO jr=jrbasdeb,jrbasfin
               basexr(:,1)= basexr(:,jr)*basexr(:,jr)+basexr(:,1)
            ENDDO
            basexr(:,1)=SQRT(basexr(:,1))
! --- write ect end
            numjr=0
            serie=1
            CALL fildirbas (vctnamout,koutbasxyo,jprbasout,numjr,serie)
            WRITE(fnameout,'("./",A,"/",A)') koutbasxyo(1:lenv(koutbasxyo)), &
     &           vctnamout(1:lenv(vctnamout))
            SELECT CASE (kflaganlxyo)
            CASE (1)
               CALL writevar (fnameout,basexr(:,1),jnxyo)
            CASE (2)
               CALL writedta (fnameout,basexr(:,1))
            CASE (3)
               CALL writeobs(fnameout,basexr(:,1),vectorms(:), &
     &              gridijkobs(:),poscoefobs(:,:))
            CASE DEFAULT
               GOTO 1000
            END SELECT
         ENDIF
!     
      ENDDO
!
! -14.- Algortho :
! -----------------
! --- orthogonality computation
      IF ((kflaggloloc.EQ.0).AND.(.NOT.(present(kinbasrefxyoz))) &
     &     .AND.(kflaganlxyo.EQ.kflagcovxyo)) THEN
         print *,'Algortho :  orthogonality computation...' 
         IF (largweight) THEN
            CALL chkortho(koutbasxyo,kflaganlxyo,basexr(:,:), &
     &           2,jprbasout,largweight, &
     &           kargxyoweight=argweight,kvectsweight=basexr(:,1))
         ELSE
            CALL chkortho(koutbasxyo,kflaganlxyo,basexr(:,:), &
     &           2,jprbasout,largweight)
         ENDIF
      ENDIF
!
! --- deallocation
      IF (allocated(vectorms)) deallocate (vectorms)
      IF (allocated(poscoefobs)) deallocate (poscoefobs)
      IF (allocated(gridijkobs)) deallocate (gridijkobs)
      IF (allocated(vectxpart)) deallocate (vectxpart)
      IF (allocated(vectxmean)) deallocate (vectxmean)
      IF (allocated(vectxect)) deallocate (vectxect)
      IF (allocated(basexr)) deallocate (basexr)
      IF (allocated(matUrrz)) deallocate (matUrrz)
      IF (allocated(covrdz)) deallocate (covrdz)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algobase','algobase')
 1001 CALL printerror2(0,1001,3,'algobase','algobase')
!
 101  WRITE (texterror,*) 'argument not valid :'
      CALL printerror2(0,101,3,'algobase','algobase',comment=texterror)
 102  WRITE (texterror,*) 'file ',tabfnin(1)(1:lenv(tabfnin(1))), &
     &     ' incoherent with ',kinzon(1:lenv(kinzon))
      CALL printerror2(0,102,3,'algobase','algobase',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algobase
