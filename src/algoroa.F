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
C ---                   ALGOROA.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 98-06 (C.E. Testut)                        ---
C --- modification : 08-09 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE calcroa
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algoroa
      use mod_main
      use mkconnect
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcroa

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcroa(
     $     kflaganlxyo,koutxyo,kinrefxyo,kinbasxyo,koutbasxyo,
     $     kflagcovyo,kinyo,kinrefyo,
     $     kflagincovyoz,kincovyoz,
     $     kflagbicovoz,koutcovoz,
     $     kflaggloloc,kinpartxyo,kinzon,kconfigo,
     $     ktextdisable)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Implement a global & local version of the ROA filter algorithm
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
      use mod_spacexyo , only : jpoend,jpitpend,jprend,jpx,jpxend,
     $     jpyend,jpz,spvaldta,poscoefobs,gridijkobs,
     $     pt1bubidx, pt2bubidx, pt3bubidx, pt4bubidx, 
     $     pt1dtalon, pt1dtalat, pt1dtadepth, pt1dtatime,
     $     pt1bublon, pt1bublat, pt1bubdepth, pt1bubtime,
     $     bubblk1, bubblk2, bubblk3, bubblk4, vectptbub, ptbub
      use hioxyo
      use hiobas
      use hiozon
      use hiocfg
      use liocrz
      use utilroa
      use utilmkh
      use utilmkto
      use utilclc
      use utilprint
      use utilvalid
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: kflaganlxyo,kflagcovyo,
     $     kflagincovyoz,kflagbicovoz,kflaggloloc
      CHARACTER(len=*), intent(in) :: koutxyo,kinrefxyo,kinbasxyo,
     $     koutbasxyo,kinyo,kinrefyo,kincovyoz,
     $     koutcovoz,kinpartxyo,kinzon,kconfigo,ktextdisable
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vectssqrdiagRi
      BIGREAL, dimension(:), allocatable, save :: vectsinnov
      BIGREAL, dimension(:), allocatable, save :: vectsweight
      BIGREAL, dimension(:), allocatable, save :: vectspart
      BIGREAL, dimension(:), allocatable, save :: vectorms
      BIGREAL, dimension(:,:), allocatable, save :: basesr
      BIGREAL, dimension(:,:), allocatable, save :: basexr
      BIGREAL, dimension(:,:), allocatable, save :: mat_yr
      BIGREAL, dimension(:), allocatable, save :: vectxpart
      BIGREAL, dimension(:), allocatable, save :: vectxweight
      BIGREAL, dimension(:), allocatable, save :: vectxa
      BIGREAL, dimension(:), allocatable, save :: vectxf
      BIGREAL, dimension(:), allocatable, save :: vecty
      BIGREAL, dimension(:), allocatable, save :: vectybub
C
      BIGREAL, dimension(:,:), allocatable :: lambda
      BIGREAL, dimension(:,:), allocatable :: delta
      BIGREAL, dimension(:,:), allocatable :: coefrz
      BIGREAL, dimension(:,:), allocatable :: coefrznolim
C
      BIGREAL, dimension(:,:), allocatable :: matBqr
C
      BIGREAL, dimension(:,:), allocatable :: iRi_zi
      BIGREAL, dimension(:,:), allocatable :: y_zi
      BIGREAL, dimension(:,:), allocatable :: vectbeta
      BIGREAL, dimension(:), allocatable :: vectalpha
      BIGREAL, dimension(:,:,:), allocatable :: deltarzi
      BIGREAL, dimension(:,:,:,:), allocatable :: matGrrzi
      BIGREAL, dimension(:,:,:), allocatable :: matCUrrz
C
      INTEGER :: allocok,jpxsize,jpssize,jpysize,jpqsize,jpusize,jpbsize
      INTEGER :: jprsize,jpitpsize,jpzsize,jpu,maxjpu,nbzonnull,jb,jitp
      INTEGER :: jprbasin,jprbasout,jpiobs
      INTEGER :: valbase,jpnbubsize
      BIGREAL :: coefrmax, coeflimite, forgfact, ratio
      INTEGER :: jnxyo,jr,js,jz,ju,jx,jq,ji,jj,jk,jt,jiobs
      INTEGER :: jbub,jdta,jz1,jz2,jnz,ios
      INTEGER :: jrbasdeb,jrbasfin,jrmatdeb,jrmatfin
      INTEGER :: jqbasdeb,jqbasfin,jqmatdeb,jqmatfin
      INTEGER :: serie,numjr,flagcfg,flagxyo
      INTEGER :: jpzalloc,jproc1,jjproc
      CHARACTER(len=bgword) :: fnamein,fnameout,vctnamin,vctnamout
      CHARACTER(len=bgword) :: titre,variable,dirnambas,text
      LOGICAL :: lmodprint,zlectinfo,lectinfo,lmoyectold
      LOGICAL :: lXa,lXaout,lPaout,linnov,lbulkave
      INTEGER :: jext,jbubfin,jbubend,jzfin,jzend,jbub1,nvpnull,jrbas
      LOGICAL :: lmodprintjnz
CC----------------------------------------------------------------------
      INTEGER, dimension(:), allocatable :: ptu,ptubub,ptlinbub
      BIGREAL, dimension(:,:), allocatable :: baseur
      BIGREAL, dimension(:,:), allocatable :: baseuq
      BIGREAL, dimension(:), allocatable :: vectuweight,bubblklin
      BIGREAL, dimension(:), allocatable :: vectuinnov
      BIGREAL, dimension(:), allocatable :: vectupart
      BIGREAL :: ensmean
      TYPE (type_poscoef), dimension(:,:), allocatable :: poscoefuobs,
     $     poscoefbubobs
      INTEGER, dimension(:), allocatable :: vectynbpoint
CC----------------------------------------------------------------------
      INTEGER :: fixjpbub
CC----------------------------------------------------------------------
      INTEGER :: jsmo,jpsmo
      CHARACTER(len=bgword), dimension (:), allocatable :: vctnaminsmo,
     $     vctnamoutsmo,dirnaminsmo,dirnamoutsmo
CC----------------------------------------------------------------------
C
      lXa=(ktextdisable(1:1).EQ.'T')
      lXaout=(ktextdisable(2:2).EQ.'T')
      IF ((.NOT.(lXa)).AND.(lXaout)) print *,' Warning : outvar = Xa-Xf'
      lPaout=(ktextdisable(3:3).EQ.'T')
      lbulkave=(ktextdisable(4:4).EQ.'T')
      linnov=(ktextdisable(5:5).EQ.'T')
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
      SELECT CASE (kflagcovyo)
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
      jpusize=jpfixjpu
C
      jpzalloc=jpzsize
#if defined MPI
      IF (MOD(jpzalloc,jpproc).NE.0) THEN
         jpzalloc=(jpzalloc/jpproc+1)*jpproc
      ENDIF
      IF (.NOT.largfixjpx) GOTO 103
#endif
C
C --- allocation basesr
      allocate ( basesr(1:jpssize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basesr(:,:) = FREAL(0.0)
C --- allocation vectsinnov
      allocate ( vectsinnov(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectsinnov(:) = FREAL(0.0)
C --- allocation vectssqrdiagRi
      allocate ( vectssqrdiagRi(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectssqrdiagRi(:) = FREAL(0.0)
C --- allocation vectsweight
      allocate ( vectsweight(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectsweight(:) = FREAL(0.0)
C --- allocation matCUrrz
      allocate ( matCUrrz(1:jprsize,1:jprsize,0:jpzalloc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      matCUrrz(:,:,0:) = FREAL(0.0)
C --- allocation lambda
      allocate ( lambda(1:jprsize,0:jpzalloc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lambda(:,0:) = FREAL(0.0)
C --- allocation delta
      allocate ( delta(1:jprsize,0:jpzalloc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      delta(:,0:) = FREAL(0.0)
C --- allocation coefrz
      allocate ( coefrz(1:jprsize,0:jpzalloc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      coefrz(:,0:) = FREAL(0.0)
C --- allocation coefrznolim
      IF (largcoefrmax) THEN
         allocate ( coefrznolim(1:jprsize,0:jpzalloc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         coefrznolim(:,0:) = FREAL(0.0)
      ENDIF
CCC---------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine algoroa &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
C
C -0.- Initialisation :
C ---------------------
C
C -0.1- Read coefrmax :
C ---------------------
C
C --- argcoefrmax
      IF (largcoefrmax) THEN
         READ(argcoefrmax,*,IOSTAT=ios) coefrmax
         IF (ios.NE.0) GOTO 101
      ELSE
         coefrmax=FREAL(0.0)
      ENDIF
C
C -0.2- Read or define config.obs :
C ---------------------------------
C
C --- allocation poscoefobs
      allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      IF (kflagcovyo.EQ.2) THEN
         DO js=1,jpssize
            poscoefobs(js,:) = type_poscoef(js,FREAL(1.0))
         ENDDO
      ELSE
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
         flagcfg=3
         CALL readcfgobs (kconfigo,flagcfg,
     $        kposcoefobs=poscoefobs(:,:))
      ENDIF
C
      IF (((lPaout).OR.(lXaout)).AND.((kflagcovyo.EQ.3).AND.
     $     ((kflagbicovoz.EQ.3).OR.(kflaganlxyo.EQ.3)))) THEN
C
C --- allocation gridijkobs
         allocate ( gridijkobs(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
         flagcfg=2
         CALL readcfgobs (kconfigo,flagcfg,
     $        kgridijkobs=gridijkobs(:))
C --- allocation vectorms
         allocate ( vectorms(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectorms(:)=FREAL(0.0)
         flagcfg=1
         CALL readcfgobs (kconfigo,flagcfg,
     $        kvectorms=vectorms(:))
C
      ENDIF
C
      serie=0
      numjr=1
      CALL fildirbas (fnamein,kinbasxyo,jprbasin,numjr,serie)
      CALL fildirbas (fnameout,koutbasxyo,jprbasout,numjr,serie)
      IF (jprbasin.GT.jprend) GOTO 1000
      IF (jprbasout.GT.jprend) GOTO 1000
      IF (jprbasin.NE.jprbasout) GOTO 1000
      CALL readinfobas(kinbasxyo,valbase)
C
C
C -0.3- Define some parameters :
C ------------------------------
C
      IF (valbase.LT.1) THEN
        PRINT *, 'Warning: input directory is an ensemble'
        jrbasdeb = 1
        jrbasfin = jprsize
        jrmatdeb = 1
        jrmatfin = jprsize
      ELSE
        PRINT *, 'Warning: input directory is a covariance square root'
        jrbasdeb = 2
        jrbasfin = jprsize
        jrmatdeb = 1
        jrmatfin = jprsize-1
      ENDIF
C
      IF (largoecorrel) THEN
         serie=0
         numjr=1
         CALL fildirbas (fnamein,argoecorrel,jpqsize,numjr,serie)
         jqbasdeb = 2
         jqbasfin = jpqsize
C --- allocation matBqr
         allocate ( matBqr(1:jpqsize,1:jprsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         matBqr(:,:) = FREAL(0.0)
      ENDIF
C
C -1.- Compute Weight*R-1 :
C -------------------------
C
C -1.1- Read weight vector :
C --------------------------
C
      IF (largweight) THEN
         IF (nprint.GE.1) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOROA : ',
     $           'Reading weight vector'
         ENDIF
C     
         IF ((validextvar(argweight)).OR.(validextdta(argweight))
     $        .OR.(validextobs(argweight))) THEN
            CALL readxyo(argweight,basesr(:,1),
     $           jnxyo,lectinfo,kflagcovyo,poscoefobs(:,:))   
         ELSE 
            GOTO 1000
         ENDIF   
C
         IF (nprint.GE.3) THEN
            titre='A few elements of the weight vector'
            variable='Wi(*)'
            CALL printtab_r (basesr(1:MIN0(jprsize,5),1),
     $           titre,variable,1) 
         ENDIF
      ENDIF
C
C -1.2- Read (or fill) and inverse observation error std vector (diag R)
C ----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ALGOROA : ',
     $          'filling observation error std vector (diag R)'
      ENDIF
C
      IF (largoestd) THEN
         IF ((validextvar(argoestd)).OR.(validextdta(argoestd))
     $        .OR.(validextobs(argoestd))) THEN
            CALL readxyo(argoestd,vectssqrdiagRi(:),
     $           jnxyo,lectinfo,kflagcovyo,poscoefobs(:,:))
         ELSE 
            GOTO 1000
         ENDIF   
      ELSE
         CALL mkyorms (vectssqrdiagRi(:),kflagcovyo)
      ENDIF
C
      vectssqrdiagRi(:) = FREAL(1.0)/(ABS(vectssqrdiagRi(:)))
C
      IF (nprint.GE.3) THEN
         titre='A few elements of sqrt ( inv(R) )'
         variable='Ri(*)'
         CALL printtab_r (vectssqrdiagRi(1:MIN0(jprsize,5)),
     $        titre,variable,1) 
      ENDIF
C
      IF (largweight) THEN
         vectssqrdiagRi(:) = vectssqrdiagRi(:)*basesr(:,1)
      ENDIF
C
C -2.- Compute the innovation vector :
C ------------------------------------
C
         IF (nprint.GE.1) print *,
     $        '--> Load and compute the innovation vector...'
C
C -2.1- Read observation vector Yt or Ot :
C ----------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ALGOROA : ',
     $                'loading observation vector (y)'
      ENDIF
C
      jnxyo=1
      lectinfo=.FALSE.
      CALL readxyo(kinyo,vectsinnov(:),
     $     jnxyo,lectinfo,kflagcovyo,poscoefobs(:,:)) 
C
C -2.2- Remove a bias from forecast (Ot -> Ot + bias) :
C -----------------------------------------------------
C
      IF (largbias) THEN
C
         IF (nprint.GE.1) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOROA: ',
     $           'loading bias state at observation points (Hbias)'
         ENDIF
         jnxyo=1
         lectinfo=.FALSE.
         CALL readxyo(argbias,basesr(:,1),
     $        jnxyo,lectinfo,kflagcovyo,poscoefobs(:,:))  
         vectsinnov(:) = vectsinnov(:) + basesr(:,1)
C
      ENDIF
C
C -2.3- Compute the weight correction to the observation correlation matrix
C -------------------------------------------------------------------------
C
      IF (largoecorrel) THEN
         IF (kflagcovyo.NE.2) THEN
C
C --- allocation vecty
            allocate ( vecty(1:jpysize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vecty(:) = FREAL(0.0)
C --- allocation vectynbpoint
            allocate ( vectynbpoint(1:jpysize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vectynbpoint(:) = 0
C
            coeflimite=FREAL(0.75)
            CALL mkhotoy(basesr(:,1),vecty(:),
     $             poscoefobs(:,:),spvaldta,kcoeflimite=coeflimite,
     $             kvectynbpoint=vectynbpoint(:))
            vecty(:) = FREAL(vectynbpoint(:))
            CALL mkhytoo(vecty(1:jpysize),vectsweight(1:jpssize),
     $                            poscoefobs(:,:))
C
            vectssqrdiagRi(:) = vectssqrdiagRi(:) 
     $                   / MAX(FREAL(1.0),SQRT(vectsweight(:)))
C
            vectsweight(:) = FREAL(0.0)
            IF (allocated(vecty)) deallocate(vecty)
            IF (allocated(vectynbpoint)) deallocate(vectynbpoint)
C
         ENDIF
      ENDIF
C
C -3.- Load forecast error covariance matrix at observation points (HS)
C ---------------------------------------------------------------------
C
      IF (kflagincovyoz.NE.4) THEN
         IF (nprint.GE.1) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOROA : ',
     $           'loading forecast error covariance matrix'
            WRITE(numout,*) '               at observation points (HS)'
         ENDIF
C --- load base.yo
         jnxyo=1
         CALL readbas(kincovyoz,basesr(:,:),jnxyo,
     $        jrbasdeb,jrbasfin,lectinfo,kflagcovyo,poscoefobs(:,:))
      ENDIF
C
C -3.1- Load forecast state at observation points (HyxXf) or (HoxXf) or Yf or (HoyYf) or Of :
C -------------------------------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ALGOROA: ',
     $        'loading forecast state at observation points (HXf)'
      ENDIF
C
      IF (valbase.LT.1) THEN
        DO js=1,jpssize
          ensmean = SUM(basesr(js,:))/jprsize
          basesr(js,:) = basesr(js,:) - ensmean
          basesr(js,:) = basesr(js,:) / SQRT(FREAL(jprsize-1))
C
          IF (linnov) THEN
             vectsinnov(js) = vectsinnov(js) - ensmean
          ENDIF
        ENDDO
      ELSE
        jnxyo=1
        lectinfo=.FALSE.
        CALL readxyo(kinrefyo,basesr(:,1),
     $       jnxyo,lectinfo,kflagcovyo,poscoefobs(:,:))
C
        IF (linnov) THEN
           vectsinnov(:) = vectsinnov(:) - basesr(:,1)
        ENDIF
      ENDIF
C
C -3a.- Read observation partition file
C -------------------------------------
C
      IF (larginpartobs) THEN
C
         IF ((.NOT.validextdta(arginpartobs)).AND.
     $       (.NOT.validextobs(arginpartobs))) GOTO 102
         IF ((kflagcovyo.EQ.2).AND.(.NOT.validextdta(arginpartobs))) GOTO 102
C --- allocation vectspart
         allocate ( vectspart(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C
         CALL readxyo(arginpartobs,vectspart(:),
     $        jnxyo,lectinfo,kflagcovyo,poscoefobs(:,:))
C
         jpiobs=NINT(MAXVAL(vectspart(:)))
C
C --- allocation deltarzi
         allocate ( deltarzi(1:jprsize,0:jpzalloc,1:jpiobs),
     $              stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C --- allocation matGrrzi
         allocate ( matGrrzi(1:jprsize,1:jprsize,0:jpzalloc,1:jpiobs),
     $              stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C
      ELSE
         jpiobs=1
      ENDIF
C
C -3b.- Read vectors of adaptive parameters
C -----------------------------------------
C
      IF (larginparadap) THEN
C --- allocation vectalpha
         allocate ( vectalpha(0:jpzalloc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectalpha(0:jpzalloc) = FREAL(1.0)
C --- allocation vectbeta
         allocate ( vectbeta(0:jpzalloc,1:jpiobs), stat=allocok )
         vectbeta(0:jpzalloc,1:jpiobs) = FREAL(1.0)
         IF (allocok.NE.0) GOTO 1001
C
C Scaling of the forecast error covariance matrix (forgetting factor)
         CALL readveccrz(arginparadap,'alpha',vectalpha(1:jpzsize))
         IF (ANY(vectalpha(1:jpzsize).LE.FREAL(0.0))) GOTO 104
C
C Scaling of the observation error covariance matrix
         CALL readmat2crz(arginparadap,'beta',vectbeta(1:jpzsize,:))
         IF (ANY(vectbeta(1:jpzsize,1:jpiobs).LE.FREAL(0.0))) GOTO 105
C
      ENDIF
C
C -3c.- Prepare output file for saving reduced rank elements
C ----------------------------------------------------------
C
      IF (largoutrz) THEN
C --- allocation iRi_zi
         allocate ( iRi_zi(0:jpzalloc,1:jpiobs), stat=allocok )
         iRi_zi(0:jpzalloc,1:jpiobs) = FREAL(0.0)
         IF (allocok.NE.0) GOTO 1001
C --- allocation y_zi
         allocate ( y_zi(0:jpzalloc,1:jpiobs), stat=allocok )
         y_zi(0:jpzalloc,1:jpiobs) = FREAL(0.0)
         IF (allocok.NE.0) GOTO 1001
C
         CALL writehdrcrz(argoutrz,jprsize-1,jpzsize,jpiobs)
      ENDIF
C
C -4.- Loop on subdomains to compute local ROA coefficients
C ----------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ALGOROA : ',
     $                'beginning loop over subdomains'
         WRITE(numout,*) '              n = ',jpzsize
      ENDIF
C
C -4.1- Prepare arrays for the management of the influence bubbles
C ----------------------------------------------------------------
C
#if defined MPI
      IF (jpfixjpz.LT.jpz) THEN
         IF (MOD(jpfixjpz,jpproc).NE.0) THEN
            jpfixjpz=(jpfixjpz/jpproc+1)*jpproc
         ENDIF
      ENDIF
#endif
C
      fixjpbub = jpfixjpz*dtaend
      limjpnz(:) = 1
      jpbubend(:)=0
      jpnbub(:)=0
      IF (kflaggloloc.NE.0) THEN
         CALL evalhdrzon(kinzon,zon_jpi,zon_jpj,zon_jpk,zon_jpt,
     $        jpbubend(1),jpz)
         IF (jpbubend(1).GT.fixjpbub) THEN
            jpnbub(1) = fixjpbub
            limjpnz(1) = 1 + jpzsize * dtaend / fixjpbub
         ELSE
            jpnbub(1) = jpbubend(1)
            limjpnz(1) = 1
         ENDIF
C
         IF (kflagincovyoz.EQ.4) THEN
            CALL readhdrzbas(kincovyoz,zon_jpi,zon_jpj,zon_jpk,zon_jpt,
     $                         jpbubend(2),jpz,jrbasdeb,jrbasfin)
            IF (jpbubend(2).GT.fixjpbub) THEN
               jpnbub(2) = fixjpbub
               limjpnz(2) = 1 + jpzsize * dtaend / fixjpbub
            ELSE
               jpnbub(2) = jpbubend(2)
               limjpnz(2) = 1
            ENDIF
         ENDIF
C
         IF (largoecorrel) THEN
            CALL readhdrzbas(argoecorrel,zon_jpi,zon_jpj,zon_jpk,zon_jpt,
     $                         jpbubend(3),jpz,jqbasdeb,jqbasfin)
            IF (jpbubend(3).GT.fixjpbub) THEN
               jpnbub(3) = fixjpbub
               limjpnz(3) = 1 + jpzsize * dtaend / fixjpbub
            ELSE
               jpnbub(3) = jpbubend(3)
               limjpnz(3) = 1
            ENDIF
         ENDIF
C
         IF (larginptzon) THEN
            CALL evalhdrzon(arginptzon,zon_jpi,zon_jpj,zon_jpk,zon_jpt,
     $           jpbubend(4),jpz)
            IF (jpbubend(4).GT.fixjpbub) THEN
               jpnbub(4) = fixjpbub
               limjpnz(4) = 1 + jpzsize * dtaend / fixjpbub
            ELSE
               jpnbub(4) = jpbubend(4)
               limjpnz(4) = 1
            ENDIF
         ENDIF
C
C --- allocation vectptbub
         allocate ( vectptbub(1:MAXVAL(jpnbub)), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectptbub(:) = 0
C
C --- allocation vectybub
         allocate ( vectybub(0:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectybub(0:) = FREAL(0.0)
C
         jpbsize=zon_jpi*zon_jpj*zon_jpk*zon_jpt*dtaend
C --- allocation ptlinbub
         allocate ( ptlinbub(1:jpbsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptlinbub(:) = 0
C --- allocation bubblklin
         allocate ( bubblklin(1:jpbsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         bubblklin(:) = FREAL(0.0)
C
C --- allocation bubblk1
         allocate ( bubblk1(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,
     $        1:jpnbub(1),1:1), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         bubblk1(:,:,:,:,:,:) = FREAL(0.0)
C --- allocation bubblk2
         IF (kflagincovyoz.EQ.4) THEN
            allocate ( bubblk2(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,
     $           1:jpnbub(2),1:jprsize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            bubblk2(:,:,:,:,:,:) = FREAL(0.0)
         ENDIF
C --- allocation bubblk3
         IF (largoecorrel) THEN
            allocate ( bubblk3(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,
     $           1:jpnbub(3),1:jpqsize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            bubblk3(:,:,:,:,:,:) = FREAL(0.0)
         ENDIF
C --- allocation bubblk4
         IF (larginptzon) THEN
            allocate ( bubblk4(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,
     $           1:jpnbub(4),1:1), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            bubblk4(:,:,:,:,:,:) = FREAL(0.0)
         ENDIF
C
C --- allocation ptbub
         allocate ( ptbub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,
     $        1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptbub(:,:,:,:,:) = 0
C
C --- allocation zone pointers
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
C
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
C
         allocate ( pt1bubidx(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         pt1bubidx(:,:) = 0
C
         IF (kflagincovyoz.EQ.4) THEN
            allocate ( pt2bubidx(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt2bubidx(:,:) = 0
         ENDIF
C
         IF (largoecorrel) THEN
            allocate ( pt3bubidx(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt3bubidx(:,:) = 0
         ENDIF
C
         IF (larginptzon) THEN
            allocate ( pt4bubidx(1:jpzsize,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            pt4bubidx(:,:) = 0
         ENDIF
C
         IF (kflagincovyoz.EQ.4) THEN
            CALL readptzbas (kincovyoz,pt2bubidx,pt1dtalon,pt1dtalat,
     $           pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth,
     $           pt1bubtime,jrbasdeb,jrbasfin)
         ENDIF
C
         IF (largoecorrel) THEN
            CALL readptzbas (argoecorrel,pt3bubidx,pt1dtalon,pt1dtalat,
     $           pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth,
     $           pt1bubtime,jqbasdeb,jqbasfin)
         ENDIF
C
         CALL readptzon (kinzon,pt1bubidx,pt1dtalon,pt1dtalat,
     $        pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth,
     $        pt1bubtime)
C
         IF (larginptzon) THEN
            CALL readptzon (arginptzon,pt4bubidx,pt1dtalon,pt1dtalat,
     $           pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth,
     $           pt1bubtime)
         ENDIF
C
      ENDIF
C
      jpu=jpssize
C
      IF (kflaggloloc.NE.0) THEN
         maxjpu=0
         nbzonnull=0
C --- allocation ptu
         allocate ( ptu(1:jpusize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptu(:) = FREAL(0.0)
C --- allocation ptubub
         IF (kflagcovyo.EQ.2) THEN
            allocate ( ptubub(1:jpusize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptubub(:) = FREAL(0.0)
         ENDIF
C --- allocation baseur
         allocate ( baseur(1:jpusize,1:jprsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         baseur(:,:) = FREAL(0.0)
C --- allocation vectuweight
         allocate ( vectuweight(1:jpusize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectuweight(:) = FREAL(0.0)
C --- allocation vectuinnov
         allocate ( vectuinnov(1:jpusize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectuinnov(:) = FREAL(0.0)
C --- allocation vectupart
         IF (larginpartobs) THEN
           allocate ( vectupart(1:jpusize), stat=allocok )
           IF (allocok.NE.0) GOTO 1001
           vectupart(:) = FREAL(0.0)
         ENDIF
C --- allocation poscoefuobs
         allocate ( poscoefuobs(1:jpusize,1:jpitpsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         poscoefuobs(:,:) = type_poscoef(0,FREAL(0.0))
C --- allocation poscoefbubobs
         IF (kflagcovyo.NE.2) THEN
            allocate (poscoefbubobs(1:jpusize,1:jpitpsize),stat=allocok)
            IF (allocok.NE.0) GOTO 1001
            poscoefbubobs(:,:) = type_poscoef(0,FREAL(0.0))
         ENDIF
C --- allocation baseuq
         IF (largoecorrel) THEN
            allocate ( baseuq(1:jpusize,1:jpqsize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            baseuq(:,:) = FREAL(0.0)
         ENDIF
      ENDIF
C
C -4.2- Loop on the blocks of subdomains and read the influence bubbles
C ---------------------------------------------------------------------
C
      DO jnz = 1,MAXVAL(limjpnz)
C
         jz1 = 1
         jz2 = jpzsize
C
         IF (kflaggloloc.NE.0) THEN
C
            IF (ANY(jpnbub(:).NE.jpbubend(:))) THEN
               jz1 = (jnz - 1) * fixjpbub / dtaend + 1
               jz2 = MIN(jnz*fixjpbub/dtaend,jpzsize)
               jpnbubsize = ( jz2 - jz1 + 1 ) * dtaend
            ENDIF
C
            zlectinfo = .FALSE.
            IF (jnz.EQ.1) zlectinfo = .TRUE.
C
            IF (jpnbub(1).EQ.jpbubend(1)) THEN
               IF (jnz.EQ.1) THEN
                  vectptbub(1:jpnbub(1)) = 
     $                 (/ (jbub, jbub = 1,jpnbub(1)) /)
                  CALL readnbubzon(kinzon,vectptbub(1:jpnbub(1)),
     $                 bubblk1(:,:,:,:,1:jpnbub(1),1),
     $                 zlectinfo)
               ENDIF
            ELSE
               vectptbub(1:jpnbubsize) = (/ ((pt1bubidx(jz,jdta),
     $              jdta=1,dtaend),jz=jz1,jz2) /)
               CALL readnbubzon(kinzon,vectptbub(1:jpnbubsize),
     $              bubblk1(:,:,:,:,1:jpnbubsize,1),
     $              zlectinfo)
            ENDIF
C
            IF (kflagincovyoz.EQ.4) THEN
               IF (jpnbub(2).EQ.jpbubend(2)) THEN
                  IF (jnz.EQ.1) THEN
                     vectptbub(1:jpnbub(2)) = 
     $                    (/ (jbub, jbub = 1,jpnbub(2)) /)
                     CALL readnbubzbas(kincovyoz,
     $                    bubblk2(:,:,:,:,1:jpnbub(2),:),
     $                    vectptbub(1:jpnbub(2)),
     $                    jrbasdeb,jrbasfin,zlectinfo)
                  ENDIF
               ELSE
                  vectptbub(1:jpnbubsize) = (/ ((pt2bubidx(jz,jdta),
     $                 jdta=1,dtaend),jz=jz1,jz2) /)
                  CALL readnbubzbas(kincovyoz,bubblk2(:,:,:,:,
     $                 1:jpnbubsize,:),
     $                 vectptbub(1:jpnbubsize),
     $                 jrbasdeb,jrbasfin,zlectinfo)
               ENDIF
            ENDIF
C
            IF (largoecorrel) THEN
               IF (jpnbub(3).EQ.jpbubend(3)) THEN
                  IF (jnz.EQ.1) THEN
                     vectptbub(1:jpnbub(3)) = 
     $                    (/ (jbub, jbub = 1,jpnbub(3)) /)
                     CALL readnbubzbas(argoecorrel,
     $                    bubblk3(:,:,:,:,1:jpnbub(3),:),
     $                    vectptbub(1:jpnbub(3)),
     $                    jqbasdeb,jqbasfin,zlectinfo)
                  ENDIF
               ELSE
                  vectptbub(1:jpnbubsize) = (/ ((pt3bubidx(jz,jdta),
     $                 jdta=1,dtaend),jz=jz1,jz2) /)
                  CALL readnbubzbas(argoecorrel,
     $                 bubblk3(:,:,:,:,1:jpnbubsize,:),
     $                 vectptbub(1:jpnbubsize),
     $                 jqbasdeb,jqbasfin,zlectinfo)
               ENDIF
            ENDIF
C
            IF (larginptzon) THEN
               IF (jpnbub(4).EQ.jpbubend(4)) THEN
                  IF (jnz.EQ.1) THEN
                     vectptbub(1:jpnbub(4)) = 
     $                    (/ (jbub, jbub = 1,jpnbub(4)) /)
                     CALL readnbubzon(arginptzon,vectptbub(1:jpnbub(4)),
     $                    bubblk4(:,:,:,:,1:jpnbub(4),1),zlectinfo)
                  ENDIF
               ELSE
                  vectptbub(1:jpnbubsize) = (/ ((pt4bubidx(jz,jdta),
     $                 jdta=1,dtaend),jz=jz1,jz2) /)
                  CALL readnbubzon(arginptzon,vectptbub(1:jpnbubsize),
     $                 bubblk4(:,:,:,:,1:jpnbubsize,1),zlectinfo)
               ENDIF
            ENDIF
         ENDIF
C
C -4.3- Loop on the subdomains
C ----------------------------
C define the obs influence zone for current subdomain :
C -----------------------------------------------------
C
         IF (nprint.GE.1) 
     $      print *,'--> Computation of ROA coefficients, nproc=',jpproc
         IF (nprint.GE.1) 
     $      print *, '  Computing for zones in the range:',jz1,jz2
         DO jz = jz1+jproc,jz2,jpproc
C
            IF (kflaggloloc.NE.0) THEN
C
               IF (larginptzon) THEN
                  DO jdta=1,dtaend
                     IF (jpnbub(4).EQ.jpbubend(4)) THEN
                        jbub = pt4bubidx(jz,jdta)
                     ELSE
                        jbub = (jz-jz1) * dtaend + jdta
                     ENDIF
                     ptbub(:,:,:,:,jdta) = NINT(bubblk4(:,:,:,:,jbub,1))
                  ENDDO
               ELSE
                  DO jdta=1,dtaend
                     CALL mkptbub (ptbub(:,:,:,:,jdta),jdta,
     $                    pt1dtalon(jz,jdta),pt1dtalat(jz,jdta),
     $                    pt1dtadepth(jz,jdta),pt1dtatime(jz,jdta),
     $                    pt1bublon(jz,jdta),pt1bublat(jz,jdta),
     $                    pt1bubdepth(jz,jdta),pt1bubtime(jz,jdta))
                  ENDDO
               ENDIF
C
               jb=0
               DO jdta=1,dtaend
                 DO jt=1,zon_jpt
                 DO jk=1,zon_jpk
                 DO jj=1,zon_jpj
                 DO ji=1,zon_jpi
                   jb=jb+1
                   ptlinbub(jb)=ptbub(ji,jj,jk,jt,jdta)
                 ENDDO
                 ENDDO
                 ENDDO
                 ENDDO
               ENDDO
C              ptlinbub(1:jpbsize)=
C    $              TRANSFER(ptbub(:,:,:,:,:),INT(1),jpbsize)
C
               jb=0
               DO jdta=1,dtaend
                 IF (jpnbub(1).EQ.jpbubend(1)) THEN
                   jbub = pt1bubidx(jz,jdta)
                 ELSE
                   jbub = (jz-jz1) * dtaend + jdta
                 ENDIF
C     
                 DO jt=1,zon_jpt
                 DO jk=1,zon_jpk
                 DO jj=1,zon_jpj
                 DO ji=1,zon_jpi
                   jb=jb+1
                   bubblklin(jb)=bubblk1(ji,jj,jk,jt,jbub,1)
                 ENDDO
                 ENDDO
                 ENDDO
                 ENDDO
C                bubblklin((1+(jdta-1)*jpbsize/dtaend)
C    $                 :(jdta*jpbsize/dtaend))=
C    $                 TRANSFER(bubblk1(:,:,:,:,jbub,1),
C    $                 FREAL(1.0),jpbsize/dtaend)
               ENDDO
 
               DO jb=1,jpbsize
                  bubblklin(jb)=
     $                 FREAL(COUNT(bubblklin(jb:jb).NE.FREAL(0.0)))
               ENDDO
 
               ptlinbub(1:jpbsize) =
     $              ptlinbub(1:jpbsize) * NINT(bubblklin(1:jpbsize))
 
               vectybub(0:jpysize) = FREAL(0.0)
               DO jb=1,jpbsize
                  vectybub(ptlinbub(jb))=bubblklin(jb)
               ENDDO
C
C              WHERE (bubblklin(:).NE.0.0_kr) vectybub(ptlinbub(:))=1.0_kr
               IF (largconnect) CALL mkconnecty(vectybub(0:jpysize))
C
C               DO jdta=1,dtaend
C                  IF (jpnbub(1).EQ.jpbubend(1)) THEN
C                     jbub = pt1bubidx(jz,jdta)
C                  ELSE
C                     jbub = (jz-jz1) * dtaend + jdta
C                  ENDIF
C
C                  DO jt=1,zon_jpt
C                  DO jk=1,zon_jpk
C                  DO jj=1,zon_jpj
C                  DO ji=1,zon_jpi
C                     vectybub(ptbub(ji,jj,jk,jt,jdta))=
C     $                    FREAL(COUNT(bubblk1(ji:ji,jj,jk,jt,jbub,1).NE.FREAL(0.0)))
C                  ENDDO
C                  ENDDO
C                  ENDDO
C                  ENDDO
C               ENDDO
C
               IF (kflagcovyo.EQ.2) THEN
                  vectsweight(1:jpssize)=vectybub(1:jpysize)
               ELSE
                  CALL mkhytoocompt(vectybub(1:jpysize),
     $                 vectsweight(1:jpssize),poscoefobs(:,:))
               ENDIF
C
               jpu = COUNT(vectsweight(1:jpssize).NE.FREAL(0.0))
C
               IF (jpu.EQ.0) THEN
                  nbzonnull = nbzonnull + 1
               ELSE
                  maxjpu=MAX(maxjpu,jpu)
                  IF (jpu.GT.jpusize) THEN
C --- desallocate some arrays
                     IF (allocated(ptu)) deallocate(ptu)
                     IF (allocated(ptubub)) deallocate(ptubub)
                     IF (allocated(baseur)) deallocate(baseur)
                     IF (allocated(baseuq)) deallocate(baseuq)
                     IF (allocated(vectuweight))  deallocate(vectuweight)
                     IF (allocated(vectuinnov)) deallocate(vectuinnov)
                     IF (allocated(poscoefuobs)) deallocate(poscoefuobs)
                     IF (allocated(poscoefbubobs))  deallocate(poscoefbubobs)
C
                     jpusize=jpu
C --- allocation ptu
                     allocate ( ptu(1:jpusize), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
C --- allocation ptubub
                     IF (kflagcovyo.EQ.2) THEN
                        allocate ( ptubub(1:jpusize), stat=allocok )
                        IF (allocok.NE.0) GOTO 1001
                     ENDIF
C --- allocation baseur
                     allocate ( baseur(1:jpusize,1:jprsize), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
C --- allocation vectuweight
                     allocate ( vectuweight(1:jpusize), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
C --- allocation vectuinnov
                     allocate ( vectuinnov(1:jpusize), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
C --- allocation poscoefuobs
                     allocate ( poscoefuobs(1:jpusize,1:jpitpsize),
     $                    stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
C --- allocation poscoefbubobs
                     IF (kflagcovyo.NE.2) THEN
                        allocate ( poscoefbubobs(1:jpusize,1:jpitpsize),
     $                       stat=allocok )
                        IF (allocok.NE.0) GOTO 1001
                     ENDIF
C --- allocation baseuq
                     IF (largoecorrel) THEN
                        allocate ( baseuq(1:jpusize,1:jpqsize),
     $                       stat=allocok )
                        IF (allocok.NE.0) GOTO 1001
                     ENDIF
                  ENDIF
C                 
                  ju=0
                  DO js=1,jpssize
                    IF (vectsweight(js).NE.FREAL(0.0)) THEN
                      ju=ju+1
                      ptu(ju)=js
                    ENDIF
                  ENDDO
                  IF (ju.NE.jpu) GOTO 1000
C                 ptu(:jpu) = PACK( (/ (js, js=1,jpssize) /) ,
C    $                        vectsweight(:) .NE. FREAL(0.0) )
                  baseur(:jpu,:) = FREAL(0.0)
                  vectuweight(:jpu) = FREAL(0.0)
                  vectuinnov(:jpu) = FREAL(0.0)
                  poscoefuobs(1:jpu,1:jpitpsize) = 
     $                 poscoefobs(ptu(1:jpu),1:jpitpsize)
                  IF (largoecorrel) THEN
                     baseuq(:jpu,:) = FREAL(0.0)
                  ENDIF
               ENDIF
            ENDIF
C
            lmodprint=(MOD(jz-1,(jpzsize/1000+1)).EQ.0)
            IF (lmodprint) print *, 'Zone=',jz,'Size=',jpu,'Proc=',jproc
C
            IF (jpu.EQ.0) THEN 
               matCUrrz(:,:,jz)=FREAL(0.0)
               DO jr=jrmatdeb,jrmatfin
                  matCUrrz(jr,jr,jz)=FREAL(1.0)
               ENDDO
               coefrz(:,jz) = FREAL(0.0)
               IF (largcoefrmax) THEN
                  coefrznolim(:,jz) = FREAL(0.0)
               ENDIF
               lambda(:,jz) = FREAL(0.0)
               delta(:,jz) = FREAL(0.0)
            ELSE
C
C -4.32- Compute pointer for line bub :
C -------------------------------------
C
               SELECT CASE (kflaggloloc)
               CASE (0)
C --- nothing
               CASE (2,3,4)
                  vectybub(0:jpysize) = FREAL(0.0)
                  DO jb=1,jpbsize
                     vectybub(ptlinbub(jb))=FREAL(jb)
                  ENDDO 
C
                  IF (largconnect) CALL mkconnecty(vectybub(0:jpysize))
C
                  IF (kflagcovyo.EQ.2) THEN
                     DO ju=1,jpusize
                        ptubub(ju) = NINT(vectybub(ptu(ju)))
                     ENDDO
                  ELSE
                     poscoefbubobs(1:jpu,1:jpitpsize)%coef = 
     $                 poscoefuobs(1:jpu,1:jpitpsize)%coef
                     DO jitp=1,jpitpsize                        
                     DO ju=1,jpusize
                        poscoefbubobs(ju,jitp)%pos = 
     $                    NINT(vectybub(poscoefuobs(ju,jitp)%pos))
                     ENDDO
                     ENDDO
                  ENDIF
               CASE DEFAULT
                  GOTO 1000
               END SELECT
C
C -4.33- Load part of uweight :
C -----------------------------
C
               SELECT CASE (kflaggloloc)
               CASE (0)
C --- nothing
               CASE (2,3,4)
C
                  jb=0
                  DO jdta=1,dtaend
                    IF (jpnbub(1).EQ.jpbubend(1)) THEN
                      jbub = pt1bubidx(jz,jdta)
                    ELSE
                      jbub = (jz-jz1) * dtaend + jdta
                    ENDIF
C     
                    DO jt=1,zon_jpt
                    DO jk=1,zon_jpk
                    DO jj=1,zon_jpj
                    DO ji=1,zon_jpi
                      jb=jb+1
                      bubblklin(jb)=bubblk1(ji,jj,jk,jt,jbub,1)
                    ENDDO
                    ENDDO
                    ENDDO
                    ENDDO
C                   bubblklin((1+(jdta-1)*jpbsize/dtaend)
C    $                    :(jdta*jpbsize/dtaend))=
C    $                    TRANSFER(bubblk1(:,:,:,:,jbub,1),
C    $                    FREAL(1.0),jpbsize/dtaend)
                  ENDDO
C
                  IF (kflagcovyo.EQ.2) THEN
                     vectuweight(1:jpu)=bubblklin(ptubub(1:jpu))
                  ELSE
                     CALL mkhytou(bubblklin(1:jpbsize),
     $                 vectuweight(:jpu),poscoefbubobs(:jpu,:))
                  ENDIF
C
               CASE DEFAULT
                  GOTO 1000
               END SELECT
C
C -4.4- Load part of A for covariance computation :
C -------------------------------------------------
C
               SELECT CASE (kflaggloloc)
               CASE (0)
C --- nothing
               CASE (2,3)
                  baseur(1:jpu,1:jprsize)=basesr(ptu(1:jpu),1:jprsize)
               CASE (4)
                  zlectinfo = .FALSE.
                  IF (kflagincovyoz.NE.4) GOTO 1000
C
                  DO jr = jrbasdeb,jrbasfin
                     jb=0
                     DO jdta=1,dtaend
                        IF (jpnbub(2).EQ.jpbubend(2)) THEN
                           jbub = pt2bubidx(jz,jdta)
                        ELSE
                           jbub = (jz-jz1) * dtaend + jdta
                        ENDIF
                        DO jt=1,zon_jpt
                        DO jk=1,zon_jpk
                        DO jj=1,zon_jpj
                        DO ji=1,zon_jpi
                          jb=jb+1
                          bubblklin(jb)=bubblk2(ji,jj,jk,jt,jbub,jr)
                        ENDDO
                        ENDDO
                        ENDDO
                        ENDDO
C                       bubblklin((1+(jdta-1)*jpbsize/dtaend)
C    $                    :(jdta*jpbsize/dtaend))=
C    $                    TRANSFER(bubblk2(:,:,:,:,jbub,jr),
C    $                    FREAL(1.0),jpbsize/dtaend)
                     ENDDO
                     IF (kflagcovyo.EQ.2) THEN
                        baseur(1:jpu,jr)=bubblklin(ptubub(1:jpu))
                     ELSE
                        CALL mkhytou(bubblklin(1:jpbsize),
     $                    baseur(:jpu,jr),poscoefbubobs(:jpu,:))
                     ENDIF
                  ENDDO
               CASE DEFAULT
                  GOTO 1000
               END SELECT
C
C -4.5- Load the inverse observation error correlation matrix 
C -----------------------------------------------------------
C
               IF (largoecorrel) THEN
                  IF (kflaggloloc.EQ.0) GOTO 1000
C 
                  DO jq = jqbasdeb,jqbasfin
                     jb=0
                     DO jdta=1,dtaend
                        IF (jpnbub(3).EQ.jpbubend(3)) THEN
                           jbub = pt3bubidx(jz,jdta)
                        ELSE
                           jbub = (jz-jz1) * dtaend + jdta
                        ENDIF
                        DO jt=1,zon_jpt
                        DO jk=1,zon_jpk
                        DO jj=1,zon_jpj
                        DO ji=1,zon_jpi
                          jb=jb+1
                          bubblklin(jb)=bubblk3(ji,jj,jk,jt,jbub,jq)
                        ENDDO
                        ENDDO
                        ENDDO
                        ENDDO
C                       bubblklin((1+(jdta-1)*jpbsize/dtaend)
C    $                    :(jdta*jpbsize/dtaend))=
C    $                    TRANSFER(bubblk3(:,:,:,:,jbub,jq),
C    $                    FREAL(1.0),jpbsize/dtaend)
                     ENDDO
                     IF (kflagcovyo.EQ.2) THEN
                        baseuq(1:jpu,jq)=bubblklin(ptubub(1:jpu))
                     ELSE
                        CALL mkhytou(bubblklin(1:jpbsize),
     $                    baseuq(:jpu,jq),poscoefbubobs(:jpu,:))
C     
                     ENDIF
                  ENDDO
               ENDIF
C
C -4.6- Weight (HSf) and inv(R) using the local bulk parameterization
C -------------------------------------------------------------------
C
               IF (kflaggloloc.EQ.0) THEN
                  vectsweight(:)=vectssqrdiagRi(:)
               ELSE
C --------------- New localization scheme (as in MWR2011 paper)
                  IF (lbulkave) THEN
C ----------------- Weight (HSf) using (Ws)^1/2
                    DO jr = jrbasdeb,jrbasfin
                      baseur(1:jpu,jr) = baseur(1:jpu,jr)
     $                      * SQRT( vectuweight(1:jpu) )
                    ENDDO
C ----------------- Compute (Wr) from (Ws) [vectuinnov = temporary storage]
                    DO ju = 1,jpu
                      vectuinnov(ju) = DOT_PRODUCT(
     $                                  baseur(ju,jrbasdeb:jrbasfin),
     $                                  baseur(ju,jrbasdeb:jrbasfin) )
     $                               / vectuweight(ju)
     $                               * vectssqrdiagRi(ptu(ju))
     $                               * vectssqrdiagRi(ptu(ju))
                    ENDDO
                    ratio = EXP ( SUM ( LOG( vectuinnov(1:jpu) ) )
     $                          / FREAL(jpu) )
C
                    DO ju = 1,jpu
                      vectuweight(ju) = FREAL(1.0) + ratio *
     $                                  ( FREAL(1.0) - vectuweight(ju) )
                    ENDDO
C ----------------- Original algorithm to compute Wr
C ----------------- [do not work with non diagonal R]
C                   DO ju = 1,jpu
C                     vectuweight(ju) = ( FREAL(1.0) - vectuweight(ju) )
C    $                   * DOT_PRODUCT( baseur(ju,jrbasdeb:jrbasfin),
C    $                                  baseur(ju,jrbasdeb:jrbasfin) )
C    $                     / vectuweight(ju)
C    $                     * vectssqrdiagRi(ptu(ju))
C    $                     * vectssqrdiagRi(ptu(ju))
C    $                     + FREAL(1.0)
C                   ENDDO
C ----------------- Weight inv(R) using (Wr)^1/2
                    vectuweight(1:jpu) = vectssqrdiagRi(ptu(1:jpu))
     $                                   / SQRT( vectuweight(1:jpu) )
C --------------- Old localization scheme (as in JGR2003 paper)
                  ELSE
                    vectuweight(1:jpu) = vectuweight(1:jpu)
     $                * vectssqrdiagRi(ptu(1:jpu))
                  ENDIF
C
                  vectuinnov(1:jpu) = vectsinnov(ptu(1:jpu))
                  IF (larginpartobs) THEN
                    vectupart(1:jpu) = vectspart(ptu(1:jpu))
                  ENDIF
               ENDIF
C
C -4.8.1- Compute SVD decomposition of the local ROA kernel matrix 
C -----------------------------------------------------------------
C         (U,lambda) = SVD ( trsp(HSf) inv(R) (HSf) )
C
               forgfact = factoubli
               IF (.NOT.larginpartobs) THEN
                 IF (kflaggloloc.EQ.0) THEN
                   CALL algoker_u(basesr(:,:),vectsweight(:),
     $                            matCUrrz(:,:,jz),lambda(:,jz))
                 ELSE
                    IF (largoecorrel) THEN
                       CALL algoker_u(baseur(:jpu,:),vectuweight(:jpu),
     $                       matCUrrz(:,:,jz),lambda(:,jz),
     $                       kbasesq=baseuq(:jpu,:),kmatBqr=matBqr(:,:))
                    ELSE
                       CALL algoker_u(baseur(:jpu,:),vectuweight(:jpu),
     $                                matCUrrz(:,:,jz),lambda(:,jz))
                    ENDIF
                 ENDIF
               ELSE
                 IF (kflaggloloc.EQ.0) THEN
                   CALL algocalcGi(basesr(:,:),vectsweight(:),
     $                             vectspart(:),matGrrzi(:,:,jz,:))
                 ELSE
                   CALL algocalcGi(baseur(:jpu,:),vectuweight(:jpu),
     $                             vectupart(:jpu),matGrrzi(:,:,jz,:))
                 ENDIF
                 IF (larginparadap) THEN
                   CALL algoker_Gi(matGrrzi(:,:,jz,:),
     $                             matCUrrz(:,:,jz),lambda(:,jz),
     $                             beta=vectbeta(jz,:))
                 ELSE
                   CALL algoker_Gi(matGrrzi(:,:,jz,:),
     $                             matCUrrz(:,:,jz),lambda(:,jz))
                 ENDIF
               ENDIF
C
C -4.8.2- Project the innovation vector on the ROA modes
C -------------------------------------------------------
C            c = U^ trsp(HSf) inv(R) (y-Hxf)
C
               IF (.NOT.larginpartobs) THEN
                 IF (kflaggloloc.EQ.0) THEN
                   CALL algocoef_u(basesr(:,:),vectsinnov(:),
     $                   vectsweight(:),matCUrrz(:,:,jz),delta(:,jz))
                   IF (largoutrz) CALL algocalciRi(vectsinnov(:),
     $                 vectsweight(:),iRi_zi(jz,:),y_zi(jz,:))
                 ELSE
                   IF (largoecorrel) THEN
                     CALL algocoef_u(baseur(:jpu,:),vectuinnov(:jpu),
     $                   vectuweight(:jpu),matCUrrz(:,:,jz),delta(:,jz),
     $                   kbasesq=baseuq(:jpu,:),kmatBqr=matBqr(:,:))
                   ELSE
                     CALL algocoef_u(baseur(:jpu,:),vectuinnov(:jpu),
     $                   vectuweight(:jpu),matCUrrz(:,:,jz),delta(:,jz))
                     IF (largoutrz) CALL algocalciRi(vectuinnov(:jpu),
     $                   vectuweight(:jpu),iRi_zi(jz,:),y_zi(jz,:))
                   ENDIF
                 ENDIF
               ELSE
                 IF (kflaggloloc.EQ.0) THEN
                   CALL algocalcDi(basesr(:,:),vectsinnov(:),
     $                             vectsweight(:),vectspart(:),
     $                             deltarzi(:,jz,:))
                   IF (largoutrz) CALL algocalciRi(vectsinnov(:),
     $                 vectsweight(:),iRi_zi(jz,:),y_zi(jz,:),
     $                 kvectspart=vectspart(:))
                 ELSE
                   CALL algocalcDi(baseur(:jpu,:),vectuinnov(:jpu),
     $                             vectuweight(:jpu),vectupart(:jpu),
     $                             deltarzi(:,jz,:))
                   IF (largoutrz) CALL algocalciRi(vectuinnov(:jpu),
     $                 vectuweight(:jpu),iRi_zi(jz,:),y_zi(jz,:),
     $                 kvectspart=vectspart(:))
                 ENDIF
                 IF (larginparadap) THEN
                   CALL algosum_Di(deltarzi(:,jz,:),delta(:,jz),
     $                            beta=vectbeta(jz,:))
                 ELSE
                   CALL algosum_Di(deltarzi(:,jz,:),delta(:,jz))
                 ENDIF
               ENDIF
C     
C -4.8.3- Scale forecast and observation error covariance matrix
C --------------------------------------------------------------
C
               IF (larginparadap) THEN
C
                 delta(:,jz) = delta(:,jz) * SQRT(vectalpha(jz))
                 lambda(:,jz) = lambda(:,jz) * vectalpha(jz)
                 IF (.NOT.larginpartobs) THEN
                   delta(:,jz) = delta(:,jz) / vectbeta(jz,1)
                   lambda(:,jz) = lambda(:,jz) / vectbeta(jz,1)
                 ENDIF
C
               ENDIF
C
C -4.8.4- Compute the ROA coefficients
C -------------------------------------
C
               CALL algocoef_u_sf(coefrz(:,jz),lambda(:,jz),delta(:,jz),
     $                            matCUrrz(:,:,jz),forgfact)
C
               IF (lmodprint) THEN
               IF (nprint.GE.3) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOROA_U: ',
     $              'computing the ROA coefficients (c)'
                  titre='ALGOROA_U: Xa-ROA coefficients'
                  variable='coefrz'
                  CALL printtab_r (coefrz(jrbasdeb:jrbasfin,jz),
     $              titre,variable,jrbasdeb)
               ENDIF
               ENDIF
C
C -4.8.5- Saturate ROA coefficients to coefrmax
C ------------------------------------------------
C
               IF ((largcoefrmax).AND.(coefrmax.GT.FREAL(0.0))) THEN
C
                  coefrznolim(:,jz)=coefrz(:,jz)
                  CALL algocoeflimit_u(coefrz(:,jz),matCUrrz(:,:,jz),
     $              coefrmax,forgfact)
C     
                  IF (lmodprint) THEN
                  IF (nprint.GE.3) THEN
                     titre='ALGOROA_U: saturated Xa-ROA coefficients'
                     variable='coefrz'
                     CALL printtab_r (coefrz(jrbasdeb:jrbasfin,jz),
     $                 titre,variable,jrbasdeb)
                  ENDIF
                  ENDIF
C
               ENDIF
C
C -4.8.6- Compute Sf to Sa transformation matrix
C ----------------------------------------------
C
               CALL algomatc_u(lambda(:,jz),matCUrrz(:,:,jz),forgfact)
C     
               IF (lmodprint) THEN
               IF (nprint.GE.3) THEN
                  titre='ALGOROA_U: Pa-ROA coefficients'
                  variable='matCUrr'
                  CALL printtab_rr (
     $              matCUrrz(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin,jz),
     $              titre,variable,jrmatdeb,jrmatdeb)
               ENDIF
               ENDIF
C
C -4.8.7- Rescale lambda and delta for output in crz file
C -------------------------------------------------------
C
               IF (largoutrz) THEN
                 IF ((larginparadap).AND.(.NOT.larginpartobs)) THEN
                   delta(2:,jz)=delta(2:,jz)
     $                      /SQRT(vectalpha(jz))*vectbeta(jz,1)
                   lambda(2:,jz)=vectalpha(jz)
     $                      /( vectbeta(jz,1)*lambda(2:,jz) )
                 ELSE
                   lambda(2:,jz)=1.0_kr/lambda(2:,jz)
                 ENDIF
                 lambda(2:,jz)=lambda(2:,jz)*forgfact
                 delta(2:,jz)=delta(2:,jz)/SQRT(forgfact)
               ENDIF
C    
            ENDIF
C
         ENDDO
      ENDDO
C
      IF (kflaggloloc.NE.0) THEN
         print *,'maxjpu=',maxjpu,' nbzonnull=',nbzonnull,' %=',
     $        FREAL(100.0*nbzonnull/jpz)
         IF (nprint.EQ.2) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOROA : maxjpu=',maxjpu
         ENDIF
C --- desallocate some arrays
         IF (allocated(ptu))  deallocate (ptu)
         IF (allocated(ptlinbub))  deallocate (ptlinbub)
         IF (allocated(bubblklin))  deallocate (bubblklin)
         IF (allocated(baseur))  deallocate (baseur)
         IF (allocated(baseuq))  deallocate (baseuq)
         IF (allocated(vectuweight))  deallocate (vectuweight)
         IF (allocated(vectuinnov))  deallocate (vectuinnov)
         IF (allocated(poscoefuobs))  deallocate (poscoefuobs)
         IF (allocated(poscoefbubobs))  deallocate (poscoefbubobs)
      ENDIF
C
#if defined MPI
      DO jz=1,jpzsize
         jproc1=mod(jz-1,jpproc)
         CALL mpi_bcast(coefrz(1,jz),jprsize,mpi_double_precision,
     $                  jproc1,mpi_comm_world,mpi_code)
         CALL mpi_bcast(matCUrrz(1,1,jz),jprsize*jprsize,
     $                  mpi_double_precision,
     $                  jproc1,mpi_comm_world,mpi_code)
         IF (largoutrz) THEN
C
           IF (larginpartobs) THEN
             DO jiobs=1,jpiobs
               CALL mpi_bcast(deltarzi(1,jz,jiobs),jprsize,
     $                        mpi_double_precision,
     $                        jproc1,mpi_comm_world,mpi_code)
               CALL mpi_bcast(matGrrzi(1,1,jz,jiobs),jprsize*jprsize,
     $                        mpi_double_precision,
     $                        jproc1,mpi_comm_world,mpi_code)
             ENDDO
           ELSE
             CALL mpi_bcast(delta(1,jz),jprsize,mpi_double_precision,
     $                      jproc1,mpi_comm_world,mpi_code)
             CALL mpi_bcast(lambda(1,jz),jprsize,mpi_double_precision,
     $                      jproc1,mpi_comm_world,mpi_code)
           ENDIF
C
           DO jiobs=1,jpiobs
             CALL mpi_bcast(iRi_zi(jz,jiobs),1,mpi_double_precision,
     $                      jproc1,mpi_comm_world,mpi_code)
             CALL mpi_bcast(y_zi(jz,jiobs),1,mpi_double_precision,
     $                      jproc1,mpi_comm_world,mpi_code)
           ENDDO
C
         ENDIF
      ENDDO
#endif
C
C -5.- Store information on local analyses
C -----------------------------------------
C
C Store information needed by adaptive scheme
      IF ((largoutrz).AND.(jproc.EQ.0)) THEN
C
         IF (nprint.GE.1) THEN
            print *,'--> Storing (coef,lambda,U) ...'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOROA : ',
     $           'storing reduced order elements'
         ENDIF
C
         IF (larginpartobs) THEN
           CALL writemat3crz(argoutrz,'Delta',
     $           deltarzi(2:jprsize,1:jpzsize,:),'rzi')
           CALL writemat4crz(argoutrz,'Gamma',
     $           matGrrzi(1:jprsize-1,1:jprsize-1,1:jpzsize,:),'rrzi')
         ELSE
           CALL writemat2crz(argoutrz,'delta',
     $           delta(2:jprsize,1:jpzsize),'rz')
           CALL writemat2crz(argoutrz,'Lambda',
     $           lambda(2:jprsize,1:jpzsize),'rz')
         ENDIF
C
         CALL writemat2crz(argoutrz,'iRi',iRi_zi(1:jpzsize,:),'zi')
         CALL writemat2crz(argoutrz,'y',y_zi(1:jpzsize,:),'zi')
C
      ENDIF
C
C Compute and write tr(HK) for each local analysis
      IF ((largtypedtadiag).AND.(jproc.EQ.0)) THEN
C
C --- allocation vecty
         allocate ( vecty(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vecty(:) = FREAL(0.0)
C
C --- read partition in vecty
         IF (kflaggloloc.NE.0) THEN
            flagxyo = 2
            lectinfo = .TRUE.
            jnxyo = 1
            lmoyectold=lmoyect
            lmoyect=.FALSE.
            CALL readxyo (kinpartxyo,vecty(:),jnxyo,
     $           lectinfo,flagxyo,poscoefobs(:,:))
            lmoyect=lmoyectold
         ELSE
            vecty(:)=FREAL(1.0)
         ENDIF
C
C --- compute tr(HK) for each local analysis
         DO jz=1,jpzsize
            lambda(2:,jz) = lambda(2:,jz) / (1. + lambda(2:,jz) )
            lambda(1,jz) = SUM(lambda(2:,jz)) / (jprsize-1)
         ENDDO
C
C --- write tr(HK)
         WRITE(fnameout,'("trHK_",A)')
     $        argtypedtadiag(1:lenv(argtypedtadiag))
         CALL mkcoeftodta(fnameout,lambda(1,:),vecty(:))
C
      ENDIF
C
C -6.- Allocation/desallocation
C -----------------------------
C
      IF ((.NOT.(lXaout)).AND.(.NOT.(lPaout))) RETURN
C --- desallocate some arrays
      IF (allocated(iRi_zi)) deallocate (iRi_zi)
      IF (allocated(y_zi)) deallocate (y_zi)
      IF (allocated(matGrrzi)) deallocate (matGrrzi)
      IF (allocated(deltarzi)) deallocate (deltarzi)
      IF (allocated(vectalpha)) deallocate (vectalpha)
      IF (allocated(vectbeta)) deallocate (vectbeta)
      IF (allocated(vectspart)) deallocate (vectspart)
      IF (allocated(vectupart)) deallocate (vectupart)
C
      IF (allocated(vectsweight)) deallocate (vectsweight)
      IF (allocated(vectsinnov))  deallocate (vectsinnov)
      IF (allocated(vectssqrdiagRi)) deallocate (vectssqrdiagRi)
      IF (allocated(lambda)) deallocate (lambda)
      IF (allocated(delta)) deallocate (delta)
      IF (allocated(matBqr)) deallocate (matBqr)
      IF (allocated(coefrznolim)) deallocate (coefrznolim)
C
      IF (allocated(pt3bubidx)) deallocate (pt3bubidx)
      IF (allocated(pt4bubidx)) deallocate (pt4bubidx)

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
      IF (allocated(vectptbub)) deallocate (vectptbub)
C
      IF (allocated(ptbub)) deallocate (ptbub)
      IF (allocated(vectybub)) deallocate (vectybub)
      IF (allocated(bubblk1)) deallocate (bubblk1)
      IF (allocated(bubblk2)) deallocate (bubblk2)
      IF (allocated(bubblk3)) deallocate (bubblk3)
      IF (allocated(bubblk4)) deallocate (bubblk4)
C
      IF (allocated(basesr))  deallocate (basesr)
C
C -9.- Compute the analysed error covariance matrix (Sa)
C       and the analysed state (Xa)
C ------------------------------------------------------
C --- allocate an additional state vector
C --- allocation mat_yr
      allocate(mat_yr(1:MIN(jpysize,jpxsize),1:jprsize),stat=allocok)
      IF (allocok.NE.0) GOTO 1001
      mat_yr(:,:) = FREAL(0.0)
C --- allocation basexr
      allocate ( basexr(1:jpxsize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basexr(:,:) = FREAL(0.0)
C --- allocation vectxa
      allocate ( vectxa(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxa(:) = FREAL(0.0)
C --- allocation vectxf
      allocate ( vectxf(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxf(:) = FREAL(0.0)
C --- allocation vectxpart
      allocate ( vectxpart(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxpart(:) = FREAL(0.0)
C --- allocation vectxweight
      allocate ( vectxweight(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectxweight(:) = FREAL(0.0)
C
      IF (nprint.GE.1) print *,
     $     '--> Computation of (Xa,Pa) ...'
      lectinfo = .FALSE.
      IF (jpproc.GT.limjpnxyo(kflaganlxyo)) GOTO 1003
      DO jnxyo=1+jproc,limjpnxyo(kflaganlxyo),jpproc
         lmodprint=(MOD(jnxyo-1,(limjpnxyo(kflaganlxyo)/5+1)).EQ.0)
         IF ((lmodprint).AND.(nprint.GE.1))
     $        print *,'Memory part number : ',
     $           jnxyo,'/',limjpnxyo(kflaganlxyo)
C
C -9.1- Read the partition in subdomains
C --------------------------------------
         IF (kflaggloloc.NE.0) THEN
C
            IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOROA: ',
     $              'reading the partition in subdomains'
            ENDIF
            IF (kflaganlxyo.NE.3) THEN
               lmoyectold=lmoyect
               lmoyect=.FALSE.
               CALL readxyo (kinpartxyo,vectxpart(:),jnxyo,
     $              lectinfo,kflaganlxyo,poscoefobs(:,:))
               lmoyect=lmoyectold
            ELSE
C --- allocation vecty
               allocate ( vecty(1:jpysize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               vecty(:) = FREAL(0.0)
               flagxyo=2
               lmoyectold=lmoyect
               lmoyect=.FALSE.
               CALL readxyo (kinpartxyo,vecty(:),jnxyo,
     $              lectinfo,flagxyo,poscoefobs(:,:))
               lmoyect=lmoyectold
               vectxpart(:)=vecty(poscoefobs(:,1)%pos)
               IF (allocated(vecty))  deallocate (vecty)
            ENDIF
C
         ENDIF
C
C -9.2- Read the input forecast covariance matrix Sf
C ---------------------------------------------------
C
         IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
C          print *,'-9.2- Read the forecast input basis Sf'
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOROA : ',
     $              'reading the forecast input basis Sf'
         ENDIF
         CALL readbas (kinbasxyo,basexr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $        lectinfo,kflaganlxyo,poscoefobs(:,:))
C
         IF (valbase.LT.1) THEN
           DO jx=1,jpxsize
             vectxf(jx) = SUM(basexr(jx,:))/jprsize
             basexr(jx,:) = basexr(jx,:) - vectxf(jx)
             basexr(jx,:) = basexr(jx,:) / SQRT(FREAL(jprsize-1))
           ENDDO
         ENDIF
C
C -9.3- Compute the analysed state (xa - xf)
C -------------------------------------------
C
         IF ((lXaout).OR.(valbase.LT.1)) THEN
            IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
C            print *,'-9.3.- Compute the analysed state (xa-xf)'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOROA : ',
     $              'computing the analysed state (xa-xf)'
            ENDIF
            vectxa(:) = FREAL(0.0)
            IF (kflaggloloc.EQ.0) THEN
               DO jr=jrbasdeb,jrbasfin
                  vectxa(:) = vectxa(:)+
     $                 basexr(:,jr)*coefrz(jr,1)
               ENDDO
            ELSE
               DO jr=jrbasdeb,jrbasfin
                  DO jx=1,jpxsize
                     vectxa(jx)=vectxa(jx)+
     $                 basexr(jx,jr)*coefrz(jr,nint(vectxpart(jx)))
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
C
C -9.4- Read the weight from the reducevar config
C -----------------------------------------------
C
         IF ((largreducevar).AND.(lXaout)) THEN
C
           IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
C              print *,' -9.4- Read the weight from the reducevar config'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOROA : ',
     $              'reading the weight from the reducevar config'
            ENDIF
            CALL readxyo (argreducevar,vectxweight(:),
     $           jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
            vectxa(:) = vectxweight(:)*vectxa(:)
C
         ENDIF
C
C -9.5- Read the forecast state vector (xf) 
C -----------------------------------------
C
         IF (((lXaout).AND.(lXa)).OR.((valbase.LT.1).AND.(lPaout))) THEN
C
            IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
C     print *,'-9.5- Read the forecast state vector (xf)'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOROA : ',
     $              'reading the forecast state vector (xf)'
            ENDIF
C
            IF (valbase.GE.1) THEN
              CALL readxyo (kinrefxyo,vectxf(:),jnxyo,
     $             lectinfo,kflaganlxyo,poscoefobs(:,:))
            ENDIF
C
            vectxa(:) = vectxf(:) + vectxa(:)
C
         ENDIF
C
C -9.6- Compute the analysed error covariance matrix (Sa)
C -------------------------------------------------------
C
         IF (lPaout) THEN
C
            IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
C            print *,'-9.6- Compute the analysed ',
C     $           'error covariance matrix (Sa)'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOROA : ',
     $             'computing the analysed error covariance matrix (Sa)'
            ENDIF
            IF (kflaggloloc.EQ.0) THEN
               CALL prodmat_zr_rr_vias (basexr(:,jrbasdeb:jrbasfin),
     $              mat_yr(:,jrbasdeb:jrbasfin),
     $              matCUrrz(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin,1),
     $              jrmatdeb,jrmatfin)
            ELSE
               CALL prodmat_wr_rrz_vias(basexr(:,jrbasdeb:jrbasfin),
     $              mat_yr(:,jrbasdeb:jrbasfin),
     $              matCUrrz(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin,0:),
     $              vectxpart(:),
     $              jrmatdeb,jrmatfin)
            ENDIF
C
         ENDIF
C
C -9.7- Write the analysed state vector (xa) or (xa - xf)
C -------------------------------------------------------
C
         IF (lXaout) THEN
C
            IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
C            print *,'-9.7- Write the analysed state vector (xa)',
C     $           ' or (xa - xf)'
               WRITE(numout,*)
               IF (lXa) WRITE(numout,*) 'ALGOROA : ',
     $              'writing the analysed state (xa)'
               IF (.NOT.(lXa)) WRITE(numout,*) 'ALGOROA : ',
     $              'writing the (Xa - Xf)'
            ENDIF
            SELECT CASE (kflaganlxyo)
            CASE (1)
              CALL writevar(koutxyo,vectxa(:),jnxyo)
            CASE (2)
              CALL writedta(koutxyo,vectxa(:))
            CASE (3)
              CALL writeobs(koutxyo,vectxa(:),vectorms(:),
     $             gridijkobs(:),poscoefobs(:,:))
            CASE DEFAULT
               GOTO 1000
            END SELECT
C
         ENDIF
C
C -9.8- Write the error covariance matrix (Sa)
C --------------------------------------------
C
         IF (lPaout) THEN
C --- affectation of outbas
            IF (jnxyo.EQ.1) CALL writeinfobas(koutbasxyo,kflaganlxyo)
C
C -9.8.1- Compute and write the diagonal of Pa into output basis directory
C ------------------------------------------------------------------------
            IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
C            print *,'-9.8.3- Compute and write  diagonal'
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOROA : ',
     $              'Compute and write the diagonal of Pa',
     $              ' into output basis directory'
            ENDIF
            vectxweight(:)=FREAL(0.0)
            DO jr=jrbasdeb,jrbasfin
               vectxweight(:)=basexr(:,jr)*basexr(:,jr)+vectxweight(:)
            ENDDO
            vectxweight(:)=SQRT(vectxweight(:))
            numjr=0
            serie=1
            CALL fildirbas (vctnamout,koutbasxyo,jprbasout,numjr,serie)
            WRITE(fnameout,'("./",A,"/",A)') 
     $           koutbasxyo(1:lenv(koutbasxyo)),
     $           vctnamout(1:lenv(vctnamout))
            SELECT CASE (kflaganlxyo)
            CASE (1)
               CALL writevar (fnameout,vectxweight(:),jnxyo)
            CASE (2)
               CALL writedta (fnameout,vectxweight(:))
            CASE (3)
               CALL writeobs(fnameout,vectxweight(:),vectorms(:),
     $              gridijkobs(:),poscoefobs(:,:))
            CASE DEFAULT
               GOTO 1000
            END SELECT
C
C -9.8.2- Compute updated ensemble if required
C --------------------------------------------
            IF (valbase.LT.1) THEN
              DO jx=1,jpxsize
                basexr(jx,:) = basexr(jx,:) * SQRT(FREAL(jprsize-1))
                basexr(jx,:) = basexr(jx,:) + vectxa(jx)
              ENDDO
            ENDIF
C
C -9.8.3- Writing output basis directory 
C --------------------------------------
            IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
C            print *,'-9.8.2- output base directory',
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOROA : writing the covariance (Sa)'
            ENDIF
            SELECT CASE (kflaganlxyo)
            CASE (1)
              CALL writebas(koutbasxyo,basexr(:,:),
     $             jnxyo,jrbasdeb,jrbasfin)
            CASE (2)
              CALL writeyobas(koutbasxyo,basexr(:,:),
     $             jrbasdeb,jrbasfin,kflaganlxyo)
            CASE (3)
              CALL writeyobas(koutbasxyo,basexr(:,:),
     $             jrbasdeb,jrbasfin,kflaganlxyo,vectorms(:),
     $             gridijkobs(:),poscoefobs(:,:))
            CASE DEFAULT
              GOTO 1000
            END SELECT
C
         ENDIF
C
      ENDDO
C
C -10.- Compute the retrospective analysis error covariance matrix (Sa)
C       and the analysis state (Xa)
C ------------------------------------------------------
C
      IF (larginsmocfg) THEN

        print *, '--> Perfom retrospective analysis... '
C
C -10.0- Read analysis filenames (states) and directory names
C        (covariance matrix)
C -----------------------------------------------------------
C
C --- read number of retrospective analysis steps: jsmo
        CALL evalhdrcfgsmo(arginsmocfg,argoutsmocfg,jsmo)
              print*, arginsmocfg(1:lenv(arginsmocfg)),' ',
     s                argoutsmocfg(1:lenv(argoutsmocfg)),' ',jsmo

        IF (jsmo.LT.1) THEN
          print*, 'jsmo<1: NO RETROSPECTIVE ANALYSIS'
        ELSE

C --- allocate vctnaminsmo,vctnamoutsmo,dirnaminsmo,dirnamoutsmo
          allocate ( vctnaminsmo(1:jsmo), stat=allocok )
          vctnaminsmo(:)=''
          allocate ( vctnamoutsmo(1:jsmo), stat=allocok )
          vctnamoutsmo(:)=''
          allocate ( dirnaminsmo(1:jsmo), stat=allocok )
          dirnaminsmo(:)=''
          allocate ( dirnamoutsmo(1:jsmo), stat=allocok )
          dirnamoutsmo(:)=''
          IF (allocok.NE.0) GOTO 1001

C Read configuration file names
          CALL readcfgsmo(arginsmocfg,argoutsmocfg,jsmo,
     $                vctnaminsmo,vctnamoutsmo,dirnaminsmo,dirnamoutsmo)

C Loop on the retrospective analysis indices
C
          DO jpsmo=1,jsmo
             print*, vctnaminsmo(jpsmo)(1:lenv(vctnaminsmo(jpsmo))),' ',
     s             vctnamoutsmo(jpsmo)(1:lenv(vctnamoutsmo(jpsmo))),' ',
     s             dirnaminsmo(jpsmo)(1:lenv(dirnaminsmo(jpsmo))),' ',
     s                dirnamoutsmo(jpsmo)(1:lenv(dirnamoutsmo(jpsmo)))

C Loop on the partition
            lectinfo = .FALSE.
            DO jnxyo=1,limjpnxyo(kflaganlxyo)
              lmodprint=(MOD(jnxyo-1,(limjpnxyo(kflaganlxyo)/5+1)).EQ.0)
              IF (lmodprint) print *,'Memory part number : ',
     $                 jnxyo,'/',limjpnxyo(kflaganlxyo)
C
C -10.1- Read the partition in subdomains
C --------------------------------------
              IF (kflaggloloc.NE.0) THEN
C
                IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
                   WRITE(numout,*)
                   WRITE(numout,*) 'ALGOROA, smoother: ',
     $              'reading the partition in subdomains'
                ENDIF
                IF (kflaganlxyo.NE.3) THEN
                  lmoyectold=lmoyect
                  lmoyect=.FALSE.
                  CALL readxyo (kinpartxyo,vectxpart(:),jnxyo,
     $                          lectinfo,kflaganlxyo,poscoefobs(:,:))
                  lmoyect=lmoyectold
                ELSE
C --- allocation vecty
                  allocate ( vecty(1:jpysize), stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  vecty(:) = FREAL(0.0)
                  flagxyo=2
                  lmoyectold=lmoyect
                  lmoyect=.FALSE.
                  CALL readxyo (kinpartxyo,vecty(:),jnxyo,
     $                         lectinfo,flagxyo,poscoefobs(:,:))
                  lmoyect=lmoyectold
                  vectxpart(:)=vecty(poscoefobs(:,1)%pos)
                  IF (allocated(vecty))  deallocate (vecty)
                ENDIF
C
              ENDIF
C
C -10.2- Read the input analysis covariance matrix Sa
C ---------------------------------------------------
C
              IF ((lXaout).OR.(lPaout)) THEN
C
                IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOROA, smoother : ',
     $                 'reading the forecast input basis Sa ',jpsmo
                ENDIF
                CALL readbas (dirnaminsmo(jpsmo),basexr(:,:),jnxyo,
     $           jrbasdeb,jrbasfin,lectinfo,kflaganlxyo,poscoefobs(:,:))
C
              ENDIF
C
C -10.3- Compute the correction (xa_after - xa_before)
C -----------------------------------------------------------------
C
              IF (lXaout) THEN
C
                IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOROA, smoother : ',
     $                'computing the update analysis state nb ',jpsmo
                ENDIF
                basexr(:,1) = FREAL(0.0)
                IF (kflaggloloc.EQ.0) THEN
                  DO jr=jrbasdeb,jrbasfin
                    basexr(:,1) = basexr(:,1)+
     $                      basexr(:,jr)*coefrz(jr,1)
                  ENDDO
                ELSE
                  DO jr=jrbasdeb,jrbasfin
                    DO jx=1,jpxsize
                      basexr(jx,1)=basexr(jx,1)+
     $                   basexr(jx,jr)*coefrz(jr,nint(vectxpart(jx)))
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
C
C -10.4- Read the weight from the reducevar config and apply on
C        correction
C -----------------------------------------------------------
C
              IF ((largreducevar).AND.(lXaout)) THEN
C
                IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOROA, smoother : ',
     $              'reading the weight from the reducevar config'
                ENDIF
                CALL readxyo (argreducevar,vectxweight(:),
     $               jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
                basexr(:,1) = vectxweight(:)*basexr(:,1)
C
              ENDIF
C
C -10.5- Read the input state vector (xa_before) and add correction
C -----------------------------------------------------------------
C
              IF ((lXaout).AND.(lXa)) THEN
C
                IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOROA, smoother : ',
     $                 'reading the input state vector ',jpsmo
                ENDIF
                CALL readxyo (vctnaminsmo(jpsmo),vectxweight(:),jnxyo,
     $              lectinfo,kflaganlxyo,poscoefobs(:,:))
                basexr(:,1) = vectxweight(:) + basexr(:,1)
C
              ENDIF
C
C -10.6- Compute the updated analysis error covariance matrix (Sa)
C -------------------------------------------------------
C
              IF (lPaout) THEN
C
                IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOROA, smoother : ',
     $      'computing the analysed error covariance matrix (Sa) ',jpsmo
                ENDIF
                IF (kflaggloloc.EQ.0) THEN
                  CALL prodmat_zr_rr_vias (basexr(:,jrbasdeb:jrbasfin),
     $                 mat_yr(:,jrbasdeb:jrbasfin),
     $                 matCUrrz(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin,1),
     $                 jrmatdeb,jrmatfin)
                ELSE
                  CALL prodmat_wr_rrz_vias(basexr(:,jrbasdeb:jrbasfin),
     $                 mat_yr(:,jrbasdeb:jrbasfin),
     $                 matCUrrz(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin,0:),
     $                 vectxpart(:),
     $                 jrmatdeb,jrmatfin)
                ENDIF
C
              ENDIF
C
C -10.7- Write the analysed state vector (xa) or (xa - xf)
C -------------------------------------------------------
C
              IF (lXaout) THEN
C
                SELECT CASE (kflaganlxyo)
                CASE (1)
                  CALL writevar(vctnamoutsmo(jpsmo),basexr(:,1),jnxyo)
                CASE (2)
                  CALL writedta(vctnamoutsmo(jpsmo),basexr(:,1))
                CASE (3)
                  CALL writeobs(vctnamoutsmo(jpsmo),basexr(:,1),
     $                 vectorms(:),gridijkobs(:),poscoefobs(:,:))
                CASE DEFAULT
                  GOTO 1000
                END SELECT
C
              ENDIF
C
C -10.8- Write the smoother error covariance matrix (Sa)
C --------------------------------------------
C
              IF (lPaout) THEN
C
C -10.8.2- Writing output basis directory
C --------------------------------------
                IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
C            print *,'-9.8.2- output base directory',
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOROA, smoother:
     $                            writing the covariance (Sa)'
                ENDIF
                SELECT CASE (kflaganlxyo)
                CASE (1)
                  CALL writebas(dirnamoutsmo(jpsmo),basexr(:,:),
     $                          jnxyo,jrbasdeb,jrbasfin)
                CASE (2)
                  CALL writeyobas(dirnamoutsmo(jpsmo),basexr(:,:),
     $                    jrbasdeb,jrbasfin,kflaganlxyo)
                CASE (3)
                  CALL writeyobas(dirnamoutsmo(jpsmo),basexr(:,:),
     $                    jrbasdeb,jrbasfin,kflaganlxyo,vectorms(:),
     $                    gridijkobs(:),poscoefobs(:,:))
                CASE DEFAULT
                  GOTO 1000
                END SELECT
C --- affectation of outbas
                IF (jnxyo.EQ.1) CALL writeinfobas(dirnamoutsmo(jpsmo),kflaganlxyo)
C
C -10.8.3- Compute and write the diagonal of Pa into output basis
C directory
C ------------------------------------------------------------------------
                IF ((jnxyo.EQ.1).AND.(nprint.GE.2)) THEN
C            print *,'-9.8.3- Compute and write  diagonal'
                  WRITE(numout,*)
                  WRITE(numout,*) 'ALGOROA, smoother: ',
     $                  'Compute and write the diagonal of Pa',
     $                  ' into output basis directory'
                ENDIF
                basexr(:,1)=FREAL(0.0)
                DO jr=jrbasdeb,jrbasfin
                  basexr(:,1)= basexr(:,jr)*basexr(:,jr)+basexr(:,1)
                ENDDO
                basexr(:,1)=SQRT(basexr(:,1))
                numjr=0
                serie=1
                CALL fildirbas (vctnamout,dirnamoutsmo(jpsmo),jprbasout,numjr,serie)
                WRITE(fnameout,'("./",A,"/",A)')
     $                dirnamoutsmo(jpsmo)(1:lenv(dirnamoutsmo(jpsmo))),
     $                vctnamout(1:lenv(vctnamout))
                SELECT CASE (kflaganlxyo)
                CASE (1)
                  CALL writevar (fnameout,basexr(:,1),jnxyo)
                CASE (2)
                  CALL writedta (fnameout,basexr(:,1))
                CASE (3)
                  CALL writeobs(fnameout,basexr(:,1),vectorms(:),
     $                 gridijkobs(:),poscoefobs(:,:))
                CASE DEFAULT
                  GOTO 1000
                END SELECT
C
C end if lPaout:
              ENDIF
C end loop partition:
            ENDDO
C end loop smoother steps:
          ENDDO
C end if jsmo>1:
        ENDIF
C end if larginsmocfg:
      ENDIF
C
C -11.- Closure
C -------------
C
C --- deallocation
      IF (allocated(vectxweight)) deallocate (vectxweight)
      IF (allocated(vectxpart)) deallocate (vectxpart)
      IF (allocated(vectorms)) deallocate (vectorms)
      IF (allocated(poscoefobs)) deallocate (poscoefobs)
      IF (allocated(gridijkobs)) deallocate (gridijkobs)
      IF (allocated(basexr)) deallocate (basexr)
      IF (allocated(vectxa)) deallocate (vectxa)
      IF (allocated(vectxf)) deallocate (vectxf)
      IF (allocated(coefrz)) deallocate (coefrz)
      IF (allocated(matCUrrz)) deallocate (matCUrrz)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algoroa','roa')
 1001 CALL printerror2(0,1001,3,'algoroa','roa')
 1003 CALL printerror2(0,1003,3,'algoroa','roa')
C
 101   WRITE (texterror,*) 'argument coefrmax not valid :',
     $     argcoefrmax(1:lenv(argcoefrmax))
      CALL printerror2(0,101,3,'algoroa','calcroa',comment=texterror)
 102  WRITE (texterror,*) 'Invalid data type for -inpartobs'
      CALL printerror2(0,102,3,'algoroa','calcroa',comment=texterror)
 103  WRITE (texterror,*) 'Parallel algorithm: -fixjpx required'
      CALL printerror2(0,103,3,'algoroa','calcroa',comment=texterror)
 104  WRITE (texterror,*) 'Invalid alpha scaling factor'
      CALL printerror2(0,104,3,'algoroa','calcroa',comment=texterror)
 105  WRITE (texterror,*) 'Invalid beta scaling factor'
      CALL printerror2(0,105,3,'algoroa','calcroa',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algoroa
