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
C ---                   ALGOSCOR.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2015-03 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE calcscor
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algoscor
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcscor

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcscor (kinxbas,kinxyo,kflagxyo,kscore,
     $                     kinpartvar,kconfigo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute probabilistic score
CCC
CCC  Method : Accumulate data region by region
CCC  ------   Compute requested score (CRPS, RCRV)
CCC
CCC  Input : kinxbas     : Cxyo directory with input ensemble
CCC  -----   kinxyo      : verification vector
CCC          kflagxyo    : Vector type (1=Vx,2=Vy,3=Vo)
CCC          kscore      : Name of probabilistic score to compute
CCC          kconfigo    : Observation operator (Vo) inputy file
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only :
     $     jpoend,jpitpend,jpx,jpxend,jpyend,jprend,jpperc,
     $     poscoefobs,gridijkobs,arraynx_jpindxend
      use hioxyo
      use hiobas
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kinxyo,kinxbas,kinpartvar
      INTEGER, intent(in) :: kflagxyo
      CHARACTER(len=*), intent(in) :: kscore
      CHARACTER(len=*), intent(in) :: kconfigo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:,:), allocatable, save :: ensemble
      BIGREAL, dimension(:), allocatable, save :: verif
      BIGREAL, dimension(:), allocatable, save :: partition
      BIGREAL, dimension(:), allocatable, save :: vectorms
C
      INTEGER :: allocok,jpssize,jpitpsize,jprsize
      INTEGER :: jnxyo,js,jr,jperc,jsend,jpregion,jregion
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbase,jrbasdeb,jrbasfin,flagcfg,eql
      BIGREAL :: eps
CC
      BIGREAL, dimension(:,:), allocatable :: aa,bb
      BIGREAL, dimension(:), allocatable :: reli,resol,crps
      BIGREAL, dimension(:), allocatable :: bias,sprea
      INTEGER, dimension(:), allocatable :: nbr
CC----------------------------------------------------------------------
C
      jprsize=jprend
      jpitpsize=1
C Maximum number of region
      jpregion=1000
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modscor/algoscor :'
         WRITE(numout,*) '         compute probabilistic score'
      ENDIF
C
      SELECT CASE (kflagxyo)
      CASE (1)
         jpssize=jpx
      CASE (2)
         jpssize=jpyend
      CASE (3)
         jpssize=jpoend
         jpitpsize=jpitpend
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Allocate Cxyo percentile array
      allocate ( ensemble(1:jpssize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ensemble(:,:) = FREAL(0.0)
C
C Allocate Vxyo arrays
      allocate ( verif(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      verif(:) = FREAL(0.0)
C
      allocate ( partition(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      partition(:) = FREAL(0.0)
C
C Allocate score arrays
      SELECT CASE (kscore)
      CASE ('crps')
C
        allocate ( aa(0:jprsize,1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        aa(:,:) = FREAL(0.0)
C
        allocate ( bb(0:jprsize,1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        bb(:,:) = FREAL(0.0)
C
        allocate ( reli(1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        reli(:) = FREAL(0.0)
C
        allocate ( resol(1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        resol(:) = FREAL(0.0)
C
        allocate ( crps(1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        crps(:) = FREAL(0.0)
C
      CASE ('rcrv')
C
        allocate ( bias(1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        bias(:) = FREAL(0.0)
C
        allocate ( sprea(1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        sprea(:) = FREAL(0.0)
C
      CASE DEFAULT
        GOTO 101
      END SELECT
C
      allocate ( nbr(1:jpregion), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      nbr(:) = 0
C
      IF (kflagxyo.EQ.3) THEN
C
C Allocate poscoefobs array
        allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
C
C Allocate gridijkobs array
        allocate ( gridijkobs(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
C
C Allocate vectorms array
        allocate ( vectorms(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectorms(:) = FREAL(0.0)
C
C Read poscoefobs, vectorms and gridijkobs arrays
        flagcfg=1
        CALL readcfgobs (kconfigo,flagcfg,
     $        kvectorms=vectorms(:))
        flagcfg=2
        CALL readcfgobs (kconfigo,flagcfg,
     $        kgridijkobs=gridijkobs(:))
        flagcfg=3
        CALL readcfgobs (kconfigo,flagcfg,
     $        kposcoefobs=poscoefobs(:,:))
      ENDIF
C
      DO jnxyo=1,limjpnxyo(MIN(3,kflagxyo))
        IF (kflagxyo.EQ.1) THEN
          jsend=arraynx_jpindxend(jnxyo)
        ELSE
          jsend=jpssize
        ENDIF
C
C -1.- Read input ensemble
C ------------------------
C
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING input ensemble'
        ENDIF
C
        jrbasdeb=1
        jrbasfin=jprsize
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxbas,ensemble(:,:),jnxyo,jrbasdeb,jrbasfin,
     $           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxbas,ensemble(:,:),jnxyo,jrbasdeb,jrbasfin,
     $           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
C     
C -2.- Read input vector
C ----------------------
C
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readxyo (kinxyo,verif(:),jnxyo,lectinfo,kflagxyo)
        CASE (3)
          CALL readxyo (kinxyo,verif(:),jnxyo,lectinfo,kflagxyo,
     $                  poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
C     
C -3.- Read partition
C -------------------
C
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readxyo (kinpartvar,partition(:),jnxyo,lectinfo,kflagxyo)
        CASE (3)
          CALL readxyo (kinpartvar,partition(:),jnxyo,lectinfo,kflagxyo,
     $                  poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
C
        IF (NINT(MINVAL(partition(:))).LT.0) GOTO 102
        IF (NINT(MAXVAL(partition(:))).GT.jpregion) GOTO 102
C
C -4.- Accumulate information from every region
C ---------------------------------------------
C
        eps = 1e-10
C
        SELECT CASE (kscore)
        CASE ('crps')
          DO js=1,jsend
            jregion = NINT(partition(js))
            IF (ABS(FREAL(jregion)-partition(js)).GT.eps) jregion=0
            IF (jregion.GT.0) THEN
              nbr(jregion) = nbr(jregion) + 1
              CALL crps_cumul(ensemble(js,:),verif(js),
     $              aa(0:,jregion),bb(0:,jregion))
            ENDIF
          ENDDO
        CASE ('rcrv')
          DO js=1,jsend
            jregion = NINT(partition(js))
            IF (ABS(FREAL(jregion)-partition(js)).GT.eps) jregion=0
            IF (jregion.GT.0) THEN
              nbr(jregion) = nbr(jregion) + 1
              CALL rcrv_cumul(ensemble(js,:),verif(js),
     $              nbr(jregion),bias(jregion),sprea(jregion))
            ENDIF
          ENDDO
        CASE DEFAULT
          GOTO 101
        END SELECT
C
      ENDDO
C
C -5.- Compute probabilistic score
C --------------------------------
C
      SELECT CASE (kscore)
      CASE ('crps')
        DO jregion=1,jpregion
          IF (nbr(jregion).GT.0) THEN
            aa(0:,jregion) = aa(0:,jregion) / nbr(jregion)
            bb(0:,jregion) = bb(0:,jregion) / nbr(jregion)
            CALL crps_score(aa(0:,jregion),bb(0:,jregion),
     $                      reli(jregion),resol(jregion),crps(jregion))
          ENDIF
        ENDDO
      CASE ('rcrv')
        DO jregion=1,jpregion
          IF (nbr(jregion).GT.1) THEN
            sprea(jregion) = sprea(jregion) / FREAL(nbr(jregion)-1)
            sprea(jregion) = SQRT(sprea(jregion))
          ENDIF
        ENDDO
      CASE DEFAULT
        GOTO 101
      END SELECT
C     
C -6.- Write probabilistic score
C ------------------------------
C
      SELECT CASE (kscore)
      CASE ('crps')
        PRINT *, 'CRPS score'
        PRINT *, '=========='
        PRINT *, 'Region  Size        Reliability    Resolution   CRPS'
        DO jregion=1,jpregion
          IF (nbr(jregion).GT.0) THEN
            PRINT '(i6.4,i10,3e15.8)',jregion,nbr(jregion),
     $           reli(jregion),resol(jregion),crps(jregion)
          ENDIF
        ENDDO
      CASE ('rcrv')
        PRINT *, 'RCRV score'
        PRINT *, '=========='
        PRINT *, 'Region  Size        Bias           Spread'
        DO jregion=1,jpregion
          IF (nbr(jregion).GT.0) THEN
            PRINT '(i6.4,i10,2e15.8)',jregion,nbr(jregion),
     $            bias(jregion),sprea(jregion)
          ENDIF
        ENDDO
      CASE DEFAULT
        GOTO 101
      END SELECT
C
C --- deallocation
      IF (allocated(verif)) deallocate(verif)
      IF (allocated(ensemble)) deallocate(ensemble)
      IF (allocated(partition)) deallocate(partition)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
C
      IF (allocated(nbr)) deallocate(nbr)
      IF (allocated(aa)) deallocate(aa)
      IF (allocated(bb)) deallocate(bb)
      IF (allocated(reli)) deallocate(reli)
      IF (allocated(resol)) deallocate(resol)
      IF (allocated(crps)) deallocate(crps)
C
      IF (allocated(bias)) deallocate(bias)
      IF (allocated(sprea)) deallocate(sprea)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algoscor','algoscor')
 1001 CALL printerror2(0,1001,3,'algoscor','algoscor')
C
C
 101  WRITE (texterror,*) 'Bad score name : ',kscore
      CALL printerror2(0,101,3,'algoscor','algoscor',
     $     comment=texterror)
 102  WRITE (texterror,*) 'Bad partition'
      CALL printerror2(0,102,3,'algoscor','algoscor',
     $     comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE crps_cumul(e,a,aa,bb)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Accumulate one more piece of information
CCC            to compute CRPS score
CCC
CCC  Method : 
CCC  ------   
CCC
CCC  Input : e     : input ensemble
CCC  -----   a     : verification (scalar)
CCC          aa    : array to update
CCC          bb    : array to update
CCC
CCC---------------------------------------------------------------------
      use mod_main
      use ensdam_anaqua
      IMPLICIT NONE
C
      BIGREAL, dimension(:), intent(inout) :: e
      BIGREAL, intent(in) :: a
      BIGREAL, dimension(0:), intent(inout) :: aa, bb
C-----------------------------------------------------------------------
      INTEGER :: n, i
C-----------------------------------------------------------------------
      n = size(e)
C
C Sort input ensemble
      CALL heapsort(e)
C
C Verfication smaller than all ensemble members
      IF(a.LT.e(1)) bb(0)=bb(0)+1.0
      IF(a.LT.e(1)) aa(0)=aa(0)+(e(1)-a)
C
C Verification inside ensemble range
      DO i=1,n-1
        bb(i)=bb(i)+MAX(MIN(a,e(i+1))-e(i),0.)
        aa(i)=aa(i)+MAX(e(i+1)-MAX(a,e(i)),0.)
      ENDDO
C
C Verfication larger than all ensemble members
      IF(a.GT.e(n)) bb(n)=bb(n)+(a-e(n))
      IF(a.GT.e(n)) aa(n)=aa(n)+1.0
C
      RETURN
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE crps_score(aa,bb,reli,resol,crps)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute CRPS score from previously
CCC            accumulated information
CCC
CCC  Method : 
CCC  ------   
CCC
CCC  Input : aa    : accumulated information
CCC  -----   bb    : accumulated information
CCC          reli  : reliability score
CCC          resol : resolution score
CCC
CCC---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
C
      BIGREAL, dimension(0:), intent(inout) :: aa, bb
      BIGREAL, intent(out) :: reli, resol, crps
C-----------------------------------------------------------------------
      INTEGER :: n, i
      BIGREAL :: gi, oi, p
C-----------------------------------------------------------------------
      n = size(aa) - 1
C
C Reinterpretation of bb(0) and aa(n)
C (Note that this does not contribute to CRPS)
      IF (bb(0).NE.0.) bb(0)=aa(0)*(1./bb(0)-1.)
      IF (aa(n).NE.0.) aa(n)=bb(n)*(1./aa(n)-1.)
C
C Compute components of CRPS
      DO i=0,n
        gi = bb(i) + aa(i)
        IF (gi.NE.0.) oi = aa(i)/gi
        p = FREAL(i)/FREAL(n)
C
        crps  = crps  + bb(i)*p*p+aa(i)*(1.-p)*(1.-p)
        reli  = reli  + gi*(oi-p)**2
        resol = resol + gi*oi*(1.-oi)
      ENDDO
C
      RETURN
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE rcrv_cumul(e,a,idx,mean,sqrs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Accumulate one more piece of information
CCC            to compute RCRV score
CCC
CCC  Method : 
CCC  ------   
CCC
CCC  Input : e     : input ensemble
CCC  -----   a     : verification (scalar)
CCC          idx   : index of accumulated piece of information
CCC          mean  : mean to update
CCC          sqrs  : square sum to update
CCC
CCC---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
C
      BIGREAL, dimension(:), intent(in) :: e
      BIGREAL, intent(in) :: a
      INTEGER, intent(in) :: idx
      BIGREAL, intent(inout) :: mean, sqrs
C-----------------------------------------------------------------------
      INTEGER :: n, i
      BIGREAL :: xmk, xsk, meanp, y
C-----------------------------------------------------------------------
      n = size(e)
C
C Compute ensemble mean (xmk) and ensemble std (xsk)
      xmk=SUM(e(:))/FREAL(n)
      xsk=SQRT( SUM((e(:)-xmk)**2)/FREAL(n-1) )
C
C Compute reduced anomaly
      y = a - xmk
      IF (xsk.NE.0.) y = y / xsk
C
C Save previous mean
      meanp = mean
C
C Update mean
      mean = mean + ( y - mean ) / FREAL(idx)
C
C Update square sum
      sqrs = sqrs + ( y - mean ) * ( y - meanp )
C
      RETURN
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algoscor
