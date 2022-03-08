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
! ---                   ALGOSCOR.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2015-03 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE calcscor
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algoscor
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcscor

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcscor (kinxbas,kinxyo,kflagxyo,kscore, &
     &                     kinpartvar,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose : Compute probabilistic score
!
!  Method : Accumulate data region by region
!  ------   Compute requested score (CRPS, RCRV)
!
!  Input : kinxbas     : Cxyo directory with input ensemble
!  -----   kinxyo      : verification vector
!          kflagxyo    : Vector type (1=Vx,2=Vy,3=Vo)
!          kscore      : Name of probabilistic score to compute
!          kconfigo    : Observation operator (Vo) inputy file
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : &
     &     jpoend,jpitpend,jpx,jpxend,jpyend,jprend,jpperc, &
     &     poscoefobs,gridijkobs,arraynx_jpindxend
      use hioxyo
      use hiobas
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinxyo,kinxbas,kinpartvar
      INTEGER, intent(in) :: kflagxyo
      CHARACTER(len=*), intent(in) :: kscore
      CHARACTER(len=*), intent(in) :: kconfigo
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:,:), allocatable, save :: ensemble
      BIGREAL, dimension(:), allocatable, save :: verif
      BIGREAL, dimension(:), allocatable, save :: partition
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      INTEGER :: allocok,jpssize,jpitpsize,jprsize
      INTEGER :: jnxyo,js,jr,jperc,jsend,jpregion,jregion
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbase,jrbasdeb,jrbasfin,flagcfg,eql
      BIGREAL :: eps
!
      BIGREAL, dimension(:,:), allocatable :: aa,bb
      BIGREAL, dimension(:), allocatable :: reli,resol,crps
      BIGREAL, dimension(:), allocatable :: bias,sprea
      INTEGER, dimension(:), allocatable :: nbr
!----------------------------------------------------------------------
!
      jprsize=jprend
      jpitpsize=1
! Maximum number of region
      jpregion=1000
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modscor/algoscor :'
         WRITE(numout,*) '         compute probabilistic score'
      ENDIF
!
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
!
! Allocate Cxyo percentile array
      allocate ( ensemble(1:jpssize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ensemble(:,:) = FREAL(0.0)
!
! Allocate Vxyo arrays
      allocate ( verif(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      verif(:) = FREAL(0.0)
!
      allocate ( partition(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      partition(:) = FREAL(0.0)
!
! Allocate score arrays
      SELECT CASE (kscore)
      CASE ('crps')
!
        allocate ( aa(0:jprsize,1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        aa(:,:) = FREAL(0.0)
!
        allocate ( bb(0:jprsize,1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        bb(:,:) = FREAL(0.0)
!
        allocate ( reli(1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        reli(:) = FREAL(0.0)
!
        allocate ( resol(1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        resol(:) = FREAL(0.0)
!
        allocate ( crps(1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        crps(:) = FREAL(0.0)
!
      CASE ('rcrv')
!
        allocate ( bias(1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        bias(:) = FREAL(0.0)
!
        allocate ( sprea(1:jpregion), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        sprea(:) = FREAL(0.0)
!
      CASE DEFAULT
        GOTO 101
      END SELECT
!
      allocate ( nbr(1:jpregion), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      nbr(:) = 0
!
      IF (kflagxyo.EQ.3) THEN
!
! Allocate poscoefobs array
        allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
!
! Allocate gridijkobs array
        allocate ( gridijkobs(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
!
! Allocate vectorms array
        allocate ( vectorms(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectorms(:) = FREAL(0.0)
!
! Read poscoefobs, vectorms and gridijkobs arrays
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
      DO jnxyo=1,limjpnxyo(MIN(3,kflagxyo))
        IF (kflagxyo.EQ.1) THEN
          jsend=arraynx_jpindxend(jnxyo)
        ELSE
          jsend=jpssize
        ENDIF
!
! -1.- Read input ensemble
! ------------------------
!
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING input ensemble'
        ENDIF
!
        jrbasdeb=1
        jrbasfin=jprsize
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxbas,ensemble(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxbas,ensemble(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
!     
! -2.- Read input vector
! ----------------------
!
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readxyo (kinxyo,verif(:),jnxyo,lectinfo,kflagxyo)
        CASE (3)
          CALL readxyo (kinxyo,verif(:),jnxyo,lectinfo,kflagxyo, &
     &                  poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
!     
! -3.- Read partition
! -------------------
!
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readxyo (kinpartvar,partition(:),jnxyo,lectinfo,kflagxyo)
        CASE (3)
          CALL readxyo (kinpartvar,partition(:),jnxyo,lectinfo,kflagxyo, &
     &                  poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
!
        IF (NINT(MINVAL(partition(:))).LT.0) GOTO 102
        IF (NINT(MAXVAL(partition(:))).GT.jpregion) GOTO 102
!
! -4.- Accumulate information from every region
! ---------------------------------------------
!
        eps = 1e-10
!
        SELECT CASE (kscore)
        CASE ('crps')
          DO js=1,jsend
            jregion = NINT(partition(js))
            !IF (ABS(FREAL(jregion)-partition(js)).GT.eps) jregion=0
            IF (jregion.GT.0) THEN
              nbr(jregion) = nbr(jregion) + 1
              CALL crps_cumul(ensemble(js,:),verif(js), &
     &              aa(0:,jregion),bb(0:,jregion))
            ENDIF
          ENDDO
        CASE ('rcrv')
          DO js=1,jsend
            jregion = NINT(partition(js))
            !IF (ABS(FREAL(jregion)-partition(js)).GT.eps) jregion=0
            IF (jregion.GT.0) THEN
              nbr(jregion) = nbr(jregion) + 1
              CALL rcrv_cumul(ensemble(js,:),verif(js), &
     &              nbr(jregion),bias(jregion),sprea(jregion))
            ENDIF
          ENDDO
        CASE DEFAULT
          GOTO 101
        END SELECT
!
      ENDDO
!
! -5.- Compute probabilistic score
! --------------------------------
!
      SELECT CASE (kscore)
      CASE ('crps')
        DO jregion=1,jpregion
          IF (nbr(jregion).GT.0) THEN
            aa(0:,jregion) = aa(0:,jregion) / nbr(jregion)
            bb(0:,jregion) = bb(0:,jregion) / nbr(jregion)
            CALL crps_score(aa(0:,jregion),bb(0:,jregion), &
     &                      reli(jregion),resol(jregion),crps(jregion))
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
!     
! -6.- Write probabilistic score
! ------------------------------
!
      SELECT CASE (kscore)
      CASE ('crps')
        PRINT *, 'CRPS score'
        PRINT *, '=========='
        PRINT *, 'Region  Size        Reliability    Resolution   CRPS'
        DO jregion=1,jpregion
          IF (nbr(jregion).GT.0) THEN
            PRINT '(i6.4,i10,3e15.8)',jregion,nbr(jregion), &
     &           reli(jregion),resol(jregion),crps(jregion)
          ENDIF
        ENDDO
      CASE ('rcrv')
        PRINT *, 'RCRV score'
        PRINT *, '=========='
        PRINT *, 'Region  Size        Bias           Spread'
        DO jregion=1,jpregion
          IF (nbr(jregion).GT.0) THEN
            PRINT '(i6.4,i10,2e15.8)',jregion,nbr(jregion), &
     &            bias(jregion),sprea(jregion)
          ENDIF
        ENDDO
      CASE DEFAULT
        GOTO 101
      END SELECT
!
! --- deallocation
      IF (allocated(verif)) deallocate(verif)
      IF (allocated(ensemble)) deallocate(ensemble)
      IF (allocated(partition)) deallocate(partition)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
!
      IF (allocated(nbr)) deallocate(nbr)
      IF (allocated(aa)) deallocate(aa)
      IF (allocated(bb)) deallocate(bb)
      IF (allocated(reli)) deallocate(reli)
      IF (allocated(resol)) deallocate(resol)
      IF (allocated(crps)) deallocate(crps)
!
      IF (allocated(bias)) deallocate(bias)
      IF (allocated(sprea)) deallocate(sprea)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algoscor','algoscor')
 1001 CALL printerror2(0,1001,3,'algoscor','algoscor')
!
!
 101  WRITE (texterror,*) 'Bad score name : ',kscore
      CALL printerror2(0,101,3,'algoscor','algoscor', &
     &     comment=texterror)
 102  WRITE (texterror,*) 'Bad partition'
      CALL printerror2(0,102,3,'algoscor','algoscor', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE crps_cumul(e,a,aa,bb)
!---------------------------------------------------------------------
!
!  Purpose : Accumulate one more piece of information
!            to compute CRPS score
!
!  Method : 
!  ------   
!
!  Input : e     : input ensemble
!  -----   a     : verification (scalar)
!          aa    : array to update
!          bb    : array to update
!
!---------------------------------------------------------------------
      use mod_main
      use ensdam_anaqua
      IMPLICIT NONE
!
      BIGREAL, dimension(:), intent(inout) :: e
      BIGREAL, intent(in) :: a
      BIGREAL, dimension(0:), intent(inout) :: aa, bb
!-----------------------------------------------------------------------
      INTEGER :: n, i
!-----------------------------------------------------------------------
      n = size(e)
!
! Sort input ensemble
      CALL heapsort(e)
!
! Verfication smaller than all ensemble members
      IF(a.LT.e(1)) bb(0)=bb(0)+1.0
      IF(a.LT.e(1)) aa(0)=aa(0)+(e(1)-a)
!
! Verification inside ensemble range
      DO i=1,n-1
        bb(i)=bb(i)+MAX(MIN(a,e(i+1))-e(i),0.)
        aa(i)=aa(i)+MAX(e(i+1)-MAX(a,e(i)),0.)
      ENDDO
!
! Verfication larger than all ensemble members
      IF(a.GT.e(n)) bb(n)=bb(n)+(a-e(n))
      IF(a.GT.e(n)) aa(n)=aa(n)+1.0
!
      RETURN
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE crps_score(aa,bb,reli,resol,crps)
!---------------------------------------------------------------------
!
!  Purpose : Compute CRPS score from previously
!            accumulated information
!
!  Method : 
!  ------   
!
!  Input : aa    : accumulated information
!  -----   bb    : accumulated information
!          reli  : reliability score
!          resol : resolution score
!
!---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
!
      BIGREAL, dimension(0:), intent(inout) :: aa, bb
      BIGREAL, intent(out) :: reli, resol, crps
!-----------------------------------------------------------------------
      INTEGER :: n, i
      BIGREAL :: gi, oi, p
!-----------------------------------------------------------------------
      n = size(aa) - 1
!
! Reinterpretation of bb(0) and aa(n)
! (Note that this does not contribute to CRPS)
      IF (bb(0).NE.0.) bb(0)=aa(0)*(1./bb(0)-1.)
      IF (aa(n).NE.0.) aa(n)=bb(n)*(1./aa(n)-1.)
!
! Compute components of CRPS
      DO i=0,n
        gi = bb(i) + aa(i)
        IF (gi.NE.0.) oi = aa(i)/gi
        p = FREAL(i)/FREAL(n)
!
        crps  = crps  + bb(i)*p*p+aa(i)*(1.-p)*(1.-p)
        reli  = reli  + gi*(oi-p)**2
        resol = resol + gi*oi*(1.-oi)
      ENDDO
!
      RETURN
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE rcrv_cumul(e,a,idx,mean,sqrs)
!---------------------------------------------------------------------
!
!  Purpose : Accumulate one more piece of information
!            to compute RCRV score
!
!  Method : 
!  ------   
!
!  Input : e     : input ensemble
!  -----   a     : verification (scalar)
!          idx   : index of accumulated piece of information
!          mean  : mean to update
!          sqrs  : square sum to update
!
!---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
!
      BIGREAL, dimension(:), intent(in) :: e
      BIGREAL, intent(in) :: a
      INTEGER, intent(in) :: idx
      BIGREAL, intent(inout) :: mean, sqrs
!-----------------------------------------------------------------------
      INTEGER :: n, i
      BIGREAL :: xmk, xsk, meanp, y
!-----------------------------------------------------------------------
      n = size(e)
!
! Compute ensemble mean (xmk) and ensemble std (xsk)
      xmk=SUM(e(:))/FREAL(n)
      xsk=SQRT( SUM((e(:)-xmk)**2)/FREAL(n-1) )
!
! Compute reduced anomaly
      y = a - xmk
      IF (xsk.NE.0.) y = y / xsk
!
! Save previous mean
      meanp = mean
!
! Update mean
      mean = mean + ( y - mean ) / FREAL(idx)
!
! Update square sum
      sqrs = sqrs + ( y - mean ) * ( y - meanp )
!
      RETURN
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algoscor
