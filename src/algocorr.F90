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
! ---                   ALGOCORR.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2004-06 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE calccorr
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algocorr
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calccorr

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calccorr (kinxyobas,koutxyo,kjxyo,kjtype,kflagxyo,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose : Compute correlation coefficients from reduced
!  -------   order covariance matrix
!
!  Method : Read reduced order covariance matrix from input
!  ------   directory (inxbas), compute correlations between
!           all state variables and the one specified by the
!           user (index kjxyo), write the result in output
!           file (koutxyo).
!
!  Input : kinxyobas   : Cxyo input directory
!  -----   koutxyo     : Vxyo output file
!          kjxyo       : Index of the state variable with which to
!                        compute correlations
!          kjtype      : Type of output vector (1=correlation,
!                        2=representer)
!          kflagxyo    : Vector type (1=Vx,2=Vy,3=Vo)
!          kconfigo    : Observation operator (Vo) inputy file
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : &
     &     jpoend,jpitpend,jpx,jpxend,jpyend,jprend, &
     &     poscoefobs,gridijkobs, &
     &     spvalvar,spvaldta,spvalobs,arraynx_jindxbeg
      use hioxyo
      use hiobas
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinxyobas,koutxyo,kconfigo
      INTEGER, intent(in) :: kflagxyo, kjxyo, kjtype
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:,:), allocatable, save :: basesr
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:), allocatable, save :: vectstd
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      BIGREAL, dimension(:), allocatable :: vectr
      INTEGER :: allocok,jpssize,jpysize,jpitpsize,jprsize
      INTEGER :: jnxyo,flagcfg,js,jr
      BIGREAL :: refstd,spval
      LOGICAL :: lectinfo,existence
      INTEGER :: valbase, jrbasdeb,jrbasfin, jr0, jr1
!----------------------------------------------------------------------
!
      jprsize=jprend
      jpysize=jpyend
      jpitpsize=1
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modcorr/algocorr :'
         WRITE(numout,*) '         compute correlations (Vxyo) from'
         WRITE(numout,*) '         reduced order cov. matrix (Cxyo)'
      ENDIF
!
      SELECT CASE (kflagxyo)
      CASE (1)
         jpssize=jpx
         spval=spvalvar
      CASE (2)
         jpssize=jpyend
         spval=spvaldta
      CASE (3)
         jpssize=jpoend
         jpitpsize=jpitpend
         spval=spvalobs
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Allocate Vxyo arrays
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
      allocate ( vectstd(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectstd(:) = FREAL(0.0)
!
! Allocate Cxyo array
      allocate ( basesr(1:jpssize,1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basesr(:,:) = FREAL(0.0)
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
! Read header information from reduced order cov. matrix directory
      CALL readinfobas(kinxyobas,valbase)
      IF (valbase.LT.1) THEN
        PRINT *, 'Warning: input directory is an ensemble'
        jr0 = 1 ; jr1 = jprsize
      ELSE
        PRINT *, 'Warning: input directory is a covariance square root'
        jr0 = 2 ; jr1 = jprsize
      ENDIF
!
! -1.- Read reduced order covariance matrix for required variable
! ---------------------------------------------------------------
! Store it in vectr vector
! (if option -fixjpx is activated)
!
      IF (limjpnxyo(MIN(3,kflagxyo)).GT.1) THEN
!
        allocate ( vectr(1:jprsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectr(:) = FREAL(0.0)
!
        jnxyo=1
        DO WHILE (kjxyo.GE.arraynx_jindxbeg(jnxyo+1))
          jnxyo=jnxyo+1 
        ENDDO
!
        jrbasdeb=1
        jrbasfin=jprsize
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &        lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &        lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
!
        vectr(:)=basesr(kjxyo-arraynx_jindxbeg(jnxyo)+1,:)
!
        IF (valbase.LT.1) THEN
          vectr(:) = vectr(:) - SUM(vectr(jr0:jr1))/jprsize
          vectr(:) = vectr(:) / SQRT(FREAL(jprsize-1))
        ENDIF
!
        refstd = FREAL(0.0)
        DO jr=jr0,jr1
          refstd = refstd + vectr(jr)*vectr(jr)
        ENDDO
        refstd = SQRT(refstd)
        IF ( refstd.LE.FREAL(0.0) ) GOTO 102
!
      ENDIF
!
      DO jnxyo=1,limjpnxyo(MIN(3,kflagxyo))
!
! -2.- Read input reduced order covariance matrix
! -----------------------------------------------
!
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING the covariance matrix'
        ENDIF
!
        jrbasdeb=1
        jrbasfin=jprsize
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
!
        IF (valbase.LT.1) THEN
          DO js=1,jpssize
            basesr(js,:) = basesr(js,:) - SUM(basesr(js,jr0:jr1))/jprsize
            basesr(js,:) = basesr(js,:) / SQRT(FREAL(jprsize-1))
          ENDDO
        ENDIF
!
! -3.- Compute correlation coefficients or representer
! ----------------------------------------------------
!
! Compute variances and covariances, and then standard deviations
        IF (limjpnxyo(MIN(3,kflagxyo)).EQ.1) THEN
          refstd = FREAL(0.0)
          DO jr=jr0,jr1
            vects(:)   = vects(:)   + basesr(:,jr)*basesr(kjxyo,jr)
            vectstd(:) = vectstd(:) + basesr(:,jr)*basesr(:,jr)
            refstd = refstd + basesr(kjxyo,jr)*basesr(kjxyo,jr)
          ENDDO
          vectstd(:) = SQRT(vectstd(:))
          refstd = SQRT(refstd)
          IF ( refstd.LE.FREAL(0.0) ) GOTO 102
        ELSE
          vects(:) = FREAL(0.0)
          vectstd(:) = FREAL(0.0)
          DO jr=jr0,jr1
            vects(:)   = vects(:)   + basesr(:,jr)*vectr(jr)
            vectstd(:) = vectstd(:) + basesr(:,jr)*basesr(:,jr)
          ENDDO
          vectstd(:) = SQRT(vectstd(:))
        ENDIF
! Check if standard deviation is positive
! and compute correlation coefficients
        WHERE ( ( vectstd(1:jpssize).LE.FREAL(0.0) ) )
          vects(:) = spval
        ENDWHERE
        SELECT CASE (kjtype)
        CASE (1)
           WHERE ( ( vectstd(1:jpssize).GT.FREAL(0.0) ) )
             vects(:) = vects(:) / vectstd(:) / refstd
           ENDWHERE
        CASE (2)
           vects(:) = vects(:) / refstd / refstd
        CASE DEFAULT
          GOTO 1000
        END SELECT
!     
! -4.- Write correlation coefficients or representer
! --------------------------------------------------
!
        SELECT CASE (kflagxyo)
        CASE (1)
           CALL writevar (koutxyo,vects(:),jnxyo)
        CASE (2)
           CALL writedta (koutxyo,vects(:))
        CASE (3)
           CALL writeobs (koutxyo,vects(:),vectorms(:), &
     &              gridijkobs(:),poscoefobs(:,:))
        CASE DEFAULT
           GOTO 1000
        END SELECT
!
      ENDDO
!
! --- deallocation
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(vectstd)) deallocate(vectstd)
      IF (allocated(basesr)) deallocate(basesr)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algocorr','calccorr')
 1001 CALL printerror2(0,1001,3,'algocorr','calccorr')
!
 102  WRITE (texterror,*) 'Variable with zero variance'
      CALL printerror2(0,102,3,'algocorr','calccorr', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algocorr
