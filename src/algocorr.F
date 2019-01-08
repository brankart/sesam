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
C ---                   ALGOCORR.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2004-06 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE calccorr
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algocorr
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calccorr

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calccorr (kinxyobas,koutxyo,kjxyo,kjtype,kflagxyo,kconfigo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute correlation coefficients from reduced
CCC  -------   order covariance matrix
CCC
CCC  Method : Read reduced order covariance matrix from input
CCC  ------   directory (inxbas), compute correlations between
CCC           all state variables and the one specified by the
CCC           user (index kjxyo), write the result in output
CCC           file (koutxyo).
CCC
CCC  Input : kinxyobas   : Cxyo input directory
CCC  -----   koutxyo     : Vxyo output file
CCC          kjxyo       : Index of the state variable with which to
CCC                        compute correlations
CCC          kjtype      : Type of output vector (1=correlation,
CCC                        2=representer)
CCC          kflagxyo    : Vector type (1=Vx,2=Vy,3=Vo)
CCC          kconfigo    : Observation operator (Vo) inputy file
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only :
     $     jpoend,jpitpend,jpx,jpxend,jpyend,jprend,
     $     poscoefobs,gridijkobs,
     $     spvalvar,spvaldta,spvalobs,arraynx_jindxbeg
      use hioxyo
      use hiobas
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kinxyobas,koutxyo,kconfigo
      INTEGER, intent(in) :: kflagxyo, kjxyo, kjtype
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:,:), allocatable, save :: basesr
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:), allocatable, save :: vectstd
      BIGREAL, dimension(:), allocatable, save :: vectorms
C
      BIGREAL, dimension(:), allocatable :: vectr
      INTEGER :: allocok,jpssize,jpysize,jpitpsize,jprsize
      INTEGER :: jnxyo,flagcfg,js,jr
      BIGREAL :: refstd,spval
      LOGICAL :: lectinfo,existence
      INTEGER :: valbase, jrbasdeb,jrbasfin, jr0, jr1
CC----------------------------------------------------------------------
C
      jprsize=jprend
      jpysize=jpyend
      jpitpsize=1
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modcorr/algocorr :'
         WRITE(numout,*) '         compute correlations (Vxyo) from'
         WRITE(numout,*) '         reduced order cov. matrix (Cxyo)'
      ENDIF
C
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
C
C Allocate Vxyo arrays
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
      allocate ( vectstd(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectstd(:) = FREAL(0.0)
C
C Allocate Cxyo array
      allocate ( basesr(1:jpssize,1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basesr(:,:) = FREAL(0.0)
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
C Read header information from reduced order cov. matrix directory
      CALL readinfobas(kinxyobas,valbase)
      IF (valbase.LT.1) THEN
        PRINT *, 'Warning: input directory is an ensemble'
        jr0 = 1 ; jr1 = jprsize
      ELSE
        PRINT *, 'Warning: input directory is a covariance square root'
        jr0 = 2 ; jr1 = jprsize
      ENDIF
C
C -1.- Read reduced order covariance matrix for required variable
C ---------------------------------------------------------------
C Store it in vectr vector
C (if option -fixjpx is activated)
C
      IF (limjpnxyo(MIN(3,kflagxyo)).GT.1) THEN
C
        allocate ( vectr(1:jprsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectr(:) = FREAL(0.0)
C
        jnxyo=1
        DO WHILE (kjxyo.GE.arraynx_jindxbeg(jnxyo+1))
          jnxyo=jnxyo+1 
        ENDDO
C
        jrbasdeb=1
        jrbasfin=jprsize
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $        lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $        lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
C
        vectr(:)=basesr(kjxyo-arraynx_jindxbeg(jnxyo)+1,:)
C
        IF (valbase.LT.1) THEN
          vectr(:) = vectr(:) - SUM(vectr(jr0:jr1))/jprsize
          vectr(:) = vectr(:) / SQRT(FREAL(jprsize-1))
        ENDIF
C
        refstd = FREAL(0.0)
        DO jr=jr0,jr1
          refstd = refstd + vectr(jr)*vectr(jr)
        ENDDO
        refstd = SQRT(refstd)
        IF ( refstd.LE.FREAL(0.0) ) GOTO 102
C
      ENDIF
C
      DO jnxyo=1,limjpnxyo(MIN(3,kflagxyo))
C
C -2.- Read input reduced order covariance matrix
C -----------------------------------------------
C
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING the covariance matrix'
        ENDIF
C
        jrbasdeb=1
        jrbasfin=jprsize
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
C
        IF (valbase.LT.1) THEN
          DO js=1,jpssize
            basesr(js,:) = basesr(js,:) - SUM(basesr(js,jr0:jr1))/jprsize
            basesr(js,:) = basesr(js,:) / SQRT(FREAL(jprsize-1))
          ENDDO
        ENDIF
C
C -3.- Compute correlation coefficients or representer
C ----------------------------------------------------
C
C Compute variances and covariances, and then standard deviations
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
C Check if standard deviation is positive
C and compute correlation coefficients
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
C     
C -4.- Write correlation coefficients or representer
C --------------------------------------------------
C
        SELECT CASE (kflagxyo)
        CASE (1)
           CALL writevar (koutxyo,vects(:),jnxyo)
        CASE (2)
           CALL writedta (koutxyo,vects(:))
        CASE (3)
           CALL writeobs (koutxyo,vects(:),vectorms(:),
     $              gridijkobs(:),poscoefobs(:,:))
        CASE DEFAULT
           GOTO 1000
        END SELECT
C
      ENDDO
C
C --- deallocation
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(vectstd)) deallocate(vectstd)
      IF (allocated(basesr)) deallocate(basesr)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algocorr','calccorr')
 1001 CALL printerror2(0,1001,3,'algocorr','calccorr')
C
 102  WRITE (texterror,*) 'Variable with zero variance'
      CALL printerror2(0,102,3,'algocorr','calccorr',
     $     comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algocorr
