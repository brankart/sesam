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
C ---                   ALGOPERC.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2008-03 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE calcperc : Compute precentiles of an ensemble of vectors
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algoperc
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcperc

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcperc(kinxyobas,koutxyobasref,kflagxyo,kconfigo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute precentiles of an ensemble of vectors
CCC
CCC  Method : Read ensemble from input directory (kinxyobas)
CCC  ------   compute percentiles, and
CCC           write the result in output directory (koutxyobasref).
CCC
CCC  Input : kinxyobas     : Cxyo input directory
CCC  -----   koutxyobasref : Cxyo output directory
CCC          kflagxyo      : Vector type (1=Vx,2=Vy,3=Vo)
CCC          kconfigo    : Observation operator (Vo) inputy file
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only :
     $     jpoend,jpitpend,jpx,jpxend,jpyend,jprend,jpperc,
     $     poscoefobs,gridijkobs,arraynx_jindxbeg
      use hioxyo
      use hiobas
      use ensdam_anaqua
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kinxyobas,koutxyobasref,kconfigo
      INTEGER, intent(in) :: kflagxyo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:,:), allocatable, save :: basesr
      BIGREAL, dimension(:,:), allocatable, save :: percsr
      BIGREAL, dimension(:), allocatable :: vectr
      BIGREAL, dimension(:), allocatable :: percdef
      BIGREAL, dimension(:), allocatable :: vectw0
      BIGREAL, dimension(:), allocatable :: vectw1
      BIGREAL, dimension(:), allocatable :: vectw
      BIGREAL, dimension(:), allocatable :: vectorms
C
      INTEGER :: allocok,jpssize,jpitpsize,jprsize
      INTEGER :: jnxyo,js,jr,jperc,il,im,ih
      LOGICAL :: lectinfo,lmodprint,lweight
      INTEGER :: valbase,jrbasdeb,jrbasfin,flagcfg
      BIGREAL :: a,b,percscal
      CHARACTER(len=bgword) :: fname
CC----------------------------------------------------------------------
C
      jprsize=jprend
      jpitpsize=1
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modanam/algoperc :'
         WRITE(numout,*) '         compute ensemble percentiles'
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
C Allocate Cxyo array
      allocate ( basesr(1:jpssize,1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basesr(:,:) = FREAL(0.0)
C
C Allocate Cxyo percentile array
      allocate ( percsr(1:jpssize,1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percsr(:,:) = FREAL(0.0)
C
C Allocate vectr (ensemble size) vector
      allocate ( vectr(1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectr(:) = FREAL(0.0)
C
C Allocate vectw (ensemble size) vector
      allocate ( vectw(1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectw(:) = FREAL(0.0)
C
C Allocate vectw0 (ensemble size) vector
      allocate ( vectw0(1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectw0(:) = FREAL(0.0)
C
C Allocate vectw1 (ensemble size) vector
      allocate ( vectw1(1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectw1(:) = FREAL(0.0)
C
C Allocate percdef (percentile definition) vector
      allocate ( percdef(1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percdef(:) = FREAL(0.0)
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
C Read percentiles definition
      IF (nprint.GE.2) THEN
          WRITE(numout,*) '    ==> READING percentile definition'
      ENDIF
C
      CALL readscalbas(koutxyobasref,'percdef',percdef)
C
      DO jnxyo=1+jproc,limjpnxyo(MIN(3,kflagxyo)),jpproc
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
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT

C       Read member weights if any
        fname=kinxyobas(1:lenv(kinxyobas)) // '/weight'
        INQUIRE(FILE=fname,EXIST=lweight)
        IF (lweight) CALL readscalbas(kinxyobas,'weight',vectw0)
C
C -2.- Compute percentiles
C ------------------------
C
        IF (traditional) THEN
C TRADITIONAL SCHEME
C ==================
        IF (lweight) THEN

          DO js=1,jpssize
C
            lmodprint=(MOD(js-1,(jpssize/5+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprint))) print *,
     $          'Variable index : ',js,'/',jpssize

C           Get ensemble members for current variable
            vectr(1:jprsize)=basesr(js,1:jprsize)
C           Get weight for each ensemble member
            vectw(:)=vectw0(:)

C           Sort ensemble members
            CALL heapsort2(vectr,vectw)

C           Compute cumulated weights
            DO jr=2,jprsize
              vectw1(jr)=(vectw(jr-1)+vectw(jr))/2
            ENDDO

            vectw(1)=0.0
            DO jr=2,jprsize
              vectw(jr)=vectw(jr-1)+vectw1(jr)
            ENDDO

C           Rescale between 0 and 1
            vectw(:)=vectw(:)/vectw(jprsize)

C           For each required quantile
            DO jperc=1,jpperc
C             Find rank by dichotomy in cumulated weights
              il = 1
              ih = jprsize
              DO WHILE (il.LT.ih-1)
                im = (il+ih)/2
                IF (percdef(jperc).GT.vectw(im)) THEN
                  il = im
                ELSE
                  ih = im
                ENDIF
              ENDDO
              jr = il
C             Compute interpolation weights (a,b)
              IF ( vectw(il+1) == vectw(il) ) THEN
                b = 0.
              ELSE
                b = ( percdef(jperc) - vectw(il) ) / ( vectw(il+1) - vectw(il) )
              ENDIF
              a  = 1.0 - b
C             Interpolate to get the required quantile
              percsr(js,jperc)=a*vectr(jr)+b*vectr(jr+1)
            ENDDO

          ENDDO

        ELSE

          DO js=1,jpssize
C
            lmodprint=(MOD(js-1,(jpssize/5+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprint))) print *,
     $          'Variable index : ',js,'/',jpssize
C
            vectr(1:jprsize)=basesr(js,1:jprsize)
            CALL heapsort(vectr)
            DO jperc=1,jpperc
              percscal = 1+percdef(jperc)*(jprsize-1)
              jr = INT(percscal)
              b  = percscal - jr
              a  = 1.0 - b
              percsr(js,jperc)=a*vectr(jr)
              IF (b.GT.0._kr) THEN
                percsr(js,jperc)=percsr(js,jperc)+b*vectr(jr+1)
              ENDIF
            ENDDO
C
          ENDDO

        ENDIF
C
        ELSE
C NEW SCHEME
C ==========
        IF (lweight) THEN
          CALL ens_quantiles(percsr,basesr,percdef,enswei=vectw0)
        ELSE
          CALL ens_quantiles(percsr,basesr,percdef)
        ENDIF
C ==========
        ENDIF
C     
C -3.- Write percentiles
C --------------------------
C
        jrbasdeb=1
        jrbasfin=jpperc
        SELECT CASE (kflagxyo)
        CASE (1)
          CALL writebas(koutxyobasref,percsr(:,:),jnxyo,jrbasdeb,jrbasfin)
        CASE (2)
          CALL writeyobas(koutxyobasref,percsr(:,:),jrbasdeb,jrbasfin,
     $                    kflagxyo)
        CASE (3)
          CALL writeyobas(koutxyobasref,percsr(:,:),jrbasdeb,jrbasfin,
     $                    kflagxyo,vectorms(:),gridijkobs(:),
     $                    poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
C
      ENDDO
C
C --- deallocation
      IF (allocated(basesr)) deallocate(basesr)
      IF (allocated(percsr)) deallocate(percsr)
      IF (allocated(vectr)) deallocate(vectr)
      IF (allocated(percdef)) deallocate(percdef)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algoanam','algoanambas')
 1001 CALL printerror2(0,1001,3,'algoanam','algoanambas')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algoperc
