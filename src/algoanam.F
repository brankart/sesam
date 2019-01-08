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
C ---                   ALGOANAM.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2008-03 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE algoanam
C --- SUBROUTINE algoanambas
C --- SUBROUTINE algoanamobs
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algoanam
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC algoanamvct, algoanambas, algoanamobs

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoanamvct (kinxyo,koutxyo,kinxbasref,
     $                        kflagxyo,kflaganam,kconfigo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute anamorphosis of a vector of variables
CCC
CCC  Method : Read ensemble percentiles from input
CCC  ------   directory (kinxbasref), read vector of variables
CCC           from input file (kinxyo)
CCC           write the result in output file (koutxyo).
CCC
CCC  Input : kinxyo      : Vxyo input file
CCC  -----   koutxyo     : Vxyo output file
CCC          kinxbasref  : Cxyo percentiles directory
CCC          kflagxyo    : Vector type (1=Vx,2=Vy,3=Vo)
CCC          kflaganam   : anamorphosis direction (T=+ or F=-)
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
      use ensdam_anatra
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kinxyo,koutxyo,kinxbasref
      CHARACTER(len=*), intent(in) :: kconfigo
      INTEGER, intent(in) :: kflagxyo
      LOGICAL, intent(in) :: kflaganam
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:,:), allocatable, save :: percsr
      BIGREAL, dimension(:), allocatable :: percref
      BIGREAL, dimension(:), allocatable :: vectorms
C
      INTEGER :: allocok,jpssize,jpitpsize,jprsize
      INTEGER :: jnxyo,js,jr,jperc,jp,jjproc
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbase,jrbasdeb,jrbasfin,flagcfg,eql
      BIGREAL :: a,b,anaval
CC----------------------------------------------------------------------
C
      jprsize=jprend
      jpitpsize=1
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modanam/algoanam :'
         WRITE(numout,*) '         compute ensemble anamorphosis'
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
      allocate ( percsr(1:jpssize,1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percsr(:,:) = FREAL(0.0)
C
C Allocate Vxyo arrays
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
C
C Allocate percref (reference pdf percentiles) vector
      allocate ( percref(1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percref(:) = FREAL(0.0)
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
C Read reference pdf percentiles
      IF (nprint.GE.2) THEN
          WRITE(numout,*) '    ==> READING reference pdf percentiles'
      ENDIF
C
      CALL readscalbas(kinxbasref,'percref',percref)
C
      IF (jpproc.GT.limjpnxyo(MIN(3,kflagxyo))) GOTO 1003
      DO jnxyo=1+jproc,limjpnxyo(MIN(3,kflagxyo)),jpproc
C
C -1.- Read percentiles
C ---------------------
C
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING percentiles'
        ENDIF
C
        jrbasdeb=1
        jrbasfin=jpperc
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxbasref,percsr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxbasref,percsr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
C     
C -4.- Read input vector
C ----------------------
C
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readxyo (kinxyo,vects(:),jnxyo,lectinfo,kflagxyo)
        CASE (3)
          CALL readxyo (kinxyo,vects(:),jnxyo,lectinfo,kflagxyo,
     $                  poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT

C
C -3.- Compute anamorphosis
C -------------------------
C
        IF (traditional) THEN
C TRADITIONAL SCHEME
C ==================
        IF (kflaganam) THEN
          DO js=1,jpssize
C
            lmodprint=(MOD(js-1,(jpssize/5+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprint))) print *,
     $        'Variable index : ',js,'/',jpssize
C
            jp=locperc_C(vects(js),percsr(js,1:jpperc),eql)
C
            IF (jp.EQ.0) THEN
              jp=1        ; a=0.0
            ELSEIF (jp.EQ.jpperc) THEN
              jp=jpperc-1 ; a=1.0
            ELSEIF (eql.EQ.0) THEN
              a=0.0
            ELSEIF (eql.EQ.1) THEN
              a=0.5
            ELSE
              a=(vects(js)      -percsr(js,jp))/
     $          (percsr(js,jp+1)-percsr(js,jp))
            ENDIF
C
            vects(js)=percref(jp)+a*(percref(jp+1)-percref(jp))
C
          ENDDO
        ELSE
          DO js=1,jpssize
C
            lmodprint=(MOD(js-1,(jpssize/5+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprint))) print *,
     $        'Variable index : ',js,'/',jpssize
C
            jp=locperc_B(vects(js),percref(1:jpperc))
C
            IF (jp.LE.0) THEN
              jp=1        ; a=0.0
            ELSEIF (jp.GE.jpperc) THEN
              jp=jpperc-1 ; a=1.0
            ELSE
              a=(vects(js)    -percref(jp))/
     $          (percref(jp+1)-percref(jp))
            ENDIF
C           
            vects(js)=percsr(js,jp)+a*(percsr(js,jp+1)-percsr(js,jp))
C
          ENDDO
        ENDIF
C NEW SCHEME
C ==================
        ELSE
        IF (kflaganam) THEN
          CALL ana_forward(vects,percsr,percref)
        ELSE
          CALL ana_backward(vects,percsr,percref)
        ENDIF
C ==================
        ENDIF
C     
C -4.- Write output vector
C ------------------------
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
      IF (allocated(percsr)) deallocate(percsr)
      IF (allocated(percref)) deallocate(percref)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algoanam','algoanam')
 1001 CALL printerror2(0,1001,3,'algoanam','algoanam')
 1003 CALL printerror2(0,1003,3,'algoanam','algoanam')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoanambas (kinxyobas,koutxyobas,kinxbasref,
     $                        kflagxyo,kflaganam,kconfigo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute anamorphosis of an ensemble of vectors
CCC
CCC  Method : Read ensemble and ensemble percentiles from input
CCC  ------   directories (kinxyobas and kinxbasref),
CCC           write the result in output directory (koutxyobas).
CCC
CCC  Input : kinxyobas   : Cxyo input directory
CCC  -----   koutxyobas  : Cxyo output directory
CCC          kinxbasref  : Cxyo percentiles directory
CCC          kflagxyo    : Vector type (1=Vx,2=Vy,3=Vo)
CCC          kflaganam   : anamorphosis direction (T=+ or F=-)
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
      use ensdam_anatra
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kinxyobas,koutxyobas,kinxbasref
      CHARACTER(len=*), intent(in) :: kconfigo
      INTEGER, intent(in) :: kflagxyo
      LOGICAL, intent(in) :: kflaganam
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:,:), allocatable, save :: basesr
      BIGREAL, dimension(:,:), allocatable, save :: percsr
      BIGREAL, dimension(:), allocatable :: percref
      BIGREAL, dimension(:), allocatable :: vectorms
C
      INTEGER :: allocok,jpssize,jpitpsize,jprsize
      INTEGER :: jnxyo,js,jr,jp,jperc,jjproc
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbase,jrbasdeb,jrbasfin,flagcfg,eql
      BIGREAL :: a,b
CC----------------------------------------------------------------------
C
      jprsize=jprend
      jpitpsize=1
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modanam/algoanambas :'
         WRITE(numout,*) '         compute ensemble anamorphosis'
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
C Allocate percref (reference pdf percentiles) vector
      allocate ( percref(1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percref(:) = FREAL(0.0)
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
C Read reference pdf percentiles
      IF (nprint.GE.2) THEN
          WRITE(numout,*) '    ==> READING reference pdf percentiles'
      ENDIF
C
      CALL readscalbas(kinxbasref,'percref',percref)
C
      IF (jpproc.GT.limjpnxyo(MIN(3,kflagxyo))) GOTO 1003
      DO jnxyo=1+jproc,limjpnxyo(MIN(3,kflagxyo)),jpproc
C
C -1.- Read percentiles
C ---------------------
C
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING percentiles'
        ENDIF
C
        jrbasdeb=1
        jrbasfin=jpperc
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxbasref,percsr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxbasref,percsr(:,:),jnxyo,jrbasdeb,jrbasfin,
     $           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
C
C -2.- Read input ensemble
C ------------------------
C
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING the ensemble'
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
C -3.- Compute anamorphosis
C -------------------------
C
        IF (traditional) THEN
C TRADITIONAL SCHEME
C ==================
        IF (kflaganam) THEN
          DO js=1,jpssize
C
          lmodprint=(MOD(js-1,(jpssize/5+1)).EQ.0)
          IF ((nprint.GE.1).AND.((lmodprint))) print *,
     $        'Variable index : ',js,'/',jpssize
C
          DO jr=1,jprsize
C
            jp=locperc_C(basesr(js,jr),percsr(js,1:jpperc),eql)
C
            IF (jp.EQ.0) THEN
              jp=1        ; a=0.0
            ELSEIF (jp.EQ.jpperc) THEN
              jp=jpperc-1 ; a=1.0
            ELSEIF (eql.EQ.0) THEN
              a=0.0
            ELSEIF (eql.EQ.1) THEN
              a=0.5
            ELSE
              a=(basesr(js,jr)  -percsr(js,jp))/
     $          (percsr(js,jp+1)-percsr(js,jp))
            ENDIF
C
            basesr(js,jr)=percref(jp)+a*(percref(jp+1)-percref(jp))
C
          ENDDO
          ENDDO
        ELSE
          DO js=1,jpssize
C
          lmodprint=(MOD(js-1,(jpssize/5+1)).EQ.0)
          IF ((nprint.GE.1).AND.((lmodprint))) print *,
     $        'Variable index : ',js,'/',jpssize
C
          DO jr=1,jprsize
C
            jp=locperc_B(basesr(js,jr),percref(1:jpperc))
C
            IF (jp.LE.0) THEN
              jp=1        ; a=0.0
            ELSEIF (jp.GE.jpperc) THEN
              jp=jpperc-1 ; a=1.0
            ELSE
              a=(basesr(js,jr)-percref(jp))/
     $          (percref(jp+1)-percref(jp))
            ENDIF
C           
            basesr(js,jr)=percsr(js,jp)+a*(percsr(js,jp+1)-percsr(js,jp))
C
          ENDDO
          ENDDO
        ENDIF
C NEW SCHEME
C ==================
        ELSE
        IF (kflaganam) THEN
          CALL ana_forward(basesr,percsr,percref)
        ELSE
          CALL ana_backward(basesr,percsr,percref)
        ENDIF
C ==================
        ENDIF
C     
C -4.- Write output ensemble
C --------------------------
C
        SELECT CASE (kflagxyo)
        CASE (1)
          CALL writebas(koutxyobas,basesr(:,:),
     $                 jnxyo,jrbasdeb,jrbasfin)
        CASE (2)
          CALL writeyobas(koutxyobas,basesr(:,:),jrbasdeb,jrbasfin,
     $                    kflagxyo)
        CASE (3)
          CALL writeyobas(koutxyobas,basesr(:,:),jrbasdeb,jrbasfin,
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
      IF (allocated(percref)) deallocate(percref)
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
 1003 CALL printerror2(0,1003,3,'algoanam','algoanambas')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoanamobs (kinxyo,kinxyobas,koutyobas,kconfigo,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Ensemble anamorphosis of observation vector
CCC
CCC  Method : Read input observation (kinxyo) and input ensemble (kinxyobas)
CCC  ------   write the sample of transformed observations (koutyobas).
CCC
CCC  Input : kinxyo      : Vo input file
CCC  -----   kinxyobas   : Co input directory
CCC          koutyobas    : Co output directory
CCC          kconfigo    : Observation operator (Vo) input file
CCC          kflagxyo    : Vector type (2=Vy,3=Vo)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilroa, only : mkyorms
      use mod_spacexyo , only :
     $     jpoend,jpyend,jprend,jpsmplend,jpitpend,jpperc,
     $     poscoefobs,gridijkobs
      use hioxyo
      use hiobas
      use utilvalid
      use ensdam_anaobs
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kinxyo,kinxyobas,koutyobas
      CHARACTER(len=*), intent(in) :: kconfigo
      INTEGER, intent(in) :: kflagxyo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:,:), allocatable, save :: enso
      BIGREAL, dimension(:,:), allocatable, save :: smplo
      BIGREAL, dimension(:), allocatable :: percref
      BIGREAL, dimension(:), allocatable :: percdef
      BIGREAL, dimension(:), allocatable, save :: oestd
      BIGREAL, dimension(:), allocatable :: vectorms
C
      INTEGER :: allocok,jpssize,jpitpsize,jprsize
      INTEGER :: js,jr,jjproc,jnxyo,jsmp,jpsmp
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbase,jrbasdeb,jrbasfin,flagcfg
      BIGREAL, DIMENSION(1) :: rjpperc
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jprbas
CC----------------------------------------------------------------------
C
      jprsize=jprend
      jpsmp=jpsmplend
C
      SELECT CASE (kflagxyo)
      CASE (2)
         jpssize=jpyend
         jpitpsize=1
      CASE (3)
         jpssize=jpoend
         jpitpsize=jpitpend
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Read number of quantiles
      CALL readscalbas(koutyobas,'percnbr',rjpperc)
      jpperc = NINT(rjpperc(1))
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modanam/algoanamobs :'
         WRITE(numout,*) '         ensemble anamorphosis of observation'
      ENDIF
C
C Allocate Vo array
      allocate ( vecto(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecto(:) = FREAL(0.0)
C
C Allocate Co array
      allocate ( enso(1:jpssize,1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      enso(:,:) = FREAL(0.0)
C
C Allocate Co array
      allocate ( smplo(1:jpssize,1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      smplo(:,:) = FREAL(0.0)
C
C Allocate percdef (percentile definition) vector
      allocate ( percdef(1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percdef(:) = FREAL(0.0)
C
C Allocate percref (reference pdf percentiles) vector
      allocate ( percref(1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percref(:) = FREAL(0.0)
C
C Allocate Vo array
      allocate ( oestd(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      oestd(:) = FREAL(0.0)
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
C
      ENDIF
C
C Read percentiles definition
      CALL readscalbas(koutyobas,'percdef',percdef)
C
C Read reference pdf percentiles
      CALL readscalbas(koutyobas,'percref',percref)
C
C Read input ensemble
      jnxyo=1
      jrbasdeb=1
      jrbasfin=jprsize
      lectinfo=.FALSE.
      SELECT CASE (kflagxyo)
      CASE (2)
        CALL readbas(kinxyobas,enso(:,:),jnxyo,jrbasdeb,jrbasfin,
     $             lectinfo,kflagxyo)
      CASE (3)
        CALL readbas(kinxyobas,enso(:,:),jnxyo,jrbasdeb,jrbasfin,
     $             lectinfo,kflagxyo,poscoefobs(:,:))
      END SELECT

C
C Read observation vector
      SELECT CASE (kflagxyo)
      CASE (2)
        CALL readxyo (kinxyo,vecto(:),jnxyo,lectinfo,kflagxyo)
      CASE (3)
        CALL readxyo (kinxyo,vecto(:),jnxyo,lectinfo,kflagxyo,
     $                  poscoefobs(:,:))
      END SELECT
C
C Read observation spread of observation error
      IF (largoestd) THEN
        SELECT CASE (kflagxyo)
        CASE (2)
          CALL readxyo(argoestd,oestd(:),
     $           jnxyo,lectinfo,kflagxyo)
        CASE (3)
          CALL readxyo(argoestd,oestd(:),
     $           jnxyo,lectinfo,kflagxyo,poscoefobs(:,:))
        END SELECT
      ELSE
        CALL mkyorms (oestd(:),kflagxyo)
      ENDIF
C
C Define type of observation error
      CALL ana_obs_init( obserror_cdf )
C
      SELECT CASE (obserror_cdf)
      CASE ('gaussian')
        print *, 'Warning: non-optimal algorithm for a Gaussian pdf'
      CASE ('gamma','beta')
      CASE DEFAULT
        print *, 'Warning: unknown observation error pdf'
      END SELECT
C
      DO jsmp=1,jpsmp
C
C Perform anamorphosis transformation
       CALL ana_obs(smplo,enso,vecto,oestd,percdef,percref)
C     
C Write output sample
        CALL fildirbas (vctnam,koutyobas,jprbas,jsmp,1)
        WRITE(fname,'("./",A,"/",A)') koutyobas(1:lenv(koutyobas)),
     $        vctnam(1:lenv(vctnam))
C
        SELECT CASE (kflagxyo)
        CASE (2)
          CALL writeyo(fname,smplo(:,1),kflagxyo)
        CASE (3)
          CALL writeyo(fname,smplo(:,1),kflagxyo,
     $         kvectsrms=vectorms,kgridijkobs=gridijkobs,
     $         kposcoefobs=poscoefobs)
       END SELECT
C
      ENDDO
C
C --- deallocation
      IF (allocated(enso)) deallocate(enso)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algoanam','algoanamobs')
 1001 CALL printerror2(0,1001,3,'algoanam','algoanamobs')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION locperc_A(v,perc)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Localize value in the percentile sequence
CCC            output = first perc if v <  min
CCC                   = last perc  if max <= v
CCC                   = perc index if inf <= v < sup
CCC
CCC  Method : search sequentially in the list of intervals
CCC  ------
CCC  Input : v    : value to localize
CCC  -----   perc : percentile sequence
CCC
CCC---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
C
      BIGREAL, intent(in) :: v
      BIGREAL, dimension(:), intent(in) :: perc
C-----------------------------------------------------------------------
      INTEGER :: jp,jperc,jpperc
C-----------------------------------------------------------------------
      jpperc=size(perc,1)
C
      jp=0
      DO jperc=1,jpperc
        IF (v.GE.perc(jperc)) THEN
         jp=jp+1
        ELSE
         EXIT
        ENDIF
      ENDDO
C
      locperc_A=jp
C
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION locperc_B(v,perc)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Localize value in the percentile sequence
CCC            output = first perc if v <  min
CCC                   = last perc  if max <= v
CCC                   = perc index if inf <= v < sup
CCC
CCC  Method : bisection in the percentile sequence
CCC  ------
CCC  Input : v    : value to localize
CCC  -----   perc : percentile sequence
CCC
CCC---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
C
      BIGREAL, intent(in) :: v
      BIGREAL, dimension(:), intent(in) :: perc
C-----------------------------------------------------------------------
      INTEGER :: jp,jp0,jp1,jperc,jpperc
C-----------------------------------------------------------------------
      jpperc=size(perc,1)
C
      jp0=0 ; jp1=jpperc+1 ; jp=(jp0+jp1)/2
      DO WHILE (jp0.NE.jp)
        IF (v.GE.perc(jp)) THEN
         jp0=jp
        ELSE
         jp1=jp
        ENDIF
        jp=(jp0+jp1)/2
      ENDDO
C
      locperc_B=jp
C
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION locperc_C(v,perc,eql)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Localize value in the percentile sequence
CCC            output = first perc if v <  min
CCC                   = last perc  if max <= v
CCC                   = median perc index such that inf <= v <= sup
CCC
CCC  Method : bisection in the percentile sequence
CCC  ------   use central interval in list of equal percentiles
CCC
CCC  Input : v    : value to localize
CCC  -----   perc : percentile sequence
CCC          eql  :-1 if inf an sup percentiles are diffrent
CCC                 0 for equal percentiles (with even rank difference)
CCC                 1 for equal percentiles (with odd rank difference)
CCC
CCC---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
C
      BIGREAL, intent(in) :: v
      BIGREAL, dimension(:), intent(in) :: perc
      INTEGER, intent(out) :: eql
C-----------------------------------------------------------------------
      INTEGER :: jp,jp0,jp1,jperc,jpperc,jpsup,jpinf
C-----------------------------------------------------------------------
      jpperc=size(perc,1)
C
      jp0=0 ; jp1=jpperc+1 ; jp=(jp0+jp1)/2
      DO WHILE (jp0.NE.jp)
        IF (v.GE.perc(jp)) THEN
         jp0=jp
        ELSE
         jp1=jp
        ENDIF
        jp=(jp0+jp1)/2
      ENDDO
      jpsup=jp
C
      jp0=0 ; jp1=jpperc+1 ; jp=(jp0+jp1)/2
      DO WHILE (jp0.NE.jp)
        IF (v.GT.perc(jp)) THEN
         jp0=jp
        ELSE
         jp1=jp
        ENDIF
        jp=(jp0+jp1)/2
      ENDDO
      jpinf=jp
C
      eql=-1
      IF (jpinf.NE.jpsup) THEN
        jpsup=MAX(jpsup,1) ; jpinf=MAX(jpinf,1)
        jp=(jpinf+jpsup)/2
        IF (jp.LT.jpperc) THEN
        IF (perc(jp).GE.perc(jp+1)) THEN
          eql=MOD(jpsup-jpinf,2)
        ENDIF
        ENDIF
      ELSE
        jp=jpsup
      ENDIF
C
      locperc_C=jp
C
      RETURN
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algoanam
