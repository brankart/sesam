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
! ---                   ALGOANAM.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2008-03 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE algoanam
! --- SUBROUTINE algoanambas
! --- SUBROUTINE algoanamobs
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algoanam
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC algoanamvct, algoanambas, algoanamobs

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoanamvct (kinxyo,koutxyo,kinxbasref, &
     &                        kflagxyo,kflaganam,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose : Compute anamorphosis of a vector of variables
!
!  Method : Read ensemble percentiles from input
!  ------   directory (kinxbasref), read vector of variables
!           from input file (kinxyo)
!           write the result in output file (koutxyo).
!
!  Input : kinxyo      : Vxyo input file
!  -----   koutxyo     : Vxyo output file
!          kinxbasref  : Cxyo percentiles directory
!          kflagxyo    : Vector type (1=Vx,2=Vy,3=Vo)
!          kflaganam   : anamorphosis direction (T=+ or F=-)
!          kconfigo    : Observation operator (Vo) inputy file
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : &
     &     jpoend,jpitpend,jpx,jpxend,jpyend,jprend,jpperc, &
     &     poscoefobs,gridijkobs,arraynx_jindxbeg
      use hioxyo
      use hiobas
      use ensdam_anatra
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinxyo,koutxyo,kinxbasref
      CHARACTER(len=*), intent(in) :: kconfigo
      INTEGER, intent(in) :: kflagxyo
      LOGICAL, intent(in) :: kflaganam
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:,:), allocatable, save :: percsr
      BIGREAL, dimension(:), allocatable :: percref
      BIGREAL, dimension(:), allocatable :: vectorms
!
      INTEGER :: allocok,jpssize,jpitpsize,jprsize
      INTEGER :: jnxyo,js,jr,jperc,jp,jjproc
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbase,jrbasdeb,jrbasfin,flagcfg,eql
      BIGREAL :: a,b,anaval
!----------------------------------------------------------------------
!
      jprsize=jprend
      jpitpsize=1
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modanam/algoanam :'
         WRITE(numout,*) '         compute ensemble anamorphosis'
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
      allocate ( percsr(1:jpssize,1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percsr(:,:) = FREAL(0.0)
!
! Allocate Vxyo arrays
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
!
! Allocate percref (reference pdf percentiles) vector
      allocate ( percref(1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percref(:) = FREAL(0.0)
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
! Read reference pdf percentiles
      IF (nprint.GE.2) THEN
          WRITE(numout,*) '    ==> READING reference pdf percentiles'
      ENDIF
!
      CALL readscalbas(kinxbasref,'percref',percref)
!
      IF (jpproc.GT.limjpnxyo(MIN(3,kflagxyo))) GOTO 1003
      DO jnxyo=1+jproc,limjpnxyo(MIN(3,kflagxyo)),jpproc
!
! -1.- Read percentiles
! ---------------------
!
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING percentiles'
        ENDIF
!
        jrbasdeb=1
        jrbasfin=jpperc
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxbasref,percsr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxbasref,percsr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
!     
! -4.- Read input vector
! ----------------------
!
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readxyo (kinxyo,vects(:),jnxyo,lectinfo,kflagxyo)
        CASE (3)
          CALL readxyo (kinxyo,vects(:),jnxyo,lectinfo,kflagxyo, &
     &                  poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT

!
! -3.- Compute anamorphosis
! -------------------------
!
        IF (traditional) THEN
! TRADITIONAL SCHEME
! ==================
        IF (kflaganam) THEN
          DO js=1,jpssize
!
            lmodprint=(MOD(js-1,(jpssize/5+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprint))) print *, &
     &        'Variable index : ',js,'/',jpssize
!
            jp=locperc_C(vects(js),percsr(js,1:jpperc),eql)
!
            IF (jp.EQ.0) THEN
              jp=1        ; a=0.0
            ELSEIF (jp.EQ.jpperc) THEN
              jp=jpperc-1 ; a=1.0
            ELSEIF (eql.EQ.0) THEN
              a=0.0
            ELSEIF (eql.EQ.1) THEN
              a=0.5
            ELSE
              a=(vects(js)      -percsr(js,jp))/ &
     &          (percsr(js,jp+1)-percsr(js,jp))
            ENDIF
!
            vects(js)=percref(jp)+a*(percref(jp+1)-percref(jp))
!
          ENDDO
        ELSE
          DO js=1,jpssize
!
            lmodprint=(MOD(js-1,(jpssize/5+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprint))) print *, &
     &        'Variable index : ',js,'/',jpssize
!
            jp=locperc_B(vects(js),percref(1:jpperc))
!
            IF (jp.LE.0) THEN
              jp=1        ; a=0.0
            ELSEIF (jp.GE.jpperc) THEN
              jp=jpperc-1 ; a=1.0
            ELSE
              a=(vects(js)    -percref(jp))/ &
     &          (percref(jp+1)-percref(jp))
            ENDIF
!           
            vects(js)=percsr(js,jp)+a*(percsr(js,jp+1)-percsr(js,jp))
!
          ENDDO
        ENDIF
! NEW SCHEME
! ==================
        ELSE
        IF (kflaganam) THEN
          WHERE(vects==special_value) vects=huge(vects)
          CALL ana_forward(vects,percsr,percref)
          WHERE(vects==huge(vects)) vects=special_value
        ELSE
          CALL ana_backward(vects,percsr,percref)
        ENDIF
! ==================
        ENDIF
!     
! -4.- Write output vector
! ------------------------
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
      IF (allocated(percsr)) deallocate(percsr)
      IF (allocated(percref)) deallocate(percref)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algoanam','algoanam')
 1001 CALL printerror2(0,1001,3,'algoanam','algoanam')
 1003 CALL printerror2(0,1003,3,'algoanam','algoanam')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoanambas (kinxyobas,koutxyobas,kinxbasref, &
     &                        kflagxyo,kflaganam,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose : Compute anamorphosis of an ensemble of vectors
!
!  Method : Read ensemble and ensemble percentiles from input
!  ------   directories (kinxyobas and kinxbasref),
!           write the result in output directory (koutxyobas).
!
!  Input : kinxyobas   : Cxyo input directory
!  -----   koutxyobas  : Cxyo output directory
!          kinxbasref  : Cxyo percentiles directory
!          kflagxyo    : Vector type (1=Vx,2=Vy,3=Vo)
!          kflaganam   : anamorphosis direction (T=+ or F=-)
!          kconfigo    : Observation operator (Vo) inputy file
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : &
     &     jpoend,jpitpend,jpx,jpxend,jpyend,jprend,jpperc, &
     &     poscoefobs,gridijkobs,arraynx_jindxbeg
      use hioxyo
      use hiobas
      use ensdam_anatra
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinxyobas,koutxyobas,kinxbasref
      CHARACTER(len=*), intent(in) :: kconfigo
      INTEGER, intent(in) :: kflagxyo
      LOGICAL, intent(in) :: kflaganam
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:,:), allocatable, save :: basesr
      BIGREAL, dimension(:,:), allocatable, save :: percsr
      BIGREAL, dimension(:), allocatable :: percref
      BIGREAL, dimension(:), allocatable :: vectorms
!
      INTEGER :: allocok,jpssize,jpitpsize,jprsize
      INTEGER :: jnxyo,js,jr,jp,jperc,jjproc
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbase,jrbasdeb,jrbasfin,flagcfg,eql
      BIGREAL :: a,b
!----------------------------------------------------------------------
!
      jprsize=jprend
      jpitpsize=1
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modanam/algoanambas :'
         WRITE(numout,*) '         compute ensemble anamorphosis'
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
! Allocate Cxyo array
      allocate ( basesr(1:jpssize,1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basesr(:,:) = FREAL(0.0)
!
! Allocate Cxyo percentile array
      allocate ( percsr(1:jpssize,1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percsr(:,:) = FREAL(0.0)
!
! Allocate percref (reference pdf percentiles) vector
      allocate ( percref(1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percref(:) = FREAL(0.0)
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
! Read reference pdf percentiles
      IF (nprint.GE.2) THEN
          WRITE(numout,*) '    ==> READING reference pdf percentiles'
      ENDIF
!
      CALL readscalbas(kinxbasref,'percref',percref)
!
      IF (jpproc.GT.limjpnxyo(MIN(3,kflagxyo))) GOTO 1003
      DO jnxyo=1+jproc,limjpnxyo(MIN(3,kflagxyo)),jpproc
!
! -1.- Read percentiles
! ---------------------
!
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING percentiles'
        ENDIF
!
        jrbasdeb=1
        jrbasfin=jpperc
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxbasref,percsr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxbasref,percsr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
!
! -2.- Read input ensemble
! ------------------------
!
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING the ensemble'
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
! -3.- Compute anamorphosis
! -------------------------
!
        IF (traditional) THEN
! TRADITIONAL SCHEME
! ==================
        IF (kflaganam) THEN
          DO js=1,jpssize
!
          lmodprint=(MOD(js-1,(jpssize/5+1)).EQ.0)
          IF ((nprint.GE.1).AND.((lmodprint))) print *, &
     &        'Variable index : ',js,'/',jpssize
!
          DO jr=1,jprsize
!
            jp=locperc_C(basesr(js,jr),percsr(js,1:jpperc),eql)
!
            IF (jp.EQ.0) THEN
              jp=1        ; a=0.0
            ELSEIF (jp.EQ.jpperc) THEN
              jp=jpperc-1 ; a=1.0
            ELSEIF (eql.EQ.0) THEN
              a=0.0
            ELSEIF (eql.EQ.1) THEN
              a=0.5
            ELSE
              a=(basesr(js,jr)  -percsr(js,jp))/ &
     &          (percsr(js,jp+1)-percsr(js,jp))
            ENDIF
!
            basesr(js,jr)=percref(jp)+a*(percref(jp+1)-percref(jp))
!
          ENDDO
          ENDDO
        ELSE
          DO js=1,jpssize
!
          lmodprint=(MOD(js-1,(jpssize/5+1)).EQ.0)
          IF ((nprint.GE.1).AND.((lmodprint))) print *, &
     &        'Variable index : ',js,'/',jpssize
!
          DO jr=1,jprsize
!
            jp=locperc_B(basesr(js,jr),percref(1:jpperc))
!
            IF (jp.LE.0) THEN
              jp=1        ; a=0.0
            ELSEIF (jp.GE.jpperc) THEN
              jp=jpperc-1 ; a=1.0
            ELSE
              a=(basesr(js,jr)-percref(jp))/ &
     &          (percref(jp+1)-percref(jp))
            ENDIF
!           
            basesr(js,jr)=percsr(js,jp)+a*(percsr(js,jp+1)-percsr(js,jp))
!
          ENDDO
          ENDDO
        ENDIF
! NEW SCHEME
! ==================
        ELSE
        IF (kflaganam) THEN
          WHERE(basesr==special_value) basesr=huge(basesr)
          CALL ana_forward(basesr,percsr,percref)
          WHERE(basesr==huge(basesr)) basesr=special_value
        ELSE
          CALL ana_backward(basesr,percsr,percref)
        ENDIF
! ==================
        ENDIF
!     
! -4.- Write output ensemble
! --------------------------
!
        SELECT CASE (kflagxyo)
        CASE (1)
          CALL writebas(koutxyobas,basesr(:,:), &
     &                 jnxyo,jrbasdeb,jrbasfin)
        CASE (2)
          CALL writeyobas(koutxyobas,basesr(:,:),jrbasdeb,jrbasfin, &
     &                    kflagxyo)
        CASE (3)
          CALL writeyobas(koutxyobas,basesr(:,:),jrbasdeb,jrbasfin, &
     &                    kflagxyo,vectorms(:),gridijkobs(:), &
     &                    poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
!
      ENDDO
!
! --- deallocation
      IF (allocated(basesr)) deallocate(basesr)
      IF (allocated(percsr)) deallocate(percsr)
      IF (allocated(percref)) deallocate(percref)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algoanam','algoanambas')
 1001 CALL printerror2(0,1001,3,'algoanam','algoanambas')
 1003 CALL printerror2(0,1003,3,'algoanam','algoanambas')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoanamobs (kinxyo,kinxyobas,koutyobas,kconfigo,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Ensemble anamorphosis of observation vector
!
!  Method : Read input observation (kinxyo) and input ensemble (kinxyobas)
!  ------   write the sample of transformed observations (koutyobas).
!
!  Input : kinxyo      : Vo input file
!  -----   kinxyobas   : Co input directory
!          koutyobas    : Co output directory
!          kconfigo    : Observation operator (Vo) input file
!          kflagxyo    : Vector type (2=Vy,3=Vo)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilroa, only : mkyorms
      use mod_spacexyo , only : &
     &     jpoend,jpyend,jprend,jpsmplend,jpitpend,jpperc, &
     &     poscoefobs,gridijkobs
      use hioxyo
      use hiobas
      use utilvalid
      use ensdam_anaobs
      use ensdam_obserror
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinxyo,kinxyobas,koutyobas
      CHARACTER(len=*), intent(in) :: kconfigo
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:,:), allocatable, save :: enso
      BIGREAL, dimension(:,:), allocatable, save :: smplo
      BIGREAL, dimension(:), allocatable :: percref
      BIGREAL, dimension(:), allocatable :: percdef
      BIGREAL, dimension(:), allocatable, save :: oestd
      BIGREAL, dimension(:), allocatable :: vectorms
!
      INTEGER :: allocok,jpssize,jpitpsize,jprsize
      INTEGER :: js,jr,jjproc,jnxyo,jsmp,jpsmp
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbase,jrbasdeb,jrbasfin,flagcfg
      BIGREAL, DIMENSION(1) :: rjpperc
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jprbas
!----------------------------------------------------------------------
!
      jprsize=jprend
      jpsmp=jpsmplend
!
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
!
! Read number of quantiles
      CALL readscalbas(koutyobas,'percnbr',rjpperc)
      jpperc = NINT(rjpperc(1))
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modanam/algoanamobs :'
         WRITE(numout,*) '         ensemble anamorphosis of observation'
      ENDIF
!
! Allocate Vo array
      allocate ( vecto(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecto(:) = FREAL(0.0)
!
! Allocate Co array
      allocate ( enso(1:jpssize,1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      enso(:,:) = FREAL(0.0)
!
! Allocate Co array
      allocate ( smplo(1:jpssize,1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      smplo(:,:) = FREAL(0.0)
!
! Allocate percdef (percentile definition) vector
      allocate ( percdef(1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percdef(:) = FREAL(0.0)
!
! Allocate percref (reference pdf percentiles) vector
      allocate ( percref(1:jpperc), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      percref(:) = FREAL(0.0)
!
! Allocate Vo array
      allocate ( oestd(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      oestd(:) = FREAL(0.0)
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
!
      ENDIF
!
! Read percentiles definition
      CALL readscalbas(koutyobas,'percdef',percdef)
!
! Read reference pdf percentiles
      CALL readscalbas(koutyobas,'percref',percref)
!
! Read input ensemble
      jnxyo=1
      jrbasdeb=1
      jrbasfin=jprsize
      lectinfo=.FALSE.
      SELECT CASE (kflagxyo)
      CASE (2)
        CALL readbas(kinxyobas,enso(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &             lectinfo,kflagxyo)
      CASE (3)
        CALL readbas(kinxyobas,enso(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &             lectinfo,kflagxyo,poscoefobs(:,:))
      END SELECT

!
! Read observation vector
      SELECT CASE (kflagxyo)
      CASE (2)
        CALL readxyo (kinxyo,vecto(:),jnxyo,lectinfo,kflagxyo)
      CASE (3)
        CALL readxyo (kinxyo,vecto(:),jnxyo,lectinfo,kflagxyo, &
     &                  poscoefobs(:,:))
      END SELECT
!
! Read observation spread of observation error
      IF (largoestd) THEN
        SELECT CASE (kflagxyo)
        CASE (2)
          CALL readxyo(argoestd,oestd(:), &
     &           jnxyo,lectinfo,kflagxyo)
        CASE (3)
          CALL readxyo(argoestd,oestd(:), &
     &           jnxyo,lectinfo,kflagxyo,poscoefobs(:,:))
        END SELECT
      ELSE
        CALL mkyorms (oestd(:),kflagxyo)
      ENDIF
!
! Define type of observation error
      obserror_type = obserror_type_sesam
!
      SELECT CASE (obserror_type_sesam)
      CASE ('gaussian')
        print *, 'Warning: non-optimal algorithm for a Gaussian pdf'
      CASE ('gamma','beta')
      CASE DEFAULT
        print *, 'Warning: unknown observation error pdf'
      END SELECT
!
      DO jsmp=1,jpsmp
!
! Perform anamorphosis transformation
       CALL ana_obs(smplo,enso,vecto,oestd,percdef,percref)
!     
! Write output sample
        CALL fildirbas (vctnam,koutyobas,jprbas,jsmp,1)
        WRITE(fname,'("./",A,"/",A)') koutyobas(1:lenv(koutyobas)), &
     &        vctnam(1:lenv(vctnam))
!
        SELECT CASE (kflagxyo)
        CASE (2)
          CALL writeyo(fname,smplo(:,1),kflagxyo)
        CASE (3)
          CALL writeyo(fname,smplo(:,1),kflagxyo, &
     &         kvectsrms=vectorms,kgridijkobs=gridijkobs, &
     &         kposcoefobs=poscoefobs)
       END SELECT
!
      ENDDO
!
! --- deallocation
      IF (allocated(enso)) deallocate(enso)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algoanam','algoanamobs')
 1001 CALL printerror2(0,1001,3,'algoanam','algoanamobs')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION locperc_A(v,perc)
!---------------------------------------------------------------------
!
!  Purpose : Localize value in the percentile sequence
!            output = first perc if v <  min
!                   = last perc  if max <= v
!                   = perc index if inf <= v < sup
!
!  Method : search sequentially in the list of intervals
!  ------
!  Input : v    : value to localize
!  -----   perc : percentile sequence
!
!---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
!
      BIGREAL, intent(in) :: v
      BIGREAL, dimension(:), intent(in) :: perc
!-----------------------------------------------------------------------
      INTEGER :: jp,jperc,jpperc
!-----------------------------------------------------------------------
      jpperc=size(perc,1)
!
      jp=0
      DO jperc=1,jpperc
        IF (v.GE.perc(jperc)) THEN
         jp=jp+1
        ELSE
         EXIT
        ENDIF
      ENDDO
!
      locperc_A=jp
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION locperc_B(v,perc)
!---------------------------------------------------------------------
!
!  Purpose : Localize value in the percentile sequence
!            output = first perc if v <  min
!                   = last perc  if max <= v
!                   = perc index if inf <= v < sup
!
!  Method : bisection in the percentile sequence
!  ------
!  Input : v    : value to localize
!  -----   perc : percentile sequence
!
!---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
!
      BIGREAL, intent(in) :: v
      BIGREAL, dimension(:), intent(in) :: perc
!-----------------------------------------------------------------------
      INTEGER :: jp,jp0,jp1,jperc,jpperc
!-----------------------------------------------------------------------
      jpperc=size(perc,1)
!
      jp0=0 ; jp1=jpperc+1 ; jp=(jp0+jp1)/2
      DO WHILE (jp0.NE.jp)
        IF (v.GE.perc(jp)) THEN
         jp0=jp
        ELSE
         jp1=jp
        ENDIF
        jp=(jp0+jp1)/2
      ENDDO
!
      locperc_B=jp
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION locperc_C(v,perc,eql)
!---------------------------------------------------------------------
!
!  Purpose : Localize value in the percentile sequence
!            output = first perc if v <  min
!                   = last perc  if max <= v
!                   = median perc index such that inf <= v <= sup
!
!  Method : bisection in the percentile sequence
!  ------   use central interval in list of equal percentiles
!
!  Input : v    : value to localize
!  -----   perc : percentile sequence
!          eql  :-1 if inf an sup percentiles are diffrent
!                 0 for equal percentiles (with even rank difference)
!                 1 for equal percentiles (with odd rank difference)
!
!---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
!
      BIGREAL, intent(in) :: v
      BIGREAL, dimension(:), intent(in) :: perc
      INTEGER, intent(out) :: eql
!-----------------------------------------------------------------------
      INTEGER :: jp,jp0,jp1,jperc,jpperc,jpsup,jpinf
!-----------------------------------------------------------------------
      jpperc=size(perc,1)
!
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
!
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
!
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
!
      locperc_C=jp
!
      RETURN
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algoanam
