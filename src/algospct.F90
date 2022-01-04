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
! ---                   ALGOSPCT.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2016-06 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE algospctvct
! --- SUBROUTINE algospctbas
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algospct
      use mod_main
      use ensdam_storfg
      use ensdam_storng
      IMPLICIT NONE
      PRIVATE

      PUBLIC algospctvct, algospctbas

      ! Public variables
      BIGREAL, PUBLIC, save :: loc_time_scale=1.    !  Localization time scale
      BIGREAL, PUBLIC, save :: loc_radius_in_deg=1. !  Localization radius

      ! Private variables/parameters
      INTEGER, save :: jpl  ! maximum degree of Legendre polynomials

      INTEGER, parameter :: jptspct=1000  ! Number of point discretizing time power spectrum
      INTEGER, parameter :: jpsmpl=1000   ! Size of sample of frequencies
      BIGREAL, parameter :: pi=3.14159265358979

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algospctvct (kinxyo,koutxyo, &
     &                        kflagxyo,kargtypeoper,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose : Compute transformation of a vector of variables
!
!  Method : Read vector of variables from input file (kinxyo),
!  ------   compute transformation for each 2D slice of variables,
!           write the result in output file (koutxyo).
!
!  Input : kinxyo      : Vxyo input file
!  -----   koutxyo     : Vxyo output file
!          kflagxyo    : Vector type (1=Vx,2=Vy)
!          kargtypeoper : transformation direction (+ or -) and
!                         maximum degree of spherical harmonics (jpl)
!          kconfigo    : Observation operator (Vo) inputy file
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_coord
      use mod_spacexyo , only : jpx,jpyend, &
     &     jpoend,jpitpend,poscoefobs,gridijkobs
      use utilroa , only : mkyorms
      use hioxyo
      use hiogrd
      use liospct
      use utilvct
      use utilvalid
      use ensdam_sphylm
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinxyo,koutxyo,kargtypeoper
      CHARACTER(len=*), intent(in) :: kconfigo
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      INTEGER :: allocok,jpssize,jpisize,jpjsize,jpitpsize,jptsize
      INTEGER :: jnxyo,js,jr,js0,js1,jsxy,jk,jt,jl,jm,jlmin
      INTEGER :: jampl,jpampl
      LOGICAL :: lectinfo
      INTEGER :: flagcfg,ios
      INTEGER :: sxyend,indsxy,nbr
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_nbr
      INTEGER, dimension(1:nbvar) :: sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt
      CHARACTER(len=varlg), dimension(1:nbvar) :: sxy_nam
      BIGREAL :: latmin, latmax, dlatmax
      BIGREAL, DIMENSION(:), allocatable :: lon, lat
      BIGREAL, DIMENSION(:), allocatable :: wei, obs, obswei
      BIGREAL, DIMENSION(:,:), allocatable :: weightij
      BIGREAL, DIMENSION(:,:), allocatable :: spct, spctstd
      BIGREAL, DIMENSION(:), allocatable :: tspct_freq, tspct_power
      BIGREAL, DIMENSION(:,:), allocatable :: ran_timeseries
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modspct/algospct :'
         WRITE(numout,*) '         compute spectral transformation'
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
      IF ( (kargtypeoper(1:1).EQ.'+').OR. &
     &     (kargtypeoper(1:1).EQ.'-').OR. &
     &     (kargtypeoper(1:1).EQ.'R').OR. &
     &     (kargtypeoper(1:1).EQ.'O') ) THEN
        jlmin = 0
        READ(kargtypeoper((posit(kargtypeoper,'_')+1): &
     &                 lenv(kargtypeoper)),*,IOSTAT=ios) jpl
        IF (ios.NE.0) THEN
          READ(kargtypeoper((posit(kargtypeoper,'_')+1): &
     &                 (posit(kargtypeoper,'~')-1)),*,IOSTAT=ios) jlmin
          IF (ios.NE.0) GOTO 1000
          READ(kargtypeoper((posit(kargtypeoper,'~')+1): &
     &                     lenv(kargtypeoper)),*,IOSTAT=ios) jpl
          IF (ios.NE.0) GOTO 1000
        ENDIF
      ELSEIF (argtypeoper(1:1).EQ.'D') THEN
        READ(kargtypeoper((posit(kargtypeoper,'_')+1): &
     &                   (posit(kargtypeoper,'-')-1)),*,IOSTAT=ios) jl
        IF (ios.NE.0) GOTO 1000
        READ(kargtypeoper((posit(kargtypeoper,'-')+1): &
     &                   lenv(kargtypeoper)),*,IOSTAT=ios) jm
        IF (ios.NE.0) GOTO 1000
        jpl = jl
      ELSE
        GOTO 1000
      ENDIF
!
! Get variable characteristics from SESAM configuration
      sxy_jpi(:) = 0 ; sxy_jpj(:) = 0
!
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         sxyend  = varend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=var_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
            sxy_nam(indsxy)=var_nam(indsxy)
            sxy_nbr(indsxy)=var_nbr(indsxy)
            sxy_jpi(indsxy)=var_jpi(indsxy)
            sxy_jpj(indsxy)=var_jpj(indsxy)
            sxy_jpk(indsxy)=var_jpk(indsxy)
            sxy_jpt(indsxy)=var_jpt(indsxy)
            IF (varngrd(indsxy).LE.2) GOTO 1000
         ENDDO
      CASE(2,3)
! --- dta
         sxyend  = dtaend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=dta_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
            sxy_nam(indsxy)=dta_nam(indsxy)
            sxy_nbr(indsxy)=dta_nbr(indsxy)
            sxy_jpi(indsxy)=dta_jpi(indsxy)
            sxy_jpj(indsxy)=dta_jpj(indsxy)
            sxy_jpk(indsxy)=dta_jpk(indsxy)
            sxy_jpt(indsxy)=dta_jpt(indsxy)
            IF (dtangrd(indsxy).LE.2) GOTO 1000
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Allocate Vxyo arrays
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
!
! Allocate spectrum array
      allocate ( spct(0:jpl,-jpl:jpl), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      spct(:,:) = FREAL(0.0)
!
      IF (kflagxyo.LT.3) THEN
!
! Allocate grid array
        jpisize=MAXVAL(sxy_jpi(:))
        jpjsize=MAXVAL(sxy_jpj(:))
!
        allocate (gridij(1:jpisize,1:jpjsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        gridij(:,:) = type_gridij(FREAL(0.0),FREAL(0.0))
!
        allocate (weightij(1:jpisize,1:jpjsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        weightij(:,:) = FREAL(0.0)
!
        allocate (lon(1:jpisize*jpjsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        lon(:) = FREAL(0.0)
!
        allocate (lat(1:jpisize*jpjsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        lat(:) = FREAL(0.0)
!
        allocate (wei(1:jpisize*jpjsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        wei(:) = FREAL(0.0)
!
      ELSE
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
! Allocate wei array
        allocate ( spctstd(0:jpl,-jpl:jpl), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        spctstd(:,:) = FREAL(0.0)
!
! Read poscoefobs and gridijkobs arrays
        flagcfg=2
        CALL readcfgobs(kconfigo,flagcfg,kgridijkobs=gridijkobs(:))
        flagcfg=3
        CALL readcfgobs(kconfigo,flagcfg,kposcoefobs=poscoefobs(:,:))
!
! Get observation error standard deviation
        jnxyo=1
        IF (largoestd) THEN
          IF ((validextvar(argoestd)).OR.(validextdta(argoestd)) &
     &        .OR.(validextobs(argoestd))) THEN
            CALL readxyo(argoestd,vectorms(:), &
     &           jnxyo,lectinfo,kflagxyo,poscoefobs(:,:))
          ELSE
            GOTO 1000
          ENDIF
        ELSE
          CALL mkyorms (vectorms(:),kflagxyo)
        ENDIF
!
        IF (ANY(vectorms(:).LE.0.0)) GOTO 101
!
        allocate (lon(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        lon(:) = FREAL(0.0)
!
        allocate (lat(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        lat(:) = FREAL(0.0)
!
        allocate (obswei(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        obswei(:) = FREAL(0.0)
!
        allocate (obs(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        obs(:) = FREAL(0.0)
!
      ENDIF
!
! -1.- Read input vector
! ----------------------
!
      jnxyo=1
      IF (kargtypeoper(1:1).EQ.'+') THEN
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readxyo(kinxyo,vects(:),jnxyo,lectinfo,kflagxyo)
        CASE DEFAULT
          GOTO 1000
        END SELECT
      ELSEIF (kargtypeoper(1:1).EQ.'O') THEN
        CALL readxyo (kinxyo,vects(:),jnxyo,lectinfo,kflagxyo, &
     &                poscoefobs(:,:))
      ENDIF
!
! -2.- Compute spectral transformation
! ------------------------------------
!
! Precomputation of Legendre polynomials
      latmin = -90. ; latmax = 90. ; dlatmax = 0.1
      CALL init_ylm( jpl, jlmin, latmin, latmax, dlatmax )
!
      IF (kargtypeoper(1:1).EQ.'+') THEN
!       Direct transformation
        js0 = 1
        DO jsxy=1,sxyend
          indsxy= sxy_ord(jsxy)
          CALL readgrd(kflagxyo,indsxy)
          CALL proj_wei(weightij)
!
          DO jt=1,sxy_jpt(indsxy)
          DO jk=1,sxy_jpk(indsxy)
!
             WRITE(*,'(2a,2i4)') 'Transforming slice: ',sxy_nam(indsxy),jt,jk
!
             CALL mk8vct(lon(1:),gridij(:,:)%longi,jk,jt, &
     &                   jsxy,nbr,kflagxyo)
             CALL mk8vct(lat(1:),gridij(:,:)%latj, jk,jt, &
     &                   jsxy,nbr,kflagxyo)
             CALL mk8vct(wei(1:),weightij(:,:), jk,jt, &
     &                   jsxy,nbr,kflagxyo)
             print *, '...fraction of the sphere',sum(wei(1:nbr))

             js1 = js0 + nbr - 1
             vects(js0:js1) = vects(js0:js1) * wei(1:nbr)
             CALL proj_ylm( spct(0:,-jpl:), vects(js0:js1), &
     &                      lon(1:nbr), lat(1:nbr) )
             CALL writespct(koutxyo,spct,kflagxyo,jsxy,jk,jt)
             js0 = js0 + nbr
!
          ENDDO
          ENDDO
        ENDDO
      ELSEIF (kargtypeoper(1:1).EQ.'-') THEN
!       Inverse transformation
        js0 = 1
        DO jsxy=1,sxyend
          indsxy= sxy_ord(jsxy)
          CALL readgrd(kflagxyo,indsxy)
!
          DO jt=1,sxy_jpt(indsxy)
          DO jk=1,sxy_jpk(indsxy)
!
             WRITE(*,'(2a,2i4)') 'Transforming slice: ',sxy_nam(indsxy),jt,jk
!
             CALL readspct(kinxyo,spct,kflagxyo,jsxy,jk,jt)
             CALL mk8vct(lon(1:),gridij(:,:)%longi,jk,jt, &
     &                   jsxy,nbr,kflagxyo)
             CALL mk8vct(lat(1:),gridij(:,:)%latj, jk,jt, &
     &                   jsxy,nbr,kflagxyo)
!
             js1 = js0 + nbr - 1
             CALL back_ylm( spct(0:,-jpl:), vects(js0:js1), &
     &                      lon(1:nbr), lat(1:nbr) )
             js0 = js0 + nbr
          ENDDO
          ENDDO
        ENDDO
      ELSEIF (kargtypeoper(1:1).EQ.'R') THEN
        IF (kflagxyo.GE.3) GOTO 1000
!       Generate random field with given spectrum in space and time

        jptsize=MAXVAL(sxy_jpt(:))
        jpampl=(jpl+1)*(jpl+1)
        allocate (ran_timeseries(1:jptsize,1:jpampl), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        allocate ( time(1:jptsize) , stat = allocok )
        IF (allocok.GT.0) GOTO 1001
        allocate ( tspct_freq(1:jptspct) , stat = allocok )
        IF (allocok.GT.0) GOTO 1001
        allocate ( tspct_power(1:jptspct) , stat = allocok )
        IF (allocok.GT.0) GOTO 1001
        allocate ( spctstd(0:jpl,-jpl:jpl), stat=allocok )
        IF (allocok.NE.0) GOTO 1001

!       Define power spectrum for spherical harmonics
!       and define time power spectrum from parameters
        CALL define_spct(spctstd,tspct_freq,tspct_power)

!       Draw random frequencies and phases for random timeseries
!       (identical for all levels, all variables)
!       separately for every spectral amplitude (jampl)
        CALL kiss_load()
        CALL def_spect_init(jptspct,jpampl,0,0)
        CALL def_sample_size(jpsmpl,0,0)
        DO jampl=1,jpampl
          CALL def_spect_power(1,jampl,tspct_freq,tspct_power)
          CALL sample_freq_1d(jampl)
        ENDDO
        CALL kiss_save()

        js0 = 1
        DO jsxy=1,sxyend
          indsxy= sxy_ord(jsxy)
          CALL readgrd(kflagxyo,indsxy)
          CALL readtime(kflagxyo,indsxy)

!         Generate random timeseries on required time grid
          DO jampl=1,jpampl
            CALL gen_field_1d(jampl,                                &
     &                  ran_timeseries(1:sxy_jpt(indsxy),jampl), &
     &                  time(1:sxy_jpt(indsxy)))
          ENDDO

!         Loop on time
          DO jt=1,sxy_jpt(indsxy)

            IF (jproc.EQ.0) THEN
              WRITE(*,'(2a,2i4)') 'Generating time slice: ',sxy_nam(indsxy),jt
            ENDIF

!           Generate amplitude of spectral harmonics for time t           
            jampl=1
            DO jl=0,jpl
            DO jm=-jl,jl
              spct(jl,jm)=ran_timeseries(jt,jampl)*spctstd(jl,jm)
              jampl=jampl+1
            ENDDO
            ENDDO

!           Loop on levels
            DO jk=1,sxy_jpk(indsxy)
!
              CALL mk8vct(lon(1:),gridij(:,:)%longi,jk,jt, &
     &                    jsxy,nbr,kflagxyo)
              CALL mk8vct(lat(1:),gridij(:,:)%latj, jk,jt, &
     &                    jsxy,nbr,kflagxyo)
!
              js1 = js0 + nbr - 1
              CALL back_ylm( spct(0:,-jpl:), vects(js0:js1), &
     &                       lon(1:nbr), lat(1:nbr) )
              js0 = js0 + nbr
            ENDDO
          ENDDO

        ENDDO

      ELSEIF (kargtypeoper(1:1).EQ.'O') THEN
!
        CALL init_regr_ylm( regr_type_sesam, regr_maxiter_sesam, regr_maxbloc_sesam, &
     &                      regr_overlap_sesam, regr_epsilon_sesam, regr_rho_sesam )
!
!       Observation transformation
        DO jsxy=1,sxyend
          indsxy= sxy_ord(jsxy)
!
          DO jt=1,sxy_jpt(indsxy)
          DO jk=1,sxy_jpk(indsxy)
!
             WRITE(*,'(2a,2i4)') 'Transforming slice: ',sxy_nam(indsxy),jt,jk
!
             spctstd(:,:) = 1.0
             IF (largweight) THEN
               CALL readspct(argweight,spctstd,kflagxyo,jsxy,jk,jt)
             ENDIF
!
             CALL prep_obs( jsxy, jt, jk, vects, &
     &                      poscoefobs, gridijkobs, vectorms, &
     &                      obs, lon, lat, obswei, nbr )
             CALL regr_ylm( spct(0:,-jpl:), spctstd(0:,-jpl:), obs(1:nbr), &
     &                      lon(1:nbr), lat(1:nbr), obswei(1:nbr) )
             CALL writespct(koutxyo,spct,kflagxyo,jsxy,jk,jt)
!
          ENDDO
          ENDDO
        ENDDO
!
      ELSEIF (kargtypeoper(1:1).EQ.'D') THEN
!       Display one single spherical harmonics
        js0 = 1
        DO jsxy=1,sxyend
          indsxy= sxy_ord(jsxy)
          CALL readgrd(kflagxyo,indsxy)
!
          DO jt=1,sxy_jpt(indsxy)
          DO jk=1,sxy_jpk(indsxy)
!
             CALL mk8vct(lon(1:),gridij(:,:)%longi,jk,jt, &
     &                   jsxy,nbr,kflagxyo)
             CALL mk8vct(lat(1:),gridij(:,:)%latj, jk,jt, &
     &                   jsxy,nbr,kflagxyo)
!
             js1 = js0 + nbr - 1
             CALL disp_ylm( vects(js0:js1), &
     &                      lon(1:nbr), lat(1:nbr), jl, jm )
             js0 = js0 + nbr
          ENDDO
          ENDDO
        ENDDO
      ENDIF
!     
! -3.- Write output data
! ----------------------
!
      IF (jproc.EQ.0) THEN
      IF ((kargtypeoper(1:1).EQ.'-').OR.(kargtypeoper(1:1).EQ.'D') &
     &                              .OR.(kargtypeoper(1:1).EQ.'R') ) THEN
        SELECT CASE (kflagxyo)
        CASE (1)
             CALL writevar (koutxyo,vects(:),jnxyo)
        CASE (2)
             CALL writedta (koutxyo,vects(:))
        CASE DEFAULT
           GOTO 1000
        END SELECT
      ENDIF
      ENDIF
!
! --- deallocation
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(spct)) deallocate(spct)
!
      IF (allocated(gridij)) deallocate(gridij)
      IF (allocated(weightij)) deallocate(weightij)
      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)
      IF (allocated(wei)) deallocate(wei)
!
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
!
      IF (allocated(ran_timeseries)) deallocate(ran_timeseries)
      IF (allocated(time)) deallocate(time)
      IF (allocated(tspct_freq)) deallocate(tspct_freq)
      IF (allocated(tspct_power)) deallocate(tspct_power)
      IF (allocated(spctstd)) deallocate(spctstd)

      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algospct','algospct')
 1001 CALL printerror2(0,1001,3,'algospct','algospct')
!
 101  WRITE (texterror,*) 'Negative observation error std'
      CALL printerror2(0,101,3,'algospct','algospct',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algospctbas (kinxyobas,koutxyobas, &
     &                        kflagxyo,kargtypeoper,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose : Compute transformation of an ensemble of vectors
!
!  Method : Read vector of variables from input file (kinxyo),
!  ------   compute transformation for each 2D slice of variables,
!           write the result in output file (koutxyo).
!
!  Input : kinxyobas   : Cxyo input directory
!  -----   koutxyobas  : Cxyo output directory
!          kflagxyo    : Vector type (1=Vx,2=Vy)
!          kargtypeoper : transformation direction (+ or -) and
!                         maximum degree of spherical harmonics (jpl)
!          kconfigo    : Observation operator (Vo) inputy file
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jprend
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinxyobas,koutxyobas,kargtypeoper
      CHARACTER(len=*), intent(in) :: kconfigo
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jprsize,jprbas,jr,serie
      LOGICAL :: lectinfo
      INTEGER :: jrbasdeb,jrbasfin
      CHARACTER(len=bgword) :: inxyo, outxyo, vctnam
!----------------------------------------------------------------------
!
      jprsize=jprend
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modspct/algospctbas :'
         WRITE(numout,*) '     compute ensemble psectral transformation'
      ENDIF
!
! Check validity of input directory
      IF (.NOT.(validextbas(kinxyobas))) GOTO 102
      SELECT CASE (kflagxyo)
      CASE (1)
         IF (.NOT.(validextvarbas(kinxyobas))) GOTO 102
      CASE (2)
         IF ((.NOT.(validextvarbas(kinxyobas))) &
     &        .AND.(.NOT.(validextdtabas(kinxyobas)))) GOTO 102
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Check validity of input directory
      IF (.NOT.(validextbas(koutxyobas))) GOTO 102
      SELECT CASE (kflagxyo)
      CASE (1)
         IF (.NOT.(validextvarbas(koutxyobas))) GOTO 102
      CASE (2)
         IF ((.NOT.(validextvarbas(koutxyobas))) &
     &        .AND.(.NOT.(validextdtabas(koutxyobas)))) GOTO 102
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -1.- Loop on ensemble files
! ---------------------------
!
      IF (nprint.GE.2) THEN
        WRITE(numout,*) '    ==> Loop on ensemble files'
      ENDIF
!
      serie=1
      jrbasdeb=1
      jrbasfin=jprsize
      lectinfo=.FALSE.
!
      DO jr=jrbasdeb,jrbasfin
! Control print
         IF (nprint.GE.3) THEN
            WRITE(numout,*) '    ==> Transforming vector number ',jr
            WRITE(*,'(a,i4)') '==> Transforming member: ',jr
         ENDIF
! Build input file name for vector number jr
         CALL fildirbas (vctnam,kinxyobas,jprbas,jr,serie)
         WRITE(inxyo,'("./",A,"/",A)') kinxyobas(1:lenv(kinxyobas)), &
     &        vctnam(1:lenv(vctnam))
! Build output file name for vector number jr
         CALL fildirbas (vctnam,koutxyobas,jprbas,jr,serie)
         WRITE(outxyo,'("./",A,"/",A)') koutxyobas(1:lenv(koutxyobas)), &
     &        vctnam(1:lenv(vctnam))
! Apply transformation to this vector
         CALL algospctvct(inxyo,outxyo,kflagxyo,kargtypeoper,kconfigo)
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algospct','algospctbas')
 1001 CALL printerror2(0,1001,3,'algospct','algospctbas')
!
 102  WRITE (texterror,*) 'Invalid ensemble directory name'
      CALL printerror2(0,102,3,'algospct','algospctbas',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE proj_wei (weightij)
!---------------------------------------------------------------------
!
!  Purpose : Compute sphercial surface element for each grid point
!
!  Method : Compute surface of each quadrangle using spherical trigonometry
!
!  Output : weightij : surface associated to each grid point (sphere = 1)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_coord
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      REAL(KIND=8), DIMENSION(:,:), INTENT( out ) :: weightij
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpisize,jpjsize,ji,jj
      BIGREAL :: lat1,lat2,lat3,lat4
      BIGREAL :: lon1,lon2,lon3,lon4
      BIGREAL :: c12, c23, c34, c41, c13, s, e, a123, a143, a
!
      REAL(KIND=8), PARAMETER :: twopi=2*3.1415926535897932384626
      REAL(KIND=8), PARAMETER :: deg2rad=twopi/360.
!----------------------------------------------------------------------
      jpisize=SIZE(weightij,1)
      jpjsize=SIZE(weightij,2)
!
! Make sure that longitude is increasing with ji
!     DO ji=2,jpisize
!     DO jj=1,jpjsize
!       IF (gridij(ji,jj)%longi.LT.gridij(ji-1,jj)%longi) THEN
!         gridij(ji,jj)%longi=gridij(ji,jj)%longi+360.
!       ENDIF
!     ENDDO
!     ENDDO
!
      DO ji=2,jpisize-1
      DO jj=2,jpjsize-1
!       coordinates of the vertices of the quadrangle
!       lon1 = 0.5 * ( gridij(ji+1,jj+1)%longi + gridij(ji,jj)%longi ) * deg2rad
!       lon2 = 0.5 * ( gridij(ji+1,jj-1)%longi + gridij(ji,jj)%longi ) * deg2rad
!       lon3 = 0.5 * ( gridij(ji-1,jj-1)%longi + gridij(ji,jj)%longi ) * deg2rad
!       lon4 = 0.5 * ( gridij(ji-1,jj+1)%longi + gridij(ji,jj)%longi ) * deg2rad
!       lat1 = ( 90. - 0.5 * ( gridij(ji+1,jj+1)%latj + gridij(ji,jj)%latj ) ) * deg2rad
!       lat2 = ( 90. - 0.5 * ( gridij(ji+1,jj-1)%latj + gridij(ji,jj)%latj ) ) * deg2rad
!       lat3 = ( 90. - 0.5 * ( gridij(ji-1,jj-1)%latj + gridij(ji,jj)%latj ) ) * deg2rad
!       lat4 = ( 90. - 0.5 * ( gridij(ji-1,jj+1)%latj + gridij(ji,jj)%latj ) ) * deg2rad
        lon1 = gridij(ji,jj)%longi * deg2rad
        lon2 = gridij(ji+1,jj)%longi * deg2rad
        lon3 = gridij(ji+1,jj+1)%longi * deg2rad
        lon4 = gridij(ji,jj+1)%longi * deg2rad
        lat1 = ( 90. - gridij(ji,jj)%latj ) * deg2rad
        lat2 = ( 90. - gridij(ji+1,jj)%latj ) * deg2rad
        lat3 = ( 90. - gridij(ji+1,jj+1)%latj ) * deg2rad
        lat4 = ( 90. - gridij(ji,jj+1)%latj ) * deg2rad
!       distance between vertices
        c12 = ACOS ( MIN( SIN(lat1)*SIN(lat2)*COS(lon1-lon2) + COS(lat1)*COS(lat2) , 1._kr) )
        c23 = ACOS ( MIN( SIN(lat2)*SIN(lat3)*COS(lon2-lon3) + COS(lat2)*COS(lat3) , 1._kr) )
        c34 = ACOS ( MIN( SIN(lat3)*SIN(lat4)*COS(lon3-lon4) + COS(lat3)*COS(lat4) , 1._kr) )
        c41 = ACOS ( MIN( SIN(lat4)*SIN(lat1)*COS(lon4-lon1) + COS(lat4)*COS(lat1) , 1._kr) )
        c13 = ACOS ( MIN( SIN(lat1)*SIN(lat3)*COS(lon1-lon3) + COS(lat1)*COS(lat3) , 1._kr) )
!       area of 1-2-3 triangle
        s = 0.5 * ( c12 + c23 + c13 )
        e = sqrt(MAX(tan(s/2)*tan((s-c12)/2)*tan((s-c23)/2)*tan((s-c13)/2),0._kr))
        a123 = 4.0 * ATAN(e)
!       area of 1-4-3 triangle
        s = 0.5 * ( c34 + c41 + c13 )
        e = sqrt(MAX(tan(s/2)*tan((s-c34)/2)*tan((s-c41)/2)*tan((s-c13)/2),0._kr))
        a143 = 4.0 * ATAN(e)
!       total area of quadrangle (sphere = 4pi)
        a = a123 + a143
!
        weightij(ji,jj) = 0.5 * a / twopi
      ENDDO
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algospct','proj_wei')
 1001 CALL printerror2(0,1001,3,'algospct','proj_wei')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE prep_obs( kjsxy, kjt, kjk, kvects, &
     &                     kposcoefobs, kgridijkobs, vectorms, &
     &                     kobs, klon, klat, kobswei, knbr )
!---------------------------------------------------------------------
!
!  Purpose : Prepare observations for current 2D slice
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo, only : obsndbs, obs_ind, obs_nbr, var_lev, dta_ord
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, INTENT( in ) :: kjsxy, kjt, kjk
      INTEGER, INTENT( out ) :: knbr
      BIGREAL, DIMENSION(:), INTENT( in ) :: kvects
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
      TYPE (type_gridijk), dimension(:), intent(in)  :: kgridijkobs
      BIGREAL, DIMENSION(:), INTENT( in ) :: vectorms
      BIGREAL, DIMENSION(:), INTENT( out ) :: kobs, kobswei
      BIGREAL, DIMENSION(:), INTENT( out ) :: klon, klat
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jposize, jpitpsize, jdbs, jo, inddta, jpksize
      LOGICAL :: accepted
!----------------------------------------------------------------------
      jposize = SIZE(kvects,1)
      jpitpsize = SIZE(kposcoefobs,2)
      IF (SIZE(kobs,1).NE.jposize) GOTO 1000
      IF (SIZE(kobswei,1).NE.jposize) GOTO 1000
      IF (SIZE(klon,1).NE.jposize) GOTO 1000
      IF (SIZE(klat,1).NE.jposize) GOTO 1000
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 1000
      IF (jposize.NE.size(kgridijkobs,1)) GOTO 1000
      jpksize = SIZE(var_lev,1)
!
! Select observations corresponding to required 2D slice
      knbr=0
      DO jdbs=1,obsndbs(kjsxy) 
        DO jo=obs_ind(kjsxy,jdbs),obs_ind(kjsxy,jdbs)+obs_nbr(kjsxy,jdbs)-1
          accepted = .TRUE.
          IF (jpitpsize.EQ.8) THEN
            inddta= dta_ord(kjsxy)
            IF (kjk.EQ.1) THEN
              accepted = kgridijkobs(jo)%levk.LT.var_lev(2,inddta)
            ELSEIF (kjk.EQ.jpksize) THEN
              accepted = kgridijkobs(jo)%levk.GT.var_lev(jpksize-1,inddta)
            ELSE
              accepted = ( kgridijkobs(jo)%levk.GT.var_lev(kjk-1,inddta) ) &
     &             .AND. ( kgridijkobs(jo)%levk.LT.var_lev(kjk+1,inddta) )
            ENDIF
          ENDIF
          IF (accepted) THEN
            knbr=knbr+1
            kobs(knbr) = kvects(jo)
            klon(knbr) = kgridijkobs(jo)%longi
            klat(knbr) = kgridijkobs(jo)%latj
            kobswei(knbr) = 1.0/vectorms(jo)
          ENDIF
        ENDDO
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algospct','prep_obs')
 1001 CALL printerror2(0,1001,3,'algospct','prep_obs')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE define_spct( spctstd, tspct_freq, tspct_power )
!---------------------------------------------------------------------
!
!  Purpose : Define spectrum for random timeseries and
!            draw random amplitudes for spherical harmonics
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, DIMENSION(0:,-jpl:), INTENT( out ) :: spctstd
      BIGREAL, DIMENSION(:), INTENT( out ) :: tspct_freq
      BIGREAL, DIMENSION(:), INTENT( out ) :: tspct_power
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jm, jl, jtspct
      BIGREAL :: wnbr, norm, freq, maxfreq
!----------------------------------------------------------------------

!     Loop on frequancies
!     to define spectrum in time
      maxfreq=2.0_kr
      DO jtspct=1,jptspct
        freq = maxfreq * REAL(jtspct,8) / REAL(jptspct,8)
        tspct_power(jtspct) = EXP(-freq*freq)
        freq = freq / loc_time_scale
        tspct_freq(jtspct) = freq
      ENDDO

!     Loop on spherical harmonics
!     to define the spectrum sph_spct
      DO jl = 0, jpl
      DO jm = -jl, jl
        wnbr = REAL(jl,8) * 2.0_kr * pi * loc_radius_in_deg / 360._kr
        spctstd(jl,jm) = EXP(-wnbr*wnbr)
      ENDDO
      ENDDO

      wnbr = REAL(jpl,8) * 2.0_kr * pi * loc_radius_in_deg / 360._kr
      PRINT *, 'Localization radius in degrees:',loc_radius_in_deg
      PRINT *, 'Length scale corresponding to lmax:',loc_radius_in_deg/wnbr
      IF (wnbr.LT.3.0_kr) THEN
        PRINT *, 'Warning: increase lmax to resolve the spectrum'
      ENDIF

!     Normalize the spectrum sph_spct
      norm = 0.
      DO jl = 0, jpl
      DO jm = -jl, jl
        norm = norm + spctstd(jl,jm)
      ENDDO
      ENDDO
      spctstd = spctstd / norm

!     Compute square root of the power spectrum
      spctstd = SQRT(spctstd)

      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algospct','prep_obs')
 1001 CALL printerror2(0,1001,3,'algospct','prep_obs')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algospct
