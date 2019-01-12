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
! ---                   LIOCRG.F90                                ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-11 (J.M. Brankart)                      ---
! --- revised      : 00-03 (J.M. Brankart)                      ---
! --- revised      : 01-06 (C.E. Testut)                        ---
! --- revised      : 03-03 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  readcrg        : Read NetCDF gridded observation 
! ---                              data base file
! --- SUBROUTINE  evalhdrcrg     : Read header of NetCDF gridded
! ---                              observation data base file
! --- SUBROUTINE  readcrgbias    : Read NetCDF gridded observation database
! ---                              file (2D constant field to add to observations)
! --- SUBROUTINE  evalhdrcrgsize : Read header of NetCDF gridded
! ---                              observation data base file (2D
! ---                              constant field to add to observations)
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liocrg
      use mod_main
      use utilfiles
      IMPLICIT NONE
      PRIVATE

      PUBLIC readcrg,evalhdrcrg,readcrgbias,evalhdrcrgsize

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcrg (kfnincrg,kvectcrg,crgobsnam, &
     &                    spvalcrgout,kgridij,kgridijk)
!---------------------------------------------------------------------
!
!  Purpose : Read NetCDF gridded observation data base file
!  -------
!  Method : Read data extraction criterion from configuration file
!  ------   Read and select observations from NetCDF database
!
!  Input :  kfnincrg    : filename
!  -----    crgobsnam   : variable name
!           spvalcrgout : special value
!
!  Output : kvectdbs    : observation values
!  ------   kgridij     : observation location (2D case)
!           kgridijk    : observation location (3D case)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilcdfvar
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnincrg
      BIGREAL, dimension(:), intent(out) :: kvectcrg
      CHARACTER(len=*), intent(in) :: crgobsnam
      BIGREAL, intent(out) :: spvalcrgout
      TYPE (type_gridij), optional, dimension (:) :: kgridij
      TYPE (type_gridijk), optional, dimension (:) :: kgridijk
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=word80) :: title,longname,unit,crgobsnamloc
      CHARACTER(len=hgword) :: line
      CHARACTER(len=word80), dimension(4) :: dimunit 
      INTEGER :: allocok,numfila
      INTEGER :: ji,jpi,jj,jpj,jk,jpk,jt,jpt,jpcrgsize,jo,jpo
      INTEGER :: ji_min,jj_min,jk_min,jt_min
      INTEGER :: ji_max,jj_max,jk_max,jt_max
      INTEGER :: jjini,jjend,jjstep
      INTEGER :: crgdim,dtadim
      BIGREAL4 :: spvalcrg,spvalxy
      BIGREAL4, dimension(:,:), allocatable :: tabcrg,lonxy,latxy,depthxy,timexy
      BIGREAL4, dimension(:), allocatable :: lon,lat,depth,time
      BIGREAL4 :: time_min,time_max,lon_min,lon_max
      BIGREAL4 :: lat_min,lat_max,depth_min,depth_max
      BIGREAL4 :: olon,olat,odepth,otime
      LOGICAL :: time_irr,depth_irr,lon_irr,lat_irr
      LOGICAL :: selected,inside,inside1,inside2
      CHARACTER(len=1) :: textexclusion
      CHARACTER(len=10) :: key,value,lon_typ,lat_typ,depth_typ,time_typ
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../mkcrgtoobs/readcrg'
         WRITE(numout,*) '    ==> READING file ',kfnincrg(1:lenv(kfnincrg))
      ENDIF
!
! Check coherence of input arguments and array sizes
      jpcrgsize=size(kvectcrg,1)
      crgobsnamloc=crgobsnam
!
      IF (.NOT.present(kgridij)) THEN
         IF (.NOT.present(kgridijk)) GOTO 1000
      ENDIF
      IF (present(kgridij)) THEN
         IF (present(kgridijk)) GOTO 1000
      ENDIF
      IF (present(kgridij)) THEN
         IF (jpcrgsize.NE.size(kgridij,1)) GOTO 102
      ENDIF
      IF (present(kgridijk)) THEN
         IF (jpcrgsize.NE.size(kgridijk,1)) GOTO 102
      ENDIF
!
      IF (present(kgridij)) THEN
         dtadim = 2
      ELSE
         dtadim = 3
      ENDIF
!
! -1.- Read optional data extraction parameters
! ---------------------------------------------
! Available options:
! dbsobsnam :  name of observed variable in NetCDF database
! lon       :  type of longitude grid (regular, xy, xt, xyt or xyzt)
! lat       :  type of latitude grid (regular, xy, xt, xyt or xyzt)
! depth     :  type of depth grid (regular, xt, xyt or xyzt)
! time      :  type of time grid (regular, xt, xyt or xyzt)
! lon_min   :  minimum observation longitude
! lon_max   :  maximum observation longitude
! lat_min   :  minimum observation latitude
! lat_max   :  maximum observation latitude
! depth_min :  minimum observation depth
! depth_max :  maximum observation depth
! time_min  :  minimum observation time
! time_max  :  maximum observation time
! ji_min    :  minimum 1st dimension grid index
! ji_max    :  maximum 1st dimension grid index
! jj_min    :  minimum 2nd dimension grid index
! jj_max    :  maximum 2nd dimension grid index
! jk_min    :  minimum 3rd dimension grid index
! jk_max    :  maximum 3rd dimension grid index
! jt_min    :  minimum 4th dimension grid index
! jt_max    :  maximum 4th dimension grid index
!
! Other variables
! lon_irr   :  irregular longitude grid (T or F)
! lat_irr   :  irregular latitude grid (T or F)
! depth_irr :  irregular depth grid (T or F)
! time_irr  :  irregular time grid (T or F)
!
      lon_irr=.FALSE.
      lat_irr=.FALSE.
      depth_irr=.FALSE.
      time_irr=.FALSE.
      lon_typ='regular'
      lat_typ='regular'
      depth_typ='regular'
      time_typ='regular'
!
      lon_min=-HUGE(lon_min)
      lon_max=HUGE(lon_max)
      lat_min=-HUGE(lat_min)
      lat_max=HUGE(lat_max)
      depth_min=-HUGE(depth_min)
      depth_max=HUGE(depth_max)
      time_min=-HUGE(time_min)
      time_max=HUGE(time_max)
!
      ji_min=-HUGE(ji_min)
      ji_max=HUGE(ji_max)
      jj_min=-HUGE(jj_min)
      jj_max=HUGE(jj_max)
      jk_min=-HUGE(jk_min)
      jk_max=HUGE(jk_max)
      jt_min=-HUGE(jt_min)
      jt_max=HUGE(jt_max)
!
      IF (larginoptcfg) THEN
        textexclusion='#'
        numfila=99
        CALL openfile(numfila,arginoptcfg)
!
        DO WHILE (line(1:3).NE.'end')
          line=readnextline(numfila,textexclusion)
          READ(line,'(2A)') key,value
!
          SELECT CASE(key(1:lenv(key)))
          CASE('dbsobsnam')
             READ(value,*) crgobsnamloc
          CASE('lon')
             READ(value,*) lon_typ
             SELECT CASE(lon_typ)
             CASE ('regular')
             CASE ('xy','xt','xyt','xyzt')
                lon_irr=.TRUE.
             CASE DEFAULT
                GOTO 104
             END SELECT
          CASE('lat')
             READ(value,*) lat_typ
             SELECT CASE(lat_typ)
             CASE ('regular')
             CASE ('xy','xt','xyt','xyzt')
                lat_irr=.TRUE.
             CASE DEFAULT
                GOTO 104
             END SELECT
          CASE('depth')
             READ(value,*) depth_typ
             SELECT CASE(depth_typ)
             CASE ('regular')
             CASE ('xt','xyt','xyzt')
                depth_irr=.TRUE.
             CASE DEFAULT
                GOTO 104
             END SELECT
          CASE('time')
             READ(value,*) time_typ
             SELECT CASE(time_typ)
             CASE ('regular')
             CASE ('xt','xyt','xyzt')
                time_irr=.TRUE.
             CASE DEFAULT
                GOTO 104
             END SELECT
          CASE('lon_min')
             READ(value,*) lon_min
          CASE('lon_max')
             READ(value,*) lon_max
          CASE('lat_min')
             READ(value,*) lat_min
          CASE('lat_max')
             READ(value,*) lat_max
          CASE('depth_min')
             READ(value,*) depth_min
          CASE('depth_max')
             READ(value,*) depth_max
          CASE('time_min')
             READ(value,*) time_min
          CASE('time_max')
             READ(value,*) time_max
          CASE('ji_min')
             READ(value,*) ji_min
          CASE('ji_max')
             READ(value,*) ji_max
          CASE('jj_min')
             READ(value,*) jj_min
          CASE('jj_max')
             READ(value,*) jj_max
          CASE('jk_min')
             READ(value,*) jk_min
          CASE('jk_max')
             READ(value,*) jk_max
          CASE('jt_min')
             READ(value,*) jt_min
          CASE('jt_max')
             READ(value,*) jt_max
          CASE('biasobsnam','dbs_max','end')
          CASE DEFAULT
             GOTO 103
          END SELECT
!         
        ENDDO
        CLOSE(numfila)
      ENDIF
!
! -2.- Read database file header
! ------------------------------
!
      CALL cdfrdim(kfnincrg,jpi,jpj,jpk,jpt,title)
      CALL cdfrvar(kfnincrg,crgobsnamloc,spvalcrg,unit,longname,crgdim)
      spvalcrgout = FREAL(spvalcrg)
!
! -3.- Allocate necessary arrays as a function of grid type
! ---------------------------------------------------------
!
! --- allocation tabcrg
      allocate ( tabcrg(1:jpi,1:jpj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabcrg(:,:) = spvalcrgout
! --- allocation lon
      allocate ( lon(1:jpi), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lon(:) = spvalcrgout
! --- allocation lat
      allocate ( lat(1:jpj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lat(:) = spvalcrgout
! --- allocation depth
      allocate ( depth(1:jpk), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      depth(:) = spvalcrgout
! --- allocation time
      allocate ( time(1:jpt), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
! --- allocation lonxy
      IF (lon_irr) THEN
         allocate ( lonxy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         lonxy(:,:) = spvalcrgout
      ENDIF
! --- allocation latxy
      IF (lat_irr) THEN
         allocate ( latxy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         latxy(:,:) = spvalcrgout
      ENDIF
! --- allocation depthxy
      IF (depth_irr) THEN
         allocate ( depthxy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         depthxy(:,:) = spvalcrgout
      ENDIF
! --- allocation timexy
      IF (time_irr) THEN
         allocate ( timexy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         timexy(:,:) = spvalcrgout
      ENDIF
!
! -4.- Read time grid
! -------------------
!
      CALL cdfrtim(kfnincrg,time)
!
! -5.- Read observation locations
! -------------------------------
! (for regular and xy grids)
!
      CALL cdfrloc(kfnincrg,lon,lat,depth,dimunit)
!
      IF (  ((lon_irr).AND.(lon_typ.EQ.'xy'))  .OR. &
     &      ((lat_irr).AND.(lat_typ.EQ.'xy'))  )   THEN
         CALL cdfrpos(kfnincrg,lonxy,latxy,spvalxy,unit)
      ENDIF
!
! -6.- Analyse database records, select reguired observations
! -----------------------------------------------------------
!
      IF ((crgdim.EQ.2).OR.(dtadim.EQ.2)) jpk=1
!
! Initialize data and grid arrays with special values
      kvectcrg(:) = spvalcrgout
      IF (dtadim.EQ.2) THEN
         kgridij(:)%longi = spvalcrgout
         kgridij(:)%latj = spvalcrgout
      ELSE
         kgridijk(:)%longi = spvalcrgout
         kgridijk(:)%latj = spvalcrgout
         kgridijk(:)%levk = spvalcrgout
      ENDIF
!
! Loop on database 4th dimension (time in regular grid)
! Exclude outside data
      jpo = 0
      DO jt=1,jpt
        inside1=(jt.GE.jt_min).AND.(jt.LE.jt_max)
        inside2=(time(jt).GE.time_min).AND.(time(jt).LE.time_max)
        selected=inside1.AND.(inside2.OR.time_irr)
        IF (selected) THEN
!
! Read time irregular grid (xy and xyt cases)
        IF (time_irr) THEN
           IF (time_typ.EQ.'xyt') THEN
             CALL cdfrtab(kfnincrg,'timexyt',jt,timexy)
           ENDIF
           IF (time_typ.EQ.'xt') THEN
             CALL cdfrtab(kfnincrg,'timext',jt,timexy(1:jpi,1))
             DO jj=2,jpj
               timexy(1:jpi,jj)=timexy(1:jpi,1)
             ENDDO
           ENDIF
        ENDIF
!
! Read depth irregular grid (xy and xyt cases)
        IF (depth_irr) THEN
           IF (depth_typ.EQ.'xyt') THEN
             CALL cdfrtab(kfnincrg,'depthxyt',jt,depthxy)
           ENDIF
           IF (depth_typ.EQ.'xt') THEN
             CALL cdfrtab(kfnincrg,'depthxt',jt,depthxy(1:jpi,1))
             DO jj=2,jpj
               depthxy(1:jpi,jj)=depthxy(1:jpi,1)
             ENDDO
           ENDIF
        ENDIF
!
! Read longitude irregular grid (xy and xyt cases)
        IF (lon_irr) THEN
           IF (lon_typ.EQ.'xyt') THEN
             CALL cdfrtab(kfnincrg,'lonxyt',jt,lonxy)
           ENDIF
           IF (lon_typ.EQ.'xt') THEN
             CALL cdfrtab(kfnincrg,'lonxt',jt,lonxy(1:jpi,1))
             DO jj=2,jpj
               lonxy(1:jpi,jj)=lonxy(1:jpi,1)
             ENDDO
           ENDIF
        ENDIF
!
! Read latitude irregular grid (xy and xyt cases)
        IF (lat_irr) THEN
           IF (lat_typ.EQ.'xyt') THEN
             CALL cdfrtab(kfnincrg,'latxyt',jt,latxy)
           ENDIF
           IF (lat_typ.EQ.'xt') THEN
             CALL cdfrtab(kfnincrg,'latxt',jt,latxy(1:jpi,1))
             DO jj=2,jpj
               latxy(1:jpi,jj)=latxy(1:jpi,1)
             ENDDO
           ENDIF
        ENDIF
!
! Loop on database 3rd dimension (depth in regular grid)
! Exclude outside data
        DO jk=1,jpk
          inside1=(jk.GE.jk_min).AND.(jk.LE.jk_max)
          inside2=(depth(jk).GE.depth_min).AND.(depth(jk).LE.depth_max)
          selected=inside1.AND.(inside2.OR.depth_irr)
          IF (selected) THEN
!
! Read time irregular grid (xyzt cases)
          IF ((time_irr).AND.(time_typ.EQ.'xyzt')) THEN
             CALL cdfrsli(kfnincrg,'timexyzt',3,jk,jt,timexy)
          ENDIF
!
! Read depth irregular grid (xyzt cases)
          IF ((depth_irr).AND.(depth_typ.EQ.'xyzt')) THEN
             CALL cdfrsli(kfnincrg,'depthxyzt',3,jk,jt,depthxy)
          ENDIF
!
! Read longitude irregular grid (xyzt cases)
          IF ((lon_irr).AND.(lon_typ.EQ.'xyzt')) THEN
             CALL cdfrsli(kfnincrg,'lonxyzt',3,jk,jt,lonxy)
          ENDIF
!
! Read latitude irregular grid (xyzt cases)
          IF ((lat_irr).AND.(lat_typ.EQ.'xyzt')) THEN
             CALL cdfrsli(kfnincrg,'latxyzt',3,jk,jt,latxy)
          ENDIF
!
! Read slice of data (along 1st and 2nd database dimensions)
          IF (crgdim.EQ.2) THEN
             CALL cdfrtab(kfnincrg,crgobsnamloc,jt,tabcrg)
          ELSE
             CALL cdfrsli(kfnincrg,crgobsnamloc,3,jk,jt,tabcrg)
          ENDIF
!
! Loop on database 1st dimension (longitude in regular grid)
! Exclude outside data
          DO ji = 1,jpi
            inside1=(ji.GE.ji_min).AND.(ji.LE.ji_max)
            inside2=(lon(ji).GE.lon_min).AND.(lon(ji).LE.lon_max)
            selected=inside1.AND.(inside2.OR.lon_irr)
            IF (selected) THEN
!
            IF (MOD(ji,2).EQ.0) THEN
              jjini = jpj
              jjend = 1
              jjstep = -1
            ELSE
              jjini = 1
              jjend = jpj
              jjstep = 1
            ENDIF
!
! Loop on database 2nd dimension (latitude in regular grid)
! Exclude outside data
            DO jj = jjini,jjend,jjstep
              inside1=(jj.GE.jj_min).AND.(jj.LE.jj_max)
              inside2=(lat(jj).GE.lat_min).AND.(lat(jj).LE.lat_max)
              selected=inside1.AND.(inside2.OR.lat_irr)
!
! Select inside data
              inside=.TRUE.
              IF (time_irr) THEN
                 otime=timexy(ji,jj)
                 inside=inside.AND.(otime.GE.time_min)
                 inside=inside.AND.(otime.LE.time_max)
              ELSE
                 otime=time(jt)
              ENDIF
              IF (depth_irr) THEN
                 odepth=depthxy(ji,jj)
                 inside=inside.AND.(odepth.GE.depth_min)
                 inside=inside.AND.(odepth.LE.depth_max)
              ELSE
                 odepth=depth(jk)
              ENDIF
              IF (lon_irr) THEN
                 olon=lonxy(ji,jj)
                 inside=inside.AND.(olon.GE.lon_min)
                 inside=inside.AND.(olon.LE.lon_max)
              ELSE
                 olat=lat(jj)
              ENDIF
              IF (lat_irr) THEN
                 olat=latxy(ji,jj)
                 inside=inside.AND.(olat.GE.lat_min)
                 inside=inside.AND.(olat.LE.lat_max)
              ELSE
                 olon=lon(ji)
              ENDIF
              selected=selected.AND.inside
!
              IF (selected) THEN
!
! Exclude special value
                IF (tabcrg(ji,jj).NE.spvalcrg) THEN
! Add observation in data vector
                  jpo = jpo + 1
                  IF (jpo.GT.jpcrgsize) GOTO 105
                  kvectcrg(jpo) = FREAL(tabcrg(ji,jj))
                  IF (dtadim.EQ.2) THEN
                     kgridij(jpo)%longi = FREAL(olon)
                     kgridij(jpo)%latj = FREAL(olat)
                  ELSE
                     kgridijk(jpo)%longi = FREAL(olon)
                     kgridijk(jpo)%latj = FREAL(olat)
                     kgridijk(jpo)%levk = FREAL(odepth)
                  ENDIF
                ENDIF
!
              ENDIF 
            ENDDO
            ENDIF 
          ENDDO
          ENDIF 
        ENDDO
        ENDIF
      ENDDO
!
! -7.- Control print
! ------------------
!
      IF (nprint.GE.4) THEN
         WRITE(numout,*) ' number of observation extracted: ',jpo
         WRITE(numout,*) ' sample of observation extracted from the database:'
         WRITE(numout,*) ' idx lon lat (depth) value'
         DO jo=1,jpo,MAX(1,jpo/25)
            IF (dtadim.EQ.2) THEN
               WRITE(numout,*) jo, kgridij(jo)%longi, &
     &            kgridij(jo)%latj, kvectcrg(jo)
            ELSE
               WRITE(numout,*) jo, kgridijk(jo)%longi, &
     &            kgridijk(jo)%latj, kgridijk(jo)%levk, &
     &            kvectcrg(jo)
            ENDIF
         ENDDO
         WRITE(numout,*) ' '
      ENDIF
!
! --- deallocate
      IF (allocated(tabcrg)) deallocate(tabcrg)
      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)
      IF (allocated(depth)) deallocate(depth)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liocrg','readcrg')
 1001 CALL printerror2(0,1001,3,'liocrg','readcrg')
!
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,3,'liocrg','readcrg',comment=texterror)
 103  WRITE (texterror,*) 'Invalid keyword in -inoptcfg file'
      CALL printerror2(0,103,3,'liocrg','readcrg',comment=texterror)
 104  WRITE (texterror,*) 'Bad grid type in -inoptcfg file'
      CALL printerror2(0,104,3,'liocrg','readcrg',comment=texterror)
 105  WRITE (texterror,*) 'Insufficient database vector size in SESAM'
      CALL printerror2(0,105,1,'liocrg','readcrg',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcrg (kfnincrg,jpcrgout)
!---------------------------------------------------------------------
!
!  Purpose : Read header of NetCDF gridded
!  -------   observation data base file
!
!  Method : Read data extraction criterion from configuration file
!  ------   Evaluate maximum number of extracted data
!
!  Input :  kfnincrg : database filename
!  -----
!  Output : jpcrgout : database size
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilcdfvar
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnincrg
      INTEGER, intent(out) :: jpcrgout
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=word80) :: title,kform
      CHARACTER(len=hgword) :: line
      INTEGER :: imax,jmax,kmax,tmax,jt,dbs_max
      INTEGER :: ji_min,jj_min,jk_min,jt_min
      INTEGER :: ji_max,jj_max,jk_max,jt_max
      BIGREAL4 :: time_min,time_max,lon_min,lon_max
      BIGREAL4 :: lat_min,lat_max,depth_min,depth_max
      BIGREAL4, allocatable, dimension(:) :: time
      INTEGER :: numfila,allocok
      LOGICAL :: time_irr,depth_irr,lon_irr,lat_irr
      CHARACTER(len=1) :: textexclusion
      CHARACTER(len=10) :: key,value,lon_typ,lat_typ,depth_typ,time_typ
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../evalhdrdbs/evalhdrcrg'
         WRITE(numout,*) '    ==> READING file ',kfnincrg(1:lenv(kfnincrg))
      ENDIF
!
! -1.- Read header of .crg file
! -----------------------------
!
      CALL cdfrdim(kfnincrg,imax,jmax,kmax,tmax,title)
!
! -2.- Read optional data extraction parameters
!      in order to evaluate the size of the database
! --------------------------------------------------
! See decription of the options in readcrg subroutine
!
      ji_min=1
      ji_max=imax
      jj_min=1
      jj_max=jmax
      jk_min=1
      jk_max=kmax
      jt_min=1
      jt_max=tmax
!
      dbs_max=0
!
      lon_min=-HUGE(lon_min)
      lon_max=HUGE(lon_max)
      lat_min=-HUGE(lat_min)
      lat_max=HUGE(lat_max)
      depth_min=-HUGE(depth_min)
      depth_max=HUGE(depth_max)
      time_min=-HUGE(time_min)
      time_max=HUGE(time_max)
!
      lon_typ='regular'
      lat_typ='regular'
      depth_typ='regular'
      time_typ='regular'
!
      IF (larginoptcfg) THEN
        textexclusion='#'
        numfila=99
        CALL openfile(numfila,arginoptcfg)
!
        DO WHILE (line(1:3).NE.'end')
          line=readnextline(numfila,textexclusion)
          READ(line,'(2A)') key,value
!
          SELECT CASE(key(1:lenv(key)))
          CASE('lon')
             READ(value,*) lon_typ
             SELECT CASE(lon_typ)
             CASE ('regular')
             CASE ('xy','xt','xyt','xyzt')
                lon_irr=.TRUE.
             CASE DEFAULT
                GOTO 104
             END SELECT
          CASE('lat')
             READ(value,*) lat_typ
             SELECT CASE(lat_typ)
             CASE ('regular')
             CASE ('xy','xt','xyt','xyzt')
                lat_irr=.TRUE.
             CASE DEFAULT
                GOTO 104
             END SELECT
          CASE('depth')
             READ(value,*) depth_typ
             SELECT CASE(depth_typ)
             CASE ('regular')
             CASE ('xt','xyt','xyzt')
                depth_irr=.TRUE.
             CASE DEFAULT
                GOTO 104
             END SELECT
          CASE('time')
             READ(value,*) time_typ
             SELECT CASE(time_typ)
             CASE ('regular')
             CASE ('xt','xyt','xyzt')
                time_irr=.TRUE.
             CASE DEFAULT
                GOTO 104
             END SELECT
          CASE('lon_min')
             READ(value,*) lon_min
          CASE('lon_max')
             READ(value,*) lon_max
          CASE('lat_min')
             READ(value,*) lat_min
          CASE('lat_max')
             READ(value,*) lat_max
          CASE('depth_min')
             READ(value,*) depth_min
          CASE('depth_max')
             READ(value,*) depth_max
          CASE('time_min')
             READ(value,*) time_min
          CASE('time_max')
             READ(value,*) time_max
          CASE('ji_min')
             READ(value,*) ji_min
          CASE('ji_max')
             READ(value,*) ji_max
          CASE('jj_min')
             READ(value,*) jj_min
          CASE('jj_max')
             READ(value,*) jj_max
          CASE('jk_min')
             READ(value,*) jk_min
          CASE('jk_max')
             READ(value,*) jk_max
          CASE('jt_min')
             READ(value,*) jt_min
          CASE('jt_max')
             READ(value,*) jt_max
          CASE('dbs_max')
             READ(value,*) dbs_max
          CASE('biasobsnam','dbsobsnam')
          CASE('end')
          CASE DEFAULT
             GOTO 103
          END SELECT
!
        ENDDO
        CLOSE(numfila)
      ENDIF
!
! -2.- Control print
! ------------------
!
      IF (nprint.GE.2) THEN
         kform='(8x,2a)'
         WRITE(numout,kform) '- Database title: ',title(1:lenv(title))
         kform='(8x,a,3i5)'
         WRITE(numout,kform) '- isiz,imin,imax: ',imax,ji_min,ji_max
         WRITE(numout,kform) '- jsiz,jmin,jmax: ',jmax,jj_min,jj_max
         WRITE(numout,kform) '- ksiz,kmin,kmax: ',kmax,jk_min,jk_max
         WRITE(numout,kform) '- tsiz,tmin,tmax: ',tmax,jt_min,jt_max
         kform='(8x,a,e12.3)'
         WRITE(numout,kform) '- lon_min:',lon_min
         WRITE(numout,kform) '- lon_max:',lon_max
         WRITE(numout,kform) '- lat_min:',lat_min
         WRITE(numout,kform) '- lat_max:',lat_max
         WRITE(numout,kform) '- depth_min:',depth_min
         WRITE(numout,kform) '- depth_max:',depth_max
         WRITE(numout,kform) '- time_min:',time_min
         WRITE(numout,kform) '- time_max:',time_max
      ENDIF
!
! -3.- Compute maximum number of data that will be extracted from database
! ------------------------------------------------------------------------
! 
      imax=ji_max-ji_min+1
      jmax=jj_max-jj_min+1
      kmax=jk_max-jk_min+1
      IF ( (.NOT.lon_irr).AND.(lon_min.EQ.lon_max) ) imax = 1
      IF ( (.NOT.lat_irr).AND.(lat_min.EQ.lat_max) ) jmax = 1
      IF ( (.NOT.depth_irr).AND.(depth_min.EQ.depth_max) ) kmax = 1
      IF (.NOT.time_irr) THEN
        IF (time_min.EQ.time_max) THEN
           tmax = 1
        ELSE
!
           allocate ( time(1:tmax), stat=allocok )
           IF (allocok.NE.0) GOTO 1001
!
           CALL cdfrtim(kfnincrg,time)
!
           DO jt=jt_min,jt_max
              IF (time(jt).LE.time_min) THEN
                 jt_min=jt
              ENDIF
              IF (time(jt).GE.time_max) THEN
                 jt_max=jt
                 exit
              ENDIF
           ENDDO
!
           tmax=jt_max-jt_min+1
!
           IF (allocated(time)) deallocate(time)
!
        ENDIF
      ELSE
        tmax=jt_max-jt_min+1
      ENDIF
!
      IF (dbs_max.GT.0) THEN
         jpcrgout = dbs_max
      ELSE
         jpcrgout = imax * jmax * kmax * tmax
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liocrg','evalhdrcrg')
 1001 CALL printerror2(0,1001,3,'liocrg','evalhdrcrg')
!
 103  WRITE (texterror,*) 'Bad key word in -inoptcfg file'
      CALL printerror2(0,103,3,'liocrg','evalhdrcrg',comment=texterror)
 104  WRITE (texterror,*) 'Bad grid type in -inoptcfg file'
      CALL printerror2(0,104,3,'liocrg','evalhdrcrg',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcrgbias (kfnincrg,crgobsnam,kmatcrg,kloncrg, &
     &                        klatcrg,spvalcrgout)
!---------------------------------------------------------------------
!
!  Purpose : Read NetCDF gridded observation database file
!  -------   (read 2D constant field to add to observations)
!
!  Method : Read data extraction criterion from configuration file
!  ------   Read and select observations from NetCDF database
!
!  Input :  kfnincrg    : filename
!  -----    crgobsnam   : variable name
!
!  Output : kmatcrg     : 2D matrix with horizontal field
!  ------   kloncrg     : position in longitude
!           klatcrg     : position in latitude
!           spvalcrgout : special value
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilcdfvar
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnincrg
      CHARACTER(len=*), intent(in) :: crgobsnam
      BIGREAL, dimension(:,:), intent(out) :: kmatcrg
      BIGREAL, dimension(:), intent(out) :: kloncrg
      BIGREAL, dimension(:), intent(out) :: klatcrg
      BIGREAL, intent(out) :: spvalcrgout
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=word80) :: title,longname,unit,crgobsnamloc
      CHARACTER(len=hgword) :: line
      CHARACTER(len=word80), dimension(4) :: dimunit
      INTEGER :: jpi,jpj,jpk,jpt,jpicrg,jpjcrg,jpkcrg
      INTEGER :: crgidx,ndim
      BIGREAL4 :: spvalcrg
      BIGREAL4, dimension(:,:), allocatable :: tabcrg
      BIGREAL4, dimension(:), allocatable :: lon,lat,depth
      INTEGER :: numfila,allocok
      CHARACTER(len=1) :: textexclusion
      CHARACTER(len=10) :: key,value
!---------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : /sesam/modobsv/mkdbstoobs/mkcrgtoobs/readcrgbias'
         WRITE(numout,*) '    ==> READING file ',kfnincrg(1:lenv(kfnincrg))
      ENDIF
!
! Check coherence of input arguments and array sizes
      jpicrg=size(kmatcrg,1)
      jpjcrg=size(kmatcrg,2)
      crgobsnamloc=crgobsnam
!
      IF (jpicrg.NE.size(kloncrg,1)) GOTO 102
      IF (jpjcrg.NE.size(klatcrg,1)) GOTO 102
!
! -1.- Read optional data extraction parameters
! ---------------------------------------------
! See option description in readcrg subroutine
      IF (larginoptcfg) THEN
        textexclusion='#'
        numfila=99
        CALL openfile(numfila,arginoptcfg)
!
        DO WHILE (line(1:3).NE.'end')
          line=readnextline(numfila,textexclusion)
          READ(line,'(2A)') key,value
!
          SELECT CASE(key(1:lenv(key)))
          CASE('biasobsnam')
             READ(value,*) crgobsnamloc
          CASE('lon','lat','depth','time')
          CASE('lon_min','lon_max','lat_min','lat_max')
          CASE('depth_min','depth_max','time_min','time_max')
          CASE('ji_min','ji_max','jj_min','jj_max')
          CASE('jk_min','jk_max','jt_min','jt_max')
          CASE('dbsobsnam','dbs_max')
          CASE('end')
          CASE DEFAULT
             GOTO 103
          END SELECT
!
        ENDDO
        CLOSE(numfila)
      ENDIF
!
! -2.- Read database file header
! ------------------------------
!
      CALL cdfrdim(kfnincrg,jpi,jpj,jpk,jpt,title)
      IF (jpi.NE.jpicrg) GOTO 1000
      IF (jpj.NE.jpjcrg) GOTO 1000
      jpkcrg = jpk
      CALL cdfrvar(kfnincrg,crgobsnamloc,spvalcrg,unit,longname,ndim)
      spvalcrgout = FREAL(spvalcrg)
!
!
! -3.- Allocate necessary arrays (grid is regular)
! ------------------------------------------------
!
! --- allocation tabcrg
      allocate ( tabcrg(1:jpicrg,1:jpjcrg), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabcrg(:,:) = spvalcrgout
! --- allocation lon
      allocate ( lon(1:jpicrg), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lon(:) = spvalcrgout
! --- allocation lat
      allocate ( lat(1:jpjcrg), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lat(:) = spvalcrgout
! --- allocation depth
      allocate ( depth(1:jpkcrg), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      depth(:) = spvalcrgout
!
! -5.- Read 2D observation locations
! ----------------------------------
!
      CALL cdfrloc(kfnincrg,lon,lat,depth,dimunit)
!
      kloncrg(:) = FREAL(lon(:))
      klatcrg(:) = FREAL(lat(:))
!
! -6.- Read 2D observation matrix
! -------------------------------
!
      crgidx = 1
      CALL cdfrsli(kfnincrg,crgobsnamloc,3,1,crgidx,tabcrg)
!
      kmatcrg(:,:) = FREAL(tabcrg(:,:))
!
! --- deallocate
      IF (allocated(tabcrg)) deallocate(tabcrg)
      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)
      IF (allocated(depth)) deallocate(depth)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liocrg','readbias')
 1001 CALL printerror2(0,1001,3,'liocrg','readbias')
!
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,3,'liocrg','readcrgbias',comment=texterror)
 103  WRITE (texterror,*) 'Bad key word in -inoptcfg file'
      CALL printerror2(0,103,3,'liocrg','readcrg',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcrgsize (kfnincrg,jpicrg,jpjcrg,jpkcrg)
!---------------------------------------------------------------------
!
!  Purpose : Read header of NetCDF gridded observation data base file
!  -------   (2D constant field to add to observations)
!
!  Method : Read data extraction criterion from configuration file
!  ------   Evaluate maximum number of extracted data
!
!  Input :  kfnincrg : database filename
!  -----
!  Output : jpicrg   : size of data base 1st dimension
!  ------   jpjcrg   : size of data base 2nd dimension
!           jpkcrg   : size of data base 3rd dimension
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilcdfvar
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnincrg
      INTEGER, intent(out) :: jpicrg,jpjcrg,jpkcrg
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=word80) :: title
      INTEGER :: imax,jmax,kmax,tmax
!----------------------------------------------------------------------
!
! Read database file header
      CALL cdfrdim(kfnincrg,imax,jmax,kmax,tmax,title)
!
! Output size of database array
      jpicrg = imax
      jpjcrg = jmax
      jpkcrg = kmax
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liocrg','evalhdrcrgsize')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liocrg
