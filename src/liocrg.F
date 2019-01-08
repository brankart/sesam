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
C ---                   LIOCRG.F                                ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-11 (J.M. Brankart)                      ---
C --- revised      : 00-03 (J.M. Brankart)                      ---
C --- revised      : 01-06 (C.E. Testut)                        ---
C --- revised      : 03-03 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  readcrg        : Read NetCDF gridded observation 
C ---                              data base file
C --- SUBROUTINE  evalhdrcrg     : Read header of NetCDF gridded
C ---                              observation data base file
C --- SUBROUTINE  readcrgbias    : Read NetCDF gridded observation database
C ---                              file (2D constant field to add to observations)
C --- SUBROUTINE  evalhdrcrgsize : Read header of NetCDF gridded
C ---                              observation data base file (2D
C ---                              constant field to add to observations)
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liocrg
      use mod_main
      use utilfiles
      IMPLICIT NONE
      PRIVATE

      PUBLIC readcrg,evalhdrcrg,readcrgbias,evalhdrcrgsize

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcrg (kfnincrg,kvectcrg,crgobsnam,
     $                    spvalcrgout,kgridij,kgridijk)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read NetCDF gridded observation data base file
CCC  -------
CCC  Method : Read data extraction criterion from configuration file
CCC  ------   Read and select observations from NetCDF database
CCC
CCC  Input :  kfnincrg    : filename
CCC  -----    crgobsnam   : variable name
CCC           spvalcrgout : special value
CCC
CCC  Output : kvectdbs    : observation values
CCC  ------   kgridij     : observation location (2D case)
CCC           kgridijk    : observation location (3D case)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilcdfvar
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrg
      BIGREAL, dimension(:), intent(out) :: kvectcrg
      CHARACTER(len=*), intent(in) :: crgobsnam
      BIGREAL, intent(out) :: spvalcrgout
      TYPE (type_gridij), optional, dimension (:) :: kgridij
      TYPE (type_gridijk), optional, dimension (:) :: kgridijk
CC----------------------------------------------------------------------
CC local declarations
CC ==================
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
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../mkcrgtoobs/readcrg'
         WRITE(numout,*) '    ==> READING file ',kfnincrg(1:lenv(kfnincrg))
      ENDIF
C
C Check coherence of input arguments and array sizes
      jpcrgsize=size(kvectcrg,1)
      crgobsnamloc=crgobsnam
C
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
C
      IF (present(kgridij)) THEN
         dtadim = 2
      ELSE
         dtadim = 3
      ENDIF
C
C -1.- Read optional data extraction parameters
C ---------------------------------------------
C Available options:
C dbsobsnam :  name of observed variable in NetCDF database
C lon       :  type of longitude grid (regular, xy, xt, xyt or xyzt)
C lat       :  type of latitude grid (regular, xy, xt, xyt or xyzt)
C depth     :  type of depth grid (regular, xt, xyt or xyzt)
C time      :  type of time grid (regular, xt, xyt or xyzt)
C lon_min   :  minimum observation longitude
C lon_max   :  maximum observation longitude
C lat_min   :  minimum observation latitude
C lat_max   :  maximum observation latitude
C depth_min :  minimum observation depth
C depth_max :  maximum observation depth
C time_min  :  minimum observation time
C time_max  :  maximum observation time
C ji_min    :  minimum 1st dimension grid index
C ji_max    :  maximum 1st dimension grid index
C jj_min    :  minimum 2nd dimension grid index
C jj_max    :  maximum 2nd dimension grid index
C jk_min    :  minimum 3rd dimension grid index
C jk_max    :  maximum 3rd dimension grid index
C jt_min    :  minimum 4th dimension grid index
C jt_max    :  maximum 4th dimension grid index
C
C Other variables
C lon_irr   :  irregular longitude grid (T or F)
C lat_irr   :  irregular latitude grid (T or F)
C depth_irr :  irregular depth grid (T or F)
C time_irr  :  irregular time grid (T or F)
C
      lon_irr=.FALSE.
      lat_irr=.FALSE.
      depth_irr=.FALSE.
      time_irr=.FALSE.
      lon_typ='regular'
      lat_typ='regular'
      depth_typ='regular'
      time_typ='regular'
C
      lon_min=-HUGE(lon_min)
      lon_max=HUGE(lon_max)
      lat_min=-HUGE(lat_min)
      lat_max=HUGE(lat_max)
      depth_min=-HUGE(depth_min)
      depth_max=HUGE(depth_max)
      time_min=-HUGE(time_min)
      time_max=HUGE(time_max)
C
      ji_min=-HUGE(ji_min)
      ji_max=HUGE(ji_max)
      jj_min=-HUGE(jj_min)
      jj_max=HUGE(jj_max)
      jk_min=-HUGE(jk_min)
      jk_max=HUGE(jk_max)
      jt_min=-HUGE(jt_min)
      jt_max=HUGE(jt_max)
C
      IF (larginoptcfg) THEN
        textexclusion='#'
        numfila=99
        CALL openfile(numfila,arginoptcfg)
C
        DO WHILE (line(1:3).NE.'end')
          line=readnextline(numfila,textexclusion)
          READ(line,'(2A)') key,value
C
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
C         
        ENDDO
        CLOSE(numfila)
      ENDIF
C
C -2.- Read database file header
C ------------------------------
C
      CALL cdfrdim(kfnincrg,jpi,jpj,jpk,jpt,title)
      CALL cdfrvar(kfnincrg,crgobsnamloc,spvalcrg,unit,longname,crgdim)
      spvalcrgout = FREAL(spvalcrg)
C
C -3.- Allocate necessary arrays as a function of grid type
C ---------------------------------------------------------
C
C --- allocation tabcrg
      allocate ( tabcrg(1:jpi,1:jpj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabcrg(:,:) = spvalcrgout
C --- allocation lon
      allocate ( lon(1:jpi), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lon(:) = spvalcrgout
C --- allocation lat
      allocate ( lat(1:jpj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lat(:) = spvalcrgout
C --- allocation depth
      allocate ( depth(1:jpk), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      depth(:) = spvalcrgout
C --- allocation time
      allocate ( time(1:jpt), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C --- allocation lonxy
      IF (lon_irr) THEN
         allocate ( lonxy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         lonxy(:,:) = spvalcrgout
      ENDIF
C --- allocation latxy
      IF (lat_irr) THEN
         allocate ( latxy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         latxy(:,:) = spvalcrgout
      ENDIF
C --- allocation depthxy
      IF (depth_irr) THEN
         allocate ( depthxy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         depthxy(:,:) = spvalcrgout
      ENDIF
C --- allocation timexy
      IF (time_irr) THEN
         allocate ( timexy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         timexy(:,:) = spvalcrgout
      ENDIF
C
C -4.- Read time grid
C -------------------
C
      CALL cdfrtim(kfnincrg,time)
C
C -5.- Read observation locations
C -------------------------------
C (for regular and xy grids)
C
      CALL cdfrloc(kfnincrg,lon,lat,depth,dimunit)
C
      IF (  ((lon_irr).AND.(lon_typ.EQ.'xy'))  .OR.
     $      ((lat_irr).AND.(lat_typ.EQ.'xy'))  )   THEN
         CALL cdfrpos(kfnincrg,lonxy,latxy,spvalxy,unit)
      ENDIF
C
C -6.- Analyse database records, select reguired observations
C -----------------------------------------------------------
C
      IF ((crgdim.EQ.2).OR.(dtadim.EQ.2)) jpk=1
C
C Initialize data and grid arrays with special values
      kvectcrg(:) = spvalcrgout
      IF (dtadim.EQ.2) THEN
         kgridij(:)%longi = spvalcrgout
         kgridij(:)%latj = spvalcrgout
      ELSE
         kgridijk(:)%longi = spvalcrgout
         kgridijk(:)%latj = spvalcrgout
         kgridijk(:)%levk = spvalcrgout
      ENDIF
C
C Loop on database 4th dimension (time in regular grid)
C Exclude outside data
      jpo = 0
      DO jt=1,jpt
        inside1=(jt.GE.jt_min).AND.(jt.LE.jt_max)
        inside2=(time(jt).GE.time_min).AND.(time(jt).LE.time_max)
        selected=inside1.AND.(inside2.OR.time_irr)
        IF (selected) THEN
C
C Read time irregular grid (xy and xyt cases)
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
C
C Read depth irregular grid (xy and xyt cases)
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
C
C Read longitude irregular grid (xy and xyt cases)
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
C
C Read latitude irregular grid (xy and xyt cases)
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
C
C Loop on database 3rd dimension (depth in regular grid)
C Exclude outside data
        DO jk=1,jpk
          inside1=(jk.GE.jk_min).AND.(jk.LE.jk_max)
          inside2=(depth(jk).GE.depth_min).AND.(depth(jk).LE.depth_max)
          selected=inside1.AND.(inside2.OR.depth_irr)
          IF (selected) THEN
C
C Read time irregular grid (xyzt cases)
          IF ((time_irr).AND.(time_typ.EQ.'xyzt')) THEN
             CALL cdfrsli(kfnincrg,'timexyzt',3,jk,jt,timexy)
          ENDIF
C
C Read depth irregular grid (xyzt cases)
          IF ((depth_irr).AND.(depth_typ.EQ.'xyzt')) THEN
             CALL cdfrsli(kfnincrg,'depthxyzt',3,jk,jt,depthxy)
          ENDIF
C
C Read longitude irregular grid (xyzt cases)
          IF ((lon_irr).AND.(lon_typ.EQ.'xyzt')) THEN
             CALL cdfrsli(kfnincrg,'lonxyzt',3,jk,jt,lonxy)
          ENDIF
C
C Read latitude irregular grid (xyzt cases)
          IF ((lat_irr).AND.(lat_typ.EQ.'xyzt')) THEN
             CALL cdfrsli(kfnincrg,'latxyzt',3,jk,jt,latxy)
          ENDIF
C
C Read slice of data (along 1st and 2nd database dimensions)
          IF (crgdim.EQ.2) THEN
             CALL cdfrtab(kfnincrg,crgobsnamloc,jt,tabcrg)
          ELSE
             CALL cdfrsli(kfnincrg,crgobsnamloc,3,jk,jt,tabcrg)
          ENDIF
C
C Loop on database 1st dimension (longitude in regular grid)
C Exclude outside data
          DO ji = 1,jpi
            inside1=(ji.GE.ji_min).AND.(ji.LE.ji_max)
            inside2=(lon(ji).GE.lon_min).AND.(lon(ji).LE.lon_max)
            selected=inside1.AND.(inside2.OR.lon_irr)
            IF (selected) THEN
C
            IF (MOD(ji,2).EQ.0) THEN
              jjini = jpj
              jjend = 1
              jjstep = -1
            ELSE
              jjini = 1
              jjend = jpj
              jjstep = 1
            ENDIF
C
C Loop on database 2nd dimension (latitude in regular grid)
C Exclude outside data
            DO jj = jjini,jjend,jjstep
              inside1=(jj.GE.jj_min).AND.(jj.LE.jj_max)
              inside2=(lat(jj).GE.lat_min).AND.(lat(jj).LE.lat_max)
              selected=inside1.AND.(inside2.OR.lat_irr)
C
C Select inside data
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
C
              IF (selected) THEN
C
C Exclude special value
                IF (tabcrg(ji,jj).NE.spvalcrg) THEN
C Add observation in data vector
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
C
              ENDIF 
            ENDDO
            ENDIF 
          ENDDO
          ENDIF 
        ENDDO
        ENDIF
      ENDDO
C
C -7.- Control print
C ------------------
C
      IF (nprint.GE.4) THEN
         WRITE(numout,*) ' number of observation extracted: ',jpo
         WRITE(numout,*) ' sample of observation extracted from the database:'
         WRITE(numout,*) ' idx lon lat (depth) value'
         DO jo=1,jpo,MAX(1,jpo/25)
            IF (dtadim.EQ.2) THEN
               WRITE(numout,*) jo, kgridij(jo)%longi,
     $            kgridij(jo)%latj, kvectcrg(jo)
            ELSE
               WRITE(numout,*) jo, kgridijk(jo)%longi,
     $            kgridijk(jo)%latj, kgridijk(jo)%levk,
     $            kvectcrg(jo)
            ENDIF
         ENDDO
         WRITE(numout,*) ' '
      ENDIF
C
C --- deallocate
      IF (allocated(tabcrg)) deallocate(tabcrg)
      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)
      IF (allocated(depth)) deallocate(depth)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrg','readcrg')
 1001 CALL printerror2(0,1001,3,'liocrg','readcrg')
C
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,3,'liocrg','readcrg',comment=texterror)
 103  WRITE (texterror,*) 'Invalid keyword in -inoptcfg file'
      CALL printerror2(0,103,3,'liocrg','readcrg',comment=texterror)
 104  WRITE (texterror,*) 'Bad grid type in -inoptcfg file'
      CALL printerror2(0,104,3,'liocrg','readcrg',comment=texterror)
 105  WRITE (texterror,*) 'Insufficient database vector size in SESAM'
      CALL printerror2(0,105,1,'liocrg','readcrg',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcrg (kfnincrg,jpcrgout)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read header of NetCDF gridded
CCC  -------   observation data base file
CCC
CCC  Method : Read data extraction criterion from configuration file
CCC  ------   Evaluate maximum number of extracted data
CCC
CCC  Input :  kfnincrg : database filename
CCC  -----
CCC  Output : jpcrgout : database size
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilcdfvar
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrg
      INTEGER, intent(out) :: jpcrgout
CC----------------------------------------------------------------------
CC local declarations
CC ==================
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
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../evalhdrdbs/evalhdrcrg'
         WRITE(numout,*) '    ==> READING file ',kfnincrg(1:lenv(kfnincrg))
      ENDIF
C
C -1.- Read header of .crg file
C -----------------------------
C
      CALL cdfrdim(kfnincrg,imax,jmax,kmax,tmax,title)
C
C -2.- Read optional data extraction parameters
C      in order to evaluate the size of the database
C --------------------------------------------------
C See decription of the options in readcrg subroutine
C
      ji_min=1
      ji_max=imax
      jj_min=1
      jj_max=jmax
      jk_min=1
      jk_max=kmax
      jt_min=1
      jt_max=tmax
C
      dbs_max=0
C
      lon_min=-HUGE(lon_min)
      lon_max=HUGE(lon_max)
      lat_min=-HUGE(lat_min)
      lat_max=HUGE(lat_max)
      depth_min=-HUGE(depth_min)
      depth_max=HUGE(depth_max)
      time_min=-HUGE(time_min)
      time_max=HUGE(time_max)
C
      lon_typ='regular'
      lat_typ='regular'
      depth_typ='regular'
      time_typ='regular'
C
      IF (larginoptcfg) THEN
        textexclusion='#'
        numfila=99
        CALL openfile(numfila,arginoptcfg)
C
        DO WHILE (line(1:3).NE.'end')
          line=readnextline(numfila,textexclusion)
          READ(line,'(2A)') key,value
C
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
C
        ENDDO
        CLOSE(numfila)
      ENDIF
C
C -2.- Control print
C ------------------
C
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
C
C -3.- Compute maximum number of data that will be extracted from database
C ------------------------------------------------------------------------
C 
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
C
           allocate ( time(1:tmax), stat=allocok )
           IF (allocok.NE.0) GOTO 1001
C
           CALL cdfrtim(kfnincrg,time)
C
           DO jt=jt_min,jt_max
              IF (time(jt).LE.time_min) THEN
                 jt_min=jt
              ENDIF
              IF (time(jt).GE.time_max) THEN
                 jt_max=jt
                 exit
              ENDIF
           ENDDO
C
           tmax=jt_max-jt_min+1
C
           IF (allocated(time)) deallocate(time)
C
        ENDIF
      ELSE
        tmax=jt_max-jt_min+1
      ENDIF
C
      IF (dbs_max.GT.0) THEN
         jpcrgout = dbs_max
      ELSE
         jpcrgout = imax * jmax * kmax * tmax
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrg','evalhdrcrg')
 1001 CALL printerror2(0,1001,3,'liocrg','evalhdrcrg')
C
 103  WRITE (texterror,*) 'Bad key word in -inoptcfg file'
      CALL printerror2(0,103,3,'liocrg','evalhdrcrg',comment=texterror)
 104  WRITE (texterror,*) 'Bad grid type in -inoptcfg file'
      CALL printerror2(0,104,3,'liocrg','evalhdrcrg',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcrgbias (kfnincrg,crgobsnam,kmatcrg,kloncrg,
     $                        klatcrg,spvalcrgout)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read NetCDF gridded observation database file
CCC  -------   (read 2D constant field to add to observations)
CCC
CCC  Method : Read data extraction criterion from configuration file
CCC  ------   Read and select observations from NetCDF database
CCC
CCC  Input :  kfnincrg    : filename
CCC  -----    crgobsnam   : variable name
CCC
CCC  Output : kmatcrg     : 2D matrix with horizontal field
CCC  ------   kloncrg     : position in longitude
CCC           klatcrg     : position in latitude
CCC           spvalcrgout : special value
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilcdfvar
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrg
      CHARACTER(len=*), intent(in) :: crgobsnam
      BIGREAL, dimension(:,:), intent(out) :: kmatcrg
      BIGREAL, dimension(:), intent(out) :: kloncrg
      BIGREAL, dimension(:), intent(out) :: klatcrg
      BIGREAL, intent(out) :: spvalcrgout
CC----------------------------------------------------------------------
CC local declarations
CC ==================
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
CCC---------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : /sesam/modobsv/mkdbstoobs/mkcrgtoobs/readcrgbias'
         WRITE(numout,*) '    ==> READING file ',kfnincrg(1:lenv(kfnincrg))
      ENDIF
C
C Check coherence of input arguments and array sizes
      jpicrg=size(kmatcrg,1)
      jpjcrg=size(kmatcrg,2)
      crgobsnamloc=crgobsnam
C
      IF (jpicrg.NE.size(kloncrg,1)) GOTO 102
      IF (jpjcrg.NE.size(klatcrg,1)) GOTO 102
C
C -1.- Read optional data extraction parameters
C ---------------------------------------------
C See option description in readcrg subroutine
      IF (larginoptcfg) THEN
        textexclusion='#'
        numfila=99
        CALL openfile(numfila,arginoptcfg)
C
        DO WHILE (line(1:3).NE.'end')
          line=readnextline(numfila,textexclusion)
          READ(line,'(2A)') key,value
C
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
C
        ENDDO
        CLOSE(numfila)
      ENDIF
C
C -2.- Read database file header
C ------------------------------
C
      CALL cdfrdim(kfnincrg,jpi,jpj,jpk,jpt,title)
      IF (jpi.NE.jpicrg) GOTO 1000
      IF (jpj.NE.jpjcrg) GOTO 1000
      jpkcrg = jpk
      CALL cdfrvar(kfnincrg,crgobsnamloc,spvalcrg,unit,longname,ndim)
      spvalcrgout = FREAL(spvalcrg)
C
C
C -3.- Allocate necessary arrays (grid is regular)
C ------------------------------------------------
C
C --- allocation tabcrg
      allocate ( tabcrg(1:jpicrg,1:jpjcrg), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabcrg(:,:) = spvalcrgout
C --- allocation lon
      allocate ( lon(1:jpicrg), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lon(:) = spvalcrgout
C --- allocation lat
      allocate ( lat(1:jpjcrg), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lat(:) = spvalcrgout
C --- allocation depth
      allocate ( depth(1:jpkcrg), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      depth(:) = spvalcrgout
C
C -5.- Read 2D observation locations
C ----------------------------------
C
      CALL cdfrloc(kfnincrg,lon,lat,depth,dimunit)
C
      kloncrg(:) = FREAL(lon(:))
      klatcrg(:) = FREAL(lat(:))
C
C -6.- Read 2D observation matrix
C -------------------------------
C
      crgidx = 1
      CALL cdfrsli(kfnincrg,crgobsnamloc,3,1,crgidx,tabcrg)
C
      kmatcrg(:,:) = FREAL(tabcrg(:,:))
C
C --- deallocate
      IF (allocated(tabcrg)) deallocate(tabcrg)
      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)
      IF (allocated(depth)) deallocate(depth)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrg','readbias')
 1001 CALL printerror2(0,1001,3,'liocrg','readbias')
C
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,3,'liocrg','readcrgbias',comment=texterror)
 103  WRITE (texterror,*) 'Bad key word in -inoptcfg file'
      CALL printerror2(0,103,3,'liocrg','readcrg',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcrgsize (kfnincrg,jpicrg,jpjcrg,jpkcrg)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read header of NetCDF gridded observation data base file
CCC  -------   (2D constant field to add to observations)
CCC
CCC  Method : Read data extraction criterion from configuration file
CCC  ------   Evaluate maximum number of extracted data
CCC
CCC  Input :  kfnincrg : database filename
CCC  -----
CCC  Output : jpicrg   : size of data base 1st dimension
CCC  ------   jpjcrg   : size of data base 2nd dimension
CCC           jpkcrg   : size of data base 3rd dimension
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilcdfvar
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnincrg
      INTEGER, intent(out) :: jpicrg,jpjcrg,jpkcrg
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80) :: title
      INTEGER :: imax,jmax,kmax,tmax
CC----------------------------------------------------------------------
C
C Read database file header
      CALL cdfrdim(kfnincrg,imax,jmax,kmax,tmax,title)
C
C Output size of database array
      jpicrg = imax
      jpjcrg = jmax
      jpkcrg = kmax
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'liocrg','evalhdrcrgsize')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liocrg
