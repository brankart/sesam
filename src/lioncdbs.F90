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
! ---                   LIONCDBS.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2010-04 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  readncdbs    : Read NetCDF gridded observation 
! ---                            data base file
! --- SUBROUTINE  evalhdrncdbs : Read header of NetCDF gridded
! ---                            observation data base file
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE lioncdbs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC readncdbs,evalhdrncdbs

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readncdbs (kfnindbs,kvectdbs,kjobs, &
     &                      spvaldbsout,kgridij,kgridijk,kgrid)
!---------------------------------------------------------------------
!
!  Purpose : Read NetCDF gridded observation data base file
!  -------
!  Method : Read data extraction criterion from configuration file
!  ------   Read and select observations from NetCDF database
!
!  Input :  kfnindbs    : filename
!  -----    kjobs       : observation index
!           spvaldbsout : special value
!
!  Output : kvectdbs    : observation values
!  ------   kgridij     : observation location (2D case)
!           kgridijk    : observation location (3D case)
!           kgrid       : observation location (4D case)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnindbs
      BIGREAL, dimension(:), intent(out) :: kvectdbs
      INTEGER, intent(in) :: kjobs
      BIGREAL, intent(out) :: spvaldbsout
      TYPE (type_gridij), optional, dimension (:) :: kgridij
      TYPE (type_gridijk), optional, dimension (:) :: kgridijk
      TYPE (type_grid4d), optional, dimension (:) :: kgrid
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=word80) :: dbsobsnam,kform
      INTEGER :: allocok,indobs,inddbs
      INTEGER :: ji,jpi,jj,jpj,jk,jpk,jt,jpt,jpdbssize,jo,jpo
      INTEGER :: ji_min,jj_min,jk_min,jt_min
      INTEGER :: ji_max,jj_max,jk_max,jt_max
      INTEGER :: jjini,jjend,jjstep
      INTEGER :: dtadim
      BIGREAL4, dimension(:,:), allocatable :: tabdbs,lonxy,latxy, &
     &                                         depthxy,timexy
      BIGREAL4, dimension(:), allocatable :: lon,lat,depth,time
      BIGREAL4 :: dbs_min,dbs_max,tim_min,tim_max,lon_min,lon_max
      BIGREAL4 :: lat_min,lat_max,dep_min,dep_max,dbsexcl
      BIGREAL4 :: lon_def,lat_def,dep_def,tim_def
      BIGREAL4 :: olon,olat,odep,otim
      LOGICAL :: tim_irr,dep_irr,lon_irr,lat_irr
      LOGICAL :: tim_var,dep_var,lon_var,lat_var
      LOGICAL :: selected,inside
!
      INTEGER, allocatable, dimension(:) :: idims,vstart,vcount,vstrid
      INTEGER, allocatable, dimension(:) :: xstart,xcount,xstrid
      INTEGER, allocatable, dimension(:) :: ystart,ycount,ystrid
      INTEGER, allocatable, dimension(:) :: zstart,zcount,zstrid
      INTEGER, allocatable, dimension(:) :: tstart,tcount,tstrid
      INTEGER :: ierr, idf, idx, idy, idz, idt, idd
      INTEGER :: idv, idvx, idvy, idvz, idvt, ndims
      INTEGER :: idtim, idtimx, idtimy, idtimz, idtimt
      INTEGER :: iddep, iddepx, iddepy, iddepz, iddept
      INTEGER :: idlat, idlatx, idlaty, idlatz, idlatt
      INTEGER :: idlon, idlonx, idlony, idlonz, idlont
      CHARACTER(len=bgword) :: xdim, ydim, zdim, tdim
      CHARACTER(len=bgword) :: xnam, ynam, znam, tnam
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../mkncdbstoobs/readncdbs'
         WRITE(numout,*) '    ==> READING file ',kfnindbs(1:lenv(kfnindbs))
      ENDIF
!
! Get database array size
      jpdbssize=size(kvectdbs,1)
!
! Get observation and database index
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
!
! Check grid type
      IF (present(kgridij)) THEN
        IF (present(kgridijk)) GOTO 1000
        IF (present(kgrid)) GOTO 1000
        IF (jpdbssize.NE.size(kgridij,1)) GOTO 102
        dtadim = 2
      ENDIF
!
      IF (present(kgridijk)) THEN
        IF (present(kgridij)) GOTO 1000
        IF (present(kgrid)) GOTO 1000
        IF (jpdbssize.NE.size(kgridijk,1)) GOTO 102
        dtadim = 3
      ENDIF
!
      IF (present(kgrid)) THEN
        IF (present(kgridij)) GOTO 1000
        IF (present(kgridijk)) GOTO 1000
        IF (jpdbssize.NE.size(kgrid,1)) GOTO 102
        dtadim = 4
      ENDIF
!
! -1.- Get data extraction parameters
! -----------------------------------
! dbsobsnam :  name of observed variable in NetCDF database
! dbsexcl   :  observation exclusion value
! dbs_min   :  minimum observation value
! dbs_max   :  maximum observation value
! lon_def   :  default longitude
! lon_min   :  minimum observation longitude
! lon_max   :  maximum observation longitude
! lat_def   :  default latitude
! lat_min   :  minimum observation latitude
! lat_max   :  maximum observation latitude
! dep_def   :  default depth
! dep_min   :  minimum observation depth
! dep_max   :  maximum observation depth
! tim_def   :  default time
! tim_min   :  minimum observation time
! tim_max   :  maximum observation time
! ji_min    :  minimum 1st dimension grid index
! ji_max    :  maximum 1st dimension grid index
! jj_min    :  minimum 2nd dimension grid index
! jj_max    :  maximum 2nd dimension grid index
! jk_min    :  minimum 3rd dimension grid index
! jk_max    :  maximum 3rd dimension grid index
! jt_min    :  minimum 4th dimension grid index
! jt_max    :  maximum 4th dimension grid index
!
      dbsobsnam=obsifil(indobs,inddbs)
      xdim=obsxdim(indobs,inddbs)
      xnam=obsxnam(indobs,inddbs) ; lon_var=xnam.NE.'none'
      ydim=obsydim(indobs,inddbs)
      ynam=obsynam(indobs,inddbs) ; lat_var=ynam.NE.'none'
      zdim=obszdim(indobs,inddbs)
      znam=obsznam(indobs,inddbs) ; dep_var=znam.NE.'none'
      tdim=obstdim(indobs,inddbs)
      tnam=obstnam(indobs,inddbs) ; tim_var=tnam.NE.'none'
!
      lon_def=obs_lon_def(indobs,inddbs)
      lat_def=obs_lat_def(indobs,inddbs)
      dep_def=obs_dep_def(indobs,inddbs)
      tim_def=obs_tim_def(indobs,inddbs)
!
      dbsexcl=obsexcl(indobs,inddbs)
      spvaldbsout = FREAL(dbsexcl)
      dbs_min=obs_min(indobs,inddbs)
      dbs_max=obs_max(indobs,inddbs)
      lon_min=obs_lon_min(indobs,inddbs)
      lon_max=obs_lon_max(indobs,inddbs)
      lat_min=obs_lat_min(indobs,inddbs)
      lat_max=obs_lat_max(indobs,inddbs)
      dep_min=obs_dep_min(indobs,inddbs)
      dep_max=obs_dep_max(indobs,inddbs)
      tim_min=obs_tim_min(indobs,inddbs)
      tim_max=obs_tim_max(indobs,inddbs)
!
      ji_min=obsimin(indobs,inddbs)
      ji_max=obsimax(indobs,inddbs)
      jj_min=obsjmin(indobs,inddbs)
      jj_max=obsjmax(indobs,inddbs)
      jk_min=obskmin(indobs,inddbs)
      jk_max=obskmax(indobs,inddbs)
      jt_min=obstmin(indobs,inddbs)
      jt_max=obstmax(indobs,inddbs)
!         
! -2.- Read database file header
! ------------------------------
!
      ierr = NF90_OPEN(kfnindbs,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 101
!
      IF (xdim.NE.'none') THEN
        ierr = NF90_INQ_DIMID(idf,xdim,idx)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=jpi)
        IF (ierr.NE.0) GOTO 103
      ELSE
        idx=-HUGE(idx) ; jpi=1
      ENDIF
!
      IF (ydim.NE.'none') THEN
        ierr = NF90_INQ_DIMID(idf,ydim,idy)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=jpj)
        IF (ierr.NE.0) GOTO 103
      ELSE
        idy=-HUGE(idy) ; jpj=1
      ENDIF
!
      IF (zdim.NE.'none') THEN
        ierr = NF90_INQ_DIMID(idf,zdim,idz)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=jpk)
        IF (ierr.NE.0) GOTO 103
      ELSE
        idz=-HUGE(idz) ; jpk=1
      ENDIF
!
      IF (tdim.NE.'none') THEN
        ierr = NF90_INQ_DIMID(idf,tdim,idt)
        IF (ierr.NE.0) GOTO 103
        ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=jpt)
        IF (ierr.NE.0) GOTO 103
      ELSE
        idt=-HUGE(idt) ; jpt=1
      ENDIF
!
! -3.- Get variable id, get the set of dimensions spanned by the variable,
!      Get indices of x, y, z, t dimensions in that set of dimensions
!      Initialize start, count and strid array for reading Netcdf variable
! ------------------------------------------------------------------------
!
      ierr = NF90_INQ_VARID(idf,dbsobsnam,idv)
      IF (ierr.NE.0) GOTO 104
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
      IF (ierr.NE.0) GOTO 104
!
      allocate(idims(ndims),vstart(ndims),vcount(ndims),vstrid(ndims))
      vstart(1:ndims)=1 ; vcount(1:ndims)=1 ; vstrid(1:ndims)=1
! 
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
      IF (ierr.NE.0) GOTO 104
!
      idvx = 0 ; idvy = 0 ; idvz = 0 ; idvt = 0
      DO idd = 1,ndims
        IF (idims(idd).EQ.idx) idvx = idd
        IF (idims(idd).EQ.idy) idvy = idd
        IF (idims(idd).EQ.idz) idvz = idd
        IF (idims(idd).EQ.idt) idvt = idd
      ENDDO
!
      IF (idvx.NE.0) vcount(idvx)=jpi
      IF (idvy.NE.0) vcount(idvy)=jpj
!
      IF (allocated(idims)) deallocate (idims)
!
!
! -4.- Get coordinate variable ids, get the set of dimensions spanned by the variable,
!      Initialize start, count and strid array for reading coordinate variables
! ------------------------------------------------------------------------------------
! --> Time
      IF (tim_var) THEN
        ierr = NF90_INQ_VARID(idf,tnam,idtim)
        IF (ierr.NE.0) GOTO 105
        ierr = NF90_INQUIRE_VARIABLE(idf,idtim,ndims=ndims)
        IF (ierr.NE.0) GOTO 105
!
        allocate(idims(ndims),tstart(ndims),tcount(ndims),tstrid(ndims))
        tstart(1:ndims)=1 ; tcount(1:ndims)=1 ; tstrid(1:ndims)=1
!
        ierr = NF90_INQUIRE_VARIABLE(idf,idtim,dimids=idims)
        IF (ierr.NE.0) GOTO 105
!
        idtimx = 0 ; idtimy = 0 ; idtimz = 0 ; idtimt = 0
        DO idd = 1,ndims
          IF (idims(idd).EQ.idx) idtimx = idd
          IF (idims(idd).EQ.idy) idtimy = idd
          IF (idims(idd).EQ.idz) idtimz = idd
          IF (idims(idd).EQ.idt) idtimt = idd
        ENDDO
!
        IF (allocated(idims)) deallocate (idims)
!
        tim_irr=(idtimx.NE.0).OR.(idtimy.NE.0).OR.(idtimz.NE.0)
!
        IF (tim_irr) THEN
          IF (idtimx.NE.0) tcount(idtimx)=jpi
          IF (idtimy.NE.0) tcount(idtimy)=jpj
        ELSE
          IF (idtimt.NE.0) tcount(idtimt)=jpt
        ENDIF
      ELSE
        tim_irr=.FALSE.
        IF (idvt.NE.0) GOTO 107
      ENDIF
!
! --> Depth
      IF (dep_var) THEN
        ierr = NF90_INQ_VARID(idf,znam,iddep)
        IF (ierr.NE.0) GOTO 105
        ierr = NF90_INQUIRE_VARIABLE(idf,iddep,ndims=ndims)
        IF (ierr.NE.0) GOTO 105
!
        allocate(idims(ndims),zstart(ndims),zcount(ndims),zstrid(ndims))
        zstart(1:ndims)=1 ; zcount(1:ndims)=1 ; zstrid(1:ndims)=1
!
        ierr = NF90_INQUIRE_VARIABLE(idf,iddep,dimids=idims)
        IF (ierr.NE.0) GOTO 105
!
        iddepx = 0 ; iddepy = 0 ; iddepz = 0 ; iddept = 0
        DO idd = 1,ndims
          IF (idims(idd).EQ.idx) iddepx = idd
          IF (idims(idd).EQ.idy) iddepy = idd
          IF (idims(idd).EQ.idz) iddepz = idd
          IF (idims(idd).EQ.idt) iddept = idd
        ENDDO
!
        IF (allocated(idims)) deallocate (idims)
!
        dep_irr=(iddepx.NE.0).OR.(iddepy.NE.0).OR.(iddept.NE.0)
        print *, 'dep_irr',dep_irr
!
        IF (dep_irr) THEN
          IF (iddepx.NE.0) zcount(iddepx)=jpi
          IF (iddepy.NE.0) zcount(iddepy)=jpj
        ELSE
          IF (iddepz.NE.0) zcount(iddepz)=jpk
        ENDIF
      ELSE
        dep_irr=.FALSE.
        IF (idvz.NE.0) GOTO 107
      ENDIF
!
! --> Latitude
      IF (lat_var) THEN
        ierr = NF90_INQ_VARID(idf,ynam,idlat)
        IF (ierr.NE.0) GOTO 105
        ierr = NF90_INQUIRE_VARIABLE(idf,idlat,ndims=ndims)
        IF (ierr.NE.0) GOTO 105
!
        allocate(idims(ndims),ystart(ndims),ycount(ndims),ystrid(ndims))
        ystart(1:ndims)=1 ; ycount(1:ndims)=1 ; ystrid(1:ndims)=1
!
        ierr = NF90_INQUIRE_VARIABLE(idf,idlat,dimids=idims)
        IF (ierr.NE.0) GOTO 105
!
        idlatx = 0 ; idlaty = 0 ; idlatz = 0 ; idlatt = 0
        DO idd = 1,ndims
          IF (idims(idd).EQ.idx) idlatx = idd
          IF (idims(idd).EQ.idy) idlaty = idd
          IF (idims(idd).EQ.idz) idlatz = idd
          IF (idims(idd).EQ.idt) idlatt = idd
        ENDDO
!
        IF (allocated(idims)) deallocate (idims)
!
        lat_irr=(idlatx.NE.0).OR.(idlatz.NE.0).OR.(idlatt.NE.0)
!
        IF (lat_irr) THEN
          IF (idlatx.NE.0) ycount(idlatx)=jpi
          IF (idlaty.NE.0) ycount(idlaty)=jpj
        ELSE
          IF (idlaty.NE.0) ycount(idlaty)=jpj
        ENDIF
      ELSE
        lat_irr=.FALSE.
        IF (idvy.NE.0) GOTO 107
      ENDIF
!
! --> Longitude
      IF (lon_var) THEN
        ierr = NF90_INQ_VARID(idf,xnam,idlon)
        IF (ierr.NE.0) GOTO 105
        ierr = NF90_INQUIRE_VARIABLE(idf,idlon,ndims=ndims)
        IF (ierr.NE.0) GOTO 105
!
        allocate(idims(ndims),xstart(ndims),xcount(ndims),xstrid(ndims))
        xstart(1:ndims)=1 ; xcount(1:ndims)=1 ; xstrid(1:ndims)=1
!
        ierr = NF90_INQUIRE_VARIABLE(idf,idlon,dimids=idims)
        IF (ierr.NE.0) GOTO 105
!
        idlonx = 0 ; idlony = 0 ; idlonz = 0 ; idlont = 0
        DO idd = 1,ndims
          IF (idims(idd).EQ.idx) idlonx = idd
          IF (idims(idd).EQ.idy) idlony = idd
          IF (idims(idd).EQ.idz) idlonz = idd
          IF (idims(idd).EQ.idt) idlont = idd
        ENDDO
!
        IF (allocated(idims)) deallocate (idims)
!
        lon_irr=(idlony.NE.0).OR.(idlonz.NE.0).OR.(idlont.NE.0)
!
        IF (lon_irr) THEN
          IF (idlonx.NE.0) xcount(idlonx)=jpi
          IF (idlony.NE.0) xcount(idlony)=jpj
        ELSE
          IF (idlonx.NE.0) xcount(idlonx)=jpi
        ENDIF
      ELSE
        lon_irr=.FALSE.
        IF (idvx.NE.0) GOTO 107
      ENDIF
!
! -5.- Allocate necessary arrays as a function of grid type
! ---------------------------------------------------------
!
! --- allocation tabdbs
      allocate ( tabdbs(1:jpi,1:jpj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabdbs(:,:) = spvaldbsout
! --- allocation longitude array
      IF (lon_irr) THEN
         allocate ( lonxy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         lonxy(:,:) = spvaldbsout
      ELSE
         allocate ( lon(1:jpi), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         lon(:) = spvaldbsout
      ENDIF
! --- allocation latitude array
      IF (lat_irr) THEN
         allocate ( latxy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         latxy(:,:) = spvaldbsout
      ELSE
         allocate ( lat(1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         lat(:) = spvaldbsout
      ENDIF
! --- allocation depth array
      IF (dep_irr) THEN
         allocate ( depthxy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         depthxy(:,:) = spvaldbsout
      ELSE
         allocate ( depth(1:jpk), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         depth(:) = spvaldbsout
      ENDIF
! --- allocation time array
      IF (tim_irr) THEN
         allocate ( timexy(1:jpi,1:jpj), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         timexy(:,:) = spvaldbsout
      ELSE
         allocate ( time(1:jpt), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         time(:) = spvaldbsout
      ENDIF
!
! -6.- Read time and depth (regular grids or no variable in file)
! ---------------------------------------------------------------
!
      IF (tim_var) THEN
        IF (.NOT.tim_irr) THEN
! grids: t
          ierr = NF90_GET_VAR(idf,idtim,time,start=tstart, &
     &                        count=tcount,stride=tstrid)
          IF (ierr.NE.0) GOTO 105
        ENDIF
      ELSE
! grids: no variable in file => default constant
        time(:) = tim_def
      ENDIF
!
      IF (dep_var) THEN
        IF (.NOT.dep_irr) THEN
! grids: z
          ierr = NF90_GET_VAR(idf,iddep,depth,start=zstart, &
     &                        count=zcount,stride=zstrid)
          IF (ierr.NE.0) GOTO 105
        ENDIF
      ELSE
! grids: no variable in file => default constant
        depth(:) = dep_def
      ENDIF
!
! -7.- Read longitudes and latitudes
!      (regular grid, xy grids or no variable in file)
! ----------------------------------------------------
!
      IF (lon_var) THEN
        IF (lon_irr) THEN
          IF ((idlonz.EQ.0).AND.(idlont.EQ.0)) THEN
! grids: xy
            ierr = NF90_GET_VAR(idf,idlon,lonxy,start=xstart, &
     &                          count=xcount,stride=xstrid)
            IF (ierr.NE.0) GOTO 105
          ENDIF
        ELSE
! grids: x
          ierr = NF90_GET_VAR(idf,idlon,lon,start=xstart, &
     &                        count=xcount,stride=xstrid)
          IF (ierr.NE.0) GOTO 105
        ENDIF
      ELSE
! grids: no variable in file => default constant
        lon(:) = lon_def
      ENDIF
!
      IF (lat_var) THEN
        IF (lat_irr) THEN
          IF ((idlatz.EQ.0).AND.(idlatt.EQ.0)) THEN
! grids: xy
            ierr = NF90_GET_VAR(idf,idlat,latxy,start=ystart, &
     &                          count=ycount,stride=ystrid)
            IF (ierr.NE.0) GOTO 105
          ENDIF
        ELSE
! grids: y
          ierr = NF90_GET_VAR(idf,idlat,lat,start=ystart, &
     &                        count=ycount,stride=ystrid)
          IF (ierr.NE.0) GOTO 105
        ENDIF
      ELSE
! grids: no variable in file => default constant
        lat(:) = lat_def
      ENDIF
!
! -8.- Analyse database records, select reguired observations
! -----------------------------------------------------------
! Initialize data and grid arrays with special values
      kvectdbs(:) = spvaldbsout
!
      SELECT CASE(dtadim)
      CASE(2)
         kgridij(:)%longi = spvaldbsout
         kgridij(:)%latj = spvaldbsout
      CASE(3)
         kgridijk(:)%longi = spvaldbsout
         kgridijk(:)%latj = spvaldbsout
         kgridijk(:)%levk = spvaldbsout
      CASE(4)
         kgrid(:)%lon = spvaldbsout
         kgrid(:)%lat = spvaldbsout
         kgrid(:)%dep = spvaldbsout
         kgrid(:)%tim = spvaldbsout
      END SELECT
!
! Loop on database 4th dimension (time in regular grid)
! Exclude outside data
      jpo = 0
      DO jt=1,jpt
        selected=(jt.GE.jt_min).AND.(jt.LE.jt_max)
        IF (.NOT.tim_irr) THEN
          selected=selected.AND.(time(jt).GE.tim_min)
          selected=selected.AND.(time(jt).LE.tim_max)
        ENDIF
        IF (selected) THEN
!
! Read time irregular grid
        IF ((tim_irr).AND.(idtimz.EQ.0)) THEN
          IF (idtimt.NE.0) tstart(idtimt)=jt
!
          IF (idtimx.EQ.0) THEN
! grids: y, yt
            ierr = NF90_GET_VAR(idf,idtim,timexy(1,1:jpj),start=tstart, &
     &                          count=tcount,stride=tstrid)
            IF (ierr.NE.0) GOTO 105
            DO ji=2,jpi
              timexy(ji,1:jpj)=timexy(1,1:jpj)
            ENDDO
          ELSEIF (idtimy.EQ.0) THEN
! grids: x, xt
            ierr = NF90_GET_VAR(idf,idtim,timexy(1:jpi,1),start=tstart, &
     &                          count=tcount,stride=tstrid)
            IF (ierr.NE.0) GOTO 105
            DO jj=2,jpj
              timexy(1:jpi,jj)=timexy(1:jpi,1)
            ENDDO
          ELSE
! grids: xy, xyt
            ierr = NF90_GET_VAR(idf,idtim,timexy,start=tstart, &
     &                          count=tcount,stride=tstrid)
            IF (ierr.NE.0) GOTO 105
          ENDIF
!
        ENDIF
!
! Read depth irregular grid
        IF ((dep_irr).AND.(iddepz.EQ.0)) THEN
          IF (iddept.NE.0) zstart(iddept)=jt
!
          IF (iddepx.EQ.0) THEN
! grids: y, yt
            ierr = NF90_GET_VAR(idf,iddep,depthxy(1,1:jpj),start=zstart, &
     &                          count=zcount,stride=zstrid)
            IF (ierr.NE.0) GOTO 105
            DO ji=2,jpi
              depthxy(ji,1:jpj)=depthxy(1,1:jpj)
            ENDDO
          ELSEIF (iddepy.EQ.0) THEN
! grids: x, xt
            ierr = NF90_GET_VAR(idf,iddep,depthxy(1:jpi,1),start=zstart, &
     &                          count=zcount,stride=zstrid)
            IF (ierr.NE.0) GOTO 105
            DO jj=2,jpj
              depthxy(1:jpi,jj)=depthxy(1:jpi,1)
            ENDDO
          ELSE
! grids: xy, xyt
            ierr = NF90_GET_VAR(idf,iddep,depthxy,start=zstart, &
     &                          count=zcount,stride=zstrid)
            IF (ierr.NE.0) GOTO 105
          ENDIF
!
        ENDIF
!
! Read latitude irregular grid
        IF ((lat_irr).AND.(idlatz.EQ.0)) THEN
          IF (idlatt.NE.0) ystart(idlatt)=jt
!
          IF (idlatx.EQ.0) THEN
! grids: yt
            ierr = NF90_GET_VAR(idf,idlat,latxy(1,1:jpj),start=ystart, &
     &                          count=ycount,stride=ystrid)
            IF (ierr.NE.0) GOTO 105
            DO ji=2,jpi
              latxy(ji,1:jpj)=latxy(1,1:jpj)
            ENDDO
          ELSEIF (idlaty.EQ.0) THEN
! grids: xt
            ierr = NF90_GET_VAR(idf,idlat,latxy(1:jpi,1),start=ystart, &
     &                          count=ycount,stride=ystrid)
            IF (ierr.NE.0) GOTO 105
            DO jj=2,jpj
              latxy(1:jpi,jj)=latxy(1:jpi,1)
            ENDDO
          ELSE
! grids: xyt
            ierr = NF90_GET_VAR(idf,idlat,latxy,start=ystart, &
     &                          count=ycount,stride=ystrid)
            IF (ierr.NE.0) GOTO 105
          ENDIF
!
        ENDIF
!
! Read longitude irregular grid
        IF ((lon_irr).AND.(idlonz.EQ.0)) THEN
          IF (idlont.NE.0) xstart(idlont)=jt
!
          IF (idlonx.EQ.0) THEN
! grids: yt
            ierr = NF90_GET_VAR(idf,idlon,lonxy(1,1:jpj),start=xstart, &
     &                          count=xcount,stride=xstrid)
            IF (ierr.NE.0) GOTO 105
            DO ji=2,jpi
              lonxy(ji,1:jpj)=lonxy(1,1:jpj)
            ENDDO
          ELSEIF (idlony.EQ.0) THEN
! grids: xt
            ierr = NF90_GET_VAR(idf,idlon,lonxy(1:jpi,1),start=xstart, &
     &                          count=xcount,stride=xstrid)
            IF (ierr.NE.0) GOTO 105
            DO jj=2,jpj
              lonxy(1:jpi,jj)=lonxy(1:jpi,1)
            ENDDO
          ELSE
! grids: xyt
            ierr = NF90_GET_VAR(idf,idlon,lonxy,start=xstart, &
     &                          count=xcount,stride=xstrid)
            IF (ierr.NE.0) GOTO 105
          ENDIF
!
        ENDIF
!
! Loop on database 3rd dimension (depth in regular grid)
! Exclude outside data
        IF (dtadim.EQ.2) jpk=1
!
        DO jk=1,jpk
          selected=(jk.GE.jk_min).AND.(jk.LE.jk_max)
          IF (.NOT.dep_irr) THEN
            selected=selected.AND.(depth(jk).GE.dep_min)
            selected=selected.AND.(depth(jk).LE.dep_max)
          ENDIF
          IF (selected) THEN
!
! Read time irregular grid
          IF ((tim_irr).AND.(idtimz.NE.0)) THEN
            tstart(idtimz)=jk
            IF (idtimt.NE.0) tstart(idtimt)=jt
!
            IF (idtimx.EQ.0) THEN
! grids: yz, yzt
              ierr = NF90_GET_VAR(idf,idtim,timexy(1,1:jpj),start=tstart, &
     &                            count=tcount,stride=tstrid)
              IF (ierr.NE.0) GOTO 105
              DO ji=2,jpi
                timexy(ji,1:jpj)=timexy(1,1:jpj)
              ENDDO
            ELSEIF (idtimy.EQ.0) THEN
! grids: xz, xzt
              ierr = NF90_GET_VAR(idf,idtim,timexy(1:jpi,1),start=tstart, &
     &                            count=tcount,stride=tstrid)
              IF (ierr.NE.0) GOTO 105
              DO jj=2,jpj
                timexy(1:jpi,jj)=timexy(1:jpi,1)
              ENDDO
            ELSE
! grids: xyz, xyzt
              ierr = NF90_GET_VAR(idf,idtim,timexy,start=tstart, &
     &                            count=tcount,stride=tstrid)
              IF (ierr.NE.0) GOTO 105
            ENDIF
!
          ENDIF
!
! Read depth irregular grid
          IF ((dep_irr).AND.(iddepz.NE.0)) THEN
            zstart(iddepz)=jk
            IF (iddept.NE.0) zstart(iddept)=jt
!
            IF (iddepx.EQ.0) THEN
! grids: yz, yzt
              ierr = NF90_GET_VAR(idf,iddep,depthxy(1,1:jpj),start=zstart, &
     &                            count=zcount,stride=zstrid)
              IF (ierr.NE.0) GOTO 105
              DO ji=2,jpi
                depthxy(ji,1:jpj)=depthxy(1,1:jpj)
              ENDDO
            ELSEIF (iddepy.EQ.0) THEN
! grids: xz, xzt
              ierr = NF90_GET_VAR(idf,iddep,depthxy(1:jpi,1),start=zstart, &
     &                            count=zcount,stride=zstrid)
              IF (ierr.NE.0) GOTO 105
              DO jj=2,jpj
                depthxy(1:jpi,jj)=depthxy(1:jpi,1)
              ENDDO
            ELSE
! grids: xyz, xyzt
              ierr = NF90_GET_VAR(idf,iddep,depthxy,start=zstart, &
     &                            count=zcount,stride=zstrid)
              IF (ierr.NE.0) GOTO 105
            ENDIF
!
          ENDIF
!
! Read latitude irregular grid
          IF ((lat_irr).AND.(idlatz.NE.0)) THEN
            ystart(idlatz)=jk
            IF (idlatt.NE.0) ystart(idlatt)=jt
!
            IF (idlatx.EQ.0) THEN
! grids: yz, yzt
              ierr = NF90_GET_VAR(idf,idlat,latxy(1,1:jpj),start=ystart, &
     &                            count=ycount,stride=ystrid)
              IF (ierr.NE.0) GOTO 105
              DO ji=2,jpi
                latxy(ji,1:jpj)=latxy(1,1:jpj)
              ENDDO
            ELSEIF (idlaty.EQ.0) THEN
! grids: xz, xzt
              ierr = NF90_GET_VAR(idf,idlat,latxy(1:jpi,1),start=ystart, &
     &                            count=ycount,stride=ystrid)
              IF (ierr.NE.0) GOTO 105
              DO jj=2,jpj
                latxy(1:jpi,jj)=latxy(1:jpi,1)
              ENDDO
            ELSE
! grids: xyz, xyzt
              ierr = NF90_GET_VAR(idf,idlat,latxy,start=ystart, &
     &                            count=ycount,stride=ystrid)
              IF (ierr.NE.0) GOTO 105
            ENDIF
!
          ENDIF
!
! Read longitude irregular grid
          IF ((lon_irr).AND.(idlonz.NE.0)) THEN
            xstart(idlonz)=jk
            IF (idlont.NE.0) xstart(idlont)=jt
!
            IF (idlonx.EQ.0) THEN
! grids: yz, yzt
              ierr = NF90_GET_VAR(idf,idlon,lonxy(1,1:jpj),start=xstart, &
     &                            count=xcount,stride=xstrid)
              IF (ierr.NE.0) GOTO 105
              DO ji=2,jpi
                lonxy(ji,1:jpj)=lonxy(1,1:jpj)
              ENDDO
            ELSEIF (idlony.EQ.0) THEN
! grids: xz, xzt
              ierr = NF90_GET_VAR(idf,idlon,lonxy(1:jpi,1),start=xstart, &
     &                            count=xcount,stride=xstrid)
              IF (ierr.NE.0) GOTO 105
              DO jj=2,jpj
                lonxy(1:jpi,jj)=lonxy(1:jpi,1)
              ENDDO
            ELSE
! grids: xyz, xyzt
              ierr = NF90_GET_VAR(idf,idlon,lonxy,start=xstart, &
     &                            count=xcount,stride=xstrid)
              IF (ierr.NE.0) GOTO 105
            ENDIF
          ENDIF
!
! Read slice of data (along 1st and 2nd database dimensions)
          IF (idvz.NE.0) vstart(idvz)=jk
          IF (idvt.NE.0) vstart(idvt)=jt
!
          IF (idvx.EQ.0) THEN
! grids: y, yz, yzt
            ierr = NF90_GET_VAR(idf,idv,tabdbs(1,1:jpj),start=vstart, &
     &                          count=vcount,stride=vstrid)
            IF (ierr.NE.0) GOTO 104
            DO ji=2,jpi
              tabdbs(ji,1:jpj)=tabdbs(1,1:jpj)
            ENDDO
          ELSEIF (idvy.EQ.0) THEN
! grids: x, xz, xzt
            ierr = NF90_GET_VAR(idf,idv,tabdbs(1:jpi,1),start=vstart, &
     &                          count=vcount,stride=vstrid)
            IF (ierr.NE.0) GOTO 104
            DO jj=2,jpj
              tabdbs(1:jpi,jj)=tabdbs(1:jpi,1)
            ENDDO
          ELSE
! grids: xy, xyz, xyzt
            ierr = NF90_GET_VAR(idf,idv,tabdbs,start=vstart, &
     &                          count=vcount,stride=vstrid)
            IF (ierr.NE.0) GOTO 104
          ENDIF
!
! Loop on database 1st dimension (longitude in regular grid)
! Exclude outside data
          DO ji = 1,jpi
            selected=(ji.GE.ji_min).AND.(ji.LE.ji_max)
            IF (.NOT.lon_irr) THEN
              selected=selected.AND.(lon(ji).GE.lon_min)
              selected=selected.AND.(lon(ji).LE.lon_max)
            ENDIF
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
              selected=(jj.GE.jj_min).AND.(jj.LE.jj_max)
              IF (.NOT.lat_irr) THEN
                selected=selected.AND.(lat(jj).GE.lat_min)
                selected=selected.AND.(lat(jj).LE.lat_max)
              ENDIF
              IF (selected) THEN
!
! Select inside data
              inside=.TRUE.
              IF (tim_irr) THEN
                 otim=timexy(ji,jj)
                 inside=inside.AND.(otim.GE.tim_min)
                 inside=inside.AND.(otim.LE.tim_max)
              ELSE
                 otim=time(jt)
              ENDIF
              IF (dep_irr) THEN
                 odep=depthxy(ji,jj)
                 inside=inside.AND.(odep.GE.dep_min)
                 inside=inside.AND.(odep.LE.dep_max)
              ELSE
                 odep=depth(jk)
              ENDIF
              IF (lon_irr) THEN
                 olon=lonxy(ji,jj)
                 inside=inside.AND.(olon.GE.lon_min)
                 inside=inside.AND.(olon.LE.lon_max)
              ELSE
                 olon=lon(ji)
              ENDIF
              IF (lat_irr) THEN
                 olat=latxy(ji,jj)
                 inside=inside.AND.(olat.GE.lat_min)
                 inside=inside.AND.(olat.LE.lat_max)
              ELSE
                 olat=lat(jj)
              ENDIF
!
              selected=inside
              selected=selected.AND.(tabdbs(ji,jj).GE.dbs_min)
              selected=selected.AND.(tabdbs(ji,jj).LE.dbs_max)
              selected=selected.AND.(tabdbs(ji,jj).NE.dbsexcl)
!
              IF (selected) THEN
!
                IF (dtangrd(indobs).EQ.4) THEN
                  DO WHILE (olon.GE.FREAL4(442.))
                    olon = olon - FREAL4(360)
                  ENDDO
                  DO WHILE (olon.LT.FREAL4(78.))
                    olon = olon + FREAL4(360)
                  ENDDO
                ENDIF
!
! Add observation in data vector
                jpo = jpo + 1
                IF (jpo.GT.jpdbssize) GOTO 106
                kvectdbs(jpo) = FREAL(tabdbs(ji,jj))
!
                SELECT CASE(dtadim)
                CASE(2)
                   kgridij(jpo)%longi = olon
                   kgridij(jpo)%latj = olat
                CASE(3)
                   kgridijk(jpo)%longi = olon
                   kgridijk(jpo)%latj = olat
                   kgridijk(jpo)%levk = odep
                CASE(4)
                   kgrid(jpo)%lon = olon
                   kgrid(jpo)%lat = olat
                   kgrid(jpo)%dep = odep
                   kgrid(jpo)%tim = otim
                END SELECT
!
              ENDIF 
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
         WRITE(numout,*) ' number of observations: ',jpo
         WRITE(numout,*) ' sample of observations:'
         SELECT CASE(dtadim)
         CASE(2)
           WRITE(numout,*) ' index lon lat value'
           kform='(i9,3e12.3)'
         CASE(3)
           WRITE(numout,*) ' index lon lat depth value'
           kform='(i9,4e12.3)'
         CASE(4)
           WRITE(numout,*) ' index lon lat depth timr value'
           kform='(i9,5e12.3)'
         END SELECT
         DO jo=1,jpo,MAX(1,jpo/25)
            SELECT CASE(dtadim)
            CASE(2)
               WRITE(numout,kform) jo, kgridij(jo)%longi, &
     &            kgridij(jo)%latj, kvectdbs(jo)
            CASE(3)
               WRITE(numout,kform) jo, kgridijk(jo)%longi, &
     &            kgridijk(jo)%latj, kgridijk(jpo)%levk, &
     &            kvectdbs(jo)
            CASE(4)
               WRITE(numout,kform) jo, kgrid(jo)%lon, &
     &            kgrid(jo)%lat, kgrid(jpo)%dep, &
     &            kgrid(jo)%tim, kvectdbs(jo)
            END SELECT
         ENDDO
         WRITE(numout,*) ' '
      ENDIF
!
! --- deallocate
      IF (allocated(tabdbs)) deallocate(tabdbs)
      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)
      IF (allocated(depth)) deallocate(depth)
      IF (allocated(time)) deallocate(time)
      IF (allocated(lonxy)) deallocate(lonxy)
      IF (allocated(latxy)) deallocate(latxy)
      IF (allocated(depthxy)) deallocate(depthxy)
      IF (allocated(timexy)) deallocate(timexy)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioncdbs','readncdbs')
 1001 CALL printerror2(0,1001,3,'lioncdbs','readncdbs')
!
 101  WRITE (texterror,*) 'Bad NetCDF file: ',kfnindbs
      CALL printerror2(0,101,3,'lioncdbs','readncdbs',comment=texterror)
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,3,'lioncdbs','readncdbs',comment=texterror)
 103  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',kfnindbs
      CALL printerror2(0,103,3,'lioncdbs','readncdbs',comment=texterror)
 104  WRITE (texterror,*) 'Error reading NetCDF variable: ',dbsobsnam
      CALL printerror2(0,104,3,'lioncdbs','readncdbs',comment=texterror)
 105  WRITE (texterror,*) 'Error reading dimension variable:'
      CALL printerror2(0,105,1,'lioncdbs','readncdbs',comment=texterror)
 106  WRITE (texterror,*) 'Insufficient database vector size'
      CALL printerror2(0,106,1,'lioncdbs','readncdbs',comment=texterror)
 107  WRITE (texterror,*) 'Incorrect dimension variables'
      CALL printerror2(0,107,1,'lioncdbs','readncdbs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrncdbs (kfnindbs,kjobs,jpdbsout)
!---------------------------------------------------------------------
!
!  Purpose : Read header of NetCDF gridded
!  -------   observation data base file
!
!  Method : Evaluate maximum number of obs to extract
!  ------
!
!  Input :  kfnindbs : database filename
!  -----
!  Output : jpdbsout : database size
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnindbs
      INTEGER, intent(in) :: kjobs
      INTEGER, intent(out) :: jpdbsout
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=word80) :: kform
      INTEGER :: jpi,jpj,jpk,jpt,dbs_max
      INTEGER :: ji_min,jj_min,jk_min,jt_min
      INTEGER :: ji_max,jj_max,jk_max,jt_max
      INTEGER :: allocok,ierr,indobs,inddbs
      INTEGER :: idx, idy, idz, idt, idf
      CHARACTER(len=bgword) :: xdim, ydim, zdim, tdim
!----------------------------------------------------------------------
! Get observation and database index
!     print *, kfnindbs,kjobs
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../evalhdrdbs/evalhdrncdbs'
         WRITE(numout,*) '    ==> READING file ',kfnindbs(1:lenv(kfnindbs))
      ENDIF
!
! -1.- Read header only if maximum database size is not in SESAM configuration
! ----------------------------------------------------------------------------
! Get database size from SESAM configuration
      dbs_max=obs_siz_max(indobs,inddbs)
!
      IF (dbs_max.GT.0) THEN
!
! 1.1 Use database size as defined in the SESAM configuration
         jpdbsout = dbs_max
!
      ELSE
!
! 1.2 Read database file header
!
! Get database info from SESAM configuration
        xdim=obsxdim(indobs,inddbs)
        ydim=obsydim(indobs,inddbs)
        zdim=obszdim(indobs,inddbs)
        tdim=obstdim(indobs,inddbs)
!
        ji_min=obsimin(indobs,inddbs)
        ji_max=obsimax(indobs,inddbs)
        jj_min=obsjmin(indobs,inddbs)
        jj_max=obsjmax(indobs,inddbs)
        jk_min=obskmin(indobs,inddbs)
        jk_max=obskmax(indobs,inddbs)
        jt_min=obstmin(indobs,inddbs)
        jt_max=obstmax(indobs,inddbs)
!
! Open database file
        ierr = NF90_OPEN(kfnindbs,NF90_NOWRITE,idf)
        IF (ierr.NE.0) GOTO 101
!
! Get x size
        IF (xdim.NE.'none') THEN
          ierr = NF90_INQ_DIMID(idf,xdim,idx)
          IF (ierr.NE.0) GOTO 103
          ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=jpi)
          IF (ierr.NE.0) GOTO 103
        ELSE
          idx=-HUGE(idx) ; jpi=1
        ENDIF
!
! Get y size
        IF (ydim.NE.'none') THEN
          ierr = NF90_INQ_DIMID(idf,ydim,idy)
          IF (ierr.NE.0) GOTO 103
          ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=jpj)
          IF (ierr.NE.0) GOTO 103
        ELSE
          idy=-HUGE(idy) ; jpj=1
        ENDIF
!
! Get z size
        IF (zdim.NE.'none') THEN
          ierr = NF90_INQ_DIMID(idf,zdim,idz)
          IF (ierr.NE.0) GOTO 103
          ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=jpk)
          IF (ierr.NE.0) GOTO 103
        ELSE
          idz=-HUGE(idz) ; jpk=1
        ENDIF
!
! Get t size
        IF (tdim.NE.'none') THEN
          ierr = NF90_INQ_DIMID(idf,tdim,idt)
          IF (ierr.NE.0) GOTO 103
          ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=jpt)
          IF (ierr.NE.0) GOTO 103
        ELSE
          idt=-HUGE(idt) ; jpt=1
        ENDIF
!
! Control print
        IF (nprint.GE.2) THEN
           kform='(8x,a,3i5)'
           WRITE(numout,kform) '- isiz,imin,imax: ',jpi,ji_min,ji_max
           WRITE(numout,kform) '- jsiz,jmin,jmax: ',jpj,jj_min,jj_max
           WRITE(numout,kform) '- ksiz,kmin,kmax: ',jpk,jk_min,jk_max
           WRITE(numout,kform) '- tsiz,tmin,tmax: ',jpt,jt_min,jt_max
        ENDIF
!
! 1.3 Compute maximum database size
        jpi=MIN(ji_max-ji_min+1,jpi)
        jpj=MIN(jj_max-jj_min+1,jpj)
        jpk=MIN(jk_max-jk_min+1,jpk)
        jpt=MIN(jt_max-jt_min+1,jpt)
        jpdbsout = jpi * jpj * jpk * jpt
!
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioncdbs','evalhdrncdbs')
 1001 CALL printerror2(0,1001,3,'lioncdbs','evalhdrncdbs')
!
 101  WRITE (texterror,*) 'Bad NetCDF file: ',kfnindbs
      CALL printerror2(0,101,3,'lioncdbs','evalhdrncdbs',comment=texterror)
 103  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',kfnindbs
      CALL printerror2(0,103,3,'lioncdbs','evalhdrncdbs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE lioncdbs
