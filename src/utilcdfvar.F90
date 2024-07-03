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
! ---                    UTILCDFVAR                             ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-03 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE cdfrdim : Read dimensions and title from CDF file
! --- SUBROUTINE cdfrsli : Read 2D slice of data from CDF file
! --- SUBROUTINE cdfrtab : Read 3D array of data from CDF file
! --- SUBROUTINE cdfwdim : Write dimensions and title in CDF file
! --- SUBROUTINE cdfwsli : Write 2D slice of data in CDF file
! --- SUBROUTINE cdfwtab : Write 3D array of data in CDF file
! --- SUBROUTINE cdfwvar : Write variable description in CDF file
! --- SUBROUTINE cdfrloc : Read regular grid location from CDF file
! --- SUBROUTINE cdfwloc : Write regular grid location in CDF file
! --- SUBROUTINE cdfrvar : Read variable description from CDF file
! --- SUBROUTINE cdfrtim : Read time values from CDF file
! --- SUBROUTINE cdfrpos : Read irregular grid location from CDF file
! --- SUBROUTINE cdfwpos : Write irregular grid location in CDF file
! ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilcdfvar
      use mod_main
      use netcdf
      IMPLICIT NONE
      include 'netcdf.inc'
      PRIVATE

      PUBLIC cdfrdim,cdfrsli,cdfrtab,cdfwdim,cdfwsli,cdfwtab
      PUBLIC cdfwvar,cdfrloc,cdfwloc,cdfrvar,cdfrtim,cdfrpos,cdfwpos

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfrdim(filename,imax,jmax,kmax,tmax,title)

      implicit none

      integer*4 imax,jmax,kmax,tmax
      character*(*) filename,title

      character*200 name
      integer*4 iid,jid,kid,tid
      integer*4 ncid, rcode, ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'lon',rcode)
      jid  = ncdid (ncid,'lat',rcode)
      kid  = ncdid (ncid,'depth',rcode)
      tid  = ncdid (ncid,'time',rcode)

      call ncdinq (ncid,iid,name,imax,rcode)
      call ncdinq (ncid,jid,name,jmax,rcode)
      call ncdinq (ncid,kid,name,kmax,rcode)
      call ncdinq (ncid,tid,name,tmax,rcode)

      ll=len(title)
      call ncagtc (ncid,ncglobal,'title',title,ll,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfrsli(filename,varname,dcut,icut,idx,c)

      implicit none

      character*(*) filename,varname
      integer*4 dcut,icut
      BIGREAL4 c(*)

      character*200 name
      integer*4 imax,jmax,kmax,idx
      integer*4 iid,jid,kid,tid,cid
      integer*4 ncid, rcode, ll
      integer*4 start(4),count(4)

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'lon',rcode)
      jid  = ncdid (ncid,'lat',rcode)
      kid  = ncdid (ncid,'depth',rcode)
      tid  = ncdid (ncid,'time',rcode)

      cid  = ncvid (ncid,varname,rcode)

      call ncdinq (ncid,iid,name,imax,rcode)
      call ncdinq (ncid,jid,name,jmax,rcode)
      call ncdinq (ncid,kid,name,kmax,rcode)

! --- get variable ---

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = idx

      count(1) = imax
      count(2) = jmax
      count(3) = kmax
      count(4) = 1

      start(dcut) = icut
      count(dcut) = 1
      call ncvgt  (ncid,cid,start,count,c,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfrtab(filename,varname,idx,c)

      implicit none

      character*(*) filename,varname
      integer*4 idx
      BIGREAL4 c(*)

      character*200 name
      integer*4 imax,jmax,kmax,tmax,ndim
      integer*4 iid,jid,kid,tid,cid
      integer*4 ncid, rcode, ll
      integer*4 start(4),count(4), vtype, cdim(4), natt

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'lon',rcode)
      jid  = ncdid (ncid,'lat',rcode)
      kid  = ncdid (ncid,'depth',rcode)
      tid  = ncdid (ncid,'time',rcode)

      cid  = ncvid (ncid,varname,rcode)

      call ncvinq (ncid,cid,name,vtype,ndim,cdim,natt,rcode)
      ndim = ndim - 1

      call ncdinq (ncid,iid,name,imax,rcode)
      call ncdinq (ncid,jid,name,jmax,rcode)
      call ncdinq (ncid,kid,name,kmax,rcode)
      call ncdinq (ncid,tid,name,tmax,rcode)

! --- get variable ---

      if (ndim.eq.1) then

        start(1) = 1
        start(2) = idx

        count (1) = imax
        count (2) = 1

      elseif (ndim.eq.2) then

        start(1) = 1
        start(2) = 1
        start(3) = idx

        count (1) = imax
        count (2) = jmax
        count (3) = 1

      else

        start(1) = 1
        start(2) = 1
        start(3) = 1
        start(4) = idx

        count (1) = imax
        count (2) = jmax
        count (3) = kmax
        count (4) = 1

      endif

      call ncvgt  (ncid,cid,start,count,c,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfwdim(filename,imax,jmax,kmax,title)

      implicit none

      integer*4 imax,jmax,kmax
      character*(*) filename,title

      integer*4 iid,jid,kid,tid
      integer*4 ncid, rcode, ll

! --- create NetCDF file ---

      ll=lenv(filename)
      ncid = nccre (filename(1:ll),ncnoclob,rcode)

! --- create dimensions ---

      iid    = ncddef (ncid,'lon',imax,rcode)
      jid    = ncddef (ncid,'lat',jmax,rcode)
      kid    = ncddef (ncid,'depth',kmax,rcode)
      tid    = ncddef (ncid,'time',ncunlim,rcode)

! --- create global attributes ---

      ll=lenv(title)
      call ncaptc (ncid,ncglobal,'title',ncchar,ll,title,rcode)

! --- Ending define mode ---

      call ncendf (ncid,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfwsli(filename,varname,dcut,icut,idx,time,c)

      implicit none

      character*(*) filename,varname
      integer*4 dcut,icut,idx
      BIGREAL4 c(*),time

      character*200 name
      integer*4 imax,jmax,kmax,tmax
      integer*4 start(4),count(4)
      integer*4 iid,jid,kid,tid,cid
      integer*4 ncid,rcode,ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'lon',rcode)
      jid  = ncdid (ncid,'lat',rcode)
      kid  = ncdid (ncid,'depth',rcode)
      tid  = ncdid (ncid,'time',rcode)

      call ncdinq (ncid,iid,name,imax,rcode)
      call ncdinq (ncid,jid,name,jmax,rcode)
      call ncdinq (ncid,kid,name,kmax,rcode)
      call ncdinq (ncid,tid,name,tmax,rcode)

      ll=lenv(varname)
      cid  = ncvid (ncid,varname(1:ll),rcode)
      tid  = ncvid (ncid,'time',rcode)

! --- if time index does not already exist

      if ( (idx.le.0) .or. (idx.gt.tmax) ) then

        idx = tmax + 1

! --- add one time index to time unlimited dimension

        start (1) = idx
        count (1) = 1
        call ncvpt (ncid,tid,start,count,time,rcode)

      endif

! --- put variable ---

      start (1) = 1
      start (2) = 1
      start (3) = 1
      start (4) = idx

      count (1) = imax
      count (2) = jmax
      count (3) = kmax
      count (4) = 1

      start(dcut) = icut
      count(dcut) = 1

      call ncvpt (ncid,cid,start,count,c,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfwtab(filename,varname,idx,time,c)

      implicit none

      character*(*) filename,varname
      integer*4 idx
      BIGREAL4 c(*),time

      character*200 name
      integer*4 imax,jmax,kmax,tmax,ndim
      integer*4 start(4),count(4)
      integer*4 iid,jid,kid,tid,cid
      integer*4 ncid,rcode,vtype,cdim(4),natt,ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'lon',rcode)
      jid  = ncdid (ncid,'lat',rcode)
      kid  = ncdid (ncid,'depth',rcode)
      tid  = ncdid (ncid,'time',rcode)

      call ncdinq (ncid,iid,name,imax,rcode)
      call ncdinq (ncid,jid,name,jmax,rcode)
      call ncdinq (ncid,kid,name,kmax,rcode)
      call ncdinq (ncid,tid,name,tmax,rcode)

      ll=lenv(varname)
      cid  = ncvid (ncid,varname(1:ll),rcode)
      tid  = ncvid (ncid,'time',rcode)

      call ncvinq (ncid,cid,name,vtype,ndim,cdim,natt,rcode)
      ndim = ndim - 1

! --- if time index does not already exist

      if ( (idx.le.0) .or. (idx.gt.tmax) ) then

        idx = tmax + 1

! --- add one time index to time unlimited dimension variable

        start (1) = idx
        count (1) = 1
        call ncvpt (ncid,tid,start,count,time,rcode)

      endif

! --- put variable ---

      if (ndim.eq.1) then

        start (1) = 1
        start (2) = idx

        count (1) = imax
        count (2) = 1

      elseif (ndim.eq.2) then

        start (1) = 1
        start (2) = 1
        start (3) = idx

        count (1) = imax
        count (2) = jmax
        count (3) = 1

      else

        start (1) = 1
        start (2) = 1
        start (3) = 1
        start (4) = idx

        count (1) = imax
        count (2) = jmax
        count (3) = kmax
        count (4) = 1

      endif

      call ncvpt (ncid,cid,start,count,c,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfwvar(filename,varname,spval,unit,name,ndim)

      implicit none

      BIGREAL4 spval
      character*(*) filename,varname,unit,name
      integer*4 ndim

      integer*4 cdim(4)
      integer*4 iid,jid,kid,tid,cid
      integer*4 ncid, rcode, ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- put file into define mode

      call ncredf (ncid,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'lon',rcode)
      jid  = ncdid (ncid,'lat',rcode)
      kid  = ncdid (ncid,'depth',rcode)
      tid  = ncdid (ncid,'time',rcode)

! --- create variable ---

      if (ndim.eq.1) then

        cdim (1) = iid
        cdim (2) = tid

        ll=lenv(varname)
        cid = ncvdef (ncid,varname(1:ll),ncfloat,2,cdim,rcode)

      elseif (ndim.eq.2) then

        cdim (1) = iid
        cdim (2) = jid
        cdim (3) = tid

        ll=lenv(varname)
        cid = ncvdef (ncid,varname(1:ll),ncfloat,3,cdim,rcode)

      else

        cdim (1) = iid
        cdim (2) = jid
        cdim (3) = kid
        cdim (4) = tid

        ll=lenv(varname)
        cid = ncvdef (ncid,varname(1:ll),ncfloat,4,cdim,rcode)

      endif

! --- create attributes ---

      call ncapt  (ncid,cid,'missing_value',ncfloat,1,spval,rcode)
      ll=lenv(unit)
      call ncaptc (ncid,cid,'units',ncchar,ll,unit,rcode)
      ll=lenv(name)
      call ncaptc (ncid,cid,'long_name',ncchar,ll,name,rcode)

! --- Ending define mode ---

      call ncendf (ncid,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfrloc(filename,lon,lat,depth,unit)

      implicit none

      character*(*) filename,unit(4)
      BIGREAL4 spval,lon(*),lat(*),depth(*)

      character*200 name
      integer*4 imax,jmax,kmax,tmax,i
      integer*4 start(1),count(1),cdim(1)
      integer*4 iid,jid,kid,tid,lid,mid,nid,oid
      integer*4 ncid, rcode,ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'lon',rcode)
      jid  = ncdid (ncid,'lat',rcode)
      kid  = ncdid (ncid,'depth',rcode)
      tid  = ncdid (ncid,'time',rcode)

      lid  = ncvid (ncid,'lon',rcode)
      mid  = ncvid (ncid,'lat',rcode)
      nid  = ncvid (ncid,'depth',rcode)
      oid  = ncvid (ncid,'time',rcode)

      call ncdinq (ncid,iid,name,imax,rcode)
      call ncdinq (ncid,jid,name,jmax,rcode)
      call ncdinq (ncid,kid,name,kmax,rcode)
      call ncdinq (ncid,tid,name,tmax,rcode)

! --- get attributes ---
      ll=len(unit(1))
      call ncagtc (ncid,lid,'units',unit(1),ll,rcode)
      ll=len(unit(2))
      call ncagtc (ncid,mid,'units',unit(2),ll,rcode)
      ll=len(unit(3))
      call ncagtc (ncid,nid,'units',unit(3),ll,rcode)
      ll=len(unit(4))
      call ncagtc (ncid,oid,'units',unit(4),ll,rcode)

! --- get variables ---

      start(1) = 1
      count (1) = imax
      call ncvgt  (ncid,lid,start,count,lon,rcode)

      IF (unit(1)(1:3).EQ.'deg') THEN
      DO i=1,imax
         DO WHILE (lon(i).GT.FREAL4(180.))
            lon(i) = lon(i) - FREAL4(360)
         ENDDO
         DO WHILE (lon(i).LE.FREAL4(-180.))
            lon(i) = lon(i) + FREAL4(360)
         ENDDO
      ENDDO
      ENDIF

      start(1) = 1
      count (1) = jmax
      call ncvgt  (ncid,mid,start,count,lat,rcode)

      start(1) = 1
      count (1) = kmax
      call ncvgt  (ncid,nid,start,count,depth,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfwloc(filename,lon,lat,depth,unit)

      implicit none

      character*(*) filename,unit(4)
      BIGREAL4 spval,lon(*),lat(*),depth(*)

      character*200 name
      integer*4 imax,jmax,kmax
      integer*4 start(1),count(1),cdim(1)
      integer*4 iid,jid,kid,tid,lid,mid,nid,oid
      integer*4 ncid, rcode,ll,i

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- put file into define mode

      call ncredf (ncid,rcode)

! --- get information on NetCDF file (dimensions) ---

      iid  = ncdid (ncid,'lon',rcode)
      jid  = ncdid (ncid,'lat',rcode)
      kid  = ncdid (ncid,'depth',rcode)
      tid  = ncdid (ncid,'time',rcode)

      call ncdinq (ncid,iid,name,imax,rcode)
      call ncdinq (ncid,jid,name,jmax,rcode)
      call ncdinq (ncid,kid,name,kmax,rcode)

! --- create dimension variable ---

      cdim (1) = iid
      lid = ncvdef (ncid,'lon',ncfloat,1,cdim,rcode)
      cdim (1) = jid
      mid = ncvdef (ncid,'lat',ncfloat,1,cdim,rcode)
      cdim (1) = kid
      nid = ncvdef (ncid,'depth',ncfloat,1,cdim,rcode)
      cdim (1) = tid
      oid = ncvdef (ncid,'time',ncfloat,1,cdim,rcode)

! --- create attributes ---

      ll=lenv(unit(1))
      call ncaptc (ncid,lid,'units',ncchar,ll,unit(1),rcode)
      ll=lenv(unit(2))
      call ncaptc (ncid,mid,'units',ncchar,ll,unit(2),rcode)
      ll=lenv(unit(3))
      call ncaptc (ncid,nid,'units',ncchar,ll,unit(3),rcode)
      ll=lenv(unit(4))
      call ncaptc (ncid,oid,'units',ncchar,ll,unit(4),rcode)

! --- Ending define mode ---

      call ncendf (ncid,rcode)

! --- put variables ---

      start (1) = 1
      count (1) = imax
      call ncvpt (ncid,lid,start,count,lon,rcode)

      start (1) = 1
      count (1) = jmax
      call ncvpt (ncid,mid,start,count,lat,rcode)

      start (1) = 1
      count (1) = kmax
      call ncvpt (ncid,nid,start,count,depth,rcode)

! --- close NetCDF file

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfrvar(filename,varname,spval,unit,name,ndim)

      implicit none

      character*(*) filename,varname,unit,name
      BIGREAL4 spval
      integer ndim

      integer ncid, rcode, vtype, cdim(4), natt, cid, ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      cid  = ncvid (ncid,varname,rcode)

      call ncvinq (ncid,cid,name,vtype,ndim,cdim,natt,rcode)
      ndim = ndim - 1

      call ncagt  (ncid,cid,'missing_value',spval,rcode)
      ll=len(unit)
      call ncagtc (ncid,cid,'units',unit,ll,rcode)
      ll=len(name)
      call ncagtc (ncid,cid,'long_name',name,ll,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfrtim(filename,time)

      implicit none

      character*(*) filename
      BIGREAL4 time(*)

      character*200 name
      integer start(1),count(1)
      integer tid,oid
      integer ncid, rcode, tmax, ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      tid  = ncdid (ncid,'time',rcode)
      oid  = ncvid (ncid,'time',rcode)

      call ncdinq (ncid,tid,name,tmax,rcode)

! --- get variables ---

      start(1) = 1
      count (1) = tmax
      call ncvgt  (ncid,oid,start,count,time,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfrpos(filename,lon,lat,spval,unit)

      implicit none

      character*(*) filename,unit
      BIGREAL4 lon(*),lat(*),spval

      character*200 name
      integer imax,jmax,i
      integer iid,jid,tid,cid,did
      integer ncid, rcode, ll
      integer start(2),count(2)

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'lon',rcode)
      jid  = ncdid (ncid,'lat',rcode)

      cid  = ncvid (ncid,'lonxy',rcode)
      did  = ncvid (ncid,'latxy',rcode)

      call ncdinq (ncid,iid,name,imax,rcode)
      call ncdinq (ncid,jid,name,jmax,rcode)

! --- get variable ---

      start(1) = 1
      start(2) = 1

      count (1) = imax
      count (2) = jmax

      call ncvgt  (ncid,cid,start,count,lon,rcode)
      call ncvgt  (ncid,did,start,count,lat,rcode)

! --- get attributes ---

      call ncagt  (ncid,cid,'missing_value',spval,rcode)
      ll=len(unit)
      call ncagtc (ncid,cid,'units',unit,ll,rcode)

      IF (unit(1:3).EQ.'deg') THEN
      DO i=1,imax*jmax
         DO WHILE (lon(i).GT.FREAL4(180.))
            lon(i) = lon(i) - FREAL4(360)
         ENDDO
         DO WHILE (lon(i).LE.FREAL4(-180.))
            lon(i) = lon(i) + FREAL4(360)
         ENDDO
      ENDDO
      ENDIF

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfwpos(filename,lon,lat,spval,unit)

      implicit none

      character*(*) filename,unit
      BIGREAL4 spval,lon(*),lat(*)

      character*200 name
      integer imax,jmax
      integer cdim(2),start(2),count(2)
      integer iid,jid,lid,mid
      integer ncid, rcode, ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'lon',rcode)
      jid  = ncdid (ncid,'lat',rcode)

      call ncdinq (ncid,iid,name,imax,rcode)
      call ncdinq (ncid,jid,name,jmax,rcode)

! --- put file into define mode

      call ncredf (ncid,rcode)

! --- create variables ---

      cdim (1) = iid
      cdim (2) = jid

      lid = ncvdef (ncid,'lonxy',ncfloat,2,cdim,rcode)
      mid = ncvdef (ncid,'latxy',ncfloat,2,cdim,rcode)

! --- create attributes ---

      call ncapt  (ncid,lid,'missing_value',ncfloat,1,spval,rcode)
      call ncapt  (ncid,mid,'missing_value',ncfloat,1,spval,rcode)

      ll=lenv(unit)
      call ncaptc (ncid,lid,'units',ncchar,ll,unit,rcode)
      call ncaptc (ncid,mid,'units',ncchar,ll,unit,rcode)

! --- Ending define mode ---

      call ncendf (ncid,rcode)

! --- put variables ---

      start (1) = 1
      start (2) = 1

      count (1) = imax
      count (2) = jmax

      call ncvpt (ncid,lid,start,count,lon,rcode)
      call ncvpt (ncid,mid,start,count,lat,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilcdfvar
