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
! ---                    UTILCDFZON                             ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-11 (J.M. Brankart)                      ---
! --- revised      : 00-03 (C.E. Testut)                        ---
! --- revised      : 03-04 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE cdfrdimzon    : Read dimensions from CZON file
! --- SUBROUTINE cdfrhdrzon    : Read header from CZON file
! --- SUBROUTINE cdfrptzon     : Read pointers from CZON file
! --- SUBROUTINE cdfrbubidxzon : Read local data section index
! ---                            from CZON file
! --- SUBROUTINE cdfrbubzon    : Read a set of local data sections
! ---                            from CZON file
! --- SUBROUTINE cdfwdimzon    : Write dimensions in CZON file
! --- SUBROUTINE cdfwhdrzon    : Write header in CZON file
! --- SUBROUTINE cdfwptzon     : Write pointers in CZON file
! --- SUBROUTINE cdfwbubzon    : Write a set of local data sections
! ---                            in CZON file
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilcdfzon
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC cdfrdimzon,cdfrhdrzon,cdfrptzon,cdfrbubidxzon
      PUBLIC cdfrbubzon,cdfwdimzon,cdfwhdrzon,cdfwptzon,cdfwbubzon

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfrdimzon(filename,jpi,jpj,jpk,jpt,jpbub,jpz, &
     &     dtaend,namelength,title)

      implicit none
      include 'netcdf.inc'

      integer jpi,jpj,jpk,jpt,jpbub,jpz,dtaend,namelength
      character*(*) filename,title

      character*200 name
      integer iid,jid,kid,tid,bubid,zid,dtaid,lid
      integer varid, ncid, rcode, ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      iid    = ncdid (ncid,'lon',rcode)
      jid    = ncdid (ncid,'lat',rcode)
      kid    = ncdid (ncid,'depth',rcode)
      tid    = ncdid (ncid,'time',rcode)
      bubid  = ncdid (ncid,'bubidx',rcode)
      zid    = ncdid (ncid,'zoneidx',rcode)
      dtaid  = ncdid (ncid,'dta',rcode)
      lid    = ncdid (ncid,'namelength',rcode)

      call ncdinq (ncid,iid,name,jpi,rcode)
      call ncdinq (ncid,jid,name,jpj,rcode)
      call ncdinq (ncid,kid,name,jpk,rcode)
      call ncdinq (ncid,tid,name,jpt,rcode)
      call ncdinq (ncid,bubid,name,jpbub,rcode)
      call ncdinq (ncid,zid,name,jpz,rcode)
      call ncdinq (ncid,dtaid,name,dtaend,rcode)
      call ncdinq (ncid,lid,name,namelength,rcode)

! --- get global attribute ---

      ll=len(title)
      call ncagtc (ncid,ncglobal,'title',title,ll,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfrhdrzon(filename,nam,dim,nbr,moy,ect)

      implicit none
      include 'netcdf.inc'

      character*(*) filename
      BIGREAL4 moy(*),ect(*)
      integer dim(*),nbr(*)
      character*(*) nam(*)

      character*200 dimname
      integer dtaid,lid
      integer moyid,ectid,dimid,nbrid,namid
      integer dtaend,namelength,namell
      integer ncid, rcode, ll
      integer start(2),count(2)

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      dtaid  = ncdid (ncid,'dta',rcode)
      lid  = ncdid (ncid,'namelength',rcode)
      call ncdinq (ncid,dtaid,dimname,dtaend,rcode)
      call ncdinq (ncid,lid,dimname,namelength,rcode)

      moyid  = ncvid (ncid,'mean',rcode)
      ectid  = ncvid (ncid,'std',rcode)
      dimid  = ncvid (ncid,'dim',rcode)
      nbrid  = ncvid (ncid,'nbr',rcode)
      namid  = ncvid (ncid,'name',rcode)

! --- get variable ---

      start(1) = 1
      count(1) = dtaend

      call ncvgt  (ncid,moyid,start,count,moy,rcode)
      call ncvgt  (ncid,ectid,start,count,ect,rcode)
      call ncvgt  (ncid,dimid,start,count,dim,rcode)
      call ncvgt  (ncid,nbrid,start,count,nbr,rcode)

      namell = len(nam)
      start(1) = 1
      count(1) = namell
      start(2) = 1
      count(2) = dtaend

      call ncvgtc (ncid,namid,start,count,nam,namell*dtaend,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfrptzon(filename,ptbubidx,ptdtalon,ptdtalat,ptdtadepth, &
     &   ptdtatime,ptbublon,ptbublat,ptbubdepth,ptbubtime,jz0,jz1,jd0,jd1)

      implicit none
      include 'netcdf.inc'

      character*(*) filename
      integer ptbubidx(*)
      integer ptdtalon(*),ptdtalat(*),ptdtadepth(*),ptdtatime(*)
      integer ptbublon(*),ptbublat(*),ptbubdepth(*),ptbubtime(*)
      integer jz0,jz1,jd0,jd1

      integer start(2),count(2)
      integer ncid, rcode, ll, varid

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- get variables ---

      start (1) = jz0
      count (1) = jz1 - jz0 + 1
      start (2) = jd0
      count (2) = jd1 - jd0 + 1

      varid  = ncvid (ncid,'ptbubidx',rcode)
      call ncvgt (ncid,varid,start,count,ptbubidx,rcode)

      varid  = ncvid (ncid,'ptdtalon',rcode)
      call ncvgt (ncid,varid,start,count,ptdtalon,rcode)

      varid  = ncvid (ncid,'ptdtalat',rcode)
      call ncvgt (ncid,varid,start,count,ptdtalat,rcode)

      varid  = ncvid (ncid,'ptdtadepth',rcode)
      call ncvgt (ncid,varid,start,count,ptdtadepth,rcode)

      varid  = ncvid (ncid,'ptdtatime',rcode)
      call ncvgt (ncid,varid,start,count,ptdtatime,rcode)

      varid  = ncvid (ncid,'ptbublon',rcode)
      call ncvgt (ncid,varid,start,count,ptbublon,rcode)

      varid  = ncvid (ncid,'ptbublat',rcode)
      call ncvgt (ncid,varid,start,count,ptbublat,rcode)

      varid  = ncvid (ncid,'ptbubdepth',rcode)
      call ncvgt (ncid,varid,start,count,ptbubdepth,rcode)

      varid  = ncvid (ncid,'ptbubtime',rcode)
      call ncvgt (ncid,varid,start,count,ptbubtime,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfrbubidxzon(filename,ptbubidx,jz0,jz1,jd0,jd1)

      implicit none
      include 'netcdf.inc'

      character*(*) filename
      integer ptbubidx(*)
      integer jz0,jz1,jd0,jd1

      integer start(2),count(2)
      integer ncid, rcode, varid, ll


! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- get variables ---

      start (1) = jz0
      count (1) = jz1 - jz0 + 1
      start (2) = jd0
      count (2) = jd1 - jd0 + 1

      varid  = ncvid (ncid,'ptbubidx',rcode)
      call ncvgt (ncid,varid,start,count,ptbubidx,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfrbubzon(filename,bub,bubidx,bubcount)

      implicit none
      include 'netcdf.inc'

      character*(*) filename
      BIGREAL4 bub(*)
      integer bubidx,bubcount

      character*200 name
      integer jpi,jpj,jpk,jpt,jpbub
      integer*4 start(5),count(5)
      integer*4 iid,jid,kid,tid,bubid,cid
      integer*4 ncid,rcode,ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- get information on NetCDF file ---

      iid    = ncdid (ncid,'lon',rcode)
      jid    = ncdid (ncid,'lat',rcode)
      kid    = ncdid (ncid,'depth',rcode)
      tid    = ncdid (ncid,'time',rcode)
      bubid  = ncdid (ncid,'bubidx',rcode)

      call ncdinq (ncid,iid,name,jpi,rcode)
      call ncdinq (ncid,jid,name,jpj,rcode)
      call ncdinq (ncid,kid,name,jpk,rcode)
      call ncdinq (ncid,tid,name,jpt,rcode)
      call ncdinq (ncid,bubid,name,jpbub,rcode)

      cid  = ncvid (ncid,'bubble',rcode)

! --- get variable ---

      start (1) = 1
      start (2) = 1
      start (3) = 1
      start (4) = 1
      start (5) = bubidx

      count (1) = jpi
      count (2) = jpj
      count (3) = jpk
      count (4) = jpt
      count (5) = bubcount

      call ncvgt (ncid,cid,start,count,bub,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfwdimzon(filename,jpi,jpj,jpk,jpt,jpbub,jpz, &
     &                   dtaend,namelength,title)

      implicit none
      include 'netcdf.inc' 

      integer jpi,jpj,jpk,jpt,jpbub,jpz,dtaend,namelength
      character*(*) filename,title

      integer iid,jid,kid,tid,bubid,zid,dtaid,lid,cdim(5)
      integer varid, ncid, rcode, ll

! --- create NetCDF file ---

      ll=lenv(filename)
      ncid = nccre (filename(1:ll),ncnoclob,rcode)

! --- create dimensions ---

      iid    = ncddef (ncid,'lon',jpi,rcode)
      jid    = ncddef (ncid,'lat',jpj,rcode)
      kid    = ncddef (ncid,'depth',jpk,rcode)
      tid    = ncddef (ncid,'time',jpt,rcode)
      bubid  = ncddef (ncid,'bubidx',ncunlim,rcode)
      zid    = ncddef (ncid,'zoneidx',jpz,rcode)
      dtaid  = ncddef (ncid,'dta',dtaend,rcode)
      lid  = ncddef (ncid,'namelength',namelength,rcode)

! --- create variable ---

      cdim (1) = iid
      cdim (2) = jid
      cdim (3) = kid
      cdim (4) = tid
      cdim (5) = bubid
      varid = ncvdef (ncid,'bubble',ncfloat,5,cdim,rcode)

      cdim (1) = lid
      cdim (2) = dtaid
      varid = ncvdef (ncid,'name',ncchar,2,cdim,rcode)

      cdim (1) = dtaid
      varid = ncvdef (ncid,'dim',nclong,1,cdim,rcode)
      varid = ncvdef (ncid,'nbr',nclong,1,cdim,rcode)
      varid = ncvdef (ncid,'mean',ncfloat,1,cdim,rcode)
      varid = ncvdef (ncid,'std',ncfloat,1,cdim,rcode)

      cdim (1) = zid
      cdim (2) = dtaid
      varid = ncvdef (ncid,'ptbubidx',nclong,2,cdim,rcode)
      varid = ncvdef (ncid,'ptdtalon',nclong,2,cdim,rcode)
      varid = ncvdef (ncid,'ptdtalat',nclong,2,cdim,rcode)
      varid = ncvdef (ncid,'ptdtadepth',nclong,2,cdim,rcode)
      varid = ncvdef (ncid,'ptdtatime',nclong,2,cdim,rcode)
      varid = ncvdef (ncid,'ptbublon',nclong,2,cdim,rcode)
      varid = ncvdef (ncid,'ptbublat',nclong,2,cdim,rcode)
      varid = ncvdef (ncid,'ptbubdepth',nclong,2,cdim,rcode)
      varid = ncvdef (ncid,'ptbubtime',nclong,2,cdim,rcode)

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
      SUBROUTINE cdfwhdrzon(filename,nam,dim,nbr,moy,ect)

      implicit none
      include 'netcdf.inc'

      character*(*) filename
      BIGREAL4 moy(*),ect(*)
      integer dim(*),nbr(*)
      character*(*) nam(*)

      character*200 dimname
      integer start(2),count(2)
      integer ncid, rcode, dtaid, lid
      integer moyid,ectid,dimid,nbrid,namid
      integer dtaend,namelength,namell,ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- get information on NetCDF file ---

      dtaid  = ncdid (ncid,'dta',rcode)
      lid  = ncdid (ncid,'namelength',rcode)
      call ncdinq (ncid,dtaid,dimname,dtaend,rcode)
      call ncdinq (ncid,lid,dimname,namelength,rcode)

      moyid  = ncvid (ncid,'mean',rcode)
      ectid  = ncvid (ncid,'std',rcode)
      dimid  = ncvid (ncid,'dim',rcode)
      nbrid  = ncvid (ncid,'nbr',rcode)
      namid  = ncvid (ncid,'name',rcode)

! --- put variables ---

      start (1) = 1
      count (1) = dtaend
      call ncvpt (ncid,moyid,start,count,moy,rcode)
      call ncvpt (ncid,ectid,start,count,ect,rcode)
      call ncvpt (ncid,dimid,start,count,dim,rcode)
      call ncvpt (ncid,nbrid,start,count,nbr,rcode)

      namell = len(nam)
      start (1) = 1
      count (1) = namell
      start (2) = 1
      count (2) = dtaend
      call ncvptc (ncid,namid,start,count,nam,namell*dtaend,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfwptzon(filename,ptbubidx,ptdtalon,ptdtalat,ptdtadepth, &
     &   ptdtatime,ptbublon,ptbublat,ptbubdepth,ptbubtime,jz0,jz1,jd0,jd1)

      implicit none
      include 'netcdf.inc'

      character*(*) filename
      integer ptbubidx(*)
      integer ptdtalon(*),ptdtalat(*),ptdtadepth(*),ptdtatime(*)
      integer ptbublon(*),ptbublat(*),ptbubdepth(*),ptbubtime(*)
      integer jz0,jz1,jd0,jd1

      integer start(2),count(2)
      integer ncid, rcode, varid, ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- put variables ---

      start (1) = jz0
      count (1) = jz1 - jz0 + 1
      start (2) = jd0
      count (2) = jd1 - jd0 + 1

      varid  = ncvid (ncid,'ptbubidx',rcode)
      call ncvpt (ncid,varid,start,count,ptbubidx,rcode)

      varid  = ncvid (ncid,'ptdtalon',rcode)
      call ncvpt (ncid,varid,start,count,ptdtalon,rcode)

      varid  = ncvid (ncid,'ptdtalat',rcode)
      call ncvpt (ncid,varid,start,count,ptdtalat,rcode)

      varid  = ncvid (ncid,'ptdtadepth',rcode)
      call ncvpt (ncid,varid,start,count,ptdtadepth,rcode)

      varid  = ncvid (ncid,'ptdtatime',rcode)
      call ncvpt (ncid,varid,start,count,ptdtatime,rcode)

      varid  = ncvid (ncid,'ptbublon',rcode)
      call ncvpt (ncid,varid,start,count,ptbublon,rcode)

      varid  = ncvid (ncid,'ptbublat',rcode)
      call ncvpt (ncid,varid,start,count,ptbublat,rcode)

      varid  = ncvid (ncid,'ptbubdepth',rcode)
      call ncvpt (ncid,varid,start,count,ptbubdepth,rcode)

      varid  = ncvid (ncid,'ptbubtime',rcode)
      call ncvpt (ncid,varid,start,count,ptbubtime,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cdfwbubzon(filename,bub,bubidx,bubcount)

      implicit none
      include 'netcdf.inc'

      character*(*) filename
      BIGREAL4 bub(*)
      integer bubidx,bubcount

      character*200 name
      integer jpi,jpj,jpk,jpt,jpbub
      integer*4 start(5),count(5)
      integer*4 iid,jid,kid,tid,bubid,cid
      integer*4 ncid,rcode,ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- get information on NetCDF file ---

      iid    = ncdid (ncid,'lon',rcode)
      jid    = ncdid (ncid,'lat',rcode)
      kid    = ncdid (ncid,'depth',rcode)
      tid    = ncdid (ncid,'time',rcode)
      bubid  = ncdid (ncid,'bubidx',rcode)

      call ncdinq (ncid,iid,name,jpi,rcode)
      call ncdinq (ncid,jid,name,jpj,rcode)
      call ncdinq (ncid,kid,name,jpk,rcode)
      call ncdinq (ncid,tid,name,jpt,rcode)
      call ncdinq (ncid,bubid,name,jpbub,rcode)

      cid  = ncvid (ncid,'bubble',rcode)

! --- if bubble index does not already exist,
! --- add one bubble index to bubidx unlimited dimension

!      if ( (bubidx.le.0) .or. (bubidx.gt.jpbub) ) then
!
!        bubidx =  jpbub + 1
!
!      endif

      if (bubidx.le.0) bubidx =  jpbub + 1

! --- put variable ---

      start (1) = 1
      start (2) = 1
      start (3) = 1
      start (4) = 1
      start (5) = bubidx

      count (1) = jpi
      count (2) = jpj
      count (3) = jpk
      count (4) = jpt
      count (5) = bubcount

      call ncvpt (ncid,cid,start,count,bub,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilcdfzon
