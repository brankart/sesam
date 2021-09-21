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
! ---                    UTILCDFPAK                             ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-10 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE cdfrdimpak : Read dimensions from CPAK file
! --- SUBROUTINE cdfrhdrpak : Read header from CPAK file
! --- SUBROUTINE cdfrpak    : Read segment of Vx vector from CPAK file
! --- SUBROUTINE cdfwdimpak : Write dimensions in CPAK file
! --- SUBROUTINE cdfwhdrpak : Write header in CPAK file
! --- SUBROUTINE cdfwpak    : Write segment of Vx vector in CPAK file
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilcdfpak
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC cdfrdimpak,cdfrhdrpak,cdfrpak
      PUBLIC cdfwdimpak,cdfwhdrpak,cdfwpak

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfrdimpak(filename,jpx,varend,namelength,title)

      implicit none
      include 'netcdf.inc' 

      integer jpx,varend,namelength
      character*(*) filename,title

      character*200 name
      integer iid,jid,kid
      integer ncid, rcode, ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'idx',rcode)
      jid  = ncdid (ncid,'var',rcode)
      kid  = ncdid (ncid,'namelength',rcode)

      call ncdinq (ncid,iid,name,jpx,rcode)
      call ncdinq (ncid,jid,name,varend,rcode)
      call ncdinq (ncid,kid,name,namelength,rcode)

! --- get global attribute ---

      ll=len(title)
      call ncagtc (ncid,ncglobal,'title',title,ll,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfrhdrpak(filename,nam,dim,nbr,moy,ect)
      use netcdf

      implicit none
      include 'netcdf.inc' 

      character*(*) filename
      BIGREAL4 moy(*),ect(*)
      integer dim(*),nbr(*)
      character*(*) nam(*)

      character*200 dimname
      integer iid,kid,moyid,ectid,dimid,nbrid,namid
      integer varend,namelength,namell
      integer ncid, rcode,ll
      integer start(2),count(2)

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      iid  = ncdid (ncid,'var',rcode)
      kid  = ncdid (ncid,'namelength',rcode)
      call ncdinq (ncid,iid,dimname,varend,rcode)
      call ncdinq (ncid,kid,dimname,namelength,rcode)

      moyid  = ncvid (ncid,'mean',rcode)
      ectid  = ncvid (ncid,'std',rcode)
      dimid  = ncvid (ncid,'dim',rcode)
      nbrid  = ncvid (ncid,'nbr',rcode)
      namid  = ncvid (ncid,'name',rcode)

! --- get variable ---

      start(1) = 1
      count(1) = varend
      call ncvgt  (ncid,moyid,start,count,moy,rcode)
      call ncvgt  (ncid,ectid,start,count,ect,rcode)
      call ncvgt  (ncid,dimid,start,count,dim,rcode)
      call ncvgt  (ncid,nbrid,start,count,nbr,rcode)

      namell = len(nam)
      start(1) = 1
      count(1) = namell
      start(2) = 1
      count(2) = varend
      !call ncvgtc (ncid,namid,start,count,nam(1:varend),namell*varend,rcode)
      rcode = NF90_GET_VAR(ncid,namid,nam(1:varend),start=start,count=count)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfrpak(filename,x,idx0,idx1)

      implicit none
      include 'netcdf.inc' 

      BIGREAL4 x(*)
      integer idx0,idx1
      character*(*) filename

      integer iid,cid
      integer ncid, rcode,ll
      integer start(1),count(1)

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncnowrit,rcode)

! --- get information on NetCDF file ---

      cid  = ncvid (ncid,'state',rcode)

! --- get variable ---

      start(1) = idx0
      count(1) = idx1 - idx0 + 1
      call ncvgt  (ncid,cid,start,count,x,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfwdimpak(filename,jpx,varend,namelength,title)

      implicit none
      include 'netcdf.inc' 

      integer jpx,varend,namelength
      character*(*) filename,title

      integer moyid,ectid,dimid,nbrid,namid
      integer iid,jid,kid,cdim(2),cid
      integer ncid, rcode, ll

! --- create NetCDF file ---

      ll=lenv(filename)
      ncid = nccre (filename(1:ll),ncnoclob,rcode)

! --- create dimensions ---

      iid    = ncddef (ncid,'idx',jpx,rcode)
      jid    = ncddef (ncid,'var',varend,rcode)
      kid    = ncddef (ncid,'namelength',namelength,rcode)

! --- create variable ---

      cdim (1) = kid
      cdim (2) = jid
      namid = ncvdef (ncid,'name',ncchar,2,cdim,rcode)

      cdim (1) = jid
      dimid = ncvdef (ncid,'dim',nclong,1,cdim,rcode)
      nbrid = ncvdef (ncid,'nbr',nclong,1,cdim,rcode)
      moyid = ncvdef (ncid,'mean',ncfloat,1,cdim,rcode)
      ectid = ncvdef (ncid,'std',ncfloat,1,cdim,rcode)

      cdim (1) = iid
      cid = ncvdef (ncid,'state',ncfloat,1,cdim,rcode)

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
      subroutine cdfwhdrpak(filename,nam,dim,nbr,moy,ect)

      implicit none
      include 'netcdf.inc' 

      character*(*) filename
      BIGREAL4 moy(*),ect(*)
      integer dim(*),nbr(*)
      character*(*) nam(*)

      character*200 dimname
      integer start(2),count(2)
      integer jid,kid,moyid,ectid,dimid,nbrid,namid
      integer varend,namelength,namell
      integer ncid, rcode, ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- get information on NetCDF file ---

      jid  = ncdid (ncid,'var',rcode)
      kid  = ncdid (ncid,'namelength',rcode)

      call ncdinq (ncid,jid,dimname,varend,rcode)
      call ncdinq (ncid,kid,dimname,namelength,rcode)

      moyid  = ncvid (ncid,'mean',rcode)
      ectid  = ncvid (ncid,'std',rcode)
      dimid  = ncvid (ncid,'dim',rcode)
      nbrid  = ncvid (ncid,'nbr',rcode)
      namid  = ncvid (ncid,'name',rcode)

! --- put variables ---

      start (1) = 1
      count (1) = varend
      call ncvpt (ncid,moyid,start,count,moy,rcode)
      call ncvpt (ncid,ectid,start,count,ect,rcode)
      call ncvpt (ncid,dimid,start,count,dim,rcode)
      call ncvpt (ncid,nbrid,start,count,nbr,rcode)

      namell = len(nam)
      start (1) = 1
      count (1) = namell
      start (2) = 1
      count (2) = varend
      call ncvptc (ncid,namid,start,count,nam,namell*varend,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine cdfwpak(filename,x,idx0,idx1)

      implicit none
      include 'netcdf.inc' 

      character*(*) filename
      BIGREAL4 x(*)
      integer idx0,idx1

      integer cdim(1),start(1),count(1)
      integer iid,cid
      integer ncid, rcode, ll

! --- open NetCDF file ---

      ll=lenv(filename)
      ncid = ncopn (filename(1:ll),ncwrite,rcode)

! --- get information on NetCDF file ---

      cid  = ncvid (ncid,'state',rcode)

! --- put variable ---

      start (1) = idx0
      count (1) = idx1 - idx0 + 1
      call ncvpt (ncid,cid,start,count,x,rcode)

! --- close NetCDF file ---

      call ncclos (ncid,rcode)

      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilcdfpak
