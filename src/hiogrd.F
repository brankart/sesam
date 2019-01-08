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
C ---                    HIOGRD.F                               ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-02  (L Parent)                          ---
C --- modification : 99-05  (C.E. Testut)                       ---
C --- modification : 99-11  (J.M. Brankart)                     ---
C --- modification : 01-06  (C.E. Testut)                       ---
C --- modification : 03-04  (J.M. Brankart)                     ---
C --- modification : 08-02  (J.M. Brankart)                     ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  readgrd : Read SESAM grids
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiogrd
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC readgrd, readlev, readtime

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readgrd (kflagxy,kindxy)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read SESAM grids (in BIMG, DIMG or CDF format)
CCC  -------
CCC  Method : Read grid arrays according to grid type and grid file format
CCC  ------
CCC  Input :  kegrd  : grid file format (1=BIMG, 2=DIMG, 3=CDF, 4=NC)
CCC  -----    kngrd  : grid type 1: lon=lon0+i*dx, lat=lat0+j*dx
CCC                              2: lon=lon(i),    lat=lat(j)
CCC                              3: lon=lon(i,j),  lat=lat(i,j)
CCC                              4: ORCA2
CCC           kjpi   : grid size (1st dimension)
CCC           kjpj   : grid size (2nd dimension)
CCC
CCC  Output : grid arrays in 'mod_coord.F'
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_coord
      use netcdf
      use utilcdfvar
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: kflagxy, kindxy
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80) :: cdrec1,cdrec2,cdrec3,cdrec4
      CHARACTER(len=word80), dimension(1:4) :: dimunit
      CHARACTER(len=word80) :: unit, kform
      CHARACTER(len=4) :: verdimg
      CHARACTER(len=bgword) :: kfname, xdim, ydim, xnam, ynam
      INTEGER :: kngrd, kegrd, kjpi, kjpj
      INTEGER :: nx, ny, nz, nt, ndim, icod
      BIGREAL4 :: lonmin, latmin, dx, dy, spval
      INTEGER :: allocok, irec
      BIGREAL4, allocatable, dimension(:,:) :: ptabij,ptablonij,ptablatij
      BIGREAL4, allocatable, dimension(:) :: ptabi,ptabj,ptabk
      INTEGER :: ji,jj,jk,jt,jdim
      INTEGER :: ierr, idf, idx, idy, idd, idxv, idyv, ndimx, ndimy
      INTEGER, allocatable, dimension(:) :: idimx,idimy,vstart,vcount,vstrid
CC----------------------------------------------------------------------
C
      SELECT CASE(kflagxy)
      CASE(1)
        kfname=varfgrd(kindxy)
        kngrd=varngrd(kindxy)
        kegrd=varegrd(kindxy)
        kjpi=var_jpi(kindxy)
        kjpj=var_jpj(kindxy)
        xdim=varxdim(kindxy)
        ydim=varydim(kindxy)
        xnam=varxnam(kindxy)
        ynam=varynam(kindxy)
      CASE(2)
        kfname=dtafgrd(kindxy)
        kngrd=dtangrd(kindxy)
        kegrd=dtaegrd(kindxy)
        kjpi=dta_jpi(kindxy)
        kjpj=dta_jpj(kindxy)
        xdim=dtaxdim(kindxy)
        ydim=dtaydim(kindxy)
        xnam=dtaxnam(kindxy)
        ynam=dtaynam(kindxy)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : ../readgrd'
         WRITE(numout,*) '    ==> READING grid file ',kfname(1:lenv(kfname))
      ENDIF
C
C -1.- Open grid file
C -------------------
C
      SELECT CASE (kegrd)
      CASE (1)
C --- BIMG format
         CALL openfile(numfil,kfname,kform=clunf)
      CASE (2)
C --- DIMG format
C Open file with incorrect record lentgh
         irecl  = ibloc*((kjpi*kjpj)/ibloc+1)*jpbyt4
         CALL openfile(numfil,kfname,clold,clunf,cldir,irecl)
C Read correct record lentgh at the beginning of the file
         READ(UNIT=numfil,REC=1,ERR=102,IOSTAT=iost) verdimg,cdrec1,irecl
         IF (verdimg.NE.'@!01') GOTO 104
         CLOSE (UNIT=numfil)
C Open file with correct record lentgh
         CALL openfile(numfil,kfname,clunk,clunf,cldir,irecl)
      CASE (3)
C --- CDF format
      CASE (4)
C --- NC format
         ierr = NF90_OPEN(kfname,NF90_NOWRITE,idf)
         IF (ierr.NE.0) GOTO 105
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -2.- Read file header
C ---------------------
C
      ndim=2
      SELECT CASE (kegrd)
      CASE (1)
C --- BIMG format
         READ(UNIT=numfil,ERR=1000,IOSTAT=iost) cdrec1
         READ(UNIT=numfil,ERR=1000,IOSTAT=iost) cdrec2
         READ(UNIT=numfil,ERR=1000,IOSTAT=iost) cdrec3
         READ(UNIT=numfil,ERR=1000,IOSTAT=iost) cdrec4
         READ(UNIT=numfil,ERR=1000,IOSTAT=iost) nx, ny, nz, nt, ndim, icod
         READ(UNIT=numfil,ERR=1000,IOSTAT=iost) lonmin,latmin,dx,dy,spval
         READ(UNIT=numfil,ERR=1000,IOSTAT=iost)
         READ(UNIT=numfil,ERR=1000,IOSTAT=iost)
      CASE (2)
C --- DIMG format
         READ(UNIT=numfil,REC=1,ERR=102,IOSTAT=iost) verdimg,cdrec1,irecl,
     $        nx, ny, nz,nt,ndim,lonmin,latmin,dx,dy,spval
      CASE (3)
C --- CDF format
         CALL cdfrdim(kfname,nx,ny,nz,nt,cdrec1)
      CASE (4)
C --- NC format
         nz=1 ; nt=1
C
         ierr = NF90_INQ_DIMID(idf,xdim,idx)
         IF (ierr.NE.0) GOTO 106
         ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=nx)
         IF (ierr.NE.0) GOTO 106
C
         ierr = NF90_INQ_DIMID(idf,ydim,idy)
         IF (ierr.NE.0) GOTO 106
         ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=ny)
         IF (ierr.NE.0) GOTO 106
C
         ierr = NF90_INQ_VARID(idf,xnam,idxv)
         IF (ierr.NE.0) GOTO 107
         ierr = NF90_INQ_VARID(idf,ynam,idyv)
         IF (ierr.NE.0) GOTO 107
C
         ierr = NF90_INQUIRE_VARIABLE(idf,idxv,ndims=ndimx)
         IF (ierr.NE.0) GOTO 107
         ierr = NF90_INQUIRE_VARIABLE(idf,idyv,ndims=ndimy)
         IF (ierr.NE.0) GOTO 107
         IF (ndimx.NE.ndimy) GOTO 107
C
         allocate(idimx(ndimx),idimy(ndimy))
         allocate(vstart(ndimx),vcount(ndimx),vstrid(ndimx))
         vstart(1:ndimx)=1 ; vcount(1:ndimx)=1 ; vstrid(1:ndimx)=1
C
         ierr = NF90_INQUIRE_VARIABLE(idf,idxv,dimids=idimx)
         IF (ierr.NE.0) GOTO 107
         ierr = NF90_INQUIRE_VARIABLE(idf,idyv,dimids=idimy)
         IF (ierr.NE.0) GOTO 107
C
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Check grid dimensions
      IF ((kjpi.NE.nx).OR.(kjpj.NE.ny).OR.(ndim.NE.2)) GOTO 103
C
C -3.- Read grid in file
C ----------------------
C
      SELECT CASE (kegrd)
      CASE (1,2)
C --- BIMG, DIMG format
C
         allocate ( ptabij(1:nx,1:ny), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptabij(:,:) = FREAL4(0.0)
C
         DO jdim=1,ndim
C	 
            SELECT CASE (kegrd)
            CASE (1)
               READ(UNIT=numfil,ERR=1000,IOSTAT=iost) ptabij(:,:)
            CASE (2)
               irec=2+(jt-1)*nz*ndim+(jk-1)*ndim+(jdim-1)
               READ(UNIT=numfil,REC=irec,ERR=1000,IOSTAT=iost) ptabij(:,:)
            CASE DEFAULT
               GOTO 1000
            END SELECT
C
            IF (jdim.EQ.1) THEN
               IF (kngrd.LE.2) THEN
                  longi(:) = FREAL(ptabij(:,1))
               ELSE
                  gridij(:,:)%longi = FREAL(ptabij(:,:))
               ENDIF
            ELSE
               IF (kngrd.LE.2) THEN
                  latj(:) = FREAL(ptabij(1,:))
               ELSE
                  gridij(:,:)%latj = FREAL(ptabij(:,:))
               ENDIF
            ENDIF
C
         ENDDO
C
      CASE (3)
C --- CDF format
C
         SELECT CASE (kngrd)
         CASE (1,2)
C Regular grid
            allocate ( ptabi(1:nx), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptabi(:) = FREAL4(0.0)
            allocate ( ptabj(1:ny), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptabj(:) = FREAL4(0.0)
            allocate ( ptabk(1:nz), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptabk(:) = FREAL4(0.0)
C
            CALL cdfrloc(kfname,ptabi,ptabj,ptabk,dimunit)
            longi(:) = FREAL(ptabi(:))
            latj(:) = FREAL(ptabj(:))
C
         CASE (3,4)
C Irregular grid
            allocate ( ptablonij(1:nx,1:ny), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptablonij(:,:) = FREAL4(0.0)
C
            allocate ( ptablatij(1:nx,1:ny), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptablatij(:,:) = FREAL4(0.0)
C
            CALL cdfrpos(kfname,ptablonij,ptablatij,spval,unit)
            gridij(:,:)%longi = FREAL(ptablonij(:,:))
            gridij(:,:)%latj = FREAL(ptablatij(:,:))
         CASE DEFAULT
            GOTO 1000
         END SELECT
C
      CASE (4)
C --- NC format
C
         SELECT CASE (kngrd)
         CASE (1,2)
C Regular grid
           allocate ( ptabi(1:nx), stat=allocok )
           IF (allocok.NE.0) GOTO 1001
           ptabi(:) = FREAL4(0.0)
           allocate ( ptabj(1:ny), stat=allocok )
           IF (allocok.NE.0) GOTO 1001
           ptabj(:) = FREAL4(0.0)
C
           vcount(1:ndimx)=1
           DO idd = 1,ndimx
             IF (idimx(idd).EQ.idx) vcount(idd)=kjpi
           ENDDO
           IF (ALL(vcount(1:ndimx).EQ.1)) GOTO 108
C
           ierr = NF90_GET_VAR(idf,idxv,ptabi,start=vstart,
     $                         count=vcount,stride=vstrid)
           IF (ierr.NE.0) GOTO 109
C
           vcount(1:ndimy)=1 
           DO idd = 1,ndimy
             IF (idimy(idd).EQ.idy) vcount(idd)=kjpj
           ENDDO
           IF (ALL(vcount(1:ndimy).EQ.1)) GOTO 108
C
           ierr = NF90_GET_VAR(idf,idyv,ptabj,start=vstart,
     $                         count=vcount,stride=vstrid)
           IF (ierr.NE.0) GOTO 109
C
           longi(:) = FREAL(ptabi(:))
           latj(:) = FREAL(ptabj(:))
C
         CASE (3,4)
C Irregular grid
           allocate ( ptabij(1:nx,1:ny), stat=allocok )
           IF (allocok.NE.0) GOTO 1001
           ptabij(:,:) = FREAL4(0.0)
C
           vcount(1:ndimx)=1
           DO idd = 1,ndimx
             IF (idimx(idd).EQ.idx) vcount(idd)=kjpi
             IF (idimx(idd).EQ.idy) vcount(idd)=kjpj
           ENDDO
           IF (ALL(vcount(1:ndimx).EQ.1)) GOTO 108
C
           ierr = NF90_GET_VAR(idf,idxv,ptabij,start=vstart,
     $                         count=vcount,stride=vstrid)
           IF (ierr.NE.0) GOTO 109
           gridij(:,:)%longi = FREAL(ptabij(:,:))
C
           ierr = NF90_GET_VAR(idf,idyv,ptabij,start=vstart,
     $                         count=vcount,stride=vstrid)
           IF (ierr.NE.0) GOTO 109
           gridij(:,:)%latj = FREAL(ptabij(:,:))
C
         CASE DEFAULT
            GOTO 1000
         END SELECT
C
         IF (kngrd.EQ.4) THEN
           DO ji=1,nx
           DO jj=1,ny
             DO WHILE (gridij(ji,jj)%longi.GE.FREAL4(442.))
               gridij(ji,jj)%longi = gridij(ji,jj)%longi - FREAL4(360)
             ENDDO
             DO WHILE (gridij(ji,jj)%longi.LT.FREAL4(78.))
               gridij(ji,jj)%longi = gridij(ji,jj)%longi + FREAL4(360)
             ENDDO
           ENDDO
           ENDDO
         ENDIF
C
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -3.- Close grid file
C --------------------
C
      SELECT CASE (kegrd)
      CASE (1,2)
         CLOSE (UNIT=numfil)
      CASE (3)
      CASE (4)
         ierr = NF90_CLOSE(idf)
         IF (ierr.NE.0) GOTO 110
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -4.- Control Print
C ------------------
C
      IF (nprint.GE.3) THEN 
C
         SELECT CASE (kegrd)
         CASE (1)
            kform='(8x,2a)'
            WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
            WRITE(numout,kform) '- Title2: ',cdrec2(1:lenv(cdrec2))
            WRITE(numout,kform) '- Title3: ',cdrec3(1:lenv(cdrec3))
            WRITE(numout,kform) '- Title4: ',cdrec4(1:lenv(cdrec4))
            kform='(8x,a,5i5)'
            WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
            kform='(8x,a,e12.3)'
            WRITE(numout,kform) '- Special value: ',spval
            kform='(8x,a,4e12.3)'
            WRITE(numout,kform) '- Regular grid: ',lonmin,latmin,dx,dy
         CASE (2)
            kform='(8x,2a)'
            WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
            kform='(8x,a,5i5)'
            WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
            kform='(8x,a,e12.3)'
            WRITE(numout,kform) '- Special value: ',spval
            kform='(8x,a,4e12.3)'
            WRITE(numout,kform) '- Regular grid: ',lonmin,latmin,dx,dy
         CASE (3)
            kform='(8x,2a)'
            WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
            kform='(8x,a,4i5)'
            WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt
            kform='(8x,a,e12.3)'
            WRITE(numout,kform) '- Special value: ',spval
         CASE (4)
            kform='(8x,a,4i5)'
            WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt
         CASE DEFAULT
            GOTO 1000
         END SELECT
C
      ENDIF
C
      IF (nprint.GE.4) THEN 
         IF (kngrd.LE.2) THEN
            WRITE(numout,*) ' - longitude regular grid:'
            WRITE(numout,*) longi(:)
            WRITE(numout,*) ' - latitude regular grid:'
            WRITE(numout,*) latj(:)
         ELSE
            WRITE(numout,*) ' - sample of longitude irregular grid:'
            DO jj = 1,kjpj,MAX(1,kjpj/5)
               WRITE(numout,*) gridij(1:kjpi:MAX(1,kjpi/5),jj)%longi
            ENDDO
            WRITE(numout,*) ' - sample of latitude irregular grid:'
            DO jj = 1,kjpj,MAX(1,kjpj/5)
               WRITE(numout,*) gridij(1:kjpi:MAX(1,kjpi/5),jj)%latj
            ENDDO
         ENDIF
      ENDIF
C
C --- deallocate arrays
      IF (allocated(ptabij)) deallocate(ptabij)
      IF (allocated(ptablonij)) deallocate(ptablonij)
      IF (allocated(ptablatij)) deallocate(ptablatij)
      IF (allocated(ptabi)) deallocate(ptabi)
      IF (allocated(ptabj)) deallocate(ptabj)
      IF (allocated(ptabk)) deallocate(ptabk)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiogrd','readgrd')
 1001 CALL printerror2(0,1001,3,'hiogrd','readgrd')
C
 102  WRITE (texterror,*) 'Error reading grid file, iost=',iost
      CALL printerror2(0,102,3,'hiogrd','readgrd',comment=texterror)
 103  WRITE (texterror,*) 'bad dimensions in grid file'
      CALL printerror2(0,103,3,'hiogrd','readgrd',comment=texterror)
 104  WRITE (texterror,*) 'bad DIMG file version'
      CALL printerror2(0,103,3,'hiogrd','readgrd',comment=texterror)
 105  WRITE (texterror,*) 'Bad NetCDF file: ',kfname
      CALL printerror2(0,105,3,'hiogrd','readgrd',comment=texterror)
 106  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',kfname
      CALL printerror2(0,106,3,'hiogrd','readgrd',comment=texterror)
 107  WRITE (texterror,*) 'Bad variable in NetCDF file: ',kfname
      CALL printerror2(0,107,3,'hiogrd','readgrd',comment=texterror)
 108  WRITE (texterror,*) 'Bad dims for NetCDF variable in: ',kfname
      CALL printerror2(0,108,3,'hiogrd','readgrd',comment=texterror)
 109  WRITE (texterror,*) 'Error reading NetCDF file: ',kfname
      CALL printerror2(0,109,3,'hiogrd','readgrd',comment=texterror)
 110  WRITE (texterror,*) 'Error closing Netcdf file: ',kfname
      CALL printerror2(0,110,3,'hiogrd','readgrd',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readlev (kflagxy,kindxy)
CCC---------------------------------------------------------------------
CCC   
CCC  Purpose : Read SESAM vertical coordinate (from every file format)
CCC  -------
CCC  Method : Read grid arrays according to grid type and grid file format
CCC  ------
CCC  Input :  kegrd  : grid file format (1=BIMG, 2=DIMG, 3=CDF, 4=NC)
CCC  -----    
CCC  Output : grid arrays in 'mod_coord.F'
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_coord
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: kflagxy, kindxy
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80), dimension(1:4) :: dimunit
      CHARACTER(len=bgword) :: kfname, zdim, znam
      CHARACTER(len=word80) :: title
      INTEGER :: kngrd, kegrd, kjpk
      INTEGER :: nx, ny, nz, nt
      INTEGER :: allocok, ierr, idf, idz, idzv, ndimz, idd
      INTEGER, allocatable, dimension(:) :: idimz,vstart,vcount,vstrid
      BIGREAL4, allocatable, dimension(:) :: ptabi,ptabj,ptabk
CC----------------------------------------------------------------------
C
      SELECT CASE(kflagxy)
      CASE(1)
        kfname=varfgrd(kindxy)
        kngrd=varngrd(kindxy)
        kjpk=var_jpk(kindxy)
        zdim=varzdim(kindxy)
        znam=varznam(kindxy)
      CASE(2)
        kfname=dtafgrd(kindxy)
        kngrd=dtangrd(kindxy)
        kjpk=dta_jpk(kindxy)
        zdim=dtazdim(kindxy)
        znam=dtaznam(kindxy)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : ../readlev'
         WRITE(numout,*) '    ==> READING level file ',kfname(1:lenv(kfname))
      ENDIF
C
C -1.- Open grid file
C -------------------
C
      SELECT CASE (kegrd)
      CASE (1)
C --- BIMG format
         GOTO 101
      CASE (2)
         GOTO 101
      CASE (3)
C --- CDF format
      CASE (4)
C --- NC format
         ierr = NF90_OPEN(kfname,NF90_NOWRITE,idf)
         IF (ierr.NE.0) GOTO 105
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -2.- Read file header
C ---------------------
C
      SELECT CASE (kegrd)
      CASE (3)
C --- CDF format
         CALL cdfrdim(kfname,nx,ny,nz,nt,title)
      CASE (4)
C --- NC format
         nx=1 ; ny=1 ; nt=1
C
         ierr = NF90_INQ_DIMID(idf,zdim,idz)
         IF (ierr.NE.0) GOTO 106
         ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=nz)
         IF (ierr.NE.0) GOTO 106
C
         ierr = NF90_INQ_VARID(idf,znam,idzv)
         IF (ierr.NE.0) GOTO 107
C
         ierr = NF90_INQUIRE_VARIABLE(idf,idzv,ndims=ndimz)
         IF (ierr.NE.0) GOTO 107
C
         allocate(idimz(ndimz))
         allocate(vstart(ndimz),vcount(ndimz),vstrid(ndimz))
         vstart(1:ndimz)=1 ; vcount(1:ndimz)=1 ; vstrid(1:ndimz)=1
C
         ierr = NF90_INQUIRE_VARIABLE(idf,idzv,dimids=idimz)
         IF (ierr.NE.0) GOTO 107
C
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Check grid dimensions
C
      IF (kjpk.NE.nz) GOTO 103
C
C -3.- Read grid in file
C ----------------------
C
      SELECT CASE (kegrd)
      CASE (3)
C --- CDF format
C
         allocate ( ptabi(1:nx), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptabi(:) = FREAL4(0.0)
         allocate ( ptabj(1:ny), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptabj(:) = FREAL4(0.0)
         allocate ( ptabk(1:nz), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptabk(:) = FREAL4(0.0)
C
         CALL cdfrloc(kfname,ptabi,ptabj,ptabk,dimunit)
         levk(:) = FREAL(ptabk(:))
C
      CASE (4)
C --- NC format
C
         allocate ( ptabk(1:nz), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         ptabk(:) = FREAL4(0.0)

         vcount(1:ndimz)=1
         DO idd = 1,ndimz
           IF (idimz(idd).EQ.idz) vcount(idd)=kjpk
         ENDDO
         IF (ALL(vcount(1:ndimz).EQ.1)) GOTO 108
C
         ierr = NF90_GET_VAR(idf,idzv,ptabk,start=vstart,
     $                       count=vcount,stride=vstrid)
         IF (ierr.NE.0) GOTO 109
C
         levk(:) = FREAL(ptabk(:))
C
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -3.- Close grid file
C --------------------
C
      SELECT CASE (kegrd)
      CASE (3)
      CASE (4)
         ierr = NF90_CLOSE(idf)
         IF (ierr.NE.0) GOTO 110
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiogrd','readgrd')
 1001 CALL printerror2(0,1001,3,'hiogrd','readgrd')
C
 101  WRITE (texterror,*) 'Unsupported file format'
      CALL printerror2(0,101,3,'hiogrd','readlev',comment=texterror)
 103  WRITE (texterror,*) 'bad dimensions in grid file'
      CALL printerror2(0,103,3,'hiogrd','readgrd',comment=texterror)
 105  WRITE (texterror,*) 'Bad NetCDF file: ',kfname
      CALL printerror2(0,105,3,'hiogrd','readgrd',comment=texterror)
 106  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',kfname
      CALL printerror2(0,106,3,'hiogrd','readgrd',comment=texterror)
 107  WRITE (texterror,*) 'Bad variable in NetCDF file: ',kfname
      CALL printerror2(0,107,3,'hiogrd','readgrd',comment=texterror)
 108  WRITE (texterror,*) 'Bad dims for NetCDF variable in: ',kfname
      CALL printerror2(0,108,3,'hiogrd','readgrd',comment=texterror)
 109  WRITE (texterror,*) 'Error reading NetCDF file: ',kfname
      CALL printerror2(0,109,3,'hiogrd','readgrd',comment=texterror)
 110  WRITE (texterror,*) 'Error closing Netcdf file: ',kfname
      CALL printerror2(0,110,3,'hiogrd','readgrd',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readtime (kflagxy,kindxy)
CCC---------------------------------------------------------------------
CCC   
CCC  Purpose : Read SESAM time coordinate (from every file format)
CCC  -------
CCC  Method : Read grid arrays according to grid type and grid file format
CCC  ------
CCC  Input :  kegrd  : grid file format (1=BIMG, 2=DIMG, 3=CDF, 4=NC)
CCC  -----    
CCC  Output : grid arrays in 'mod_coord.F'
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_coord
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: kflagxy, kindxy
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=word80), dimension(1:4) :: dimunit
      CHARACTER(len=bgword) :: kfname, tdim, tnam
      CHARACTER(len=word80) :: title
      INTEGER :: kngrd, kegrd, kjpt
      INTEGER :: nx, ny, nz, nt
      INTEGER :: allocok, ierr, idf, idt, idtv, ndimt, idd
      INTEGER, allocatable, dimension(:) :: idimt,vstart,vcount,vstrid
      BIGREAL4, allocatable, dimension(:) :: ptabt
CC----------------------------------------------------------------------
C
      SELECT CASE(kflagxy)
      CASE(1)
        kfname=varfgrd(kindxy)
        kngrd=varngrd(kindxy)
        kjpt=var_jpt(kindxy)
        tdim=vartdim(kindxy)
        tnam=vartnam(kindxy)
      CASE(2)
        kfname=dtafgrd(kindxy)
        kngrd=dtangrd(kindxy)
        kjpt=dta_jpt(kindxy)
        tdim=dtatdim(kindxy)
        tnam=dtatnam(kindxy)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : ../readtime'
         WRITE(numout,*) '    ==> READING time file ',kfname(1:lenv(kfname))
      ENDIF
C
C -1.- Open grid file
C -------------------
C
      SELECT CASE (kegrd)
      CASE (1)
         GOTO 101
      CASE (2)
         GOTO 101
      CASE (3)
C --- CDF format
C
      CASE (4)
C --- NC format
         ierr = NF90_OPEN(kfname,NF90_NOWRITE,idf)
         IF (ierr.NE.0) GOTO 105
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -2.- Read file header
C ---------------------
C
      SELECT CASE (kegrd)
      CASE (3)
C --- CDF format
         CALL cdfrdim(kfname,nx,ny,nz,nt,title)
      CASE (4)
C --- NC format
         nx=1 ; ny=1 ; nz=1
C
         ierr = NF90_INQ_DIMID(idf,tdim,idt)
         IF (ierr.NE.0) GOTO 106
         ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=nt)
         IF (ierr.NE.0) GOTO 106
C
         ierr = NF90_INQ_VARID(idf,tnam,idtv)
         IF (ierr.NE.0) GOTO 107
C
         ierr = NF90_INQUIRE_VARIABLE(idf,idtv,ndims=ndimt)
         IF (ierr.NE.0) GOTO 107
C
         allocate(idimt(ndimt))
         allocate(vstart(ndimt),vcount(ndimt),vstrid(ndimt))
         vstart(1:ndimt)=1 ; vcount(1:ndimt)=1 ; vstrid(1:ndimt)=1
C
         ierr = NF90_INQUIRE_VARIABLE(idf,idtv,dimids=idimt)
         IF (ierr.NE.0) GOTO 107
C
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Check grid dimensions
C
      IF (kjpt.NE.nt) GOTO 103
C
C -3.- Read grid in file
C ----------------------
C
      allocate ( ptabt(1:nt), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabt(:) = FREAL4(0.0)
C
      SELECT CASE (kegrd)
      CASE (3)
C --- CDF format
C
         CALL cdfrtim(kfname,ptabt)
C
      CASE (4)
C --- NC format
C
         vcount(1:ndimt)=1
         DO idd = 1,ndimt
           IF (idimt(idd).EQ.idt) vcount(idd)=kjpt
         ENDDO
         IF (ALL(vcount(1:ndimt).EQ.1)) GOTO 108
C
         ierr = NF90_GET_VAR(idf,idtv,ptabt,start=vstart,
     $                       count=vcount,stride=vstrid)
         IF (ierr.NE.0) GOTO 109
C
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      time(:) = FREAL(ptabt(:))
C
C -3.- Close grid file
C --------------------
C
      SELECT CASE (kegrd)
      CASE (3)
      CASE (4)
         ierr = NF90_CLOSE(idf)
         IF (ierr.NE.0) GOTO 110
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiogrd','readtime')
 1001 CALL printerror2(0,1001,3,'hiogrd','readtime')
C
 101  WRITE (texterror,*) 'Unsupported file format'
      CALL printerror2(0,101,3,'hiogrd','readtime',comment=texterror)
 103  WRITE (texterror,*) 'bad dimensions in grid file'
      CALL printerror2(0,103,3,'hiogrd','readgrd',comment=texterror)
 105  WRITE (texterror,*) 'Bad NetCDF file: ',kfname
      CALL printerror2(0,105,3,'hiogrd','readgrd',comment=texterror)
 106  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',kfname
      CALL printerror2(0,106,3,'hiogrd','readgrd',comment=texterror)
 107  WRITE (texterror,*) 'Bad variable in NetCDF file: ',kfname
      CALL printerror2(0,107,3,'hiogrd','readgrd',comment=texterror)
 108  WRITE (texterror,*) 'Bad dims for NetCDF variable in: ',kfname
      CALL printerror2(0,108,3,'hiogrd','readgrd',comment=texterror)
 109  WRITE (texterror,*) 'Error reading NetCDF file: ',kfname
      CALL printerror2(0,109,3,'hiogrd','readgrd',comment=texterror)
 110  WRITE (texterror,*) 'Error closing Netcdf file: ',kfname
      CALL printerror2(0,110,3,'hiogrd','readgrd',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiogrd
