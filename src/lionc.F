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
C ---                  LIONC.F                                  ---
C ---                                                           ---
C --- original     : 07-11 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- Routines to read and write from generic NetCDF files
C ---
C --- SUBROUTINE  readnc       : Read variable field from NetCDF file
C --- SUBROUTINE  writenc      : Write variable field in NetCDF file
C --- SUBROUTINE  evalhdrmsknc : Get variable dimensions from mask files
C --- SUBROUTINE  readmsknc    : Read mask arrays from NetCDF mask files
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE lionc
      use mod_main
      use hiogrd
      use utilvct
      IMPLICIT NONE
      PRIVATE

      PUBLIC readnc,writenc,evalhdrmsknc,readmsknc

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readnc(kfname,kvectsout,kjsxy,ksomsxynbr,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read variable field (of Vx or Vy object) from NetCDF file
CCC  -------
CCC  Method : Open, read and close NetCDF file (.nc or .ncdta)
CCC  ------
CCC  Input :  kfname     : Filename
CCC  -----    kjsxy      : Index of variable field to read
CCC           kflagxyo   : vector type (Vx,Vy)
CCC  Output : kvectsout  : variable field 1D vector (jpx or jpy)
CCC  ------   ksomsxynbr : number of values loaded
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfname
      BIGREAL, dimension(:), intent(out) :: kvectsout
      INTEGER, intent(in) :: kjsxy,kflagxyo
      INTEGER, intent(out) :: ksomsxynbr
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpssize
      BIGREAL4, allocatable, dimension(:,:) :: ptabij
      INTEGER :: indsxy,js,sompartsxynbr,jextdta,jextvar
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt
      CHARACTER(len=bgword) :: sxy_nam, kform
      CHARACTER(len=bgword), dimension(1:4) :: dimunit
C
      INTEGER, allocatable, dimension(:) :: idims,vstart,vcount,vstrid
      INTEGER :: ierr, idf, idx, idy, idz, idt, idv, idd
      INTEGER :: idvx, idvy, idvz, idvt
      CHARACTER(len=bgword) :: xdim, ydim, zdim, tdim
      INTEGER :: nx, ny, nz, nt, ndims
      INTEGER :: ji, jj, jk, jt
      INTEGER :: jpiend, jpjend, jpkend, jptend
      LOGICAL :: varfile
CC----------------------------------------------------------------------
C Size of output vector
      jpssize=size(kvectsout,1)
C
C -1.- Get variable characteristics from SESAM configuration
C ----------------------------------------------------------
C
C Get file type (Vx or Vy) and variable index
      SELECT CASE(kflagxyo)
      CASE(1)
C --- var
         indsxy  = var_ord(kjsxy)
         varfile = .TRUE.
      CASE(2)
C --- dta
         indsxy  = dta_ord(kjsxy)
         jextvar = indext(kfname,extvartab,nbextvar)
         jextdta = indext(kfname,extdtatab,nbextdta)
         IF (jextdta.EQ.4) THEN
            varfile=.FALSE.
         ELSEIF (jextvar.EQ.4) THEN
            varfile=.TRUE.
         ELSE
            GOTO 1000
         ENDIF
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Set file dimensions
      sxy_jpi = 1
      sxy_jpj = 1
      sxy_jpk = 1
      sxy_jpt = 1
      IF (varfile) THEN
C --- read Vx or Vy variable field from var file
         sxy_nam = varifil(indsxy)
         sxy_dim = var_dim(indsxy)
         IF (sxy_dim.GE.1) sxy_jpi = var_jpi(indsxy)
         IF (sxy_dim.GE.2) sxy_jpj = var_jpj(indsxy)
         IF (sxy_dim.GE.3) sxy_jpk = var_jpk(indsxy)
         IF (sxy_dim.GE.4) sxy_jpt = var_jpt(indsxy)
         xdim = varxdim(indsxy)
         ydim = varydim(indsxy)
         zdim = varzdim(indsxy)
         tdim = vartdim(indsxy)
      ELSE
C --- read Vy variable field from dta file
         sxy_nam = dtaifil(indsxy)
         sxy_dim = dta_dim(indsxy)
         IF (sxy_dim.GE.1) sxy_jpi = dta_jpi(indsxy)
         IF (sxy_dim.GE.2) sxy_jpj = dta_jpj(indsxy)
         IF (sxy_dim.GE.3) sxy_jpk = dta_jpk(indsxy)
         IF (sxy_dim.GE.4) sxy_jpt = dta_jpt(indsxy)
         xdim = dtaxdim(indsxy)
         ydim = dtaydim(indsxy)
         zdim = dtazdim(indsxy)
         tdim = dtatdim(indsxy)
      ENDIF
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readfil/readnc'
         WRITE(numout,*) '    ==> READING file ',kfname(1:lenv(kfname))
         WRITE(numout,*) '    ==> READING variable: ',sxy_nam
      ENDIF
C
C Set size of array to read
      jpiend = 1
      jpjend = 1
      jpkend = 1
      jptend = 1
      SELECT CASE(kflagxyo)
      CASE(1)
         IF (sxy_dim.GE.1) jpiend = var_jpi(indsxy)
         IF (sxy_dim.GE.2) jpjend = var_jpj(indsxy)
         IF (sxy_dim.GE.3) jpkend = var_jpk(indsxy)
         IF (sxy_dim.GE.4) jptend = var_jpt(indsxy)
      CASE(2)
         IF (sxy_dim.GE.1) jpiend = dta_jpi(indsxy)
         IF (sxy_dim.GE.2) jpjend = dta_jpj(indsxy)
         IF (sxy_dim.GE.3) jpkend = dta_jpk(indsxy)
         IF (sxy_dim.GE.4) jptend = dta_jpt(indsxy)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -2.- Open NetCDF file, get and check dimensions
C -----------------------------------------------
C
      ierr = NF90_OPEN(kfname,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 101
C
      IF (sxy_dim.GE.1) THEN
        ierr = NF90_INQ_DIMID(idf,xdim,idx)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=nx)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpi.NE.nx) GOTO 102
      ENDIF
C
      IF (sxy_dim.GE.2) THEN
        ierr = NF90_INQ_DIMID(idf,ydim,idy)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=ny)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpj.NE.ny) GOTO 102
      ENDIF
C
      IF (sxy_dim.GE.3) THEN
        ierr = NF90_INQ_DIMID(idf,zdim,idz)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=nz)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpk.NE.nz) GOTO 102
      ENDIF
C
      IF (sxy_dim.GE.4) THEN
        ierr = NF90_INQ_DIMID(idf,tdim,idt)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=nt)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpt.NE.nt) GOTO 102
      ENDIF
C
C -3.- Get variable id, get the set of dimensions spanned by the variable,
C      Get indices of x, y, z, t dimensions in that set of dimensions
C      Initialize start, count and strid array for reading Netcdf variable
C ------------------------------------------------------------------------
C
      ierr = NF90_INQ_VARID(idf,sxy_nam,idv)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
      IF (ierr.NE.0) GOTO 103
      IF (sxy_dim.GT.ndims) GOTO 103
C
      allocate(idims(ndims),vstart(ndims),vcount(ndims),vstrid(ndims))
      vstart(1:ndims)=1 ; vcount(1:ndims)=1 ; vstrid(1:ndims)=1
C 
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
      IF (ierr.NE.0) GOTO 103
C
      idvx = 0 ; idvy = 0 ; idvz = 0 ; idvt = 0
      DO idd = 1,ndims
        IF (idims(idd).EQ.idx) idvx = idd
        IF (idims(idd).EQ.idy) idvy = idd
        IF (idims(idd).EQ.idz) idvz = idd
        IF (idims(idd).EQ.idt) idvt = idd
      ENDDO
C
      IF ( (sxy_dim.GE.1).AND.(idvx.EQ.0) ) GOTO 105
      IF ( (sxy_dim.GE.2).AND.(idvy.EQ.0) ) GOTO 105
      IF ( (sxy_dim.GE.3).AND.(idvz.EQ.0) ) GOTO 105
      IF ( (sxy_dim.GE.4).AND.(idvt.EQ.0) ) GOTO 105
C
      DO idd = 1,ndims
        IF (idd.EQ.idvx) vcount(idd)=sxy_jpi
        IF (idd.EQ.idvy) vcount(idd)=sxy_jpj
      ENDDO
C
C -4.- Read 2D arrays in NetCDF file
C ----------------------------------
C
C --- allocation ptabij
      allocate ( ptabij(1:sxy_jpi,1:sxy_jpj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL4(0.0)
C
      js=1
      DO jt=1,jptend
      DO jk=1,jpkend
C
         IF (idvz.NE.0) vstart(idvz)=jk
         IF (idvt.NE.0) vstart(idvt)=jt
C
         ierr = NF90_GET_VAR(idf,idv,ptabij,start=vstart,count=vcount,
     $                       stride=vstrid)
         IF (ierr.NE.0) GOTO 104
C
         CALL mk4vct(kvectsout(js:),ptabij(:,:),
     $               jk,jt,kjsxy,sompartsxynbr,kflagxyo)
         js  = js + sompartsxynbr
C
      ENDDO
      ENDDO
C
      ksomsxynbr=js-1
C
      IF (allocated(ptabij)) deallocate(ptabij)
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 106
C
C -5.- Control Print
C ------------------
C
      IF (nprint.GE.3) THEN
         kform='(8x,a,4i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt
         IF (  (jpkend.NE.nz).OR.(jptend.NE.nt) ) THEN
            kform='(8x,a,4i5)'
            WRITE(numout,kform) '- Size of loaded array : ',jpiend,
     $                                  jpjend,jpkend,jptend
         ENDIF
         kform='(8x,a,i9)'
         WRITE(numout,kform) '- Number of loaded values : ',ksomsxynbr
      ENDIF
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'lionc','readnc')
 1001 CALL printerror2(0,1001,3,'lionc','readnc')
C
 101  WRITE (texterror,*) 'Bad NetCDF file: ',kfname
      CALL printerror2(0,101,3,'lionc','readnc',comment=texterror)
 102  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',kfname
      CALL printerror2(0,102,3,'lionc','readnc',comment=texterror)
 103  WRITE (texterror,*) 'Bad variable in NetCDF file: ',sxy_nam
      CALL printerror2(0,103,3,'lionc','readnc',comment=texterror)
 104  WRITE (texterror,*) 'Error reading NetCDF variable: ',sxy_nam
      CALL printerror2(0,104,3,'lionc','readnc',comment=texterror)
 105  WRITE (texterror,*) 'Bad dimensions for NetCDF variable: ',sxy_nam
      CALL printerror2(0,105,3,'lionc','readnc',comment=texterror)
 106  WRITE (texterror,*) 'Error closing Netcdf file: ',kfname
      CALL printerror2(0,106,3,'lionc','readnc',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writenc(kfnoutnc,kvectsin,
     $           ktxy_indvct,ktxy_indmsk,txyend,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write variable field (of Vx or Vy object) in NetCDF file
CCC  -------
CCC  Method : Open, write and close NetCDF file (.nc or .ncdta)
CCC  ------
CCC  Input : kfnoutnc    : Filename
CCC  -----   kvectsin    : 1D SESAM vector object (Vx or Vy)
CCC          ktxy_indvct : array with variable field beginning indices
CCC                                  (in vector object Vx or Vy)
CCC          ktxy_indmsk : table with variable field indices
CCC                                  (in variable list)
CCC          txyend      : number of variable fields to write
CCC          kflagxyo    : vector type (Vx,Vy)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_coord
      use mod_spacexyo , only : spvaldta,spvalvar
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutnc
      BIGREAL, dimension(:), intent(in) :: kvectsin
      INTEGER, dimension(:), intent(in) :: ktxy_indvct,ktxy_indmsk
      INTEGER, intent(in) :: txyend,kflagxyo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: nx, ny, nz, nt
      BIGREAL4 :: spval
      BIGREAL4, dimension(1:nbvar) :: sxyvmsk
      INTEGER  ::allocok
      BIGREAL4, allocatable, dimension(:,:) :: ptabij
      INTEGER :: sxyend, sxy_dim_max
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_dim,sxyemsk,
     $     sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxy_nbr,sxyngrd,sxyegrd
      INTEGER :: jsxy,indsxy,sompartxynbr,jtxy
      INTEGER, dimension(1:nbvar) :: txy_indtab
      CHARACTER(len=bgword), dimension(1:nbvar) :: sxy_nam
      CHARACTER(len=bgword), dimension(1:nbvar) :: sxyfgrd,sxyfmsk
      INTEGER :: ji, jj, jk, jt
      INTEGER :: jpiend, jpjend, jpkend, jptend
      LOGICAL :: filexists
      INTEGER, allocatable, dimension(:) :: idims,vstart,vcount,vstrid
      INTEGER :: ierr, idf, idx, idy, idz, idt, idv, idd
      INTEGER :: idvx, idvy, idvz, idvt, ndims
      CHARACTER(len=bgword) :: xdim, ydim, zdim, tdim
      CHARACTER(len=bgword) :: cpcom
CC----------------------------------------------------------------------
C
C -1.- Get variable characteristics from SESAM configuration
C ----------------------------------------------------------
C
      SELECT CASE(kflagxyo)
      CASE(1)
C --- var
         sxyend  = varend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=var_ord(jsxy)
            indsxy= sxy_ord(jsxy)
C
            sxy_nam(indsxy)=varofil(indsxy)
            sxy_dim(indsxy)=var_dim(indsxy)
            sxy_jpi(indsxy)=var_jpi(indsxy)
            sxy_jpj(indsxy)=var_jpj(indsxy)
            sxy_jpk(indsxy)=var_jpk(indsxy)
            sxy_jpt(indsxy)=var_jpt(indsxy)
            sxy_nbr(indsxy)=var_nbr(indsxy)
            sxyemsk(indsxy)=varemsk(indsxy)
            sxyfmsk(indsxy)=varfmsk(indsxy)
            sxyvmsk(indsxy)=FREAL4(varvmsk(indsxy))
            sxyngrd(indsxy)=dtangrd(indsxy)
            sxyfgrd(indsxy)=dtafgrd(indsxy)
            sxyegrd(indsxy)=dtaegrd(indsxy)
         ENDDO
C
         jtxy=1
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         xdim = varxdim(indsxy)
         ydim = varydim(indsxy)
         zdim = varzdim(indsxy)
         tdim = vartdim(indsxy)
C
      CASE(2)
C --- dta
         sxyend  = dtaend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=dta_ord(jsxy)
            indsxy= sxy_ord(jsxy)
C
            sxy_nam(indsxy)=dtaofil(indsxy)
            sxy_dim(indsxy)=dta_dim(indsxy)
            sxy_jpi(indsxy)=dta_jpi(indsxy)
            sxy_jpj(indsxy)=dta_jpj(indsxy)
            sxy_jpk(indsxy)=dta_jpk(indsxy)
            sxy_jpt(indsxy)=dta_jpt(indsxy)
            sxy_nbr(indsxy)=dta_nbr(indsxy)
            sxyemsk(indsxy)=dtaemsk(indsxy)
            sxyfmsk(indsxy)=dtafmsk(indsxy)
            sxyvmsk(indsxy)=FREAL4(dtavmsk(indsxy))
            sxyngrd(indsxy)=dtangrd(indsxy)
            sxyfgrd(indsxy)=dtafgrd(indsxy)
            sxyegrd(indsxy)=dtaegrd(indsxy)
         ENDDO
C
         jtxy=1
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         xdim = dtaxdim(indsxy)
         ydim = dtaydim(indsxy)
         zdim = dtazdim(indsxy)
         tdim = dtatdim(indsxy)
C
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Compute dimensions of NetCDF file
      jtxy=1
      jsxy=ktxy_indmsk(jtxy)
      indsxy=sxy_ord(jsxy)
C
      sxy_dim_max = sxy_dim(indsxy)
      jpiend = sxy_jpi(indsxy)
      jpjend = sxy_jpj(indsxy)
      jpkend = sxy_jpk(indsxy)
      jptend = sxy_jpt(indsxy)
C
      DO jtxy=2,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
C
         sxy_dim_max = max(sxy_dim(indsxy),sxy_dim_max)
         jpiend = max(jpiend,sxy_jpi(indsxy))
         jpjend = max(jpjend,sxy_jpj(indsxy))
         jpkend = max(jpkend,sxy_jpk(indsxy))
         jptend = max(jptend,sxy_jpt(indsxy))
      ENDDO
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writenfil/writenc'
         WRITE(numout,*) '    ==> WRITING file: ',
     $                            kfnoutnc(1:lenv(kfnoutnc))
      ENDIF
C
C -2.- Check if output file exists
C      If yes, open NetCDF file, check dimensions
C      If not, create output Netcdf file, and write header
C -------------------------------------------------------------
C
      filexists=.FALSE.
      INQUIRE (FILE=kfnoutnc,EXIST=filexists)
      IF (filexists) THEN
C
        ierr = NF90_OPEN(kfnoutnc,NF90_WRITE,idf)
        IF (ierr.NE.0) GOTO 101
C
      ELSE
C
        IF (ALL(sxyemsk(sxy_ord(1:sxyend)).EQ.4)) THEN
C
C Copy mask file as a prototype .nc output file (if in .nc format)
C
          jsxy=ktxy_indmsk(1)
          indsxy=sxy_ord(jsxy)
          cpcom='cp '//sxyfmsk(indsxy)(1:lenv(sxyfmsk(indsxy)))
     $               //' '//kfnoutnc(1:lenv(kfnoutnc))
          CALL system(cpcom)
C
          ierr = NF90_OPEN(kfnoutnc,NF90_WRITE,idf)
          IF (ierr.NE.0) GOTO 101
C
        ELSE
C
C Create .nc file from scratch using SESAM configuration
C         CALL createnc(kfnoutnc,ktxy_indvct,ktxy_indmsk,txyend,kflagxyo)
          GOTO 106
C
        ENDIF
C
      ENDIF
C
C Inquire about dimensions
      idx = -1
      IF (sxy_dim_max.GE.1) THEN
        ierr = NF90_INQ_DIMID(idf,xdim,idx)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=nx)
        IF (ierr.NE.0) GOTO 102
        IF (jpiend.NE.nx) GOTO 102
      ENDIF
C
      idy = -1
      IF (sxy_dim_max.GE.2) THEN
        ierr = NF90_INQ_DIMID(idf,ydim,idy)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=ny)
        IF (ierr.NE.0) GOTO 102
        IF (jpjend.NE.ny) GOTO 102
      ENDIF
C
      idz = -1 ; nz = 1
      IF (sxy_dim_max.GE.3) THEN
        ierr = NF90_INQ_DIMID(idf,zdim,idz)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=nz)
        IF (ierr.NE.0) GOTO 102
        IF (jpkend.NE.nz) GOTO 102
      ENDIF
C
      idt = -1 ; nt = 1
      IF (sxy_dim_max.GE.4) THEN
        ierr = NF90_INQ_DIMID(idf,tdim,idt)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=nt)
        IF (ierr.NE.0) GOTO 102
        IF (jptend.NE.nt) GOTO 102
      ENDIF
C
C -3.- Loop on variables to write in NetCDF file
C      Get variable id, get the set of dimensions spanned by the variable,
C      Get indices of x, y, z, t dimensions in that set of dimensions
C      Initialize start, count and strid array for writing Netcdf variable
C ------------------------------------------------------------------------
C
C --- allocation ptabij
      allocate ( ptabij(1:jpiend,1:jpjend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL4(0.0)
C
      DO jtxy=1,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         spval=sxyvmsk(indsxy)
         txy_indtab(jtxy)=ktxy_indvct(jtxy)
C
         ierr = NF90_INQ_VARID(idf,sxy_nam(indsxy),idv)
         IF (ierr.NE.0) GOTO 103
         ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
         IF (ierr.NE.0) GOTO 103
         IF (sxy_dim(indsxy).GT.ndims) GOTO 103
C
         allocate(idims(ndims),vstart(ndims),vcount(ndims),vstrid(ndims))
C
         ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
         IF (ierr.NE.0) GOTO 103
C
         idvx = 0 ; idvy = 0 ; idvz = 0 ; idvt = 0
         DO idd = 1,ndims
           IF (idims(idd).EQ.idx) idvx = idd
           IF (idims(idd).EQ.idy) idvy = idd
           IF (idims(idd).EQ.idz) idvz = idd
           IF (idims(idd).EQ.idt) idvt = idd
         ENDDO
C
         IF ( (sxy_dim(indsxy).GE.1).AND.(idvx.EQ.0) ) GOTO 105
         IF ( (sxy_dim(indsxy).GE.2).AND.(idvy.EQ.0) ) GOTO 105
         IF ( (sxy_dim(indsxy).GE.3).AND.(idvz.EQ.0) ) GOTO 105
         IF ( (sxy_dim(indsxy).GE.4).AND.(idvt.EQ.0) ) GOTO 105
C
         vstart(1:ndims)=1 ; vcount(1:ndims)=1 ; vstrid(1:ndims)=1
         DO idd = 1,ndims
           IF (idd.EQ.idvx) vcount(idd)=sxy_jpi(indsxy)
           IF (idd.EQ.idvy) vcount(idd)=sxy_jpj(indsxy)
         ENDDO
C
C -4.- Write 2D arrays in NetCDF file
C -----------------------------------
C
         DO jt=1,sxy_jpt(indsxy)
         DO jk=1,sxy_jpk(indsxy)
C
            sompartxynbr=0
            IF (txy_indtab(jtxy).LE.size(kvectsin,1)) THEN
               CALL unmk4vct(kvectsin(txy_indtab(jtxy):),
     $              ptabij(:,:),jk,jt,jsxy,
     $              sompartxynbr,spval,kflagxyo) 
            ELSE
               ptabij(:,:) = spval
            ENDIF
C
            IF (idvz.NE.0) vstart(idvz)=jk
            IF (idvt.NE.0) vstart(idvt)=jt
C
            ierr = NF90_PUT_VAR(idf,idv,ptabij,start=vstart,count=vcount,
     $                          stride=vstrid)
            IF (ierr.NE.0) GOTO 104
C
            txy_indtab(jtxy)=txy_indtab(jtxy)+sompartxynbr
C
         ENDDO
         ENDDO
C
         IF (allocated(vstart)) deallocate (vstart)
         IF (allocated(vcount)) deallocate (vcount)
         IF (allocated(vstrid)) deallocate (vstrid)
         IF (allocated(idims)) deallocate (idims)
C
      ENDDO
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
C
C Coherence test
      DO jtxy=1,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         IF (sxy_nbr(indsxy).NE.(txy_indtab(jtxy)-ktxy_indvct(jtxy)))
     $                   GOTO 1000
      ENDDO
C
C --- deallocation
      IF (allocated(ptabij)) deallocate (ptabij)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lionc','writenc')
 1001 CALL printerror2(0,1001,3,'lionc','writenc')
C
 101  WRITE (texterror,*) 'Bad NetCDF file: ',kfnoutnc
      CALL printerror2(0,101,3,'lionc','writenc',comment=texterror)
 102  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',kfnoutnc
      CALL printerror2(0,102,3,'lionc','writenc',comment=texterror)
 103  WRITE (texterror,*) 'Bad variable in NetCDF file: ',sxy_nam
      CALL printerror2(0,103,3,'lionc','writenc',comment=texterror)
 104  WRITE (texterror,*) 'Error writing NetCDF variable'
      CALL printerror2(0,104,3,'lionc','writenc',comment=texterror)
 105  WRITE (texterror,*) 'Bad dimensions for NetCDF variable: ',sxy_nam
      CALL printerror2(0,105,3,'lionc','writenc',comment=texterror)
 106  WRITE (texterror,*) 'Output .nc file does not exist: ',kfnoutnc
      CALL printerror2(0,106,3,'lionc','writenc',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnoutnc
      CALL printerror2(0,107,3,'lionc','writenc',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrmsknc(kfname,kjpi,kjpj,kjpk,kjpt,kindxy,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Get variable dimensions from mask files
CCC  -------
CCC  Method : Open, read and close NetCDF file header
CCC  ------
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfname
      INTEGER, intent(in) :: kindxy,kflagxyo
      INTEGER, intent(out) :: kjpi,kjpj,kjpk,kjpt
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80) :: kform
      CHARACTER(len=bgword) :: xdim, ydim, zdim, tdim
      INTEGER :: ierr, idf, idx, idy, idz, idt, indsxy
CC----------------------------------------------------------------------
      SELECT CASE (kflagxyo)
      CASE(1)
         indsxy = kindxy
         xdim = varxdim(indsxy)
         ydim = varydim(indsxy)
         zdim = varzdim(indsxy)
         tdim = vartdim(indsxy)
      CASE(2)
         indsxy = kindxy
         xdim = dtaxdim(indsxy)
         ydim = dtaydim(indsxy)
         zdim = dtazdim(indsxy)
         tdim = dtatdim(indsxy)
      END SELECT
C
C -1.- Read dimensions in NetCDF
C ------------------------------
C
      ierr = NF90_OPEN(kfname,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 101
C
      ierr = NF90_INQ_DIMID(idf,xdim,idx)
      IF (ierr.NE.0) THEN
        kjpi = 1
      ELSE
        ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=kjpi)
        IF (ierr.NE.0) GOTO 102
      ENDIF
C
      ierr = NF90_INQ_DIMID(idf,ydim,idy)
      IF (ierr.NE.0) THEN
        kjpj = 1
      ELSE
        ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=kjpj)
        IF (ierr.NE.0) GOTO 102
      ENDIF
C
      ierr = NF90_INQ_DIMID(idf,zdim,idz)
      IF (ierr.NE.0) THEN
        kjpk = 1
      ELSE
        ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=kjpk)
        IF (ierr.NE.0) GOTO 102
      ENDIF
C
      ierr = NF90_INQ_DIMID(idf,tdim,idt)
      IF (ierr.NE.0) THEN
        kjpt = 1
      ELSE
        ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=kjpt)
        IF (ierr.NE.0) GOTO 102
      ENDIF
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 106
C
C -2.- Control Print
C ------------------
C
      IF (nprint.GE.3) THEN 
         kform='(8x,a,4i5)'
         WRITE(numout,kform) '- Dimensions: ',kjpi,kjpj,kjpk,kjpt
      ENDIF
C
      RETURN
C
C --- error management
C
 101  WRITE (texterror,*) 'Bad NetCDF file: ',kfname
      CALL printerror2(0,101,3,'lionc','evalhdrmsknc',comment=texterror)
 102  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',kfname
      CALL printerror2(0,102,3,'lionc','evalhdrmsknc',comment=texterror)
 106  WRITE (texterror,*) 'Error closing Netcdf file: ',kfname
      CALL printerror2(0,106,3,'lionc','evalhdrmsknc',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmsknc(jsxy,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read mask arrays from NetCDF mask files (.cdta, .cdf)
CCC  -------
CCC  Method : Open, read and close NetCDF file
CCC  ------
CCC  Input : jsxy     : index variable to load
CCC  -----   kflagxyo : vector type (Vx,Vy)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use netcdf
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: jsxy,kflagxyo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: nx, ny, nz, nt, ndims
      BIGREAL4 :: sxyvmsk
      BIGREAL4, dimension(:), allocatable :: lev
      BIGREAL4, dimension(:,:), allocatable :: ptab
      CHARACTER(len=bgword) :: sxy_nam, sxyznam
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxypmsk,sxydmsk 
      INTEGER :: indsxy, ind_msk, js, sxynbr, allocok
      LOGICAL :: sxymsea
      CHARACTER(len=bgword) :: varunit, longname, kform, sxyfmsk
      CHARACTER(len=bgword), dimension(1:4) :: dimunit
      INTEGER :: ji, jj, jk, jt, jpiend, jpjend, jpkend, jptend
      INTEGER, allocatable, dimension(:) :: idims,vstart,vcount,vstrid
      INTEGER, allocatable, dimension(:) :: zstart,zcount
      INTEGER :: ierr, idf, idx, idy, idz, idt, idv, idd
      INTEGER :: idvx, idvy, idvz, idvt, idvzz
      CHARACTER(len=bgword) :: xdim, ydim, zdim, tdim
CC----------------------------------------------------------------------
C
C -1.- Get variable characteristics from SESAM configuration
C ----------------------------------------------------------
C
      SELECT CASE(kflagxyo)
      CASE(1)
C --- var
         indsxy  = var_ord(jsxy)
         sxy_nam = varifil(indsxy)
         sxy_dim = var_dim(indsxy)
         sxy_jpi = var_jpi(indsxy)
         sxy_jpj = var_jpj(indsxy)
         sxy_jpk = var_jpk(indsxy)
         sxy_jpt = var_jpt(indsxy)
         sxyfmsk = varfmsk(indsxy)
         sxypmsk = varpmsk(indsxy)
         sxydmsk = vardmsk(indsxy)
         sxyvmsk = FREAL4(varvmsk(indsxy))
         sxymsea = varmsea(indsxy)
         xdim = varxdim(indsxy)
         ydim = varydim(indsxy)
         zdim = varzdim(indsxy)
         tdim = vartdim(indsxy)
         ind_msk = jsxy-1
      CASE(2)
C --- dta
         indsxy  = dta_ord(jsxy)
         sxy_nam = dtaifil(indsxy)
         sxy_dim = dta_dim(indsxy)
         sxy_jpi = dta_jpi(indsxy)
         sxy_jpj = dta_jpj(indsxy)
         sxy_jpk = dta_jpk(indsxy)
         sxy_jpt = dta_jpt(indsxy)
         sxyfmsk = dtafmsk(indsxy)
         sxypmsk = dtapmsk(indsxy)
         sxydmsk = dtadmsk(indsxy)
         sxyvmsk = FREAL4(dtavmsk(indsxy))
         sxymsea = dtamsea(indsxy)
         xdim = dtaxdim(indsxy)
         ydim = dtaydim(indsxy)
         zdim = dtazdim(indsxy)
         tdim = dtatdim(indsxy)
         sxyznam = dtaznam(indsxy)
         ind_msk = jsxy-1+varend
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Control print
      IF (nprint.GE.2) THEN 
         WRITE(numout,*) '*** ROUTINE : sesam/readmsk/readmsknc'
         WRITE(numout,*) '    ==> READING ',sxyfmsk(1:lenv(sxyfmsk))
      ENDIF
C
C Set size of array to read
      jpiend = 1 ; jpjend = 1 ; jpkend = 1 ; jptend = 1
      IF (sxydmsk.GE.1) jpiend = sxy_jpi
      IF (sxydmsk.GE.2) jpjend = sxy_jpj
      IF (sxydmsk.GE.3) jpkend = sxy_jpk
      IF (sxydmsk.GE.4) jptend = sxy_jpt
C
C -2.- Open NetCDF file, get and check dimensions
C -----------------------------------------------
C
      ierr = NF90_OPEN(sxyfmsk,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 101
C
      idx = -1
      IF (sxy_dim.GE.1) THEN
        WRITE(numout,*)
        ierr = NF90_INQ_DIMID(idf,xdim,idx)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=nx)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpi.NE.nx) GOTO 102
      ENDIF
C
      idy = -1
      IF (sxy_dim.GE.2) THEN
        ierr = NF90_INQ_DIMID(idf,ydim,idy)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=ny)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpj.NE.ny) GOTO 102
      ENDIF
C
      idz = -1 ; nz = 1
      IF (sxy_dim.GE.3) THEN
        ierr = NF90_INQ_DIMID(idf,zdim,idz)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=nz)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpk.NE.nz) GOTO 102
        IF (size(mask,3).GT.nz) GOTO 102
      ENDIF
C
      idt = -1 ; nt = 1
      IF (sxy_dim.GE.4) THEN
        ierr = NF90_INQ_DIMID(idf,tdim,idt)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=nt)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpt.NE.nt) GOTO 102
      ENDIF
C
C -3.- Get variable id, get the set of dimensions spanned by the variable,
C      Get indices of x, y, z, t dimensions in that set of dimensions
C      Initialize start, count and strid array for reading Netcdf variable
C ------------------------------------------------------------------------
C
      ierr = NF90_INQ_VARID(idf,sxy_nam,idv)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
      IF (ierr.NE.0) GOTO 103
      IF (sxy_dim.GT.ndims) GOTO 103
C
      allocate(idims(ndims),vstart(ndims),vcount(ndims),vstrid(ndims))
      vstart(1:ndims)=1 ; vcount(1:ndims)=1 ; vstrid(1:ndims)=1
C 
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
      IF (ierr.NE.0) GOTO 103
C
      idvx = 0 ; idvy = 0 ; idvz = 0 ; idvt = 0
      DO idd = 1,ndims
        IF (idims(idd).EQ.idx) idvx = idd
        IF (idims(idd).EQ.idy) idvy = idd
        IF (idims(idd).EQ.idz) idvz = idd
        IF (idims(idd).EQ.idt) idvt = idd
      ENDDO
C
      IF ( (sxy_dim.GE.1).AND.(idvx.EQ.0) ) GOTO 105
      IF ( (sxy_dim.GE.2).AND.(idvy.EQ.0) ) GOTO 105
      IF ( (sxy_dim.GE.3).AND.(idvz.EQ.0) ) GOTO 105
      IF ( (sxy_dim.GE.4).AND.(idvt.EQ.0) ) GOTO 105
C
      DO idd = 1,ndims
        IF (idd.EQ.idvx) vcount(idd)=sxy_jpi
        IF (idd.EQ.idvy) vcount(idd)=sxy_jpj
      ENDDO
C
C -2.- Read vertical levels (for Vx object only)
C ----------------------------------------------
C Separate this operations from the reading of the mask
C Create for this a new routine readlev, only called in obsv module
C
      IF ((kflagxyo.EQ.2).AND.(sxy_dim.GE.3)) THEN
C allocation lev
        allocate ( lev(1:nz), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        lev(:) = FREAL4(0.0)
C
C       Read vertical levels from files
        ierr = NF90_INQ_VARID(idf,sxyznam,idvzz)
        IF (ierr.NE.0) THEN
C
          IF (sxydmsk.EQ.2) THEN
            var_lev(1:sxy_jpk,indsxy) = FREAL(0.0)
          ELSEIF (sxydmsk.EQ.sxy_dim) THEN
            GOTO 103
          ENDIF
C
        ELSE
C
          allocate(zstart(1),zcount(1))
          zstart(1)=1 ; zcount(1)=nz
          ierr = NF90_GET_VAR(idf,idvzz,lev,start=zstart,count=zcount)
          deallocate(zstart,zcount)
          IF (ierr.NE.0) GOTO 104
C
          IF (sxydmsk.EQ.2) THEN
            var_lev(1:sxy_jpk,indsxy) = FREAL(lev(1))
          ELSEIF (sxydmsk.EQ.sxy_dim) THEN
            var_lev(1:jpkend,indsxy)=FREAL(lev(1:jpkend))
          ENDIF
C
        ENDIF
C
        IF (allocated(lev)) deallocate(lev)
C
      ENDIF
C
C -3.- Read CDF mask arrays
C -------------------------
C
C --- allocation ptab
      allocate ( ptab(1:jpiend,1:jpjend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptab(:,:) = FREAL4(0.0)
C
      js=1
      DO jt=1,jptend
      DO jk=1,jpkend
C
         IF (idvz.NE.0) vstart(idvz)=jk
         IF (idvt.NE.0) vstart(idvt)=jt
C
         ierr = NF90_GET_VAR(idf,idv,ptab,start=vstart,count=vcount,
     $                       stride=vstrid)
         IF (ierr.NE.0) GOTO 104
C
         DO jj=1,jpjend
         DO ji=1,jpiend
            IF (sxymsea.EQV.(sxyvmsk.EQ.ptab(ji,jj))) THEN
               mask(ji,jj,jk,jt)=mask(ji,jj,jk,jt)+(2**(ind_msk))
               js=js+1
            ENDIF
         ENDDO
         ENDDO
C
      ENDDO
      ENDDO
      sxynbr=js-1
C
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 106
C
      IF (allocated(ptab)) deallocate(ptab)
C
C -5.- Control Print :
C --------------------
C
      IF (nprint.GE.3) THEN 
         kform='(8x,a,5i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndims
         kform='(8x,a,e12.3)'
         IF (  (jpkend.NE.nz).OR.(jptend.NE.nt) ) THEN
            kform='(8x,a,4i5)'
            WRITE(numout,kform) '- Size of loaded array : ',jpiend,
     $                                  jpjend,jpkend,jptend
         ENDIF
         kform='(8x,a,i9)'
         WRITE(numout,kform) '- Number of loaded variables : ',sxynbr
         kform='(8x,a)'
         WRITE(numout,kform) '- Vertical levels : '
         kform='(8x,7e12.3)'
         WRITE(numout,kform) var_lev(1:jpkend,indsxy)
      ENDIF
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'lionc','readmsknc')
 1001 CALL printerror2(0,1001,3,'lionc','readmsknc')
C
 101  WRITE (texterror,*) 'Bad NetCDF file: ',sxyfmsk
      CALL printerror2(0,101,3,'lionc','readmsknc',comment=texterror)
 102  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',sxyfmsk
      CALL printerror2(0,102,3,'lionc','readmsknc',comment=texterror)
 103  WRITE (texterror,*) 'Bad variable in NetCDF file: ',sxy_nam
      CALL printerror2(0,103,3,'lionc','readmsknc',comment=texterror)
 104  WRITE (texterror,*) 'Error reading NetCDF variable: ',sxy_nam
      CALL printerror2(0,104,3,'lionc','readmsknc',comment=texterror)
 105  WRITE (texterror,*) 'Bad dimensions for NetCDF variable: ',sxy_nam
      CALL printerror2(0,105,3,'lionc','readmsknc',comment=texterror)
 106  WRITE (texterror,*) 'Error closing Netcdf file: ',sxyfmsk
      CALL printerror2(0,106,3,'lionc','readmsknc',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE lionc
