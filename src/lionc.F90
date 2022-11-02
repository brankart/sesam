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
! ---                  LIONC.F90                                  ---
! ---                                                           ---
! --- original     : 07-11 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- Routines to read and write from generic NetCDF files
! ---
! --- SUBROUTINE  readnc       : Read variable field from NetCDF file
! --- SUBROUTINE  writenc      : Write variable field in NetCDF file
! --- SUBROUTINE  evalhdrmsknc : Get variable dimensions from mask files
! --- SUBROUTINE  readmsknc    : Read mask arrays from NetCDF mask files
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE lionc
      use mod_main
      use hiogrd
      use utilvct
      IMPLICIT NONE
      PRIVATE

      PUBLIC readnc,writenc,evalhdrmsknc,readmsknc

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readnc(kfname,kvectsout,kjsxy,ksomsxynbr,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Read variable field (of Vx or Vy object) from NetCDF file
!  -------
!  Method : Open, read and close NetCDF file (.nc or .ncdta)
!  ------
!  Input :  kfname     : Filename
!  -----    kjsxy      : Index of variable field to read
!           kflagxyo   : vector type (Vx,Vy)
!  Output : kvectsout  : variable field 1D vector (jpx or jpy)
!  ------   ksomsxynbr : number of values loaded
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
      CHARACTER(len=*), intent(in) :: kfname
      BIGREAL, dimension(:), intent(out) :: kvectsout
      INTEGER, intent(in) :: kjsxy,kflagxyo
      INTEGER, intent(out) :: ksomsxynbr
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jpssize
      BIGREAL4, allocatable, dimension(:,:) :: ptabij
      INTEGER :: indsxy,js,sompartsxynbr,jextdta,jextvar
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt
      CHARACTER(len=bgword) :: sxy_nam, kform
      CHARACTER(len=bgword), dimension(1:4) :: dimunit
!
      INTEGER, allocatable, dimension(:) :: idims,vstart,vcount,vstrid
      INTEGER :: ierr, idf, idx, idy, idz, idt, idv, idd
      INTEGER :: idvx, idvy, idvz, idvt
      CHARACTER(len=bgword) :: xdim, ydim, zdim, tdim
      INTEGER :: nx, ny, nz, nt, ndims
      INTEGER :: ji, jj, jk, jt
      INTEGER :: jpiend, jpjend, jpkend, jptend
      LOGICAL :: varfile
!----------------------------------------------------------------------
! Size of output vector
      jpssize=size(kvectsout,1)
!
! -1.- Get variable characteristics from SESAM configuration
! ----------------------------------------------------------
!
! Get file type (Vx or Vy) and variable index
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         indsxy  = var_ord(kjsxy)
         varfile = .TRUE.
      CASE(2)
! --- dta
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
!
! Set file dimensions
      sxy_jpi = 1
      sxy_jpj = 1
      sxy_jpk = 1
      sxy_jpt = 1
      IF (varfile) THEN
! --- read Vx or Vy variable field from var file
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
! --- read Vy variable field from dta file
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
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readfil/readnc'
         WRITE(numout,*) '    ==> READING file ',kfname(1:lenv(kfname))
         WRITE(numout,*) '    ==> READING variable: ',sxy_nam
      ENDIF
!
! Set size of array to read
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
!
! -2.- Open NetCDF file, get and check dimensions
! -----------------------------------------------
!
      ierr = NF90_OPEN(kfname,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 101
!
      IF (sxy_dim.GE.1) THEN
        ierr = NF90_INQ_DIMID(idf,xdim,idx)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=nx)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpi.NE.nx) GOTO 102
      ENDIF
!
      IF (sxy_dim.GE.2) THEN
        ierr = NF90_INQ_DIMID(idf,ydim,idy)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=ny)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpj.NE.ny) GOTO 102
      ENDIF
!
      IF (sxy_dim.GE.3) THEN
        IF (zdim.NE.'none') THEN
          ierr = NF90_INQ_DIMID(idf,zdim,idz)
          IF (ierr.NE.0) GOTO 102
          ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=nz)
          IF (ierr.NE.0) GOTO 102
        ELSE
          nz = 1
        ENDIF
        IF (sxy_jpk.NE.nz) GOTO 102
      ENDIF
!
      IF (sxy_dim.GE.4) THEN
        ierr = NF90_INQ_DIMID(idf,tdim,idt)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=nt)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpt.NE.nt) GOTO 102
      ENDIF
!
! -3.- Get variable id, get the set of dimensions spanned by the variable,
!      Get indices of x, y, z, t dimensions in that set of dimensions
!      Initialize start, count and strid array for reading Netcdf variable
! ------------------------------------------------------------------------
!
      ierr = NF90_INQ_VARID(idf,sxy_nam,idv)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
      IF (ierr.NE.0) GOTO 103
      IF (zdim.NE.'none') THEN
        IF (sxy_dim.GT.ndims) GOTO 103
      ELSE
        IF (sxy_dim.GT.ndims+1) GOTO 103
      ENDIF
!
      allocate(idims(ndims),vstart(ndims),vcount(ndims),vstrid(ndims))
      vstart(1:ndims)=1 ; vcount(1:ndims)=1 ; vstrid(1:ndims)=1
! 
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
      IF (ierr.NE.0) GOTO 103
!
      idvx = 0 ; idvy = 0 ; idvz = 0 ; idvt = 0
      DO idd = 1,ndims
        IF (idims(idd).EQ.idx) idvx = idd
        IF (idims(idd).EQ.idy) idvy = idd
        IF (idims(idd).EQ.idz) idvz = idd
        IF (idims(idd).EQ.idt) idvt = idd
      ENDDO
!
      IF ( (sxy_dim.GE.1).AND.(idvx.EQ.0) ) GOTO 105
      IF ( (sxy_dim.GE.2).AND.(idvy.EQ.0) ) GOTO 105
      IF (zdim.NE.'none') THEN
        IF ( (sxy_dim.GE.3).AND.(idvz.EQ.0) ) GOTO 105
      ENDIF
      IF ( (sxy_dim.GE.4).AND.(idvt.EQ.0) ) GOTO 105
!
      DO idd = 1,ndims
        IF (idd.EQ.idvx) vcount(idd)=sxy_jpi
        IF (idd.EQ.idvy) vcount(idd)=sxy_jpj
      ENDDO
!
! -4.- Read 2D arrays in NetCDF file
! ----------------------------------
!
! --- allocation ptabij
      allocate ( ptabij(1:sxy_jpi,1:sxy_jpj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL4(0.0)
!
      js=1
      DO jt=1,jptend
      DO jk=1,jpkend
!
         IF (idvz.NE.0) vstart(idvz)=jk
         IF (idvt.NE.0) vstart(idvt)=jt
!
         ierr = NF90_GET_VAR(idf,idv,ptabij,start=vstart,count=vcount, &
     &                       stride=vstrid)
         IF (ierr.NE.0) GOTO 104
!
         CALL mk4vct(kvectsout(js:),ptabij(:,:), &
     &               jk,jt,kjsxy,sompartsxynbr,kflagxyo)
         js  = js + sompartsxynbr
!
      ENDDO
      ENDDO
!
      ksomsxynbr=js-1
!
      IF (allocated(ptabij)) deallocate(ptabij)
!
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 106
!
! -5.- Control Print
! ------------------
!
      IF (nprint.GE.3) THEN
         kform='(8x,a,4i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt
         IF (  (jpkend.NE.nz).OR.(jptend.NE.nt) ) THEN
            kform='(8x,a,4i5)'
            WRITE(numout,kform) '- Size of loaded array : ',jpiend, &
     &                                  jpjend,jpkend,jptend
         ENDIF
         kform='(8x,a,i9)'
         WRITE(numout,kform) '- Number of loaded values : ',ksomsxynbr
      ENDIF
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'lionc','readnc')
 1001 CALL printerror2(0,1001,3,'lionc','readnc')
!
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
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writenc(kfnoutnc,kvectsin, &
     &           ktxy_indvct,ktxy_indmsk,txyend,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Write variable field (of Vx or Vy object) in NetCDF file
!  -------
!  Method : Open, write and close NetCDF file (.nc or .ncdta)
!  ------
!  Input : kfnoutnc    : Filename
!  -----   kvectsin    : 1D SESAM vector object (Vx or Vy)
!          ktxy_indvct : array with variable field beginning indices
!                                  (in vector object Vx or Vy)
!          ktxy_indmsk : table with variable field indices
!                                  (in variable list)
!          txyend      : number of variable fields to write
!          kflagxyo    : vector type (Vx,Vy)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_coord
      use mod_spacexyo , only : spvaldta,spvalvar
      use netcdf
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutnc
      BIGREAL, dimension(:), intent(in) :: kvectsin
      INTEGER, dimension(:), intent(in) :: ktxy_indvct,ktxy_indmsk
      INTEGER, intent(in) :: txyend,kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: nx, ny, nz, nt
      BIGREAL4 :: spval
      BIGREAL4, dimension(1:nbvar) :: sxyvmsk
      INTEGER  ::allocok
      BIGREAL4, allocatable, dimension(:,:) :: ptabij
      INTEGER :: sxyend, sxy_dim_max
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_dim,sxyemsk, &
     &     sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxy_nbr,sxyngrd,sxyegrd
      INTEGER :: jsxy,indsxy,sompartxynbr,jtxy
      INTEGER, dimension(1:nbvar) :: txy_indtab
      CHARACTER(len=bgword), dimension(1:nbvar) :: sxy_nam,sxyzdim
      CHARACTER(len=bgword), dimension(1:nbvar) :: sxyfgrd,sxyfmsk
      INTEGER :: ji, jj, jk, jt
      INTEGER :: jpiend, jpjend, jpkend, jptend
      LOGICAL :: filexists
      INTEGER, allocatable, dimension(:) :: idims,vstart,vcount,vstrid
      INTEGER :: ierr, idf, idx, idy, idz, idt, idv, idd
      INTEGER :: idvx, idvy, idvz, idvt, ndims
      CHARACTER(len=bgword) :: xdim, ydim, zdim, tdim
      CHARACTER(len=bgword) :: cpcom
!----------------------------------------------------------------------
!
! -1.- Get variable characteristics from SESAM configuration
! ----------------------------------------------------------
!
      zdim = 'none'

      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         sxyend  = varend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=var_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
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
            sxyzdim(indsxy)=varzdim(indsxy)
         ENDDO
!
         jtxy=1
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         xdim = varxdim(indsxy)
         ydim = varydim(indsxy)
         tdim = vartdim(indsxy)
         DO jtxy = 1,txyend
           jsxy=ktxy_indmsk(jtxy)
           indsxy=sxy_ord(jsxy)
           IF (varzdim(indsxy).NE.'none') THEN
              zdim = varzdim(indsxy)
           ENDIF
         ENDDO
!
      CASE(2)
! --- dta
         sxyend  = dtaend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=dta_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
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
            sxyzdim(indsxy)=dtazdim(indsxy)
         ENDDO
!
         jtxy=1
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         xdim = dtaxdim(indsxy)
         ydim = dtaydim(indsxy)
         tdim = dtatdim(indsxy)
         DO jtxy = 1,txyend
           jsxy=ktxy_indmsk(jtxy)
           indsxy=sxy_ord(jsxy)
           IF (dtazdim(indsxy).NE.'none') THEN
              zdim = dtazdim(indsxy)
           ENDIF
         ENDDO
!
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Compute dimensions of NetCDF file
      jtxy=1
      jsxy=ktxy_indmsk(jtxy)
      indsxy=sxy_ord(jsxy)
!
      sxy_dim_max = sxy_dim(indsxy)
      jpiend = sxy_jpi(indsxy)
      jpjend = sxy_jpj(indsxy)
      jpkend = sxy_jpk(indsxy)
      jptend = sxy_jpt(indsxy)
!
      DO jtxy=2,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
!
         sxy_dim_max = max(sxy_dim(indsxy),sxy_dim_max)
         jpiend = max(jpiend,sxy_jpi(indsxy))
         jpjend = max(jpjend,sxy_jpj(indsxy))
         jpkend = max(jpkend,sxy_jpk(indsxy))
         jptend = max(jptend,sxy_jpt(indsxy))
      ENDDO
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writenfil/writenc'
         WRITE(numout,*) '    ==> WRITING file: ', &
     &                            kfnoutnc(1:lenv(kfnoutnc))
      ENDIF
!
! -2.- Check if output file exists
!      If yes, open NetCDF file, check dimensions
!      If not, create output Netcdf file, and write header
! -------------------------------------------------------------
!
      filexists=.FALSE.
      INQUIRE (FILE=kfnoutnc,EXIST=filexists)
      IF (filexists) THEN
!
        ierr = NF90_OPEN(kfnoutnc,NF90_WRITE,idf)
        IF (ierr.NE.0) GOTO 101
!
      ELSE
!
        IF (ALL(sxyemsk(sxy_ord(1:sxyend)).EQ.4)) THEN
!
! Copy mask file as a prototype .nc output file (if in .nc format)
!
          jsxy=ktxy_indmsk(1)
          indsxy=sxy_ord(jsxy)
          cpcom='cp '//sxyfmsk(indsxy)(1:lenv(sxyfmsk(indsxy))) &
     &               //' '//kfnoutnc(1:lenv(kfnoutnc))
          CALL system(cpcom)
!
          ierr = NF90_OPEN(kfnoutnc,NF90_WRITE,idf)
          IF (ierr.NE.0) GOTO 101
!
        ELSE
!
! Create .nc file from scratch using SESAM configuration
!         CALL createnc(kfnoutnc,ktxy_indvct,ktxy_indmsk,txyend,kflagxyo)
          GOTO 106
!
        ENDIF
!
      ENDIF
!
! Inquire about dimensions
      idx = -1
      IF (sxy_dim_max.GE.1) THEN
        ierr = NF90_INQ_DIMID(idf,xdim,idx)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=nx)
        IF (ierr.NE.0) GOTO 102
        IF (jpiend.NE.nx) GOTO 102
      ENDIF
!
      idy = -1
      IF (sxy_dim_max.GE.2) THEN
        ierr = NF90_INQ_DIMID(idf,ydim,idy)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=ny)
        IF (ierr.NE.0) GOTO 102
        IF (jpjend.NE.ny) GOTO 102
      ENDIF
!
      idz = -1 ; nz = 1
      IF (sxy_dim_max.GE.3) THEN
        IF (zdim.NE.'none') THEN
          ierr = NF90_INQ_DIMID(idf,zdim,idz)
          IF (ierr.NE.0) GOTO 102
          ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=nz)
          IF (ierr.NE.0) GOTO 102
        ELSE
          nz = 1
        ENDIF
        IF (jpkend.NE.nz) GOTO 102
      ENDIF
!
      idt = -1 ; nt = 1
      IF (sxy_dim_max.GE.4) THEN
        ierr = NF90_INQ_DIMID(idf,tdim,idt)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=nt)
        IF (ierr.NE.0) GOTO 102
        IF (jptend.NE.nt) GOTO 102
      ENDIF
!
! -3.- Loop on variables to write in NetCDF file
!      Get variable id, get the set of dimensions spanned by the variable,
!      Get indices of x, y, z, t dimensions in that set of dimensions
!      Initialize start, count and strid array for writing Netcdf variable
! ------------------------------------------------------------------------
!
! --- allocation ptabij
      allocate ( ptabij(1:jpiend,1:jpjend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL4(0.0)
!
      DO jtxy=1,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         spval=sxyvmsk(indsxy)
         txy_indtab(jtxy)=ktxy_indvct(jtxy)
!
         ierr = NF90_INQ_VARID(idf,sxy_nam(indsxy),idv)
         IF (ierr.NE.0) GOTO 103
         ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
         IF (ierr.NE.0) GOTO 103
         IF (sxyzdim(indsxy).NE.'none') THEN
           IF (sxy_dim(indsxy).GT.ndims) GOTO 103
         ELSE
           IF (sxy_dim(indsxy).GT.ndims+1) GOTO 103
         ENDIF
!
         allocate(idims(ndims),vstart(ndims),vcount(ndims),vstrid(ndims))
!
         ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
         IF (ierr.NE.0) GOTO 103
!
         idvx = 0 ; idvy = 0 ; idvz = 0 ; idvt = 0
         DO idd = 1,ndims
           IF (idims(idd).EQ.idx) idvx = idd
           IF (idims(idd).EQ.idy) idvy = idd
           IF (idims(idd).EQ.idz) idvz = idd
           IF (idims(idd).EQ.idt) idvt = idd
         ENDDO
!
         IF ( (sxy_dim(indsxy).GE.1).AND.(idvx.EQ.0) ) GOTO 105
         IF ( (sxy_dim(indsxy).GE.2).AND.(idvy.EQ.0) ) GOTO 105
         IF (sxyzdim(indsxy).NE.'none') THEN
           IF ( (sxy_dim(indsxy).GE.3).AND.(idvz.EQ.0) ) GOTO 105
         ENDIF
         IF ( (sxy_dim(indsxy).GE.4).AND.(idvt.EQ.0) ) GOTO 105
!
         vstart(1:ndims)=1 ; vcount(1:ndims)=1 ; vstrid(1:ndims)=1
         DO idd = 1,ndims
           IF (idd.EQ.idvx) vcount(idd)=sxy_jpi(indsxy)
           IF (idd.EQ.idvy) vcount(idd)=sxy_jpj(indsxy)
         ENDDO
!
! -4.- Write 2D arrays in NetCDF file
! -----------------------------------
!
         DO jt=1,sxy_jpt(indsxy)
         DO jk=1,sxy_jpk(indsxy)
!
            sompartxynbr=0
            IF (txy_indtab(jtxy).LE.size(kvectsin,1)) THEN
               CALL unmk4vct(kvectsin(txy_indtab(jtxy):), &
     &              ptabij(:,:),jk,jt,jsxy, &
     &              sompartxynbr,spval,kflagxyo) 
            ELSE
               ptabij(:,:) = spval
            ENDIF
!
            IF (idvz.NE.0) vstart(idvz)=jk
            IF (idvt.NE.0) vstart(idvt)=jt
!
            ierr = NF90_PUT_VAR(idf,idv,ptabij,start=vstart,count=vcount, &
     &                          stride=vstrid)
            IF (ierr.NE.0) GOTO 104
!
            txy_indtab(jtxy)=txy_indtab(jtxy)+sompartxynbr
!
         ENDDO
         ENDDO
!
         IF (allocated(vstart)) deallocate (vstart)
         IF (allocated(vcount)) deallocate (vcount)
         IF (allocated(vstrid)) deallocate (vstrid)
         IF (allocated(idims)) deallocate (idims)
!
      ENDDO
!
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 107
!
! Coherence test
      DO jtxy=1,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         IF (sxy_nbr(indsxy).NE.(txy_indtab(jtxy)-ktxy_indvct(jtxy))) &
     &                   GOTO 1000
      ENDDO
!
! --- deallocation
      IF (allocated(ptabij)) deallocate (ptabij)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lionc','writenc')
 1001 CALL printerror2(0,1001,3,'lionc','writenc')
!
 101  WRITE (texterror,*) 'Bad NetCDF file: ',kfnoutnc
      CALL printerror2(0,101,3,'lionc','writenc',comment=texterror)
 102  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',kfnoutnc
      CALL printerror2(0,102,3,'lionc','writenc',comment=texterror)
 103  WRITE (texterror,*) 'Bad variable in NetCDF file: ',sxy_nam(indsxy)
      CALL printerror2(0,103,3,'lionc','writenc',comment=texterror)
 104  WRITE (texterror,*) 'Error writing NetCDF variable'
      CALL printerror2(0,104,3,'lionc','writenc',comment=texterror)
 105  WRITE (texterror,*) 'Bad dimensions for NetCDF variable: ',sxy_nam(indsxy)
      CALL printerror2(0,105,3,'lionc','writenc',comment=texterror)
 106  WRITE (texterror,*) 'Output .nc file does not exist: ',kfnoutnc
      CALL printerror2(0,106,3,'lionc','writenc',comment=texterror)
 107  WRITE (texterror,*) 'Error closing Netcdf file: ',kfnoutnc
      CALL printerror2(0,107,3,'lionc','writenc',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrmsknc(kfname,kjpi,kjpj,kjpk,kjpt,kindxy,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Get variable dimensions from mask files
!  -------
!  Method : Open, read and close NetCDF file header
!  ------
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
      CHARACTER(len=*), intent(in) :: kfname
      INTEGER, intent(in) :: kindxy,kflagxyo
      INTEGER, intent(out) :: kjpi,kjpj,kjpk,kjpt
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=word80) :: kform
      CHARACTER(len=bgword) :: xdim, ydim, zdim, tdim
      INTEGER :: ierr, idf, idx, idy, idz, idt, indsxy
!----------------------------------------------------------------------
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
!
! -1.- Read dimensions in NetCDF
! ------------------------------
!
      ierr = NF90_OPEN(kfname,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 101
!
      ierr = NF90_INQ_DIMID(idf,xdim,idx)
      IF (ierr.NE.0) THEN
        kjpi = 1
      ELSE
        ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=kjpi)
        IF (ierr.NE.0) GOTO 102
      ENDIF
!
      ierr = NF90_INQ_DIMID(idf,ydim,idy)
      IF (ierr.NE.0) THEN
        kjpj = 1
      ELSE
        ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=kjpj)
        IF (ierr.NE.0) GOTO 102
      ENDIF
!
      ierr = NF90_INQ_DIMID(idf,zdim,idz)
      IF (ierr.NE.0) THEN
        kjpk = 1
      ELSE
        ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=kjpk)
        IF (ierr.NE.0) GOTO 102
      ENDIF
!
      ierr = NF90_INQ_DIMID(idf,tdim,idt)
      IF (ierr.NE.0) THEN
        kjpt = 1
      ELSE
        ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=kjpt)
        IF (ierr.NE.0) GOTO 102
      ENDIF
!
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 106
!
! -2.- Control Print
! ------------------
!
      IF (nprint.GE.3) THEN 
         kform='(8x,a,4i5)'
         WRITE(numout,kform) '- Dimensions: ',kjpi,kjpj,kjpk,kjpt
      ENDIF
!
      RETURN
!
! --- error management
!
 101  WRITE (texterror,*) 'Bad NetCDF file: ',kfname
      CALL printerror2(0,101,3,'lionc','evalhdrmsknc',comment=texterror)
 102  WRITE (texterror,*) 'Bad dimension in NetCDF file: ',kfname
      CALL printerror2(0,102,3,'lionc','evalhdrmsknc',comment=texterror)
 106  WRITE (texterror,*) 'Error closing Netcdf file: ',kfname
      CALL printerror2(0,106,3,'lionc','evalhdrmsknc',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmsknc(jsxy,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Read mask arrays from NetCDF mask files (.cdta, .cdf)
!  -------
!  Method : Open, read and close NetCDF file
!  ------
!  Input : jsxy     : index variable to load
!  -----   kflagxyo : vector type (Vx,Vy)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use netcdf
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(in) :: jsxy,kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
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
!----------------------------------------------------------------------
!
! -1.- Get variable characteristics from SESAM configuration
! ----------------------------------------------------------
!
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
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
! --- dta
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
!
! Control print
      IF (nprint.GE.2) THEN 
         WRITE(numout,*) '*** ROUTINE : sesam/readmsk/readmsknc'
         WRITE(numout,*) '    ==> READING ',sxyfmsk(1:lenv(sxyfmsk))
      ENDIF
!
! Set size of array to read
      jpiend = 1 ; jpjend = 1 ; jpkend = 1 ; jptend = 1
      IF (sxydmsk.GE.1) jpiend = sxy_jpi
      IF (sxydmsk.GE.2) jpjend = sxy_jpj
      IF (sxydmsk.GE.3) jpkend = sxy_jpk
      IF (sxydmsk.GE.4) jptend = sxy_jpt
!
! -2.- Open NetCDF file, get and check dimensions
! -----------------------------------------------
!
      ierr = NF90_OPEN(sxyfmsk,NF90_NOWRITE,idf)
      IF (ierr.NE.0) GOTO 101
!
      idx = -1
      IF (sxy_dim.GE.1) THEN
        WRITE(numout,*)
        ierr = NF90_INQ_DIMID(idf,xdim,idx)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idx,len=nx)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpi.NE.nx) GOTO 102
      ENDIF
!
      idy = -1
      IF (sxy_dim.GE.2) THEN
        ierr = NF90_INQ_DIMID(idf,ydim,idy)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idy,len=ny)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpj.NE.ny) GOTO 102
      ENDIF
!
      idz = -1 ; nz = 1
      IF (sxy_dim.GE.3) THEN
        IF (zdim.NE.'none') THEN
          ierr = NF90_INQ_DIMID(idf,zdim,idz)
          IF (ierr.NE.0) GOTO 102
          ierr = NF90_INQUIRE_DIMENSION(idf,idz,len=nz)
          IF (ierr.NE.0) GOTO 102
        ELSE
          nz=1
        ENDIF
        IF (sxy_jpk.NE.nz) GOTO 102
        !IF (size(mask,3).GT.nz) GOTO 102
      ENDIF
!
      idt = -1 ; nt = 1
      IF (sxy_dim.GE.4) THEN
        ierr = NF90_INQ_DIMID(idf,tdim,idt)
        IF (ierr.NE.0) GOTO 102
        ierr = NF90_INQUIRE_DIMENSION(idf,idt,len=nt)
        IF (ierr.NE.0) GOTO 102
        IF (sxy_jpt.NE.nt) GOTO 102
      ENDIF
!
! -3.- Get variable id, get the set of dimensions spanned by the variable,
!      Get indices of x, y, z, t dimensions in that set of dimensions
!      Initialize start, count and strid array for reading Netcdf variable
! ------------------------------------------------------------------------
!
      ierr = NF90_INQ_VARID(idf,sxy_nam,idv)
      IF (ierr.NE.0) GOTO 103
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,ndims=ndims)
      IF (ierr.NE.0) GOTO 103
      IF (zdim.NE.'none') THEN
        IF (sxy_dim.GT.ndims) GOTO 103
      ELSE
        IF (sxy_dim.GT.ndims+1) GOTO 103
      ENDIF
!
      allocate(idims(ndims),vstart(ndims),vcount(ndims),vstrid(ndims))
      vstart(1:ndims)=1 ; vcount(1:ndims)=1 ; vstrid(1:ndims)=1
! 
      ierr = NF90_INQUIRE_VARIABLE(idf,idv,dimids=idims)
      IF (ierr.NE.0) GOTO 103
!
      idvx = 0 ; idvy = 0 ; idvz = 0 ; idvt = 0
      DO idd = 1,ndims
        IF (idims(idd).EQ.idx) idvx = idd
        IF (idims(idd).EQ.idy) idvy = idd
        IF (idims(idd).EQ.idz) idvz = idd
        IF (idims(idd).EQ.idt) idvt = idd
      ENDDO
!
      IF ( (sxy_dim.GE.1).AND.(idvx.EQ.0) ) GOTO 105
      IF ( (sxy_dim.GE.2).AND.(idvy.EQ.0) ) GOTO 105
      IF (zdim.NE.'none') THEN
        IF ( (sxy_dim.GE.3).AND.(idvz.EQ.0) ) GOTO 105
      ENDIF
      IF ( (sxy_dim.GE.4).AND.(idvt.EQ.0) ) GOTO 105
!
      DO idd = 1,ndims
        IF (idd.EQ.idvx) vcount(idd)=sxy_jpi
        IF (idd.EQ.idvy) vcount(idd)=sxy_jpj
      ENDDO
!
! -2.- Read vertical levels (for Vx object only)
! ----------------------------------------------
! Separate this operations from the reading of the mask
! Create for this a new routine readlev, only called in obsv module
!
      IF (zdim.NE.'none') THEN
      IF ((kflagxyo.EQ.2).AND.(sxy_dim.GE.3)) THEN
! allocation lev
        allocate ( lev(1:nz), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        lev(:) = FREAL4(0.0)
!
!       Read vertical levels from files
        ierr = NF90_INQ_VARID(idf,sxyznam,idvzz)
        IF (ierr.NE.0) THEN
!
          IF (sxydmsk.EQ.2) THEN
            var_lev(1:sxy_jpk,indsxy) = FREAL(0.0)
          ELSEIF (sxydmsk.EQ.sxy_dim) THEN
            GOTO 103
          ENDIF
!
        ELSE
!
          allocate(zstart(1),zcount(1))
          zstart(1)=1 ; zcount(1)=nz
          ierr = NF90_GET_VAR(idf,idvzz,lev,start=zstart,count=zcount)
          deallocate(zstart,zcount)
          IF (ierr.NE.0) GOTO 104
!
          IF (sxydmsk.EQ.2) THEN
            var_lev(1:sxy_jpk,indsxy) = FREAL(lev(1))
          ELSEIF (sxydmsk.EQ.sxy_dim) THEN
            var_lev(1:jpkend,indsxy)=FREAL(lev(1:jpkend))
          ENDIF
!
        ENDIF
!
        IF (allocated(lev)) deallocate(lev)
!
      ENDIF
      ENDIF
!
! -3.- Read CDF mask arrays
! -------------------------
!
! --- allocation ptab
      allocate ( ptab(1:jpiend,1:jpjend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptab(:,:) = FREAL4(0.0)
!
      js=1
      DO jt=1,jptend
      DO jk=1,jpkend
!
         IF (idvz.NE.0) vstart(idvz)=jk
         IF (idvt.NE.0) vstart(idvt)=jt
!
         ierr = NF90_GET_VAR(idf,idv,ptab,start=vstart,count=vcount, &
     &                       stride=vstrid)
         IF (ierr.NE.0) GOTO 104
!
         DO jj=1,jpjend
         DO ji=1,jpiend
            IF (sxymsea.EQV.(sxyvmsk.EQ.ptab(ji,jj))) THEN
               mask(ji,jj,jk,jt)=mask(ji,jj,jk,jt)+(2**(ind_msk))
               js=js+1
            ENDIF
         ENDDO
         ENDDO
!
      ENDDO
      ENDDO
      sxynbr=js-1
!
      ierr = NF90_CLOSE(idf)
      IF (ierr.NE.0) GOTO 106
!
      IF (allocated(ptab)) deallocate(ptab)
!
! -5.- Control Print :
! --------------------
!
      IF (nprint.GE.3) THEN 
         kform='(8x,a,5i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndims
         kform='(8x,a,e12.3)'
         IF (  (jpkend.NE.nz).OR.(jptend.NE.nt) ) THEN
            kform='(8x,a,4i5)'
            WRITE(numout,kform) '- Size of loaded array : ',jpiend, &
     &                                  jpjend,jpkend,jptend
         ENDIF
         kform='(8x,a,i9)'
         WRITE(numout,kform) '- Number of loaded variables : ',sxynbr
         kform='(8x,a)'
         WRITE(numout,kform) '- Vertical levels : '
         kform='(8x,7e12.3)'
         WRITE(numout,kform) var_lev(1:jpkend,indsxy)
      ENDIF
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'lionc','readmsknc')
 1001 CALL printerror2(0,1001,3,'lionc','readmsknc')
!
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
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE lionc
