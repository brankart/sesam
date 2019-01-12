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
! ---                  LIOCDF.F90                                 ---
! ---                                                           ---
! --- original     : 99-03 (J.M. Brankart)                      ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 03-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  readcdf       : Read variable field from NetCDF file
! --- SUBROUTINE  writencdf     : Write variable field in NetCDF file
! --- SUBROUTINE  evalhdrmskcdf : Get variable dimensions from mask files
! --- SUBROUTINE  readmskcdf    : Read mask arrays from NetCDF mask files
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liocdf
      use mod_main
      use hiogrd
      use utilvct
      use utilcdfvar
      IMPLICIT NONE
      PRIVATE

      PUBLIC readcdf,writencdf,evalhdrmskcdf,readmskcdf

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcdf(kfname,kvectsout,kjsxy,ksomsxynbr,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Read variable field (of Vx or Vy object) from NetCDF file
!  -------
!  Method : Open, read and close NetCDF file (.cdf or .cdta)
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
      use mod_cfgxyo
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
      BIGREAL4, allocatable, dimension(:,:,:) :: ptabijk
      INTEGER :: indsxy,js,sompartsxynbr,jextdta,jextvar
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt
      CHARACTER(len=varlg) :: sxy_nam
      CHARACTER(len=bgword), dimension(1:4) :: dimunit
      CHARACTER(len=bgword) :: ktitle, kform
!
      INTEGER :: nx, ny, nz, nt
      INTEGER :: ji, jj, jk, jt
      INTEGER :: jpiend, jpjend, jpkend, jptend
      LOGICAL :: varfile
!----------------------------------------------------------------------
! Size of output vector
      jpssize=size(kvectsout,1)
!
! Get variable characteristics from SESAM configuration
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
         IF (jextdta.EQ.3) THEN
            varfile=.FALSE.
         ELSEIF (jextvar.EQ.3) THEN
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
         sxy_nam = var_nam(indsxy)
         sxy_dim = var_dim(indsxy)
         IF (sxy_dim.GE.1) sxy_jpi = var_jpi(indsxy)
         IF (sxy_dim.GE.2) sxy_jpj = var_jpj(indsxy)
         IF (sxy_dim.GE.3) sxy_jpk = var_jpk(indsxy)
         IF (sxy_dim.GE.4) sxy_jpt = var_jpt(indsxy)
      ELSE
! --- read Vy variable field from dta file
         sxy_nam = dta_nam(indsxy)
         sxy_dim = dta_dim(indsxy)
         IF (sxy_dim.GE.1) sxy_jpi = dta_jpi(indsxy)
         IF (sxy_dim.GE.2) sxy_jpj = dta_jpj(indsxy)
         IF (sxy_dim.GE.3) sxy_jpk = dta_jpk(indsxy)
         IF (sxy_dim.GE.4) sxy_jpt = dta_jpt(indsxy)
      ENDIF
      IF (sxy_nam.EQ.'T') THEN
         sxy_nam='TEM'
      ENDIF
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readfil/readcdf'
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
! -1.- Read NetCDF file header
! ----------------------------
!
      CALL cdfrdim (kfname,nx,ny,nz,nt,ktitle)
!
      IF (   (sxy_jpi.NE.nx) &
     &    .OR.(sxy_jpj.NE.ny) &
     &    .OR.(sxy_jpk.GT.nz) &
     &    .OR.(sxy_jpt.GT.nt) ) GOTO 102
!
! -2.- Read 2D arrays in NetCDF file
! ----------------------------------
!
      SELECT CASE (sxy_dim)
      CASE (2,3)
! --- allocation ptabij
        allocate ( ptabij(1:sxy_jpi,1:sxy_jpj), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        ptabij(:,:) = FREAL4(0.0)
      CASE (4)
! --- allocation ptabij
        allocate ( ptabijk(1:sxy_jpi,1:sxy_jpj,1:sxy_jpk), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        ptabij(:,:) = FREAL4(0.0)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      js=1
      SELECT CASE (sxy_dim)
      CASE (2)
! ==>  2D
         jk=1
         jt=1
         IF (js.LE.jpssize) THEN
            CALL cdfrtab(kfname,sxy_nam,jt,ptabij)
            CALL mk4vct(kvectsout(js:),ptabij(:,:), &
     &             jk,jt,kjsxy,sompartsxynbr,kflagxyo)
            js  = js + sompartsxynbr
         ENDIF
      CASE (3)
! ==>  3D
         DO jt=1,jptend
         DO jk=1,jpkend
            IF (js.LE.jpssize) THEN
               CALL cdfrsli(kfname,sxy_nam,3,jk,jt,ptabij)
               CALL mk4vct(kvectsout(js:),ptabij(:,:), &
     &              jk,jt,kjsxy,sompartsxynbr,kflagxyo)
               js  = js + sompartsxynbr
            ENDIF
         ENDDO
         ENDDO
      CASE (4)
! ==>  4D
         DO jt=1,jptend
            IF (js.LE.jpssize) THEN
               CALL cdfrtab(kfname,sxy_nam,jt,ptabijk)
               DO jk=1,jpkend
                  CALL mk4vct(kvectsout(js:),ptabijk(:,:,jk), &
     &                 jk,jt,kjsxy,sompartsxynbr,kflagxyo)
                  js  = js + sompartsxynbr
               ENDDO
            ENDIF
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      ksomsxynbr=js-1
!
      IF (allocated(ptabij)) deallocate(ptabij)
      IF (allocated(ptabijk)) deallocate(ptabijk)
!
! -3.- Control Print
! ------------------
!
      IF (nprint.GE.3) THEN
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title: ',ktitle(1:lenv(ktitle))
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
 1000 CALL printerror2(0,1000,1,'liocdf','readcdf')
 1001 CALL printerror2(0,1001,3,'liocdf','readcdf')
!
 102  WRITE (texterror,*) 'bad dimensions in NetCDF file: ',kfname
      CALL printerror2(0,102,3,'liocdf','readcdf',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writencdf(kfnoutncdf,ktitle,idast,kvectsin, &
     &        ktxy_indvct,ktxy_indmsk,txyend,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Write variable field (of Vx or Vy object) in NetCDF file
!  -------
!  Method : Open, write and close NetCDF file (.cdf or .cdta)
!  ------
!  Input : kfnoutncdf  : Filename
!  -----   ktitle : File title
!          idast  : date (obsolete)
!          kvectsin    : 1D SESAM vector object (Vx or Vy)
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
      use mod_cfgxyo
      use mod_coord
      use mod_spacexyo , only : spvaldta,spvalvar
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutncdf
      CHARACTER(len=*), intent(in) :: ktitle
      BIGREAL, intent(in) :: idast
      BIGREAL, dimension(:), intent(in) :: kvectsin
      INTEGER, dimension(:), intent(in) :: ktxy_indvct,ktxy_indmsk
      INTEGER, intent(in) :: txyend,kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: nx, ny, nz, nt
      BIGREAL4 :: lonmin, latmin, dx, dy, spval
      BIGREAL4 :: date
      INTEGER  ::allocok
      BIGREAL4, allocatable, dimension(:) :: lev,lon,lat
      BIGREAL4, allocatable, dimension(:,:) :: ptabij,lonxy,latxy
      BIGREAL4, allocatable, dimension(:,:,:) :: ptabijk
      INTEGER :: sxyend
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_dim, &
     &     sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxy_nbr,sxyngrd,sxyegrd
      INTEGER :: jsxy,indsxy,sompartxynbr
      INTEGER :: jtxy,txy_dimmax,indsxy_dimmax
      INTEGER, dimension(1:nbvar) :: txy_indtab
      CHARACTER(len=varlg), dimension(1:nbvar) :: sxy_nam
      CHARACTER(len=bgword), dimension(1:nbvar) :: sxyfgrd
      CHARACTER*100 unitnam,dimunit(4),llunit
      INTEGER :: ji, jj, jk, jt
      INTEGER :: jpiend, jpjend, jpkend, jptend
      LOGICAL :: existence
!----------------------------------------------------------------------
!
! Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         sxyend  = varend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=var_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
            sxy_nam(indsxy)=var_nam(indsxy)
            sxy_dim(indsxy)=var_dim(indsxy)
            sxy_jpi(indsxy)=var_jpi(indsxy)
            sxy_jpj(indsxy)=var_jpj(indsxy)
            sxy_jpk(indsxy)=var_jpk(indsxy)
            sxy_jpt(indsxy)=var_jpt(indsxy)
            sxy_nbr(indsxy)=var_nbr(indsxy)
            sxyngrd(indsxy)=varngrd(indsxy)
            sxyfgrd(indsxy)=varfgrd(indsxy)
            sxyegrd(indsxy)=varegrd(indsxy)
            IF (sxy_nam(indsxy).EQ.'T') THEN
               sxy_nam(indsxy)='TEM'
            ENDIF
         ENDDO
      CASE(2)
! --- dta
         sxyend  = dtaend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=dta_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
            sxy_nam(indsxy)=dta_nam(indsxy)
            sxy_dim(indsxy)=dta_dim(indsxy)
            sxy_jpi(indsxy)=dta_jpi(indsxy)
            sxy_jpj(indsxy)=dta_jpj(indsxy)
            sxy_jpk(indsxy)=dta_jpk(indsxy)
            sxy_jpt(indsxy)=dta_jpt(indsxy)
            sxy_nbr(indsxy)=dta_nbr(indsxy)
            sxyngrd(indsxy)=dtangrd(indsxy)
            sxyfgrd(indsxy)=dtafgrd(indsxy)
            sxyegrd(indsxy)=dtaegrd(indsxy)
            IF (sxy_nam(indsxy).EQ.'T') THEN
               sxy_nam(indsxy)='TEM'
            ENDIF
         ENDDO
      CASE(3)
         GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Compute dimensions of NetCDF file
      jtxy=1
      jsxy=ktxy_indmsk(jtxy)
      indsxy=sxy_ord(jsxy)
!
      jpiend = sxy_jpi(indsxy)
      jpjend = sxy_jpj(indsxy)
      jpkend = sxy_jpk(indsxy)        
      jptend = sxy_jpt(indsxy)       
!
      DO jtxy=2,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
!
         jpiend = max(jpiend,sxy_jpi(indsxy))
         jpjend = max(jpjend,sxy_jpj(indsxy))
         jpkend = max(jpkend,sxy_jpk(indsxy))
         jptend = max(jptend,sxy_jpt(indsxy))
      ENDDO
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writenfil/writencdf'
         WRITE(numout,*) '    ==> WRITING file: ', &
     &                            kfnoutncdf(1:lenv(kfnoutncdf))
      ENDIF
!
! allocate all necessary arrays
! --- allocation lev
      allocate ( lev(1:jpkend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lev(:) = FREAL4(0.0)
!
! --- allocation levk
      allocate ( levk(1:jpkend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      levk(:) = FREAL8(0.0)
!
! --- allocation lon
      allocate ( lon(1:jpiend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lon(:) = FREAL4(0.0)
!
! --- allocation longi
      allocate ( longi(1:jpiend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      longi(:) = FREAL8(0.0)
!
! --- allocation lat
      allocate ( lat(1:jpjend) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lat(:) = FREAL4(0.0)
!
! --- allocation latj
      allocate ( latj(1:jpjend) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      latj(:) = FREAL8(0.0)
!
! --- allocation gridij
      allocate ( gridij(1:jpiend,1:jpjend) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      gridij(:,:) = type_gridij(FREAL(0.0),FREAL(0.0))
!
! --- allocation ptabij
      allocate ( ptabij(1:jpiend,1:jpjend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL4(0.0)
!
! --- allocation ptabijk
      allocate ( ptabijk(1:jpiend,1:jpjend,1:jpkend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabijk(:,:,:) = FREAL4(0.0)
!
! --- allocation lonxy
      allocate ( lonxy(1:jpiend,1:jpjend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lonxy(:,:) = FREAL4(0.0)
!
! --- allocation latxy
      allocate ( latxy(1:jpiend,1:jpjend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      latxy(:,:) = FREAL4(0.0)
!
! -1.- Set parameters to write in NetCDF file
! -------------------------------------------
! dimensions
      nx     = jpiend
      ny     = jpjend
      nz     = jpkend
      nt     = jptend
      lonmin = FREAL4(1.)
      latmin = FREAL4(1.)
      dx     = FREAL4(1.)
      dy     = FREAL4(1.)
! special value
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         spval  = FREAL4(spvalvar)
      CASE(2)
! --- dta
         spval  = FREAL4(spvaldta)
      CASE(3)
         GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
! units
      unitnam = 'none'
      dimunit(1) = 'degrees_east'
      dimunit(2) = 'degrees_north'
      dimunit(3) = 'meters'
      dimunit(4) = 'days'
      llunit = 'degrees'
! maximum variable dimension
      txy_dimmax = 0
      DO jtxy=1,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy= sxy_ord(jsxy)
         IF (sxy_dim(indsxy).GT.txy_dimmax) THEN
           txy_dimmax=sxy_dim(indsxy)
           indsxy_dimmax=indsxy
         ENDIF
      ENDDO
! default grid locations
      DO ji=1,jpiend
         longi(ji) = lonmin + (ji-1) * dx
      ENDDO
      DO jj=1,jpjend
         latj(jj) = latmin + (jj-1) * dy
      ENDDO
      date = FREAL4(0.)
!
      levk(1:jpkend) = FREAL4(var_lev(1:jpkend,indsxy_dimmax))
!
! read grid from grid file (if such file exists)
!
      existence=.FALSE.
      INQUIRE (FILE=sxyfgrd(indsxy),EXIST=existence)
      IF (existence) THEN
         CALL readgrd (kflagxyo,indsxy)
      ENDIF
      lon(1:jpiend) = FREAL4(longi(1:jpiend))
      lat(1:jpjend) = FREAL4(latj(1:jpjend))
      lev(1:jpkend) = FREAL4(levk(1:jpkend))
      IF (existence.AND.(sxyngrd(1).GE.3)) THEN
        lonxy(1:jpiend,1:jpjend) = FREAL4(gridij(1:jpiend,1:jpjend)%longi)
        latxy(1:jpiend,1:jpjend) = FREAL4(gridij(1:jpiend,1:jpjend)%latj)
      ENDIF
!
! -2.- Write NetCDF file header
! -----------------------------
!
      CALL cdfwdim(kfnoutncdf,nx,ny,nz,ktitle)
      CALL cdfwloc(kfnoutncdf,lon,lat,lev,dimunit)
      IF (existence.AND.(sxyngrd(1).GE.3)) THEN
         CALL cdfwpos(kfnoutncdf,lonxy,latxy,spval,llunit)
      ENDIF
!
! -3.- Create variables
! ---------------------
!
      DO jtxy=1,txyend
         txy_indtab(jtxy)=ktxy_indvct(jtxy)
      ENDDO
!
      DO jtxy=1,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy= sxy_ord(jsxy)
         CALL cdfwvar(kfnoutncdf,sxy_nam(indsxy),spval, &
     &        unitnam,sxy_nam(indsxy),sxy_dim(indsxy))
         IF (nprint.GE.2) THEN 
            WRITE(numout,*) '    ==> WRITING variable: ', &
     &                                    sxy_nam(indsxy)
         ENDIF
      ENDDO
!
! -4.- Write 2D arrays in NetCDF file
! -----------------------------------
!
      DO jtxy=1,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         SELECT CASE (sxy_dim(indsxy))
         CASE (2)
! ==>  2D
            jk=1
            jt=1
            sompartxynbr=0
            IF (txy_indtab(jtxy).LE.size(kvectsin,1)) THEN
               CALL unmk4vct(kvectsin(txy_indtab(jtxy):), &
     &              ptabij(:,:),jk,jt,jsxy, &
     &              sompartxynbr,spval,kflagxyo) 
            ELSE
               ptabij(:,:) = spval
            ENDIF
            CALL cdfwtab(kfnoutncdf,sxy_nam(indsxy), &
     &           jt,date,ptabij)
            txy_indtab(jtxy)=txy_indtab(jtxy)+sompartxynbr
         CASE (3)
! ==>  3D
            DO jt=1,sxy_jpt(indsxy)
            DO jk=1,sxy_jpk(indsxy)
               sompartxynbr=0
               IF (txy_indtab(jtxy).LE.size(kvectsin,1)) THEN
                  CALL unmk4vct(kvectsin(txy_indtab(jtxy):), &
     &                 ptabij(:,:),jk,jt,jsxy, &
     &                 sompartxynbr,spval,kflagxyo) 
               ELSE
                  ptabij(:,:) = spval
               ENDIF
               CALL cdfwsli(kfnoutncdf,sxy_nam(indsxy), &
     &              3,jk,jt,date,ptabij)
               txy_indtab(jtxy)=txy_indtab(jtxy)+sompartxynbr
            ENDDO
            ENDDO
         CASE (4)
! ==>  4D
            DO jt=1,sxy_jpt(indsxy)
               DO jk=1,sxy_jpk(indsxy)
                  sompartxynbr=0
                  IF (txy_indtab(jtxy).LE.size(kvectsin,1)) THEN
                     CALL unmk4vct(kvectsin(txy_indtab(jtxy):), &
     &                    ptabijk(:,:,jk),jk,jt,jsxy, &
     &                    sompartxynbr,spval,kflagxyo)
                  ELSE
                     ptabijk(:,:,jk) = spval
                  ENDIF
                  txy_indtab(jtxy)=txy_indtab(jtxy)+sompartxynbr
               ENDDO
               CALL cdfwtab(kfnoutncdf,sxy_nam(indsxy), &
     &              jt,date,ptabijk)
            ENDDO
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
      ENDDO
!
! coherence test
      DO jtxy=1,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         IF (sxy_nbr(indsxy).NE.(txy_indtab(jtxy)-ktxy_indvct(jtxy))) &
     &                   GOTO 1000
      ENDDO
!
! --- deallocation
      IF (allocated(levk)) deallocate(levk)
      IF (allocated(lev)) deallocate(lev)
      IF (allocated(longi))  deallocate(longi)
      IF (allocated(lon))  deallocate(lon)
      IF (allocated(latj)) deallocate (latj)
      IF (allocated(lat)) deallocate (lat)
      IF (allocated(lonxy)) deallocate (lonxy)
      IF (allocated(latxy)) deallocate (latxy)
      IF (allocated(ptabij)) deallocate (ptabij)
      IF (allocated(ptabijk)) deallocate (ptabijk)
      IF (allocated(gridij)) deallocate (gridij)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'liocdf','writencdf')
 1001 CALL printerror2(0,1001,3,'liocdf','writencdf')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrmskcdf(kfname,kjpi,kjpj,kjpk,kjpt)
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
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfname
      INTEGER, intent(out) :: kjpi,kjpj,kjpk,kjpt
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=word80) :: ktitle, kform
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN 
         WRITE(numout,*)  &
     &    '*** ROUTINE : sesam/readmsk/evalhdrmsk/evalhdrmskcdf'
         WRITE(numout,*) '    ==> READING ',kfname(1:lenv(kfname))
      ENDIF
!
! -1.- Read NetCDF file header
! ----------------------------
!
      CALL cdfrdim (kfname,kjpi,kjpj,kjpk,kjpt,ktitle)
!
! -2.- Control Print
! ------------------
!
      IF (nprint.GE.3) THEN 
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title: ',ktitle(1:lenv(ktitle))
         kform='(8x,a,4i5)'
         WRITE(numout,kform) '- Dimensions: ',kjpi,kjpj,kjpk,kjpt
      ENDIF
!
      RETURN
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmskcdf(jsxy,kflagxyo)
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
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(in) :: jsxy,kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=word80) :: ktitle
      INTEGER :: nx, ny, nz, nt, ndim
      BIGREAL4 :: spval, sxyvmsk
      BIGREAL4, dimension(:), allocatable :: lon,lat,lev
      BIGREAL4, dimension(:,:), allocatable :: ptab
      BIGREAL4, dimension(:,:,:), allocatable :: ptabk
      CHARACTER(len=varlg) :: sxy_nam
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxypmsk,sxydmsk 
      INTEGER :: indsxy, ind_msk, js, sxynbr, allocok
      LOGICAL :: sxymsea
      CHARACTER(len=bgword) :: varunit, longname, kform, sxyfmsk
      CHARACTER(len=bgword), dimension(1:4) :: dimunit
      INTEGER :: ji, jj, jk, jt, jpiend, jpjend, jpkend, jptend
!----------------------------------------------------------------------
!
! Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         indsxy  = var_ord(jsxy)
         sxy_nam = var_nam(indsxy)
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
         ind_msk = jsxy-1
      CASE(2)
! --- dta
         indsxy  = dta_ord(jsxy)
         sxy_nam = dta_nam(indsxy)
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
         ind_msk = jsxy-1+varend
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Control print
      IF (nprint.GE.2) THEN 
         WRITE(numout,*) '*** ROUTINE : sesam/readmsk/readmskcdf'
         WRITE(numout,*) '    ==> READING ',sxyfmsk(1:lenv(sxyfmsk))
      ENDIF
!
! Set size of array to read
      SELECT CASE (sxydmsk)
      CASE (1)
! ==>  1D
         jpiend = sxy_jpi
         jpjend = 1
         jpkend = 1
         jptend = 1
      CASE (2)
! ==>  2D
         jpiend =  sxy_jpi
         jpjend =  sxy_jpj
         jpkend = 1
         jptend = 1
      CASE (3)
! ==>  3D
         jpiend =  sxy_jpi
         jpjend =  sxy_jpj
         jpkend =  sxy_jpk     
         jptend = 1
      CASE (4)
! ==>  4D
         jpiend =  sxy_jpi
         jpjend =  sxy_jpj
         jpkend =  sxy_jpk    
         jptend =  sxy_jpt
      CASE DEFAULT
! ==> ERROR
         GOTO 102
      END SELECT
!
! -1.- Read NetCDF file header
! ----------------------------
!
      CALL cdfrdim(sxyfmsk,nx,ny,nz,nt,ktitle)
      CALL cdfrvar(sxyfmsk,sxy_nam,spval,varunit,longname,ndim)
!
      IF (    (jpiend.NE.nx) &
     &    .OR.(jpjend.NE.ny) &
     &    .OR.(jpkend.GT.nz) &
     &    .OR.(jptend.GT.nt) ) GOTO 102
!
      IF (size(mask,3).GT.nz) GOTO 102
!
! -2.- Set special value if optional mask is used
! -----------------------------------------------
!
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var 
         IF (largvarmsk) THEN
            indsxy  = var_ord(jsxy)
            varvmsk(indsxy) = spval
            sxyvmsk = FREAL4(varvmsk(indsxy))
         ENDIF
      CASE(2)
! --- dta
         IF (largdtamsk) THEN
            indsxy  = dta_ord(jsxy)
            dtavmsk(indsxy) = spval
            sxyvmsk = FREAL4(dtavmsk(indsxy))
         ENDIF
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -3.- Read vertical levels
! -------------------------
! --- allocation lon
      allocate ( lon(1:nx), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lon(:) = FREAL4(0.0)
! --- allocation lat
      allocate ( lat(1:ny), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lat(:) = FREAL4(0.0)
! --- allocation lev
      allocate ( lev(1:nz), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lev(:) = FREAL4(0.0)
!
      CALL cdfrloc (sxyfmsk,lon,lat,lev,dimunit)
!
      IF (allocated(lon))  deallocate(lon)
      IF (allocated(lat)) deallocate (lat)
!
      IF (kflagxyo.EQ.1) THEN
         IF (sxydmsk.EQ.2) THEN
            indsxy=var_ord(jsxy)
            var_lev(1:sxy_jpk,indsxy) = FREAL(lev(1))
         ELSEIF (sxydmsk.EQ.sxy_dim) THEN
            indsxy=var_ord(jsxy)
            var_lev(1:jpkend,indsxy)=FREAL(lev(1:jpkend))
         ELSE
            IF (sxydmsk.NE.sxy_dim) GOTO 103
         ENDIF
      ENDIF
!
      IF (allocated(lev)) deallocate(lev)
!
! -4.- Read CDF mask arrays
! -------------------------
!
      allocate ( ptabk(1:jpiend,1:jpjend,1:jpkend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabk(:,:,:) = FREAL4(0.0)
!
      js=1
      DO jt=1,jptend
         CALL cdfrtab(sxyfmsk,sxy_nam,jt,ptabk)
         DO jk=1,jpkend
            DO jj=1,jpjend
               DO ji=1,jpiend
                  IF (sxymsea.EQV.(sxyvmsk.EQ.ptabk(ji,jj,jk))) THEN
                     mask(ji,jj,jk,jt)=mask(ji,jj,jk,jt)+(2**(ind_msk))
                     js=js+1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      sxynbr=js-1
!
      IF (allocated(ptabk)) deallocate(ptabk)
!
! -5.- Control Print :
! --------------------
!
      IF (nprint.GE.3) THEN 
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title: ',ktitle(1:lenv(ktitle))
         kform='(8x,a,5i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
         kform='(8x,a,e12.3)'
         WRITE(numout,kform) '- Special value: ',spval
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
 1000 CALL printerror2(0,1000,1,'liocdf','readmskcdf')
 1001 CALL printerror2(0,1001,3,'liocdf','readmskcdf')
!
 102  WRITE (texterror,*) 'bad dimensions in mask file: ',sxyfmsk
      CALL printerror2(0,102,3,'liocdf','readmskcdf', &
     &     comment=texterror)
 103  WRITE (texterror,*) 'incoherence between mask files and ', &
     &                    'SESAM configuration'
      CALL printerror2(0,103,3,'liocdf','readmskcdf', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liocdf
