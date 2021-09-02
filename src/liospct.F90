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
! ---                  LIOSPCT.F90                                ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 16-06  (J.M. Brankart)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  readspct      : Read spectrum from NetCDF file
! --- SUBROUTINE  writespct     : Write spectrum in NetCDF file
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liospct
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC readspct,writespct

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readspct(kinspct,kspct,kflagxyo,kjsxy,kjk,kjt)
!---------------------------------------------------------------------
!
!  Purpose : Read spectrum from NetCDF file
!  -------
!  Input :  kinspct   : input file
!  -----    kflagxyo  : vector type (Vx,Vy)
!           kjsxy     : variable index
!           kjk       : horizontal slice index
!           kjt       : time slice index
!  Output : kspct     : spectrum to read
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo
      use utilcdfvar
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinspct
      BIGREAL, dimension(:,:), intent(out) :: kspct
      INTEGER, intent(in) :: kflagxyo,kjsxy,kjk,kjt
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpisize,jpjsize,jpksize,jptsize
      INTEGER :: jpi1,jpj1,jpk1,jpt1
      INTEGER :: indsxy,jsxy,sxy_dim
      CHARACTER(len=varlg) :: sxy_nam
      INTEGER :: allocok
      LOGICAL :: incompatible
      BIGREAL4, allocatable, dimension(:,:) :: ptabij
      BIGREAL :: spval
      CHARACTER(len=word80) :: cdrec1
!----------------------------------------------------------------------
!
! Check input arguments
      jpisize = size(kspct,1)
      jpjsize = size(kspct,2)
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../algospct/readspct'
         WRITE(numout,*) '    ==> READING file ',kinspct(1:lenv(kinspct))
      ENDIF
!
! Get characteristics of variables
      jpksize=0 ; jptsize=0
      SELECT CASE (kflagxyo)
      CASE(1)
        indsxy=var_ord(kjsxy)
        sxy_nam=var_nam(indsxy)
        sxy_dim=var_dim(indsxy)
        DO jsxy = 1,varend
          indsxy  = var_ord(jsxy)
          jpksize = MAX(jpksize,var_jpk(indsxy))
          jptsize = MAX(jptsize,var_jpt(indsxy))
        ENDDO
        spval=spvalvar
      CASE(2,3)
        indsxy=dta_ord(kjsxy)
        sxy_nam=dta_nam(indsxy)
        sxy_dim=dta_dim(indsxy)
        DO jsxy = 1,dtaend
          indsxy  = dta_ord(jsxy)
          jpksize = MAX(jpksize,dta_jpk(indsxy))
          jptsize = MAX(jptsize,dta_jpt(indsxy))
        ENDDO
        spval=spvaldta
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -1.- Read and check file dimensions
! -----------------------------------
!
      CALL cdfrdim(kinspct,jpi1,jpj1,jpk1,jpt1,cdrec1)
      IF (jpi1.LT.jpisize) GOTO 102
      IF (jpj1.LT.jpjsize) GOTO 102
      IF (jpk1.NE.jpksize) GOTO 102
      IF (jpt1.NE.jptsize) GOTO 102
!
! -2.- Read spectrum from file
! ----------------------------
!
! --- allocation ptabij
      allocate ( ptabij(1:jpi1,1:jpj1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
      SELECT CASE (sxy_dim)
      CASE (2)
         CALL cdfrtab(kinspct,sxy_nam,kjt,ptabij)
      CASE (3,4)
         CALL cdfrsli(kinspct,sxy_nam,3,kjk,kjt,ptabij)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      kspct(:,:) = ptabij(1:jpisize,jpj1/2+1-jpjsize/2:jpj1/2+1+jpjsize/2)
      WHERE(kspct==spval) kspct=0.0
!
! --- deallocation
      IF (allocated(ptabij)) deallocate(ptabij)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liospct','readspct')
 1001 CALL printerror2(0,1001,1,'liospct','readspct')
!
 102  WRITE (texterror,*) 'Incompatible spectrum file'
      CALL printerror2(0,102,3,'liospct','readspct',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writespct(koutspct,kspct,kflagxyo,kjsxy,kjk,kjt)
!---------------------------------------------------------------------
!
!  Purpose : Write spectrum in NetCDF file
!  -------
!  Input : koutspct   : output file
!  -----   kspct      : spectrum to write
!          kflagxyo   : vector type (Vx,Vy)
!          kjsxy      : variable index
!          kjk        : horizontal slice index
!          kjt        : time slice index
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
      CHARACTER(len=*), intent(in) :: koutspct
      BIGREAL, dimension(:,:), intent(in) :: kspct
      INTEGER, intent(in) :: kflagxyo,kjsxy,kjk,kjt
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpisize,jpjsize,jpksize,jptsize
      INTEGER :: jpi1,jpj1,jpk1,jpt1,ji,jj,jk
      INTEGER :: jsxy, indsxy, sxyend, indsxy_dimmax
      INTEGER, dimension(1:nbvar) :: sxy_ord, sxy_dim
      CHARACTER(len=varlg), dimension(1:nbvar) :: sxy_nam
      CHARACTER(len=varlg) :: sxy_nam1
      INTEGER :: sxy_dim1
      INTEGER allocok
      LOGICAL :: existence
      BIGREAL4, allocatable, dimension(:,:) :: ptabij
      BIGREAL4, allocatable, dimension(:) :: lon, lat, lev
      BIGREAL4 :: time, spval
      CHARACTER(len=word80) :: cdrec1
      CHARACTER(len=bgword), dimension(1:4) :: dimunit
      CHARACTER(len=bgword) :: unitnam
!----------------------------------------------------------------------
!
! Check input arguments
      jpisize = size(kspct,1)
      jpjsize = size(kspct,2)
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../algospct/writespct'
         WRITE(numout,*) '    ==> WRITING file ',koutspct(1:lenv(koutspct))
      ENDIF
!
! Get characteristics of variables
      SELECT CASE (kflagxyo)
      CASE(1)
        sxyend=varend
        indsxy_dimmax=1 ; jpksize=0 ; jptsize=0
        DO jsxy = 1,sxyend
          sxy_ord(jsxy)=var_ord(jsxy)
          indsxy=sxy_ord(jsxy)
          sxy_nam(indsxy)=var_nam(indsxy)
          sxy_dim(indsxy)=var_dim(indsxy)
          IF (sxy_dim(indsxy).GT.sxy_dim(indsxy_dimmax)) THEN
            indsxy_dimmax=indsxy
          ENDIF
          jpksize = MAX(jpksize,var_jpk(indsxy))
          jptsize = MAX(jptsize,var_jpt(indsxy))
        ENDDO
        indsxy=sxy_ord(kjsxy)
        sxy_nam1=sxy_nam(indsxy)
        sxy_dim1=sxy_dim(indsxy)
      CASE(2,3)
        sxyend=dtaend
        indsxy_dimmax=1 ; jpksize=0 ; jptsize=0
        DO jsxy = 1,sxyend
          sxy_ord(jsxy)=dta_ord(jsxy)
          indsxy= sxy_ord(jsxy)
          sxy_nam(indsxy)=dta_nam(indsxy)
          sxy_dim(indsxy)=dta_dim(indsxy)
          IF (sxy_dim(indsxy).GT.sxy_dim(indsxy_dimmax)) THEN
            indsxy_dimmax=indsxy
          ENDIF
          jpksize = MAX(jpksize,dta_jpk(indsxy))
          jptsize = MAX(jptsize,dta_jpt(indsxy))
        ENDDO
        indsxy=sxy_ord(kjsxy)
        sxy_nam1=sxy_nam(indsxy)
        sxy_dim1=sxy_dim(indsxy)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -1.- Check or write spectrum file dimensions
! --------------------------------------------
!
! Check if spectrum file already exists
      INQUIRE (FILE=koutspct,EXIST=existence)
!
      IF (existence) THEN
         CALL cdfrdim(koutspct,jpi1,jpj1,jpk1,jpt1,cdrec1)
         IF (jpi1.NE.jpisize) GOTO 102
         IF (jpj1.NE.jpjsize) GOTO 102
         IF (jpk1.NE.jpksize) GOTO 102
      ELSE
         cdrec1='SESAM spectrum file'
         CALL cdfwdim(koutspct,jpisize,jpjsize,jpksize,cdrec1)

         allocate ( lon(1:jpisize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         allocate ( lat(1:jpjsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         allocate ( lev(1:jpksize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001

         DO ji=1,jpisize
           lon(ji) = FREAL4(ji-1)
         ENDDO
         DO jj=1,jpjsize
           lat(jj) = FREAL4(jj-1-jpjsize/2)
         ENDDO
         DO jk=1,jpksize
           lev(jk) = FREAL4(var_lev(jk,indsxy_dimmax))
         ENDDO

         dimunit(1) = 'none'
         dimunit(2) = 'none'
         dimunit(3) = 'meters'
         dimunit(4) = 'days'

         CALL cdfwloc(koutspct,lon,lat,lev,dimunit)

         unitnam = 'none'
         spval = FREAL4(-9999.)
         DO jsxy = 1,sxyend
           indsxy=sxy_ord(jsxy)
           CALL cdfwvar(koutspct,sxy_nam(indsxy),spval, &
     &                  unitnam,sxy_nam(indsxy),sxy_dim(indsxy))
         ENDDO
      ENDIF
!
! -2.- Write spectrum in file
! ---------------------------
!
! --- allocation ptabij
      allocate ( ptabij(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
      ptabij(:,:) = kspct(:,:) 
!
      time = 0.0
      SELECT CASE (sxy_dim1)
      CASE (2)
         CALL cdfwtab(koutspct,sxy_nam1,kjt,time,ptabij)
      CASE (3,4)
         CALL cdfwsli(koutspct,sxy_nam1,3,kjk,kjt,time,ptabij)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! --- deallocation
      IF (allocated(ptabij)) deallocate(ptabij)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liospct','writespct')
 1001 CALL printerror2(0,1001,3,'liospct','writespct')
!
 102  WRITE (texterror,*) 'inconsistent dimensions in spectrum file:', &
     &                        koutspct(1:lenv(koutspct))
      CALL printerror2(0,102,3,'liospct','writespct',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liospct
