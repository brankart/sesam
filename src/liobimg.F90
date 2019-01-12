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
! ---                  LIOBIMG.F90                                ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12  (C.E. Testut)                       ---
! --- modification : 99-05  (C.E. Testut)                       ---
! --- modification : 01-06  (C.E. Testut)                       ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE readbimg       : Read variable field from BIMG file
! --- SUBROUTINE writenbimg     : Write variable field in BIMG file
! --- SUBROUTINE evalhdrmskbimg : Get variable dimensions from mask files
! --- SUBROUTINE readmskbimg    : Read mask arrays from BIMG mask files
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liobimg
      use mod_main
      use utilvct
      use utilfiles
      IMPLICIT NONE
      PRIVATE

      PUBLIC readbimg,writenbimg,evalhdrmskbimg,readmskbimg

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readbimg(kfname,kvectsout,kjsxy,ksomsxynbr,kflagxyo)
!-----------------------------------------------------------------
!
!  Purpose : Read variable field (of Vx or Vy object) from BIMG file
!  -------
!  Method : Open, read and close DIMG file (.dimg or .dta)
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
      CHARACTER(len=word80) :: cdrec1, cdrec2, cdrec3, cdrec4, kform
      INTEGER :: nx, ny, nz, nt, ndim, icod
      BIGREAL4 :: lonmin, latmin, dx, dy, spval
      BIGREAL4 :: date
      INTEGER :: allocok,jpssize
      BIGREAL4, allocatable, dimension(:) :: lev
      BIGREAL4, allocatable, dimension(:,:) :: ptabij
      INTEGER :: indsxy,js,sompartsxynbr,jextdta,jextvar
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,ksxypos
      INTEGER :: ji, jj, jk, jt, jdim
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
         varfile=.TRUE.
      CASE(2)
! --- dta
         indsxy  = dta_ord(kjsxy)
         jextvar=indext(kfname,extvartab,nbextvar)
         jextdta=indext(kfname,extdtatab,nbextdta)
         IF (jextdta.EQ.1) THEN
            varfile=.FALSE.
         ELSEIF (jextvar.EQ.1) THEN
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
         sxy_dim = var_dim(indsxy)
         IF (sxy_dim.GE.1) sxy_jpi = var_jpi(indsxy)
         IF (sxy_dim.GE.2) sxy_jpj = var_jpj(indsxy)
         IF (sxy_dim.GE.3) sxy_jpk = var_jpk(indsxy)
         IF (sxy_dim.GE.4) sxy_jpt = var_jpt(indsxy)
         ksxypos = mod(varipos(indsxy),10)
      ELSE
! --- read Vy variable field from dta file
         sxy_dim = dta_dim(indsxy)
         IF (sxy_dim.GE.1) sxy_jpi = dta_jpi(indsxy)
         IF (sxy_dim.GE.2) sxy_jpj = dta_jpj(indsxy)
         IF (sxy_dim.GE.3) sxy_jpk = dta_jpk(indsxy)
         IF (sxy_dim.GE.4) sxy_jpt = dta_jpt(indsxy)
         ksxypos = mod(dtaipos(indsxy),10)
      ENDIF
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readfil/readbimg'
         WRITE(numout,*) '    ==> READING file ',kfname(1:lenv(kfname))
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
! --- allocation lev
      allocate ( lev(1:sxy_jpk), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lev(:) = FREAL4(0.0)
!
! --- allocation ptabij
      allocate ( ptabij(1:sxy_jpi,1:sxy_jpj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL4(0.0)
!
! -1.- Open BIMG file
! -------------------
!
      CALL openfile(numfil,kfname,kform=clunf)
!
! -2.- Read BIMG file header
! --------------------------
!
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec1
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec2
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec3
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec4
!
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) nx, ny, nz, nt, ndim, icod
!
      IF (   (sxy_jpi.NE.nx) &
     &    .OR.(sxy_jpj.NE.ny) &
     &    .OR.(sxy_jpk.GT.nz) &
     &    .OR.(sxy_jpt.GT.nt) &
     &    .OR.(ksxypos.GT.ndim) ) GOTO 102
!
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) lonmin,latmin,dx,dy,spval
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) lev(1:sxy_jpk)
!
! -3.- Read 2D arrays in BIMG file
! --------------------------------
!
      js=1
      ksomsxynbr=0
      sompartsxynbr=0
      DO jt=1,sxy_jpt
         READ(UNIT=numfil,ERR=101,IOSTAT=iost) date
!         
         IF (jt.LE.jptend) THEN
            DO jk=1,sxy_jpk
               DO jdim=1,(ksxypos-1)
                  READ(UNIT=numfil,ERR=101,IOSTAT=iost)
               ENDDO
               IF (jk.LE.jpkend) THEN
                  IF (js.LE.jpssize) THEN
                     READ(UNIT=numfil,ERR=101,IOSTAT=iost) ptabij(:,:)
                     CALL mk4vct(kvectsout(js:),ptabij(:,:), &
     &                    jk,jt,kjsxy,sompartsxynbr,kflagxyo)
                  ELSE
                     sompartsxynbr=0
                  ENDIF
                  js  = js + sompartsxynbr
                  ksomsxynbr= ksomsxynbr + sompartsxynbr
                  DO jdim=(ksxypos+1),ndim
                     READ(UNIT=numfil,ERR=101,IOSTAT=iost)
                  ENDDO
               ELSE
                  DO jdim=(ksxypos),ndim
                     READ(UNIT=numfil,ERR=101,IOSTAT=iost)
                  ENDDO               
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!
! -3.- Close BIMG file 
! --------------------
!
      CLOSE (UNIT=numfil)
!
! -4.- Control Print
! ------------------
!
      IF (nprint.GE.3) THEN 
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title1: ',cdrec1(1:lenv(cdrec1))
         WRITE(numout,kform) '- Title2: ',cdrec2(1:lenv(cdrec2))
         WRITE(numout,kform) '- Title3: ',cdrec3(1:lenv(cdrec3))
         WRITE(numout,kform) '- Title4: ',cdrec4(1:lenv(cdrec4))
         kform='(8x,a,5i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
         kform='(8x,a,e12.3)'
         WRITE(numout,kform) '- Special value: ',spval
         kform='(8x,a,4e12.3)'
         WRITE(numout,kform) '- Regular grid: ',lonmin,latmin,dx,dy
         IF ( (jpkend.NE.nz).OR.(jptend.NE.nt).OR.(ksxypos.NE.1) ) THEN
            kform='(8x,a,4i5)'
            WRITE(numout,kform) '- Size of loaded array : ',jpiend, &
     &                                  jpjend,jpkend,jptend
            kform='(8x,a,i5,a,i5)'
            WRITE(numout,kform) '- Variable : ',ksxypos,'/',ndim
         ENDIF
      ENDIF
!
! --- deallocation
      IF (allocated( lev)) deallocate(lev)
      IF (allocated( ptabij)) deallocate(ptabij)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liobimg','readbimg')
 1001 CALL printerror2(0,1001,3,'liobimg','readbimg')
!
 101  WRITE (texterror,*) 'error reading BIMG file, iost=',iost
      CALL printerror2(0,101,3,'liobimg','readbimg',comment=texterror)
 102  WRITE (texterror,*) 'bad dimensions in BIMG file: ',kfname
      CALL printerror2(0,102,3,'liobimg','readbimg',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writenbimg(kfnoutnbimg,lgcdrec1,lgcdrec2, &
     &     lgcdrec3,lgcdrec4,idast,kvectsin, &
     &     ktxy_indvct,ktxy_indmsk,txyend,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Write variable field (of Vx or Vy object) in BIMG file
!  -------
!  Method : Open, write and close DIMG file (.bimg or .bdta)
!  ------
!  Input : kfnoutnbimg : Filename
!  -----   lgcdrec1    : 1st record of BIMG file header
!          lgcdrec2    : 2nd record of BIMG file header
!          lgcdrec3    : 3rd record of BIMG file header
!          lgcdrec4    : 4th record of BIMG file header
!          idast       : date (obsolete)
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
      use mod_spacexyo , only : spvaldta,spvalvar
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutnbimg
      CHARACTER(len=word80), intent(in) :: lgcdrec1,lgcdrec2
      CHARACTER(len=word80), intent(in) :: lgcdrec3,lgcdrec4
      BIGREAL, intent(in) :: idast
      BIGREAL, dimension(:), intent(in) :: kvectsin
      INTEGER, dimension(:), intent(in) :: ktxy_indvct,ktxy_indmsk
      INTEGER, intent(in) :: txyend,kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: nx, ny, nz, nt, ndim, icod
      BIGREAL4 :: lonmin, latmin, dx, dy, spval
      BIGREAL4 ::  date
      CHARACTER(len=word80) :: cdrec1, cdrec2, cdrec3, cdrec4, kform
      INTEGER  ::allocok
      BIGREAL4, allocatable, dimension(:) :: lev
      BIGREAL4, allocatable, dimension(:,:) :: ptabij
      INTEGER :: sxyend
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_dim &
     &     ,sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxy_nbr
      INTEGER :: jsxy,indsxy,sompartxynbr
      INTEGER :: jtxy
      INTEGER, dimension(1:nbvar) :: txy_indtab
      INTEGER :: ji, jj, jk, jt
      INTEGER :: jpiend, jpjend, jpkend, jptend
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
            sxy_dim(indsxy)=var_dim(indsxy)
            sxy_jpi(indsxy)=var_jpi(indsxy)
            sxy_jpj(indsxy)=var_jpj(indsxy)
            sxy_jpk(indsxy)=var_jpk(indsxy)
            sxy_jpt(indsxy)=var_jpt(indsxy)
            sxy_nbr(indsxy)=var_nbr(indsxy)
         ENDDO
      CASE(2)
! --- dta
         sxyend  = dtaend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=dta_ord(jsxy)
            indsxy= sxy_ord(jsxy)
!
            sxy_dim(indsxy)=dta_dim(indsxy)
            sxy_jpi(indsxy)=dta_jpi(indsxy)
            sxy_jpj(indsxy)=dta_jpj(indsxy)
            sxy_jpk(indsxy)=dta_jpk(indsxy)
            sxy_jpt(indsxy)=dta_jpt(indsxy)
            sxy_nbr(indsxy)=dta_nbr(indsxy)
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writenfil/writenbimg'
         WRITE(numout,*) '    ==> WRITING file: ', &
     &                           kfnoutnbimg(1:lenv(kfnoutnbimg))
      ENDIF
!
! Compute dimensions of BIMG file
      jtxy=1
      jsxy=ktxy_indmsk(jtxy)
      indsxy=sxy_ord(jsxy)
      SELECT CASE (sxy_dim(indsxy))
      CASE (1)
! ==>  1D
         jpiend = sxy_jpi(indsxy)
         jpjend = 1
         jpkend = 1
         jptend = 1
         DO jtxy=2,txyend
            jsxy=ktxy_indmsk(jtxy)
            indsxy=sxy_ord(jsxy)
            IF (jpiend.NE.sxy_jpi(indsxy)) GOTO 102
         ENDDO
      CASE (2)
! ==>  2D
         jpiend = sxy_jpi(indsxy)
         jpjend = sxy_jpj(indsxy)
         jpkend = 1
         jptend = 1
         DO jtxy=2,txyend
            jsxy=ktxy_indmsk(jtxy)
            indsxy=sxy_ord(jsxy)
            IF (jpiend.NE.sxy_jpi(indsxy)) GOTO 102
            IF (jpjend.NE.sxy_jpj(indsxy)) GOTO 102
         ENDDO
      CASE (3)
! ==>  3D
         jpiend = sxy_jpi(indsxy)
         jpjend = sxy_jpj(indsxy)
         jpkend = sxy_jpk(indsxy)      
         jptend = 1
         DO jtxy=2,txyend
            jsxy=ktxy_indmsk(jtxy)
            indsxy=sxy_ord(jsxy)
            IF (jpiend.NE.sxy_jpi(indsxy)) GOTO 102
            IF (jpjend.NE.sxy_jpj(indsxy)) GOTO 102
            IF (jpkend.NE.sxy_jpk(indsxy)) GOTO 102
         ENDDO
      CASE (4)
! ==>  4D
         jpiend = sxy_jpi(indsxy)   
         jpjend = sxy_jpj(indsxy)   
         jpkend = sxy_jpk(indsxy)        
         jptend = sxy_jpt(indsxy)       
         DO jtxy=2,txyend
            jsxy=ktxy_indmsk(jtxy)
            indsxy=sxy_ord(jsxy)
            IF (jpiend.NE.sxy_jpi(indsxy)) GOTO 102
            IF (jpjend.NE.sxy_jpj(indsxy)) GOTO 102
            IF (jpkend.NE.sxy_jpk(indsxy)) GOTO 102
            IF (jptend.NE.sxy_jpt(indsxy)) GOTO 102
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! --- allocation lev
      allocate ( lev(1:jpkend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lev(:) = FREAL4(0.0)
!
! --- allocation ptabij
      allocate ( ptabij(1:jpiend,1:jpjend), stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      ptabij(:,:) = FREAL4(0.0)
!
! -1.- Set parameters to write in DIMG file
! -----------------------------------------
!
      nx     = jpiend
      ny     = jpjend
      nz     = jpkend
      nt     = jptend
      ndim   = txyend
      icod   = 1001
      lonmin = FREAL4(1.)
      latmin = FREAL4(1.)
      dx     = FREAL4(1.)
      dy     = FREAL4(1.)
      SELECT CASE(kflagxyo)
      CASE(1)
! --- var
         spval  = FREAL4(spvalvar)
      CASE(2)
! --- dta
         spval  = FREAL4(spvaldta)
      CASE DEFAULT
         GOTO 1000
      END SELECT
      lev(1:jpkend) =  FREAL4(var_lev(1:jpkend,indsxy))
      date   = FREAL4(idast)
!
      txy_indtab(1:txyend)=ktxy_indvct(1:txyend)
!
! -1.- Open BIMG file
! --------------------
!
      CALL openfile(numfil,kfnoutnbimg,kstatus=clunk,kform=clunf)
!
! -2.- Write BIMG file header
! ---------------------------
!
      REWIND(UNIT=numfil)
!
      WRITE(cdrec1,'(a)') lgcdrec1(1:80)
      WRITE(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec1
      WRITE(cdrec2,'(a)') lgcdrec2(1:80)
      WRITE(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec2
      WRITE(cdrec3,'(a)') lgcdrec3(1:80)
      WRITE(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec3
      WRITE(cdrec4,'(a)') lgcdrec4(1:80)
      WRITE(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec4
!
      WRITE(UNIT=numfil,ERR=101,IOSTAT=iost) nx, ny, nz, nt, ndim, icod
      WRITE(UNIT=numfil,ERR=101,IOSTAT=iost) lonmin,latmin,dx,dy,spval
      WRITE(UNIT=numfil,ERR=101,IOSTAT=iost) lev(1:jpkend)
!
! -3.- Write 2D arrays in DIMG file
! ---------------------------------
!
      DO jt=1,jptend
         IF (jptend.EQ.1) THEN
            WRITE(UNIT=numfil,ERR=101,IOSTAT=iost) date
         ELSE
            WRITE(UNIT=numfil,ERR=101,IOSTAT=iost) jt
         ENDIF
!
         DO jk=1,jpkend
            DO jtxy=1,txyend
               jsxy=ktxy_indmsk(jtxy)
               sompartxynbr=0
               IF (txy_indtab(jtxy).LE.size(kvectsin,1)) THEN
                  CALL unmk4vct(kvectsin(txy_indtab(jtxy):), &
     &                 ptabij(:,:),jk,jt,jsxy, &
     &                 sompartxynbr,spval,kflagxyo) 
               ELSE
                  ptabij(:,:) = spval
               ENDIF
               WRITE(UNIT=numfil,ERR=101,IOSTAT=iost) ptabij(:,:)
               txy_indtab(jtxy)=txy_indtab(jtxy)+sompartxynbr
            ENDDO
         ENDDO
      ENDDO
!
! coherence test
      DO jtxy=1,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         IF (sxy_nbr(indsxy).NE.(txy_indtab(jtxy)-ktxy_indvct(jtxy))) &
     &                   GOTO 101
      ENDDO
!
! -3.- Close BIMG file
! --------------------
!
      CLOSE (UNIT=numfil)
!
! -4.- Control Print
! ------------------
!
      IF (nprint.GE.3) THEN 
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title1: ',cdrec1(1:lenv(cdrec1))
         WRITE(numout,kform) '- Title2: ',cdrec2(1:lenv(cdrec2))
         WRITE(numout,kform) '- Title3: ',cdrec3(1:lenv(cdrec3))
         WRITE(numout,kform) '- Title4: ',cdrec4(1:lenv(cdrec4))
         kform='(8x,a,5i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
         kform='(8x,a,e12.3)'
         WRITE(numout,kform) '- Special value: ',spval
         kform='(8x,a,4e12.3)'
         WRITE(numout,kform) '- Regular grid: ',lonmin,latmin,dx,dy
      ENDIF
!
! --- deallocation
      IF (allocated(lev)) deallocate(lev)
      IF (allocated(ptabij)) deallocate (ptabij)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liobimg','writenbimg')
 1001 CALL printerror2(0,1001,3,'liobimg','writenbimg')
!
 101  WRITE (texterror,*) 'error writing BIMG file, iost=',iost
      CALL printerror2(0,101,3,'liobimg','writenbimg',comment=texterror)
 102  WRITE (texterror,*) 'inconsistent data dimensions'
      CALL printerror2(0,102,3,'liobimg','writenbimg',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrmskbimg(kfname,kjpi,kjpj,kjpk,kjpt)
!---------------------------------------------------------------------
!
!  Purpose : Get variable dimensions from mask files
!  -------
!  Method : Open, read and close BIMG file header
!  ------
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
      CHARACTER(len=word80) :: cdrec1, cdrec2, cdrec3, cdrec4, kform
      INTEGER :: nx, ny, nz, nt, ndim, icod
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/readmsk/evalhdrmsk/evalhdrmskbimg'
         WRITE(numout,*) '    ==> READING ',kfname(1:lenv(kfname))
      ENDIF
!
! -1.- Open BIMG file
! -------------------
!
      CALL openfile(numfil,kfname,kform=clunf)
!
! -2.- Read BIMG file header
! --------------------------
!
      REWIND(UNIT=numfil)
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec1
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec2
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec3
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec4
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) nx, ny, nz, nt, ndim, icod
!
! -3.- Close BIMG file
! --------------------
!
      CLOSE (UNIT=numfil)
!
! -4.- Control Print
! ------------------
!
      IF (nprint.GE.3) THEN 
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title1: ',cdrec1(1:lenv(cdrec1))
         WRITE(numout,kform) '- Title2: ',cdrec2(1:lenv(cdrec2))
         WRITE(numout,kform) '- Title3: ',cdrec3(1:lenv(cdrec3))
         WRITE(numout,kform) '- Title4: ',cdrec4(1:lenv(cdrec4))
         kform='(8x,a,5i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
      ENDIF
!
! -5.- Set output variables
! -------------------------
!
      kjpi=nx
      kjpj=ny
      kjpk=nz
      kjpt=nt
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liobimg','evalhdrmskbimg')
 1001 CALL printerror2(0,1001,3,'liobimg','evalhdrmskbimg')
!
 101  WRITE (texterror,*) 'error reading BIMG file, iost=',iost
      CALL printerror2(0,101,3,'liobimg','evalhdrmskbimg',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmskbimg(jsxy,kflagxyo)
!---------------------------------------------------------------------
!
!  Purpose : Read mask arrays from BIMG mask files (.dta, .dimg)
!  -------
!  Method : Open, read and close BIMG file
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
      CHARACTER(len=word80) :: cdrec1, cdrec2, cdrec3, cdrec4, kform
      INTEGER :: nx, ny, nz, nt, ndim, icod
      BIGREAL4 :: lonmin, latmin, dx, dy, spval
      BIGREAL4 :: date
      BIGREAL4, dimension(:), allocatable :: lev
      BIGREAL4, dimension(:,:), allocatable :: ptab
      CHARACTER(len=varlg) :: sxy_nam
      CHARACTER(len=bgword) :: sxyfmsk
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxypmsk,sxydmsk 
      INTEGER :: indsxy
      BIGREAL4 :: sxyvmsk
      LOGICAl :: sxymsea
      INTEGER :: sxy_pos,ind_msk,js,sxynbr
      INTEGER :: jjref,jkref,jtref
      INTEGER :: allocok,jpisize,jpjsize,jpksize
      INTEGER :: ji, jj, jk, jt, jdim
      INTEGER :: jpiend, jpjend, jpkend, jptend
!----------------------------------------------------------------------
      jpksize=size(mask,3)
! --- allocation lev
      allocate ( lev(1:jpksize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lev(:) = FREAL4(0.0)
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
      jpisize=sxy_jpi
      jpjsize=sxy_jpj
!
! --- allocation ptab
      allocate ( ptab(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptab(:,:) = FREAL4(0.0)
!
! Control print
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/readmsk/readmskbimg'
         WRITE(numout,*) '    ==> READING ',sxyfmsk(1:lenv(sxyfmsk))
      ENDIF
!
! Set size of array to read
      sxy_pos =mod(sxypmsk,10)
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
! -1.- Open BIMG file
! -------------------
!
      CALL openfile(numfil,sxyfmsk,kform=clunf)
!
! -2.- Read BIMG file header
! --------------------------
!
      REWIND(UNIT=numfil)
!
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec1
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec2
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec3
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) cdrec4
!
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) nx, ny, nz, nt, ndim, icod
!
      IF (   (jpiend.NE.nx) &
     &    .OR.(jpjend.NE.ny) &
     &    .OR.(jpkend.GT.nz) &
     &    .OR.(jptend.GT.nt) &
     &    .OR.(sxy_pos.GT.ndim) ) GOTO 102
!
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) lonmin,latmin,dx,dy,spval
      READ(UNIT=numfil,ERR=101,IOSTAT=iost) lev(1:nz)
!
! -3.- Set special value if optional mask is used
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
!         print *,sxyvmsk
      CASE(2)
! --- dta
         IF (largdtamsk) THEN
            indsxy  = dta_ord(jsxy)
            dtavmsk(indsxy)= spval
            sxyvmsk = FREAL4(dtavmsk(indsxy))
         ENDIF
!         print *,sxyvmsk
      CASE(3)
         GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -4.- Read BIMG mask arrays
! --------------------------
!
      js=1
      DO jt=1,jptend
         READ(UNIT=numfil,ERR=101,IOSTAT=iost) date
!
         DO jk=1,jpkend
            DO jdim=1,(sxy_pos-1)
               READ(UNIT=numfil,ERR=101,IOSTAT=iost)
            ENDDO
            READ(UNIT=numfil,ERR=101,IOSTAT=iost) ptab(:jpiend,:jpjend)
            DO jdim=(sxy_pos+1),ndim
               READ(UNIT=numfil,ERR=101,IOSTAT=iost)
            ENDDO
            DO jj=1,jpjend
               DO ji=1,jpiend
                  IF (sxymsea.EQV.(sxyvmsk.EQ.ptab(ji,jj))) THEN
                     mask(ji,jj,jk,jt)=mask(ji,jj,jk,jt)+(2**(ind_msk))
                     js=js+1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      sxynbr=js-1
!
! -5.- Close BIMG file 
! --------------------
!
      CLOSE (UNIT=numfil)
!
! -6.- Set vertical levels (Vx masks only)
! ----------------------------------------
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
! -7.- Control Print
! ------------------
!
      IF (nprint.GE.3) THEN 
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title1: ',cdrec1(1:lenv(cdrec1))
         WRITE(numout,kform) '- Title2: ',cdrec2(1:lenv(cdrec2))
         WRITE(numout,kform) '- Title3: ',cdrec3(1:lenv(cdrec3))
         WRITE(numout,kform) '- Title4: ',cdrec4(1:lenv(cdrec4))
         kform='(8x,a,5i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
         kform='(8x,a,e12.3)'
         WRITE(numout,kform) '- Special value: ',spval
         kform='(8x,a,4e12.3)'
         WRITE(numout,kform) '- Regular grid: ',lonmin,latmin,dx,dy
         IF ( (jpkend.NE.nz).OR.(jptend.NE.nt).OR.(sxy_pos.NE.1) ) THEN
            kform='(8x,a,4i5)'
            WRITE(numout,kform) '- Size of loaded array : ',jpiend, &
     &                                  jpjend,jpkend,jptend
            kform='(8x,a,i5,a,i5)'
            WRITE(numout,kform) '- Variable : ',sxy_pos,'/',ndim
         ENDIF
      ENDIF
!
! --- deallocation
      IF (allocated(lev)) deallocate(lev)
      IF (allocated(ptab)) deallocate(ptab)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'liobimg','readmskbimg')
 1001 CALL printerror2(0,1001,3,'liobimg','readmskbimg')
!
 101  WRITE (texterror,*) 'error reading BIMG file, iost=',iost
      CALL printerror2(0,101,3,'liobimg','readmskbimg',comment=texterror)
 102  WRITE (texterror,*) 'bad dimensions in BIMG file'
      CALL printerror2(0,102,3,'liobimg','readmskbimg',comment=texterror)
 103  WRITE (texterror,*) 'incoherence between mask files and ', &
     &                    'SESAM configuration'
      CALL printerror2(0,103,3,'liobimg','readmskbimg',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liobimg
