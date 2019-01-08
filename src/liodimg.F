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
C ---                  LIODIMG.F                                ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12  (C.E. Testut)                       ---
C --- modification : 99-05  (C.E. Testut)                       ---
C --- modification : 01-06  (C.E. Testut)                       ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE readdimg       : Read variable field from DIMG file
C --- SUBROUTINE writendimg     : Write variable field in DIMG file
C --- SUBROUTINE evalhdrmskdimg : Get variable dimensions from mask files
C --- SUBROUTINE readmskdimg    : Read mask arrays from DIMG mask files
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE liodimg
      use mod_main
      use utilvct
      use utilfiles
      IMPLICIT NONE
      PRIVATE

      PUBLIC readdimg,writendimg,evalhdrmskdimg,readmskdimg

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readdimg(kfname,kvectsout,kjsxy,ksomsxynbr,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read variable field (of Vx or Vy object) from DIMG file
CCC  -------
CCC  Method : Open, read and close DIMG file (.dimg or .dta)
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
      CHARACTER(len=word80) :: cdrec1,kform
      CHARACTER(len=4) :: verdimg
      INTEGER :: irec
      INTEGER :: nx, ny, nz, nt, ndim
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
CC----------------------------------------------------------------------
C Size of output vector
      jpssize=size(kvectsout,1)
C
C Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE(1)
C --- var
         indsxy  = var_ord(kjsxy)
         varfile=.TRUE.
      CASE(2)
C --- dta
         indsxy  = dta_ord(kjsxy)
         jextvar=indext(kfname,extvartab,nbextvar)
         jextdta=indext(kfname,extdtatab,nbextdta)
         IF (jextdta.EQ.2) THEN
            varfile=.FALSE.
         ELSEIF (jextvar.EQ.2) THEN
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
         sxy_dim = var_dim(indsxy)
         IF (sxy_dim.GE.1) sxy_jpi = var_jpi(indsxy)
         IF (sxy_dim.GE.2) sxy_jpj = var_jpj(indsxy)
         IF (sxy_dim.GE.3) sxy_jpk = var_jpk(indsxy)
         IF (sxy_dim.GE.4) sxy_jpt = var_jpt(indsxy)
         ksxypos = mod(varipos(indsxy),10)
      ELSE
C --- read Vy variable field from dta file
         sxy_dim = dta_dim(indsxy)
         IF (sxy_dim.GE.1) sxy_jpi = dta_jpi(indsxy)
         IF (sxy_dim.GE.2) sxy_jpj = dta_jpj(indsxy)
         IF (sxy_dim.GE.3) sxy_jpk = dta_jpk(indsxy)
         IF (sxy_dim.GE.4) sxy_jpt = dta_jpt(indsxy)
         ksxypos = mod(dtaipos(indsxy),10)
      ENDIF
C
C Control print
      IF (nprint.GE.2) THEN 
         WRITE(numout,*) '*** ROUTINE : ../readfil/readdimg'
         WRITE(numout,*) '    ==> READING file ',kfname(1:lenv(kfname))
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
C --- allocation lev
      allocate ( lev(1:sxy_jpk), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lev(:) = FREAL4(0.0)
C
C --- allocation ptabij
      allocate ( ptabij(1:sxy_jpi,1:sxy_jpj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL4(0.0)
C
C -1.- Open DIMG file
C -------------------
C
C Open file with incorrect record lentgh
      irecl  = ibloc*((sxy_jpi*sxy_jpj)/ibloc+1)*jpbyt4
      CALL openfile (numfil,kfname,clunk,clunf,cldir,irecl)
C Read correct record lentgh at the beginning of the file
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) verdimg,cdrec1,irecl
      IF (verdimg.NE.'@!01') GOTO 103
      CLOSE (UNIT=numfil)
C Open file with correct record lentgh
      CALL openfile (numfil,kfname,clunk,clunf,cldir,irecl)
C
C -2.- Read DIMG file header (first record)
C -----------------------------------------
C
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) verdimg,cdrec1,irecl,
     $                                       nx, ny, nz,nt,ndim,
     $                                       lonmin,latmin,dx,dy,
     $                                       spval,
     $                                       lev(1:sxy_jpk),
     $                                       (date,jt=1,nt)
C
      IF (   (sxy_jpi.NE.nx)
     $    .OR.(sxy_jpj.NE.ny)
     $    .OR.(sxy_jpk.GT.nz)
     $    .OR.(sxy_jpt.GT.nt)
     $    .OR.(ksxypos.GT.ndim) ) GOTO 102
C
C -3.- Read 2D arrays in DIMG file
C --------------------------------
C
      js=1
      ksomsxynbr=0
      jdim=ksxypos
      DO jt=1,jptend
         DO jk=1,jpkend
            sompartsxynbr=0
            IF (js.LE.jpssize) THEN
               irec=2+(jt-1)*nz*ndim+(jk-1)*ndim+(jdim-1)
               READ(UNIT=numfil,REC=irec,ERR=101,IOSTAT=iost) ptabij(:,:)
               CALL mk4vct(kvectsout(js:),ptabij(:,:),
     $              jk,jt,kjsxy,sompartsxynbr,kflagxyo)
            ENDIF
            js  = js + sompartsxynbr
            ksomsxynbr= ksomsxynbr + sompartsxynbr
         ENDDO
      ENDDO
C
C -4.- Close DIMG file
C --------------------
C
      CLOSE (UNIT=numfil)
C
C -5.- Control Print
C ------------------
C
      IF (nprint.GE.3) THEN 
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
         kform='(8x,a,5i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
         kform='(8x,a,e12.3)'
         WRITE(numout,kform) '- Special value: ',spval
         kform='(8x,a,i8)'
         WRITE(numout,kform) '- Record length: ',irecl
         kform='(8x,a,4e12.3)'
         WRITE(numout,kform) '- Regular grid: ',lonmin,latmin,dx,dy
         IF ( (jpkend.NE.nz).OR.(jptend.NE.nt).OR.(ksxypos.NE.1) ) THEN
            kform='(8x,a,4i5)'
            WRITE(numout,kform) '- Size of loaded array : ',jpiend,
     $                                  jpjend,jpkend,jptend
            kform='(8x,a,i5,a,i5)'
            WRITE(numout,kform) '- Variable : ',ksxypos,'/',ndim
         ENDIF
      ENDIF
C
C --- deallocation
      IF (allocated( lev)) deallocate(lev)
      IF (allocated( ptabij)) deallocate(ptabij)
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'liodimg','readdimg')
 1001 CALL printerror2(0,1001,3,'liodimg','readdimg')
C
 101  WRITE (texterror,*) 'error reading DIMG file, iost=',iost
      CALL printerror2(0,101,3,'liodimg','readdimg',comment=texterror)
 102  WRITE (texterror,*) 'bad dimensions in DIMG file: ',kfname
      CALL printerror2(0,102,3,'liodimg','readdimg',comment=texterror)
 103  WRITE (texterror,*) 'bad DIMG file version'
      CALL printerror2(0,103,3,'liodimg','readdimg',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writendimg(kfnoutndimg,lgcdrec1,idast,kvectsin,
     $                ktxy_indvct,ktxy_indmsk,txyend,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write variable field (of Vx or Vy object) in DIMG file
CCC  -------
CCC  Method : Open, write and close DIMG file (.cdf or .cdta)
CCC  ------
CCC  Input : kfnoutndimg : Filename
CCC  -----   cdrec1      : File title
CCC          idast       : date (obsolete)
CCC          kvectsin    : 1D SESAM vector object (Vx or Vy)
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
      use mod_spacexyo , only : spvaldta,spvalvar
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutndimg
      CHARACTER(len=*), intent(in) :: lgcdrec1
      BIGREAL, intent(in) :: idast
      BIGREAL, dimension(:), intent(in) :: kvectsin
      INTEGER, dimension(:), intent(in) :: ktxy_indvct,ktxy_indmsk
      INTEGER, intent(in) :: txyend, kflagxyo
CC----------------------------------------------------------------------
      CHARACTER(len=4) :: verdimg
      CHARACTER(len=word80) :: cdrec1, kform
      INTEGER :: irec
      INTEGER :: nx, ny, nz, nt, ndim
      BIGREAL4 :: lonmin, latmin, dx, dy, spval
      BIGREAL4 :: date
      INTEGER :: allocok
      BIGREAL4, allocatable, dimension(:) :: lev
      BIGREAL4, allocatable, dimension(:,:) :: ptabij
      INTEGER :: sxyend
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_dim,
     $     sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt,sxy_nbr
      INTEGER :: jsxy,indsxy,sompartxynbr
      INTEGER :: jtxy
      INTEGER, dimension(1:nbvar) :: txy_indtab
      INTEGER :: ji, jj, jk, jt
      INTEGER :: jpiend, jpjend, jpkend, jptend
CC----------------------------------------------------------------------
C
C Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE(1)
C --- var
         sxyend  = varend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=var_ord(jsxy)
            indsxy= sxy_ord(jsxy)
C
            sxy_dim(indsxy)=var_dim(indsxy)
            sxy_jpi(indsxy)=var_jpi(indsxy)
            sxy_jpj(indsxy)=var_jpj(indsxy)
            sxy_jpk(indsxy)=var_jpk(indsxy)
            sxy_jpt(indsxy)=var_jpt(indsxy)
            sxy_nbr(indsxy)=var_nbr(indsxy)
         ENDDO
      CASE(2)
C --- dta
         sxyend  = dtaend
         DO jsxy = 1,sxyend
            sxy_ord(jsxy)=dta_ord(jsxy)
            indsxy= sxy_ord(jsxy)
C
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
C
C Control print
      IF (nprint.GE.2) THEN 
         WRITE(numout,*) '*** ROUTINE : ../writenfil/writendimg'
         WRITE(numout,*) '    ==> WRITING file: ',
     $                           kfnoutndimg(1:lenv(kfnoutndimg))
      ENDIF
C
C Compute dimensions of DIMG file
      jtxy=1
      jsxy=ktxy_indmsk(jtxy)
      indsxy=sxy_ord(jsxy)
      SELECT CASE (sxy_dim(indsxy))
      CASE (1)
C ==>  1D
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
C ==>  2D
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
C ==>  3D
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
C ==>  4D
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
C
C --- allocation lev
      allocate ( lev(1:jpkend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lev(:) = FREAL4(0.0)
C
C --- allocation ptabij
      allocate ( ptabij(1:jpiend,1:jpjend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptabij(:,:) = FREAL4(0.0)
C
C -1.- Set parameters to write in DIMG file
C -----------------------------------------
C
#if defined _NEC
      irecl   = (jpiend*jpjend)*jpbyt4
#else
      irecl   = ibloc*((jpiend*jpjend)/ibloc+1)*jpbyt4
#endif
C
      verdimg = '@!01'
      nx      = jpiend
      ny      = jpjend
      nz      = jpkend
      nt      = jptend
      ndim    = txyend
      lonmin  = FREAL4(1.)
      latmin  = FREAL4(1.)
      dx      = FREAL4(1.)
      dy      = FREAL4(1.)
C
      SELECT CASE(kflagxyo)
      CASE(1)
C --- var
         spval  = FREAL4(spvalvar)
      CASE(2)
C --- dta
         spval  = FREAL4(spvaldta)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      lev(1:jpkend) =  FREAL4(var_lev(1:jpkend,indsxy))
      date    = FREAL4(idast)
C
      txy_indtab(1:txyend)=ktxy_indvct(1:txyend)
C
C -1.- Open DIMG file
C -------------------
C
      CALL openfile (numfil,kfnoutndimg,clunk,clunf,cldir,irecl)
C
C -2.- Write DIMG file header (first record)
C ------------------------------------------
C
      WRITE(cdrec1,'(a)') lgcdrec1(1:80)
      IF (jptend.EQ.1) THEN
         WRITE(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) 
     $        verdimg,cdrec1,irecl,
     $        nx, ny, nz,nt,ndim,
     $        lonmin,latmin,dx,dy,
     $        spval,
     $        lev(1:jpkend),
     $        (date,jt=1,jptend)
      ELSE
         WRITE(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) 
     $        verdimg,cdrec1,irecl,
     $        nx, ny, nz,nt,ndim,
     $        lonmin,latmin,dx,dy,
     $        spval,
     $        lev(1:jpkend),
     $        (jt,jt=1,jptend)
      ENDIF
C
C -4.- Write 2D arrays in DIMG file
C ---------------------------------
C
      DO jt=1,jptend
      DO jk=1,jpkend
         DO jtxy=1,txyend
            jsxy=ktxy_indmsk(jtxy)
            sompartxynbr=0
            IF (txy_indtab(jtxy).LE.size(kvectsin,1)) THEN
               CALL unmk4vct(kvectsin(txy_indtab(jtxy):),
     $              ptabij(:,:),jk,jt,jsxy,
     $              sompartxynbr,spval,kflagxyo) 
            ELSE
               ptabij(:,:) = spval
            ENDIF
            irec=2+(jt-1)*nz*ndim+(jk-1)*ndim+(jtxy-1)
            WRITE(UNIT=numfil,REC=irec,ERR=101,IOSTAT=iost) ptabij(:,:)
            txy_indtab(jtxy)=txy_indtab(jtxy)+sompartxynbr
         ENDDO
      ENDDO
      ENDDO
C
C coherence test
      DO jtxy=1,txyend
         jsxy=ktxy_indmsk(jtxy)
         indsxy=sxy_ord(jsxy)
         IF (sxy_nbr(indsxy).NE.(txy_indtab(jtxy)-ktxy_indvct(jtxy)))
     $                   GOTO 102
      ENDDO
C
C -5.- Close DIMG file
C --------------------
C
      CLOSE (UNIT=numfil)
C
C -6.- Control Print
C ------------------
C
      IF (nprint.GE.3) THEN 
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
         kform='(8x,a,5i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
         kform='(8x,a,e12.3)'
         WRITE(numout,kform) '- Special value: ',spval
         kform='(8x,a,i8)'
         WRITE(numout,kform) '- Record length: ',irecl
         kform='(8x,a,4e12.3)'
         WRITE(numout,kform) '- Regular grid: ',lonmin,latmin,dx,dy
      ENDIF
C
C --- deallocation
      IF (allocated(lev)) deallocate(lev)
      IF (allocated(ptabij)) deallocate (ptabij)
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'liodimg','writendimg')
 1001 CALL printerror2(0,1001,3,'liodimg','writendimg')
C
 101  WRITE (texterror,*) 'error writing DIMG file, iost=',iost
      CALL printerror2(0,101,3,'liodimg','writendimg',comment=texterror)
 102  WRITE (texterror,*) 'inconsistent data dimensions'
      CALL printerror2(0,102,3,'liodimg','writendimg',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrmskdimg(kfname,kjpi,kjpj,kjpk,kjpt)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Get variable dimensions from mask files
CCC  -------
CCC  Method : Open, read and close DIMG file header
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfname
      INTEGER, intent(out) :: kjpi,kjpj,kjpk,kjpt
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80) :: cdrec1, kform
      CHARACTER(len=4) :: verdimg
      INTEGER :: nx, ny, nz, nt, ndim
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN 
         WRITE(numout,*) '*** ROUTINE : sesam/readmsk/evalhdrmsk/evalhdrmskdimg'
         WRITE(numout,*) '    ==> READING ',kfname(1:lenv(kfname))
      ENDIF
C
C -1.- Open DIMG file
C -------------------
C
C Open file with incorrect record lentgh
      irecl   = 1000
      CALL openfile (numfil,kfname,clunk,clunf,cldir,irecl)
C Read correct record lentgh at the beginning of the file
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) verdimg,cdrec1,irecl
      IF (verdimg.NE.'@!01') GOTO 103
      CLOSE (UNIT=numfil)
C Open file with correct record lentgh
      CALL openfile (numfil,kfname,clunk,clunf,cldir,irecl)
C
C -2.- Read dimg file header
C --------------------------
C
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) verdimg,cdrec1,irecl,
     $                                       nx, ny, nz,nt,ndim
C
C -3.- Close DIMG file
C --------------------
C
      CLOSE (UNIT=numfil)
C
C -4.- Control Print
C ------------------
C
      IF (nprint.GE.3) THEN 
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
         kform='(8x,a,5i5)'
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
         kform='(8x,a,i8)'
         WRITE(numout,kform) '- Record length: ',irecl
      ENDIF
C
C -5.- Set output variables
C -------------------------
C
      kjpi=nx
      kjpj=ny
      kjpk=nz
      kjpt=nt
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'liodimg','evalhdrmskdimg')
 1001 CALL printerror2(0,1001,3,'liodimg','evalhdrmskdimg')
C
 101  WRITE (texterror,*) 'error reading DIMG file, iost=',iost
      CALL printerror2(0,101,3,'liodimg','evalhdrmskdimg',
     $     comment=texterror)
 103  WRITE (texterror,*) 'bad DIMG file version'
      CALL printerror2(0,103,3,'liodimg','evalhdrmsdimg',
     $     comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmskdimg(jsxy,kflagxyo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read mask arrays from DIMG mask files (.dta, .dimg)
CCC  -------
CCC  Method : Open, read and close DIMG file
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
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      INTEGER, intent(in) :: jsxy,kflagxyo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=word80) :: cdrec1, kform
      CHARACTER(len=4) :: verdimg
      INTEGER :: irec
      INTEGER :: nx, ny, nz, nt, ndim
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
CC----------------------------------------------------------------------
      jpksize=size(mask,3)
C --- allocation lev
      allocate ( lev(1:jpksize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lev(:) = FREAL4(0.0)
CC----------------------------------------------------------------------
C
C Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE(1)
C --- var
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
C --- dta
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
C
      jpisize=sxy_jpi
      jpjsize=sxy_jpj
C
C --- allocation ptab
      allocate ( ptab(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptab(:,:) = FREAL4(0.0)
C
C Control print
      IF (nprint.GE.1) THEN 
         WRITE(numout,*) '*** ROUTINE : sesam/readmsk/readmskdimg'
         WRITE(numout,*) '    ==> READING ',sxyfmsk(1:lenv(sxyfmsk))
      ENDIF
C
C Set size of array to read
      sxy_pos =mod(sxypmsk,10)
      SELECT CASE (sxydmsk)
      CASE (1)
C ==>  1D
         jpiend = sxy_jpi
         jpjend = 1
         jpkend = 1
         jptend = 1
      CASE (2)
C ==>  2D
         jpiend =  sxy_jpi
         jpjend =  sxy_jpj
         jpkend = 1
         jptend = 1
      CASE (3)
C ==>  3D
         jpiend =  sxy_jpi
         jpjend =  sxy_jpj
         jpkend =  sxy_jpk     
         jptend = 1
      CASE (4)
C ==>  4D
         jpiend =  sxy_jpi
         jpjend =  sxy_jpj
         jpkend =  sxy_jpk    
         jptend =  sxy_jpt
      CASE DEFAULT
C ==> ERROR
         GOTO 102
      END SELECT
C
C -1.- Open DIMG file
C -------------------
C Open file with incorrect record lentgh
      irecl   = 1000
      CALL openfile (numfil,sxyfmsk,clunk,clunf,cldir,irecl)
C Read correct record lentgh at the beginning of the file
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) verdimg,cdrec1,irecl
      IF (verdimg.NE.'@!01') GOTO 104
      CLOSE (UNIT=numfil)
C Open file with correct record lentgh
      CALL openfile (numfil,sxyfmsk,clunk,clunf,cldir,irecl)
C
C -2.- Read DIMG file header
C --------------------------
C
      REWIND(UNIT=numfil)
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) verdimg,cdrec1,irecl,
     $                                       nx, ny, nz,nt,ndim,
     $                                       lonmin,latmin,dx,dy,
     $                                       spval,
     $                                       lev(1:nz),
     $                                       (date,jt=1,nt)
C
      IF (   (jpiend.NE.nx)
     $    .OR.(jpjend.NE.ny)
     $    .OR.(jpkend.GT.nz)
     $    .OR.(jptend.GT.nt)
     $    .OR.(sxy_pos.GT.ndim) ) GOTO 101
C
C -3.- Set special value if optional mask is used
C -----------------------------------------------
C
      SELECT CASE(kflagxyo)
      CASE(1)
C --- var
         IF (largvarmsk) THEN
            indsxy  = var_ord(jsxy)
            varvmsk(indsxy) = spval
            sxyvmsk = FREAL4(varvmsk(indsxy))
         ENDIF
      CASE(2)
C --- dta
         IF (largdtamsk) THEN
            indsxy  = dta_ord(jsxy)
            dtavmsk(indsxy) = spval
            sxyvmsk = FREAL4(dtavmsk(indsxy))
         ENDIF
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -4.- Read DIMG mask arrays
C --------------------------
C
      js=1
      DO jt=1,jptend
         DO jk=1,jpkend
            jdim=sxy_pos
            irec=2+(jt-1)*nz*ndim+(jk-1)*ndim+(jdim-1)
            READ(UNIT=numfil,REC=irec,ERR=101,IOSTAT=iost) 
     $           ptab(:jpiend,:jpjend)
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
C
C -5.- Close DIMG file
C --------------------
C
      CLOSE (UNIT=numfil)
C
C -6.- Set vertical levels (Vx masks only)
C ----------------------------------------
C
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
C
C -7.- Control Print
C ------------------
C
      IF (nprint.GE.3) THEN 
         kform='(8x,2a)'
         WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
         kform='(8x,a,5i5)'             
         WRITE(numout,kform) '- Dimensions: ',nx,ny,nz,nt,ndim
         kform='(8x,a,e12.3)'
         WRITE(numout,kform) '- Special value: ',spval
         kform='(8x,a,i8)'
         WRITE(numout,kform) '- Record length: ',irecl
         kform='(8x,a,4e12.3)'
         WRITE(numout,kform) '- Regular grid: ',lonmin,latmin,dx,dy
         IF ( (jpkend.NE.nz).OR.(jptend.NE.nt).OR.(sxy_pos.NE.1) ) THEN
            kform='(8x,a,4i5)'
            WRITE(numout,kform) '- Size of loaded array : ',jpiend,
     $                                  jpjend,jpkend,jptend
            kform='(8x,a,i5,a,i5)'
            WRITE(numout,kform) '- Variable : ',sxy_pos,'/',ndim
         ENDIF
      ENDIF
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'liodimg','readmskdimg')
 1001 CALL printerror2(0,1001,3,'liodimg','readmskdimg')
C
 101  WRITE (texterror,*) 'error reading DIMG file, iost=',iost
      CALL printerror2(0,101,3,'liodimg','readmskdimg',comment=texterror)
 102  WRITE (texterror,*) 'bad dimensions in DIMG file'
      CALL printerror2(0,102,3,'liodimg','readmskbdmg',comment=texterror)
 103  WRITE (texterror,*) 'incoherence between mask files and ',
     $                    'SESAM configuration'
      CALL printerror2(0,103,3,'liodimg','readmskdimg',comment=texterror)
 104  WRITE (texterror,*) 'bad DIMG file version'
      CALL printerror2(0,104,3,'liodimg','readmskdimg',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE liodimg
