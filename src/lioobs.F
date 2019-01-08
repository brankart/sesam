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
C ---                  LIOOBS.F                                 ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12  (C.E. Testut)                       ---
C --- modification : 99-05  (C.E. Testut)                       ---
C --- modification : 01-06  (C.E. Testut)                       ---
C --- modification : 03-03  (J.M. Brankart)                     ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE evalhdrfileobs  : Read individual '.obs' file header
C --- SUBROUTINE readvalfileobs  : Read individual '.obs' file
C ---                              (observed values only)
C --- SUBROUTINE readcfgfileobs  : Read individual '.obs' file
C ---                              (observation configuration only)
C --- SUBROUTINE writehdrfileobs : Write empty individual '.obs' file
C --- SUBROUTINE writefileobs    : Write individual '.obs' file
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE lioobs
      use mod_main
      use utilfiles
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrfileobs,readvalfileobs,readcfgfileobs
      PUBLIC writehdrfileobs,writefileobs

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrfileobs (kfninobs,indobs,jpoloc,jpitploc)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read individual '.obs' file header
CCC  -------
CCC  Method : Read first record of direct access binary file
CCC  ------
CCC  Input :  kfninobs : filename
CCC  -----    indobs   : observed variable index
CCC
CCC  Output : jpoloc   : number of observations in file
CCC  ------   jpitploc : number of interpolation point in obs operator
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend, jpyend
      IMPLICIT NONE
C----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      INTEGER, intent(in) :: indobs
      INTEGER, intent(out) :: jpoloc,jpitploc
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=4) :: versobs
      CHARACTER(len=word80) :: cdrec1,kform
      INTEGER :: jpiobs,jpjobs,jpkobs,jptobs
      INTEGER :: jpxendobs,jpyendobs
      INTEGER :: allocok
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./evalhdrobs/evalhdrfileobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
C
C -1.- Open direct access observation file
C ----------------------------------------
C
      irecl  = ibloc*((1000)/ibloc+1)*jpbyt4
      CALL openfile(numfil,kfninobs,clold,clunf,cldir,irecl)
C
C -2.- Read observation file header
C ---------------------------------
C
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) 
     $     versobs,cdrec1,irecl,jpiobs,jpjobs,jpkobs,jptobs,
     $     jpxendobs,jpyendobs,jpoloc,jpitploc
      IF (versobs.NE.'@obs') GOTO 102
C
C -3.- Close observation file
C ---------------------------
C
      CLOSE (UNIT=numfil)
C
C -4.- Control print
C ------------------
C
      IF (nprint.GE.3) THEN
         kform='(8x,2a)'
         WRITE(numout,kform) '- Version: ',versobs
         WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
         kform='(8x,a,i8)'
         WRITE(numout,kform) '- Record length: ',irecl
         kform='(8x,a,3i5)'
         WRITE(numout,kform) '- Vy grid dimensions: ',jpiobs,jpjobs,jpkobs
         kform='(8x,a,2i5)'
         WRITE(numout,kform) '- Vx and Vy sizes: ',jpxendobs,jpyendobs
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpoloc
         WRITE(numout,kform) '- Number of interpolation points: ',jpitploc
      ENDIF
C
C -5.- Coherence test
C -------------------
C
      IF ((jpoloc.NE.0).AND.(jpitploc.NE.0)) THEN
         IF (jpiobs.NE.dta_jpi(indobs)) GOTO 103
         IF (jpjobs.NE.dta_jpj(indobs)) GOTO 103
         IF (jpkobs.NE.dta_jpk(indobs)) GOTO 103
         IF (jptobs.NE.dta_jpt(indobs)) GOTO 103
         IF (jpxendobs.NE.jpxend) GOTO 103
         IF (jpyendobs.NE.jpyend) GOTO 103
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lioobs','evalhdrfileobs')
 1001 CALL printerror2(0,1001,3,'lioobs','evalhdrfileobs')
C
 101  WRITE (texterror,*) 'Error reading observation file, iost=',iost
      CALL printerror2(0,101,3,'lioobs','evalhdrfileobs',comment=texterror)
 102  WRITE (texterror,*) 'Bad observation file'
      CALL printerror2(0,102,3,'lioobs','evalhdrfileobs',comment=texterror)
 103  WRITE (texterror,*) 'Configuration in .obs file header',
     $     ' incompatible with SESAM configuration'
      CALL printerror2(0,103,3,'lioobs','evalhdrfileobs',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readvalfileobs(kfninobs,kvecto,kjobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read individual '.obs' file (observed values only)
CCC  -------
CCC  Method : Read direct access binary file
CCC  ------
CCC  Input :  kfninobs    : filename
CCC  -----    kjobs       : observed variable index
CCC
CCC  Output : kvecto      : vector of observed values in file
CCC  ------ 
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpyend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      BIGREAL, dimension(:), intent(out) :: kvecto
      INTEGER, intent(in) :: kjobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jposize,jpitpsize
      INTEGER :: indobs,inddbs
      CHARACTER(len=word80) :: cdrec1, kform
      INTEGER :: jrec
      CHARACTER(len=4) :: versobs
      BIGREAL4, dimension(:), allocatable :: ptabo
      BIGREAL :: sxyo_moy,sxyo_ect
      INTEGER :: jpiobs,jpjobs,jpkobs,jptobs
      INTEGER :: jpxendobs,jpyendobs,jpolocobs,jpitplocobs
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readvalobs/readvalfileobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
C
C Set input array sizes and observation indices
      jposize = size(kvecto,1)
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
      jpitpsize=obs_itp(indobs,inddbs)
C
C -1.- Open direct access observation file
C ----------------------------------------
C
      irecl  = ibloc*(1000/ibloc+1)*jpbyt4
      CALL openfile (numfil,kfninobs,clold,clunf,cldir,irecl)
C
C -2.- Read observation file header
C ---------------------------------
C
      versobs='@obs'
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) 
     $     versobs,cdrec1,irecl,jpiobs,jpjobs,jpkobs,jptobs,
     $     jpxendobs,jpyendobs,jpolocobs,jpitplocobs
      IF (versobs.NE.'@obs') GOTO 104
C
C -3.- Close observation file
C ---------------------------
C
      CLOSE (UNIT=numfil)
C
C -4.- Control print
C ------------------
C
      IF (nprint.GE.4) THEN
         kform='(8x,2a)'
         WRITE(numout,kform) '- Version: ',versobs
         WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
         kform='(8x,a,i8)'
         WRITE(numout,kform) '- Record length: ',irecl
         kform='(8x,a,3i5)'
         WRITE(numout,kform) '- Vy grid dimensions: ',jpiobs,jpjobs,jpkobs
         kform='(8x,a,2i5)'
         WRITE(numout,kform) '- Vx and Vy sizes: ',jpxendobs,jpyendobs
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpolocobs
         WRITE(numout,kform) '- Number of interpolation points: ',jpitplocobs
      ENDIF
C
C -5.- Coherence test
C -------------------
C
      IF (jpiobs.NE.dta_jpi(indobs)) GOTO 105
      IF (jpjobs.NE.dta_jpj(indobs)) GOTO 105
      IF (jpkobs.NE.dta_jpk(indobs)) GOTO 105
      IF (jpkobs.NE.dta_jpk(indobs)) GOTO 105
      IF (jpxendobs.NE.jpxend) GOTO 105
      IF (jpyendobs.NE.jpyend) GOTO 105
      IF (jpolocobs.NE.jposize) GOTO 105
      IF (jpitplocobs.NE.jpitpsize) GOTO 105
C
C -6.- Re-open file with the right record length
C ----------------------------------------------
C
      CALL openfile (numfil,kfninobs,clold,clunf,cldir,irecl)
C
C -7.- Read observation file
C --------------------------
C
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
C Read observation values
      jrec=5
      READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
      IF (lmoyect) THEN
         sxyo_moy=FREAL4(obs_moy(indobs,inddbs))
         sxyo_ect=FREAL4(obs_ect(indobs,inddbs))
         kvecto(:) = (FREAL(ptabo(:))-sxyo_moy)/sxyo_ect
      ELSE
         kvecto(:) = FREAL(ptabo(:))
      ENDIF
      kvecto(:) = FREAL(ptabo(:))
C
C -8.- Close observation file
C ---------------------------
C
      CLOSE (UNIT=numfil)
C
C --- deallocation ptabo
      IF (allocated(ptabo)) deallocate (ptabo)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lioobs','readvalfileobs')
 1001 CALL printerror2(0,1001,3,'lioobs','readvalfileobs')
C
 101  WRITE (texterror,*) 'Error reading observation file, iost=',iost
      CALL printerror2(0,101,3,'lioobs','readvalfileobs',comment=texterror)
 104  WRITE (texterror,*) 'Bad observation file'
      CALL printerror2(0,104,3,'lioobs','readvalfileobs',comment=texterror)
 105  WRITE (texterror,*) 'Configuration in .obs file header',
     $     ' incompatible with SESAM configuration'
      CALL printerror2(0,105,3,'lioobs','readvalfileobs',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgfileobs(kfninobs,kjobs,kflagcfg,
     $    kvectorms,kgridijkobs,kposcoefobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read individual '.obs' file (observation configuration only)
CCC  -------
CCC  Method : Read direct access binary file
CCC  ------
CCC  Input :  kfninobs    : filename
CCC  -----    kjobs       : observed variable index
CCC
CCC  Output : kgridijkobs : observation location (x,y,z)
CCC  ------   kposcoefobs : observation operator (interpolation points
CCC                         and interpolation coefficients)
CCC           kvectorms   : associated error value (obsolete)
CCC           kflagcfg    : what element of the configuration to read
CCC                         (1=kvectorms, 2=kgridijkobs, 3=kposcoefobs)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpyend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      INTEGER, intent(in) :: kjobs,kflagcfg
      BIGREAL, dimension(:), optional, intent(out) :: kvectorms
      TYPE (type_gridijk), dimension(:), optional, intent(out) :: 
     $     kgridijkobs
      TYPE (type_poscoef), dimension(:,:), optional, intent(out) :: 
     $     kposcoefobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jposize,jpitpsize
      INTEGER :: indobs,inddbs
      CHARACTER(len=word80) :: cdrec1,kform
      INTEGER :: jrec,jitp
      CHARACTER(len=4) :: versobs
      BIGREAL4, dimension(:), allocatable :: ptabo
      INTEGER :: jpiobs,jpjobs,jpkobs,jptobs
      INTEGER :: jpxendobs,jpyendobs,jpolocobs,jpitplocobs
      INTEGER :: allocok
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readcfgobs/readcfgfileobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
C
C Set input array sizes and observation indices
C Check coherence of input array sizes
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
      SELECT CASE (kflagcfg)
      CASE(1)
         IF (.NOT.(present(kvectorms))) GOTO 1000
         jposize = size(kvectorms,1)
         jpitpsize = obs_itp(indobs,inddbs)
      CASE(2)
         IF (.NOT.(present(kgridijkobs))) GOTO 1000
         jposize = size(kgridijkobs,1)
         jpitpsize = obs_itp(indobs,inddbs)
      CASE(3)
         IF (.NOT.(present(kposcoefobs))) GOTO 1000
         jposize = size(kposcoefobs,1)
         jpitpsize = size(kposcoefobs,2)
         IF (jpitpsize.LT.obs_itp(indobs,inddbs)) GOTO 102
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C -1.- Open direct access observation file
C ----------------------------------------
C
      irecl  = ibloc*(1000/ibloc+1)*jpbyt4
      CALL openfile (numfil,kfninobs,clold,clunf,cldir,irecl)
C
C -2.- Read observation file header
C ---------------------------------
C
      versobs='@obs'
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost) 
     $     versobs,cdrec1,irecl,jpiobs,jpjobs,jpkobs,jptobs,
     $     jpxendobs,jpyendobs,jpolocobs,jpitplocobs
      IF (versobs.NE.'@obs') GOTO 104
C
C -3.- Close observation file
C ---------------------------
C
      CLOSE (UNIT=numfil)
C
C -4.- Control print
C ------------------
C
      IF (nprint.GE.4) THEN
         kform='(8x,2a)'
         WRITE(numout,kform) '- Version: ',versobs
         WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
         kform='(8x,a,i8)'
         WRITE(numout,kform) '- Record length: ',irecl
         kform='(8x,a,3i5)'
         WRITE(numout,kform) '- Vy grid dimensions: ',jpiobs,jpjobs,jpkobs
         kform='(8x,a,2i5)'
         WRITE(numout,kform) '- Vx and Vy sizes: ',jpxendobs,jpyendobs
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpolocobs
         WRITE(numout,kform) '- Number of interpolation points: ',jpitplocobs
      ENDIF
C
C -5.- Coherence test
C -------------------
C
      IF (jpiobs.NE.dta_jpi(indobs)) GOTO 105
      IF (jpjobs.NE.dta_jpj(indobs)) GOTO 105
      IF (jpkobs.NE.dta_jpk(indobs)) GOTO 105
      IF (jpkobs.NE.dta_jpk(indobs)) GOTO 105
      IF (jpxendobs.NE.jpxend) GOTO 105
      IF (jpyendobs.NE.jpyend) GOTO 105
      IF (jpolocobs.NE.jposize) GOTO 105
      IF (jpitplocobs.NE.jpitpsize) GOTO 105
C
C -6.- Re-open file with the right record length
C ----------------------------------------------
C
      CALL openfile (numfil,kfninobs,clold,clunf,cldir,irecl)
C
C -7.- Read observation file
C --------------------------
C
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
C Read observation longitudes
      jrec=2
      IF (kflagcfg.EQ.2) THEN
         READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
         kgridijkobs(:)%longi = FREAL(ptabo(:))
      ENDIF
C
C Read observation latitudes
      jrec=jrec+1
      IF (kflagcfg.EQ.2) THEN
         READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
         kgridijkobs(:)%latj = FREAL(ptabo(:))
      ENDIF
C
C Read observation depths
      jrec=jrec+1
      IF (kflagcfg.EQ.2) THEN
         READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
         kgridijkobs(:)%levk = FREAL(ptabo(:))
      ENDIF
C
C Do not read observation values
      jrec=jrec+1
C
C Read associated observation error (obsolete)
      jrec=jrec+1
      IF (kflagcfg.EQ.1) THEN
         READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
         kvectorms(:) = FREAL(ptabo(:))
      ENDIF
C
C Read observation operator
      IF (kflagcfg.EQ.3) THEN
         DO jitp=1,jpitpsize
C
            jrec=jrec+1
            READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
            kposcoefobs(:,jitp)%pos = FREAL(ptabo(:))
C
            jrec=jrec+1
            READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
            kposcoefobs(:,jitp)%coef = FREAL(ptabo(:))
C     
         ENDDO
      ENDIF
C
C -8.- Close observation file
C ---------------------------
C
      CLOSE (UNIT=numfil)
C
C --- deallocation
      IF (allocated(ptabo)) deallocate(ptabo)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lioobs','readcfgfileobs')
 1001 CALL printerror2(0,1001,3,'lioobs','readcfgfileobs')
C
 101  WRITE (texterror,*) 'Error reading observation file, iost=',iost
      CALL printerror2(0,101,3,'lioobs','readcfgfileobs',comment=texterror)
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,1,'lioobs','readcfgfileobs',comment=texterror)
 104  WRITE (texterror,*) 'Bad observation file'
      CALL printerror2(0,104,3,'lioobs','readcfgfileobs',comment=texterror)
 105  WRITE (texterror,*) 'Configuration in .obs file header',
     $     ' incompatible with SESAM configuration'
      CALL printerror2(0,105,3,'lioobs','readcfgfileobs',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrfileobs(kfnoutobs,kjobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write empty '.obs' file
CCC  -------
CCC  Method : Write first record of direct access binary file
CCC  ------
CCC  Input :  kfnoutobs : filename
CCC  -----    kjobs    : observed variable index
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpyend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutobs
      INTEGER, intent(in) :: kjobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jposize,jpitpsize
      INTEGER :: indobs,inddbs
      CHARACTER(len=word80) :: cdrec1,kform
      INTEGER :: jrec
      CHARACTER(len=4) :: versobs
      INTEGER :: jpiobs,jpjobs,jpkobs,jptobs
      INTEGER :: jpxendobs,jpyendobs,jpolocobs,jpitplocobs
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writeobs/writehdrfileobs'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutobs(1:lenv(kfnoutobs))
      ENDIF
C
C Set output array sizes and observation indices
      jposize = 0
      jpitpsize = 0
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
C
C -1.- Open direct access observation file
C ----------------------------------------
C
      irecl  = ibloc*((1000)/ibloc+1)*jpbyt4
      CALL openfile (numfil,kfnoutobs,clunk,clunf,cldir,irecl)
C
C -2.- Write observation file header
C ----------------------------------
C
      WRITE (cdrec1,'(2a)') 'Observation vector from database: ',
     $     obs_nam(indobs,inddbs)(1:lenv(obs_nam(indobs,inddbs)))
C
      versobs = '@obs'
      jpiobs = dta_jpi(indobs)
      jpjobs = dta_jpj(indobs)
      jpkobs = dta_jpk(indobs)
      jptobs = dta_jpt(indobs)
      jpxendobs = jpxend
      jpyendobs = jpyend
      jpolocobs = jposize
      jpitplocobs = jpitpsize
C
      WRITE(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost)  
     $     versobs,cdrec1,irecl,jpiobs,jpjobs,jpkobs,jptobs,
     $     jpxendobs,jpyendobs,jpolocobs,jpitplocobs
C
C -3.- Close observation file
C ---------------------------
C
      CLOSE (UNIT=numfil)
C
C -4.- Control print
C ------------------
C
      IF (nprint.GE.4) THEN
         kform='(8x,2a)'
         WRITE(numout,kform) '- Version: ',versobs
         WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
         kform='(8x,a,i8)'
         WRITE(numout,kform) '- Record length: ',irecl
         kform='(8x,a,3i5)'
         WRITE(numout,kform) '- Vy grid dimensions: ',jpiobs,jpjobs,jpkobs
         kform='(8x,a,2i5)'
         WRITE(numout,kform) '- Vx and Vy sizes: ',jpxendobs,jpyendobs
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpolocobs
         WRITE(numout,kform) '- Number of interpolation points: ',jpitplocobs
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lioobs','writehdrfileobs')
C
 101  WRITE (texterror,*) 'Error writing observation file, iost=',iost
      CALL printerror2(0,101,3,'lioobs','writehdrfileobs',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writefileobs(kfnoutobs,kvecto,kvectorms,kgridijkobs,
     $     kposcoefobs,kjobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write individual '.obs' file
CCC  -------
CCC  Method : Write direct access binary file
CCC  ------
CCC  Input : kfnoutobs   : filename
CCC  -----   kjobs       : observed variable index
CCC          kvecto      : vector of observed values in file
CCC          kvectorms   : associated error value (obsolete)
CCC          kgridijkobs : observation location (x,y,z)
CCC          kposcoefobs : observation operator (interpolation points
CCC                         and interpolation coefficients)
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpyend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutobs
      BIGREAL, dimension(:), intent(in) :: kvecto,kvectorms
      TYPE (type_gridijk), dimension(:), intent(in)  :: kgridijkobs
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
      INTEGER, intent(in) :: kjobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jposize,jpitpsize
      INTEGER :: indobs,inddbs
      CHARACTER(len=word80) :: cdrec1,kform
      INTEGER :: jrec,jitp,jo
      CHARACTER(len=4) :: versobs
      BIGREAL4, dimension(:), allocatable :: ptabo
      BIGREAL4 :: sxyo_moy,sxyo_ect
      INTEGER :: jpiobs,jpjobs,jpkobs,jptobs
      INTEGER :: jpxendobs,jpyendobs,jpolocobs,jpitplocobs
CC----------------------------------------------------------------------
C
C Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writeobs/writehdrfileobs'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutobs(1:lenv(kfnoutobs))
      ENDIF
C
C Set input array sizes and observation indices
      jposize = size(kvecto,1)
      jpitpsize = size(kposcoefobs,2)
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
C
C Check coherence of input array sizes
      IF (jposize.NE.size(kvectorms,1)) GOTO 102
      IF (jposize.NE.size(kgridijkobs,1)) GOTO 102
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 102
      IF (jpitpsize.NE.obs_itp(indobs,inddbs)) GOTO 102
C
C -1.- Open direct access observation file
C ----------------------------------------
C
#if defined _NEC
      irecl = max(jposize*jpbyt4,2048)
#else
      irecl = ibloc*((jposize)/ibloc+1)*jpbyt4
#endif
      CALL openfile (numfil,kfnoutobs,clunk,clunf,cldir,irecl)
C
C -2.- Write observation file header
C ----------------------------------
C
      WRITE (cdrec1,'(2a)') 'Observation vector from database: ',
     $     obs_nam(indobs,inddbs)(1:lenv(obs_nam(indobs,inddbs)))
C
      versobs = '@obs'
      jpiobs = dta_jpi(indobs)
      jpjobs = dta_jpj(indobs)
      jpkobs = dta_jpk(indobs)
      jptobs = dta_jpt(indobs)
      jpxendobs = jpxend
      jpyendobs = jpyend
      jpolocobs = jposize
      jpitplocobs = jpitpsize
C
      WRITE(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost)  
     $     versobs,cdrec1,irecl,jpiobs,jpjobs,jpkobs,jptobs,
     $     jpxendobs,jpyendobs,jpolocobs,jpitplocobs
C
C -3.- Write observation file
C ---------------------------
C
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
C Write observation longitudes
      jrec=2
      ptabo(:) = FREAL4(kgridijkobs(:)%longi)
      WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
C
C Write observation latitudes
      jrec=jrec+1
      ptabo(:) = FREAL4(kgridijkobs(:)%latj)
      WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
C
C Write observation depths
      jrec=jrec+1
      ptabo(:) = FREAL4(kgridijkobs(:)%levk)
      WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
C
C Write observation values
      jrec=jrec+1
      IF (lmoyect) THEN
         sxyo_moy=FREAL4(obs_moy(indobs,inddbs))
         sxyo_ect=FREAL4(obs_ect(indobs,inddbs))
         ptabo(:) = FREAL4(kvecto(:))*sxyo_ect+sxyo_moy
      ELSE
         ptabo(:) = FREAL4(kvecto(:))
      ENDIF
      WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
C
C Write associated observation error (obsolete)
      jrec=jrec+1
      IF (lmoyect) THEN
         sxyo_moy=FREAL4(obs_moy(indobs,inddbs))
         sxyo_ect=FREAL4(obs_ect(indobs,inddbs))
         ptabo(:) = FREAL4(kvectorms(:))*sxyo_ect+sxyo_moy
      ELSE
         ptabo(:) = FREAL4(kvectorms(:))
      ENDIF
      WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
C
C Write observation operator
      DO jitp=1,obs_itp(indobs,inddbs)
C
         jrec=jrec+1
         ptabo(:) = FREAL4(kposcoefobs(:,jitp)%pos)
         WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
C
         jrec=jrec+1
         ptabo(:) = FREAL4(kposcoefobs(:,jitp)%coef)
         WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
C
      ENDDO
C
C -4.- Close observation file
C ---------------------------
C
      CLOSE (UNIT=numfil)
C
C -5.- Control print
C ------------------
C
      IF (nprint.GE.4) THEN
         kform='(8x,2a)'
         WRITE(numout,kform) '- Version: ',versobs
         WRITE(numout,kform) '- Title: ',cdrec1(1:lenv(cdrec1))
         kform='(8x,a,i8)'
         WRITE(numout,kform) '- Record length: ',irecl
         kform='(8x,a,3i5)'
         WRITE(numout,kform) '- Vy grid dimensions: ',jpiobs,jpjobs,jpkobs
         kform='(8x,a,2i5)'
         WRITE(numout,kform) '- Vx and Vy sizes: ',jpxendobs,jpyendobs
         kform='(8x,a,i5)'
         WRITE(numout,kform) '- Number of observations in file: ',jpolocobs
         WRITE(numout,kform) '- Number of interpolation points: ',jpitplocobs
      ENDIF

      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lioobs','writefileobs')
 1001 CALL printerror2(0,1001,3,'lioobs','writefileobs')
C
 101  WRITE (texterror,*) 'Error writing observation file, iost=',iost
      CALL printerror2(0,101,3,'lioobs','writefileobs',comment=texterror)
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,1,'lioobs','writefileobs',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE lioobs
