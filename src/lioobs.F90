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
! ---                  LIOOBS.F90                                 ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12  (C.E. Testut)                       ---
! --- modification : 99-05  (C.E. Testut)                       ---
! --- modification : 01-06  (C.E. Testut)                       ---
! --- modification : 03-03  (J.M. Brankart)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE evalhdrfileobs  : Read individual '.obs' file header
! --- SUBROUTINE readvalfileobs  : Read individual '.obs' file
! ---                              (observed values only)
! --- SUBROUTINE readcfgfileobs  : Read individual '.obs' file
! ---                              (observation configuration only)
! --- SUBROUTINE writehdrfileobs : Write empty individual '.obs' file
! --- SUBROUTINE writefileobs    : Write individual '.obs' file
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE lioobs
      use mod_main
      use utilfiles
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrfileobs,readvalfileobs,readcfgfileobs
      PUBLIC writehdrfileobs,writefileobs

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrfileobs (kfninobs,indobs,jpoloc,jpitploc)
!---------------------------------------------------------------------
!
!  Purpose : Read individual '.obs' file header
!  -------
!  Method : Read first record of direct access binary file
!  ------
!  Input :  kfninobs : filename
!  -----    indobs   : observed variable index
!
!  Output : jpoloc   : number of observations in file
!  ------   jpitploc : number of interpolation point in obs operator
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend, jpyend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      INTEGER, intent(in) :: indobs
      INTEGER, intent(out) :: jpoloc,jpitploc
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=4) :: versobs
      CHARACTER(len=word80) :: cdrec1,kform
      INTEGER :: jpiobs,jpjobs,jpkobs,jptobs
      INTEGER :: jpxendobs,jpyendobs
      INTEGER :: allocok
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./evalhdrobs/evalhdrfileobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
!
! -1.- Open direct access observation file
! ----------------------------------------
!
      irecl  = ibloc*((1000)/ibloc+1)*jpbyt4
      CALL openfile(numfil,kfninobs,clold,clunf,cldir,irecl)
!
! -2.- Read observation file header
! ---------------------------------
!
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost)  &
     &     versobs,cdrec1,irecl,jpiobs,jpjobs,jpkobs,jptobs, &
     &     jpxendobs,jpyendobs,jpoloc,jpitploc
      IF (versobs.NE.'@obs') GOTO 102
!
! -3.- Close observation file
! ---------------------------
!
      CLOSE (UNIT=numfil)
!
! -4.- Control print
! ------------------
!
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
!
! -5.- Coherence test
! -------------------
!
      IF ((jpoloc.NE.0).AND.(jpitploc.NE.0)) THEN
         IF (jpiobs.NE.dta_jpi(indobs)) GOTO 103
         IF (jpjobs.NE.dta_jpj(indobs)) GOTO 103
         IF (jpkobs.NE.dta_jpk(indobs)) GOTO 103
         IF (jptobs.NE.dta_jpt(indobs)) GOTO 103
         IF (jpxendobs.NE.jpxend) GOTO 103
         IF (jpyendobs.NE.jpyend) GOTO 103
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioobs','evalhdrfileobs')
 1001 CALL printerror2(0,1001,3,'lioobs','evalhdrfileobs')
!
 101  WRITE (texterror,*) 'Error reading observation file, iost=',iost
      CALL printerror2(0,101,3,'lioobs','evalhdrfileobs',comment=texterror)
 102  WRITE (texterror,*) 'Bad observation file'
      CALL printerror2(0,102,3,'lioobs','evalhdrfileobs',comment=texterror)
 103  WRITE (texterror,*) 'Configuration in .obs file header', &
     &     ' incompatible with SESAM configuration'
      CALL printerror2(0,103,3,'lioobs','evalhdrfileobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readvalfileobs(kfninobs,kvecto,kjobs)
!---------------------------------------------------------------------
!
!  Purpose : Read individual '.obs' file (observed values only)
!  -------
!  Method : Read direct access binary file
!  ------
!  Input :  kfninobs    : filename
!  -----    kjobs       : observed variable index
!
!  Output : kvecto      : vector of observed values in file
!  ------ 
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpyend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      BIGREAL, dimension(:), intent(out) :: kvecto
      INTEGER, intent(in) :: kjobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jposize,jpitpsize
      INTEGER :: indobs,inddbs
      CHARACTER(len=word80) :: cdrec1, kform
      INTEGER :: jrec
      CHARACTER(len=4) :: versobs
      BIGREAL4, dimension(:), allocatable :: ptabo
      BIGREAL :: sxyo_moy,sxyo_ect
      INTEGER :: jpiobs,jpjobs,jpkobs,jptobs
      INTEGER :: jpxendobs,jpyendobs,jpolocobs,jpitplocobs
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readvalobs/readvalfileobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
!
! Set input array sizes and observation indices
      jposize = size(kvecto,1)
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
      jpitpsize=obs_itp(indobs,inddbs)
!
! -1.- Open direct access observation file
! ----------------------------------------
!
      irecl  = ibloc*(1000/ibloc+1)*jpbyt4
      CALL openfile (numfil,kfninobs,clold,clunf,cldir,irecl)
!
! -2.- Read observation file header
! ---------------------------------
!
      versobs='@obs'
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost)  &
     &     versobs,cdrec1,irecl,jpiobs,jpjobs,jpkobs,jptobs, &
     &     jpxendobs,jpyendobs,jpolocobs,jpitplocobs
      IF (versobs.NE.'@obs') GOTO 104
!
! -3.- Close observation file
! ---------------------------
!
      CLOSE (UNIT=numfil)
!
! -4.- Control print
! ------------------
!
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
!
! -5.- Coherence test
! -------------------
!
      IF (jpiobs.NE.dta_jpi(indobs)) GOTO 105
      IF (jpjobs.NE.dta_jpj(indobs)) GOTO 105
      IF (jpkobs.NE.dta_jpk(indobs)) GOTO 105
      IF (jpkobs.NE.dta_jpk(indobs)) GOTO 105
      IF (jpxendobs.NE.jpxend) GOTO 105
      IF (jpyendobs.NE.jpyend) GOTO 105
      IF (jpolocobs.NE.jposize) GOTO 105
      IF (jpitplocobs.NE.jpitpsize) GOTO 105
!
! -6.- Re-open file with the right record length
! ----------------------------------------------
!
      CALL openfile (numfil,kfninobs,clold,clunf,cldir,irecl)
!
! -7.- Read observation file
! --------------------------
!
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
! Read observation values
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
!
! -8.- Close observation file
! ---------------------------
!
      CLOSE (UNIT=numfil)
!
! --- deallocation ptabo
      IF (allocated(ptabo)) deallocate (ptabo)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioobs','readvalfileobs')
 1001 CALL printerror2(0,1001,3,'lioobs','readvalfileobs')
!
 101  WRITE (texterror,*) 'Error reading observation file, iost=',iost
      CALL printerror2(0,101,3,'lioobs','readvalfileobs',comment=texterror)
 104  WRITE (texterror,*) 'Bad observation file'
      CALL printerror2(0,104,3,'lioobs','readvalfileobs',comment=texterror)
 105  WRITE (texterror,*) 'Configuration in .obs file header', &
     &     ' incompatible with SESAM configuration'
      CALL printerror2(0,105,3,'lioobs','readvalfileobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcfgfileobs(kfninobs,kjobs,kflagcfg, &
     &    kvectorms,kgridijkobs,kposcoefobs)
!---------------------------------------------------------------------
!
!  Purpose : Read individual '.obs' file (observation configuration only)
!  -------
!  Method : Read direct access binary file
!  ------
!  Input :  kfninobs    : filename
!  -----    kjobs       : observed variable index
!
!  Output : kgridijkobs : observation location (x,y,z)
!  ------   kposcoefobs : observation operator (interpolation points
!                         and interpolation coefficients)
!           kvectorms   : associated error value (obsolete)
!           kflagcfg    : what element of the configuration to read
!                         (1=kvectorms, 2=kgridijkobs, 3=kposcoefobs)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpyend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobs
      INTEGER, intent(in) :: kjobs,kflagcfg
      BIGREAL, dimension(:), optional, intent(out) :: kvectorms
      TYPE (type_gridijk), dimension(:), optional, intent(out) ::  &
     &     kgridijkobs
      TYPE (type_poscoef), dimension(:,:), optional, intent(out) ::  &
     &     kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jposize,jpitpsize
      INTEGER :: indobs,inddbs
      CHARACTER(len=word80) :: cdrec1,kform
      INTEGER :: jrec,jitp
      CHARACTER(len=4) :: versobs
      BIGREAL4, dimension(:), allocatable :: ptabo
      INTEGER :: jpiobs,jpjobs,jpkobs,jptobs
      INTEGER :: jpxendobs,jpyendobs,jpolocobs,jpitplocobs
      INTEGER :: allocok
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../readcfgobs/readcfgfileobs'
         WRITE(numout,*) '    ==> READING file ',kfninobs(1:lenv(kfninobs))
      ENDIF
!
! Set input array sizes and observation indices
! Check coherence of input array sizes
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
!
! -1.- Open direct access observation file
! ----------------------------------------
!
      irecl  = ibloc*(1000/ibloc+1)*jpbyt4
      CALL openfile (numfil,kfninobs,clold,clunf,cldir,irecl)
!
! -2.- Read observation file header
! ---------------------------------
!
      versobs='@obs'
      READ(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost)  &
     &     versobs,cdrec1,irecl,jpiobs,jpjobs,jpkobs,jptobs, &
     &     jpxendobs,jpyendobs,jpolocobs,jpitplocobs
      IF (versobs.NE.'@obs') GOTO 104
!
! -3.- Close observation file
! ---------------------------
!
      CLOSE (UNIT=numfil)
!
! -4.- Control print
! ------------------
!
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
!
! -5.- Coherence test
! -------------------
!
      IF (jpiobs.NE.dta_jpi(indobs)) GOTO 105
      IF (jpjobs.NE.dta_jpj(indobs)) GOTO 105
      IF (jpkobs.NE.dta_jpk(indobs)) GOTO 105
      IF (jpkobs.NE.dta_jpk(indobs)) GOTO 105
      IF (jpxendobs.NE.jpxend) GOTO 105
      IF (jpyendobs.NE.jpyend) GOTO 105
      IF (jpolocobs.NE.jposize) GOTO 105
      IF (jpitplocobs.NE.jpitpsize) GOTO 105
!
! -6.- Re-open file with the right record length
! ----------------------------------------------
!
      CALL openfile (numfil,kfninobs,clold,clunf,cldir,irecl)
!
! -7.- Read observation file
! --------------------------
!
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
! Read observation longitudes
      jrec=2
      IF (kflagcfg.EQ.2) THEN
         READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
         kgridijkobs(:)%longi = FREAL(ptabo(:))
      ENDIF
!
! Read observation latitudes
      jrec=jrec+1
      IF (kflagcfg.EQ.2) THEN
         READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
         kgridijkobs(:)%latj = FREAL(ptabo(:))
      ENDIF
!
! Read observation depths
      jrec=jrec+1
      IF (kflagcfg.EQ.2) THEN
         READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
         kgridijkobs(:)%levk = FREAL(ptabo(:))
      ENDIF
!
! Do not read observation values
      jrec=jrec+1
!
! Read associated observation error (obsolete)
      jrec=jrec+1
      IF (kflagcfg.EQ.1) THEN
         READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
         kvectorms(:) = FREAL(ptabo(:))
      ENDIF
!
! Read observation operator
      IF (kflagcfg.EQ.3) THEN
         DO jitp=1,jpitpsize
!
            jrec=jrec+1
            READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
            kposcoefobs(:,jitp)%pos = FREAL(ptabo(:))
!
            jrec=jrec+1
            READ(numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
            kposcoefobs(:,jitp)%coef = FREAL(ptabo(:))
!     
         ENDDO
      ENDIF
!
! -8.- Close observation file
! ---------------------------
!
      CLOSE (UNIT=numfil)
!
! --- deallocation
      IF (allocated(ptabo)) deallocate(ptabo)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioobs','readcfgfileobs')
 1001 CALL printerror2(0,1001,3,'lioobs','readcfgfileobs')
!
 101  WRITE (texterror,*) 'Error reading observation file, iost=',iost
      CALL printerror2(0,101,3,'lioobs','readcfgfileobs',comment=texterror)
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,1,'lioobs','readcfgfileobs',comment=texterror)
 104  WRITE (texterror,*) 'Bad observation file'
      CALL printerror2(0,104,3,'lioobs','readcfgfileobs',comment=texterror)
 105  WRITE (texterror,*) 'Configuration in .obs file header', &
     &     ' incompatible with SESAM configuration'
      CALL printerror2(0,105,3,'lioobs','readcfgfileobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrfileobs(kfnoutobs,kjobs)
!---------------------------------------------------------------------
!
!  Purpose : Write empty '.obs' file
!  -------
!  Method : Write first record of direct access binary file
!  ------
!  Input :  kfnoutobs : filename
!  -----    kjobs    : observed variable index
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpyend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutobs
      INTEGER, intent(in) :: kjobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jposize,jpitpsize
      INTEGER :: indobs,inddbs
      CHARACTER(len=word80) :: cdrec1,kform
      INTEGER :: jrec
      CHARACTER(len=4) :: versobs
      INTEGER :: jpiobs,jpjobs,jpkobs,jptobs
      INTEGER :: jpxendobs,jpyendobs,jpolocobs,jpitplocobs
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writeobs/writehdrfileobs'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutobs(1:lenv(kfnoutobs))
      ENDIF
!
! Set output array sizes and observation indices
      jposize = 0
      jpitpsize = 0
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
!
! -1.- Open direct access observation file
! ----------------------------------------
!
      irecl  = ibloc*((1000)/ibloc+1)*jpbyt4
      CALL openfile (numfil,kfnoutobs,clunk,clunf,cldir,irecl)
!
! -2.- Write observation file header
! ----------------------------------
!
      WRITE (cdrec1,'(2a)') 'Observation vector from database: ', &
     &     obs_nam(indobs,inddbs)(1:lenv(obs_nam(indobs,inddbs)))
!
      versobs = '@obs'
      jpiobs = dta_jpi(indobs)
      jpjobs = dta_jpj(indobs)
      jpkobs = dta_jpk(indobs)
      jptobs = dta_jpt(indobs)
      jpxendobs = jpxend
      jpyendobs = jpyend
      jpolocobs = jposize
      jpitplocobs = jpitpsize
!
      WRITE(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost)   &
     &     versobs,cdrec1,irecl,jpiobs,jpjobs,jpkobs,jptobs, &
     &     jpxendobs,jpyendobs,jpolocobs,jpitplocobs
!
! -3.- Close observation file
! ---------------------------
!
      CLOSE (UNIT=numfil)
!
! -4.- Control print
! ------------------
!
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
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioobs','writehdrfileobs')
!
 101  WRITE (texterror,*) 'Error writing observation file, iost=',iost
      CALL printerror2(0,101,3,'lioobs','writehdrfileobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writefileobs(kfnoutobs,kvecto,kvectorms,kgridijkobs, &
     &     kposcoefobs,kjobs)
!---------------------------------------------------------------------
!
!  Purpose : Write individual '.obs' file
!  -------
!  Method : Write direct access binary file
!  ------
!  Input : kfnoutobs   : filename
!  -----   kjobs       : observed variable index
!          kvecto      : vector of observed values in file
!          kvectorms   : associated error value (obsolete)
!          kgridijkobs : observation location (x,y,z)
!          kposcoefobs : observation operator (interpolation points
!                         and interpolation coefficients)
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpyend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutobs
      BIGREAL, dimension(:), intent(in) :: kvecto,kvectorms
      TYPE (type_gridijk), dimension(:), intent(in)  :: kgridijkobs
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
      INTEGER, intent(in) :: kjobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jposize,jpitpsize
      INTEGER :: indobs,inddbs
      CHARACTER(len=word80) :: cdrec1,kform
      INTEGER :: jrec,jitp,jo
      CHARACTER(len=4) :: versobs
      BIGREAL4, dimension(:), allocatable :: ptabo
      BIGREAL4 :: sxyo_moy,sxyo_ect
      INTEGER :: jpiobs,jpjobs,jpkobs,jptobs
      INTEGER :: jpxendobs,jpyendobs,jpolocobs,jpitplocobs
!----------------------------------------------------------------------
!
! Control print
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ../writeobs/writehdrfileobs'
         WRITE(numout,*) '    ==> WRITING file ',kfnoutobs(1:lenv(kfnoutobs))
      ENDIF
!
! Set input array sizes and observation indices
      jposize = size(kvecto,1)
      jpitpsize = size(kposcoefobs,2)
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
!
! Check coherence of input array sizes
      IF (jposize.NE.size(kvectorms,1)) GOTO 102
      IF (jposize.NE.size(kgridijkobs,1)) GOTO 102
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 102
      IF (jpitpsize.NE.obs_itp(indobs,inddbs)) GOTO 102
!
! -1.- Open direct access observation file
! ----------------------------------------
!
#if defined _NEC
      irecl = max(jposize*jpbyt4,2048)
#else
      irecl = ibloc*((jposize)/ibloc+1)*jpbyt4
#endif
      CALL openfile (numfil,kfnoutobs,clunk,clunf,cldir,irecl)
!
! -2.- Write observation file header
! ----------------------------------
!
      WRITE (cdrec1,'(2a)') 'Observation vector from database: ', &
     &     obs_nam(indobs,inddbs)(1:lenv(obs_nam(indobs,inddbs)))
!
      versobs = '@obs'
      jpiobs = dta_jpi(indobs)
      jpjobs = dta_jpj(indobs)
      jpkobs = dta_jpk(indobs)
      jptobs = dta_jpt(indobs)
      jpxendobs = jpxend
      jpyendobs = jpyend
      jpolocobs = jposize
      jpitplocobs = jpitpsize
!
      WRITE(UNIT=numfil,REC=1,ERR=101,IOSTAT=iost)   &
     &     versobs,cdrec1,irecl,jpiobs,jpjobs,jpkobs,jptobs, &
     &     jpxendobs,jpyendobs,jpolocobs,jpitplocobs
!
! -3.- Write observation file
! ---------------------------
!
      allocate ( ptabo(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
! Write observation longitudes
      jrec=2
      ptabo(:) = FREAL4(kgridijkobs(:)%longi)
      WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
!
! Write observation latitudes
      jrec=jrec+1
      ptabo(:) = FREAL4(kgridijkobs(:)%latj)
      WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
!
! Write observation depths
      jrec=jrec+1
      ptabo(:) = FREAL4(kgridijkobs(:)%levk)
      WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
!
! Write observation values
      jrec=jrec+1
      IF (lmoyect) THEN
         sxyo_moy=FREAL4(obs_moy(indobs,inddbs))
         sxyo_ect=FREAL4(obs_ect(indobs,inddbs))
         ptabo(:) = FREAL4(kvecto(:))*sxyo_ect+sxyo_moy
      ELSE
         ptabo(:) = FREAL4(kvecto(:))
      ENDIF
      WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
!
! Write associated observation error (obsolete)
      jrec=jrec+1
      IF (lmoyect) THEN
         sxyo_moy=FREAL4(obs_moy(indobs,inddbs))
         sxyo_ect=FREAL4(obs_ect(indobs,inddbs))
         ptabo(:) = FREAL4(kvectorms(:))*sxyo_ect+sxyo_moy
      ELSE
         ptabo(:) = FREAL4(kvectorms(:))
      ENDIF
      WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
!
! Write observation operator
      DO jitp=1,obs_itp(indobs,inddbs)
!
         jrec=jrec+1
         ptabo(:) = FREAL4(kposcoefobs(:,jitp)%pos)
         WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
!
         jrec=jrec+1
         ptabo(:) = FREAL4(kposcoefobs(:,jitp)%coef)
         WRITE(UNIT=numfil,REC=jrec,ERR=101,IOSTAT=iost) ptabo(:)
!
      ENDDO
!
! -4.- Close observation file
! ---------------------------
!
      CLOSE (UNIT=numfil)
!
! -5.- Control print
! ------------------
!
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
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioobs','writefileobs')
 1001 CALL printerror2(0,1001,3,'lioobs','writefileobs')
!
 101  WRITE (texterror,*) 'Error writing observation file, iost=',iost
      CALL printerror2(0,101,3,'lioobs','writefileobs',comment=texterror)
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,1,'lioobs','writefileobs',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE lioobs
