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
! ---                   HIOBAS.F90                                ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12  (C.E. Testut)                       ---
! --- modification : 99-05  (C.E. Testut)                       ---
! --- modification : 99-11  (J.M. Brankart)                     ---
! --- modification : 01-06  (C.E. Testut)                       ---
! --- modification : 03-03  (J.M. Brankart)                     ---
! --- modification : 07-09  (F. Castruccio)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE readbas      : Read Cx, Cy or Co object from input directory
! --- SUBROUTINE writebas     : Write Cx object to output directory
! --- SUBROUTINE writeyobas   : Write Cy or Co object to output directory
!
! --- SUBROUTINE readhdrzbas  : Read Cz configuration
! --- SUBROUTINE readptzbas   : Read Cz pointers
! --- SUBROUTINE readnbubzbas : Read list of Cz local covariance matrices
! --- SUBROUTINE writehdrzbas : Write Cz file headers
! --- SUBROUTINE writeptzbas  : Write Cz pointers
! --- 
! --- SUBROUTINE readinfobas  : Read information files from covariance directory
! --- SUBROUTINE writeinfobas : Write information files to covariance directory
!
! --- SUBROUTINE readvectb    : Read constraint b vector
! --- SUBROUTINE readmatzmb   : Read local constraint b vectors
!
! --- SUBROUTINE readscalbas  : Read scalar for every ensemble member
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiobas
      use mod_main
      use hioxyo
      use utilvalid
      use utilfiles
      IMPLICIT NONE
      PRIVATE

      PUBLIC readbas,writebas,writeyobas
      PUBLIC readhdrzbas,readptzbas,readnbubzbas,writehdrzbas
      PUBLIC writeptzbas,readinfobas,writeinfobas
      PUBLIC readvectb,readmatzmb,readscalbas

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readbas(dirnambas,kbasesr,kjnxyo,kjrbasdeb,kjrbasfin, &
     &     klectinfo,kflagxyo,kposcoefobs)
!---------------------------------------------------------------------
!
!  Purpose : Read Cx, Cy or Co reduced order covariance matrix objects
!  -------   from input directory
!
!  Method : Read ensemble of vectors (error modes) defining
!  ------   the reduced order covariance matrix
!
!  Input : dirnambas   : name of the directory
!  -----   kjnx        : block index
!          kjrbasdeb   : first vector to read (in the ensemble)
!          kjrbasfin   : last vector to read (in the ensemble)
!          klectinfo   : read or not header of input files
!          kflagxyo    : covariance type (1=Cx,2=Cy,3=Co)
!          kposcoefobs : observation operator (for Co objects only)
!
!  Output : kbasesr  : covariance object (s*r array)
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      BIGREAL, dimension(:,:), intent(out) :: kbasesr
      INTEGER , intent(in) :: kjrbasdeb,kjrbasfin,kjnxyo,kflagxyo
      LOGICAL, intent(in) :: klectinfo
      TYPE (type_poscoef), dimension(:,:), optional, intent(in) ::  &
     &     kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpssize,jprsize
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jr,jprbas
      LOGICAL :: lectinfo
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readbas :'
         WRITE(numout,*) '         read Cx, Cy or Co covariance matrix'
         WRITE(numout,*) '    ==> READING directory ', &
     &                                dirnambas(1:lenv(dirnambas))
      ENDIF
!
! Check validity of input directory
      jpssize = size(kbasesr,1)
      jprsize = size(kbasesr,2)
      IF (.NOT.(validextbas(dirnambas))) GOTO 102
      SELECT CASE (kflagxyo)
      CASE (1)
         IF (.NOT.(validextvarbas(dirnambas))) GOTO 102
      CASE (2)
         IF ((.NOT.(validextvarbas(dirnambas))) &
     &        .AND.(.NOT.(validextdtabas(dirnambas)))) GOTO 102
      CASE (3)
         IF (.NOT.(present(kposcoefobs))) GOTO 1000
         IF (jpssize.NE.size(kposcoefobs,1))  GOTO 1000
         IF ((.NOT.(validextvarbas(dirnambas))) &
     &        .AND.(.NOT.(validextdtabas(dirnambas))) &
     &        .AND.(.NOT.(validextobsbas(dirnambas)))) GOTO 102
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Loop on vectors in the ensemble
      serie=1
      lectinfo = klectinfo
      DO jr=kjrbasdeb,kjrbasfin
! Control print
         IF (nprint.GE.3) THEN
            WRITE(numout,*) '    ==> LOADING vector number ',jr
         ENDIF
! Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &        vctnam(1:lenv(vctnam))
! Load vector number jr
         SELECT CASE(kflagxyo)
         CASE (1,2)
            CALL readxyo(fname,kbasesr(:,jr), &
     &           kjnxyo,lectinfo,kflagxyo)
         CASE (3)
            CALL readxyo(fname,kbasesr(:,jr),kjnxyo,lectinfo, &
     &           kflagxyo,kposcoefobs(:,:))
         CASE DEFAULT
            GOTO 1000
         END SELECT
         lectinfo = .FALSE.
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','readbas')
!
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','readbas',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writebas(dirnambas,kbasexr,kjnx,kjrbasdeb, &
     &        kjrbasfin)
!---------------------------------------------------------------------
!
!  Purpose : Write Cx reduced order covariance
!  -------   matrix objects to input directory
!
!  Method : Write ensemble of vectors (error modes) defining
!  ------   the reduced order covariance matrix
!
!  Input : dirnambas   : name of the directory
!  -----   kbasexr     : covariance object (x*r array)
!          kjnx        : block index
!          kjrbasdeb   : first vector to read (in the ensemble)
!          kjrbasfin   : last vector to read (in the ensemble)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      BIGREAL, dimension(:,:), intent(in) :: kbasexr
      INTEGER , intent(in) :: kjrbasdeb,kjrbasfin,kjnx
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: lectinfo,serie,jr,jprbas
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writebas :'
         WRITE(numout,*) '         write Cx covariance matrix'
         WRITE(numout,*) '    ==> WRITING directory ', &
     &                                dirnambas(1:lenv(dirnambas))
      ENDIF
!
! Check validity of output directory
      IF (.NOT.(validextbas(dirnambas))) GOTO 102
      IF (.NOT.(validextvarbas(dirnambas))) GOTO 102
! Check validity of vector indices
      IF (kjrbasdeb.LT.0) GOTO 104
      IF (kjrbasfin.GT.size(kbasexr,2)) GOTO 104
!
! Loop on vectors in the ensemble
      serie=1
      DO jr=kjrbasdeb,kjrbasfin
! Control print
         IF (nprint.GE.3) THEN
            WRITE(numout,*) '    ==> WRITING vector number ',jr
         ENDIF
! Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &        vctnam(1:lenv(vctnam))
! Writing vector number jr
         CALL writevar (fname,kbasexr(:,jr),kjnx)
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','writebas')
!
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','writebas',comment=texterror)
 104  WRITE (texterror,*) 'Vector index out of ensemble: size=',size(kbasexr,2), &
     &     ' jr_first=',kjrbasdeb,' jr_last=',kjrbasfin
      CALL printerror2(0,104,1,'hiobas','writebas',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeyobas(dirnambas,kbasesr,kjrbasdeb, &
     &     kjrbasfin,kflagxyo,kvectsrms,kgridijkobs, &
     &     kposcoefobs)
!---------------------------------------------------------------------
!
!  Purpose : Write Cy or Co reduced order covariance
!  -------   matrix objects to output directory
!
!  Method : Write ensemble of vectors (error modes) defining
!  ------   the reduced order covariance matrix
!
!  Input : dirnambas   : name of the directory
!  -----   kbasesr     : covariance object (s*r array)
!          kjrbasdeb   : first vector to read (in the ensemble)
!          kjrbasfin   : last vector to read (in the ensemble)
!          kflagxyo    : covariance type (2=Cy,3=Co)
!          kvectsrms   : observation error (for Co vectors only) [obsolete]
!          kgridijkobs : observation location (for Co vectors only)
!          kposcoefobs : observation operator (for Co objects only)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      BIGREAL, dimension(:,:), intent(in) :: kbasesr
      INTEGER , intent(in) :: kjrbasdeb,kjrbasfin
      INTEGER, intent(in) :: kflagxyo
      BIGREAL, dimension(:), optional, intent(in) :: kvectsrms
      TYPE (type_gridijk), dimension(:), optional,  intent(in)  ::  &
     &     kgridijkobs
      TYPE (type_poscoef), dimension(:,:), optional,  intent(in) ::  &
     &     kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jr,jprbas
      INTEGER :: jpssize,jprsize
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writeyobas :'
         WRITE(numout,*) '         write Cy or Co covariance matrix'
         WRITE(numout,*) '    ==> WRITING directory ', &
     &                                dirnambas(1:lenv(dirnambas))
      ENDIF
!
! Check validity of output directory and routine arguments
      jpssize = size(kbasesr,1)
      jprsize = size(kbasesr,2)
      IF (.NOT.(validextbas(dirnambas))) GOTO 102
      SELECT CASE (kflagxyo)
      CASE (1)
         GOTO 1000
      CASE (2)
         IF (.NOT.(validextdtabas(dirnambas))) GOTO 102
      CASE (3)
         IF (.NOT.(validextobsbas(dirnambas))) GOTO 102
         IF (.NOT.(present(kvectsrms))) GOTO 1000
         IF (.NOT.(present(kgridijkobs))) GOTO 1000
         IF (.NOT.(present(kposcoefobs))) GOTO 1000
         IF (jpssize.NE.size(kvectsrms,1))  GOTO 1000
         IF (jpssize.NE.size(kgridijkobs,1))  GOTO 1000
         IF (jpssize.NE.size(kposcoefobs,1))  GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Loop on vectors in the ensemble
      serie=1
      DO jr=kjrbasdeb,kjrbasfin
! Control print
         IF (nprint.GE.3) THEN
            WRITE(numout,*) '    ==> WRITING vector number ',jr
         ENDIF
! Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &        vctnam(1:lenv(vctnam))
! Writing vector number jr
         SELECT CASE (kflagxyo)
         CASE (1)
            GOTO 1000
         CASE (2)
            CALL writeyo(fname,kbasesr(:,jr),kflagxyo)
         CASE (3)
            CALL writeyo(fname,kbasesr(:,jr),kflagxyo, &
     &           kvectsrms=kvectsrms,kgridijkobs=kgridijkobs, &
     &           kposcoefobs=kposcoefobs)
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','writeyobas')
!
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','writeyobas',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readhdrzbas(dirnambas,kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz,kjrbasdeb,kjrbasfin)
!---------------------------------------------------------------------
!
!  Purpose : Read Cz configuration
!  -------
!  Method : Read header from the first vector of the ensemble
!  ------   
!  Input :  dirnambas   : name of the directory
!  -----    kjrbasdeb   : first vector to read (in the ensemble)
!           kjrbasfin   : last vector to read (in the ensemble)
!
!  Output : kjpz        : number of local data sections
!  ------   kjpbub      : number of different local data sections
!           kzon_jpi    : size of local data sections (along X dimension)
!           kzon_jpj    : size of local data sections (along Y dimension)
!           kzon_jpk    : size of local data sections (along Z dimension)
!           kzon_jpt    : size of local data sections (along T dimension)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use hiozon , only : evalhdrzon
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, intent(out) :: kzon_jpi,kzon_jpj,kzon_jpk, &
     &                        kzon_jpt,kjpbub,kjpz
      INTEGER, intent(in) :: kjrbasdeb,kjrbasfin
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jprbas
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readhdrzbas :'
         WRITE(numout,*) '         read Cz configuration'
         WRITE(numout,*) '    ==> READING in directory ', &
     &                                dirnambas(1:lenv(dirnambas))
      ENDIF
!
! Check validity of input directory
      IF (.NOT.(validextzonbas(dirnambas))) GOTO 102
!
! Build file name for vector number kjrbasdeb
      serie=1
      CALL fildirbas (vctnam,dirnambas,jprbas,kjrbasdeb,serie)
      WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &     vctnam(1:lenv(vctnam))
! Read header from vector number kjrbasdeb
      CALL evalhdrzon(fname,kzon_jpi,kzon_jpj,kzon_jpk, &
     &        kzon_jpt,kjpbub,kjpz)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','readhdrzbas')
!
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','readhdrzbas',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readptzbas(dirnambas,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime,kjrbasdeb,kjrbasfin)
!---------------------------------------------------------------------
!  Purpose : Read Cz pointers
!  -------
!  Method : Read pointers from the first vector of the ensemble
!  ------   
!  Input :  dirnambas   : name of the directory
!  -----    kjrbasdeb   : first vector to read (in the ensemble)
!           kjrbasfin   : last vector to read (in the ensemble)
!
!  Output : ptbubidx   :
!  ------   ptdtalon   :
!           ptdtalat   :
!           ptdtadepth :    Cz pointers, see description
!           ptdtatime  :    in 'mod_spacexyo.F90'
!           ptbublon   :
!           ptbublat   :
!           ptbubdepth :
!           ptbubtime  :
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use hiozon , only : readptzon
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, dimension (:,:), intent(out) :: ptbubidx
      INTEGER, dimension (:,:), intent(out) :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), intent(out) :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), intent(out) :: ptbublon, ptbublat
      INTEGER, dimension (:,:), intent(out) :: ptbubdepth, ptbubtime
      INTEGER , intent(in) :: kjrbasdeb,kjrbasfin
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jprbas
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readptzbas :'
         WRITE(numout,*) '         read Cz pointers'
         WRITE(numout,*) '    ==> READING in directory ', &
     &                                dirnambas(1:lenv(dirnambas))
      ENDIF
!
! Check validity of input directory
      IF (.NOT.(validextzonbas(dirnambas))) GOTO 102
!
! Build file name for vector number kjrbasdeb
      serie=1
      CALL fildirbas (vctnam,dirnambas,jprbas,kjrbasdeb,serie)
      WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &     vctnam(1:lenv(vctnam))
! Read pointers from vector number kjrbasdeb
      CALL readptzon(fname,ptbubidx,ptdtalon,ptdtalat,ptdtadepth, &
     &              ptdtatime,ptbublon,ptbublat,ptbubdepth,ptbubtime)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','readptzbas')
!
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','readptzbas',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readnbubzbas(dirnambas,kbasenbubr,kptbubidx, &
     &           kjrbasdeb,kjrbasfin,klectinfo)
!---------------------------------------------------------------------
!
!  Purpose : Read Cz reduced order covariance matrices
!  -------   for a list of local data sections
!
!  Method : Read list of local data sections from all
!  ------   vectors of the ensemble
!
!  Input :  dirnambas  : name of the directory
!  -----    kjrbasdeb  : first vector to read (in the ensemble)
!           kjrbasfin  : last vector to read (in the ensemble)
!           klectinfo  : read or not header of input files
!           kptbubidx  : list of local data sections to read
!
!  Output : kbasenbubr : ensemble of list of local data sections
!  ------                (zon_jpi*zon_jpi*zon_jpk*zon_jpt*
!                        nbub*r)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use hiozon , only : readnbubzon
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      CHARACTER(len=*), intent(in) :: dirnambas
      BIGREAL, dimension(:,:,:,:,:,:), intent(out) :: kbasenbubr
      INTEGER, dimension(:), intent(in) :: kptbubidx
      INTEGER , intent(in) :: kjrbasdeb,kjrbasfin
      LOGICAL, intent(in) :: klectinfo
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jr,jprbas
!----------------------------------------------------------------------
!
      IF (nprint.GE.3) THEN
         WRITE(numout,*) '*** ROUTINE : ./readnbubzbas :'
         WRITE(numout,*) '         read list Cz local data sections'
         WRITE(numout,*) '    ==> READING in directory ', &
     &                                dirnambas(1:lenv(dirnambas))
      ENDIF
!
! Check validity of input directory and routine arguments
      IF (klectinfo) THEN
         IF (size(kbasenbubr,1).NE.zon_jpi) GOTO 1000
         IF (size(kbasenbubr,2).NE.zon_jpj) GOTO 1000
         IF (size(kbasenbubr,3).NE.zon_jpk) GOTO 1000
         IF (size(kbasenbubr,4).NE.zon_jpt) GOTO 1000
         IF (size(kptbubidx,1).NE.size(kbasenbubr,5)) GOTO 1000
         IF (kjrbasfin.GT.size(kbasenbubr,6)) GOTO 1000
         IF (.NOT.(validextzonbas(dirnambas))) GOTO 102
      ENDIF
!
! Loop on vectors in the ensemble
      serie=1
      DO jr=kjrbasdeb,kjrbasfin
! Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &        vctnam(1:lenv(vctnam))
! Load list of local data sections for vector number jr
         CALL readnbubzon(fname,kptbubidx,kbasenbubr(:,:,:,:,:,jr), &
     &                    klectinfo)
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','readnbubzbas')
!
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','readnbubzbas',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrzbas(dirnambas,kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz,kjrbasdeb,kjrbasfin)
!---------------------------------------------------------------------
!
!  Purpose : Write Cz file headers
!  -------
!  Method : Loop on Cz members and write header of every file
!  ------   in Cz directory
!
!  Input : dirnambas   : name of the directory
!  -----   kjrbasdeb   : first vector to write (in the ensemble)
!          kjrbasfin   : last vector to read (in the ensemble)
!          kjpz        : number of local data sections
!          kjpbub      : number of different local data sections
!          kzon_jpi    : size of local data sections (along X dimension)
!          kzon_jpj    : size of local data sections (along Y dimension)
!          kzon_jpk    : size of local data sections (along Z dimension)
!          kzon_jpt    : size of local data sections (along T dimension)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use hiozon , only : writehdrzon
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, intent(in) :: kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz
      INTEGER, intent(in) :: kjrbasdeb,kjrbasfin
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jr,jprbas
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writehdrzbas :'
         WRITE(numout,*) '         write Cz file headers'
         WRITE(numout,*) '    ==> WRITING in directory ', &
     &                                dirnambas(1:lenv(dirnambas))
      ENDIF
!
! Check validity of output directory
      IF (.NOT.(validextbas(dirnambas))) GOTO 102
      IF (.NOT.(validextzonbas(dirnambas))) GOTO 102
!
! Loop on vectors in the ensemble

      serie=1
      DO jr=kjrbasdeb,kjrbasfin
! Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &        vctnam(1:lenv(vctnam))
! Write header of vector number jr
         CALL writehdrzon (fname,kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz)
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','writehdrzbas')
!
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','writehdrzbas',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeptzbas(dirnambas,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime,kjrbasdeb,kjrbasfin)
!---------------------------------------------------------------------
!
!  Purpose : Write Cz file pointers
!  -------
!  Method : Loop on Cz members and write pointers of every file
!  ------   of Cz directory
!
!  Input : dirnambas  : name of the directory
!  -----   kjrbasdeb  : first vector to write (in the ensemble)
!          kjrbasfin  : last vector to read (in the ensemble)
!          ptbubidx   :
!          ptdtalon   :
!          ptdtalat   :
!          ptdtadepth :    Cz pointers, see description
!          ptdtatime  :    in 'mod_spacexyo.F90'
!          ptbublon   :
!          ptbublat   :
!          ptbubdepth :
!          ptbubtime  :
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use hiozon , only : writeptzon
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, dimension (:,:), intent(in) :: ptbubidx
      INTEGER, dimension (:,:), intent(in) :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), intent(in) :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), intent(in) :: ptbublon, ptbublat
      INTEGER, dimension (:,:), intent(in) :: ptbubdepth, ptbubtime
      INTEGER, intent(in) :: kjrbasdeb,kjrbasfin
!----------------------------------------------------------------------
! local declarations
! ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jr,jprbas
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readptzbas :'
         WRITE(numout,*) '         write Cz pointers'
         WRITE(numout,*) '    ==> WRITING in directory ', &
     &                                dirnambas(1:lenv(dirnambas))
      ENDIF
!
! Check validity of input directory
      IF (.NOT.(validextbas(dirnambas))) GOTO 102
      IF (.NOT.(validextzonbas(dirnambas))) GOTO 102
!
! Loop on vectors in the ensemble
      serie=1
      DO jr=kjrbasdeb,kjrbasfin
! Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &        vctnam(1:lenv(vctnam))
! Write pointers for vector number jr
         CALL writeptzon (fname,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime)
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','writeptzbas')
!
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','writeptzbas',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readinfobas(dirnambas,kvalbase,kvalp,kvctp)
!---------------------------------------------------------------------
!
!  Purpose : Read information files from covariance directory
!  --------
!  Method : Open information files and read covariance type,
!  ------   eigenvalues and eigenvectors
!
!  Input :  dirnambas  : name of the directory
!  -----
!  Output : kvalbase   : covariance type
!  ------   kvalp      : eigenvalues
!           kvctp      : eigenvectors
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, intent(out) :: kvalbase
      BIGREAL, dimension(:), optional, intent(out) :: kvalp
      BIGREAL, dimension(:,:), optional, intent(out) :: kvctp
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpr1size,jpr2size
      INTEGER :: jr,jr1,serie,numjr,jprbas
      CHARACTER(len=bgword) :: infonam,fname,kform
      LOGICAL :: existence
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readinfobas :'
         WRITE(numout,*) '         read Cx, Cy or Co information files'
         WRITE(numout,*) '    ==> READING in directory ', &
     &                                dirnambas(1:lenv(dirnambas))
      ENDIF
!
! Check input arguments
      IF (present(kvctp)) THEN
         jpr1size = size(kvctp,1)
         jpr2size = size(kvctp,2)
         IF ((present(kvalp)).AND. &
     &        (jpr2size.NE.size(kvalp,1))) GOTO 1000
      ENDIF
!
! -1.- Read the 'valbase' file (covariance type)
! ----------------------------------------------
! Build file name
      serie=0
      numjr=1
      CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
      WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &                        infonam(1:lenv(infonam))
      INQUIRE (FILE=fname,EXIST=existence)
!
      IF (existence) THEN
!
! Read informations
         CALL openfile(numfil,fname)
         READ(UNIT=numfil,FMT='(i3.3)',ERR=101) kvalbase
         CLOSE (UNIT=numfil)
!
         IF (kvalbase.GT.0) THEN
!
! -2.- Read the 'valp' file (eigenvalues)
! ---------------------------------------
!
            IF (present(kvalp)) THEN
! Build file name
               numjr=2
               CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
               WRITE(fname,'("./",A,"/",A)')  &
     &              dirnambas(1:lenv(dirnambas)), &
     &              infonam(1:lenv(infonam))
! Read eigenvalues
               CALL openfile(numfil,fname)
               DO jr=1,jpr1size
                  READ(UNIT=numfil,FMT='(E12.6E2)',ERR=101) kvalp(jr)
               ENDDO
               CLOSE (UNIT=numfil)
            ENDIF
!
! -2.- Read the 'vctp' file (eigenvectors)
! ----------------------------------------
!
            IF (present(kvctp)) THEN
! Build file name
               numjr=3
               CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
               WRITE(fname,'("./",A,"/",A)')  &
     &              dirnambas(1:lenv(dirnambas)), &
     &              infonam(1:lenv(infonam))
! Read eigenvectors
               CALL openfile(numfil,fname)
               DO jr=1,jpr2size
               DO jr1=1,jpr1size
                  READ(UNIT=numfil,FMT='(E12.6E2)',ERR=101)  &
     &                 kvctp(jr1,jr)
               ENDDO
               ENDDO
               CLOSE (UNIT=numfil)
            ENDIF
         ENDIF
      ELSE
         kvalbase=0
         IF (present(kvalp)) kvalp(:) = FREAL(0.0)
         IF (present(kvctp)) kvctp(:,:) = FREAL(0.0)
      ENDIF
!
      IF (nprint.GE.2) THEN
         kform='(8x,a,i3)'
         WRITE(numout,kform) ' - Covariance type: ',kvalbase
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','readinfobas')
!
 101  WRITE (texterror,*) 'error reading ASCII file, iost=',iost
      CALL printerror2(0,101,3,'hiobas','readinfobas', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeinfobas(dirnambas,kvalbase,kvalp,kvctp)
!---------------------------------------------------------------------
!
!  Purpose : Write information files to covariance directory
!  -------
!  Method : Open information files and write covariance type,
!  ------   eigenvalues and eigenvectors
!
!  Input :  dirnambas  : name of the directory
!  -----    kvalbase   : covariance type
!           kvalp      : eigenvalues
!           kvctp      : eigenvectors
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, intent(in) :: kvalbase
      BIGREAL, dimension(:), optional, intent(in) :: kvalp
      BIGREAL, dimension(:,:), optional, intent(in) :: kvctp
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpr1size,jpr2size,jr,jr1,serie,numjr,jprbas
      CHARACTER(len=bgword) :: infonam,fname
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writeinfobas :'
         WRITE(numout,*) '         write Cx, Cy or Co information files'
         WRITE(numout,*) '    ==> WRITING in directory ', &
     &                                dirnambas(1:lenv(dirnambas))
      ENDIF
!
! Check input arguments
      IF (present(kvctp)) THEN
         jpr1size = size(kvctp,1)
         jpr2size = size(kvctp,2)
         IF ((present(kvalp)).AND.(jpr2size.NE.size(kvalp,1))) GOTO 1000
      ENDIF
!
! -1.- Write the 'valbase' file (covariance type)
! -----------------------------------------------
! Build file name
      serie=0
      numjr=1
      CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
      WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &                        infonam(1:lenv(infonam))
! Write informations
      CALL openfile(numfil,fname,kstatus=clunk)
      WRITE(UNIT=numfil,FMT='(i3.3)') kvalbase
      CLOSE (UNIT=numfil)
!
! -2.- Write the 'valp' file (eigenvalues)
! ----------------------------------------
!
      IF (present(kvalp)) THEN
! Build file name
         numjr=2
         CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &        infonam(1:lenv(infonam))
! Write eigenvalues
         CALL openfile(numfil,fname,kstatus=clunk)
         DO jr=1,jpr2size
            WRITE(UNIT=numfil,FMT='(E12.6E2)') kvalp(jr)
         ENDDO
         CLOSE (UNIT=numfil)
      ENDIF
!
! -3.- Write the 'vctp' file (eigenvectors)
! -----------------------------------------
!
      IF (present(kvctp)) THEN
! Build file name
         numjr=3
         CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)), &
     &        infonam(1:lenv(infonam))
! Write eigenvectors
         CALL openfile(numfil,fname,kstatus=clunk)
         DO jr=1,jpr2size
         DO jr1=1,jpr1size
            WRITE(UNIT=numfil,FMT='(E12.6E2)') kvctp(jr1,jr)
         ENDDO
         ENDDO
         CLOSE (UNIT=numfil)
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','writeinfobas')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readvectb(dirnambas,kvectb)
!---------------------------------------------------------------------
!
!  Purpose : Read constraint b vector
!  --------
!  Method : Build b vector file name and
!  ------   read b vector
!
!  Input :  dirnambas  : name of the directory
!  -----
!  Output : kvectb     : b vector
!  ------
!---------------------------------------------------------------------
! modules
! =======
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(out) :: kvectb
      CHARACTER(len=*), intent(in) :: dirnambas
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpmsize, jm
      CHARACTER(len=bgword) :: vectbname
!----------------------------------------------------------------------
!
      jpmsize=size(kvectb,1)
!
      vectbname='vectb'
      vectbname=dirnambas(1:lenv(dirnambas)) // '/' &
     &                  //  vectbname(1:lenv(vectbname))
!
      OPEN (unit=10,file=vectbname)
      DO jm=1,jpmsize
        READ(10,*) kvectb(jm)
      ENDDO
      CLOSE(10)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','readvectb')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmatzmb(dirnambas,kmatzmb)
!---------------------------------------------------------------------
!
!  Purpose : Read local constraint b vectors
!  --------
!  Method : Build b vector file name and
!  ------   read b vector
!
!  Input :  dirnambas  : name of the directory
!  -----
!  Output : kmatzmb    : local b vectors
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo, only : texterror
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:,:), intent(out) :: kmatzmb
      CHARACTER(len=*), intent(in) :: dirnambas
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpmsize, jpzsize, jm, jz
      CHARACTER(len=bgword) :: vectbname
!----------------------------------------------------------------------
!
      jpzsize=size(kmatzmb,1)
      jpmsize=size(kmatzmb,2)
!
      vectbname='matzmb'
      vectbname=dirnambas(1:lenv(dirnambas)) // '/' &
     &                  //  vectbname(1:lenv(vectbname))
!
      OPEN (unit=10,file=vectbname)
      DO jz=1,jpzsize
        READ(10,*) ( kmatzmb(jz,jm), jm=1,jpmsize )
      ENDDO
      CLOSE(10)
!
      RETURN
!
! --- error management
!
 101  WRITE (texterror,*) 'Error reading file: ',vectbname
      CALL printerror2(0,102,3,'hiobas','readmatzmb',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readscalbas(dirnambas,filename,scal)
!---------------------------------------------------------------------
!
!  Purpose : Read scalar for every ensemble member
!  --------
!  Method : Read from file in ensemble directory
!  ------
!
!  Input :  dirnambas  : name of the directory
!  -----    filename   : name of the file
!
!  Output : scal       : scalar array
!  ------
!---------------------------------------------------------------------
! modules
! =======
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(out) :: scal
      CHARACTER(len=*), intent(in) :: dirnambas
      CHARACTER(len=*), intent(in) :: filename
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jppsize, jp
      CHARACTER(len=bgword) :: fname
!----------------------------------------------------------------------
!
      jppsize=size(scal,1)
!
      fname=dirnambas(1:lenv(dirnambas)) // '/' &
     &                  //  filename(1:lenv(filename))
!
      OPEN (unit=10,file=fname)
      DO jp=1,jppsize
        READ(10,*) scal(jp)
      ENDDO
      CLOSE(10)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiobas','readscalbas')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiobas
