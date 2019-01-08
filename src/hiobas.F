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
C ---                   HIOBAS.F                                ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12  (C.E. Testut)                       ---
C --- modification : 99-05  (C.E. Testut)                       ---
C --- modification : 99-11  (J.M. Brankart)                     ---
C --- modification : 01-06  (C.E. Testut)                       ---
C --- modification : 03-03  (J.M. Brankart)                     ---
C --- modification : 07-09  (F. Castruccio)                     ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE readbas      : Read Cx, Cy or Co object from input directory
C --- SUBROUTINE writebas     : Write Cx object to output directory
C --- SUBROUTINE writeyobas   : Write Cy or Co object to output directory
C
C --- SUBROUTINE readhdrzbas  : Read Cz configuration
C --- SUBROUTINE readptzbas   : Read Cz pointers
C --- SUBROUTINE readnbubzbas : Read list of Cz local covariance matrices
C --- SUBROUTINE writehdrzbas : Write Cz file headers
C --- SUBROUTINE writeptzbas  : Write Cz pointers
C --- 
C --- SUBROUTINE readinfobas  : Read information files from covariance directory
C --- SUBROUTINE writeinfobas : Write information files to covariance directory
C
C --- SUBROUTINE readvectb    : Read constraint b vector
C --- SUBROUTINE readmatzmb   : Read local constraint b vectors
C
C --- SUBROUTINE readscalbas  : Read scalar for every ensemble member
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readbas(dirnambas,kbasesr,kjnxyo,kjrbasdeb,kjrbasfin,
     $     klectinfo,kflagxyo,kposcoefobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read Cx, Cy or Co reduced order covariance matrix objects
CCC  -------   from input directory
CCC
CCC  Method : Read ensemble of vectors (error modes) defining
CCC  ------   the reduced order covariance matrix
CCC
CCC  Input : dirnambas   : name of the directory
CCC  -----   kjnx        : block index
CCC          kjrbasdeb   : first vector to read (in the ensemble)
CCC          kjrbasfin   : last vector to read (in the ensemble)
CCC          klectinfo   : read or not header of input files
CCC          kflagxyo    : covariance type (1=Cx,2=Cy,3=Co)
CCC          kposcoefobs : observation operator (for Co objects only)
CCC
CCC  Output : kbasesr  : covariance object (s*r array)
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      BIGREAL, dimension(:,:), intent(out) :: kbasesr
      INTEGER , intent(in) :: kjrbasdeb,kjrbasfin,kjnxyo,kflagxyo
      LOGICAL, intent(in) :: klectinfo
      TYPE (type_poscoef), dimension(:,:), optional, intent(in) :: 
     $     kposcoefobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpssize,jprsize
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jr,jprbas
      LOGICAL :: lectinfo
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readbas :'
         WRITE(numout,*) '         read Cx, Cy or Co covariance matrix'
         WRITE(numout,*) '    ==> READING directory ',
     $                                dirnambas(1:lenv(dirnambas))
      ENDIF
C
C Check validity of input directory
      jpssize = size(kbasesr,1)
      jprsize = size(kbasesr,2)
      IF (.NOT.(validextbas(dirnambas))) GOTO 102
      SELECT CASE (kflagxyo)
      CASE (1)
         IF (.NOT.(validextvarbas(dirnambas))) GOTO 102
      CASE (2)
         IF ((.NOT.(validextvarbas(dirnambas)))
     $        .AND.(.NOT.(validextdtabas(dirnambas)))) GOTO 102
      CASE (3)
         IF (.NOT.(present(kposcoefobs))) GOTO 1000
         IF (jpssize.NE.size(kposcoefobs,1))  GOTO 1000
         IF ((.NOT.(validextvarbas(dirnambas)))
     $        .AND.(.NOT.(validextdtabas(dirnambas)))
     $        .AND.(.NOT.(validextobsbas(dirnambas)))) GOTO 102
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Loop on vectors in the ensemble
      serie=1
      lectinfo = klectinfo
      DO jr=kjrbasdeb,kjrbasfin
C Control print
         IF (nprint.GE.3) THEN
            WRITE(numout,*) '    ==> LOADING vector number ',jr
         ENDIF
C Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $        vctnam(1:lenv(vctnam))
C Load vector number jr
         SELECT CASE(kflagxyo)
         CASE (1,2)
            CALL readxyo(fname,kbasesr(:,jr),
     $           kjnxyo,lectinfo,kflagxyo)
         CASE (3)
            CALL readxyo(fname,kbasesr(:,jr),kjnxyo,lectinfo,
     $           kflagxyo,kposcoefobs(:,:))
         CASE DEFAULT
            GOTO 1000
         END SELECT
         lectinfo = .FALSE.
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','readbas')
C
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','readbas',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writebas(dirnambas,kbasexr,kjnx,kjrbasdeb,
     $        kjrbasfin)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write Cx reduced order covariance
CCC  -------   matrix objects to input directory
CCC
CCC  Method : Write ensemble of vectors (error modes) defining
CCC  ------   the reduced order covariance matrix
CCC
CCC  Input : dirnambas   : name of the directory
CCC  -----   kbasexr     : covariance object (x*r array)
CCC          kjnx        : block index
CCC          kjrbasdeb   : first vector to read (in the ensemble)
CCC          kjrbasfin   : last vector to read (in the ensemble)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      BIGREAL, dimension(:,:), intent(in) :: kbasexr
      INTEGER , intent(in) :: kjrbasdeb,kjrbasfin,kjnx
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: lectinfo,serie,jr,jprbas
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writebas :'
         WRITE(numout,*) '         write Cx covariance matrix'
         WRITE(numout,*) '    ==> WRITING directory ',
     $                                dirnambas(1:lenv(dirnambas))
      ENDIF
C
C Check validity of output directory
      IF (.NOT.(validextbas(dirnambas))) GOTO 102
      IF (.NOT.(validextvarbas(dirnambas))) GOTO 102
C Check validity of vector indices
      IF (kjrbasdeb.LT.0) GOTO 104
      IF (kjrbasfin.GT.size(kbasexr,2)) GOTO 104
C
C Loop on vectors in the ensemble
      serie=1
      DO jr=kjrbasdeb,kjrbasfin
C Control print
         IF (nprint.GE.3) THEN
            WRITE(numout,*) '    ==> WRITING vector number ',jr
         ENDIF
C Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $        vctnam(1:lenv(vctnam))
C Writing vector number jr
         CALL writevar (fname,kbasexr(:,jr),kjnx)
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','writebas')
C
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','writebas',comment=texterror)
 104  WRITE (texterror,*) 'Vector index out of ensemble: size=',size(kbasexr,2),
     $     ' jr_first=',kjrbasdeb,' jr_last=',kjrbasfin
      CALL printerror2(0,104,1,'hiobas','writebas',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeyobas(dirnambas,kbasesr,kjrbasdeb,
     $     kjrbasfin,kflagxyo,kvectsrms,kgridijkobs,
     $     kposcoefobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write Cy or Co reduced order covariance
CCC  -------   matrix objects to output directory
CCC
CCC  Method : Write ensemble of vectors (error modes) defining
CCC  ------   the reduced order covariance matrix
CCC
CCC  Input : dirnambas   : name of the directory
CCC  -----   kbasesr     : covariance object (s*r array)
CCC          kjrbasdeb   : first vector to read (in the ensemble)
CCC          kjrbasfin   : last vector to read (in the ensemble)
CCC          kflagxyo    : covariance type (2=Cy,3=Co)
CCC          kvectsrms   : observation error (for Co vectors only) [obsolete]
CCC          kgridijkobs : observation location (for Co vectors only)
CCC          kposcoefobs : observation operator (for Co objects only)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      BIGREAL, dimension(:,:), intent(in) :: kbasesr
      INTEGER , intent(in) :: kjrbasdeb,kjrbasfin
      INTEGER, intent(in) :: kflagxyo
      BIGREAL, dimension(:), optional, intent(in) :: kvectsrms
      TYPE (type_gridijk), dimension(:), optional,  intent(in)  :: 
     $     kgridijkobs
      TYPE (type_poscoef), dimension(:,:), optional,  intent(in) :: 
     $     kposcoefobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jr,jprbas
      INTEGER :: jpssize,jprsize
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writeyobas :'
         WRITE(numout,*) '         write Cy or Co covariance matrix'
         WRITE(numout,*) '    ==> WRITING directory ',
     $                                dirnambas(1:lenv(dirnambas))
      ENDIF
C
C Check validity of output directory and routine arguments
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
C
C Loop on vectors in the ensemble
      serie=1
      DO jr=kjrbasdeb,kjrbasfin
C Control print
         IF (nprint.GE.3) THEN
            WRITE(numout,*) '    ==> WRITING vector number ',jr
         ENDIF
C Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $        vctnam(1:lenv(vctnam))
C Writing vector number jr
         SELECT CASE (kflagxyo)
         CASE (1)
            GOTO 1000
         CASE (2)
            CALL writeyo(fname,kbasesr(:,jr),kflagxyo)
         CASE (3)
            CALL writeyo(fname,kbasesr(:,jr),kflagxyo,
     $           kvectsrms=kvectsrms,kgridijkobs=kgridijkobs,
     $           kposcoefobs=kposcoefobs)
         CASE DEFAULT
            GOTO 1000
         END SELECT
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','writeyobas')
C
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','writeyobas',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readhdrzbas(dirnambas,kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz,kjrbasdeb,kjrbasfin)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read Cz configuration
CCC  -------
CCC  Method : Read header from the first vector of the ensemble
CCC  ------   
CCC  Input :  dirnambas   : name of the directory
CCC  -----    kjrbasdeb   : first vector to read (in the ensemble)
CCC           kjrbasfin   : last vector to read (in the ensemble)
CCC
CCC  Output : kjpz        : number of local data sections
CCC  ------   kjpbub      : number of different local data sections
CCC           kzon_jpi    : size of local data sections (along X dimension)
CCC           kzon_jpj    : size of local data sections (along Y dimension)
CCC           kzon_jpk    : size of local data sections (along Z dimension)
CCC           kzon_jpt    : size of local data sections (along T dimension)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use hiozon , only : evalhdrzon
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, intent(out) :: kzon_jpi,kzon_jpj,kzon_jpk,
     $                        kzon_jpt,kjpbub,kjpz
      INTEGER, intent(in) :: kjrbasdeb,kjrbasfin
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jprbas
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readhdrzbas :'
         WRITE(numout,*) '         read Cz configuration'
         WRITE(numout,*) '    ==> READING in directory ',
     $                                dirnambas(1:lenv(dirnambas))
      ENDIF
C
C Check validity of input directory
      IF (.NOT.(validextzonbas(dirnambas))) GOTO 102
C
C Build file name for vector number kjrbasdeb
      serie=1
      CALL fildirbas (vctnam,dirnambas,jprbas,kjrbasdeb,serie)
      WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $     vctnam(1:lenv(vctnam))
C Read header from vector number kjrbasdeb
      CALL evalhdrzon(fname,kzon_jpi,kzon_jpj,kzon_jpk,
     $        kzon_jpt,kjpbub,kjpz)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','readhdrzbas')
C
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','readhdrzbas',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readptzbas(dirnambas,ptbubidx,ptdtalon,ptdtalat,
     $           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth,
     $           ptbubtime,kjrbasdeb,kjrbasfin)
CCC---------------------------------------------------------------------
CCC  Purpose : Read Cz pointers
CCC  -------
CCC  Method : Read pointers from the first vector of the ensemble
CCC  ------   
CCC  Input :  dirnambas   : name of the directory
CCC  -----    kjrbasdeb   : first vector to read (in the ensemble)
CCC           kjrbasfin   : last vector to read (in the ensemble)
CCC
CCC  Output : ptbubidx   :
CCC  ------   ptdtalon   :
CCC           ptdtalat   :
CCC           ptdtadepth :    Cz pointers, see description
CCC           ptdtatime  :    in 'mod_spacexyo.F'
CCC           ptbublon   :
CCC           ptbublat   :
CCC           ptbubdepth :
CCC           ptbubtime  :
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use hiozon , only : readptzon
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, dimension (:,:), intent(out) :: ptbubidx
      INTEGER, dimension (:,:), intent(out) :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), intent(out) :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), intent(out) :: ptbublon, ptbublat
      INTEGER, dimension (:,:), intent(out) :: ptbubdepth, ptbubtime
      INTEGER , intent(in) :: kjrbasdeb,kjrbasfin
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jprbas
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readptzbas :'
         WRITE(numout,*) '         read Cz pointers'
         WRITE(numout,*) '    ==> READING in directory ',
     $                                dirnambas(1:lenv(dirnambas))
      ENDIF
C
C Check validity of input directory
      IF (.NOT.(validextzonbas(dirnambas))) GOTO 102
C
C Build file name for vector number kjrbasdeb
      serie=1
      CALL fildirbas (vctnam,dirnambas,jprbas,kjrbasdeb,serie)
      WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $     vctnam(1:lenv(vctnam))
C Read pointers from vector number kjrbasdeb
      CALL readptzon(fname,ptbubidx,ptdtalon,ptdtalat,ptdtadepth,
     $              ptdtatime,ptbublon,ptbublat,ptbubdepth,ptbubtime)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','readptzbas')
C
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','readptzbas',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readnbubzbas(dirnambas,kbasenbubr,kptbubidx,
     $           kjrbasdeb,kjrbasfin,klectinfo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read Cz reduced order covariance matrices
CCC  -------   for a list of local data sections
CCC
CCC  Method : Read list of local data sections from all
CCC  ------   vectors of the ensemble
CCC
CCC  Input :  dirnambas  : name of the directory
CCC  -----    kjrbasdeb  : first vector to read (in the ensemble)
CCC           kjrbasfin  : last vector to read (in the ensemble)
CCC           klectinfo  : read or not header of input files
CCC           kptbubidx  : list of local data sections to read
CCC
CCC  Output : kbasenbubr : ensemble of list of local data sections
CCC  ------                (zon_jpi*zon_jpi*zon_jpk*zon_jpt*
CCC                        nbub*r)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use hiozon , only : readnbubzon
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ==================
      CHARACTER(len=*), intent(in) :: dirnambas
      BIGREAL, dimension(:,:,:,:,:,:), intent(out) :: kbasenbubr
      INTEGER, dimension(:), intent(in) :: kptbubidx
      INTEGER , intent(in) :: kjrbasdeb,kjrbasfin
      LOGICAL, intent(in) :: klectinfo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jr,jprbas
CC----------------------------------------------------------------------
C
      IF (nprint.GE.3) THEN
         WRITE(numout,*) '*** ROUTINE : ./readnbubzbas :'
         WRITE(numout,*) '         read list Cz local data sections'
         WRITE(numout,*) '    ==> READING in directory ',
     $                                dirnambas(1:lenv(dirnambas))
      ENDIF
C
C Check validity of input directory and routine arguments
      IF (klectinfo) THEN
         IF (size(kbasenbubr,1).NE.zon_jpi) GOTO 1000
         IF (size(kbasenbubr,2).NE.zon_jpj) GOTO 1000
         IF (size(kbasenbubr,3).NE.zon_jpk) GOTO 1000
         IF (size(kbasenbubr,4).NE.zon_jpt) GOTO 1000
         IF (size(kptbubidx,1).NE.size(kbasenbubr,5)) GOTO 1000
         IF (kjrbasfin.GT.size(kbasenbubr,6)) GOTO 1000
         IF (.NOT.(validextzonbas(dirnambas))) GOTO 102
      ENDIF
C
C Loop on vectors in the ensemble
      serie=1
      DO jr=kjrbasdeb,kjrbasfin
C Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $        vctnam(1:lenv(vctnam))
C Load list of local data sections for vector number jr
         CALL readnbubzon(fname,kptbubidx,kbasenbubr(:,:,:,:,:,jr),
     $                    klectinfo)
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','readnbubzbas')
C
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','readnbubzbas',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrzbas(dirnambas,kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz,kjrbasdeb,kjrbasfin)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write Cz file headers
CCC  -------
CCC  Method : Loop on Cz members and write header of every file
CCC  ------   in Cz directory
CCC
CCC  Input : dirnambas   : name of the directory
CCC  -----   kjrbasdeb   : first vector to write (in the ensemble)
CCC          kjrbasfin   : last vector to read (in the ensemble)
CCC          kjpz        : number of local data sections
CCC          kjpbub      : number of different local data sections
CCC          kzon_jpi    : size of local data sections (along X dimension)
CCC          kzon_jpj    : size of local data sections (along Y dimension)
CCC          kzon_jpk    : size of local data sections (along Z dimension)
CCC          kzon_jpt    : size of local data sections (along T dimension)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use hiozon , only : writehdrzon
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, intent(in) :: kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz
      INTEGER, intent(in) :: kjrbasdeb,kjrbasfin
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jr,jprbas
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writehdrzbas :'
         WRITE(numout,*) '         write Cz file headers'
         WRITE(numout,*) '    ==> WRITING in directory ',
     $                                dirnambas(1:lenv(dirnambas))
      ENDIF
C
C Check validity of output directory
      IF (.NOT.(validextbas(dirnambas))) GOTO 102
      IF (.NOT.(validextzonbas(dirnambas))) GOTO 102
C
C Loop on vectors in the ensemble

      serie=1
      DO jr=kjrbasdeb,kjrbasfin
C Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $        vctnam(1:lenv(vctnam))
C Write header of vector number jr
         CALL writehdrzon (fname,kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz)
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','writehdrzbas')
C
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','writehdrzbas',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeptzbas(dirnambas,ptbubidx,ptdtalon,ptdtalat,
     $           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth,
     $           ptbubtime,kjrbasdeb,kjrbasfin)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write Cz file pointers
CCC  -------
CCC  Method : Loop on Cz members and write pointers of every file
CCC  ------   of Cz directory
CCC
CCC  Input : dirnambas  : name of the directory
CCC  -----   kjrbasdeb  : first vector to write (in the ensemble)
CCC          kjrbasfin  : last vector to read (in the ensemble)
CCC          ptbubidx   :
CCC          ptdtalon   :
CCC          ptdtalat   :
CCC          ptdtadepth :    Cz pointers, see description
CCC          ptdtatime  :    in 'mod_spacexyo.F'
CCC          ptbublon   :
CCC          ptbublat   :
CCC          ptbubdepth :
CCC          ptbubtime  :
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use hiozon , only : writeptzon
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, dimension (:,:), intent(in) :: ptbubidx
      INTEGER, dimension (:,:), intent(in) :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), intent(in) :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), intent(in) :: ptbublon, ptbublat
      INTEGER, dimension (:,:), intent(in) :: ptbubdepth, ptbubtime
      INTEGER, intent(in) :: kjrbasdeb,kjrbasfin
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      CHARACTER(len=bgword) :: vctnam,fname
      INTEGER :: serie,jr,jprbas
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readptzbas :'
         WRITE(numout,*) '         write Cz pointers'
         WRITE(numout,*) '    ==> WRITING in directory ',
     $                                dirnambas(1:lenv(dirnambas))
      ENDIF
C
C Check validity of input directory
      IF (.NOT.(validextbas(dirnambas))) GOTO 102
      IF (.NOT.(validextzonbas(dirnambas))) GOTO 102
C
C Loop on vectors in the ensemble
      serie=1
      DO jr=kjrbasdeb,kjrbasfin
C Build file name for vector number jr
         CALL fildirbas (vctnam,dirnambas,jprbas,jr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $        vctnam(1:lenv(vctnam))
C Write pointers for vector number jr
         CALL writeptzon (fname,ptbubidx,ptdtalon,ptdtalat,
     $           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth,
     $           ptbubtime)
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','writeptzbas')
C
 102  WRITE (texterror,*) 'Invalid covariance directory name'
      CALL printerror2(0,102,3,'hiobas','writeptzbas',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readinfobas(dirnambas,kvalbase,kvalp,kvctp)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read information files from covariance directory
CCC  --------
CCC  Method : Open information files and read covariance type,
CCC  ------   eigenvalues and eigenvectors
CCC
CCC  Input :  dirnambas  : name of the directory
CCC  -----
CCC  Output : kvalbase   : covariance type
CCC  ------   kvalp      : eigenvalues
CCC           kvctp      : eigenvectors
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, intent(out) :: kvalbase
      BIGREAL, dimension(:), optional, intent(out) :: kvalp
      BIGREAL, dimension(:,:), optional, intent(out) :: kvctp
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpr1size,jpr2size
      INTEGER :: jr,jr1,serie,numjr,jprbas
      CHARACTER(len=bgword) :: infonam,fname,kform
      LOGICAL :: existence
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./readinfobas :'
         WRITE(numout,*) '         read Cx, Cy or Co information files'
         WRITE(numout,*) '    ==> READING in directory ',
     $                                dirnambas(1:lenv(dirnambas))
      ENDIF
C
C Check input arguments
      IF (present(kvctp)) THEN
         jpr1size = size(kvctp,1)
         jpr2size = size(kvctp,2)
         IF ((present(kvalp)).AND.
     $        (jpr2size.NE.size(kvalp,1))) GOTO 1000
      ENDIF
C
C -1.- Read the 'valbase' file (covariance type)
C ----------------------------------------------
C Build file name
      serie=0
      numjr=1
      CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
      WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $                        infonam(1:lenv(infonam))
      INQUIRE (FILE=fname,EXIST=existence)
C
      IF (existence) THEN
C
C Read informations
         CALL openfile(numfil,fname)
         READ(UNIT=numfil,FMT='(i3.3)',ERR=101) kvalbase
         CLOSE (UNIT=numfil)
C
         IF (kvalbase.GT.0) THEN
C
C -2.- Read the 'valp' file (eigenvalues)
C ---------------------------------------
C
            IF (present(kvalp)) THEN
C Build file name
               numjr=2
               CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
               WRITE(fname,'("./",A,"/",A)') 
     $              dirnambas(1:lenv(dirnambas)),
     $              infonam(1:lenv(infonam))
C Read eigenvalues
               CALL openfile(numfil,fname)
               DO jr=1,jpr1size
                  READ(UNIT=numfil,FMT='(E12.6E2)',ERR=101) kvalp(jr)
               ENDDO
               CLOSE (UNIT=numfil)
            ENDIF
C
C -2.- Read the 'vctp' file (eigenvectors)
C ----------------------------------------
C
            IF (present(kvctp)) THEN
C Build file name
               numjr=3
               CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
               WRITE(fname,'("./",A,"/",A)') 
     $              dirnambas(1:lenv(dirnambas)),
     $              infonam(1:lenv(infonam))
C Read eigenvectors
               CALL openfile(numfil,fname)
               DO jr=1,jpr2size
               DO jr1=1,jpr1size
                  READ(UNIT=numfil,FMT='(E12.6E2)',ERR=101) 
     $                 kvctp(jr1,jr)
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
C
      IF (nprint.GE.2) THEN
         kform='(8x,a,i3)'
         WRITE(numout,kform) ' - Covariance type: ',kvalbase
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','readinfobas')
C
 101  WRITE (texterror,*) 'error reading ASCII file, iost=',iost
      CALL printerror2(0,101,3,'hiobas','readinfobas',
     $     comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeinfobas(dirnambas,kvalbase,kvalp,kvctp)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write information files to covariance directory
CCC  -------
CCC  Method : Open information files and write covariance type,
CCC  ------   eigenvalues and eigenvectors
CCC
CCC  Input :  dirnambas  : name of the directory
CCC  -----    kvalbase   : covariance type
CCC           kvalp      : eigenvalues
CCC           kvctp      : eigenvectors
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      INTEGER, intent(in) :: kvalbase
      BIGREAL, dimension(:), optional, intent(in) :: kvalp
      BIGREAL, dimension(:,:), optional, intent(in) :: kvctp
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpr1size,jpr2size,jr,jr1,serie,numjr,jprbas
      CHARACTER(len=bgword) :: infonam,fname
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) '*** ROUTINE : ./writeinfobas :'
         WRITE(numout,*) '         write Cx, Cy or Co information files'
         WRITE(numout,*) '    ==> WRITING in directory ',
     $                                dirnambas(1:lenv(dirnambas))
      ENDIF
C
C Check input arguments
      IF (present(kvctp)) THEN
         jpr1size = size(kvctp,1)
         jpr2size = size(kvctp,2)
         IF ((present(kvalp)).AND.(jpr2size.NE.size(kvalp,1))) GOTO 1000
      ENDIF
C
C -1.- Write the 'valbase' file (covariance type)
C -----------------------------------------------
C Build file name
      serie=0
      numjr=1
      CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
      WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $                        infonam(1:lenv(infonam))
C Write informations
      CALL openfile(numfil,fname,kstatus=clunk)
      WRITE(UNIT=numfil,FMT='(i3.3)') kvalbase
      CLOSE (UNIT=numfil)
C
C -2.- Write the 'valp' file (eigenvalues)
C ----------------------------------------
C
      IF (present(kvalp)) THEN
C Build file name
         numjr=2
         CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $        infonam(1:lenv(infonam))
C Write eigenvalues
         CALL openfile(numfil,fname,kstatus=clunk)
         DO jr=1,jpr2size
            WRITE(UNIT=numfil,FMT='(E12.6E2)') kvalp(jr)
         ENDDO
         CLOSE (UNIT=numfil)
      ENDIF
C
C -3.- Write the 'vctp' file (eigenvectors)
C -----------------------------------------
C
      IF (present(kvctp)) THEN
C Build file name
         numjr=3
         CALL fildirbas (infonam,dirnambas,jprbas,numjr,serie)
         WRITE(fname,'("./",A,"/",A)') dirnambas(1:lenv(dirnambas)),
     $        infonam(1:lenv(infonam))
C Write eigenvectors
         CALL openfile(numfil,fname,kstatus=clunk)
         DO jr=1,jpr2size
         DO jr1=1,jpr1size
            WRITE(UNIT=numfil,FMT='(E12.6E2)') kvctp(jr1,jr)
         ENDDO
         ENDDO
         CLOSE (UNIT=numfil)
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','writeinfobas')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readvectb(dirnambas,kvectb)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read constraint b vector
CCC  --------
CCC  Method : Build b vector file name and
CCC  ------   read b vector
CCC
CCC  Input :  dirnambas  : name of the directory
CCC  -----
CCC  Output : kvectb     : b vector
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(out) :: kvectb
      CHARACTER(len=*), intent(in) :: dirnambas
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpmsize, jm
      CHARACTER(len=bgword) :: vectbname
CC----------------------------------------------------------------------
C
      jpmsize=size(kvectb,1)
C
      vectbname='vectb'
      vectbname=dirnambas(1:lenv(dirnambas)) // '/'
     $                  //  vectbname(1:lenv(vectbname))
C
      OPEN (unit=10,file=vectbname)
      DO jm=1,jpmsize
        READ(10,*) kvectb(jm)
      ENDDO
      CLOSE(10)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','readvectb')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readmatzmb(dirnambas,kmatzmb)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read local constraint b vectors
CCC  --------
CCC  Method : Build b vector file name and
CCC  ------   read b vector
CCC
CCC  Input :  dirnambas  : name of the directory
CCC  -----
CCC  Output : kmatzmb    : local b vectors
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo, only : texterror
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:,:), intent(out) :: kmatzmb
      CHARACTER(len=*), intent(in) :: dirnambas
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpmsize, jpzsize, jm, jz
      CHARACTER(len=bgword) :: vectbname
CC----------------------------------------------------------------------
C
      jpzsize=size(kmatzmb,1)
      jpmsize=size(kmatzmb,2)
C
      vectbname='matzmb'
      vectbname=dirnambas(1:lenv(dirnambas)) // '/'
     $                  //  vectbname(1:lenv(vectbname))
C
      OPEN (unit=10,file=vectbname)
      DO jz=1,jpzsize
        READ(10,*) ( kmatzmb(jz,jm), jm=1,jpmsize )
      ENDDO
      CLOSE(10)
C
      RETURN
C
C --- error management
C
 101  WRITE (texterror,*) 'Error reading file: ',vectbname
      CALL printerror2(0,102,3,'hiobas','readmatzmb',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readscalbas(dirnambas,filename,scal)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read scalar for every ensemble member
CCC  --------
CCC  Method : Read from file in ensemble directory
CCC  ------
CCC
CCC  Input :  dirnambas  : name of the directory
CCC  -----    filename   : name of the file
CCC
CCC  Output : scal       : scalar array
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(out) :: scal
      CHARACTER(len=*), intent(in) :: dirnambas
      CHARACTER(len=*), intent(in) :: filename
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jppsize, jp
      CHARACTER(len=bgword) :: fname
CC----------------------------------------------------------------------
C
      jppsize=size(scal,1)
C
      fname=dirnambas(1:lenv(dirnambas)) // '/'
     $                  //  filename(1:lenv(filename))
C
      OPEN (unit=10,file=fname)
      DO jp=1,jppsize
        READ(10,*) scal(jp)
      ENDDO
      CLOSE(10)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiobas','readscalbas')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiobas
