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
C ---                  HIOZON.F                                 ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-11  (C.E. Testut)                       ---
C --- modification : 99-11  (J.M. Brankart)                     ---
C --- modification : 00-03  (C.E. Testut)                       ---
C --- modification : 03-03  (J.M. Brankart)                     ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE evalhdrzon   : Read Vz file header
C --- SUBROUTINE readptzon    : Read Vz pointers
C --- SUBROUTINE readnbubzon  : Read list of Vz local data sections
C
C --- SUBROUTINE writehdrzon  : Write Vz file header
C --- SUBROUTINE writenbubzon : Write list of Vz local data sections
C --- SUBROUTINE writeptzon   : Write Vz pointers
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiozon
      use mod_main
      use lioczon
      use utilvalid
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrzon,readptzon,readnbubzon
      PUBLIC writehdrzon,writenbubzon,writeptzon

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrzon (kfninzon,kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read Vz file header
CCC  -------
CCC  Method : Call appropriate routine to read the right 'zon' format
CCC  ------
CCC  Input :   kfninzon : filename
CCC  -----
CCC  Output :  kzon_jpi : size of local data sections (along X dimension)
CCC  ------    kzon_jpj : size of local data sections (along Y dimension)
CCC            kzon_jpk : size of local data sections (along Z dimension)
CCC            kzon_jpt : size of local data sections (along T dimension)
CCC            kjpbub   : number of different local data sections
CCC            kjpz     : number of local data sections
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
       CHARACTER(len=*), intent(in) :: kfninzon
       INTEGER, intent(out) :: kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jextzon
CC----------------------------------------------------------------------
C
C Check validity of input file
      IF (.NOT.(validextzon(kfninzon))) GOTO 101
      jextzon=indext(kfninzon,extzontab,nbextzon)
C
C Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
C --- czon file format
            CALL evalhdrczon (kfninzon,kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiozon','evalhdrzon')
 1001 CALL printerror2(0,1001,3,'hiozon','evalhdrzon')
C
 101  WRITE (texterror,*) 'Invalid file extension: ',
     $          kfninzon(1:lenv(kfninzon))
      CALL printerror2(0,101,3,'hiozon','evalhdrzon',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readptzon (kfninzon,ptbubidx,ptdtalon,ptdtalat,
     $           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth,
     $           ptbubtime)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read Vz pointers
CCC  -------
CCC  Method : Call appropriate routine to read the right 'zon' format
CCC  ------
CCC  Input :  kfninzon   : filename
CCC  -----
CCC  Output : ptbubidx   :
CCC  ------   ptdtalon   :
CCC           ptdtalat   :
CCC           ptdtadepth :  Vz pointers, see description
CCC           ptdtatime  :  in 'mod_spacexyo.F'
CCC           ptbublon   :
CCC           ptbublat   :
CCC           ptbubdepth :
CCC           ptbubtime  :
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninzon
      INTEGER, dimension (:,:), intent(out) :: ptbubidx
      INTEGER, dimension (:,:), intent(out) :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), intent(out) :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), intent(out) :: ptbublon, ptbublat
      INTEGER, dimension (:,:), intent(out) :: ptbubdepth, ptbubtime
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jextzon
CC----------------------------------------------------------------------
C
C Check validity of input directory
      IF (.NOT.(validextzon(kfninzon))) GOTO 101
      jextzon = indext(kfninzon,extzontab,nbextzon)
C
C Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
C --- czon file format
         CALL readptczon(kfninzon,ptbubidx,ptdtalon,ptdtalat,ptdtadepth,
     $           ptdtatime,ptbublon,ptbublat,ptbubdepth,ptbubtime)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiozon','readptzon')
 1001 CALL printerror2(0,1001,3,'hiozon','readptzon')
C
 101  WRITE (texterror,*) 'Invalid file extension: ',
     $          kfninzon(1:lenv(kfninzon))
      CALL printerror2(0,101,3,'hiozon','readptzon',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readnbubzon(kfninzon,kptbubidx,kbub,klectinfo) 
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read list of Vz local data sections
CCC  -------
CCC  Method : Call appropriate routine to read the right 'zon' format
CCC  ------
CCC  Input :  kfninzon  : filename
CCC  -----    kptbubidx : list of local data sections to read
CCC           klectinfo : check or not zone file header
CCC
CCC  Output : kbub : list of local data sections
CCC  ------          (zon_jpi*zon_jpi*zon_jpk*zon_jpt*nbub)
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninzon
      INTEGER, dimension(:), intent(in) :: kptbubidx
      BIGREAL, dimension(:,:,:,:,:), intent(out) :: kbub
      LOGICAL, intent(in) :: klectinfo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jextzon
CC----------------------------------------------------------------------
C
      IF (klectinfo) THEN
C Check size of input arrays
         IF (size(kbub,1).NE.zon_jpi) GOTO 1000
         IF (size(kbub,2).NE.zon_jpj) GOTO 1000
         IF (size(kbub,3).NE.zon_jpk) GOTO 1000
         IF (size(kbub,4).NE.zon_jpt) GOTO 1000
         IF (size(kptbubidx,1).NE.size(kbub,5)) GOTO 1000
C Check validity of input file
         IF (.NOT.(validextzon(kfninzon))) GOTO 101
      ENDIF
      jextzon=indext(kfninzon,extzontab,nbextzon)
C
C Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
C --- czon file format
         CALL readnbubczon(kfninzon,kptbubidx(:),
     $        kbub(:,:,:,:,:),klectinfo)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiozon','readnbubzon')
C
 101  WRITE (texterror,*) 'Invalid file extension: ',
     $          kfninzon(1:lenv(kfninzon))
      CALL printerror2(0,101,3,'hiozon','readnbubzon',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrzon (kfnoutzon,kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write Vz file header
CCC  -------
CCC  Method : Call appropriate routine to write the right 'zon' format
CCC  ------
CCC  Input : kfnoutzon : filename
CCC  -----   kzon_jpi  : size of local data sections (along X dimension)
CCC          kzon_jpj  : size of local data sections (along Y dimension)
CCC          kzon_jpk  : size of local data sections (along Z dimension)
CCC          kzon_jpt  : size of local data sections (along T dimension)
CCC          kjpbub    : number of different local data sections
CCC          kjpz      : number of local data sections
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutzon
      INTEGER, intent(in) :: kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jextzon
CC----------------------------------------------------------------------
C
C Check validity of output file
      IF (.NOT.(validextzon(kfnoutzon))) GOTO 101
      jextzon = indext(kfnoutzon,extzontab,nbextzon)
C
C Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
C --- czon file format
         CALL writehdrczon(kfnoutzon,kzon_jpi,kzon_jpj,kzon_jpk,
     $             kzon_jpt,kjpbub,kjpz)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiozon','writehdrzon')
C
 101  WRITE (texterror,*) 'Invalid file extension: ',
     $          kfnoutzon(1:lenv(kfnoutzon))
      CALL printerror2(0,101,3,'hiozon','writehdrzon',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writenbubzon (kfnoutzon,kbub,kjbub)
CCC---------------------------------------------------------------------
CCC   
CCC  Purpose : Write list of Vz local data sections
CCC  -------
CCC  Method : Call appropriate routine to write the right 'zon' format
CCC  ------
CCC  Input : kfnoutzon : filename
CCC  -----   kjbub     : index of first local data section to write
CCC          kbub      : array of local data sections to write
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ==================
      CHARACTER(len=*), intent(in) :: kfnoutzon
      BIGREAL, dimension(:,:,:,:,:), intent(in) :: kbub
      INTEGER, intent(in) :: kjbub
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jextzon
CC----------------------------------------------------------------------
C
C Check validity of output file
      IF (.NOT.(validextzon(kfnoutzon))) GOTO 101
      jextzon = indext(kfnoutzon,extzontab,nbextzon)
C
C Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
C --- czon file format
         CALL writenbubczon(kfnoutzon,kbub,kjbub)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiozon','writenbubzon')
C
 101  WRITE (texterror,*) 'Invalid file extension: ',
     $          kfnoutzon(1:lenv(kfnoutzon))
      CALL printerror2(0,101,3,'hiozon','writenbubzon',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeptzon (kfnoutzon,ptbubidx,ptdtalon,ptdtalat,
     $           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth,
     $           ptbubtime)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write Vz pointers
CCC  -------
CCC  Method : Call appropriate routine to write the right 'zon' format
CCC  ------
CCC  Input : kfnoutzon   : filename
CCC  -----   ptbubidx   :
CCC          ptdtalon   :
CCC          ptdtalat   :
CCC          ptdtadepth :  Vz pointers, see description
CCC          ptdtatime  :  in 'mod_spacexyo.F'
CCC          ptbublon   :
CCC          ptbublat   :
CCC          ptbubdepth :
CCC          ptbubtime  :
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutzon
      INTEGER, dimension (:,:), intent(in) :: ptbubidx
      INTEGER, dimension (:,:), intent(in) :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), intent(in) :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), intent(in) :: ptbublon, ptbublat
      INTEGER, dimension (:,:), intent(in) :: ptbubdepth, ptbubtime
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jextzon
CC----------------------------------------------------------------------
C
C Check validity of input directory
      IF (.NOT.(validextzon(kfnoutzon))) GOTO 101
      jextzon = indext(kfnoutzon,extzontab,nbextzon)
C
C Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
C --- czon file format
         CALL writeptczon(kfnoutzon,ptbubidx,ptdtalon,ptdtalat,
     $           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth,
     $           ptbubtime)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'hiozon','writeptzon')
 1001 CALL printerror2(0,1001,3,'hiozon','writeptzon')
C
 101  WRITE (texterror,*) 'Invalid file extension: ',
     $          kfnoutzon(1:lenv(kfnoutzon))
      CALL printerror2(0,101,3,'hiozon','writeptzon',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiozon
