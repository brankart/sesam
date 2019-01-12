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
! ---                  HIOZON.F90                                 ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-11  (C.E. Testut)                       ---
! --- modification : 99-11  (J.M. Brankart)                     ---
! --- modification : 00-03  (C.E. Testut)                       ---
! --- modification : 03-03  (J.M. Brankart)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE evalhdrzon   : Read Vz file header
! --- SUBROUTINE readptzon    : Read Vz pointers
! --- SUBROUTINE readnbubzon  : Read list of Vz local data sections
!
! --- SUBROUTINE writehdrzon  : Write Vz file header
! --- SUBROUTINE writenbubzon : Write list of Vz local data sections
! --- SUBROUTINE writeptzon   : Write Vz pointers
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiozon
      use mod_main
      use lioczon
      use utilvalid
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrzon,readptzon,readnbubzon
      PUBLIC writehdrzon,writenbubzon,writeptzon

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrzon (kfninzon,kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz)
!---------------------------------------------------------------------
!
!  Purpose : Read Vz file header
!  -------
!  Method : Call appropriate routine to read the right 'zon' format
!  ------
!  Input :   kfninzon : filename
!  -----
!  Output :  kzon_jpi : size of local data sections (along X dimension)
!  ------    kzon_jpj : size of local data sections (along Y dimension)
!            kzon_jpk : size of local data sections (along Z dimension)
!            kzon_jpt : size of local data sections (along T dimension)
!            kjpbub   : number of different local data sections
!            kjpz     : number of local data sections
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
       CHARACTER(len=*), intent(in) :: kfninzon
       INTEGER, intent(out) :: kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextzon
!----------------------------------------------------------------------
!
! Check validity of input file
      IF (.NOT.(validextzon(kfninzon))) GOTO 101
      jextzon=indext(kfninzon,extzontab,nbextzon)
!
! Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
! --- czon file format
            CALL evalhdrczon (kfninzon,kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiozon','evalhdrzon')
 1001 CALL printerror2(0,1001,3,'hiozon','evalhdrzon')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &          kfninzon(1:lenv(kfninzon))
      CALL printerror2(0,101,3,'hiozon','evalhdrzon',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readptzon (kfninzon,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime)
!---------------------------------------------------------------------
!
!  Purpose : Read Vz pointers
!  -------
!  Method : Call appropriate routine to read the right 'zon' format
!  ------
!  Input :  kfninzon   : filename
!  -----
!  Output : ptbubidx   :
!  ------   ptdtalon   :
!           ptdtalat   :
!           ptdtadepth :  Vz pointers, see description
!           ptdtatime  :  in 'mod_spacexyo.F90'
!           ptbublon   :
!           ptbublat   :
!           ptbubdepth :
!           ptbubtime  :
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninzon
      INTEGER, dimension (:,:), intent(out) :: ptbubidx
      INTEGER, dimension (:,:), intent(out) :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), intent(out) :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), intent(out) :: ptbublon, ptbublat
      INTEGER, dimension (:,:), intent(out) :: ptbubdepth, ptbubtime
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextzon
!----------------------------------------------------------------------
!
! Check validity of input directory
      IF (.NOT.(validextzon(kfninzon))) GOTO 101
      jextzon = indext(kfninzon,extzontab,nbextzon)
!
! Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
! --- czon file format
         CALL readptczon(kfninzon,ptbubidx,ptdtalon,ptdtalat,ptdtadepth, &
     &           ptdtatime,ptbublon,ptbublat,ptbubdepth,ptbubtime)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiozon','readptzon')
 1001 CALL printerror2(0,1001,3,'hiozon','readptzon')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &          kfninzon(1:lenv(kfninzon))
      CALL printerror2(0,101,3,'hiozon','readptzon',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readnbubzon(kfninzon,kptbubidx,kbub,klectinfo) 
!---------------------------------------------------------------------
!
!  Purpose : Read list of Vz local data sections
!  -------
!  Method : Call appropriate routine to read the right 'zon' format
!  ------
!  Input :  kfninzon  : filename
!  -----    kptbubidx : list of local data sections to read
!           klectinfo : check or not zone file header
!
!  Output : kbub : list of local data sections
!  ------          (zon_jpi*zon_jpi*zon_jpk*zon_jpt*nbub)
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninzon
      INTEGER, dimension(:), intent(in) :: kptbubidx
      BIGREAL, dimension(:,:,:,:,:), intent(out) :: kbub
      LOGICAL, intent(in) :: klectinfo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextzon
!----------------------------------------------------------------------
!
      IF (klectinfo) THEN
! Check size of input arrays
         IF (size(kbub,1).NE.zon_jpi) GOTO 1000
         IF (size(kbub,2).NE.zon_jpj) GOTO 1000
         IF (size(kbub,3).NE.zon_jpk) GOTO 1000
         IF (size(kbub,4).NE.zon_jpt) GOTO 1000
         IF (size(kptbubidx,1).NE.size(kbub,5)) GOTO 1000
! Check validity of input file
         IF (.NOT.(validextzon(kfninzon))) GOTO 101
      ENDIF
      jextzon=indext(kfninzon,extzontab,nbextzon)
!
! Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
! --- czon file format
         CALL readnbubczon(kfninzon,kptbubidx(:), &
     &        kbub(:,:,:,:,:),klectinfo)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiozon','readnbubzon')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &          kfninzon(1:lenv(kfninzon))
      CALL printerror2(0,101,3,'hiozon','readnbubzon',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrzon (kfnoutzon,kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz)
!---------------------------------------------------------------------
!
!  Purpose : Write Vz file header
!  -------
!  Method : Call appropriate routine to write the right 'zon' format
!  ------
!  Input : kfnoutzon : filename
!  -----   kzon_jpi  : size of local data sections (along X dimension)
!          kzon_jpj  : size of local data sections (along Y dimension)
!          kzon_jpk  : size of local data sections (along Z dimension)
!          kzon_jpt  : size of local data sections (along T dimension)
!          kjpbub    : number of different local data sections
!          kjpz      : number of local data sections
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutzon
      INTEGER, intent(in) :: kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextzon
!----------------------------------------------------------------------
!
! Check validity of output file
      IF (.NOT.(validextzon(kfnoutzon))) GOTO 101
      jextzon = indext(kfnoutzon,extzontab,nbextzon)
!
! Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
! --- czon file format
         CALL writehdrczon(kfnoutzon,kzon_jpi,kzon_jpj,kzon_jpk, &
     &             kzon_jpt,kjpbub,kjpz)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiozon','writehdrzon')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &          kfnoutzon(1:lenv(kfnoutzon))
      CALL printerror2(0,101,3,'hiozon','writehdrzon',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writenbubzon (kfnoutzon,kbub,kjbub)
!---------------------------------------------------------------------
!   
!  Purpose : Write list of Vz local data sections
!  -------
!  Method : Call appropriate routine to write the right 'zon' format
!  ------
!  Input : kfnoutzon : filename
!  -----   kjbub     : index of first local data section to write
!          kbub      : array of local data sections to write
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      CHARACTER(len=*), intent(in) :: kfnoutzon
      BIGREAL, dimension(:,:,:,:,:), intent(in) :: kbub
      INTEGER, intent(in) :: kjbub
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextzon
!----------------------------------------------------------------------
!
! Check validity of output file
      IF (.NOT.(validextzon(kfnoutzon))) GOTO 101
      jextzon = indext(kfnoutzon,extzontab,nbextzon)
!
! Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
! --- czon file format
         CALL writenbubczon(kfnoutzon,kbub,kjbub)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiozon','writenbubzon')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &          kfnoutzon(1:lenv(kfnoutzon))
      CALL printerror2(0,101,3,'hiozon','writenbubzon',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeptzon (kfnoutzon,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime)
!---------------------------------------------------------------------
!
!  Purpose : Write Vz pointers
!  -------
!  Method : Call appropriate routine to write the right 'zon' format
!  ------
!  Input : kfnoutzon   : filename
!  -----   ptbubidx   :
!          ptdtalon   :
!          ptdtalat   :
!          ptdtadepth :  Vz pointers, see description
!          ptdtatime  :  in 'mod_spacexyo.F90'
!          ptbublon   :
!          ptbublat   :
!          ptbubdepth :
!          ptbubtime  :
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutzon
      INTEGER, dimension (:,:), intent(in) :: ptbubidx
      INTEGER, dimension (:,:), intent(in) :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), intent(in) :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), intent(in) :: ptbublon, ptbublat
      INTEGER, dimension (:,:), intent(in) :: ptbubdepth, ptbubtime
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jextzon
!----------------------------------------------------------------------
!
! Check validity of input directory
      IF (.NOT.(validextzon(kfnoutzon))) GOTO 101
      jextzon = indext(kfnoutzon,extzontab,nbextzon)
!
! Select the right 'zon' file format
      SELECT CASE (jextzon)
      CASE (1)
! --- czon file format
         CALL writeptczon(kfnoutzon,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiozon','writeptzon')
 1001 CALL printerror2(0,1001,3,'hiozon','writeptzon')
!
 101  WRITE (texterror,*) 'Invalid file extension: ', &
     &          kfnoutzon(1:lenv(kfnoutzon))
      CALL printerror2(0,101,3,'hiozon','writeptzon',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiozon
