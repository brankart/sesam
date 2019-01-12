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
! ---                  LIOCZON.F90                                ---
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
! --- SUBROUTINE evalhdrczon   : Read 'czon' file header
! --- SUBROUTINE readptczon    : Read Vz pointers from 'czon' file
! --- SUBROUTINE readnbubczon  : Read list of local data sections from 'czon' file
!
! --- SUBROUTINE writehdrczon  : Write 'czon' file header
! --- SUBROUTINE writenbubczon : Write list of local data sections in 'czon' file
! --- SUBROUTINE writeptczon   : Write Vz pointers in 'czon' file
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE lioczon
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrczon,readptczon,readnbubczon
      PUBLIC writehdrczon,writenbubczon,writeptczon

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrczon (kfninzon,kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz)
!---------------------------------------------------------------------
!
!  Purpose : Read 'czon' file dimensions
!  -------
!  Method :  Call low level routine defining czon format (NetCDF)
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
      use mod_main
      use mod_cfgxyo
      use utilcdfzon
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
      CHARACTER(len=word80) :: cdrec1
      INTEGER :: dtaend1,varlg1
!----------------------------------------------------------------------
!
! Read 'czon' file dimensions
      CALL cdfrdimzon(kfninzon,kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz,dtaend1,varlg1,cdrec1)
!
      RETURN
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readptczon (kfninzon,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime)
!---------------------------------------------------------------------
!
!  Purpose : Read Vz pointers from 'czon' file
!  -------
!  Method :  Call low level routine defining czon format (NetCDF)
!  ------
!  Input :   kfninzon   : filename
!  -----
!  Output :  ptbubidx   :
!  ------    ptdtalon   :
!            ptdtalat   :
!            ptdtadepth :  Vz pointers, see description
!            ptdtatime  :  in 'mod_spacexyo.F90'
!            ptbublon   :
!            ptbublat   :
!            ptbubdepth :
!            ptbubtime  :
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilcdfzon
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
      INTEGER :: jzstart,dtastart,jzend
!----------------------------------------------------------------------
!
! Check size of input arrays
!
      dtastart = 1
!
      IF (size(ptbubidx,2).NE.dtaend) GOTO 1000
      IF (size(ptbublon,2).NE.dtaend) GOTO 1000
      IF (size(ptbublat,2).NE.dtaend) GOTO 1000
      IF (size(ptbubtime,2).NE.dtaend) GOTO 1000
      IF (size(ptbubdepth,2).NE.dtaend) GOTO 1000
      IF (size(ptdtalon,2).NE.dtaend) GOTO 1000
      IF (size(ptdtalat,2).NE.dtaend) GOTO 1000
      IF (size(ptdtatime,2).NE.dtaend) GOTO 1000
      IF (size(ptdtadepth,2).NE.dtaend) GOTO 1000
!
      jzstart = 1
      jzend = size(ptbubidx,1)
!
      IF (size(ptbublon,1).NE.jzend) GOTO 1000
      IF (size(ptbublat,1).NE.jzend) GOTO 1000
      IF (size(ptbubtime,1).NE.jzend) GOTO 1000
      IF (size(ptbubdepth,1).NE.jzend) GOTO 1000
      IF (size(ptdtalon,1).NE.jzend) GOTO 1000
      IF (size(ptdtalat,1).NE.jzend) GOTO 1000
      IF (size(ptdtatime,1).NE.jzend) GOTO 1000
      IF (size(ptdtadepth,1).NE.jzend) GOTO 1000
!
! Read pointers from .czon file
!
      CALL cdfrptzon(kfninzon,ptbubidx,ptdtalon, &
     &       ptdtalat,ptdtadepth,ptdtatime, &
     &       ptbublon,ptbublat,ptbubdepth, &
     &       ptbubtime,jzstart,jzend,dtastart,dtaend)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioczon','readptczon')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readnbubczon(kfninzon,kptbubidx,kbub,klectinfo)
!---------------------------------------------------------------------
!
!  Purpose : Read list of local data sections from 'czon' file
!  -------
!  Method :  Call low level routine defining czon format (NetCDF)
!  ------    to read selected local data sections
!
!  Input :   kfninzon  : filename
!  -----     kptbubidx : index of local data sections to read
!            klectinfo : check or not 'czon' file header
!
!  Output :  kbub : list of local data sections
!  ------           (zon_jpi*zon_jpi*zon_jpk*zon_jpt*nbub)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpz,jpyend
      use utilcdfzon
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
      CHARACTER(len=word80) :: cdrec1
      INTEGER :: dtaend1,jpyend1
      INTEGER :: varlg1
      INTEGER, dimension(1:nbvar) :: dta_ord1, &
     &                 dta_dim1,dta_nbr1,dta_ind1
      CHARACTER(len=varlg), dimension(1:nbvar) :: dta_nam1
      BIGREAL4, dimension(1:nbvar) :: dta_moy1,dta_ect1
      INTEGER :: allocok,jzstart,dtastart,jz,jbubmin,jbubmax, &
     &     jbubdeb,jbubfin,jbub,jpnbubsize,bubcount
      LOGICAL :: incompatible,sameend,sameord,affect
      INTEGER :: indvar,jdta,jdta1,inddta
      INTEGER :: kzon_jpi,kzon_jpj,kzon_jpk,kzon_jpt,kjpbub,kjpz
      LOGICAL, dimension(1:nbvar) :: transfert
      BIGREAL4, allocatable, dimension(:,:,:,:,:) :: bub
!----------------------------------------------------------------------
      jpnbubsize=size(kptbubidx,1)
!
! -1.- Check 'czon' file dimensions
! ---------------------------------
!
      CALL cdfrdimzon(kfninzon,kzon_jpi,kzon_jpj,kzon_jpk, &
     &     kzon_jpt,kjpbub,kjpz,dtaend1,varlg1,cdrec1)
!
      IF (klectinfo) THEN
         IF (kzon_jpi.NE.zon_jpi) GOTO 1000
         IF (kzon_jpj.NE.zon_jpj) GOTO 1000
         IF (kzon_jpk.NE.zon_jpk) GOTO 1000
         IF (kzon_jpt.NE.zon_jpt) GOTO 1000
         IF (kjpz.NE.jpz) GOTO 1000
      ENDIF
!
! -2.- Read 'czon' file header
! ----------------------------
!
      CALL cdfrhdrzon(kfninzon,dta_nam1,dta_dim1, &
     &     dta_nbr1,dta_moy1,dta_ect1)
!
! -3.- Check 'czon' file compatibilty if required
! -----------------------------------------------
!
      IF (klectinfo) THEN 
!
! Compute size of Vy object from 'czon' file informations
         dta_ind1(1)=1
         DO jdta1=2,dtaend1
            dta_ind1(jdta1)=dta_ind1(jdta1-1)+dta_nbr1(jdta1-1)
         ENDDO
         jpyend1=SUM(dta_nbr1(1:dtaend1))
!
! Compute order of variables fields in 'czon' file
         DO jdta1=1,dtaend1
            affect = .FALSE.
            LOOP1 : DO inddta=1,nbvar
               IF (dta_nam1(jdta1).EQ.dta_nam(inddta)) THEN
                  dta_ord1(jdta1)=inddta
                  affect = .TRUE.
                  EXIT LOOP1
               ENDIF
            ENDDO LOOP1  
            IF (.NOT.(affect)) GOTO 103
         ENDDO
         dta_ord1((dtaend1+1):nbvar) = 0
!
! Check if 'czon' file configuration is identical to SESAM configuration
         sameend = .TRUE.
         sameend = ((dtaend1.EQ.dtaend).AND.sameend)
         sameend = ((jpyend1.EQ.jpyend).AND.sameend)
         sameend = ((varlg1.LE.varlg).AND.sameend)
         sameord = .TRUE.
         DO jdta1=1,dtaend1
            sameord = ((dta_ord1(jdta1).EQ.dta_ord(jdta1)).AND.sameord)
         ENDDO
!
         incompatible = (.NOT.(sameord.AND.sameend))
         DO jdta=1,dtaend
            inddta=dta_ord(jdta)
            incompatible = ((dta_nam(inddta).NE.dta_nam1(jdta)) &
     &           .OR.incompatible)
            incompatible = ((dta_dim(inddta).NE.dta_dim1(jdta)) &
     &           .OR.incompatible)
            incompatible = ((dta_nbr(inddta).NE.dta_nbr1(jdta)) &
     &           .OR.incompatible)
         ENDDO
!
         IF (incompatible) GOTO 102
      ENDIF
!
! -4.- Take decision for centering/reducing input variables
! ---------------------------------------------------------
!
      transfert(:)=.FALSE.
      IF (ANY(dta_ect1(1:dtaend1).NE.0.0)) THEN
         DO jdta=1,dtaend1
            inddta=dta_ord(jdta)
            transfert(jdta) = ( &
     &        ((FREAL4(dta_moy(inddta)).NE.dta_moy1(jdta)) &
     &        .AND.(lmoyect)) &
     &        .OR. ((dta_moy1(jdta).NE.(FREAL4(0.0))) &
     &        .AND.(.NOT.lmoyect)) &
     &        .OR. transfert(jdta) )
            transfert(jdta) = ( &
     &           ((FREAL4(dta_ect(inddta)).NE.dta_ect1(jdta)) &
     &           .AND.(lmoyect)) &
     &           .OR. ((dta_ect1(jdta).NE.(FREAL4(1.0))) &
     &           .AND.(.NOT.lmoyect)) &
     &           .OR. transfert(jdta) )
            IF (transfert(jdta)) THEN
               GOTO 104
            ENDIF
         ENDDO

      ENDIF
!
! -6.- Read local data sections from 'czon' file
! ----------------------------------------------
!
! Allocate current local data section array
      allocate ( bub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &     1:jpnbubsize),stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      bub(:,:,:,:,:) = FREAL4 (0.0)
!
! Load list of local data sections from file
      jbubmin=MINVAL(kptbubidx(:))
      jbubmax=MAXVAL(kptbubidx(:))
      IF (jbubmin.LE.0) GOTO 1000
      IF (jbubmax.GT.kjpbub) GOTO 1000
      DO jbubdeb=jbubmin,jbubmax,jpnbubsize
         jbubfin=MIN(jbubdeb-1+jpnbubsize,jbubmax)
!
         bubcount=jbubfin-(jbubdeb-1)
         CALL cdfrbubzon(kfninzon,bub,jbubdeb,bubcount)
      ENDDO
!
      kbub(:,:,:,:,:)=FREAL(bub(:,:,:,:,:))
!
      IF (allocated(bub)) deallocate(bub)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioczon','readnbubczon')
 1001 CALL printerror2(0,1001,3,'lioczon','readnbubczon')
!
 102  WRITE (texterror,*) 'Parameters in .czon file header', &
     &     ' incompatible with SESAM configuration'
      CALL printerror2(0,102,3,'lioczon','readnbubczon',comment=texterror)
 103  WRITE (texterror,*) 'Variable names in .czon file', &
     &     ' incompatible with SESAM configuration'
      CALL printerror2(0,103,3,'lioczon','readnbubczon',comment=texterror)
 104  WRITE (texterror,*) 'No centering/reducing', &
     &     ' is possible with readnbubzon'
      CALL printerror2(0,104,3,'lioczon','readnbubczon',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrczon (kfnoutzon,kzon_jpi,kzon_jpj,kzon_jpk, &
     &           kzon_jpt,kjpbub,kjpz)
!---------------------------------------------------------------------
!
!  Purpose : Write 'czon' file header
!  -------
!  Method :  Call low level routine defining czon format (NetCDF)
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
      use mod_main
      use mod_cfgxyo
      use utilcdfzon
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
      INTEGER :: jdta,inddta
      INTEGER, dimension(1:nbvar) :: dta_dim1,dta_nbr1
      CHARACTER(len=varlg), dimension(1:nbvar) :: dta_nam1
      BIGREAL4, dimension(1:nbvar) :: dta_moy1,dta_ect1
      CHARACTER(len=word80) :: title
!----------------------------------------------------------------------
!
! Set header variables
!
      title = 'SESAM Vz object'
!
      DO jdta = 1,dtaend
         inddta = dta_ord(jdta)
         dta_nam1(jdta) = dta_nam(inddta)
         dta_dim1(jdta) = dta_dim(inddta)
         dta_nbr1(jdta) = dta_nbr(inddta)
         IF (lmoyect) THEN
            dta_moy1(jdta) = FREAL4(dta_moy(inddta))
            dta_ect1(jdta) = FREAL4(dta_ect(inddta))
         ELSE
            dta_moy1(jdta) = FREAL4(0.0)
            dta_ect1(jdta) = FREAL4(1.0)
         ENDIF
      ENDDO
!
! Write 'czon' file header
!
      CALL cdfwdimzon(kfnoutzon,kzon_jpi,kzon_jpj,kzon_jpk,kzon_jpt, &
     &               kjpbub,kjpz,dtaend,varlg,title)
!
      CALL cdfwhdrzon(kfnoutzon,dta_nam1,dta_dim1,dta_nbr1, &
     &     dta_moy1,dta_ect1)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioczon','writehdrczon')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writenbubczon (kfnoutzon,kbub,kjbub)
!---------------------------------------------------------------------
!   
!  Purpose : Write list of local data sections in 'czon' file
!  -------
!  Method : Call low level routine defining czon format (NetCDF)
!  ------   to write selected local data sections
!
!  Input : kfnoutzon : filename
!  -----   kjbub     : index of first local data section to write
!          kbub      : array of local data sections to write
!
!---------------------------------------------------------------------
! modules
! ======= 
      use mod_main
      use mod_cfgxyo
      use utilcdfzon
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnoutzon
      BIGREAL, dimension(:,:,:,:,:), intent(in) :: kbub
      INTEGER, intent(in) :: kjbub
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL4, dimension(:,:,:,:,:), allocatable :: kbub4
      INTEGER :: jpisize,jpjsize,jpksize,jptsize, &
     &     jpzsize,jpnbubsize,bubcount
      INTEGER :: allocok
!----------------------------------------------------------------------
      jpisize = size(kbub,1)
      jpjsize = size(kbub,2)
      jpksize = size(kbub,3)
      jptsize = size(kbub,4)
      jpnbubsize=size(kbub,5)
! allocate local data section array to write (kr4 real kind)
      allocate ( kbub4(1:jpisize,1:jpjsize,1:jpksize, &
     &     1:jptsize,1:jpnbubsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
      kbub4(:,:,:,:,:) = FREAL4(kbub(:,:,:,:,:))
!
! Write block of local data sections in 'czon' file
!
      bubcount=jpnbubsize
      CALL cdfwbubzon(kfnoutzon,kbub4,kjbub,bubcount)
!
      IF (allocated(kbub4)) deallocate (kbub4)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioczon','writenbubczon')
 1001 CALL printerror2(0,1001,3,'lioczon','writenbubczon')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeptczon (kfnoutzon,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime)
!---------------------------------------------------------------------
!
!  Purpose : Write Vz pointers in 'czon' file
!  -------
!  Method : Call low level routine defining czon format (NetCDF)
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
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use utilcdfzon
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
      INTEGER :: dtastart,jzstart,jzend
!----------------------------------------------------------------------
!
! Check size of input arrays
!
      dtastart = 1
!
      IF (size(ptbubidx,2).NE.dtaend) GOTO 1000
      IF (size(ptbublon,2).NE.dtaend) GOTO 1000
      IF (size(ptbublat,2).NE.dtaend) GOTO 1000
      IF (size(ptbubtime,2).NE.dtaend) GOTO 1000
      IF (size(ptbubdepth,2).NE.dtaend) GOTO 1000
      IF (size(ptdtalon,2).NE.dtaend) GOTO 1000
      IF (size(ptdtalat,2).NE.dtaend) GOTO 1000
      IF (size(ptdtatime,2).NE.dtaend) GOTO 1000
      IF (size(ptdtadepth,2).NE.dtaend) GOTO 1000
!
      jzstart = 1
      jzend = size(ptbubidx,1)
!
      IF (size(ptbublon,1).NE.jzend) GOTO 1000
      IF (size(ptbublat,1).NE.jzend) GOTO 1000
      IF (size(ptbubtime,1).NE.jzend) GOTO 1000
      IF (size(ptbubdepth,1).NE.jzend) GOTO 1000
      IF (size(ptdtalon,1).NE.jzend) GOTO 1000
      IF (size(ptdtalat,1).NE.jzend) GOTO 1000
      IF (size(ptdtatime,1).NE.jzend) GOTO 1000
      IF (size(ptdtadepth,1).NE.jzend) GOTO 1000
!
! Write pointers in .czon file
!
      CALL cdfwptzon(kfnoutzon,ptbubidx,ptdtalon,ptdtalat, &
     &       ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &       ptbubtime,jzstart,jzend,dtastart,dtaend)

!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'lioczon','writeptczon')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE lioczon
