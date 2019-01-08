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
C ---                  LIOCZON.F                                ---
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
C --- SUBROUTINE evalhdrczon   : Read 'czon' file header
C --- SUBROUTINE readptczon    : Read Vz pointers from 'czon' file
C --- SUBROUTINE readnbubczon  : Read list of local data sections from 'czon' file
C
C --- SUBROUTINE writehdrczon  : Write 'czon' file header
C --- SUBROUTINE writenbubczon : Write list of local data sections in 'czon' file
C --- SUBROUTINE writeptczon   : Write Vz pointers in 'czon' file
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE lioczon
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrczon,readptczon,readnbubczon
      PUBLIC writehdrczon,writenbubczon,writeptczon

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrczon (kfninzon,kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read 'czon' file dimensions
CCC  -------
CCC  Method :  Call low level routine defining czon format (NetCDF)
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
      use mod_main
      use mod_cfgxyo
      use utilcdfzon
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
      CHARACTER(len=word80) :: cdrec1
      INTEGER :: dtaend1,varlg1
CC----------------------------------------------------------------------
C
C Read 'czon' file dimensions
      CALL cdfrdimzon(kfninzon,kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz,dtaend1,varlg1,cdrec1)
C
      RETURN
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readptczon (kfninzon,ptbubidx,ptdtalon,ptdtalat,
     $           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth,
     $           ptbubtime)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read Vz pointers from 'czon' file
CCC  -------
CCC  Method :  Call low level routine defining czon format (NetCDF)
CCC  ------
CCC  Input :   kfninzon   : filename
CCC  -----
CCC  Output :  ptbubidx   :
CCC  ------    ptdtalon   :
CCC            ptdtalat   :
CCC            ptdtadepth :  Vz pointers, see description
CCC            ptdtatime  :  in 'mod_spacexyo.F'
CCC            ptbublon   :
CCC            ptbublat   :
CCC            ptbubdepth :
CCC            ptbubtime  :
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilcdfzon
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
      INTEGER :: jzstart,dtastart,jzend
CC----------------------------------------------------------------------
C
C Check size of input arrays
C
      dtastart = 1
C
      IF (size(ptbubidx,2).NE.dtaend) GOTO 1000
      IF (size(ptbublon,2).NE.dtaend) GOTO 1000
      IF (size(ptbublat,2).NE.dtaend) GOTO 1000
      IF (size(ptbubtime,2).NE.dtaend) GOTO 1000
      IF (size(ptbubdepth,2).NE.dtaend) GOTO 1000
      IF (size(ptdtalon,2).NE.dtaend) GOTO 1000
      IF (size(ptdtalat,2).NE.dtaend) GOTO 1000
      IF (size(ptdtatime,2).NE.dtaend) GOTO 1000
      IF (size(ptdtadepth,2).NE.dtaend) GOTO 1000
C
      jzstart = 1
      jzend = size(ptbubidx,1)
C
      IF (size(ptbublon,1).NE.jzend) GOTO 1000
      IF (size(ptbublat,1).NE.jzend) GOTO 1000
      IF (size(ptbubtime,1).NE.jzend) GOTO 1000
      IF (size(ptbubdepth,1).NE.jzend) GOTO 1000
      IF (size(ptdtalon,1).NE.jzend) GOTO 1000
      IF (size(ptdtalat,1).NE.jzend) GOTO 1000
      IF (size(ptdtatime,1).NE.jzend) GOTO 1000
      IF (size(ptdtadepth,1).NE.jzend) GOTO 1000
C
C Read pointers from .czon file
C
      CALL cdfrptzon(kfninzon,ptbubidx,ptdtalon,
     $       ptdtalat,ptdtadepth,ptdtatime,
     $       ptbublon,ptbublat,ptbubdepth,
     $       ptbubtime,jzstart,jzend,dtastart,dtaend)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lioczon','readptczon')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readnbubczon(kfninzon,kptbubidx,kbub,klectinfo)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Read list of local data sections from 'czon' file
CCC  -------
CCC  Method :  Call low level routine defining czon format (NetCDF)
CCC  ------    to read selected local data sections
CCC
CCC  Input :   kfninzon  : filename
CCC  -----     kptbubidx : index of local data sections to read
CCC            klectinfo : check or not 'czon' file header
CCC
CCC  Output :  kbub : list of local data sections
CCC  ------           (zon_jpi*zon_jpi*zon_jpk*zon_jpt*nbub)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpz,jpyend
      use utilcdfzon
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
      CHARACTER(len=word80) :: cdrec1
      INTEGER :: dtaend1,jpyend1
      INTEGER :: varlg1
      INTEGER, dimension(1:nbvar) :: dta_ord1,
     $                 dta_dim1,dta_nbr1,dta_ind1
      CHARACTER(len=varlg), dimension(1:nbvar) :: dta_nam1
      BIGREAL4, dimension(1:nbvar) :: dta_moy1,dta_ect1
      INTEGER :: allocok,jzstart,dtastart,jz,jbubmin,jbubmax,
     $     jbubdeb,jbubfin,jbub,jpnbubsize,bubcount
      LOGICAL :: incompatible,sameend,sameord,affect
      INTEGER :: indvar,jdta,jdta1,inddta
      INTEGER :: kzon_jpi,kzon_jpj,kzon_jpk,kzon_jpt,kjpbub,kjpz
      LOGICAL, dimension(1:nbvar) :: transfert
      BIGREAL4, allocatable, dimension(:,:,:,:,:) :: bub
CC----------------------------------------------------------------------
      jpnbubsize=size(kptbubidx,1)
C
C -1.- Check 'czon' file dimensions
C ---------------------------------
C
      CALL cdfrdimzon(kfninzon,kzon_jpi,kzon_jpj,kzon_jpk,
     $     kzon_jpt,kjpbub,kjpz,dtaend1,varlg1,cdrec1)
C
      IF (klectinfo) THEN
         IF (kzon_jpi.NE.zon_jpi) GOTO 1000
         IF (kzon_jpj.NE.zon_jpj) GOTO 1000
         IF (kzon_jpk.NE.zon_jpk) GOTO 1000
         IF (kzon_jpt.NE.zon_jpt) GOTO 1000
         IF (kjpz.NE.jpz) GOTO 1000
      ENDIF
C
C -2.- Read 'czon' file header
C ----------------------------
C
      CALL cdfrhdrzon(kfninzon,dta_nam1,dta_dim1,
     $     dta_nbr1,dta_moy1,dta_ect1)
C
C -3.- Check 'czon' file compatibilty if required
C -----------------------------------------------
C
      IF (klectinfo) THEN 
C
C Compute size of Vy object from 'czon' file informations
         dta_ind1(1)=1
         DO jdta1=2,dtaend1
            dta_ind1(jdta1)=dta_ind1(jdta1-1)+dta_nbr1(jdta1-1)
         ENDDO
         jpyend1=SUM(dta_nbr1(1:dtaend1))
C
C Compute order of variables fields in 'czon' file
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
C
C Check if 'czon' file configuration is identical to SESAM configuration
         sameend = .TRUE.
         sameend = ((dtaend1.EQ.dtaend).AND.sameend)
         sameend = ((jpyend1.EQ.jpyend).AND.sameend)
         sameend = ((varlg1.LE.varlg).AND.sameend)
         sameord = .TRUE.
         DO jdta1=1,dtaend1
            sameord = ((dta_ord1(jdta1).EQ.dta_ord(jdta1)).AND.sameord)
         ENDDO
C
         incompatible = (.NOT.(sameord.AND.sameend))
         DO jdta=1,dtaend
            inddta=dta_ord(jdta)
            incompatible = ((dta_nam(inddta).NE.dta_nam1(jdta))
     $           .OR.incompatible)
            incompatible = ((dta_dim(inddta).NE.dta_dim1(jdta))
     $           .OR.incompatible)
            incompatible = ((dta_nbr(inddta).NE.dta_nbr1(jdta))
     $           .OR.incompatible)
         ENDDO
C
         IF (incompatible) GOTO 102
      ENDIF
C
C -4.- Take decision for centering/reducing input variables
C ---------------------------------------------------------
C
      transfert(:)=.FALSE.
      IF (ANY(dta_ect1(1:dtaend1).NE.0.0)) THEN
         DO jdta=1,dtaend1
            inddta=dta_ord(jdta)
            transfert(jdta) = (
     $        ((FREAL4(dta_moy(inddta)).NE.dta_moy1(jdta))
     $        .AND.(lmoyect))
     $        .OR. ((dta_moy1(jdta).NE.(FREAL4(0.0)))
     $        .AND.(.NOT.lmoyect))
     $        .OR. transfert(jdta) )
            transfert(jdta) = (
     $           ((FREAL4(dta_ect(inddta)).NE.dta_ect1(jdta))
     $           .AND.(lmoyect))
     $           .OR. ((dta_ect1(jdta).NE.(FREAL4(1.0)))
     $           .AND.(.NOT.lmoyect))
     $           .OR. transfert(jdta) )
            IF (transfert(jdta)) THEN
               GOTO 104
            ENDIF
         ENDDO

      ENDIF
C
C -6.- Read local data sections from 'czon' file
C ----------------------------------------------
C
C Allocate current local data section array
      allocate ( bub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,
     $     1:jpnbubsize),stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      bub(:,:,:,:,:) = FREAL4 (0.0)
C
C Load list of local data sections from file
      jbubmin=MINVAL(kptbubidx(:))
      jbubmax=MAXVAL(kptbubidx(:))
      IF (jbubmin.LE.0) GOTO 1000
      IF (jbubmax.GT.kjpbub) GOTO 1000
      DO jbubdeb=jbubmin,jbubmax,jpnbubsize
         jbubfin=MIN(jbubdeb-1+jpnbubsize,jbubmax)
C
         bubcount=jbubfin-(jbubdeb-1)
         CALL cdfrbubzon(kfninzon,bub,jbubdeb,bubcount)
      ENDDO
C
      kbub(:,:,:,:,:)=FREAL(bub(:,:,:,:,:))
C
      IF (allocated(bub)) deallocate(bub)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lioczon','readnbubczon')
 1001 CALL printerror2(0,1001,3,'lioczon','readnbubczon')
C
 102  WRITE (texterror,*) 'Parameters in .czon file header',
     $     ' incompatible with SESAM configuration'
      CALL printerror2(0,102,3,'lioczon','readnbubczon',comment=texterror)
 103  WRITE (texterror,*) 'Variable names in .czon file',
     $     ' incompatible with SESAM configuration'
      CALL printerror2(0,103,3,'lioczon','readnbubczon',comment=texterror)
 104  WRITE (texterror,*) 'No centering/reducing',
     $     ' is possible with readnbubzon'
      CALL printerror2(0,104,3,'lioczon','readnbubczon',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writehdrczon (kfnoutzon,kzon_jpi,kzon_jpj,kzon_jpk,
     $           kzon_jpt,kjpbub,kjpz)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write 'czon' file header
CCC  -------
CCC  Method :  Call low level routine defining czon format (NetCDF)
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
      use mod_main
      use mod_cfgxyo
      use utilcdfzon
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
      INTEGER :: jdta,inddta
      INTEGER, dimension(1:nbvar) :: dta_dim1,dta_nbr1
      CHARACTER(len=varlg), dimension(1:nbvar) :: dta_nam1
      BIGREAL4, dimension(1:nbvar) :: dta_moy1,dta_ect1
      CHARACTER(len=word80) :: title
CC----------------------------------------------------------------------
C
C Set header variables
C
      title = 'SESAM Vz object'
C
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
C
C Write 'czon' file header
C
      CALL cdfwdimzon(kfnoutzon,kzon_jpi,kzon_jpj,kzon_jpk,kzon_jpt,
     $               kjpbub,kjpz,dtaend,varlg,title)
C
      CALL cdfwhdrzon(kfnoutzon,dta_nam1,dta_dim1,dta_nbr1,
     $     dta_moy1,dta_ect1)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lioczon','writehdrczon')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writenbubczon (kfnoutzon,kbub,kjbub)
CCC---------------------------------------------------------------------
CCC   
CCC  Purpose : Write list of local data sections in 'czon' file
CCC  -------
CCC  Method : Call low level routine defining czon format (NetCDF)
CCC  ------   to write selected local data sections
CCC
CCC  Input : kfnoutzon : filename
CCC  -----   kjbub     : index of first local data section to write
CCC          kbub      : array of local data sections to write
CCC
CCC---------------------------------------------------------------------
CC modules
CC ======= 
      use mod_main
      use mod_cfgxyo
      use utilcdfzon
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfnoutzon
      BIGREAL, dimension(:,:,:,:,:), intent(in) :: kbub
      INTEGER, intent(in) :: kjbub
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL4, dimension(:,:,:,:,:), allocatable :: kbub4
      INTEGER :: jpisize,jpjsize,jpksize,jptsize,
     $     jpzsize,jpnbubsize,bubcount
      INTEGER :: allocok
CC----------------------------------------------------------------------
      jpisize = size(kbub,1)
      jpjsize = size(kbub,2)
      jpksize = size(kbub,3)
      jptsize = size(kbub,4)
      jpnbubsize=size(kbub,5)
C allocate local data section array to write (kr4 real kind)
      allocate ( kbub4(1:jpisize,1:jpjsize,1:jpksize,
     $     1:jptsize,1:jpnbubsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
      kbub4(:,:,:,:,:) = FREAL4(kbub(:,:,:,:,:))
C
C Write block of local data sections in 'czon' file
C
      bubcount=jpnbubsize
      CALL cdfwbubzon(kfnoutzon,kbub4,kjbub,bubcount)
C
      IF (allocated(kbub4)) deallocate (kbub4)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lioczon','writenbubczon')
 1001 CALL printerror2(0,1001,3,'lioczon','writenbubczon')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE writeptczon (kfnoutzon,ptbubidx,ptdtalon,ptdtalat,
     $           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth,
     $           ptbubtime)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Write Vz pointers in 'czon' file
CCC  -------
CCC  Method : Call low level routine defining czon format (NetCDF)
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
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use utilcdfzon
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
      INTEGER :: dtastart,jzstart,jzend
CC----------------------------------------------------------------------
C
C Check size of input arrays
C
      dtastart = 1
C
      IF (size(ptbubidx,2).NE.dtaend) GOTO 1000
      IF (size(ptbublon,2).NE.dtaend) GOTO 1000
      IF (size(ptbublat,2).NE.dtaend) GOTO 1000
      IF (size(ptbubtime,2).NE.dtaend) GOTO 1000
      IF (size(ptbubdepth,2).NE.dtaend) GOTO 1000
      IF (size(ptdtalon,2).NE.dtaend) GOTO 1000
      IF (size(ptdtalat,2).NE.dtaend) GOTO 1000
      IF (size(ptdtatime,2).NE.dtaend) GOTO 1000
      IF (size(ptdtadepth,2).NE.dtaend) GOTO 1000
C
      jzstart = 1
      jzend = size(ptbubidx,1)
C
      IF (size(ptbublon,1).NE.jzend) GOTO 1000
      IF (size(ptbublat,1).NE.jzend) GOTO 1000
      IF (size(ptbubtime,1).NE.jzend) GOTO 1000
      IF (size(ptbubdepth,1).NE.jzend) GOTO 1000
      IF (size(ptdtalon,1).NE.jzend) GOTO 1000
      IF (size(ptdtalat,1).NE.jzend) GOTO 1000
      IF (size(ptdtatime,1).NE.jzend) GOTO 1000
      IF (size(ptdtadepth,1).NE.jzend) GOTO 1000
C
C Write pointers in .czon file
C
      CALL cdfwptzon(kfnoutzon,ptbubidx,ptdtalon,ptdtalat,
     $       ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth,
     $       ptbubtime,jzstart,jzend,dtastart,dtaend)

C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'lioczon','writeptczon')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE lioczon
