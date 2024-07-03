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
! ---                   ALGOINTF.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-11 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE algointfxyoz
! --- SUBROUTINE algointfxyozbas
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algointf
      use mod_main
      use mkdtatozon
      use utilconstraint
      use hiogrd
      IMPLICIT NONE
      PRIVATE

      PUBLIC algointfxyoz,algointfxyozbas

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algointfxyoz (kinxyoz,koutxyoz, &
     &     kflagxyoz,kconfigoz)
!---------------------------------------------------------------------
!
!  Purpose : Interface management of xyoz vector
!  -------
!  Method :
!  ------
!  Input :
!  -----
!  Output :
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_coord
      use mod_spacexyo , only : jpx,jpxend,jpyend,jpy,jpz, &
     &     gridijkobs,poscoefobs,jpoend,jpitpend
      use hioxyo
      use hiozon
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      CHARACTER(len=*), intent(in) :: kinxyoz,koutxyoz
      INTEGER, intent(in) :: kflagxyoz
      CHARACTER(len=*), intent(in), optional :: kconfigoz
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:), allocatable, save :: vecty
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      INTEGER :: allocok,jpssize,jpysize,jpitpsize
      INTEGER :: jnxyo,flagcfg,flagxyo,jpisize,jpjsize
      LOGICAL :: lectinfo,lmodprint,allmemjpbub
      INTEGER :: jbub,jbub1,jdta,kjdta,jz,kjz,jprbasout1
      INTEGER :: jdtafin,jpnzsize,jzdeb,jzfin,jz1,fixjpbub
      INTEGER :: zon_jpi1,zon_jpj1,zon_jpk1,zon_jpt1, &
     &     jpbub1,jpbubsize,jpz1,jbubfin
      INTEGER, dimension (:,:), allocatable :: ptbubidx,ptbubidx1
      INTEGER, dimension (:), allocatable :: vectptbub
      INTEGER, dimension (:,:), allocatable :: ptdtalon, ptdtalat, &
     &     ptdtalon1, ptdtalat1
      INTEGER, dimension (:,:), allocatable :: ptdtadepth, ptdtatime, &
     &     ptdtadepth1, ptdtatime1
      INTEGER, dimension (:,:), allocatable :: ptbublon, ptbublat, &
     &     ptbublon1, ptbublat1
      INTEGER, dimension (:,:), allocatable :: ptbubdepth, ptbubtime, &
     &     ptbubdepth1, ptbubtime1
      BIGREAL, dimension (:,:,:,:,:), allocatable :: bub
!----------------------------------------------------------------------
!
      jpysize=jpyend
      jpitpsize=1
      SELECT CASE (kflagxyoz)
      CASE (1)
         jpssize=jpx
      CASE (2)
         jpssize=jpyend
      CASE (3)
         jpssize=jpoend
         jpitpsize=jpitpend
      CASE (4)
         jpssize=jpyend
      CASE DEFAULT
         GOTO 1000
      END SELECT
! --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
!---------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&         routine algointfxyoz             &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      lectinfo = .TRUE.
!
      IF (((kflagxyoz.EQ.3).OR.(kflagxyoz.EQ.4)) &
     &     .AND.(.NOT.(present(kconfigoz)))) GOTO 1000
!
      SELECT CASE (kflagxyoz)
      CASE (1,2,4)
! ---
      CASE (3)
! --- allocation vectorms
         allocate ( vectorms(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectorms(:) = FREAL(0.0)
! --- allocation gridijkobs
         allocate ( gridijkobs(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
! --- allocation poscoefobs
         allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
!
! --- Reading config.obs
         flagcfg=1
         CALL readcfgobs (kconfigoz,flagcfg, &
     &        kvectorms=vectorms(:))
         flagcfg=2
         CALL readcfgobs (kconfigoz,flagcfg, &
     &        kgridijkobs=gridijkobs(:))
         flagcfg=3
         CALL readcfgobs (kconfigoz,flagcfg, &
     &        kposcoefobs=poscoefobs(:,:))
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -1.- Reading and writing xyoz object :
! --------------------------------------
!
      SELECT CASE (kflagxyoz)
      CASE (1,2,3)
         DO jnxyo=1,limjpnxyo(MIN(3,kflagxyoz))
         lmodprint=(MOD(jnxyo-1,(limjpnxyo(MIN(3,kflagxyoz))/5+1)).EQ.0)
         IF ((lmodprint).AND.(nprint.GE.1)) &
     &        print *,'Memory part number : ', &
     &           jnxyo,'/',limjpnxyo(MIN(3,kflagxyoz))
!
! -1.1- Reading xyoz object :
! --------------------------
!
         SELECT CASE (kflagxyoz)
         CASE (1,2)
            CALL readxyo(kinxyoz,vects(:), &
     &           jnxyo,lectinfo,kflagxyoz)
         CASE (3)
            CALL readxyo(kinxyoz,vects(:), &
     &           jnxyo,lectinfo,kflagxyoz,poscoefobs(:,:))
         CASE DEFAULT
            GOTO 1000
         END SELECT
!     
! -1.2- Apply constraint if requested
! -----------------------------------
!
#ifdef FLOWSAMPLER
         IF (dyn_constraint.AND.(kflagxyoz.EQ.1)) THEN
           IF (MAXVAL(varngrd(1:varend)).GT.1) GOTO 1000

           jpisize=MAXVAL(var_jpi(1:varend))
           jpjsize=MAXVAL(var_jpj(1:varend))

           allocate (longi(1:jpisize), stat=allocok )
           IF (allocok.NE.0) GOTO 1001
           longi(:) = FREAL(0.0)

           allocate (latj(1:jpjsize), stat=allocok )
           IF (allocok.NE.0) GOTO 1001
           latj(:) = FREAL(0.0)

           CALL readgrd (kflagxyoz,1)
           CALL apply_constraint(vects)

           IF (allocated(longi)) deallocate(longi)
           IF (allocated(latj)) deallocate(latj)

         ENDIF
#endif
!     
! -1.3- Writing xyoz object :
! --------------------------
!
         SELECT CASE (kflagxyoz)
         CASE (1)
            CALL writevar (koutxyoz,vects(:),jnxyo)
         CASE (2)
            CALL writedta (koutxyoz,vects(:))
         CASE (3)
            CALL writeobs(koutxyoz,vects(:),vectorms(:),gridijkobs(:), &
     &           poscoefobs(:,:))
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
         lectinfo = .FALSE.
!
         ENDDO
      CASE (4)
! ---
         IF (validextvar(kinxyoz).OR.validextdta(kinxyoz)) THEN
            CALL dtatozon (kinxyoz,kconfigoz,koutxyoz)
         ELSE
!
            CALL evalhdrzon (kconfigoz,zon_jpi,zon_jpj,zon_jpk, &
     &           zon_jpt,jpbub,jpz)
            fixjpbub=jpbub
! --- allocation zone pointers
            allocate ( ptbubidx(1:jpz,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptbubidx(:,:) = 0
!
            allocate ( ptdtalon(1:jpz,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptdtalon(:,:) = 0
            allocate ( ptdtalat(1:jpz,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptdtalat(:,:) = 0
            allocate ( ptdtadepth(1:jpz,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptdtadepth(:,:) = 0
            allocate ( ptdtatime(1:jpz,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptdtatime(:,:) = 0
!
            allocate ( ptbublon(1:jpz,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptbublon(:,:) = 0
            allocate ( ptbublat(1:jpz,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptbublat(:,:) = 0
            allocate ( ptbubdepth(1:jpz,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptbubdepth(:,:) = 0
            allocate ( ptbubtime(1:jpz,1:dtaend), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            ptbubtime(:,:) = 0
!
            CALL readptzon (kconfigoz,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime)
            CALL writehdrzon (koutxyoz,zon_jpi,zon_jpj,zon_jpk, &
     &           zon_jpt,jpbub,jpz)
            CALL writeptzon (koutxyoz,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime)
!
            IF (kconfigoz.EQ.kinxyoz) THEN
               allmemjpbub=(jpbub.LE.fixjpbub)
               IF (allmemjpbub) THEN
                  jpbubsize=jpbub
               ELSE
                  jpbubsize=((fixjpbub/dtaend)+1)*dtaend
               ENDIF
! --- allocation bub
               allocate ( bub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &              1:jpbubsize), &
     &              stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               bub(:,:,:,:,:) = FREAL(0.0)
               lectinfo=.TRUE.
! --- allocation vectptbub
               allocate ( vectptbub(1:jpbubsize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               IF (allmemjpbub) THEN
                  vectptbub(:) = (/ ( jbub, jbub=1,jpbubsize) /)
                  CALL readnbubzon(kinxyoz,vectptbub(:), &
     &                 bub(:,:,:,:,:),lectinfo)
                  jbub =1
                  CALL writenbubzon(koutxyoz, &
     &                 bub(:,:,:,:,:),jbub)
               ELSE
                  DO jbub=1,jpbub,jpbubsize
                     jbubfin=MIN(jpbubsize,jpbub-(jbub-1))
                     vectptbub(1:jbubfin) =  &
     &                    (/ ( jbub+(jbub1-1) , jbub1=1,jbubfin) /)
                     CALL readnbubzon(kinxyoz, &
     &                    vectptbub(1:jbubfin), &
     &                    bub(:,:,:,:,1:jbubfin),lectinfo)
                     lectinfo=.FALSE.
                     CALL writenbubzon(koutxyoz, &
     &                    bub(:,:,:,:,1:jbubfin),jbub)
                 ENDDO
               ENDIF
            ELSE
               IF (jpbub.NE.(jpz*dtaend)) GOTO 103
               jbub=0
               DO jz=1,jpz
               DO jdta=1,dtaend
                 jbub=jbub+1
                 ptbubidx(jz,jdta)=jbub
               ENDDO
               ENDDO
               CALL evalhdrzon (kinxyoz,zon_jpi1,zon_jpj1,zon_jpk1, &
     &              zon_jpt1,jpbub1,jpz1)
               IF (zon_jpi1.NE.zon_jpi) GOTO 102
               IF (zon_jpj1.NE.zon_jpj) GOTO 102
               IF (zon_jpk1.NE.zon_jpk) GOTO 102
               IF (zon_jpt1.NE.zon_jpt) GOTO 102
               IF (jpz1.NE.jpz) GOTO 102
! --- allocation zone pointers
               allocate ( ptbubidx1(1:jpz,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               ptbubidx1(:,:) = 0
              allocate ( ptdtalon1(1:jpz,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               ptdtalon1(:,:) = 0
               allocate ( ptdtalat1(1:jpz,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               ptdtalat1(:,:) = 0
               allocate ( ptdtadepth1(1:jpz,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               ptdtadepth1(:,:) = 0
               allocate ( ptdtatime1(1:jpz,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               ptdtatime1(:,:) = 0
               allocate ( ptbublon1(1:jpz,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               ptbublon1(:,:) = 0
               allocate ( ptbublat1(1:jpz,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               ptbublat1(:,:) = 0
               allocate ( ptbubdepth1(1:jpz,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               ptbubdepth1(:,:) = 0
               allocate ( ptbubtime1(1:jpz,1:dtaend), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               ptbubtime1(:,:) = 0
               CALL readptzon (kinxyoz,ptbubidx1,ptdtalon1,ptdtalat1, &
     &              ptdtadepth1,ptdtatime1,ptbublon1,ptbublat1, &
     &              ptbubdepth1,ptbubtime1)
               IF (ANY(ptdtalon1(:,:).NE.ptdtalon(:,:))) GOTO 102
               IF (ANY(ptdtalat1(:,:).NE.ptdtalat(:,:))) GOTO 102
               IF (ANY(ptdtadepth1(:,:).NE.ptdtadepth(:,:))) GOTO 102
               IF (ANY(ptdtatime1(:,:).NE.ptdtatime(:,:))) GOTO 102
               IF (ANY(ptbublon1(:,:).NE.ptbublon(:,:))) GOTO 102
               IF (ANY(ptbublat1(:,:).NE.ptbublat(:,:))) GOTO 102
               IF (ANY(ptbubdepth1(:,:).NE.ptbubdepth(:,:))) GOTO 102
               IF (ANY(ptbubtime1(:,:).NE.ptbubtime(:,:))) GOTO 102
!
               IF (allocated(ptdtalon1)) deallocate(ptdtalon1)
               IF (allocated(ptdtalat1)) deallocate(ptdtalat1)
               IF (allocated(ptdtadepth1)) deallocate(ptdtadepth1)
               IF (allocated(ptdtatime1)) deallocate(ptdtatime1)
               IF (allocated(ptbublon1)) deallocate(ptbublon1)
               IF (allocated(ptbublat1)) deallocate(ptbublat1)
               IF (allocated(ptbubdepth1)) deallocate(ptbubdepth1)
               IF (allocated(ptbubtime1)) deallocate(ptbubtime1)
!
               allmemjpbub=(jpbub.LE.fixjpbub)
               IF (allmemjpbub) THEN
                  jpbubsize=jpbub
               ELSE
                  jpbubsize=((fixjpbub/dtaend)+1)*dtaend
               ENDIF
! --- allocation bub
               allocate ( bub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &              1:jpbubsize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               bub(:,:,:,:,:) = FREAL(0.0)
               lectinfo=.TRUE.
! --- allocation vectptbub
               allocate ( vectptbub(1:jpbubsize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               IF (allmemjpbub) THEN
                  vectptbub(:) =  &
     &               (/ ((jdta+(jz-1)*dtaend, jdta=1,dtaend),jz=1,jpz) /)
                  CALL readnbubzon(kinxyoz,vectptbub(:), &
     &                 bub(:,:,:,:,:),lectinfo)
                  jbub =1
                  CALL writenbubzon(koutxyoz, &
     &                 bub(:,:,:,:,:),jbub)
               ELSE
                  DO jbub=1,jpbub,jpbubsize
                     jbubfin=MIN(jpbubsize,jpbub-(jbub-1))
                     jzdeb=(jbub-1)/dtaend
                     jzfin=(jbubfin-1)/dtaend
                     vectptbub(1:jbubfin) =  &
     &                    (/ ((ptbubidx1(jz,jdta), &
     &                    jdta=1,dtaend),jz=jzdeb,jzfin) /)
                     CALL readnbubzon(kinxyoz, &
     &                    vectptbub(1:jbubfin), &
     &                    bub(:,:,:,:,1:jbubfin),lectinfo)
                     lectinfo=.FALSE.
                      CALL writenbubzon(koutxyoz, &
     &                    bub(:,:,:,:,1:jbubfin),jbub)
                  ENDDO
               ENDIF
            ENDIF               
         ENDIF               
!
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! --- deallocate
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(bub)) deallocate(bub)
      IF (allocated(ptbubidx)) deallocate(ptbubidx)
      IF (allocated(ptdtalon)) deallocate(ptdtalon)
      IF (allocated(ptdtalat)) deallocate(ptdtalat)
      IF (allocated(ptdtadepth)) deallocate(ptdtadepth)
      IF (allocated(ptdtatime)) deallocate(ptdtatime)
      IF (allocated(ptbublon)) deallocate(ptbublon)
      IF (allocated(ptbublat)) deallocate(ptbublat)
      IF (allocated(ptbubdepth)) deallocate(ptbubdepth)
      IF (allocated(ptbubtime)) deallocate(ptbubtime)
      IF (allocated(ptbubidx1)) deallocate(ptbubidx1)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algointf','algointfxyoz')
 1001 CALL printerror2(0,1001,3,'algointf','algointfxyoz')
!
 102  WRITE (texterror,*) 'file ',kinxyoz(1:lenv(kinxyoz)), &
     &     ' incoherent with ',kconfigoz(1:lenv(kconfigoz))
      CALL printerror2(0,102,3,'algointf','algointfxyoz', &
     &     comment=texterror)
 103  WRITE (texterror,*) 'file ',kconfigoz(1:lenv(kconfigoz)), &
     &     ' not allowed for transfert on ',kinxyoz(1:lenv(kinxyoz))
      CALL printerror2(0,103,3,'algointf','algointfxyoz', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algointfxyozbas (kinxyozbas,koutxyozbas, &
     &     kflagxyoz,kconfigoz)
!---------------------------------------------------------------------
!
!  Purpose : Interface management of xyoz vector
!  -------
!  Method :
!  ------
!  Input :
!  -----
!  Output :
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_coord
      use mod_spacexyo , only : &
     &     gridijkobs,poscoefobs,jpoend,jpitpend, &
     &     jpx,jpxend,jpyend,jpy,jpz,jprend
      use hioxyo
      use hiobas
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinxyozbas,koutxyozbas
      INTEGER, intent(in) :: kflagxyoz
      CHARACTER(len=*), intent(in), optional :: kconfigoz
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      INTEGER :: allocok,jpssize,jpysize,jpitpsize,jprsize
      INTEGER :: jnxyo,flagcfg
      LOGICAL :: lectinfo
      INTEGER :: valbase
      INTEGER :: serie,numjr,jprbas,jprbasin,jprbasout
      CHARACTER(len=bgword) :: fnamein,fnameout,vctnamin,vctnamout
!----------------------------------------------------------------------
!
      jprsize=jprend
      jpysize=jpyend
      jpitpsize=1
      SELECT CASE (kflagxyoz)
      CASE (1)
         jpssize=jpx
      CASE (2)
         jpssize=jpyend
      CASE (3)
         jpssize=jpoend
         jpitpsize=jpitpend
      CASE (4)
         jpssize=jpz
      CASE DEFAULT
         GOTO 1000
      END SELECT
      IF (kflagxyoz.NE.4) THEN
! --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
      ENDIF
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&       routine algointfxyozbas            &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      IF (((kflagxyoz.EQ.3).OR.(kflagxyoz.EQ.4)) &
     &     .AND.(.NOT.(present(kconfigoz)))) GOTO 1000
!
      SELECT CASE (kflagxyoz)
      CASE (1,2,4)
! ---
      CASE (3)
! --- allocation vectorms
         allocate ( vectorms(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectorms(:) = FREAL(0.0)
! --- allocation gridijkobs
         allocate ( gridijkobs(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
! --- allocation poscoefobs
         allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
!
! --- Reading config.obs
         flagcfg=1
         CALL readcfgobs (kconfigoz,flagcfg, &
     &        kvectorms=vectorms(:))
         flagcfg=2
         CALL readcfgobs (kconfigoz,flagcfg, &
     &        kgridijkobs=gridijkobs(:))
         flagcfg=3
         CALL readcfgobs (kconfigoz,flagcfg, &
     &        kposcoefobs=poscoefobs(:,:))
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      lectinfo = .TRUE.
      valbase=0
!
! -0.1- Control test :
! --------------------
!
      serie=0
      numjr=1
      CALL fildirbas (fnamein,kinxyozbas,jprbasin,numjr,serie)
      CALL fildirbas (fnameout,koutxyozbas,jprbasout,numjr,serie)
      IF (jprbasin.GT.jprend) GOTO 1000
      IF (jprbasout.GT.jprbasin) GOTO 101
      IF (jprbasout.NE.jprbasin) THEN
         IF (nprint.GE.1) THEN
           WRITE(numout,*) ' '
           WRITE(numout,*) ' Attention :' 
           WRITE(numout,*) ' la taille de la base en sortie'
           WRITE(numout,*) ' est inferieur a celle de la base en entree'
           WRITE(numout,*) ' jprbas entree :',jprbasin
           WRITE(numout,*) ' jprbas sortie :',jprbasout
           WRITE(numout,*) ' ==> troncature de la base'
           WRITE(numout,*) ' '
         ENDIF 
      ENDIF
!
! -1.1- Read the infos :
! ----------------------
!
      CALL readinfobas(kinxyozbas,valbase)
!
! -1.2- Write the infos :
! -----------------------
!
      CALL writeinfobas(koutxyozbas,valbase)
!
! -1.- Reading and writing xyoz object :
! --------------------------------------
!
      serie=1
      DO numjr=1,jprbasout
!
         CALL fildirbas (vctnamin,kinxyozbas,jprbasin,numjr,serie)
         WRITE(fnamein,'("./",A,"/",A)') kinxyozbas(1:lenv(kinxyozbas)), &
     &        vctnamin(1:lenv(vctnamin))
!
         CALL fildirbas (vctnamout,koutxyozbas,jprbas,numjr,serie)
         WRITE(fnameout,'("./",A,"/",A)')  &
     &        koutxyozbas(1:lenv(koutxyozbas)), &
     &        vctnamout(1:lenv(vctnamout))
!     
         SELECT CASE (kflagxyoz)
         CASE (1,2,3)
! ---
         DO jnxyo=1,limjpnxyo(MIN(3,kflagxyoz))
!
! -1.1- Reading xyoz object :
! --------------------------
!
            SELECT CASE (kflagxyoz)
            CASE (1,2)
               CALL readxyo(fnamein,vects(:), &
     &              jnxyo,lectinfo,kflagxyoz)
            CASE (3)
               CALL readxyo(fnamein,vects(:), &
     &              jnxyo,lectinfo,kflagxyoz,poscoefobs(:,:))
            CASE (4)
! ---
            CASE DEFAULT
               GOTO 1000
            END SELECT
!     
! -1.1- Writing xyoz object :
! --------------------------
!
            SELECT CASE (kflagxyoz)
            CASE (1)
               CALL writevar (fnameout,vects(:),jnxyo)
            CASE (2)
               CALL writedta (fnameout,vects(:))
            CASE (3)
               CALL writeobs(fnameout,vects(:),vectorms(:), &
     &              gridijkobs(:),poscoefobs(:,:))
            CASE (4)
! ---
            CASE DEFAULT
               GOTO 1000
            END SELECT
!
            lectinfo = .FALSE.
!     
         ENDDO
         CASE (4)
! ---
            CALL algointfxyoz (fnamein,fnameout, &
     &           kflagxyoz,kconfigoz)
!
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
!     ==> control print
         IF (nprint.GE.2) THEN
            WRITE(numout,*) ' ==>  vectin  colonne :',vctnamin
            WRITE(numout,*) ' ==>  vectout colonne :',vctnamout
         ENDIF
!
      ENDDO
!
! --- deallocate
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algointf','algointfxyozbas')
 1001 CALL printerror2(0,1001,3,'algointf','algointfxyozbas')
!
 101  WRITE (texterror,*) 'bad dimension'
      CALL printerror2(0,101,3,'algointf','algointfxyozbas', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algointf
