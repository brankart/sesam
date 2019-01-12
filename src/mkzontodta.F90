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
! ---                    MKZONTODTA.F90                           ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-11 (J.M. Brankart)                      ---
! --- modification : 03-03 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! ---
! --- SUBROUTINE mkzontodta
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkzontodta
      use mod_main
      use mkconnect
      IMPLICIT NONE
      PRIVATE

      PUBLIC zontodta

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE zontodta (kfninzon,kfnoutdta,kjz)
!---------------------------------------------------------------------
!
!  Purpose : Interface local data sections corresponding to one
!  -------   subsystem into Vy vector
!
!  Method : Read local data sections from zon file, paste it
!  ------   on Vy vector, and write Vy vector in dta file
!
!  Input :  kfninzon  : Vz input filename
!  -----    kfnoutdta : Vy output filename
!           kjz       : subsystem index
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpyend,jpz,ptbub,bubblk1
      use hioxyo
      use hiozon
      use utilmkto
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninzon,kfnoutdta
      INTEGER, intent(in) :: kjz
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vecty
!
      INTEGER :: allocok,jpysize
      INTEGER :: ji,jj,jk,jt,jdta,jbub
      LOGICAL :: lectinfo
      INTEGER, dimension (:,:), allocatable :: ptbubidx
      INTEGER, dimension (:,:), allocatable :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), allocatable :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), allocatable :: ptbublon, ptbublat
      INTEGER, dimension (:,:), allocatable :: ptbubdepth, ptbubtime
      INTEGER, dimension (1:1) :: bublist
!----------------------------------------------------------------------
      jpysize=jpyend
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modzone/mkzontodta :'
         WRITE(numout,*) '         interface Vz vector to Vy vector'
         WRITE(numout,*) '         (one of the subsystems)'
      ENDIF
!
! Allocate Vy array
      allocate ( vecty(0:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecty(0:jpysize) = FREAL(0.0)
!
! Allocate pointers arrays
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
      ptbubidx(:,:) = 0
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
! Allocate local data section arrays
      allocate ( bubblk1(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &     1:1,1:1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      allocate ( ptbub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &     1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptbub(:,:,:,:,:) = 0
!
! Read pointers from Vz file
      CALL readptzon(kfninzon,ptbubidx,ptdtalon,ptdtalat, &
     &        ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &        ptbubtime)
!
! Loop on variables in Vy object, read corresponding
! local data section and paste it on Vy object
      lectinfo=.TRUE.
      DO jdta=1,dtaend
         CALL mkptbub (ptbub(:,:,:,:,jdta),jdta, &
     &                 ptdtalon(kjz,jdta),ptdtalat(kjz,jdta), &
     &                 ptdtadepth(kjz,jdta),ptdtatime(kjz,jdta), &
     &                 ptbublon(kjz,jdta),ptbublat(kjz,jdta), &
     &                 ptbubdepth(kjz,jdta),ptbubtime(kjz,jdta))
         bublist(1) = ptbubidx(kjz,jdta)
         CALL readnbubzon(kfninzon,bublist(:),bubblk1(:,:,:,:,:,1),lectinfo)
         DO jt=1,zon_jpt
         DO jk=1,zon_jpk
         DO jj=1,zon_jpj
         DO ji=1,zon_jpi
            vecty(ptbub(ji,jj,jk,jt,jdta)) = bubblk1(ji,jj,jk,jt,1,1)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         lectinfo=.FALSE.
      ENDDO
!
! Take into account possible connection lines
      IF (largconnect) CALL mkconnecty(vecty(0:jpysize))
!
! Write Vy vector in dta file
      CALL writedta (kfnoutdta,vecty(1:jpysize))
!
! --- deallocation
      IF (allocated(vecty)) deallocate(vecty)
      IF (allocated(ptbubidx)) deallocate(ptbubidx)
      IF (allocated(ptdtalon)) deallocate(ptdtalon)
      IF (allocated(ptdtalat)) deallocate(ptdtalat)
      IF (allocated(ptdtadepth)) deallocate(ptdtadepth)
      IF (allocated(ptdtatime)) deallocate(ptdtatime)
      IF (allocated(ptbublon)) deallocate(ptbublon)
      IF (allocated(ptbublat)) deallocate(ptbublat)
      IF (allocated(ptbubdepth)) deallocate(ptbubdepth)
      IF (allocated(ptbubtime)) deallocate(ptbubtime)
      IF (allocated(ptbub)) deallocate(ptbub)
      IF (allocated(bubblk1)) deallocate(bubblk1)
!
      RETURN
!
! --- error management
!
 1001 CALL printerror2(0,1001,3,'mkzontodta','zontodta')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkzontodta
