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
! ---                    MKPTTOZON.F90                            ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 00-01 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! ---
! --- SUBROUTINE mkpttozon
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkpttozon
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC pttozon

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE pttozon (kfninzon,kfnoutptzon)
!---------------------------------------------------------------------
!
!  Purpose :
!  -------
!  Method :
!  ------
!  Input :           : no
!  -----
!  Output :          : no
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpz, pt1bubidx, &
     &     pt1dtalon, pt1dtalat, pt1dtadepth, pt1dtatime, &
     &     pt1bublon, pt1bublat, pt1bubdepth, pt1bubtime
      use hiozon
      use utilmkto
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      CHARACTER(len=*), intent(in) :: kfninzon,kfnoutptzon
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jpbubsize
      LOGICAL :: lectinfo,lmoyectold
      BIGREAL, dimension (:,:,:,:,:), allocatable :: bub
      INTEGER, dimension (:,:,:,:,:), allocatable :: ptbub
      INTEGER :: jz,jbub,jdta
!----------------------------------------------------------------------
!
      jpbubsize = jpz * dtaend
! --- allocation bub
      allocate ( bub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &     1:jpbubsize),stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      bub(:,:,:,:,:) = FREAL(0.0)
! --- allocation ptbub
      allocate ( ptbub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &     1:jpbubsize),stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptbub(:,:,:,:,:) = 0
! --- allocation zone pointers
      allocate ( pt1bubidx(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bubidx(:,:) = 0
!
      allocate ( pt1dtalon(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1dtalon(:,:) = 0
      allocate ( pt1dtalat(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1dtalat(:,:) = 0
      allocate ( pt1dtadepth(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1dtadepth(:,:) = 0
      allocate ( pt1dtatime(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1dtatime(:,:) = 0
!
      allocate ( pt1bublon(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bubidx(:,:) = 0
      allocate ( pt1bublat(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bublat(:,:) = 0
      allocate ( pt1bubdepth(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bubdepth(:,:) = 0
      allocate ( pt1bubtime(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bubtime(:,:) = 0
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mkpttozon &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
      lmoyectold=lmoyect
      lmoyect=.FALSE.
!
! -1.- Write the pt.zon file header :
! -----------------------------------
!
      CALL writehdrzon (kfnoutptzon,zon_jpi,zon_jpj,zon_jpk, &
     &           zon_jpt,jpbubsize,jpz)
!
! -2.- Read the .zon pointers :
! -----------------------------
      CALL readptzon  (kfninzon,pt1bubidx,pt1dtalon,pt1dtalat, &
     &           pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth, &
     &           pt1bubtime)
!
! -3.- Create the ptzon :
! -----------------------
!
      jbub = 0
      DO jz = 1,jpz
      DO jdta=1,dtaend
         jbub = jbub + 1
!
         CALL mkptbub (ptbub(:,:,:,:,jbub),jdta, &
     &        pt1dtalon(jz,jdta),pt1dtalat(jz,jdta), &
     &        pt1dtadepth(jz,jdta),pt1dtatime(jz,jdta), &
     &        pt1bublon(jz,jdta),pt1bublat(jz,jdta), &
     &        pt1bubdepth(jz,jdta),pt1bubtime(jz,jdta))
         pt1bubidx(jz,jdta) = jbub
!
      ENDDO
      ENDDO
      bub(:,:,:,:,:) = FREAL(ptbub(:,:,:,:,:))
!
! -4.- Write the pt.czon file :
! -----------------------------
!
      jbub=1
      CALL writenbubzon (kfnoutptzon,bub(:,:,:,:,1:jpbubsize),jbub)
      CALL writeptzon (kfnoutptzon,pt1bubidx,pt1dtalon,pt1dtalat, &
     &           pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth, &
     &           pt1bubtime)
!
      lmoyect=lmoyectold
!
! --- deallocation
      IF (allocated(bub)) deallocate(bub)
      IF (allocated(pt1bubidx)) deallocate(pt1bubidx)
      IF (allocated(pt1dtalon)) deallocate(pt1dtalon)
      IF (allocated(pt1dtalat)) deallocate(pt1dtalat)
      IF (allocated(pt1dtadepth)) deallocate(pt1dtadepth)
      IF (allocated(pt1dtatime)) deallocate(pt1dtatime)
      IF (allocated(pt1bublon)) deallocate(pt1bublon)
      IF (allocated(pt1bublat)) deallocate(pt1bublat)
      IF (allocated(pt1bubdepth)) deallocate(pt1bubdepth)
      IF (allocated(pt1bubtime)) deallocate(pt1bubtime)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkpttozon','pttozon')
 1001 CALL printerror2(0,1001,3,'mkpttozon','pttozon')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkpttozon
