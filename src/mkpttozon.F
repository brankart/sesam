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
C ---                    MKPTTOZON.F                            ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 00-01 (C.E. Testut)                        ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C ---
C --- SUBROUTINE mkpttozon
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkpttozon
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC pttozon

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE pttozon (kfninzon,kfnoutptzon)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
CCC  -------
CCC  Method :
CCC  ------
CCC  Input :           : no
CCC  -----
CCC  Output :          : no
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpz, pt1bubidx,
     $     pt1dtalon, pt1dtalat, pt1dtadepth, pt1dtatime,
     $     pt1bublon, pt1bublat, pt1bubdepth, pt1bubtime
      use hiozon
      use utilmkto
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ==================
      CHARACTER(len=*), intent(in) :: kfninzon,kfnoutptzon
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpbubsize
      LOGICAL :: lectinfo,lmoyectold
      BIGREAL, dimension (:,:,:,:,:), allocatable :: bub
      INTEGER, dimension (:,:,:,:,:), allocatable :: ptbub
      INTEGER :: jz,jbub,jdta
CC----------------------------------------------------------------------
C
      jpbubsize = jpz * dtaend
C --- allocation bub
      allocate ( bub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,
     $     1:jpbubsize),stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      bub(:,:,:,:,:) = FREAL(0.0)
C --- allocation ptbub
      allocate ( ptbub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,
     $     1:jpbubsize),stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptbub(:,:,:,:,:) = 0
C --- allocation zone pointers
      allocate ( pt1bubidx(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bubidx(:,:) = 0
C
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
C
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
CCC---------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mkpttozon &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
C
      lmoyectold=lmoyect
      lmoyect=.FALSE.
C
C -1.- Write the pt.zon file header :
C -----------------------------------
C
      CALL writehdrzon (kfnoutptzon,zon_jpi,zon_jpj,zon_jpk,
     $           zon_jpt,jpbubsize,jpz)
C
C -2.- Read the .zon pointers :
C -----------------------------
      CALL readptzon  (kfninzon,pt1bubidx,pt1dtalon,pt1dtalat,
     $           pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth,
     $           pt1bubtime)
C
C -3.- Create the ptzon :
C -----------------------
C
      jbub = 0
      DO jz = 1,jpz
      DO jdta=1,dtaend
         jbub = jbub + 1
C
         CALL mkptbub (ptbub(:,:,:,:,jbub),jdta,
     $        pt1dtalon(jz,jdta),pt1dtalat(jz,jdta),
     $        pt1dtadepth(jz,jdta),pt1dtatime(jz,jdta),
     $        pt1bublon(jz,jdta),pt1bublat(jz,jdta),
     $        pt1bubdepth(jz,jdta),pt1bubtime(jz,jdta))
         pt1bubidx(jz,jdta) = jbub
C
      ENDDO
      ENDDO
      bub(:,:,:,:,:) = FREAL(ptbub(:,:,:,:,:))
C
C -4.- Write the pt.czon file :
C -----------------------------
C
      jbub=1
      CALL writenbubzon (kfnoutptzon,bub(:,:,:,:,1:jpbubsize),jbub)
      CALL writeptzon (kfnoutptzon,pt1bubidx,pt1dtalon,pt1dtalat,
     $           pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth,
     $           pt1bubtime)
C
      lmoyect=lmoyectold
C
C --- deallocation
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
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mkpttozon','pttozon')
 1001 CALL printerror2(0,1001,3,'mkpttozon','pttozon')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkpttozon
