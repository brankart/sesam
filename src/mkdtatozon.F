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
C ---                    MKDTATOZON.F                           ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-11 (J.M. Brankart)                      ---
C --- modification : 00-03 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                        ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C ---
C --- SUBROUTINE dtatozon
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkdtatozon
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC dtatozon

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE dtatozon (kfnindta,kfninzon,kfnoutzon)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
CCC  -------
CCC     Interface a .dta file into a .zon file
CCC           using a predefined zone configuration
CCC  Method :
CCC  ------
CCC  Input :
CCC  -----
CCC  Output :
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpz, jpyend, pt1bubidx,
     $     pt1dtalon, pt1dtalat, pt1dtadepth, pt1dtatime,
     $     pt1bublon, pt1bublat, pt1bubdepth, pt1bubtime
      use hioxyo
      use hiozon
      use utilmkto
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninzon,kfnindta,kfnoutzon
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vecty
C
      INTEGER :: allocok,jpysize,jnxyo,flagxyo,jpbubsize
      LOGICAL :: lectinfo,lmoyectold
      BIGREAL, dimension (:,:,:,:,:), allocatable :: bub
C     $     ,maty
      INTEGER :: ji,jj,jk,jt,jz,jbub,jdta
C     $     ,inddta,inddtamsk
C      INTEGER :: jpi,jpj,jpk,jpt,jy
C      INTEGER :: zshift_ji1,zshift_jj1,zshift_jk1,zshift_jt1
CC----------------------------------------------------------------------
C
      jpysize=jpyend
      jpbubsize = jpz * dtaend
C --- allocation vecty
      allocate ( vecty(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecty(:) = FREAL(0.0)
C --- allocation bub
      allocate ( bub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,
     $     1:jpbubsize),stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      bub(:,:,:,:,:) = FREAL(0.0)
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
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&            routine mkdtatozon            &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
C
C -1.- Read the .dta file
C -----------------------
C
      lmoyectold=lmoyect
      lmoyect=.FALSE.
      lectinfo = .TRUE.
      jnxyo=1
      flagxyo=2
C      CALL readdta(kfnindta,vecty,lectinfo)
      CALL readxyo(kfnindta,vecty(:),jnxyo,lectinfo,flagxyo)
C
C -2.- Write the .zon file header
C -------------------------------
C
      CALL writehdrzon (kfnoutzon,zon_jpi,zon_jpj,zon_jpk,
     $           zon_jpt,jpbubsize,jpz)
C
C -3.- Read the .zon pointers
C ---------------------------
      CALL readptzon  (kfninzon,pt1bubidx,pt1dtalon,pt1dtalat,
     $           pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth,
     $           pt1bubtime)
C
C -4.- Select the y portion for each zone to write
C ------------------------------------------------
      jbub = 0
      DO jz = 1,jpz
      DO jdta = 1,dtaend
         jbub = jbub + 1
         CALL mkdtatobub (vecty(:),bub(:,:,:,:,jbub),jdta,
     $        pt1dtalon(jz,jdta),pt1dtalat(jz,jdta),
     $        pt1dtadepth(jz,jdta),pt1dtatime(jz,jdta),
     $        pt1bublon(jz,jdta),pt1bublat(jz,jdta),
     $        pt1bubdepth(jz,jdta),pt1bubtime(jz,jdta))
         pt1bubidx(jz,jdta) = jbub
      ENDDO
      ENDDO
C
C -5.- Write the .zon pointers
C ----------------------------
      jbub=1
      CALL writenbubzon (kfnoutzon,bub(:,:,:,:,1:jpbubsize),jbub)
      CALL writeptzon (kfnoutzon,pt1bubidx,pt1dtalon,pt1dtalat,
     $           pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth,
     $           pt1bubtime)
C
      lmoyect=lmoyectold
C
C --- desallocation
      IF (allocated(vecty)) deallocate(vecty)
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
 1000 CALL printerror2(0,1000,1,'mkdtatozon','dtatozon')
 1001 CALL printerror2(0,1001,3,'mkdtatozon','dtatozon')
C
 101  WRITE (texterror,*) 'argument not valid :'
      CALL printerror2(0,101,3,'mkdtatozon','dtatozon',
     $     comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkdtatozon
