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
! ---                   ALGOFILT.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2000-05 (J.M. Brankart)                    ---
! --- modification : 2001-11 (C.E. Testut)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! ---
! --- SUBROUTINE calcfilt
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algofilt
      use mod_main
      use mkconnect
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcfilt

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcfilt(kouty,kflagyo,kinyo,kinzon,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose : Filter a dta file 
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
      use mod_spacexyo , only : jpyend, jpz, jpitpend, jpoend, &
     &     pt1bubidx, pt4bubidx, &
     &     pt1dtalon, pt1dtalat, pt1dtadepth, pt1dtatime, &
     &     pt1bublon, pt1bublat, pt1bubdepth, pt1bubtime, &
     &     ptbub, bubblk1, bubblk4, vectptbub, &
     &     poscoefobs, spvaldta
      use hioxyo
      use hiozon
      use utilmkh
      use utilmkto
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(in) :: kflagyo
      CHARACTER(len=*), intent(in) :: kinyo,kinzon,kouty,kconfigo
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:), allocatable, save :: vectsfilt
      BIGREAL, dimension(:), allocatable, save :: vecty
      BIGREAL, dimension(:), allocatable, save :: vectybub
!
      INTEGER :: jy ,jz, jpu, jbub, ji, jj, jk, jt, jdta, jns, jnxyo
      INTEGER :: jpysize, jpzsize, jpssize, jpitpsize, js, flagcfg
      LOGICAL :: lectinfo, zlectinfo, lmoyectold, lmodprint
      INTEGER, dimension(:), allocatable :: ptu,jputab
      INTEGER :: allocok, flagxyo
      TYPE (type_poscoef), dimension(:,:), allocatable :: poscoefuobs
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine algofilt &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Initialisation
! -------------------
!
      SELECT CASE (kflagyo)
      CASE(2)
         jpssize=jpyend
         jpitpsize=1
      CASE(3)
         jpssize=jpoend
         jpitpsize=jpitpend
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      jpysize = jpyend
      jpzsize = jpz
      lectinfo = .FALSE.
      flagxyo = kflagyo
      jns = 1
!
! --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
!
! --- allocation vectsfilt
      allocate ( vectsfilt(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectsfilt(:) = FREAL(0.0)
!
! --- allocation vecty
      allocate ( vecty(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecty(:) = FREAL(0.0)
!
! --- allocation poscoefobs
      allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      IF (kflagyo.EQ.2) THEN
         DO js=1,jpssize
            poscoefobs(js,:) = type_poscoef(js,FREAL(1.0))
         ENDDO
      ELSE
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
         flagcfg=3
         CALL readcfgobs (kconfigo,flagcfg, &
     &        kposcoefobs=poscoefobs(:,:))
      ENDIF
!
! --- allocation jputab
      allocate ( jputab(1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      jputab(:) = 0
!
! -2.- Read input dta file
! ------------------------
!
      jnxyo=1
      CALL readxyo(kinyo,vects(:), &
     &     jnxyo,lectinfo,kflagyo,poscoefobs(:,:))
!
! -3.- Prepare the filter bubbles
! -------------------------------
!
      jpbubend(:)=0
      CALL evalhdrzon(kinzon,zon_jpi,zon_jpj,zon_jpk,zon_jpt, &
     &     jpbubend(1),jpz)
!
      IF (larginptzon) THEN
! --- jpbubend(4) by symmetry with algoroa
         CALL evalhdrzon(arginptzon,zon_jpi,zon_jpj,zon_jpk,zon_jpt, &
     &        jpbubend(4),jpz)
      ENDIF
!
! --- allocation zone pointers
      allocate ( pt1dtalon(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1dtalon(:,:) = 0
      allocate ( pt1dtalat(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1dtalat(:,:) = 0
      allocate ( pt1dtadepth(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1dtadepth(:,:) = 0
      allocate ( pt1dtatime(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1dtatime(:,:) = 0
!
      allocate ( pt1bublon(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bublon(:,:) = 0
      allocate ( pt1bublat(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bublat(:,:) = 0
      allocate ( pt1bubdepth(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bubdepth(:,:) = 0
      allocate ( pt1bubtime(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bubtime(:,:) = 0
!
      allocate ( pt1bubidx(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bubidx(:,:) = 0
!
      CALL readptzon (kinzon,pt1bubidx,pt1dtalon,pt1dtalat, &
     &      pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth, &
     &      pt1bubtime)
!
      IF (larginptzon) THEN
         allocate ( pt4bubidx(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         pt4bubidx(:,:) = 0
!
         CALL readptzon (arginptzon,pt4bubidx,pt1dtalon,pt1dtalat, &
     &        pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth, &
     &        pt1bubtime)
      ENDIF
!
! --- allocation vectptbub
      allocate ( vectptbub(1:MAXVAL(jpbubend)), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectptbub(:) = 0
!
! --- allocation bubblk1
      allocate ( bubblk1(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &     1:jpbubend(1),1:1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      bubblk1(:,:,:,:,:,:) = FREAL(0.0)
!
      zlectinfo = .TRUE.
      vectptbub(1:jpbubend(1)) = (/ (jbub, jbub = 1,jpbubend(1)) /)
      CALL readnbubzon(kinzon,vectptbub(1:jpbubend(1)), &
     &                    bubblk1(:,:,:,:,1:jpbubend(1),1),zlectinfo)
!
! --- allocation bubblk4
      IF (larginptzon) THEN
         allocate ( bubblk4(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &        1:jpbubend(4),1:1), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         bubblk4(:,:,:,:,:,:) = FREAL(0.0)
!
         zlectinfo = .TRUE.
         vectptbub(1:jpbubend(4)) = (/ (jbub, jbub = 1,jpbubend(4)) /)
         CALL readnbubzon(arginptzon,vectptbub(1:jpbubend(4)), &
     &        bubblk4(:,:,:,:,1:jpbubend(4),1),zlectinfo)
      ENDIF
!
! --- allocation ptbub
      allocate ( ptbub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &     1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptbub(:,:,:,:,:) = 0
!
! --- allocation vectybub
      allocate ( vectybub(0:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectybub(0:) = FREAL(0.0)
!
! --- allocation ptu
      allocate ( ptu(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptu(:) = 0
!
! -4.- Run the filter 
! -------------------
!
      DO jz=1,jpzsize
         jy = jz
!
         DO jdta=1,dtaend
            jbub = pt1bubidx(jz,jdta)
            jputab(jdta) = COUNT(bubblk1(:,:,:,:,jbub,1).NE.FREAL(0.0))
         ENDDO
         jpu = SUM(jputab(:))
!
         IF ((jpu.NE.0).OR.(kflagyo.NE.2)) THEN
!
            vectybub(0:) = FREAL(0.0)
!
            DO jdta=1,dtaend
               IF (jputab(jdta).NE.0) THEN
                  IF (larginptzon) THEN
                     jbub = pt4bubidx(jz,jdta)
                     ptbub(:,:,:,:,jdta) = NINT(bubblk4(:,:,:,:,jbub,1))
                  ELSE
                     CALL mkptbub (ptbub(:,:,:,:,jdta),jdta, &
     &                    pt1dtalon(jz,jdta),pt1dtalat(jz,jdta), &
     &                    pt1dtadepth(jz,jdta),pt1dtatime(jz,jdta), &
     &                    pt1bublon(jz,jdta),pt1bublat(jz,jdta), &
     &                    pt1bubdepth(jz,jdta),pt1bubtime(jz,jdta))
                  ENDIF
                  jbub = pt1bubidx(jz,jdta)
                  DO jt=1,zon_jpt
                  DO jk=1,zon_jpk
                  DO jj=1,zon_jpj
                  DO ji=1,zon_jpi
                     vectybub(ptbub(ji,jj,jk,jt,jdta)) =  &
     &                    bubblk1(ji,jj,jk,jt,jbub,1)
                  ENDDO
                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF
            ENDDO
!
            IF (largconnect) CALL mkconnecty(vectybub(0:jpysize))
!
            IF (kflagyo.EQ.2) THEN
               IF (jpysize.NE.jpssize) GOTO 1000
               vectsfilt(1:jpysize) = vectybub(1:jpysize)
            ELSE
               CALL mkhytoo(vectybub(1:jpysize), &
     &              vectsfilt(1:jpssize),poscoefobs(:,:))
            ENDIF
!
            jpu = COUNT(vectsfilt(1:jpssize) .NE. FREAL(0.0))
            ptu(1:jpu) = PACK( (/ (js, js=1,jpssize) /) ,  &
     &           vectsfilt(1:jpssize) .NE. FREAL(0.0) )
!     
            IF (jpu.GT.0) THEN
               vecty(jy) = DOT_PRODUCT(vects(ptu(1:jpu)), &
     &              vectsfilt(ptu(1:jpu))) &
     &              / SUM(vectsfilt(ptu(1:jpu)))
            ELSE
               vecty(jy) = spvaldta
            ENDIF
         ELSE
            vecty(jy) = vects(jy)
         ENDIF
!
         lmodprint=(MOD(jz-1,(jpzsize/10+1)).EQ.0)
         IF ((nprint.GE.1).AND.((lmodprint))) print *, &
     &        'Zone number : ',jz,'/',jpzsize,'( Size =',jpu,')'
!
      ENDDO
!
! -5.- Write output dta file
! --------------------------
!
      CALL writedta(kouty,vecty)
!
      IF (allocated(pt1bubidx)) deallocate (pt1bubidx)
      IF (allocated(pt1bublon)) deallocate (pt1bublon)
      IF (allocated(pt1bublat)) deallocate (pt1bublat)
      IF (allocated(pt1bubdepth)) deallocate (pt1bubdepth)
      IF (allocated(pt1bubtime)) deallocate (pt1bubtime)
      IF (allocated(pt1dtalon)) deallocate (pt1dtalon)
      IF (allocated(pt1dtalat)) deallocate (pt1dtalat)
      IF (allocated(pt1dtadepth)) deallocate (pt1dtadepth)
      IF (allocated(pt1dtatime)) deallocate (pt1dtatime)
!
      IF (allocated(vectptbub)) deallocate (vectptbub)
      IF (allocated(vectybub)) deallocate (vectybub)
      IF (allocated(ptbub)) deallocate (ptbub)
      IF (allocated(bubblk1)) deallocate (bubblk1)
      IF (allocated(ptu)) deallocate (ptu)
!
      IF (allocated(poscoefobs)) deallocate (poscoefobs)
      IF (allocated(vects)) deallocate (vects)
      IF (allocated(vectsfilt)) deallocate (vectsfilt)
      IF (allocated(vecty)) deallocate (vecty)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'algofilt','calcfilt')
 1001 CALL printerror2(0,1001,3,'algofilt','calcfilt')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algofilt
