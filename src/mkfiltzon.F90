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
! ---                   MKFILTZON.F90                             ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2000-05 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! ---
! --- SUBROUTINE filtzon
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkfiltzon
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC filtzon

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE filtzon(kargincfg,kargoutzon,kargindta)
!---------------------------------------------------------------------
!
!  Purpose : Generate a filter convolution matrix
!  -------        from an  isotropic and homogeneous
!                 filter spectral function
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
      use mod_mask
      use mod_cont
      use mod_spacexyo , only : jpy, jpyend, jpz
      use hioxyo
      use hiozon
      use utilfilt
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kargoutzon,kargincfg
      CHARACTER(len=*), intent(in), optional :: kargindta
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vectylength
      BIGREAL, dimension(:), allocatable, save :: vectbublength
!
      INTEGER :: allocok, numfila
      INTEGER :: jdta, jdta1, ji, jj, jk, jt, jy, jbub, jpysize, jnxyo
      INTEGER :: inddta, inddtamsk, flagxyo
      INTEGER :: zon_li,zon_lj
      INTEGER :: deli, delj
      BIGREAL :: delr,filt_lengthold
      INTEGER, dimension (:,:), allocatable :: ptbubidx
      INTEGER, dimension (:,:), allocatable :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), allocatable :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), allocatable :: ptbublon, ptbublat
      INTEGER, dimension (:,:), allocatable :: ptbubdepth, ptbubtime
      LOGICAL :: lmoyectold, homogeneous, lectinfo, lmodprint
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mkfiltzon &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Initialisation
! -------------------
!
      jpysize = jpy
      lectinfo = .FALSE.
      flagxyo = 2
!
      homogeneous = .NOT. present(kargindta)
!
      numfila=10
      CALL openfile(numfila,kargincfg)
      READ(numfila,*) zon_li, zon_lj
      READ(numfila,*) filt_typ
      SELECT CASE (filt_typ(1:lenv(filt_typ)))
      CASE ('spline')
         READ(numfila,*) filt_ord, filt_length, filt_threshold
! --- allocation filt_coef
         allocate ( filt_coef(0:filt_ord), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         READ(numfila,*) filt_coef(0:filt_ord)
      CASE('constant')
!
         filt_length=FREAL(1000000000000.0)
      CASE DEFAULT
         GOTO 102
      END SELECT
      CLOSE(numfila)
!
! -2.- Read the non-homogeneous filter length
! -------------------------------------------
!
      IF (.NOT.homogeneous) THEN
!
! --- allocation vectylength
         allocate ( vectylength(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectylength(:) = FREAL(0.0)
!
         jnxyo=1
         lmoyectold=lmoyect
         lmoyect=.FALSE.
         CALL readxyo(kargindta,vectylength(:), &
     &                   jnxyo,lectinfo,flagxyo)
         lmoyect=lmoyectold
!
         filt_lengthold=vectylength(1)
         jpbub=1
         DO jy = 2, jpyend
            IF (vectylength(jy).NE.filt_lengthold) THEN
               filt_lengthold=vectylength(jy)
               jpbub=jpbub+1
            ENDIF
         ENDDO
!
! --- allocation vectbublength
         allocate ( vectbublength(1:jpbub), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectbublength(:) = FREAL(0.0)
!
         jbub=1
         vectbublength(jbub)=vectylength(1)
         vectylength(1)=FREAL(jbub)
         DO jy = 2, jpyend
            IF (vectylength(jy).NE.vectbublength(jbub)) THEN
               jbub=jbub+1
               vectbublength(jbub)=vectylength(jy)
            ENDIF
            vectylength(jy)=FREAL(jbub)
         ENDDO
!
      ENDIF
!
! -3.- Write the zone file header
! -------------------------------
!
      zon_jpi = 1 + 2 * zon_li
      zon_jpj = 1 + 2 * zon_lj
      zon_jpk = 1
      zon_jpt = 1
      jpz = jpysize
!
      IF (homogeneous) THEN
         jpbub = 2
      ELSE
         jpbub = jpbub + 1
      ENDIF
!
      lmoyectold=lmoyect
      lmoyect=.FALSE.
      CALL writehdrzon(kargoutzon,zon_jpi,zon_jpj,zon_jpk,zon_jpt, &
     &     jpbub,jpz)
      lmoyect=lmoyectold
!
! -4.- Generate and write the filter convolution matrix
! -----------------------------------------------------
!
! --- allocation filt_bub
      allocate ( filt_bub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,1:1), &
     &     stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      filt_bub(:,:,:,:,:) = FREAL(0.0)
!
      DO jbub = 1,jpbub-1
!
         IF (.NOT.homogeneous) THEN
!
            lmodprint=(MOD(jbub-1,((jpbub-1)/10+1)).EQ.0)
            IF ((nprint.GE.1).AND.((lmodprint))) print *, &
     &           'Bubble index : ',jbub,'/',jpbub
!
            filt_length = vectbublength(jbub)
         ENDIF
         IF (filt_length.GT.FREAL(0.0)) THEN
!
            DO jt = 1,zon_jpt
            DO jk = 1,zon_jpk
            DO jj = 1,zon_jpj
            DO ji = 1,zon_jpi
               deli = ABS(ji - zon_li - 1)
               delj = ABS(jj - zon_lj - 1)
               delr = SQRT( FREAL( deli * deli + delj * delj ) )
               filt_bub(ji,jj,jk,jt,1) = calcfilt(delr)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
!
            filt_bub(:,:,:,:,1) =  &
     &           filt_bub(:,:,:,:,1) / filt_bub(zon_li+1,zon_lj+1,1,1,1)
!
            WHERE ( filt_bub(:,:,:,:,1) .LT. filt_threshold )  &
     &           filt_bub(:,:,:,:,1) = FREAL(0.0)
         ELSE
!
            IF (filt_length.EQ.FREAL(0.0)) THEN
               filt_bub(:,:,:,:,1) = FREAL(0.0)
               filt_bub(1 + zon_li,1 + zon_lj,1,1,1) = FREAL(1.0)
            ELSE
               filt_bub(:,:,:,:,1) = FREAL(0.0)
            ENDIF
!
         ENDIF
!
         CALL writenbubzon(kargoutzon,filt_bub,jbub)
!
      ENDDO
!
! --- Write the empty bubble at the last position
!
      IF (homogeneous) THEN
         jbub = 2
      ELSE
         jbub = jpbub
      ENDIF
!
      filt_bub(:,:,:,:,1) = FREAL(0.0)
      CALL writenbubzon(kargoutzon,filt_bub,jbub)
!
      IF (allocated(filt_bub))  deallocate (filt_bub) 
      IF (allocated(filt_coef))  deallocate (filt_coef) 
!
! -5.- Generate and write zone pointers
! -------------------------------------
!-
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
      ptbublon(:,:) = zon_li + 1
      ptbublat(:,:) = zon_lj + 1
      ptbubdepth(:,:) = 1
      ptbubtime(:,:) = 1
!
      jy = 0
      DO jdta=1,dtaend
         inddta = dta_ord(jdta)
         inddtamsk = jdta - 1 + varend
!
         DO jt=1,dta_jpt(inddta)
         DO jk=1,dta_jpk(inddta)
            DO jj=1,dta_jpj(inddta)
            DO ji=1,dta_jpi(inddta)
               IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
                  jy=jy+1
                  DO jdta1 = 1,dtaend
                     ptdtalon(jy,jdta1) = ji
                     ptdtalat(jy,jdta1) = jj
                     ptdtadepth(jy,jdta1) = jk
                     ptdtatime(jy,jdta1) = jt
!
                     IF (homogeneous) THEN
                        IF (jdta1.EQ.jdta) THEN
                           ptbubidx(jy,jdta1) = 1
                        ELSE
                           ptbubidx(jy,jdta1) = 2
                        ENDIF
                     ELSE
                        IF (jdta1.EQ.jdta) THEN
                           ptbubidx(jy,jdta1) = NINT(vectylength(jy))
                        ELSE
                           ptbubidx(jy,jdta1) = jpbub
                        ENDIF
                     ENDIF
!
                  ENDDO
               ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
!
      ENDDO
!
      CALL writeptzon (kargoutzon,ptbubidx, &
     &           ptdtalon,ptdtalat,ptdtadepth,ptdtatime, &
     &           ptbublon,ptbublat,ptbubdepth,ptbubtime)
!
      IF (allocated(ptbubidx)) deallocate (ptbubidx)
      IF (allocated(ptbublon)) deallocate (ptbublon)
      IF (allocated(ptbublat)) deallocate (ptbublat)
      IF (allocated(ptbubdepth)) deallocate (ptbubdepth)
      IF (allocated(ptbubtime)) deallocate (ptbubtime)
      IF (allocated(ptdtalon)) deallocate (ptdtalon)
      IF (allocated(ptdtalat)) deallocate (ptdtalat)
      IF (allocated(ptdtadepth)) deallocate (ptdtadepth)
      IF (allocated(ptdtatime)) deallocate (ptdtatime)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkfiltzon','filtzon')
 1001 CALL printerror2(0,1001,3,'mkfiltzon','filtzon')
!
 102  WRITE (texterror,*) 'incorrect filter type : ', &
     &                       filt_typ(1:lenv(filt_typ))
      CALL printerror2(0,102,3,'mkfiltzon','filtzon', &
     &      comment=texterror)
 103  WRITE (texterror,*) 'non homogeneous filter only ', &
     &                       'available for the spline type'
      CALL printerror2(0,103,3,'mkfiltzon','filtzon', &
     &      comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkfiltzon
