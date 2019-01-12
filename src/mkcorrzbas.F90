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
! ---                   MKCORRZBAS.F90                            ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2000-01 (J.M. Brankart)                    ---
! --- modification : 2000-03 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! ---
! --- SUBROUTINE mkcorrzbas
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkcorrzbas
      use mod_main
      use utilzone
      IMPLICIT NONE
      PRIVATE

      PUBLIC corrzbas

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE corrzbas(karginzon,kargincfg,kargoutzbas)
!---------------------------------------------------------------------
!
!  Purpose : Generate a local error correlation matrix 
!  -------        from an  isotropic and homogeneous 
!                 analytical correlation function
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
      use mod_spacexyo , only : jprend, jpz
      use hiozon
      use hiobas
      use utilroa
      use utilvalid
      use utilprint
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: karginzon,kargincfg, &
     &                     kargoutzbas
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jpbsize,jprsize,jbub,nvpnull,numjr,serie
      INTEGER :: valbase,jprbasout,jdta,jdta1
      INTEGER, dimension (:,:), allocatable :: ptbubidx
      INTEGER, dimension (:,:), allocatable :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), allocatable :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), allocatable :: ptbublon, ptbublat
      INTEGER, dimension (:,:), allocatable :: ptbubdepth, ptbubtime
      BIGREAL, dimension (:,:,:,:,:), allocatable :: bub
      BIGREAL, dimension (:,:), allocatable :: matbb, matUbb
      BIGREAL, dimension (:), allocatable :: lambda
      INTEGER :: jil,jjl,jkl,jtl,jic,jjc,jkc,jtc,jbl,jbc,jb
      INTEGER :: jrbasdeb, jrbasfin, jrmatdeb, jrmatfin, jrmat, jr
      BIGREAL :: clgt_i, clgt_j, clgt_k, clgt_t
      INTEGER :: ctyp_i, ctyp_j, ctyp_k, ctyp_t
      CHARACTER(len=bgword) :: titre,variable,fnameout,vctnamout
      LOGICAL :: inverse,normalize
!----------------------------------------------------------------------
!
      jpbsize=zon_jpj*zon_jpj*zon_jpk*zon_jpt
      jprsize=jprend
      IF (jpbsize*dtaend.LT.jprsize-1) GOTO 102
      IF (MOD(jprsize-1,dtaend).NE.0) GOTO 103
! --- allocation matbb
      allocate ( matbb(1:jpbsize,1:jpbsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      matbb(:,:) = FREAL(0.0)
! --- allocation matUbb
      allocate ( matUbb(1:jpbsize,1:jpbsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      matUbb(:,:) = FREAL(0.0)
! --- allocation lambda
      allocate ( lambda(1:jpbsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      lambda(:) = FREAL(0.0)
! --- allocation bub
      allocate ( bub(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,1:1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      bub(:,:,:,:,:) = FREAL(0.0)
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
! ----------------------------------------------------------------------
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&           routine mkcorrzbas             &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
!
! -1.- Initialisation
! -------------------
!
      jrbasdeb = 2
      jrbasfin = jprsize
      jrmatdeb = jrbasdeb - 1
      jrmatfin = jrbasfin - 1
!
      CALL openfile(10,kargincfg)
      READ(10,*) inverse,normalize
      READ(10,*) clgt_i, clgt_j, clgt_k, clgt_t
      READ(10,*) ctyp_i, ctyp_j, ctyp_k, ctyp_t
      CLOSE(10)
!
! -2.- Transfer zone configuration to zbas
! ----------------------------------------
!
      jpbub = dtaend
      CALL writehdrzbas(kargoutzbas,zon_jpi,zon_jpj,zon_jpk, &
     &           zon_jpt,jpbub,jpz,jrbasdeb,jrbasfin)
      CALL readptzon (karginzon,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime)
!
      DO jdta = 1,dtaend
         ptbubidx(:,jdta) = jdta
      ENDDO
!
      CALL writeptzbas(kargoutzbas,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime,jrbasdeb,jrbasfin)
!
      jr = 0
      CALL writehdrzbas(kargoutzbas,zon_jpi,zon_jpj,zon_jpk, &
     &           zon_jpt,jpbub,jpz,jr,jr)
      CALL writeptzbas(kargoutzbas,ptbubidx,ptdtalon,ptdtalat, &
     &           ptdtadepth,ptdtatime,ptbublon,ptbublat,ptbubdepth, &
     &           ptbubtime,jr,jr)
!
! -3.- Compute correlation matrix: matbb
! --------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKCORRZBAS : ', &
     &         'Computing the correlation matrix'
      ENDIF
!
      jbl = 0
      DO jtl = 1,zon_jpt
         DO jkl = 1,zon_jpk
            DO jjl = 1,zon_jpj
               DO jil = 1,zon_jpi
                  jbl = jbl + 1
                  jbc = 0
                  DO jtc = 1,zon_jpt
                     DO jkc = 1,zon_jpk
                        DO jjc = 1,zon_jpj
                           DO jic = 1,zon_jpi
                              jbc = jbc + 1
                              IF (jbc.GE.jbl) THEN
                                 matbb(jbl,jbc) = fctcorrel &
     &                            (jil,jjl,jkl,jtl,jic,jjc,jkc,jtc, &
     &                             clgt_i,clgt_j,clgt_k,clgt_t, &
     &                             ctyp_i,ctyp_j,ctyp_k,ctyp_t)
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                  matbb(jbl:jpbsize,jbl)=matbb(jbl,jbl:jpbsize)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
! -4.- Compute SVD decomposition of matbb
! ---------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKCORRZBAS : ', &
     &     'Computing SVD decomposition of the correlation matrix'
      ENDIF
!
      nvpnull=0
      CALL mkrs (matbb(1:jpbsize,1:jpbsize), &
     &     matUbb(1:jpbsize,1:jpbsize), &
     &     lambda(1:jpbsize),nvpnull)
!
      IF (nprint.GE.2) THEN
         titre='MKCORRZBAS: eigenvalues of correlation matrix'
         variable='SQRT(lambda)'
         CALL printtab_r (SQRT(lambda(jrmatdeb:jrmatfin)), &
     &              titre,variable,jrmatdeb)
      ENDIF
!
! -5.- Write diagonal of reduced order correlation matrix
! -------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKCORRZBAS : ', &
     &             'Writing the diagonal of the reduced ', &
     &             'rank correlation matrix'
      ENDIF
!
      jbub = 1
      bub(:,:,:,:,jbub) = FREAL(0.0)
      DO jrmat = 1,(jprsize-1)/dtaend
         jbl = 0
         DO jtl = 1,zon_jpt
            DO jkl = 1,zon_jpk
               DO jjl = 1,zon_jpj
                  DO jil = 1,zon_jpi
                     jbl = jbl + 1
                     bub(jil,jjl,jkl,jtl,jbub) = &
     &                      bub(jil,jjl,jkl,jtl,jbub) + &
     &                      matUbb(jbl,jrmat) * matUbb(jbl,jrmat) * lambda(jrmat)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      bub(:,:,:,:,jbub) = SQRT(bub(:,:,:,:,jbub))
!
      numjr=0
      serie=1
      CALL fildirbas (vctnamout,kargoutzbas,jprbasout,numjr,serie)
      WRITE(fnameout,'("./",A,"/",A)') kargoutzbas(1:lenv(kargoutzbas)), &
     &           vctnamout(1:lenv(vctnamout))
!
      DO jdta = 1,dtaend
         CALL writenbubzon(fnameout,bub(:,:,:,:,jbub:jbub),jdta)
      ENDDO
!
! -6.- Compute the modes of the correlation matrix : GAMMA = THETA transp(THETA)
! ------------------------------------------------------------------------------
!
      lambda(1:jpbsize) = SQRT(lambda(1:jpbsize))
!
      IF (inverse) THEN
         DO jb = 1,jpbsize
            matUbb(1:jpbsize,jb) = matUbb(1:jpbsize,jb) / lambda(jb)
         ENDDO
      ELSE
         DO jb = 1,jpbsize
            matUbb(1:jpbsize,jb) = matUbb(1:jpbsize,jb) * lambda(jb)
         ENDDO
      ENDIF
!
      IF (normalize) THEN
         IF (inverse) THEN
            jbub = 1
            DO jrmat = 1,(jprsize-1)/dtaend
               jbl = 0
               DO jtl = 1,zon_jpt
               DO jkl = 1,zon_jpk
                  DO jjl = 1,zon_jpj
                  DO jil = 1,zon_jpi
                     jbl = jbl + 1
                     matUbb(jbl,jrmat) = matUbb(jbl,jrmat) * bub(jil,jjl,jkl,jtl,jbub)
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
            ENDDO
         ELSE
            jbub = 1
            DO jrmat = 1,(jprsize-1)/dtaend
               jbl = 0
               DO jtl = 1,zon_jpt
               DO jkl = 1,zon_jpk
                  DO jjl = 1,zon_jpj
                  DO jil = 1,zon_jpi
                     jbl = jbl + 1
                     matUbb(jbl,jrmat) = matUbb(jbl,jrmat) / bub(jil,jjl,jkl,jtl,jbub)
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!
! -7.- Write output modes
! -----------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKCORRZBAS : ', &
     &     'Writing the reduced rank correlation matrix'
      ENDIF
!
      valbase = 4
      CALL writeinfobas(kargoutzbas,valbase)
!
      jbub = 1
      jr = jrbasdeb - 1
      DO jrmat = 1,(jprsize-1)/dtaend
         DO jdta = 1,dtaend
            jr = jr + 1
!
            numjr=jr
            serie=1
            CALL fildirbas (vctnamout,kargoutzbas,jprbasout,numjr,serie)
            WRITE(fnameout,'("./",A,"/",A)')  &
     &           kargoutzbas(1:lenv(kargoutzbas)), &
     &              vctnamout(1:lenv(vctnamout))
!
            DO jdta1 = 1,dtaend
               IF (jdta1.EQ.jdta) THEN
                  jbl = 0
                  DO jtl = 1,zon_jpt
                     DO jkl = 1,zon_jpk
                        DO jjl = 1,zon_jpj
                           DO jil = 1,zon_jpi
                              jbl = jbl + 1
                              bub(jil,jjl,jkl,jtl,jbub) = matUbb(jbl,jrmat)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ELSE
                  bub(:,:,:,:,jbub) = FREAL(0.0)
               ENDIF
!
               CALL writenbubzon(fnameout,bub(:,:,:,:,jbub:jbub),jdta1)
            ENDDO
!
         ENDDO
      ENDDO
!
! --- desallocate arrays
      IF (allocated(matbb)) deallocate(matbb)
      IF (allocated(matUbb)) deallocate(matUbb)
      IF (allocated(bub)) deallocate(bub)
!
      IF (allocated(ptbubidx)) deallocate(ptbubidx)
      IF (allocated(ptdtalon)) deallocate(ptdtalon)
      IF (allocated(ptdtalat)) deallocate(ptdtalat)
      IF (allocated(ptdtadepth)) deallocate(ptdtadepth)
      IF (allocated(ptdtatime)) deallocate(ptdtatime)
      IF (allocated(ptbublon)) deallocate(ptbublon)
      IF (allocated(ptbublat)) deallocate(ptbublat)
      IF (allocated(ptbubdepth)) deallocate(ptbubdepth)
      IF (allocated(ptbubtime)) deallocate(ptbubtime)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkcorrzbas','corrzbas')
 1001 CALL printerror2(0,1001,3,'mkcorrzbas','corrzbas')
!
 102  WRITE (texterror,*) 'base rank larger than bubble size'
      CALL printerror2(0,102,3,'mkcorrzbas','corrzbas', &
     &     comment=texterror)
 103  WRITE (texterror,*) 'base rank must be a multiple of dtaend'
      CALL printerror2(0,103,3,'mkcorrzbas','corrzbas', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkcorrzbas
