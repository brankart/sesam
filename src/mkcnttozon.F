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
C ---                   MKCNTTOZON.F                            ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2000-02 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C ---
C --- SUBROUTINE mkcnttozon
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkcnttozon
      use mod_main
      use utilzone
      IMPLICIT NONE
      PRIVATE

      PUBLIC cnttozon

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE cnttozon(kargincfg,kargoutzon)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Generate influence bubbles from a set of contours
CCC  -------        and function types
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
      use mod_coord
      use mod_cont
      use mod_spacexyo , only : jpy, jpz
      use hiozon
      use hiocnt
      use hiogrd
      use utilmkto
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kargincfg,kargoutzon
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:,:), allocatable, save :: baseyr
C
      INTEGER :: allocok,jpysize
      INTEGER :: ji,jj,jk,jt,jy,jdta,jdta1,jc,jp,jnxyo,jlay,dimcnt,jbub,jz
      INTEGER :: dtaecnt1, flagxy
      INTEGER :: inddta1, inddta, inddtamsk, typcnt
      INTEGER, dimension(:), allocatable  :: dtaecnt
      CHARACTER(len=bgword) :: dtafcnt1
      CHARACTER(len=bgword), dimension(:), allocatable  :: dtafcnt
      CHARACTER(len=varlg) :: dta_nam1
      LOGICAL :: lmoyectold
      TYPE (type_gridij) :: grdpt
      BIGREAL :: levpt
      BIGREAL, dimension(4) :: jix,jjy
      BIGREAL :: r, s
      BIGREAL, dimension (:,:,:,:,:), allocatable :: zon_vect
      INTEGER, dimension (:,:), allocatable :: ptbubidx
      INTEGER, dimension (:,:), allocatable :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), allocatable :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), allocatable :: ptbublon, ptbublat
      INTEGER, dimension (:,:), allocatable :: ptbubdepth, ptbubtime
CC----------------------------------------------------------------------
C
      jpysize=jpy
C --- allocation dtaecnt
      allocate ( dtaecnt(1:varend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      dtaecnt(:) = 0
C --- allocation dtafcnt
      allocate ( dtafcnt(1:varend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      dtafcnt(:) = ''
C ----------------------------------------------------------------------
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&           routine mkcnttozon             &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
C
C -1.- Initialisation
C -------------------
C
      CALL openfile(10,kargincfg)
C
C --- read grid description
C
      DO jdta1 = 1,dtaend
         READ(10,*) dta_nam1, dtafcnt1, dtaecnt1
         DO jdta = 1,dtaend
            inddta = dta_ord(jdta)
            IF (dta_nam1(1:lenv(dta_nam1))
     $          .EQ.dta_nam(inddta)(1:lenv(dta_nam(inddta))) ) THEN
               dtafcnt(inddta) = dtafcnt1
               dtaecnt(inddta) = dtaecnt1
            ENDIF
         ENDDO
      ENDDO
      CLOSE(10)
C
C -2.- Determine the total number of zones
C ----------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKCNTTOZON : ',
     $         'Determining the total number of zones...'
      ENDIF
C
C --- Loop on dta variable fields
C
      jpz = 0
      DO jdta = 1,dtaend
         inddta = dta_ord(jdta)
C
C --- allocate and read contour file
C
         CALL evalhdrcnt(dtafcnt(inddta))
C
C --- allocation jpp
         allocate ( jpp(1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         jpp(:) = 0
C --- allocation contij
         allocate ( contij(1:jppend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contij(:,:)%longi = FREAL(0.0)
         contij(:,:)%latj = FREAL(0.0)
C --- allocation jplay
         allocate ( jplay(1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         jplay(:) = 0
C --- allocation contlevmin
         allocate ( contlevmin(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contlevmin(:,:) = FREAL(0.0)
C --- allocation contlevmax
         allocate ( contlevmax(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contlevmax(:,:) = FREAL(0.0)
C --- allocation contidx
         allocate ( contidx(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contidx(:,:) = 0
C
         typcnt = 0
         CALL readcnt(dtafcnt(inddta),typcnt)
C
         jpz = MAX(jpz,ABS(MAXVAL(contidx(:,:))))
C
         IF (allocated(contij)) deallocate(contij)
         IF (allocated(contlevmin)) deallocate(contlevmin)
         IF (allocated(contlevmax)) deallocate(contlevmax)
         IF (allocated(jpp)) deallocate(jpp)
         IF (allocated(jplay)) deallocate(jplay)
         IF (allocated(contidx)) deallocate(contidx)
C
      ENDDO
C
C --- allocation baseyr
      allocate ( baseyr(1:jpysize,1:jpz), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      baseyr(:,:) = FREAL(0.0)
C
C -2.- Generate the output Vy vector
C ----------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKCNTTOZON : ',
     $         'Generating the Vz object...'
      ENDIF
C
C --- Loop on dta variable fields
C
      jy = 0
      DO jdta = 1,dtaend
         inddta = dta_ord(jdta)
         inddtamsk = jdta - 1 + varend
C
C --- allocate and read contour file
C
         CALL evalhdrcnt(dtafcnt(inddta))
C
C --- allocation jpp
         allocate ( jpp(1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         jpp(:) = 0
C --- allocation contij
         allocate ( contij(1:jppend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contij(:,:)%longi = FREAL(0.0)
         contij(:,:)%latj = FREAL(0.0)
C --- allocation jplay
         allocate ( jplay(1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         jplay(:) = 0
C --- allocation contlevmin
         allocate ( contlevmin(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contlevmin(:,:) = FREAL(0.0)
C --- allocation contlevmax
         allocate ( contlevmax(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contlevmax(:,:) = FREAL(0.0)
C --- allocation contidx
         allocate ( contidx(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contidx(:,:) = 0
C --- allocation contf0
         allocate ( contf0(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contf0(:,:) = FREAL(0.0)
C --- allocation contfx
         allocate ( contfx(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contfx(:,:) = 0
C --- allocation contfy
         allocate ( contfy(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contfy(:,:) = 0
C --- allocation contfz
         allocate ( contfz(1:jplayend,1:jpc), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         contfz(:,:) = 0
C
         typcnt = 1
         CALL readcnt(dtafcnt(inddta),typcnt)
C
         IF (MOD(dtaecnt(inddta),10).EQ.0) THEN
C
C --- allocate and read grid file
C
            IF (dtangrd(inddta).LE.2) THEN
C --- allocation longi
               allocate ( longi(1:dta_jpi(inddta)) , stat = allocok )
               IF (allocok.GT.0) GOTO 1001
               longi(:) = FREAL(0.0)
C --- allocation latj
               allocate ( latj(1:dta_jpj(inddta)) , stat = allocok )
               IF (allocok.GT.0) GOTO 1001
               latj(:) = FREAL(0.0)
            ELSE
C --- allocation gridij
               allocate ( gridij(1:dta_jpi(inddta),1:dta_jpj(inddta)) ,
     $              stat = allocok )
               IF (allocok.GT.0) GOTO 1001
               gridij(:,:) = type_gridij(FREAL(0.0),FREAL(0.0))
            ENDIF
C
            flagxy=2
            CALL readgrd (flagxy,inddta)
C
         ENDIF
C
C --- Loop on grid indices
C
         DO jt=1,dta_jpt(inddta)
         DO jk=1,dta_jpk(inddta)
            DO jj=1,dta_jpj(inddta)
            DO ji=1,dta_jpi(inddta)
               IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
                  jy = jy + 1
                  IF (jy.GT.jpy) GOTO 1000
C
                  IF (MOD(dtaecnt(inddta),10).EQ.0) THEN
                     IF (dtangrd(inddta).LE.2) THEN
                        grdpt%longi = longi(ji)
                        grdpt%latj  = latj(jj)
                     ELSE
                        grdpt%longi = gridij(ji,jj)%longi
                        grdpt%latj  = gridij(ji,jj)%latj
                     ENDIF
                  ELSE
                     grdpt%longi = FREAL(ji)
                     grdpt%latj  = FREAL(jj)
                  ENDIF
C
                  IF (dtaecnt(inddta)/10.EQ.0) THEN
                     levpt = var_lev(jk,inddta)
                  ELSE
                     levpt = FREAL(jk)
                  ENDIF
C
                  DO jc = 1,jpc
                     IF (inarea(grdpt,jc)) THEN
                        DO jlay = 1,jplay(jc)
                           IF ( ( levpt.GE.contlevmin(jlay,jc) ) .AND.
     $                        ( levpt.LE.contlevmax(jlay,jc) ) ) THEN
C
                              baseyr(jy,contidx(jlay,jc)) = contf0(jlay,jc)
C
                              IF ( (contfx(jlay,jc).NE.0) .OR. 
     $                                (contfy(jlay,jc).NE.0) ) THEN
C
                                 IF (jpp(jc).NE.5) GOTO 102
C
                                 jix (1) = contij(1,jc)%longi
                                 jix (2) = contij(2,jc)%longi
                                 jix (3) = contij(3,jc)%longi
                                 jix (4) = contij(4,jc)%longi
                                 jjy (1) = contij(1,jc)%latj
                                 jjy (2) = contij(2,jc)%latj
                                 jjy (3) = contij(3,jc)%latj
                                 jjy (4) = contij(4,jc)%latj
C
                                 CALL mkxytors(grdpt%longi,grdpt%latj,jix,jjy,r,s)
C
                                 dimcnt = 1
                                 baseyr(jy,contidx(jlay,jc)) = baseyr(jy,contidx(jlay,jc))
     $                                     * fxyzcnt(r,dimcnt,jlay,jc)
                                 dimcnt = 2
                                 baseyr(jy,contidx(jlay,jc)) = baseyr(jy,contidx(jlay,jc))
     $                                     * fxyzcnt(s,dimcnt,jlay,jc)
C
                              ENDIF
C
                              IF ( (contfz(jlay,jc).NE.0) ) THEN
                                 IF ( contlevmax(jlay,jc) .NE. contlevmin(jlay,jc) ) THEN
C
                                    r = (levpt - contlevmin(jlay,jc))
                                    r = r / (contlevmax(jlay,jc) - contlevmin(jlay,jc))
C
                                    dimcnt = 3
                                    baseyr(jy,contidx(jlay,jc)) = baseyr(jy,contidx(jlay,jc))
     $                                        * fxyzcnt(r,dimcnt,jlay,jc)
C
                                 ENDIF
                              ENDIF
C
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
C
C --- deallocate grid and contour arrays
C
         IF (allocated(longi)) deallocate(longi)
         IF (allocated(latj)) deallocate(latj)
         IF (allocated(gridij)) deallocate(gridij)
         IF (allocated(contij)) deallocate(contij)
         IF (allocated(contlevmin)) deallocate(contlevmin)
         IF (allocated(contlevmax)) deallocate(contlevmax)
         IF (allocated(jpp)) deallocate(jpp)
         IF (allocated(jplay)) deallocate(jplay)
         IF (allocated(contidx)) deallocate(contidx)
         IF (allocated(contf0)) deallocate(contf0)
         IF (allocated(contfx)) deallocate(contfx)
         IF (allocated(contfy)) deallocate(contfy)
         IF (allocated(contfz)) deallocate(contfz)
C
      ENDDO
C
      IF (jy.NE.jpy) GOTO 1000
C
C -3.- Writing the Vz object
C --------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'MKCNTTOZON : ',
     $             'Writing the Vz object... '
      ENDIF
C
      zon_jpi = MAXVAL(dta_jpi(1:dtaend))
      zon_jpj = MAXVAL(dta_jpj(1:dtaend))
      zon_jpk = MAXVAL(dta_jpk(1:dtaend))
      zon_jpt = MAXVAL(dta_jpt(1:dtaend))
      jpbub = jpz * dtaend
C
      CALL writehdrzon(kargoutzon,zon_jpi,zon_jpj,zon_jpk,zon_jpt,jpbub,jpz)
C
C --- allocation zon_vect
      allocate ( zon_vect(1:zon_jpi,1:zon_jpj,1:zon_jpk,
     $                    1:zon_jpt,1:1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      zon_vect(:,:,:,:,:) = FREAL(0.0)
C --- allocation zone pointers
      allocate ( ptbubidx(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptbubidx(:,:) = 0
C
      jbub = 0
      DO jz = 1,jpz
         jy = 0
         DO jdta = 1,dtaend
            inddta=dta_ord(jdta)
            inddtamsk=jdta-1+varend
C
            DO jt = 1,dta_jpt(inddta)
            DO jk = 1,dta_jpk(inddta)
               DO jj = 1,dta_jpj(inddta)
               DO ji = 1,dta_jpi(inddta)
                  IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
                     jy = jy + 1
                     zon_vect(ji,jj,jk,jt,1) = baseyr(jy,jz)
                  ENDIF
               ENDDO
               ENDDO
            ENDDO
            ENDDO
C
            jbub = jbub + 1
            ptbubidx(jz,jdta) = jbub
            CALL writenbubzon (kargoutzon,zon_vect,jbub)
C
         ENDDO
         IF (jy.NE.jpy) GOTO 1000
      ENDDO
C
      IF (jbub.NE.jpbub) GOTO 1000
C
C --- desallocate arrays
      IF (allocated(baseyr)) deallocate(baseyr)
      IF (allocated(dtaecnt)) deallocate(dtaecnt)
      IF (allocated(dtafcnt)) deallocate(dtafcnt)
      IF (allocated(zon_vect)) deallocate(zon_vect)
C --- allocation zone pointers
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
C
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
C
      ptdtalon(:,:) = 1
      ptdtalat(:,:) = 1
      ptdtadepth(:,:) = 1
      ptdtatime(:,:) = 1
      ptbublon(:,:) = 1
      ptbublat(:,:) = 1
      ptbubdepth(:,:) = 1
      ptbubtime(:,:) = 1
C
      CALL writeptzon (kargoutzon,ptbubidx,
     $           ptdtalon,ptdtalat,ptdtadepth,ptdtatime,
     $           ptbublon,ptbublat,ptbubdepth,ptbubtime)
C
C --- deallocation
      IF (allocated(ptbubidx)) deallocate(ptbubidx)
      IF (allocated(ptdtalon)) deallocate(ptdtalon)
      IF (allocated(ptdtalat)) deallocate(ptdtalat)
      IF (allocated(ptdtadepth)) deallocate(ptdtadepth)
      IF (allocated(ptdtatime)) deallocate(ptdtatime)
      IF (allocated(ptbublon)) deallocate(ptbublon)
      IF (allocated(ptbublat)) deallocate(ptbublat)
      IF (allocated(ptbubdepth)) deallocate(ptbubdepth)
      IF (allocated(ptbubtime)) deallocate(ptbubtime)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mkcnttozon','cnttozon')
 1001 CALL printerror2(0,1001,3,'mkcnttozon','cnttozon')
C
 102  WRITE (texterror,*) 'contour should be quadrilateral'
      CALL printerror2(0,102,3,'mkcnttozon','cnttozon',
     $     comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkcnttozon
