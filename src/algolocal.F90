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
! ---                   ALGOLOCAL.F90                           ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2021-06 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! ---
! --- SUBROUTINE calclocal
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algolocal
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calclocal

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calclocal(kargoutpartvar,kargoutzon,kargincfg)
!---------------------------------------------------------------------
!
!  Purpose : Generate a partition of the domain and the corresponding influence bubbles
!  -------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : jpx, jpz
      use hioxyo
      use hiozon
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kargoutpartvar,kargoutzon, &
     &     kargincfg
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vectx
!
      INTEGER :: jpxsize
      INTEGER :: subdom_size_i,dom_size_i,part_nbr_i,zon_li
      INTEGER :: subdom_size_j,dom_size_j,part_nbr_j,zon_lj
      INTEGER :: subdom_size_k,dom_size_k,part_nbr_k,zon_lk
      INTEGER :: subdom_size_t,dom_size_t,part_nbr_t,zon_lt
      INTEGER, dimension (:,:,:,:), allocatable :: part_index
      INTEGER, dimension (:), allocatable :: pini_i,pini_j,pini_k,pini_t
      INTEGER :: ji,jj,jk,jt,jx,jz,jvar,jbub,jdta,jnxyo,jzidx
      INTEGER :: jimin,jimax,jjmin,jjmax,jkmin,jkmax,jtmin,jtmax
      INTEGER :: indvarmsk,indvar
      BIGREAL, dimension (:,:,:,:), allocatable :: part_vect
      BIGREAL, dimension (:,:,:,:,:), allocatable :: zon_vect
      INTEGER :: zon_lcut_i,zon_lcut_j,zon_lcut_k,zon_lcut_t
      BIGREAL :: zon_lw_i,zon_lw_j,zon_lw_k,zon_lw_t
      INTEGER, dimension (:,:), allocatable :: ptbubidx
      INTEGER, dimension (:,:), allocatable :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), allocatable :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), allocatable :: ptbublon, ptbublat
      INTEGER, dimension (:,:), allocatable :: ptbubdepth, ptbubtime
      LOGICAL :: inmask,lmoyectold,lectinfo,ingrid
      INTEGER ::  allocok,numfila,flagxyo
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
      BIGREAL :: zjic, zjis, zjid, zjjc, zjjs, zjjd
      BIGREAL :: zjkc, zjks, zjkd, zjtc, zjts, zjtd
      BIGREAL :: shifti, shiftj, shiftk, shiftt
      BIGREAL :: expi, expj, expk, expt
!----------------------------------------------------------------------
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine calclocal &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
! ----------------------------------------------------------------------
      jpxsize=jpx
! --- allocation vectx
      allocate ( vectx(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectx(:) = FREAL(0.0)
!
! -0.- Initialisation
! -------------------
!
      textexclusion='#'
      numfila=10
      CALL openfile(numfila,kargincfg)
!
      text=readnextline(numfila,textexclusion)
      READ(text,*) subdom_size_i, subdom_size_j, subdom_size_k, subdom_size_t
      text=readnextline(numfila,textexclusion)
      READ(text,*) zon_lcut_i, zon_lcut_j, zon_lcut_k, zon_lcut_t
      text=readnextline(numfila,textexclusion)
      READ(text,*) zon_lw_i, zon_lw_j, zon_lw_k, zon_lw_t
!
      CLOSE(numfila)
!
      IF (subdom_size_i==0) zon_lcut_i=0
      IF (subdom_size_j==0) zon_lcut_j=0
      IF (subdom_size_k==0) zon_lcut_k=0
      IF (subdom_size_t==0) zon_lcut_t=0
!
      jpbub=1
!
! -1.- Give subdomain index (jz) to every subdomain (ji,jj)
! ---------------------------------------------------------
! Compute domain size for coarse grid
      dom_size_i = 0
      dom_size_j = 0
      dom_size_k = 0
      dom_size_t = 0
      DO jvar = 1,varend
        indvar = var_ord(jvar)
        IF (var_jpi(indvar).GT.dom_size_i) dom_size_i = var_jpi(indvar)
        IF (var_jpj(indvar).GT.dom_size_j) dom_size_j = var_jpj(indvar)
        IF (var_jpk(indvar).GT.dom_size_k) dom_size_k = var_jpk(indvar)
        IF (var_jpt(indvar).GT.dom_size_t) dom_size_t = var_jpt(indvar)
      ENDDO     
      IF (subdom_size_i==0) subdom_size_i=dom_size_i
      IF (subdom_size_j==0) subdom_size_j=dom_size_j
      IF (subdom_size_k==0) subdom_size_k=dom_size_k
      IF (subdom_size_t==0) subdom_size_t=dom_size_t
!
! Compute number of subdomains to cover the full grid
      part_nbr_i = ( dom_size_i - 1 ) / subdom_size_i + 1
      part_nbr_j = ( dom_size_j - 1 ) / subdom_size_j + 1
      part_nbr_k = ( dom_size_k - 1 ) / subdom_size_k + 1
      part_nbr_t = ( dom_size_t - 1 ) / subdom_size_t + 1

! --- allocation part_index
      allocate ( part_index(1:part_nbr_i,1:part_nbr_j, &
                     &      1:part_nbr_k,1:part_nbr_t), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      part_index(:,:,:,:) = 0
!
! Compute subdomain index (jz) for each subdomain (ji,jj)
      jz = 0
      DO jt = 1,part_nbr_t
      DO jk = 1,part_nbr_k
      DO jj = 1,part_nbr_j
      DO ji = 1,part_nbr_i
        jz = jz + 1
        part_index(ji,jj,jk,jt) = jz
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      jpz = jz
!
! -2.- Store initial grid point (pini_i,pini_j) of each subdomain (jz)
! --------------------------------------------------------------------
! --- allocation pini_i
      allocate ( pini_i(1:jpz), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pini_i(:) = 0
! --- allocation pini_j
      allocate ( pini_j(1:jpz), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pini_j(:) = 0
! --- allocation pini_k
      allocate ( pini_k(1:jpz), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pini_k(:) = 0
! --- allocation pini_t
      allocate ( pini_t(1:jpz), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pini_t(:) = 0
!
      DO ji = 1,part_nbr_i
      DO jj = 1,part_nbr_j
      DO jk = 1,part_nbr_k
      DO jt = 1,part_nbr_t
        pini_i(part_index(ji,jj,jk,jt)) = 1 + (ji-1) * subdom_size_i
        pini_j(part_index(ji,jj,jk,jt)) = 1 + (jj-1) * subdom_size_j
        pini_k(part_index(ji,jj,jk,jt)) = 1 + (jk-1) * subdom_size_k
        pini_t(part_index(ji,jj,jk,jt)) = 1 + (jt-1) * subdom_size_t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!
! --- desallocation of arrays
      IF (allocated(part_index)) deallocate (part_index)
!
! -3.2- Eliminate subdomains outside model mask
! ---------------------------------------------
!
      jz = 1
      DO WHILE (jz .LE. jpz)
!
        inmask = .FALSE.
        DO jvar = 1,varend
          indvarmsk = jvar - 1
          indvar = var_ord(jvar)
!
! Set the boundaries of the subdomain
          jimin=pini_i(jz) 
          jjmin=pini_j(jz) 
          jkmin=pini_k(jz) 
          jtmin=pini_t(jz) 
          jimax=MIN( pini_i(jz)+subdom_size_i-1, var_jpi(indvar) )
          jjmax=MIN( pini_j(jz)+subdom_size_j-1, var_jpj(indvar) )
          jkmax=MIN( pini_k(jz)+subdom_size_k-1, var_jpk(indvar) )
          jtmax=MIN( pini_t(jz)+subdom_size_t-1, var_jpt(indvar) )
!
! Check the lower boundary of the subdomain
          ingrid = jimin.LE.var_jpi(indvar)
          ingrid = ingrid .AND. ( jjmin.LE.var_jpj(indvar) ) 
          ingrid = ingrid .AND. ( jkmin.LE.var_jpk(indvar) ) 
          ingrid = ingrid .AND. ( jtmin.LE.var_jpt(indvar) ) 
!
! Check if there is at least one value in mask
          IF (ingrid) THEN
            DO jt=  jtmin,jtmax
            DO jk = jkmin,jkmax
            DO jj = jjmin,jjmax
            DO ji = jimin,jimax
              inmask = inmask .OR.  &
     &         (IBITS(mask(ji,jj,jk,jt),indvarmsk,1).NE.0) 
            ENDDO
            ENDDO
            ENDDO
            ENDDO
          ENDIF
!
        ENDDO
!
        IF (.NOT.inmask) THEN
          DO jzidx = jz+1,jpz
             pini_i(jzidx-1) = pini_i(jzidx)
             pini_j(jzidx-1) = pini_j(jzidx)
             pini_k(jzidx-1) = pini_k(jzidx)
             pini_t(jzidx-1) = pini_t(jzidx)
          ENDDO
          jpz = jpz - 1
        ELSE
          jz = jz + 1
        ENDIF
!
      ENDDO
!
! -4.- Generate partition vector
! ------------------------------
!
! --- allocation part_vect
      jimin=1 ; jimax=dom_size_i
      jjmin=1 ; jjmax=dom_size_j
      jkmin=1 ; jkmax=dom_size_k
      jtmin=1 ; jtmax=dom_size_t
      allocate(part_vect(jimin:jimax,jjmin:jjmax,jkmin:jkmax,jtmin:jtmax),stat=allocok)
      IF (allocok.NE.0) GOTO 1001
      part_vect(:,:,:,:) = FREAL(0.0)
!
      DO jz = 1,jpz
        DO jt=pini_t(jz),MIN(pini_t(jz)+subdom_size_t-1,jtmax)
        DO jk=pini_k(jz),MIN(pini_k(jz)+subdom_size_k-1,jkmax)
        DO jj=pini_j(jz),MIN(pini_j(jz)+subdom_size_j-1,jjmax)
        DO ji=pini_i(jz),MIN(pini_i(jz)+subdom_size_i-1,jimax)
          part_vect(ji,jj,jk,jt)=FREAL(jz)
        ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDDO
!
      jx = 0
      DO jvar = 1,varend
         indvar = var_ord(jvar)
         indvarmsk = jvar - 1
         DO jt=1,var_jpt(indvar)
         DO jk=1,var_jpk(indvar)
            DO jj=1,var_jpj(indvar)
            DO ji=1,var_jpi(indvar)
               IF (IBITS(mask(ji,jj,jk,jt),indvarmsk,1).NE.0) THEN
                  jx = jx + 1
                  vectx(jx) = part_vect(ji,jj,jk,jt)
               ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      IF (jx.NE.jpx) GOTO 1000
!
! -5.- Write partition vector
! ---------------------------
!
      lmoyectold=lmoyect
      lmoyect=.FALSE.
      jnxyo = 1
      CALL writevar (kargoutpartvar,vectx(:),jnxyo)
      lmoyect=lmoyectold
!
! --- desallocation of arrays
      IF (allocated(vectx)) deallocate (vectx)
      IF (allocated(part_vect)) deallocate (part_vect)
!
! -6.- Write the zone file header
! -------------------------------
!
      zon_li = MIN(zon_lcut_i,dom_size_i)
      zon_lj = MIN(zon_lcut_j,dom_size_j)
      zon_lk = MIN(zon_lcut_k,dom_size_k)
      zon_lt = MIN(zon_lcut_t,dom_size_t)
      zon_jpi = subdom_size_i + 2 * zon_li
      zon_jpj = subdom_size_j + 2 * zon_lj
      zon_jpk = subdom_size_k + 2 * zon_lk
      zon_jpt = subdom_size_t + 2 * zon_lt
!
      lmoyectold=lmoyect
      lmoyect=.FALSE.
      CALL writehdrzon(kargoutzon, &
     &     zon_jpi,zon_jpj,zon_jpk,zon_jpt,jpbub,jpz)
      lmoyect=lmoyectold
!
! -7.- Generate and write influence bubbles
! -----------------------------------------
!
! --- allocation zon_vect
      allocate ( zon_vect(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt,1:1), &
     &     stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      zon_vect(:,:,:,:,:) = FREAL(0.0)

! --- Position and size of the subdomain in the bubble
      zjic = FREAL(zon_jpi+1) / 2
      zjis = FREAL(subdom_size_i) / 2
      shifti = zjic / FREAL(zon_lw_i)
      shifti = EXP ( - shifti * shifti )

      zjjc = FREAL(zon_jpj+1) / 2
      zjjs = FREAL(subdom_size_j) / 2
      shiftj = zjjc / FREAL(zon_lw_j)
      shiftj = EXP ( - shiftj * shiftj )

      zjkc = FREAL(zon_jpk+1) / 2
      zjks = FREAL(subdom_size_k) / 2
      shiftk = zjkc / FREAL(zon_lw_k)
      shiftk = EXP ( - shiftk * shiftk )

      zjtc = FREAL(zon_jpt+1) / 2
      zjts = FREAL(subdom_size_t) / 2
      shiftt = zjtc / FREAL(zon_lw_t)
      shiftt = EXP ( - shiftt * shiftt )

      DO jt = 1,zon_jpt
      DO jk = 1,zon_jpk
      DO jj = 1,zon_jpj
      DO ji = 1,zon_jpi

        ! distance to subdomain in i, j, k, t
        zjid = MAX( ABS( FREAL(ji) - zjic ) - zjis , 0._kr ) / FREAL(zon_lw_i)
        zjjd = MAX( ABS( FREAL(jj) - zjjc ) - zjjs , 0._kr ) / FREAL(zon_lw_j)
        zjkd = MAX( ABS( FREAL(jk) - zjkc ) - zjks , 0._kr ) / FREAL(zon_lw_k)
        zjtd = MAX( ABS( FREAL(jt) - zjtc ) - zjts , 0._kr ) / FREAL(zon_lw_t)

        ! decrease to zero at the external edges
        expi = EXP ( - zjid * zjid )
        expi = ( expi - shifti ) / ( 1.0_kr - shifti )
        expi = MAX ( expi , 0._kr )

        expj = EXP ( - zjjd * zjjd )
        expj = ( expj - shiftj ) / ( 1.0_kr - shiftj )
        expj = MAX ( expj , 0._kr )

        expk = EXP ( - zjkd * zjkd )
        expk = ( expk - shiftk ) / ( 1.0_kr - shiftk )
        expk = MAX ( expk , 0._kr )

        expt = EXP ( - zjtd * zjtd )
        expt = ( expt - shiftt ) / ( 1.0_kr - shiftt )
        expt = MAX ( expt , 0._kr )

        ! combine expontential decrease in all directions
        zon_vect(ji,jj,jk,jt,1) = expi * expj * expk * expt

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!
      DO jbub = 1, jpbub
        CALL writenbubzon(kargoutzon,zon_vect(:,:,:,:,:),jbub)
      ENDDO
!
! --- desallocation of arrays
      IF (allocated(zon_vect)) deallocate (zon_vect)
!
! -8.- Generate and write zone pointers
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
      ptbubdepth(:,:) = zon_lk + 1
      ptbubtime(:,:) = zon_lt + 1
!
      DO jdta = 1,dtaend
        ptbubidx(1:jpz,jdta) = 1
        ptdtalon(1:jpz,jdta) = pini_i(1:jpz)
        ptdtalat(1:jpz,jdta) = pini_j(1:jpz)
        ptdtadepth(1:jpz,jdta) = pini_k(1:jpz)
        ptdtatime(1:jpz,jdta) = pini_t(1:jpz)
      ENDDO
!
      CALL writeptzon (kargoutzon,ptbubidx, &
     &           ptdtalon,ptdtalat,ptdtadepth,ptdtatime, &
     &           ptbublon,ptbublat,ptbubdepth,ptbubtime)
!
! --- deallocation
      IF (allocated(ptbubtime)) deallocate (ptbubtime)
      IF (allocated(ptbubdepth)) deallocate (ptbubdepth)
      IF (allocated(ptbublat)) deallocate (ptbublat)
      IF (allocated(ptbublon)) deallocate (ptbublon)
      IF (allocated(ptdtatime)) deallocate (ptdtatime)
      IF (allocated(ptdtadepth)) deallocate (ptdtadepth)
      IF (allocated(ptdtalat)) deallocate (ptdtalat)
      IF (allocated(ptdtalon)) deallocate (ptdtalon)
      IF (allocated(ptbubidx)) deallocate (ptbubidx)
      IF (allocated(pini_i)) deallocate (pini_i)
      IF (allocated(pini_j)) deallocate (pini_j)
      IF (allocated(pini_k)) deallocate (pini_k)
      IF (allocated(pini_t)) deallocate (pini_t)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algolocal','calclocal')
 1001 CALL printerror2(0,1001,3,'algolocal','calclocal')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algolocal
