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
C ---                   ALGOZONE.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-11 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C ---
C --- SUBROUTINE calczone
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algozone
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calczone

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calczone(kargoutpartvar,kargoutzon,kargincfg)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Generate a partition of the domain 
CCC                      and the corresponding influence bubbles
CCC  -------
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
      use mod_spacexyo , only : jpx, jpz
      use hioxyo
      use hiozon
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kargoutpartvar,kargoutzon,
     $     kargincfg
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vectx
C
      INTEGER :: jpxsize,jpnest
      INTEGER :: subdomain_size_i,domain_size_i,part_nbr_i,shift_i
      INTEGER :: subdomain_size_j,domain_size_j,part_nbr_j,shift_j
      INTEGER :: zon_li,zon_lj
      INTEGER :: jpsbd, jvar1, jdta1, var_sbd1, jpzh
      INTEGER, dimension (:), allocatable :: var_sbd, dta_bub1
      INTEGER, dimension (:,:), allocatable :: dta_bub
      BIGREAL, dimension (:), allocatable :: bubcor,bubmin
      INTEGER, dimension (:,:), allocatable :: part_index
      INTEGER, dimension (:), allocatable :: pini_i,pini_j
      INTEGER :: ji,jj,jk,jt,jx,jz,jvar,jbub,jdta,jnxyo,jzidx,jnest
      INTEGER :: jimin,jimax,jjmin,jjmax,jpj
      INTEGER :: indvarmsk,indvar,inddta,jz1,jz2,jsbd,jipart,jjpart
      BIGREAL, dimension (:,:), allocatable :: part_vect
      BIGREAL, dimension (:,:,:,:), allocatable :: zon_vect
      BIGREAL, dimension (:,:,:,:,:), allocatable :: zon_vect1
      BIGREAL :: zon_lcut_i,zon_lcut_j,zon_lw_i,zon_lw_j,zon_ji,zon_jj
      BIGREAL, dimension (:), allocatable :: zon_weight_i,zon_weight_j
      BIGREAL, dimension (:,:), allocatable :: zon_weight_ij
      INTEGER, dimension (:,:), allocatable :: ptbubidx
      INTEGER, dimension (:,:), allocatable :: ptdtalon, ptdtalat
      INTEGER, dimension (:,:), allocatable :: ptdtadepth, ptdtatime
      INTEGER, dimension (:,:), allocatable :: ptbublon, ptbublat
      INTEGER, dimension (:,:), allocatable :: ptbubdepth, ptbubtime
      INTEGER :: zji1,zji2,zji3,zji4,zji5,zji6
      INTEGER :: zjj1,zjj2,zjj3,zjj4,zjj5,zjj6
      CHARACTER(len=word80) :: title
      CHARACTER(len=varlg) :: dta_nam1, var_nam1
      LOGICAL :: inmask,lmoyectold,lectinfo,ingrid
      INTEGER ::  allocok,numfila,flagxyo
      INTEGER, dimension (:,:,:,:), allocatable :: reducemask 
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
      BIGREAL :: shift, factor, exparg
      LOGICAL :: flagspct
      INTEGER :: zon_jpj_all,jj1,jj_ini,jj_med
CC----------------------------------------------------------------------
C
      jpnest=nestend
      jpxsize=jpx
C --- allocation vectx
      allocate ( vectx(1:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectx(:) = FREAL(0.0)
C --- allocation var_sbd
      allocate ( var_sbd(1:varend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      var_sbd(:) = 1
C ----------------------------------------------------------------------
      IF (nprint.GE.1) THEN
         WRITE(numout,*)
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine algozone &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*)
      ENDIF
C
C -0.- Initialisation
C -------------------
C
      textexclusion='#'
      numfila=10
      CALL openfile(numfila,kargincfg)
C
      text=readnextline(numfila,textexclusion)
      READ(text,*) subdomain_size_i, subdomain_size_j
      text=readnextline(numfila,textexclusion)
      READ(text,*) zon_lcut_i, zon_lcut_j
      text=readnextline(numfila,textexclusion)
      READ(text,*) zon_lw_i, zon_lw_j
C
      text=readnextline(numfila,textexclusion)
      flagspct = (text(1:4).EQ.'spct')
      IF ((text(1:3).EQ.'end').OR.flagspct) THEN
C
C Default subdomain configuration (all variables together)
C
        jpsbd=1
        var_sbd(1:varend)=1
        jpbub=1
C
C --- allocation bubcor
        allocate ( bubcor(1:jpbub), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        bubcor(:) = FREAL(1.0)
C --- allocation bubmin
        allocate ( bubmin(1:jpbub), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        bubmin(:) = FREAL(0.0)
C --- allocation dta_bub
        allocate ( dta_bub(1:dtaend,1:jpsbd), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        dta_bub(:,:) = 1
C
      ELSE
C
C Read subdomain configuration from config file
C
        READ(text,*) jpsbd
        DO jvar1 = 1,varend
          text=readnextline(numfila,textexclusion)
          READ(text,*) var_nam1, var_sbd1
          IF (var_sbd1.GT.jpsbd) GOTO 1000
          IF (var_sbd1.LE.0) GOTO 1000
          DO jvar = 1,varend
            indvar = var_ord(jvar)
            IF (var_nam1(1:lenv(var_nam1))
     $        .EQ.var_nam(indvar)(1:lenv(var_nam(indvar))) ) THEN
              var_sbd(jvar) = var_sbd1
            ENDIF
          ENDDO
        ENDDO
        text=readnextline(numfila,textexclusion)
        READ(text,*) jpbub
C
C --- allocation bubcor
        allocate ( bubcor(1:jpbub), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        bubcor(:) = FREAL(0.0)
C --- allocation bubmin
        allocate ( bubmin(1:jpbub), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        bubmin(:) = FREAL(0.0)
C --- allocation dta_bub
        allocate ( dta_bub(1:dtaend,1:jpsbd), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        dta_bub(:,:) = 1
C --- allocation dta_bub1
        allocate ( dta_bub1(1:jpsbd), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        dta_bub1(:) = 1
C
        DO jbub = 1,jpbub
          text=readnextline(numfila,textexclusion)
          READ(text,*) bubcor(jbub),bubmin(jbub)
          IF (bubmin(jbub).LT.FREAL(0.0)) GOTO 1000
        ENDDO
C
        DO jdta1 = 1,dtaend
          text=readnextline(numfila,textexclusion)
          READ(text,*) dta_nam1, dta_bub1(1:jpsbd)
          IF (ANY(dta_bub1(1:jpsbd).GT.jpbub)) GOTO 1000
          IF (ANY(dta_bub1(1:jpsbd).LE.0))     GOTO 1000
          DO jdta = 1,dtaend
            inddta = dta_ord(jdta)
            IF (dta_nam1(1:lenv(dta_nam1))
     $        .EQ.dta_nam(inddta)(1:lenv(dta_nam(inddta))) ) THEN
              dta_bub(jdta,1:jpsbd) = dta_bub1(1:jpsbd)
            ENDIF
          ENDDO
        ENDDO
C
        IF (allocated(dta_bub1)) deallocate (dta_bub1)
      ENDIF
C
C Close config file
      CLOSE(10)
C
C -1.- Give subdomain index (jz) to every subdomain (ji,jj)
C ---------------------------------------------------------
C Compute domain size for coarse grid
      domain_size_i = 0
      domain_size_j = 0
      DO jvar = 1,varend
        indvar = var_ord(jvar)
        IF (varnest(indvar).EQ.0) THEN 
          IF (var_jpi(indvar).GT.domain_size_i) 
     $        domain_size_i = var_jpi(indvar)
          IF (var_jpj(indvar).GT.domain_size_j) 
     $        domain_size_j = var_jpj(indvar)
        ENDIF
      ENDDO     
C
C Convert it to finest grid resolution
      DO jnest=1,jpnest
        domain_size_i = domain_size_i * nestxres(jnest)
        domain_size_j = domain_size_j * nestyres(jnest)
      ENDDO
C
C Compute number of subdomains to cover the full grid
      part_nbr_i = ( domain_size_i - 1 ) / subdomain_size_i + 1
      part_nbr_j = ( domain_size_j - 1 ) / subdomain_size_j + 1

C --- allocation part_index
      allocate ( part_index(1:part_nbr_i,1:part_nbr_j), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      part_index(:,:) = 0
C
C Compute subdomain index (jz) for each subdomain (ji,jj)
      jz = 0
      DO jj = part_nbr_j,1,-1
         DO ji = 1,part_nbr_i
           jz = jz + 1
           part_index(ji,jj) = jz
         ENDDO
      ENDDO
      jpzh = jz
C
C -2.- Store south-west grid point (pini_i,pini_j) of each subdomain (jz)
C -----------------------------------------------------------------------
C --- allocation pini_i
      allocate ( pini_i(1:jpzh), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pini_i(:) = 0
C --- allocation pini_j
      allocate ( pini_j(1:jpzh), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pini_j(:) = 0
C
      DO ji = 1,part_nbr_i
      DO jj = 1,part_nbr_j
        pini_i(part_index(ji,jj)) = 1 + (ji-1) * subdomain_size_i
        pini_j(part_index(ji,jj)) = 1 + (jj-1) * subdomain_size_j
      ENDDO
      ENDDO
C
C --- reference to south-west corner of finest grid
      shift_i=0 ; shift_j=0
      DO jnest=1,jpnest
        shift_i = (shift_i + nestxori(jnest) - 1 ) * nestxres(jnest)
        shift_j = (shift_j + nestyori(jnest) - 1 ) * nestyres(jnest)
      ENDDO
      pini_i(:) = pini_i(:) - shift_i
      pini_j(:) = pini_j(:) - shift_j
C
C --- desallocation of arrays
      IF (allocated(part_index)) deallocate (part_index)
C
C -3.1- Create mask from user file (argreducevar)
C -----------------------------------------------
C
      IF (largreducevar) THEN
         lmoyectold=lmoyect
         lmoyect=.FALSE.
         flagxyo=1
         jnxyo=1
         lectinfo=.TRUE.
         CALL readxyo (argreducevar,vectx(:),jnxyo,lectinfo,flagxyo)
         lmoyect=lmoyectold
C --- allocation reducemask
         allocate ( reducemask(1:size(mask,1),1:size(mask,2),
     $        1:size(mask,3),1:size(mask,4)), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         reducemask(:,:,:,:) = 0
C
         DO jvar=1,varend
            indvar=var_ord(jvar)
            indvarmsk=jvar-1
            jx=var_ind(indvar)
            DO jt=1,var_jpt(indvar)
            DO jk=1,var_jpk(indvar)
            DO jj=1,var_jpj(indvar)
            DO ji=1,var_jpi(indvar)
               IF (IBITS(mask(ji,jj,jk,jt),indvarmsk,1).NE.0) THEN
                  IF (vectx(jx).NE.FREAL(0.0)) THEN
                     reducemask(ji,jj,jk,jt) = reducemask(ji,jj,jk,jt) 
     $                    + (2**(indvarmsk))
                  ENDIF
                  jx=jx+1
               ENDIF
            ENDDO
            ENDDO
            ENDDO
            ENDDO
            IF (jx.NE.(var_ind(indvar)+var_nbr(indvar))) GOTO 1000
         ENDDO
C
         vectx(:)=FREAL(0.0)
      ENDIF
C
C -3.2- Eliminate subdomains outside model mask
C ---------------------------------------------
C
      jz = 1
      DO WHILE (jz .LE. jpzh)
C
        inmask = .FALSE.
        DO jvar = 1,varend
          indvarmsk = jvar - 1
          indvar = var_ord(jvar)
C
C Set the boundaries of the subdomain
          jimin=pini_i(jz) 
          jjmin=pini_j(jz) 
          jimax=pini_i(jz)+subdomain_size_i-1
          jjmax=pini_j(jz)+subdomain_size_j-1
C Convert to the right grid level
          DO jnest=jpnest,varnest(indvar)+1,-1
C Round to next coarse grid point
            IF (jimin.LE.1) THEN
              jimin=(jimin-1)/nestxres(jnest)+nestxori(jnest)
            ELSE
              jimin=(jimin-2)/nestxres(jnest)+nestxori(jnest)+1
            ENDIF
            IF (jjmin.LE.1) THEN
              jjmin=(jjmin-1)/nestyres(jnest)+nestyori(jnest)
            ELSE
              jjmin=(jjmin-2)/nestyres(jnest)+nestyori(jnest)+1
            ENDIF
C Round to previous coarse grid point
            IF (jimax.LT.1) THEN
              jimax=(jimax)/nestxres(jnest)+nestxori(jnest)-1   !ok
            ELSE
              jimax=(jimax-1)/nestxres(jnest)+nestxori(jnest)   !ok
            ENDIF
            IF (jjmax.LT.1) THEN
              jjmax=(jjmax)/nestyres(jnest)+nestyori(jnest)-1
            ELSE
              jjmax=(jjmax-1)/nestyres(jnest)+nestyori(jnest)
            ENDIF
          ENDDO
C Check if at least one grid point is inlcuded in the subdomain
C (subdomian maybe empty when converted to coarser grid)
          ingrid=((jimin.LE.jimax).AND.(jjmin.LE.jjmax))
C
C Check the boundaries of the subdomain
          IF (jimin.LT.1) THEN
            IF (jimax.GE.1) THEN
              jimin=1
            ELSE
              ingrid=.FALSE.
            ENDIF
          ENDIF
C
          IF (jimax.GT.var_jpi(indvar)) THEN
            IF (jimin.LE.var_jpi(indvar)) THEN
              jimax=var_jpi(indvar)
            ELSE
              ingrid=.FALSE.
            ENDIF
          ENDIF
C
          IF (jjmin.LT.1) THEN
            IF (jjmax.GE.1) THEN
              jjmin=1
            ELSE
              ingrid=.FALSE.
            ENDIF
          ENDIF
C
          IF (jjmax.GT.var_jpj(indvar)) THEN
            IF (jjmin.LE.var_jpj(indvar)) THEN
              jjmax=var_jpj(indvar)
            ELSE
              ingrid=.FALSE.
            ENDIF
          ENDIF
C
C Check if there is at least one value in mask
          IF (ingrid) THEN
            DO jt=  1,var_jpt(indvar)
            DO jk = 1,var_jpk(indvar)
            DO ji = jimin,jimax
            DO jj = jjmin,jjmax
C
               IF (largreducevar) THEN
                 inmask = inmask .OR. 
     $           (IBITS(reducemask(ji,jj,jk,jt),indvarmsk,1).NE.0) 
               ELSE
                 inmask = inmask .OR. 
     $           (IBITS(mask(ji,jj,jk,jt),indvarmsk,1).NE.0) 
               ENDIF
C
            ENDDO
            ENDDO
            ENDDO
            ENDDO
          ENDIF
C
        ENDDO
C
        IF (.NOT.inmask) THEN
          DO jzidx = jz+1,jpzh
             pini_i(jzidx-1) = pini_i(jzidx)
             pini_j(jzidx-1) = pini_j(jzidx)
          ENDDO
          jpzh = jpzh - 1
        ELSE
          jz = jz + 1
        ENDIF
C
      ENDDO
C
C -4.- Generate partition vector
C ------------------------------
C
C --- desallocation of arrays
      IF (allocated(reducemask)) deallocate (reducemask)
C --- allocation part_vect
      jimin=1-shift_i ; jimax=domain_size_i-shift_i
      jjmin=1-shift_j ; jjmax=domain_size_j-shift_j
      allocate(part_vect(jimin:jimax,jjmin:jjmax),stat=allocok)
      IF (allocok.NE.0) GOTO 1001
      part_vect(:,:) = FREAL(0.0)
C
      DO jz = 1,jpzh
        DO ji=pini_i(jz),MIN(pini_i(jz)+subdomain_size_i-1,jimax)
        DO jj=pini_j(jz),MIN(pini_j(jz)+subdomain_size_j-1,jjmax)
          part_vect(ji,jj)=FREAL(jz)
        ENDDO
        ENDDO
      ENDDO
C
      jx = 0
      DO jvar = 1,varend
         indvar = var_ord(jvar)
         indvarmsk = jvar - 1
         DO jt=1,var_jpt(indvar)
         DO jk=1,var_jpk(indvar)
            DO jj=1,var_jpj(indvar)
            DO ji=1,var_jpi(indvar)
               jipart=ji ; jjpart=jj
               DO jnest=varnest(indvar),jpnest-1
                 jipart=(jipart-nestxori(jnest+1))*nestxres(jnest+1)+1
                 jjpart=(jjpart-nestyori(jnest+1))*nestyres(jnest+1)+1
               ENDDO
               IF (IBITS(mask(ji,jj,jk,jt),indvarmsk,1).NE.0) THEN
                  jx = jx + 1
                  vectx(jx) = part_vect(jipart,jjpart) 
     $                       + ( var_sbd(jvar) - 1) * jpzh
               ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      IF (jx.NE.jpx) GOTO 1000
C
C -5.- Write partition vector
C ---------------------------
C
      lmoyectold=lmoyect
      lmoyect=.FALSE.
      jnxyo = 1
      CALL writevar (kargoutpartvar,vectx(:),jnxyo)
      lmoyect=lmoyectold
C
C --- desallocation of arrays
      IF (allocated(vectx)) deallocate (vectx)
      IF (allocated(var_sbd)) deallocate (var_sbd)
      IF (allocated(part_vect)) deallocate (part_vect)
C
C -6.- Write the zone file header
C -------------------------------
C
      zon_li = MIN(NINT(zon_lcut_i),domain_size_i)
      zon_lj = MIN(NINT(zon_lcut_j),domain_size_j)
      zon_jpi = subdomain_size_i + 2 * zon_li
      zon_jpj = subdomain_size_j + 2 * zon_lj
      zon_jpk = MAXVAL(dta_jpk(dta_ord(1:dtaend)))
      zon_jpt = MAXVAL(dta_jpt(dta_ord(1:dtaend)))
      jpz = jpzh * jpsbd
C
      IF (flagspct) THEN
        jpj = MAXVAL(var_jpj(:))
        zon_jpj_all = jpj
      ELSE
        zon_jpj_all = zon_jpj
      ENDIF
C
      lmoyectold=lmoyect
      lmoyect=.FALSE.
      CALL writehdrzon(kargoutzon,
     $     zon_jpi,zon_jpj_all,zon_jpk,zon_jpt,jpbub,jpz)
      lmoyect=lmoyectold
C
C -7.- Generate and write influence bubbles
C -----------------------------------------
C
C --- allocation zon_weight_i
      allocate ( zon_weight_i(1:zon_li), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      zon_weight_i(:) = FREAL(0.0)
C
C --- allocation zon_weight_j
      allocate ( zon_weight_j(1:zon_lj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      zon_weight_j(:) = FREAL(0.0)
C
C --- allocation zon_weight_ij
      allocate ( zon_weight_ij(1:zon_li,1:zon_lj), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      zon_weight_ij(:,:) = FREAL(0.0)
C
      DO ji = 1,zon_li
         IF (zon_lw_i.EQ.FREAL(0.)) THEN
            zon_weight_i(ji) = FREAL(1.0)
         ELSE
            zon_weight_i(ji) = 
     $           EXP( - FREAL(ji*ji)/FREAL(zon_lw_i*zon_lw_i) )
         ENDIF
      ENDDO
C
      DO jj = 1,zon_lj
         IF (zon_lw_j.EQ.FREAL(0.)) THEN
            zon_weight_j(jj) = FREAL(1.0)
         ELSE
            zon_weight_j(jj) = 
     $           EXP( - FREAL(jj*jj)/FREAL(zon_lw_j*zon_lw_j) )
         ENDIF
      ENDDO
C
      DO ji = 1,zon_li
         DO jj = 1,zon_lj
            IF (zon_lw_i.EQ.FREAL(0.)) THEN
               zon_ji = FREAL(0.)
            ELSE
               zon_ji = FREAL(ji) / zon_lw_i
            ENDIF
            IF (zon_lw_j.EQ.FREAL(0.)) THEN
               zon_jj = FREAL(0.)
            ELSE
               zon_jj = FREAL(jj) / zon_lw_j
            ENDIF
            exparg = FREAL(zon_ji*zon_ji+zon_jj*zon_jj)
            zon_weight_ij(ji,jj) = EXP( - exparg )
C
            IF ((zon_lw_i.EQ.FREAL(0.)).AND.(zon_lw_j.EQ.FREAL(0.))) THEN
              zon_ji = FREAL(ji) / FREAL(zon_li)
              zon_jj = FREAL(jj) / FREAL(zon_lj)
              exparg = FREAL(zon_ji*zon_ji+zon_jj*zon_jj)
              IF (exparg.GT.1.0_kr) THEN
                zon_weight_ij(ji,jj) = 0.0_kr
              ENDIF
            ENDIF
C
         ENDDO
      ENDDO
C
C Decrease to zero at the external edges
C
      IF ((zon_lw_i.EQ.FREAL(0.)).OR.(zon_lw_j.EQ.FREAL(0.))) THEN
        shift = 0.0_kr
        factor = 1.0_kr
      ELSE
        shift = MIN ( FREAL(zon_li)/zon_lw_i , FREAL(zon_lj)/zon_lw_j )
        shift = EXP ( - shift * shift )
        factor = 1.0_kr / ( 1.0_kr - shift )
      ENDIF
C
      zon_weight_i(:) = ( zon_weight_i(:) - shift ) * factor
      zon_weight_j(:) = ( zon_weight_j(:) - shift ) * factor
      zon_weight_ij(:,:) = ( zon_weight_ij(:,:) - shift ) * factor
      WHERE (zon_weight_i(:).LT.0.0_kr) zon_weight_i(:) = 0.0_kr
      WHERE (zon_weight_j(:).LT.0.0_kr) zon_weight_j(:) = 0.0_kr
      WHERE (zon_weight_ij(:,:).LT.0.0_kr) zon_weight_ij(:,:) = 0.0_kr
C
C --- allocation zon_vect
      allocate ( zon_vect(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt),
     $     stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      zon_vect(:,:,:,:) = FREAL(0.0)
C
C --- allocation zon_vect1
      allocate ( zon_vect1(1:zon_jpi,1:zon_jpj_all,1:zon_jpk,1:zon_jpt,1:1),
     $     stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      zon_vect1(:,:,:,:,:) = FREAL(0.0)
C
C --- fill the zon_vect array
      zji1 = 1
      zji2 = zon_li
      zji3 = zon_li + 1
      zji4 = zon_jpi - zon_li
      zji5 = zon_jpi - zon_li + 1
      zji6 = zon_jpi
C
      zjj1 = 1
      zjj2 = zon_lj
      zjj3 = zon_lj + 1
      zjj4 = zon_jpj - zon_lj
      zjj5 = zon_jpj - zon_lj + 1
      zjj6 = zon_jpj
C
      zon_vect(zji3:zji4,zjj3:zjj4,1:zon_jpk,1:zon_jpt) = FREAL(1.0)
C
      DO jt = 1,zon_jpt
      DO jk = 1,zon_jpk

         DO jj = zjj3,zjj4
            zon_vect(zji1:zji2,jj,jk,jt) = zon_weight_i (zon_li:1:-1)
            zon_vect(zji5:zji6,jj,jk,jt) = zon_weight_i (1:zon_li)
         ENDDO
C
         DO ji = zji3,zji4
            zon_vect(ji,zjj1:zjj2,jk,jt) = zon_weight_j (zon_lj:1:-1)
            zon_vect(ji,zjj5:zjj6,jk,jt) = zon_weight_j (1:zon_lj)
         ENDDO
C
         zon_vect(zji1:zji2,zjj1:zjj2,jk,jt)
     $                    = zon_weight_ij(zon_li:1:-1,zon_lj:1:-1)
         zon_vect(zji5:zji6,zjj5:zjj6,jk,jt)
     $                    = zon_weight_ij(1:zon_li,1:zon_lj)
         zon_vect(zji5:zji6,zjj1:zjj2,jk,jt)
     $                    = zon_weight_ij(1:zon_li,zon_lj:1:-1)
         zon_vect(zji1:zji2,zjj5:zjj6,jk,jt)
     $                    = zon_weight_ij(zon_li:1:-1,1:zon_lj)
C
      ENDDO
      ENDDO
C
      IF (flagspct) THEN
        DO jbub = 1, part_nbr_j
          jj_ini = 1 + (jbub-1) * subdomain_size_j
          jj_med = 1 + (jbub-1) * subdomain_size_j + subdomain_size_j/2
          zon_vect1(:,:,:,:,1)=0.0
          DO jj=1,zon_jpj
            jj1 = jj_ini + jj - zon_lj - 1
            IF ((jj1.GE.1).AND.(jj1.LE.zon_jpj_all)) THEN
              zon_vect1(:,jj1,:,:,1)=zon_vect(:,jj,:,:)*bubcor(1)
            ENDIF
          ENDDO
          IF (jj_med.GE.jpj/2+1) THEN
            DO jj=1,jpj/2
              zon_vect1(:,jj,:,:,1) = zon_vect1(:,jpj-jj+1,:,:,1)
            ENDDO
          ELSEIF (jj_med.LT.jpj/2+1) THEN
            DO jj=1,jpj/2
              zon_vect1(:,jpj-jj+1,:,:,1) = zon_vect1(:,jj,:,:,1)
            ENDDO
          ENDIF
          CALL writenbubzon(kargoutzon,zon_vect1(:,:,:,:,:),jbub)
        ENDDO
      ELSE
        DO jbub = 1, jpbub
          zon_vect1(:,:,:,:,1)=zon_vect(:,:,:,:)*bubcor(jbub)
          WHERE (zon_vect1(:,:,:,:,1).LT.bubmin(jbub))
     $       zon_vect1(:,:,:,:,1) = FREAL(0.0)
          CALL writenbubzon(kargoutzon,zon_vect1(:,:,:,:,:),jbub)
        ENDDO
      ENDIF
C
C --- desallocation of arrays
      IF (allocated(zon_vect)) deallocate (zon_vect)
      IF (allocated(zon_vect1)) deallocate (zon_vect1)
      IF (allocated(bubcor)) deallocate (bubcor)
      IF (allocated(bubmin)) deallocate (bubmin)
      IF (allocated(zon_weight_i)) deallocate (zon_weight_i)
      IF (allocated(zon_weight_j)) deallocate (zon_weight_j)
      IF (allocated(zon_weight_ij)) deallocate (zon_weight_ij)
C
C -8.- Generate and write zone pointers
C -------------------------------------
C-
C --- allocation zone pointers
      allocate ( ptbubidx(1:jpz,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ptbubidx(:,:) = 0
C
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
C
      ptdtadepth(:,:) = 1
      ptdtatime(:,:) = 1
      ptbublon(:,:) = zon_li + 1
      ptbublat(:,:) = zon_lj + 1
      ptbubdepth(:,:) = 1
      ptbubtime(:,:) = 1
C
      DO jdta = 1,dtaend
         DO jsbd = 1,jpsbd
            jz1 = 1 + jpzh * (jsbd-1)
            jz2 = jpzh * jsbd
            ptbubidx(jz1:jz2,jdta) = dta_bub(jdta,jsbd)
            ptdtalon(jz1:jz2,jdta) = pini_i(1:jpzh)
            ptdtalat(jz1:jz2,jdta) = pini_j(1:jpzh)
         ENDDO
      ENDDO
C
      IF (flagspct) THEN
        ptbublat(:,:) = 1
        ptdtalat(:,:) = 1
        DO jdta = 1,dtaend
        DO jz = 1,jpzh
          jbub = ( pini_j(jz) - 1 ) / subdomain_size_j + 1
          ptbubidx(jz,jdta) = jbub
        ENDDO
        ENDDO
      ENDIF
C
      CALL writeptzon (kargoutzon,ptbubidx,
     $           ptdtalon,ptdtalat,ptdtadepth,ptdtatime,
     $           ptbublon,ptbublat,ptbubdepth,ptbubtime)
C
C --- deallocation
      IF (allocated(ptbubtime)) deallocate (ptbubtime)
      IF (allocated(ptbubdepth)) deallocate (ptbubdepth)
      IF (allocated(ptbublat)) deallocate (ptbublat)
      IF (allocated(ptbublon)) deallocate (ptbublon)
      IF (allocated(ptdtatime)) deallocate (ptdtatime)
      IF (allocated(ptdtadepth)) deallocate (ptdtadepth)
      IF (allocated(ptdtalat)) deallocate (ptdtalat)
      IF (allocated(ptdtalon)) deallocate (ptdtalon)
      IF (allocated(ptbubidx)) deallocate (ptbubidx)
      IF (allocated(dta_bub)) deallocate (dta_bub)
      IF (allocated(pini_i)) deallocate (pini_i)
      IF (allocated(pini_j)) deallocate (pini_j)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algozone','calczone')
 1001 CALL printerror2(0,1001,3,'algozone','calczone')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algozone
