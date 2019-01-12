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
! ---                    UTILMKTO.F90                             ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-11 (J.M. Brankart)                      ---
! --- revised      : 07-11 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE mkdtatobub : Extract one bubble from a Vy vector
! --- SUBROUTINE mkptbub    : Generate pointers on the Vy object
! ---                         for one bubble
! --- SUBROUTINE mkxytors   : Inverse the quadrangular mapping
! --- SUBROUTINE mkrstoxy   : Quadrangular mapping
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilmkto
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC mkdtatobub,mkptbub,mkxytors,mkrstoxy

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkdtatobub (kvecty,kbub,kjdta, &
     &                 kptdtalon,kptdtalat, &
     &                 kptdtadepth,kptdtatime, &
     &                 kptbublon,kptbublat, &
     &                 kptbubdepth,kptbubtime)
!---------------------------------------------------------------------
!
!  Purpose : Extract one bubble from a Vy vector
!  -------
!  Method :  Find ccordinates limits of the bubble in 4D mask
!  ------    Loop on the mask 4 dimensions,
!                incrementing index (jy) in Vy vector
!            If we are inside the bubble limits,
!                affect Vy values to bubble array
!
!  Input :   kvecty      : Vy vector
!  -----     kjdta       : dta variable index
!            kptdtalon   : dta longitude index
!            kptdtalat   : dta latitude index
!            kptdtadepth : dta depth index
!            kptdtatime  : dta time index
!            kptbublon   : bub longitude index
!            kptbublat   : bub latitude index
!            kptbubdepth : bub depth index
!            kptbubtime  : bub time index
!
!  Output :  kbub      : extracted bubble
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(in) :: kvecty
      BIGREAL, dimension(:,:,:,:), intent(out) :: kbub
      INTEGER, intent(in) :: kjdta, &
     &     kptdtalon,kptdtalat, &
     &     kptdtadepth,kptdtatime, &
     &     kptbublon,kptbublat, &
     &     kptbubdepth,kptbubtime
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: inddta,inddtamsk
      INTEGER :: jpysize,jpisize,jpjsize,jpksize,jptsize
      INTEGER :: jy,ji,jj,jk,jt
      INTEGER :: siz_ji1,siz_ji2, siz_jj1,siz_jj2, &
     &     siz_jk1,siz_jk2, siz_jt1,siz_jt2
      INTEGER :: zlim_ji1,zlim_ji2, zlim_jj1,zlim_jj2, &
     &     zlim_jk1,zlim_jk2, zlim_jt1,zlim_jt2
!----------------------------------------------------------------------
!
! Get size of input arrays
      jpysize=size(kvecty,1)
      jpisize=size(kbub,1)
      jpjsize=size(kbub,2)
      jpksize=size(kbub,3)
      jptsize=size(kbub,4)
!----------------------------------------------------------------------
!
! Initialize bubble to zero and
! get first index of the kjdta variable in Vy vector
      kbub(:,:,:,:)=FREAL(0.0)
      inddta=dta_ord(kjdta)
      inddtamsk=kjdta-1+varend
      jy=dta_ind(inddta)
!
! Find ccordinates limits of the bubble in 4D mask
      zlim_jt1=kptdtatime-(kptbubtime-1)
      zlim_jt2=zlim_jt1-1+jptsize
      zlim_jk1=kptdtadepth-(kptbubdepth-1)
      zlim_jk2=zlim_jk1-1+jpksize
      zlim_jj1=kptdtalat-(kptbublat-1)
      zlim_jj2=zlim_jj1-1+jpjsize
      zlim_ji1=kptdtalon-(kptbublon-1)
      zlim_ji2=zlim_ji1-1+jpisize
!
      siz_jt1=max(zlim_jt1,1)
      siz_jt2=min(zlim_jt2,dta_jpt(inddta))
      siz_jk1=max(zlim_jk1,1)
      siz_jk2=min(zlim_jk2,dta_jpk(inddta))
      siz_jj1=max(zlim_jj1,1)
      siz_jj2=min(zlim_jj2,dta_jpj(inddta))
      siz_ji1=max(zlim_ji1,1)
      siz_ji2=min(zlim_ji2,dta_jpi(inddta))
!
! Loop on the mask 4 dimensions, incrementing index (jy) in Vy vector
      DO jt=1,siz_jt1-1
         DO jk=1,dta_jpk(inddta)
         DO jj=1,dta_jpj(inddta)
         DO ji=1,dta_jpi(inddta)
            IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) jy=jy+1
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!
      DO jt=siz_jt1,siz_jt2
         DO jk=1,siz_jk1-1
            DO jj=1,dta_jpj(inddta)
            DO ji=1,dta_jpi(inddta)
               IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) jy=jy+1
            ENDDO
            ENDDO
         ENDDO
!
         DO jk=siz_jk1,siz_jk2
            DO jj=1,siz_jj1-1
               DO ji=1,dta_jpi(inddta)
                  IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0)  &
     &                 jy=jy+1
               ENDDO
            ENDDO
!
            DO jj=siz_jj1,siz_jj2
               DO ji=1,siz_ji1-1
                  IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0)  &
     &                 jy=jy+1
               ENDDO
               DO ji=siz_ji1,siz_ji2
                  IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
!     
!                    Affect Vy values to bubble array
                     kbub(ji+1-zlim_ji1,jj+1-zlim_jj1, &
     &                    jk+1-zlim_jk1,jt+1-zlim_jt1)=kvecty(jy)
                     jy=jy+1
!
                  ENDIF
               ENDDO
               DO ji=siz_ji2+1,dta_jpi(inddta)
                  IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0)  &
     &                 jy=jy+1
               ENDDO
            ENDDO
!     
            DO jj=siz_jj2+1,dta_jpj(inddta)
               DO ji=1,dta_jpi(inddta)
                  IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0)  &
     &                 jy=jy+1
               ENDDO
            ENDDO
         ENDDO
!
         DO jk=siz_jk2+1,dta_jpk(inddta)
            DO jj=1,dta_jpj(inddta)
            DO ji=1,dta_jpi(inddta)
               IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) jy=jy+1
            ENDDO
            ENDDO
         ENDDO
      ENDDO
!
      DO jt=siz_jt2+1,dta_jpt(inddta)
         DO jk=1,dta_jpk(inddta)
         DO jj=1,dta_jpj(inddta)
         DO ji=1,dta_jpi(inddta)
            IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) jy=jy+1
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      IF (dta_nbr(inddta).NE.(jy-dta_ind(inddta))) GOTO 1000
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilmkto','mkdtatobub')
 1001 CALL printerror2(0,1001,3,'utilmkto','mkdtatobub')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkptbub (kptbub,kjdta, &
     &                 kptdtalon,kptdtalat, &
     &                 kptdtadepth,kptdtatime, &
     &                 kptbublon,kptbublat, &
     &                 kptbubdepth,kptbubtime)
!---------------------------------------------------------------------
!
!  Purpose : Generate pointers on the Vy object for one bubble
!  -------
!  Method :  In presence of connections (largconnect=.TRUE.):
!  ------      Compute 4D bubbles indices (in 1D arrays)
!                  taking the connections into account
!              Use these 4D indices to select the bubble pointers
!                  from de global pointers (pty)
!                  precomputed in the mkconnect routine
!            Else:
!              Find ccordinates limits of the bubble in 4D mask
!              Loop on the mask 4 dimensions,
!                  incrementing index (jy) in Vy vector
!              If we are inside the bubble limits,
!                  affect Vy indices to pointer array
!
!  Input :   kjdta       : dta variable index
!  -----     kptdtalon   : dta longitude index
!            kptdtalat   : dta latitude index
!            kptdtadepth : dta depth index
!            kptdtatime  : dta time index
!            kptbublon   : bub longitude index
!            kptbublat   : bub latitude index
!            kptbubdepth : bub depth index
!            kptbubtime  : bub time index
!
!  Output :  kptbub      : bubble with pointers
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_coord
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, dimension(:,:,:,:), intent(out) :: kptbub
      INTEGER, intent(in) :: kjdta, &
     &     kptdtalon,kptdtalat, &
     &     kptdtadepth,kptdtatime, &
     &     kptbublon,kptbublat, &
     &     kptbubdepth,kptbubtime
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: inddta,inddtamsk
      INTEGER :: jpisize,jpjsize,jpksize,jptsize,jpbsize
      INTEGER :: jy,ji,jj,jk,jt,jb,jcon,jline1,jline2
      INTEGER :: siz_ji1,siz_ji2, siz_jj1,siz_jj2, &
     &     siz_jk1,siz_jk2, siz_jt1,siz_jt2
      INTEGER :: zlim_ji1,zlim_ji2, zlim_jj1,zlim_jj2, &
     &     zlim_jk1,zlim_jk2, zlim_jt1,zlim_jt2
!----------------------------------------------------------------------
!
! Get size of input arrays
      jpisize=size(kptbub,1)
      jpjsize=size(kptbub,2)
      jpksize=size(kptbub,3)
      jptsize=size(kptbub,4)
!---------------------------------------------------------------------
!
! Initialize bubble to zero and
! get index of the kjdta variable
      kptbub(:,:,:,:)=0
      inddta=dta_ord(kjdta)
!
      IF (largconnect) THEN
!
! Transform bubble 4D indices to four 1D array
        jb=0
        DO jt=1,jptsize
        DO jk=1,jpksize
          DO jj=1,jpjsize
          DO ji=1,jpisize
            jb=jb+1
            ptbubcon(jb,1)=kptdtalon-kptbublon+ji
            ptbubcon(jb,2)=kptdtalat-kptbublat+jj
            ptbubcon(jb,3)=kptdtadepth-kptbubdepth+jk
            ptbubcon(jb,4)=kptdtatime-kptbubtime+jt
          ENDDO
          ENDDO
        ENDDO
        ENDDO
        jpbsize=jb
!
! Modify bubble 4D indices in 1D array according to grid connections
        DO jcon=1,jpcon
          DO jline1=1,2
          DO jb=1,jpbsize
            IF (signcon(jcon,jline1)*ptbubcon(jb,dir1con(jcon,jline1)) &
     &                     .GT.  signcon(jcon,jline1)*idxcon(jcon,jline1)) THEN
              IF (ptbubcon(jb,dir2con(jcon,jline1)).GE.mincon(jcon,jline1)) THEN
              IF (ptbubcon(jb,dir2con(jcon,jline1)).LE.maxcon(jcon,jline1)) THEN
                 jline2=3-jline1
                 ptbubcon(jb,dir1con(jcon,jline2)) = idxcon(jcon,jline2) &
     &                  - signcon(jcon,jline2) * signcon(jcon,jline1) &
     &                  * ( ptbubcon(jb,dir1con(jcon,jline1)) - idxcon(jcon,jline1) )
                 IF (dir1con(jcon,jline1).EQ.dir1con(jcon,jline2)) THEN
                    ptbubcon(jb,dir2con(jcon,jline2)) = mincon(jcon,jline2) &
     &                  - signcon(jcon,jline2) * signcon(jcon,jline1) &
     &                  * ( ptbubcon(jb,dir2con(jcon,jline1)) - mincon(jcon,jline1) )
                 ELSE
                    ptbubcon(jb,dir2con(jcon,jline2)) = mincon(jcon,jline2) &
     &                  + signcon(jcon,jline2) * signcon(jcon,jline1) &
     &                  * ( ptbubcon(jb,dir2con(jcon,jline1)) - mincon(jcon,jline1) )
                 ENDIF
              ENDIF
              ENDIF
            ENDIF
          ENDDO
          ENDDO
        ENDDO
!
! Mask modified bubble 4D indices according to grid limits
        ptbubmask(:)=.TRUE.
        WHERE(ptbubcon(:,1).GT.dta_jpi(inddta)) ptbubmask(:)=.FALSE.
        WHERE(ptbubcon(:,2).GT.dta_jpj(inddta)) ptbubmask(:)=.FALSE.
        WHERE(ptbubcon(:,3).GT.dta_jpk(inddta)) ptbubmask(:)=.FALSE.
        WHERE(ptbubcon(:,4).GT.dta_jpt(inddta)) ptbubmask(:)=.FALSE.
        WHERE(ptbubcon(:,1).LT.1) ptbubmask(:)=.FALSE.
        WHERE(ptbubcon(:,2).LT.1) ptbubmask(:)=.FALSE.
        WHERE(ptbubcon(:,3).LT.1) ptbubmask(:)=.FALSE.
        WHERE(ptbubcon(:,4).LT.1) ptbubmask(:)=.FALSE.
!
! Use these 4D indices to select the bubble pointers
! from de global pointers (pty) precomputed in the mkconnect routine
        DO jb=1,jpbsize
           IF (ptbubmask(jb)) THEN
              ptbubcon(jb,0)=   &
     &          pty(ptbubcon(jb,1),ptbubcon(jb,2),ptbubcon(jb,3),ptbubcon(jb,4),kjdta)
           ELSE
              ptbubcon(jb,0)=0
           ENDIF
        ENDDO
!
! Transform the 4D indices 1D array back to 4D pointer array
        jb=0
        DO jt=1,jptsize
        DO jk=1,jpksize
          DO jj=1,jpjsize
          DO ji=1,jpisize
            jb=jb+1
            kptbub(ji,jj,jk,jt)=ptbubcon(jb,0)
           ENDDO
           ENDDO
        ENDDO
        ENDDO
!
      ELSE
!
! Get first index of the kjdta variable in Vy vector
        inddtamsk=kjdta-1+varend
        jy=dta_ind(inddta)
!
! Find ccordinates limits of the bubble in 4D mask
        zlim_jt1=kptdtatime-(kptbubtime-1)
        zlim_jt2=zlim_jt1-1+jptsize
        zlim_jk1=kptdtadepth-(kptbubdepth-1)
        zlim_jk2=zlim_jk1-1+jpksize
        zlim_jj1=kptdtalat-(kptbublat-1)
        zlim_jj2=zlim_jj1-1+jpjsize
        zlim_ji1=kptdtalon-(kptbublon-1)
        zlim_ji2=zlim_ji1-1+jpisize
!
        siz_jt1=max(min(zlim_jt1,dta_jpt(inddta)+1),1)
        siz_jt2=min(max(zlim_jt2,0),dta_jpt(inddta))
        siz_jk1=max(min(zlim_jk1,dta_jpk(inddta)+1),1)
        siz_jk2=min(max(zlim_jk2,0),dta_jpk(inddta))
        siz_jj1=max(min(zlim_jj1,dta_jpj(inddta)+1),1)
        siz_jj2=min(max(zlim_jj2,0),dta_jpj(inddta))
        siz_ji1=max(min(zlim_ji1,dta_jpi(inddta)+1),1)
        siz_ji2=min(max(zlim_ji2,0),dta_jpi(inddta))
!
! Loop on the mask 4 dimensions, incrementing index (jy) in Vy vector
        DO jt=1,siz_jt1-1
           DO jk=1,dta_jpk(inddta)
           DO jj=1,dta_jpj(inddta)
           DO ji=1,dta_jpi(inddta)
              IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) jy=jy+1
           ENDDO
           ENDDO
           ENDDO
        ENDDO
!
        DO jt=siz_jt1,siz_jt2
           DO jk=1,siz_jk1-1
              DO jj=1,dta_jpj(inddta)
              DO ji=1,dta_jpi(inddta)
                 IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) jy=jy+1
              ENDDO
              ENDDO
           ENDDO
!
           DO jk=siz_jk1,siz_jk2
              DO jj=1,siz_jj1-1
                 DO ji=1,dta_jpi(inddta)
                    IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0)  &
     &                 jy=jy+1
                 ENDDO
              ENDDO
!
              DO jj=siz_jj1,siz_jj2
                 DO ji=1,siz_ji1-1
                    IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0)  &
     &                 jy=jy+1
                 ENDDO
                 DO ji=siz_ji1,siz_ji2
                    IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
!     
!                      Affect Vy indices to pointer array
                       kptbub(ji+1-zlim_ji1,jj+1-zlim_jj1, &
     &                    jk+1-zlim_jk1,jt+1-zlim_jt1)=jy
                       jy=jy+1
!
                    ENDIF
                 ENDDO
                 DO ji=siz_ji2+1,dta_jpi(inddta)
                    IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0)  &
     &                 jy=jy+1
                 ENDDO
              ENDDO
!     
              DO jj=siz_jj2+1,dta_jpj(inddta)
                 DO ji=1,dta_jpi(inddta)
                    IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0)  &
     &                 jy=jy+1
                 ENDDO
              ENDDO
           ENDDO
!
           DO jk=siz_jk2+1,dta_jpk(inddta)
              DO jj=1,dta_jpj(inddta)
              DO ji=1,dta_jpi(inddta)
                 IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) jy=jy+1
              ENDDO
              ENDDO
           ENDDO
        ENDDO
!
        DO jt=siz_jt2+1,dta_jpt(inddta)
           DO jk=1,dta_jpk(inddta)
           DO jj=1,dta_jpj(inddta)
           DO ji=1,dta_jpi(inddta)
              IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) jy=jy+1
           ENDDO
           ENDDO
           ENDDO
        ENDDO
!
        IF (dta_nbr(inddta).NE.(jy-dta_ind(inddta))) GOTO 1000
!
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilmkto','mkptbub')
 1001 CALL printerror2(0,1001,3,'utilmkto','mkptbub')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkxytors(jox,joy,jix,jjy,r,s)
!
!     Inverse the nonlinear quadrangular mapping 
!     using a bidimensional Newton-Raphson method
!
      IMPLICIT NONE
!
      BIGREAL, intent(in) :: jox,joy
      BIGREAL, dimension(4), intent(in) :: jix,jjy
      BIGREAL, intent(out) :: r,s
!
      BIGREAL :: rho,epsilon,rr,ss
      BIGREAL, dimension(2,2) :: aa,aai,aaj
      BIGREAL, dimension(2) :: bb,tt,cc
!
      INTEGER :: iter
!
      epsilon = FREAL(0.0001)
!
      aa(1,1) = FREAL(0.25) * ( jix(2) + jix(3) - jix(1) - jix(4) )
      aa(1,2) = FREAL(0.25) * ( jix(3) + jix(4) - jix(1) - jix(2) )
      aa(2,1) = FREAL(0.25) * ( jjy(2) + jjy(3) - jjy(1) - jjy(4) )
      aa(2,2) = FREAL(0.25) * ( jjy(3) + jjy(4) - jjy(1) - jjy(2) )
!
      bb(1) = jox - FREAL(0.25) * ( jix(1) + jix(2) + jix(3) + jix(4) )
      bb(2) = joy - FREAL(0.25) * ( jjy(1) + jjy(2) + jjy(3) + jjy(4) )
!
      tt(1) = FREAL(0.25) * ( jix(1) + jix(3) - jix(2) - jix(4) )
      tt(2) = FREAL(0.25) * ( jjy(1) + jjy(3) - jjy(2) - jjy(4) )
!
      rr = FREAL(0.) 
      ss = FREAL(0.)
      LOOP1 : DO iter = 1,100
         aaj(1,1) = aa(1,1) + ss * tt(1)
         aaj(1,2) = aa(1,2) + rr * tt(1)
         aaj(2,1) = aa(2,1) + ss * tt(2)
         aaj(2,2) = aa(2,2) + rr * tt(2)
!
         rho = aaj(1,1) * aaj(2,2) - aaj(1,2) * aaj(2,1)
!
         aai(1,1) = aaj(2,2) / rho
         aai(1,2) = - aaj(1,2) / rho
         aai(2,1) = - aaj(2,1) / rho
         aai(2,2) = aaj(1,1) / rho
!
         cc(1) = rr * ss * tt(1) - bb(1)
         cc(2) = rr * ss * tt(2) - bb(2)
         cc(1) = cc(1) + aa(1,1) * rr + aa(1,2) * ss
         cc(2) = cc(2) + aa(2,1) * rr + aa(2,2) * ss
!
         r = rr - aai(1,1) * cc(1) - aai(1,2) * cc(2)
         s = ss - aai(2,1) * cc(1) - aai(2,2) * cc(2)
!
         IF ( (ABS(r-rr).LT.epsilon) .AND. (ABS(s-ss).LT.epsilon) ) THEN
            exit LOOP1
         ELSE
            rr = r
            ss = s
         ENDIF
      ENDDO LOOP1
!
      IF (r.GT.FREAL(1.)) r = FREAL(1.)
      IF (s.GT.FREAL(1.)) s = FREAL(1.)
      IF (r.LT.FREAL(-1.)) r = FREAL(-1.)
      IF (s.LT.FREAL(-1.)) s = FREAL(-1.)
!
      RETURN
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkrstoxy(jox,joy,jix,jjy,r,s)
!
!     Quadrangular mapping 
!
      IMPLICIT NONE
!
      BIGREAL, intent(out) :: jox,joy
      BIGREAL, dimension(4), intent(in) :: jix,jjy
      BIGREAL, intent(in) :: r,s
!
      BIGREAL :: fn1,fn2,fn3,fn4
!
      fn1 = FREAL(0.25) * ( 1 - r ) * ( 1 - s )
      fn2 = FREAL(0.25) * ( 1 + r ) * ( 1 - s )
      fn3 = FREAL(0.25) * ( 1 + r ) * ( 1 + s )
      fn4 = FREAL(0.25) * ( 1 - r ) * ( 1 + s )
!
      jox = fn1 * jix(1) + fn2 * jix(2) + fn3 * jix(3) + fn4 * jix(4)
      joy = fn1 * jjy(1) + fn2 * jjy(2) + fn3 * jjy(3) + fn4 * jjy(4)
!
      RETURN
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilmkto
