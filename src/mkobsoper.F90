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
! ---                    MKOBSOPER.F90                            ---
! ---                                                           ---
! --- original     : 2018-06 JM Brankart                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  obsoper
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkobsoper
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC obsoper

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE obsoper (kvectdbs,kgriddbs, &
     &           kvecto,kgridobs,kposcoefobs, &
     &           kjpoend,kspvaldbs,kjobs,kvectreducedta)
!---------------------------------------------------------------------
!
!  Purpose :
!  -------
!  This routine build the observation operator
!      by localizaing observations in the model grid and
!      by computing the interpolation coefficients
!
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
      use mod_coord , only : longi,latj,levk,time,gridij
      use ensdam_interp
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      BIGREAL, dimension(:), intent(in) :: kvectdbs
      TYPE (type_grid4d), dimension(:), intent(inout) :: kgriddbs
      BIGREAL, dimension(:), intent(out) :: kvecto
      TYPE (type_gridijk), dimension(:), intent(out) :: kgridobs
      TYPE (type_poscoef), dimension(:,:), intent(out) :: kposcoefobs
      INTEGER, intent(out) :: kjpoend
      BIGREAL, intent(in) :: kspvaldbs
      INTEGER, intent(in) :: kjobs
      BIGREAL, intent(in), dimension(:), optional :: kvectreducedta
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,kjpdbssize,kjpitpsize
      INTEGER :: indobs,inddbs,inddtamsk,jdta,jdbs,jo,jy,jpitp,jitp
      INTEGER :: jpi,jpj,jpk,jpt,ji,jj,jk,jt,i,j,k,l,dimi,dimj,dimk,dimt
      INTEGER, dimension(:,:,:,:), allocatable :: tabind
      BIGREAL, dimension(2,2,2,2) :: w
      BIGREAL :: w1d
      LOGICAL :: inmask, located
      CHARACTER(len=20) :: gtype
!----------------------------------------------------------------------
! allocation dynamique
! =====================
      kjpdbssize=size(kposcoefobs,1)
      kjpitpsize=size(kposcoefobs,2)
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&            routine mkobsoper             &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! Initial size of observation vector
! ----------------------------------
      kjpoend=0
!          
! Set dimensions of model grid array
! ----------------------------------
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
!
      inddtamsk=0
      LOOP1 : DO jdta=1,dtaend
         IF (dta_ord(jdta).EQ.indobs) THEN
            inddtamsk=(varend-1)+jdta
            EXIT LOOP1
         ENDIF
      ENDDO LOOP1
      IF (inddtamsk.EQ.0) GOTO 1000
!
      jpi=dta_jpi(indobs)
      jpj=dta_jpj(indobs)
      jpk=dta_jpk(indobs)
      jpt=dta_jpt(indobs)
!
! Identify degenerate dimensions
! ------------------------------
      dimi = 0 ; dimj = 0 ; dimk = 0 ; dimt = 0
      If (jpi.GT.1) dimi = 1
      If (jpj.GT.1) dimj = 1
      If (jpk.GT.1) dimk = 1
      If (jpt.GT.1) dimt = 1
!
      IF (dtangrd(indobs).GT.2) THEN
        IF ((dimi.EQ.0).OR.(dimj.EQ.0)) GOTO 1000
      ENDIF
!
      gtype='cartesian'
      IF (dtangrd(indobs).GT.3) gtype='spherical'
!
! Check dimension of observation operator
! ---------------------------------------
      jpitp=1
      IF (jpi.GT.1) jpitp = jpitp * 2
      IF (jpj.GT.1) jpitp = jpitp * 2
      IF (jpk.GT.1) jpitp = jpitp * 2
      IF (jpt.GT.1) jpitp = jpitp * 2
!
      IF (jpitp.NE.kjpitpsize) GOTO 1000
!
! Control print
! -------------
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' routine mkobsoper'
         WRITE(numout,*) 'inddtamsk  =',inddtamsk
         WRITE(numout,*) 'dta_nam()  =',dta_nam(indobs) &
     &        (1:lenv(dta_nam(indobs)))
         WRITE(numout,*) 'obs_nam()  =',obs_nam(indobs,inddbs) &
     &        (1:lenv(obs_nam(indobs,inddbs)))
         WRITE(numout,*) 'jpitpdbs   =',kjpitpsize
         WRITE(numout,*) 'reducedta  =',present(kvectreducedta)
         WRITE(numout,*) ' '
      ENDIF
!
! -1.- Set indices of grid points (set to 0, if masked)
! -----------------------------------------------------
! --- allocation tabind
      allocate ( tabind(1:jpi,1:jpj,1:jpk,1:jpt), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabind (:,:,:,:) = 0
!
      jy=dta_ind(indobs)
      DO jt=1,jpt
      DO jk=1,jpk
      DO jj=1,jpj
      DO ji=1,jpi
         tabind(ji,jj,jk,jt)=IBITS(mask(ji,jj,jk,jt),inddtamsk,1)
         IF (tabind(ji,jj,jk,jt).NE.0) THEN
           tabind(ji,jj,jk,jt)=jy
           IF (present(kvectreducedta)) THEN
             IF (kvectreducedta(jy).EQ.FREAL(0.0)) THEN
               tabind(ji,jj,jk,jt)=0
             ENDIF
           ENDIF
           jy=jy+1
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO  
!
      IF (dta_nbr(indobs).NE.(jy-dta_ind(indobs))) GOTO 1000
!
! -2.- Preprocessing for irregular 2D grid
! ----------------------------------------
      IF (dtangrd(indobs).GT.2) THEN
        CALL grid2D_init(gridij(:,:)%longi,gridij(:,:)%latj,gtype)
      ENDIF
!
! -3.- Loop on all observations extracted from the database
! ---------------------------------------------------------
      jo=0 ; ji=1 ; jj=1 ; jk=1 ; jt=1
      DO jdbs=1,kjpdbssize
        IF (kgriddbs(jdbs)%lat.NE.kspvaldbs) THEN
          IF (nprint.GE.2) THEN
            IF (MOD(jdbs-1,MAX(1,kjpdbssize/20)).EQ.0) THEN
              print *, ' ------ Observation ',jdbs, '/',kjpdbssize, &
     &                 '   has been considered...'
            ENDIF
          ENDIF
!
! 3a - Locate observation in model grid
! -------------------------------------
          located = .TRUE.
! --- locate observation in model horizontal grid
          IF (dtangrd(indobs).GT.2) THEN
            located=located.AND.grid2D_locate &
     &           (kgriddbs(jdbs)%lon,kgriddbs(jdbs)%lat,ji,jj)
          ELSE
            IF (dimi.EQ.1) THEN
              located=located.AND.grid1D_locate &
     &           (longi,kgriddbs(jdbs)%lon,ji)
            ENDIF
            IF (dimj.EQ.1) THEN
              located=located.AND.grid1D_locate &
     &           (latj,kgriddbs(jdbs)%lat,jj)
            ENDIF
          ENDIF
!
! --- locate observation in model vertical grid
          IF (dimk.EQ.1) THEN
            located=located.AND.grid1D_locate &
     &         (levk,kgriddbs(jdbs)%dep,jk)
          ENDIF
!
! --- locate observation in model time grid
          IF (dimt.EQ.1) THEN
            located=located.AND.grid1D_locate &
     &         (time,kgriddbs(jdbs)%tim,jt)
          ENDIF
!
          IF (located) THEN
!
            IF (nprint.GE.3) THEN
              WRITE(numout,*) ' Observation ',jdbs, &
     &                        ' located : ',ji,jj,jk,jt
            ENDIF
!
! 3b - check if grid cell is inside observation mask
! --------------------------------------------------
            jitp=0 ; inmask = .TRUE.
            DO i=ji,ji+dimi
            DO j=jj,jj+dimj
            DO k=jk,jk+dimk
            DO l=jt,jt+dimt
              jitp=jitp+1
              inmask = inmask .AND. (tabind(i,j,k,l).NE.0)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
!
            IF (inmask) THEN
              jo= jo+1

              IF (nprint.GE.2) THEN
                IF (MOD(jo,1000).EQ.0) THEN
                  print *, ' Observation ',jo, &
     &                     '   in mask: ',ji,jj,jk
                ENDIF
              ENDIF
!
! 3c - compute interpolation weights
! ----------------------------------
! --- initialize interpolation weights
              w(:,:,:,:) = 1.0
!
! --- compute interpolation weights in model horizontal grid
              IF (dtangrd(indobs).GT.2) THEN
                CALL grid2D_interp(kgriddbs(jdbs)%lon, &
     &                     kgriddbs(jdbs)%lat,ji,jj, &
     &                     w(:,:,1,1))
                w(:,:,2,1) = w(:,:,1,1)
                w(:,:,:,2) = w(:,:,:,1)
              ELSE
                IF (dimi.EQ.1) THEN
                  CALL grid1D_interp(longi,kgriddbs(jdbs)%lon,ji,w1d)
                  w(1,:,:,:) = w(1,:,:,:) * (1-w1d)
                  w(2,:,:,:) = w(2,:,:,:) * w1d
                ENDIF
                IF (dimj.EQ.1) THEN
                  CALL grid1D_interp(latj,kgriddbs(jdbs)%lat,jj,w1d)
                  w(:,1,:,:) = w(:,1,:,:) * (1-w1d)
                  w(:,2,:,:) = w(:,2,:,:) * w1d
                ENDIF
              ENDIF
!
! --- compute interpolation weights in model vertical grid
              IF (dimk.EQ.1) THEN
                CALL grid1D_interp(levk,kgriddbs(jdbs)%dep,jk,w1d)
                w(:,:,1,:) = w(:,:,1,:) * (1-w1d)
                w(:,:,2,:) = w(:,:,2,:) * w1d
              ENDIF
!
! --- compute interpolation weights in model time grid
              IF (dimt.EQ.1) THEN
                CALL grid1D_interp(time,kgriddbs(jdbs)%tim,jt,w1d)
                w(:,:,:,1) = w(:,:,:,1) * (1-w1d)
                w(:,:,:,2) = w(:,:,:,2) * w1d
              ENDIF
!
! 3d - fill the routine output arrays
! -----------------------------------
! --- fill the observation arrays
              kgridobs(jo)%longi=kgriddbs(jdbs)%lon
              kgridobs(jo)%latj=kgriddbs(jdbs)%lat
              kgridobs(jo)%levk=kgriddbs(jdbs)%dep
              kvecto(jo)=kvectdbs(jdbs)
!
! --- fill the observation operator arrays
              jitp=0
              DO i=0,dimi
              DO j=0,dimj
              DO k=0,dimk
              DO l=0,dimt
                jitp=jitp+1
                kposcoefobs(jo,jitp)%pos=tabind(ji+i,jj+j,jk+k,jt+l)
                kposcoefobs(jo,jitp)%coef=w(1+i,1+j,1+k,1+l)
              ENDDO
              ENDDO
              ENDDO
              ENDDO
!
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!
! Output final number of observations
! -----------------------------------
      kjpoend=jo
      PRINT *,'Number of observations:',kjpoend
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Number of observations:',kjpoend
         WRITE(numout,*) ' '
      ENDIF
!
      IF (kjpoend.EQ.0) PRINT *, 'Warning : Empty observation vector'
!
! --- deallocation
      IF (allocated(tabind)) deallocate(tabind)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkobsoper','obsoper')
 1001 CALL printerror2(0,1001,3,'mkobsoper','obsoper')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkobsoper
