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
! ---                   MKINTERP2DH.F90                           ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-02 (L. Parent)                          ---
! --- modification : 99-05 (C.E. Testut)                        ---
! --- modification : 99-11 (J.M. Brankart)                      ---
! --- modification : 01-06 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE mkinterp2dh
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkinterp2dh
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC interp2dh

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE interp2dh (kvectdbs,kvectrmsdbs,kgridijdbs, &
     &     kvecto,kvectorms,kgridijkobs,kposcoefobs, &
     &     kjpoend,ldbsrms,kspvaldbs,kjobs,kvectreducedta)
!---------------------------------------------------------------------
!
!  Purpose :
!  -------
!  This routine build the observation operator
!      by computing the interpolation coeffcients
!  Method :
!  ------
! --- technique 1 (CET): for an irregular grid of type 1,2
! --- technique 2 (LP) : for an irregular grid of type 1,2
!                            but constant and imposed dx interval
! --- technique 3 (CET): for an irregular grid of type 1,2
!                            (much quicker than technique 1 i
!                            and slightly quicker than 2)
! --- technique 4 (JMB): for an irregular grid of type 1,2,3
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
      use mod_coord , only : longi,latj,gridij
      use utilmkto
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(in) :: kvectdbs
      BIGREAL, dimension(:), intent(in) :: kvectrmsdbs
      TYPE (type_gridij), dimension(:), intent(inout) :: kgridijdbs
      BIGREAL, dimension(:), intent(out) :: kvecto
      BIGREAL, dimension(:), intent(out) :: kvectorms
      TYPE (type_gridijk), dimension(:), intent(out) :: kgridijkobs
      TYPE (type_poscoef), dimension(:,:), intent(out) :: kposcoefobs
      INTEGER, intent(out) :: kjpoend
      LOGICAL, intent(in) :: ldbsrms
      BIGREAL, intent(in) :: kspvaldbs
      INTEGER, intent(in) :: kjobs
      BIGREAL, intent(in), dimension(:), optional :: kvectreducedta
!----------------------------------------------------------------------
! local declarations (for technique 4)
! ====================================
      LOGICAL, dimension(:,:),allocatable :: viewed
      LOGICAL :: localized,stopped
      INTEGER, dimension(4) :: ivic,jvic
      INTEGER :: jiloc,jjloc,icell
      INTEGER :: ngrow,kngrd,nobsloc,nobsrej1,nobsrej2
      BIGREAL :: dcell,dcellmin,dcellold
      BIGREAL, dimension(4) :: jix,jjy
      BIGREAL :: jox,joy,r,s
! local declarations
! ==================
      INTEGER :: allocok,kjpdbssize,kjpitpsize,kjpisize,kjpjsize
      BIGREAL :: a,b,c,d,w1,w2,w3,w4
      INTEGER :: indobs,inddbs
      INTEGER, dimension(:,:),allocatable :: tabind
      INTEGER, dimension(:,:),allocatable :: tabnbobs
      TYPE type_inddbs
          INTEGER :: indi
          INTEGER :: indj
      END TYPE type_inddbs
      TYPE (type_inddbs), dimension(:),allocatable :: tabinddbs
      INTEGER :: inddtamsk,jdta
      INTEGER :: ji,jj,jk,jt,jy
      INTEGER :: jdbs,jo,jbas,jmil,jhaut
      INTEGER :: indice,nbpoint
      LOGICAL :: ptvalid
      INTEGER :: technique,increment
      BIGREAL :: latmin,latmax,lonmin,lonmax
! --- CETBOX
!      PARAMETER (dx=0.2544_kr)
! --- CETATL1
!      PARAMETER (dx=1.0_kr)
! --- LPTDH8
!      PARAMETER (dx=1.0_kr)
      BIGREAL, parameter :: dx=1.0_kr
!----------------------------------------------------------------------
! statement functions (technique 4)
! =================================
      BIGREAL :: fside,eucl_distance,square_distance,cell_distance
      LOGICAL :: sameside,insquare,incell
      BIGREAL :: fx1,fx2,fy1,fy2,fx,fy,fxa,fxb,fya,fyb
      BIGREAL :: fn1,fn2,fn3,fn4,fr,fs,baryc,calcr,calcs
      INTEGER :: jfdbs,jfi,jfj
      fn1(fr,fs) = FREAL(0.25) * ( 1 - r ) * ( 1 - s )
      fn2(fr,fs) = FREAL(0.25) * ( 1 + r ) * ( 1 - s )
      fn3(fr,fs) = FREAL(0.25) * ( 1 + r ) * ( 1 + s )
      fn4(fr,fs) = FREAL(0.25) * ( 1 - r ) * ( 1 + s )
      fside(fx,fx1,fx2,fy,fy1,fy2)=(fy-fy1)*(fx2-fx1)-(fx-fx1)*(fy2-fy1)
      sameside(fxa,fxb,fx1,fx2,fya,fyb,fy1,fy2) =  &
     &     ( ( fside(fxa,fx1,fx2,fya,fy1,fy2) *  &
     &                   fside(fxb,fx1,fx2,fyb,fy1,fy2) ) .GE. FREAL(0.))
      eucl_distance(fx1,fx2,fy1,fy2) =  &
     &          (fx1-fx2)*(fx1-fx2) + (fy1-fy2)*(fy1-fy2)
      baryc(fxa,fxb,fx1,fx2) = FREAL(0.25) * ( fxa + fxb + fx1 + fx2 )
      calcr(jfi,jfj,jfdbs) = &
     &            ( kgridijdbs(jfdbs)%longi - longi(jfi) ) / &
     &            ( longi(jfi+1)            - longi(jfi) )
      calcs(jfi,jfj,jfdbs) = &
     &            ( kgridijdbs(jfdbs)%latj  - latj (jfj) ) / &
     &            ( latj (jfj+1)            - latj (jfj) )
      insquare(jfi,jfj,jfdbs) =  &
     &            ( longi(jfi)   .LE. kgridijdbs(jfdbs)%longi ) .AND. &
     &            ( longi(jfi+1) .GT. kgridijdbs(jfdbs)%longi ) .AND. &
     &            ( latj(jfj)    .LE. kgridijdbs(jfdbs)%latj  ) .AND. &
     &            ( latj(jfj+1)  .GT. kgridijdbs(jfdbs)%latj  )
      incell(jfi,jfj,jfdbs) =  &
     &   sameside ( kgridijdbs(jfdbs)%longi   , gridij(jfi+1,jfj+1)%longi, &
     &              gridij(jfi,jfj)%longi     , gridij(jfi+1,jfj)%longi, &
     &              kgridijdbs(jfdbs)%latj    , gridij(jfi+1,jfj+1)%latj, &
     &              gridij(jfi,jfj)%latj      , gridij(jfi+1,jfj)%latj )    .AND. &
     &   sameside ( kgridijdbs(jfdbs)%longi   , gridij(jfi,jfj+1)%longi, &
     &              gridij(jfi+1,jfj)%longi   , gridij(jfi+1,jfj+1)%longi, &
     &              kgridijdbs(jfdbs)%latj    , gridij(jfi,jfj+1)%latj, &
     &              gridij(jfi+1,jfj)%latj    , gridij(jfi+1,jfj+1)%latj )  .AND. &
     &   sameside ( kgridijdbs(jfdbs)%longi   , gridij(jfi,jfj)%longi, &
     &              gridij(jfi+1,jfj+1)%longi , gridij(jfi,jfj+1)%longi, &
     &              kgridijdbs(jfdbs)%latj    , gridij(jfi,jfj)%latj, &
     &              gridij(jfi+1,jfj+1)%latj  , gridij(jfi,jfj+1)%latj )    .AND. &
     &   sameside ( kgridijdbs(jfdbs)%longi   , gridij(jfi+1,jfj)%longi, &
     &              gridij(jfi,jfj+1)%longi   , gridij(jfi,jfj)%longi, &
     &              kgridijdbs(jfdbs)%latj    , gridij(jfi+1,jfj)%latj, &
     &              gridij(jfi,jfj+1)%latj    , gridij(jfi,jfj)%latj )
      square_distance(jfi,jfj,jfdbs) =  &
     &   eucl_distance ( kgridijdbs(jfdbs)%longi                 , &
     &                   ( longi(jfi) + longi(jfi+1)  ) * FREAL(0.5)    , &
     &                   kgridijdbs(jfdbs)%latj                  , &
     &                   ( latj(jfj)  + latj(jfj+1)   ) * FREAL(0.5)    )
      cell_distance(jfi,jfj,jfdbs) =  &
     &   eucl_distance ( kgridijdbs(jfdbs)%longi                          , &
     &   baryc ( gridij(jfi,jfj)%longi     , gridij(jfi+1,jfj)%longi , &
     &           gridij(jfi+1,jfj+1)%longi , gridij(jfi,jfj+1)%longi   )  , &
     &                   kgridijdbs(jfdbs)%latj                           , &
     &   baryc ( gridij(jfi,jfj)%latj      , gridij(jfi+1,jfj)%latj  , &
     &           gridij(jfi+1,jfj+1)%latj  , gridij(jfi,jfj+1)%latj    )  )
!----------------------------------------------------------------------
! statement functions
! ===================
      BIGREAL :: fcoefa,fcoefb,fcoefc,fcoefd
      BIGREAL :: fw1,fw2,fw3,fw4,jfa,jfb,jfc,jfd
      fcoefa(jfdbs,jfi)=kgridijdbs(jfdbs)%longi-longi(jfi)
      fcoefb(jfi)=longi(jfi+1)-longi(jfi)
      fcoefc(jfdbs,jfj)=kgridijdbs(jfdbs)%latj-latj(jfj)
      fcoefd(jfj)=latj(jfj+1)-latj(jfj)
      fw1(jfa,jfb,jfc,jfd)=(FREAL(1.0)-jfa/jfb)*(FREAL(1.0)-jfc/jfd)
      fw2(jfa,jfb,jfc,jfd)=(jfa/jfb)*(FREAL(1.0)-jfc/jfd)
      fw3(jfa,jfb,jfc,jfd)=(jfc/jfd)*(FREAL(1.0)-jfa/jfb)
      fw4(jfa,jfb,jfc,jfd)=(jfa/jfb)*(jfc/jfd)
!----------------------------------------------------------------------
! dynamic allocation
! ==================
      kjpdbssize=size(kposcoefobs,1)
      kjpitpsize=size(kposcoefobs,2)
! --- allocation tabinddbs
      allocate ( tabinddbs(1:kjpdbssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabinddbs(:) = type_inddbs(0,0)
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&          routine mkinterp2dh             &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&      Interpolation of observations       &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation
! -------------------
!
! -0.1- Initialisation of parameters
! ----------------------------------
!
      kjpoend=0
!
! -0.2- Initialisation of interpolation points
! --------------------------------------------
!
      IF (kjpitpsize.NE.4) GOTO 1000
!
! 4 points interpolation scheme
! the 4 points closest to observation are used
!
! index of the grid points (1,2,3,4):
! a) for techniques 1,2,3
! 3-------------4
! |     .obs    |
! 1-------------2
! b) for techniques 4
! 4-------------3
! |     .obs    |
! 1-------------2
!
      a=FREAL(0.0)
      b=FREAL(0.0)
      c=FREAL(0.0)
      d=FREAL(0.0)
      w1=FREAL(0.0)
      w2=FREAL(0.0)
      w3=FREAL(0.0)
      w4=FREAL(0.0)
!          
! -0.3- Initialisation of mask parameter
! --------------------------------------
!
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
      IF (obs_dim(indobs,inddbs).LT.2) GOTO 1000
      IF (dta_dim(indobs).NE.2) GOTO 1000
      inddtamsk=0
      LOOP1 : DO jdta=1,dtaend
         IF (dta_ord(jdta).EQ.indobs) THEN
            inddtamsk=(varend-1)+jdta
            EXIT LOOP1
         ENDIF
      ENDDO LOOP1
      IF (inddtamsk.EQ.0) GOTO 1000
!
      kjpisize=dta_jpi(indobs)
      kjpjsize=dta_jpj(indobs)
! --- allocation tabind
      allocate ( tabind(1:kjpisize,1:kjpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabind (:,:) = 0
! --- allocation tabnbobs
      allocate ( tabnbobs(1:kjpisize,1:kjpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabnbobs (:,:) = 0
!
! -0.4- Initialisation of indices array 
! -------------------------------------
!
      jk=1
      jt=1
      jy=dta_ind(indobs)
      DO jj=1,kjpjsize
      DO ji=1,kjpisize
         increment = ABS(IBITS(mask(ji,jj,jk,jt),inddtamsk,1))
         tabind(ji,jj)=jy*increment
         jy=jy+increment
      ENDDO
      ENDDO  
      IF ((dta_nbr(indobs).NE.(jy-dta_ind(indobs))) &
     &     .AND.dta_dim(indobs).EQ.2) GOTO 1000
!
! -0.5- Selection of interpolation technique
! ------------------------------------------
!
      SELECT CASE(dtangrd(indobs))
      CASE(1)
         technique=3
      CASE(2)
         technique=3
      CASE(3,4)
         technique=4
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -0.6- Initialisation control print
! ----------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' routine mkinterp2dh'
         WRITE(numout,*) 'inddtamsk  =',inddtamsk
         WRITE(numout,*) 'dta_nam()  =',dta_nam(indobs) &
     &        (1:lenv(dta_nam(indobs)))
         WRITE(numout,*) 'obs_nam()  =',obs_nam(indobs,inddbs) &
     &        (1:lenv(obs_nam(indobs,inddbs)))
         WRITE(numout,*) 'jpitpdbs   =',kjpitpsize	
         WRITE(numout,*) 'ldbsrms    =',ldbsrms	
         WRITE(numout,*) 'technique  =',technique	
         WRITE(numout,*) 'dx         =',dx	
         WRITE(numout,*) 'reducedta  =',present(kvectreducedta)	
         WRITE(numout,*) ' '
      ENDIF
      IF ((technique.GT.4).OR.(technique.LT.1)) GOTO 1000
!
      SELECT CASE(dtangrd(indobs))
      CASE(1,2)
         lonmin=MINVAL(longi(:))
         lonmax=MAXVAL(longi(:))
         latmin=MINVAL(latj(:))
         latmax=MAXVAL(latj(:))
      CASE(3,4)
         lonmin=MINVAL(gridij(:,:)%longi, &
     &        mask=(IBITS(mask(:,:,jk,jt),inddtamsk,1).NE.0))
         lonmax=MAXVAL(gridij(:,:)%longi, &
     &        mask=(IBITS(mask(:,:,jk,jt),inddtamsk,1).NE.0))
         latmin=MINVAL(gridij(:,:)%latj, &
     &        mask=(IBITS(mask(:,:,jk,jt),inddtamsk,1).NE.0))
         latmax=MAXVAL(gridij(:,:)%latj, &
     &        mask=(IBITS(mask(:,:,jk,jt),inddtamsk,1).NE.0))
      CASE DEFAULT
         GOTO 1000
      END SELECT
      PRINT '(4(A,F7.2))','removing observations which are not in ', &
     &             lonmin,' ',latmin,' ',lonmax,' ',latmax
      ji=0
      jj=0
      DO jdbs=1,kjpdbssize
         IF (kgridijdbs(jdbs)%latj.NE.kspvaldbs) THEN                    
            IF ((kgridijdbs(jdbs)%longi.LT.lonmin).OR. &
     &           (kgridijdbs(jdbs)%longi.GT.lonmax).OR. &
     &           (kgridijdbs(jdbs)%latj.LT.latmin).OR. &
     &           (kgridijdbs(jdbs)%latj.GT.latmax)) THEN
               kgridijdbs(jdbs)%latj=kspvaldbs
               ji=ji+1
            ENDIF
            jj=jj+1
         ENDIF
      ENDDO
      print *,'number of observations removed: ',ji,'out of ',jj
!
! -1.- Localize observation and set the observation operator
! ----------------------------------------------------------
!
      SELECT CASE (technique)
      CASE (1)
! --- on fait le tri suivant kjpisize
         jo=0
         DO jj=1,kjpjsize-1
         DO ji=1,kjpisize-1
            ptvalid = .TRUE.
            ptvalid = ( ptvalid .AND. &
     &           (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) )
            ptvalid = ( ptvalid .AND. &
     &           (IBITS(mask(ji+1,jj,jk,jt),inddtamsk,1).NE.0) )
            ptvalid = ( ptvalid .AND. &
     &           (IBITS(mask(ji,jj+1,jk,jt),inddtamsk,1).NE.0) )
            ptvalid = ( ptvalid .AND. &
     &           (IBITS(mask(ji+1,jj+1,jk,jt),inddtamsk,1).NE.0) )
            IF (ptvalid.AND.present(kvectreducedta)) THEN
               ptvalid = ( ptvalid .AND. &
     &              (kvectreducedta(tabind(ji,jj)).NE.FREAL(0.0)) )
               ptvalid = ( ptvalid .AND. &
     &              (kvectreducedta(tabind(ji+1,jj)).NE.FREAL(0.0)) )
               ptvalid = ( ptvalid .AND. &
     &              (kvectreducedta(tabind(ji,jj+1)).NE.FREAL(0.0)) )
               ptvalid = ( ptvalid .AND. &
     &              (kvectreducedta(tabind(ji+1,jj+1)).NE.FREAL(0.0)) )
            ENDIF
!                 
            IF (ptvalid) THEN
               DO jdbs=1,kjpdbssize
                  IF (kgridijdbs(jdbs)%latj.NE.kspvaldbs) THEN
!                     
                     IF ( (longi(ji).LE.kgridijdbs(jdbs)%longi) .AND. &
     &                    (longi(ji+1).GT.kgridijdbs(jdbs)%longi) .AND. &
     &                    (latj(jj).LE.kgridijdbs(jdbs)%latj) .AND. &
     &                    (latj(jj+1).GT.kgridijdbs(jdbs)%latj ) ) THEN
! --- on calcule les poids
!                           PRINT *,'on n''en tient un :'
!     $                          ,jdbs,ji,jj,longi(ji),longi(ji+1),latj(jj),latj(jj+1)
                        a=fcoefa(jdbs,ji)
                        b=fcoefb(ji)
                        c=fcoefc(jdbs,jj)
                        d=fcoefd(jj)
                        w1=fw1(a,b,c,d)
                        w2=fw2(a,b,c,d)
                        w3=fw3(a,b,c,d)
                        w4=fw4(a,b,c,d)
!                           PRINT *,'a,b,c,d,w1,w2,w3,w4'
!     $                          ,a,b,c,d,w1,w2,w3,w4
! --- on incremente                       
                        jo= jo+1
! --- on affecte les valeurs
                        kgridijkobs(jo)%longi=kgridijdbs(jdbs)%longi
                        kgridijkobs(jo)%latj=kgridijdbs(jdbs)%latj
                        kvecto(jo)=kvectdbs(jdbs)
                        IF (ldbsrms) THEN
                           kvectorms(jo)=kvectrmsdbs(jdbs)
                        ENDIF
                        kposcoefobs(jo,1)%pos=tabind(ji,jj)
                        kposcoefobs(jo,1)%coef=w1
                        kposcoefobs(jo,2)%pos=tabind(ji+1,jj)
                        kposcoefobs(jo,2)%coef=w2
                        kposcoefobs(jo,3)%pos=tabind(ji,jj+1)
                        kposcoefobs(jo,3)%coef=w3
                        kposcoefobs(jo,4)%pos=tabind(ji+1,jj+1)
                        kposcoefobs(jo,4)%coef=w4
! --- on desactive
                        kgridijdbs(jdbs)%latj=kspvaldbs
!     
                     ENDIF
                  ENDIF	
               ENDDO
            ENDIF	
         ENDDO
         ENDDO
         kjpoend=jo
         PRINT *,'kjpoend=',kjpoend
!         
      CASE(2,3)
! --- on fait le tri suivant les obs
!
! --- compter le nombre de points par carre
         jo=0
         IF (technique.EQ.2) THEN
           DO jdbs=1,kjpdbssize
               IF (kgridijdbs(jdbs)%latj.NE.kspvaldbs) THEN
                  ji=int( (kgridijdbs(jdbs)%longi-longi(1))/dx + 1 )
                  jj=1
                  IF ((ji.LE.(kjpisize-1)).AND.(ji.GE.1)) THEN
                  DO WHILE ( jj.LE.(kjpjsize-1) )
                     ptvalid = .TRUE.
                     ptvalid = ( ptvalid .AND. &
     &                    (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) )
                     ptvalid = ( ptvalid .AND. &
     &                    (IBITS(mask(ji+1,jj,jk,jt),inddtamsk,1).NE.0) )
                     ptvalid = ( ptvalid .AND. &
     &                    (IBITS(mask(ji,jj+1,jk,jt),inddtamsk,1).NE.0) )
                     ptvalid = ( ptvalid .AND. &
     &                    (IBITS(mask(ji+1,jj+1,jk,jt),inddtamsk,1).NE.0) )
                     IF (ptvalid.AND.present(kvectreducedta)) THEN
                        ptvalid = ( ptvalid .AND. &
     &                       (kvectreducedta(tabind(ji,jj)).NE.FREAL(0.0)) )
                        ptvalid = ( ptvalid .AND. &
     &                       (kvectreducedta(tabind(ji+1,jj)).NE.FREAL(0.0)) )
                        ptvalid = ( ptvalid .AND. &
     &                       (kvectreducedta(tabind(ji,jj+1)).NE.FREAL(0.0)) )
                        ptvalid = ( ptvalid .AND. &
     &                       (kvectreducedta(tabind(ji+1,jj+1)).NE.FREAL(0.0)) )
                     ENDIF
                     ptvalid = ( ptvalid .AND. &
     &                    (latj(jj).LE.kgridijdbs(jdbs)%latj) )
                     ptvalid = ( ptvalid .AND. &
     &                    (latj(jj+1).GT.kgridijdbs(jdbs)%latj ))
                     IF (ptvalid) THEN
!     --- on incremente             
                        jo= jo+1
! --- on incremente tabnbobs
                        tabnbobs(ji,jj)=tabnbobs(ji,jj) + 1
!
                        tabinddbs(jdbs) = type_inddbs ( ji , jj )
                        jj=kjpjsize                     
                     ENDIF               
                     jj=jj+1                              
                  ENDDO
                  ENDIF
               ENDIF                
            ENDDO
         ELSE
            DO jdbs=1,kjpdbssize
               IF (kgridijdbs(jdbs)%latj.NE.kspvaldbs) THEN
                  ptvalid = .TRUE.
! --- selection du carre en i par dichotomie
                  jbas=1
                  jmil=kjpisize/2
                  jhaut=kjpisize
                  ptvalid = ( ptvalid .AND. &
     &                 (longi(jbas).LE.kgridijdbs(jdbs)%longi) )
                  ptvalid = ( ptvalid .AND. &
     &                 (longi(jhaut).GT.kgridijdbs(jdbs)%longi ))
                  IF (ptvalid) THEN
                     DO WHILE ( jbas.NE.jmil )
                        IF (kgridijdbs(jdbs)%longi.GE.longi(jmil)) THEN
                           jbas=jmil
                           jhaut=jhaut
                           jmil=jbas+(jhaut-jbas)/2
                        ELSE
                           jbas=jbas
                           jhaut=jmil
                           jmil=jbas+(jhaut-jbas)/2
                        ENDIF
                     ENDDO
                     ji=jmil
! --- selection du carre en j par dichotomie
                     jbas=1
                     jmil=kjpjsize/2
                     jhaut=kjpjsize
                     ptvalid = ( ptvalid .AND. &
     &                    (latj(jbas).LE.kgridijdbs(jdbs)%latj) )
                     ptvalid = ( ptvalid .AND. &
     &                    (latj(jhaut).GT.kgridijdbs(jdbs)%latj ))
                     IF (ptvalid) THEN
                        DO WHILE ( jbas.NE.jmil )
                           IF (kgridijdbs(jdbs)%latj.GE.latj(jmil)) THEN
                              jbas=jmil
                              jhaut=jhaut
                              jmil=jbas+(jhaut-jbas)/2
                           ELSE
                              jbas=jbas
                              jhaut=jmil
                              jmil=jbas+(jhaut-jbas)/2
                           ENDIF
                        ENDDO
                        jj=jmil
!
                        ptvalid = ( ptvalid .AND. &
     &                       (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) )
                        ptvalid = ( ptvalid .AND. &
     &                       (IBITS(mask(ji+1,jj,jk,jt),inddtamsk,1).NE.0) )
                        ptvalid = ( ptvalid .AND. &
     &                       (IBITS(mask(ji,jj+1,jk,jt),inddtamsk,1).NE.0) )
                        ptvalid = ( ptvalid .AND. &
     &                       (IBITS(mask(ji+1,jj+1,jk,jt),inddtamsk,1).NE.0) )
                        IF (ptvalid.AND.present(kvectreducedta)) THEN
                           ptvalid = ( ptvalid .AND. &
     &                          (kvectreducedta(tabind(ji,jj)).NE.FREAL(0.0)) )
                           ptvalid = ( ptvalid .AND. &
     &                          (kvectreducedta(tabind(ji+1,jj)).NE.FREAL(0.0)) )
                           ptvalid = ( ptvalid .AND. &
     &                          (kvectreducedta(tabind(ji,jj+1)).NE.FREAL(0.0)) )
                           ptvalid = ( ptvalid .AND. &
     &                          (kvectreducedta(tabind(ji+1,jj+1)).NE.FREAL(0.0)) )
                        ENDIF
                        IF (ptvalid) THEN
! --- on incremente             
                           jo= jo+1
! --- on incremente tabnbobs
                           tabnbobs(ji,jj)=tabnbobs(ji,jj) + 1
!
                           tabinddbs(jdbs) = type_inddbs ( ji , jj )
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF                
            ENDDO
         ENDIF                
         kjpoend=jo
         PRINT *,'kjpoend=',kjpoend
! --- reindexation de tabnbobs
         indice=1
         nbpoint=0
         DO jj=1,kjpjsize
         DO ji=1,kjpisize
            nbpoint=tabnbobs(ji,jj)
            tabnbobs(ji,jj)=indice
            indice=indice+nbpoint
         ENDDO        
         ENDDO        
! --- Affectation :
!         print *,kjpdbssize,kjpisize,kjpjsize,kjpoend
         DO jdbs=1,kjpdbssize      
            IF ((tabinddbs(jdbs)%indi.NE.0) &
     &           .AND.(tabinddbs(jdbs)%indj.NE.0)) THEN
! --- calcul de jo
               jo=tabnbobs(tabinddbs(jdbs)%indi,tabinddbs(jdbs)%indj)
! --- on incremente tabnbobs
               tabnbobs(tabinddbs(jdbs)%indi,tabinddbs(jdbs)%indj)=jo+1
! --- on calcule les poids                       
               a=fcoefa(jdbs,tabinddbs(jdbs)%indi)
               b=fcoefb(tabinddbs(jdbs)%indi)
               c=fcoefc(jdbs,tabinddbs(jdbs)%indj)
               d=fcoefd(tabinddbs(jdbs)%indj)
               w1=fw1(a,b,c,d)
               w2=fw2(a,b,c,d)
               w3=fw3(a,b,c,d)
               w4=fw4(a,b,c,d)
! --- on affecte les valeurs
               kgridijkobs(jo)%longi=kgridijdbs(jdbs)%longi
               kgridijkobs(jo)%latj=kgridijdbs(jdbs)%latj
               kvecto(jo)=kvectdbs(jdbs)
               IF (ldbsrms) THEN
                  kvectorms(jo)=kvectrmsdbs(jdbs)
               ENDIF
               kposcoefobs(jo,1)%pos=tabind(tabinddbs(jdbs)%indi, &
     &              tabinddbs(jdbs)%indj)
               kposcoefobs(jo,1)%coef=w1
               kposcoefobs(jo,2)%pos=tabind(tabinddbs(jdbs)%indi+1, &
     &              tabinddbs(jdbs)%indj)
               kposcoefobs(jo,2)%coef=w2
               kposcoefobs(jo,3)%pos=tabind(tabinddbs(jdbs)%indi, &
     &              tabinddbs(jdbs)%indj+1)
               kposcoefobs(jo,3)%coef=w3
               kposcoefobs(jo,4)%pos=tabind(tabinddbs(jdbs)%indi+1, &
     &              tabinddbs(jdbs)%indj+1)
               kposcoefobs(jo,4)%coef=w4
            ENDIF
         ENDDO
      CASE(4)
         kngrd = dtangrd(indobs)
         nobsloc = 0
         nobsrej1 = 0
         nobsrej2 = 0
!
! --- dynamic allocation
         allocate ( viewed(1:kjpisize,1:kjpjsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         viewed (:,:) = .FALSE.
!
! --- set starting research cell for 1st observation
         jiloc = kjpisize/2
         jjloc = kjpjsize/2
!
! --- loop on valid database elements
         jo=0
         DO jdbs=1,kjpdbssize
            IF (kgridijdbs(jdbs)%latj.NE.kspvaldbs) THEN
!
! !! ADDED FOR CREG4
! --- set starting research cell
               jiloc = kjpisize/2
               jjloc = kjpjsize/2
! !! END ADDED FOR CREG4
!
               IF (nprint.GE.2) THEN
                  IF (MOD(jdbs-1,MAX(1,kjpdbssize/5)).EQ.0) THEN
                     print *, ' ------ Observation ',jdbs, '/',kjpdbssize, &
     &                        '   has been considered...'
                  ENDIF
               ENDIF
!
! --- initialize state of research description parameters
               viewed(:,:) = .FALSE.
               stopped  = .FALSE.
               ngrow = 0
               dcellold = -1.
!
               SELECT CASE(kngrd)
                  CASE(1,2)
                     localized = insquare(jiloc,jjloc,jdbs)
                  CASE(3,4)
                     localized = incell(jiloc,jjloc,jdbs)
                  CASE DEFAULT
                     GOTO 1000
               END SELECT
!
! --- follow a path accross the mesh towards the observation cell
               DO WHILE ( (.NOT.localized) .AND. (.NOT.stopped) )
!
                  viewed(jiloc,jjloc) = .TRUE.
                  ivic(1) = MAX(1,jiloc-1)
                  jvic(1) = jjloc
                  ivic(2) = jiloc
                  jvic(2) = MAX(1,jjloc-1)
                  ivic(3) = MIN(jiloc+1,kjpisize-1)
                  jvic(3) = jjloc
                  ivic(4) = jiloc
                  jvic(4) = MIN(jjloc+1,kjpjsize-1)
!
! --- select cell nearest to obseration
                  dcellmin = FREAL(-1.)
                  DO icell = 1,4
                     IF ( .NOT.viewed(ivic(icell),jvic(icell)) ) THEN
                        SELECT CASE(kngrd)
                           CASE(1,2)
                        dcell = square_distance(ivic(icell),jvic(icell),jdbs)
                           CASE(3,4)
                        dcell = cell_distance(ivic(icell),jvic(icell),jdbs)
                           CASE DEFAULT
                              GOTO 1000
                        END SELECT
                        IF ( (dcell.LT.dcellmin) .OR. (dcellmin.EQ.FREAL(-1.)) ) THEN
                           dcellmin = dcell
                           jiloc = ivic(icell)
                           jjloc = jvic(icell)
                        ENDIF
                     ENDIF
                  ENDDO
!
! --- detect abnormal end of path
                  IF (dcellmin.EQ.FREAL(-1.)) THEN
                     stopped = .TRUE.
                     nobsrej1 = nobsrej1 + 1
                     IF (nprint.GE.3) THEN
                        WRITE(numout,*) ' Observation ',jdbs, &
     &                                  ' out of bounds (no more path)'
                     ENDIF
                  ENDIF
!
                  IF (ngrow.GT.3) THEN
                     stopped = .TRUE.
                     nobsrej2 = nobsrej2 + 1
                     IF (nprint.GE.3) THEN
                        WRITE(numout,*) ' Observation ',jdbs, &
     &                               ' out of bounds (path going away)'
                     ENDIF
                  ENDIF
!
! --- store information when path is going away from observation
                  IF (ngrow.GT.0) THEN
                     ngrow = ngrow + 1
                  ENDIF
!
                  IF ( (dcellmin.GT.dcellold) .AND. (dcellold.GT.0.) ) THEN
                     ngrow = 1
                  ENDIF
!
                  dcellold = dcellmin
!
! --- test current cell
                  SELECT CASE(kngrd)
                     CASE(1,2)
                        localized = insquare(jiloc,jjloc,jdbs)
                     CASE(3,4)
                        localized = incell(jiloc,jjloc,jdbs)
                     CASE DEFAULT
                        GOTO 1000
                  END SELECT
               ENDDO
!
! --- save localized observation in .obs file
               IF (localized) THEN
                  ji = jiloc
                  jj = jjloc
!
                  nobsloc = nobsloc + 1
                  IF (nprint.GE.3) THEN
                     WRITE(numout,*) ' Observation ',jdbs, &
     &                               ' localized : ',ji,jj
                  ENDIF
!
! --- save only if grid cell in not masked
                  ptvalid = .TRUE.
                  ptvalid = ( ptvalid .AND. &
     &               (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) )
                  ptvalid = ( ptvalid .AND. &
     &               (IBITS(mask(ji+1,jj,jk,jt),inddtamsk,1).NE.0) )
                  ptvalid = ( ptvalid .AND. &
     &               (IBITS(mask(ji,jj+1,jk,jt),inddtamsk,1).NE.0) )
                  ptvalid = ( ptvalid .AND. &
     &               (IBITS(mask(ji+1,jj+1,jk,jt),inddtamsk,1).NE.0) )
                  IF (ptvalid.AND.present(kvectreducedta)) THEN
                     ptvalid = ( ptvalid .AND. &
     &                    (kvectreducedta(tabind(ji,jj)).NE.FREAL(0.0)) )
                     ptvalid = ( ptvalid .AND. &
     &                    (kvectreducedta(tabind(ji+1,jj)).NE.FREAL(0.0)) )
                     ptvalid = ( ptvalid .AND. &
     &                    (kvectreducedta(tabind(ji,jj+1)).NE.FREAL(0.0)) )
                     ptvalid = ( ptvalid .AND. &
     &                    (kvectreducedta(tabind(ji+1,jj+1)).NE.FREAL(0.0)) )
                  ENDIF
                  IF (ptvalid) THEN
                     jo= jo+1
!
                     IF (nprint.GE.3) THEN
                        IF (MOD(jo,1000).EQ.0) THEN
                            print *, ' Observation ',jo, &
     &                               '   localized: ',ji,jj
                        ENDIF
                     ENDIF
!
! --- compute reduced coordinates
!
                     SELECT CASE(kngrd)
                        CASE(1,2)
!
                           r = FREAL(2.0)*calcr(ji,jj,jdbs)-FREAL(1.0)
                           s = FREAL(2.0)*calcs(ji,jj,jdbs)-FREAL(1.0)
!
                        CASE(3,4)
!
                           jox = kgridijdbs(jdbs)%longi
                           joy = kgridijdbs(jdbs)%latj
!
                           jix (1) = gridij(ji,jj)%longi
                           jix (2) = gridij(ji+1,jj)%longi
                           jix (3) = gridij(ji+1,jj+1)%longi
                           jix (4) = gridij(ji,jj+1)%longi
                           jjy (1) = gridij(ji,jj)%latj
                           jjy (2) = gridij(ji+1,jj)%latj
                           jjy (3) = gridij(ji+1,jj+1)%latj
                           jjy (4) = gridij(ji,jj+1)%latj
!
                           CALL mkxytors(jox,joy,jix,jjy,r,s)
!
                        CASE DEFAULT
                           GOTO 1000
                     END SELECT
!
! --- set the interpolation weights
!
                     w1=fn1(r,s)
                     w2=fn2(r,s)
                     w3=fn3(r,s)
                     w4=fn4(r,s)
!
! --- fill the observation arrays
!
                     kgridijkobs(jo)%longi=kgridijdbs(jdbs)%longi
                     kgridijkobs(jo)%latj=kgridijdbs(jdbs)%latj
                     kvecto(jo)=kvectdbs(jdbs)
                     IF (ldbsrms) THEN
                        kvectorms(jo)=kvectrmsdbs(jdbs)
                     ENDIF
                     kposcoefobs(jo,1)%pos=tabind(ji,jj)
                     kposcoefobs(jo,1)%coef=w1
                     kposcoefobs(jo,2)%pos=tabind(ji+1,jj)
                     kposcoefobs(jo,2)%coef=w2
                     kposcoefobs(jo,3)%pos=tabind(ji,jj+1)
                     kposcoefobs(jo,3)%coef=w4
                     kposcoefobs(jo,4)%pos=tabind(ji+1,jj+1)
                     kposcoefobs(jo,4)%coef=w3
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         kjpoend=jo
         PRINT *,'kjpoend=',kjpoend
!         
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -2.- Control print:
! -------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Maximum size of data set                 = ',kjpdbssize
         WRITE(numout,*) ' Total number of observations             = ', &
     &                                         nobsloc + nobsrej1 + nobsrej2
         WRITE(numout,*) ' Total number of rejected observations    = ',  &
     &                                         nobsrej1 + nobsrej2
         WRITE(numout,*) ' Total number of localized observations   = ',nobsloc
         WRITE(numout,*) ' Total number of selected observations    = ',kjpoend
         WRITE(numout,*) ' '
      ENDIF
!      IF (kjpoend.EQ.0) GOTO 102
      IF (kjpoend.EQ.0) PRINT *, 'Warning : Empty observation vector'
!
! --- desallocate arrays
      IF (allocated(tabind)) deallocate(tabind)
      IF (allocated(tabnbobs)) deallocate(tabnbobs)
      IF (allocated(tabinddbs)) deallocate(tabinddbs)
      IF (allocated(viewed)) deallocate(viewed)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkinterp2dh','mkinterp2dh')
 1001 CALL printerror2(0,1001,3,'mkinterp2dh','mkinterp2dh')
!
 102  WRITE (texterror,*) 'the domain is probably incompatible'
      CALL printerror2(0,102,3,'mkinterp2dh','mkinterp2dh', &
     &      comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkinterp2dh
