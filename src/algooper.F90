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
! ---                   ALGOOPER.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 00-02 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 01-11 (C.E. Testut)                        ---
! --- modification : 02-02 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE algoopervct
! --- SUBROUTINE algooperzon
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algooper
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC algoopervct,algooperzon

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoopervct(kflaganlxyo,koutxyo, &
     &     kflagnbxyo,kincfg,kinxyo,kinrefxyo, &
     &     ktextoper,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose : Tools for arithmetic operation on state and basis :
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
      use mod_coord
      use mod_spacexyo , only :  &
     &     poscoefobs,gridijkobs, &
     &     jpoend,jpitpend,jpx,jpxend, &
     &     jpyend,spvalvar,spvaldta,spvalobs
      use hioxyo
      use hiocfg
      use hiogrd
      use ensdam_storng
      use ensdam_stoanam
      use ensdam_anaqua
      use ensdam_interp
      use ensdam_obserror
      use utilvct
      use utilvalid
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(in) :: kflaganlxyo,kflagnbxyo
      CHARACTER(len=*), intent(in) :: koutxyo, &
     &     kincfg,kinxyo,kinrefxyo, &
     &     ktextoper,kconfigo
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vectsin
      BIGREAL, dimension(:), allocatable, save :: vectsinref
      BIGREAL, dimension(:), allocatable, save :: vectsout
!
      BIGREAL, dimension(:), allocatable, save :: vectonew
      BIGREAL, dimension(:), allocatable, save :: vectorms
      BIGREAL, dimension(:), allocatable, save :: vectormsnew
!
      TYPE (type_gridijk), dimension(:), allocatable :: gridijkobsnew
      TYPE (type_poscoef), dimension(:,:), allocatable :: poscoefobsnew
!
      INTEGER :: allocok,jpssize,jpisize,jpjsize,jpksize
      INTEGER :: jpitpsize,sxyend,indsxy,nbr,js0,js1
      INTEGER, dimension(1:nbvar) :: sxy_ord,sxy_jpi,sxy_jpj,sxy_jpk,sxy_jpt
      INTEGER :: jnxyo,js,ji,jk,jt,jpi,jextvar,jjproc,jqua,jsxy
      INTEGER :: flagcfg,flagxyo,ios,jpgroup,jgroup
      INTEGER :: jpoendold,jpitpsizeold,jpoendnew,jpitpsizenew,indobs, &
     &     inddbs,indobs1,inddbs1,jodeb,jofin,jodebnew,jofinnew,jobs,jitp
      LOGICAL :: lectinfo,lmodprint
      CHARACTER(len=bgword), dimension(:), allocatable :: nam_groupoper
      BIGREAL :: spvals,spvalsin,spvalsinref,cst,norm,x,delta,delta2
      BIGREAL :: gpmin,gpmax,gpval,ri
!      TYPE (type_gridijk), dimension(:), allocatable :: gridijkobs
!      TYPE (type_poscoef), dimension(:,:), allocatable :: poscoefobs
      BIGREAL, dimension(:), allocatable, save :: garray
      INTEGER, dimension(:), allocatable, save :: jpoendtab,jpitpendtab
      INTEGER, dimension(:,:), allocatable, save :: obs_nbrold, &
     &     obs_indold,obs_nbrnew,obs_nbrtot,obs_indnew, &
     &     obs_itpold,obs_itpnew
      INTEGER, dimension(:,:,:), allocatable, save :: obs_nbrjgrp, &
     &     obs_indjgrp,obs_itpjgrp
      REAL(KIND=8) :: gran, xcoef, mu, nu, gammak, gammath, betaa, betab
      REAL(KIND=8) :: xmin, xmax, qua, distmin, dist
      BIGREAL, DIMENSION(:), allocatable :: lon, lat
      REAL(KIND=8), PARAMETER :: twopi=2*3.1415926535897932384626
      REAL(KIND=8), PARAMETER :: deg2rad=twopi/360.
      REAL(KIND=8), PARAMETER :: earthrad=6.371e3
      INTEGER :: jday,jdaytarget,year,month,day
      BIGREAL, dimension(:), allocatable, save :: dum1d
      BIGREAL, dimension(:,:), allocatable, save :: dum2d
      BIGREAL, dimension(:,:), allocatable, save :: zvalues
!----------------------------------------------------------------------
      SELECT CASE (kflaganlxyo)
      CASE(1)
         jpssize=jpx
         jpitpsize=1
         spvals=spvalvar
         sxyend=varend
         DO jsxy = 1,sxyend
           sxy_ord(jsxy)=var_ord(jsxy)
           indsxy= sxy_ord(jsxy)
           sxy_jpi(indsxy)=var_jpi(indsxy)
           sxy_jpj(indsxy)=var_jpj(indsxy)
           sxy_jpk(indsxy)=var_jpk(indsxy)
           sxy_jpt(indsxy)=var_jpt(indsxy)
         ENDDO
      CASE(2)
         jpssize=jpyend
         jpitpsize=1
         spvals=spvaldta
         sxyend=dtaend
         DO jsxy = 1,sxyend
           sxy_ord(jsxy)=dta_ord(jsxy)
           indsxy= sxy_ord(jsxy)
           sxy_jpi(indsxy)=dta_jpi(indsxy)
           sxy_jpj(indsxy)=dta_jpj(indsxy)
           sxy_jpk(indsxy)=dta_jpk(indsxy)
           sxy_jpt(indsxy)=dta_jpt(indsxy)
         ENDDO
      CASE(3)
         jpssize=jpoend
         jpitpsize=jpitpend
         spvals=spvalobs
         sxyend=0
      CASE DEFAULT
         GOTO 1000
      END SELECT
! --- allocation vectsout
      allocate ( vectsout(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectsout(:) = FREAL(0.0)
      IF (kflagnbxyo.GE.0) THEN
! --- allocation vectsin
         allocate ( vectsin(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectsin(:) = FREAL(0.0)
      ENDIF
      IF (kflagnbxyo.EQ.2) THEN
! --- allocation vectsinref
         allocate ( vectsinref(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectsinref(:) = FREAL(0.0)
      ENDIF
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine algooper &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
! -0.1- reading or define config.obs :
! ------------------------------------
!
! --- allocation poscoefobs
      allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      IF (kflaganlxyo.LE.2) THEN
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
      IF (kflaganlxyo.EQ.3) THEN
!
! --- allocation gridijkobs
         allocate ( gridijkobs(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
         flagcfg=2
         CALL readcfgobs (kconfigo,flagcfg, &
     &        kgridijkobs=gridijkobs(:))
! --- allocation vectorms
         allocate ( vectorms(1:jpssize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectorms(:)=FREAL(0.0)
         flagcfg=1
         CALL readcfgobs (kconfigo,flagcfg, &
     &        kvectorms=vectorms(:))
!
      ENDIF
!
! -1.- Operation computation :
! -----------------------------
!
      IF (nprint.GE.1) print *, &
     &     '--> Computation of operation ...'
      lectinfo = .FALSE.
      IF (jpproc.GT.limjpnxyo(kflaganlxyo)) GOTO 1003
      DO jnxyo=1+jproc,limjpnxyo(kflaganlxyo),jpproc
         lmodprint=(MOD(jnxyo-1,(limjpnxyo(kflaganlxyo)/5+1)).EQ.0)
         IF ((lmodprint).AND.(nprint.GE.1)) &
     &        print *,'Memory part number : ', &
     &           jnxyo,'/',limjpnxyo(kflaganlxyo)
!
! -1.1- Read in vector :
! ----------------------
!
         IF (kflagnbxyo.GE.1) THEN
            IF (nprint.GE.2) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOOPER : ', &
     &              'Reading in vector'
            ENDIF
!     
            spvalsin=spvalvar
            IF (validextdta(kinxyo)) spvalsin=spvaldta
            IF (validextobs(kinxyo)) spvalsin=spvalobs
            CALL readxyo(kinxyo,vectsin(:), &
     &           jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))   
         ENDIF
!
! -1.2- Read inref vector :
! -------------------------
!
         IF (kflagnbxyo.EQ.2) THEN
            IF (nprint.GE.2) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOOPER : ', &
     &              'Reading inref vector'
            ENDIF
!
            spvalsinref=spvalvar
            IF (validextdta(kinrefxyo)) spvalsinref=spvaldta
            IF (validextobs(kinrefxyo)) spvalsinref=spvalobs
            CALL readxyo(kinrefxyo,vectsinref(:), &
     &           jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))   
         ENDIF
!
! -2- Operation :
! -----------------
!
         SELECT CASE (kflagnbxyo)
         CASE (0)
            SELECT CASE (posit(ktextoper,'_'))
            CASE (2:bgword)        
               READ(ktextoper((posit(ktextoper,'_')+1): &
     &                 lenv(ktextoper)),*,IOSTAT=ios) cst
               IF (ios.NE.0) GOTO 101
               GOTO 101
            CASE (0)        
               SELECT CASE (ktextoper(1:lenv(ktextoper)))
               CASE ('mean')
                  IF (jnxyo.EQ.1+jproc) THEN
                     CALL evalhdrcfgoper(kincfg,jpgroup)
! --- allocation nam_grouparea
                     allocate ( nam_groupoper(1:jpgroup), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     nam_groupoper(:) = 'file#name.dimg'
                     CALL readcfgoper(kincfg,nam_groupoper)
                     DO jgroup=1,jpgroup
                        SELECT CASE (kflaganlxyo)
                        CASE (1)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (2)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (3)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          .OR.validextobs(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE DEFAULT
                           GOTO 1000
                        END SELECT
                     ENDDO
                  ENDIF
                  vectsout(:)=0.0_kr
                  DO jgroup=1,jpgroup
                     spvalsin=spvalvar
                     IF (validextdta(nam_groupoper(jgroup))) &
     &                    spvalsin=spvaldta
                     IF (validextobs(nam_groupoper(jgroup))) &
     &                    spvalsin=spvalobs
                     CALL readxyo(nam_groupoper(jgroup),vectsin(:), &
     &                    jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))  
                     vectsout(:)=vectsout(:)+vectsin(:)
                  ENDDO   
                  vectsout(:)=vectsout(:)/FREAL(jpgroup)
               CASE ('std')
                  allocate ( vectsinref(1:jpssize), stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
!
                  IF (jnxyo.EQ.1+jproc) THEN
                     CALL evalhdrcfgoper(kincfg,jpgroup)
! --- allocation nam_grouparea
                     allocate ( nam_groupoper(1:jpgroup), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     nam_groupoper(:) = 'file#name.dimg'
                     CALL readcfgoper(kincfg,nam_groupoper)
                     DO jgroup=1,jpgroup
                        SELECT CASE (kflaganlxyo)
                        CASE (1)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (2)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (3)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          .OR.validextobs(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE DEFAULT
                           GOTO 1000
                        END SELECT
                     ENDDO
                  ENDIF
                  vectsinref(:)=0.0_kr
                  vectsout(:)=0.0_kr
                  DO jgroup=1,jpgroup
                     spvalsin=spvalvar
                     IF (validextdta(nam_groupoper(jgroup))) &
     &                    spvalsin=spvaldta
                     IF (validextobs(nam_groupoper(jgroup))) &
     &                    spvalsin=spvalobs
                     CALL readxyo(nam_groupoper(jgroup),vectsin(:), &
     &                    jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
                     DO js=1,jpssize
                       delta=vectsin(js)-vectsinref(js)
                       vectsinref(js)=vectsinref(js)+delta/jgroup
                       delta2=vectsin(js)-vectsinref(js)
                       vectsout(js)=vectsout(js)+delta*delta2
                     ENDDO
                  ENDDO
                  vectsout(:)=SQRT(vectsout(:)/FREAL(jpgroup-1))
               CASE ('lincomb')
                  IF (jnxyo.EQ.1+jproc) THEN
                     CALL evalhdrcfgoper(kincfg,jpgroup)
! --- allocation nam_grouparea
                     allocate ( nam_groupoper(1:jpgroup), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     nam_groupoper(:) = 'file#name.dimg'
                     CALL readcfgoper(kincfg,nam_groupoper)
                     DO jgroup=1,jpgroup
                        SELECT CASE (kflaganlxyo)
                        CASE (1)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (2)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (3)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          .OR.validextobs(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE DEFAULT
                           GOTO 1000
                        END SELECT
                     ENDDO
                  ENDIF
                  OPEN(UNIT=11,file='coef.txt')
                  DO jgroup=1,jpgroup
                     READ(11,*) xcoef
                     spvalsin=spvalvar
                     IF (validextdta(nam_groupoper(jgroup))) &
     &                    spvalsin=spvaldta
                     IF (validextobs(nam_groupoper(jgroup))) &
     &                    spvalsin=spvalobs
                     CALL readxyo(nam_groupoper(jgroup),vectsin(:), &
     &                    jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
                     print *, 'coef,minmax',xcoef,MINVAL(vectsin(:)), MAXVAL(vectsin(:))
                     vectsout(:)=vectsout(:)+vectsin(:)*xcoef
                  ENDDO
                  CLOSE(11)
                  print *, 'vectsout',MINVAL(vectsout(:)), MAXVAL(vectsout(:))
               CASE ('zfile')
                  IF (kflaganlxyo.EQ.3) GOTO 101

                  ! Get max number of levels
                  jpisize = MAXVAL(sxy_jpi(1:sxyend))
                  jpjsize = MAXVAL(sxy_jpj(1:sxyend))
                  jpksize = MAXVAL(sxy_jpk(1:sxyend))

                  ! Allocate zvalues
                  allocate ( zvalues(1:jpksize,1:sxyend), stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  ! Allocate dummy arrays
                  allocate ( dum1d(1:jpisize*jpjsize), stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  allocate ( dum2d(1:jpisize,1:jpjsize), stat=allocok )
                  IF (allocok.NE.0) GOTO 1001

                  ! Read z values for each variable
                  CALL openfile(11,kincfg)
                  DO jsxy=1,sxyend
                    indsxy= sxy_ord(jsxy)
                    DO jk=1,sxy_jpk(jsxy)
                      READ(11,*) zvalues(jk,jsxy)
                    ENDDO
                  ENDDO
                  CLOSE(11)

                  ! Loop on variables, time slices and levels
                  js0 = 0
                  DO jsxy=1,sxyend  ! loop on variables
                    indsxy= sxy_ord(jsxy)
                    DO jt=1,sxy_jpt(indsxy)  ! loop on time slices
                    DO jk=1,sxy_jpk(indsxy)  ! loop on levels

                      ! Get size of 2D slice -> nbr
                      CALL mk8vct(dum1d,dum2d,jk,jt,jsxy,nbr,kflaganlxyo)

                      ! Put constant zvalue in the corresponding segment of vector
                      js1 = js0 + nbr - 1
                      vectsout(js0:js1) = zvalues(jk,jsxy)
                      js0 = js0 + nbr
!
                    ENDDO
                    ENDDO
                  ENDDO

               CASE ('probmap')
                  IF (jnxyo.EQ.1+jproc) THEN
                     CALL evalhdrcfgoper(kincfg,jpgroup)
! --- allocation nam_grouparea
                     allocate ( nam_groupoper(1:jpgroup), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     nam_groupoper(:) = 'file#name.dimg'
                     CALL readcfgoper(kincfg,nam_groupoper)
                     DO jgroup=1,jpgroup
                        SELECT CASE (kflaganlxyo)
                        CASE (1)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (2)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (3)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          .OR.validextobs(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE DEFAULT
                           GOTO 1000
                        END SELECT
                     ENDDO
                  ENDIF
                  OPEN(UNIT=11,file='minmax.txt')
                  READ(11,*) xmin, xmax
                  CLOSE(11)
                  DO jgroup=1,jpgroup
                     spvalsin=spvalvar
                     IF (validextdta(nam_groupoper(jgroup))) &
     &                    spvalsin=spvaldta
                     IF (validextobs(nam_groupoper(jgroup))) &
     &                    spvalsin=spvalobs
                     CALL readxyo(nam_groupoper(jgroup),vectsin(:), &
     &                    jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
                     WHERE ((vectsin(:).GE.xmin).AND.(vectsin(:).LE.xmax)) &
     &                    vectsout(:)=vectsout(:)+1.0
                  ENDDO
                  vectsout(:)=vectsout(:)/jpgroup
               CASE ('condave')
                  IF (jnxyo.EQ.1+jproc) THEN
                     CALL evalhdrcfgoper(kincfg,jpgroup)
! --- allocation nam_grouparea
                     allocate ( nam_groupoper(1:jpgroup), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     nam_groupoper(:) = 'file#name.dimg'
                     CALL readcfgoper(kincfg,nam_groupoper)
                     DO jgroup=1,jpgroup
                        SELECT CASE (kflaganlxyo)
                        CASE (1)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (2)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (3)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          .OR.validextobs(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE DEFAULT
                           GOTO 1000
                        END SELECT
                     ENDDO
                  ENDIF
                  OPEN(UNIT=11,file='minmax.txt')
                  READ(11,*) xmin, xmax
                  CLOSE(11)
                  DO jgroup=1,jpgroup
                     spvalsin=spvalvar
                     IF (validextdta(nam_groupoper(jgroup))) &
     &                    spvalsin=spvaldta
                     IF (validextobs(nam_groupoper(jgroup))) &
     &                    spvalsin=spvalobs
                     CALL readxyo(nam_groupoper(jgroup),vectsin(:), &
     &                    jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
                     WHERE ((vectsin(:).GE.xmin).AND.(vectsin(:).LE.xmax)) &
     &                    vectsout(:)=vectsout(:)+vectsin(:)
                  ENDDO
                  vectsout(:)=vectsout(:)/jpgroup
               CASE ('maxabs')
                  IF (jnxyo.EQ.1+jproc) THEN
                     CALL evalhdrcfgoper(kincfg,jpgroup)
! --- allocation nam_grouparea
                     allocate ( nam_groupoper(1:jpgroup), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     nam_groupoper(:) = 'file#name.dimg'
                     CALL readcfgoper(kincfg,nam_groupoper)
                     DO jgroup=1,jpgroup
                        SELECT CASE (kflaganlxyo)
                        CASE (1)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (2)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (3)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          .OR.validextobs(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE DEFAULT
                           GOTO 1000
                        END SELECT
                     ENDDO
                  ENDIF
                  DO jgroup=1,jpgroup
                     spvalsin=spvalvar
                     IF (validextdta(nam_groupoper(jgroup))) &
     &                    spvalsin=spvaldta
                     IF (validextobs(nam_groupoper(jgroup))) &
     &                    spvalsin=spvalobs
                     CALL readxyo(nam_groupoper(jgroup),vectsin(:), &
     &                    jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
                     vectsout(:)=MAX(vectsout(:),ABS(vectsin(:)))
                  ENDDO
               CASE ('meanabs')
                  IF (jnxyo.EQ.1+jproc) THEN
                     CALL evalhdrcfgoper(kincfg,jpgroup)
! --- allocation nam_grouparea
                     allocate ( nam_groupoper(1:jpgroup), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     nam_groupoper(:) = 'file#name.dimg'
                     CALL readcfgoper(kincfg,nam_groupoper)
                     DO jgroup=1,jpgroup
                        SELECT CASE (kflaganlxyo)
                        CASE (1)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (2)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (3)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          .OR.validextobs(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE DEFAULT
                           GOTO 1000
                        END SELECT
                     ENDDO
                  ENDIF
                  DO jgroup=1,jpgroup
                     spvalsin=spvalvar
                     IF (validextdta(nam_groupoper(jgroup))) &
     &                    spvalsin=spvaldta
                     IF (validextobs(nam_groupoper(jgroup))) &
     &                    spvalsin=spvalobs
                     CALL readxyo(nam_groupoper(jgroup),vectsin(:), &
     &                    jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
                     vectsout(:)=vectsout(:)+ABS(vectsin(:))
                  ENDDO
                  vectsout(:)=vectsout(:)/FREAL(jpgroup)
               CASE ('meansqr')
                  IF (jnxyo.EQ.1+jproc) THEN
                     CALL evalhdrcfgoper(kincfg,jpgroup)
! --- allocation nam_grouparea
                     allocate ( nam_groupoper(1:jpgroup), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     nam_groupoper(:) = 'file#name.dimg'
                     CALL readcfgoper(kincfg,nam_groupoper)
                     DO jgroup=1,jpgroup
                        SELECT CASE (kflaganlxyo)
                        CASE (1)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (2)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE (3)
                           IF (.NOT.(validextvar(nam_groupoper(jgroup)) &
     &                          .OR.validextdta(nam_groupoper(jgroup)) &
     &                          .OR.validextobs(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE DEFAULT
                           GOTO 1000
                        END SELECT
                     ENDDO
                  ENDIF
                  DO jgroup=1,jpgroup
                     spvalsin=spvalvar
                     IF (validextdta(nam_groupoper(jgroup))) &
     &                    spvalsin=spvaldta
                     IF (validextobs(nam_groupoper(jgroup))) &
     &                    spvalsin=spvalobs
                     CALL readxyo(nam_groupoper(jgroup),vectsin(:), &
     &                    jnxyo,lectinfo,kflaganlxyo,poscoefobs(:,:))
                     vectsout(:)=vectsout(:)+vectsin(:)*vectsin(:)
                  ENDDO
                  vectsout(:)=vectsout(:)/FREAL(jpgroup)
               CASE ('grpobs')
                  IF (jnxyo.EQ.1+jproc) THEN
                     CALL evalhdrcfgoper(kincfg,jpgroup)
! --- allocation nam_grouparea
                     allocate ( nam_groupoper(1:jpgroup), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     nam_groupoper(:) = 'file#name.obs'
                     CALL readcfgoper(kincfg,nam_groupoper)
                     DO jgroup=1,jpgroup
                        SELECT CASE (kflaganlxyo)
                        CASE (1,2)
                           GOTO 1000
                        CASE (3)
                           IF (.NOT.(validextobs(nam_groupoper(jgroup)) &
     &                          )) GOTO 105
                        CASE DEFAULT
                           GOTO 1000
                        END SELECT
                     ENDDO
!
                  jpoendold=jpoend
                  jpitpsizeold=jpitpend
!
                  allocate ( obs_nbrold(1:nbvar,1:jpndbs), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  obs_nbrold(:,:)=obs_nbr(:,:)
                  allocate ( obs_indold(1:nbvar,1:jpndbs), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  obs_indold(:,:)=obs_ind(:,:)
                  allocate ( obs_itpold(1:nbvar,1:jpndbs), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  obs_itpold(:,:)=obs_itp(:,:)
                  allocate ( obs_nbrnew(1:nbvar,1:jpndbs), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  obs_nbrnew(:,:) = 0
                  allocate ( obs_nbrtot(1:nbvar,1:jpndbs), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  obs_nbrtot(:,:) = 0
                  allocate ( obs_indnew(1:nbvar,1:jpndbs), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  obs_indnew(:,:) = 0
                  allocate ( obs_itpnew(1:nbvar,1:jpndbs), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  obs_itpnew(:,:) = 0
!
                  allocate ( obs_nbrjgrp(1:nbvar,1:jpndbs,1:jpgroup), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  obs_nbrjgrp(:,:,:)=0
                  allocate ( obs_indjgrp(1:nbvar,1:jpndbs,1:jpgroup), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  obs_indjgrp(:,:,:)=0
                  allocate ( obs_itpjgrp(1:nbvar,1:jpndbs,1:jpgroup), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  obs_itpjgrp(:,:,:)=0
!
                  allocate ( jpoendtab(1:jpgroup), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  jpoendtab(:)=0
                  allocate ( jpitpendtab(1:jpgroup), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  jpitpendtab(:)=0
!
                  DO jgroup=1,jpgroup
                     CALL evalhdrobs(nam_groupoper(jgroup), &
     &                    jpoendtab(jgroup),jpitpendtab(jgroup))
                     obs_nbrjgrp(:,:,jgroup)=obs_nbr(:,:)
                     obs_indjgrp(:,:,jgroup)=obs_ind(:,:)
                     obs_itpjgrp(:,:,jgroup)=obs_itp(:,:)
                     obs_nbrnew(:,:) = obs_nbrnew(:,:) &
     &                    + obs_nbrjgrp(:,:,jgroup)
                     obs_itpnew(:,:) = MAX(obs_itpnew(:,:),obs_itp(:,:))
                  ENDDO
                  jpoendnew=SUM(jpoendtab(1:jpgroup))
                  jpitpsizenew = MAX(jpitpend, &
     &                 MAXVAL(jpitpendtab(1:jpgroup)))
!
                  indobs = obs_ord(1)
                  inddbs = obsnord(1)
                  obs_indnew(indobs,inddbs) = 1
                  DO jobs=2,obsend
                     indobs = obs_ord(jobs)
                     inddbs=obsnord(jobs)
                     indobs1 = obs_ord(jobs-1)
                     inddbs1=obsnord(jobs-1)
                     obs_indnew(indobs,inddbs) = obs_indnew(indobs1, &
     &                    inddbs1) + obs_nbrnew(indobs1,inddbs1)
                  ENDDO
!
! --- allocation poscoefobsnew
                  allocate ( poscoefobsnew(1:jpoendnew,1:jpitpsizenew), &
     &                 stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  poscoefobsnew(:,:) = type_poscoef(0,FREAL(0.0))
! --- allocation gridijkobsnew
                  allocate ( gridijkobsnew(1:jpoendnew), stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  gridijkobsnew(:)=type_gridijk(FREAL(0.0),FREAL(0.0), &
     &                 FREAL(0.0))
! --- allocation vectormsnew
                  allocate ( vectormsnew(1:jpoendnew), stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  vectormsnew(:)=FREAL(0.0)
! --- allocation vectonew
                  allocate ( vectonew(1:jpoendnew), stat=allocok )
                  IF (allocok.NE.0) GOTO 1001
                  vectonew(:)=FREAL(0.0)
!
                  IF (allocated(vectsin)) deallocate (vectsin)
                  IF (allocated(vectorms)) deallocate(vectorms)
                  IF (allocated(poscoefobs)) deallocate(poscoefobs)
                  IF (allocated(gridijkobs)) deallocate(gridijkobs)
!
                  DO jgroup=1,jpgroup
!     
                     jpoend=jpoendtab(jgroup)
!                     jpitpsize=jpitpendtab(jgroup)
                     obs_nbr(:,:)=obs_nbrjgrp(:,:,jgroup)
                     obs_ind(:,:)=obs_indjgrp(:,:,jgroup)  
                     obs_itp(:,:)=obs_itpjgrp(:,:,jgroup)                     
! --- allocation poscoefobs
                     allocate ( poscoefobs(1:jpoend,1:jpitpend), &
     &                    stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
! --- allocation gridijkobs
                     allocate ( gridijkobs(1:jpoend), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0), &
     &                    FREAL(0.0))
! --- allocation vectorms
                     allocate ( vectorms(1:jpoend), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     vectorms(:)=FREAL(0.0)
! --- allocation vectsin
                     allocate ( vectsin(1:jpoend), stat=allocok )
                     IF (allocok.NE.0) GOTO 1001
                     vectsin(:)=FREAL(0.0)
!
                     flagcfg=3
                     CALL readcfgobs (nam_groupoper(jgroup),flagcfg, &
     &                    kposcoefobs=poscoefobs(:,:))
!                     print *,poscoefobs(1:10,:)
                     flagcfg=2
                     CALL readcfgobs (nam_groupoper(jgroup),flagcfg, &
     &                    kgridijkobs=gridijkobs(:))
                     flagcfg=1
                     CALL readcfgobs (nam_groupoper(jgroup),flagcfg, &
     &                    kvectorms=vectorms(:))
                     CALL readvalobs(nam_groupoper(jgroup),vectsin(:))
!
                     DO jobs=1,obsend
                        indobs = obs_ord(jobs)
                        inddbs = obsnord(jobs)
                        jodeb = obs_ind(indobs,inddbs)
                        jofin = jodeb - 1 + obs_nbr(indobs,inddbs)
                        jodebnew = obs_indnew(indobs,inddbs) +  &
     &                       obs_nbrtot(indobs,inddbs)
                        jofinnew = jodebnew - 1 +  &
     &                          obs_nbr(indobs,inddbs)
!
                        IF (obs_nbr(indobs,inddbs).NE.0) THEN
                           vectonew(jodebnew:jofinnew) = &
     &                          vectsin(jodeb:jofin)
                           vectormsnew(jodebnew:jofinnew) = &
     &                          vectorms(jodeb:jofin)
                           DO jitp=1,jpitpend
                              poscoefobsnew(jodebnew:jofinnew,jitp) = &
     &                             poscoefobs(jodeb:jofin,jitp)
                           ENDDO
                           gridijkobsnew(jodebnew:jofinnew) = &
     &                          gridijkobs(jodeb:jofin)
                        ENDIF
                        obs_nbrtot(indobs,inddbs) = obs_nbrtot(indobs, &
     &                       inddbs)+obs_nbr(indobs,inddbs)
                     ENDDO
!
                     IF (allocated(vectsin)) deallocate (vectsin)
                     IF (allocated(vectorms)) deallocate(vectorms)
                     IF (allocated(poscoefobs)) deallocate(poscoefobs)
                     IF (allocated(gridijkobs)) deallocate(gridijkobs)
!
                  ENDDO
!
                  jpoend=jpoendnew
                  jpitpend=jpitpsizenew
                  obs_nbr(:,:)=obs_nbrnew(:,:)
                  obs_ind(:,:)=obs_indnew(:,:)
                  obs_itp(:,:)=obs_itpnew(:,:) 
                  CALL writeobs(koutxyo,vectonew(:),vectormsnew(:), &
     &              gridijkobsnew(:),poscoefobsnew(:,:))
!
                  obs_nbr(:,:)=obs_nbrold(:,:)
                  obs_ind(:,:)=obs_indold(:,:)
                  jpoend=jpoendold
                  jpitpsize=jpitpsizeold
!
                  IF (allocated(obs_nbrnew)) deallocate(obs_nbrnew)
                  IF (allocated(obs_nbrold)) deallocate(obs_nbrold)
                  IF (allocated(obs_nbrjgrp)) deallocate(obs_nbrjgrp)
                  IF (allocated(vectonew)) deallocate (vectonew)
                  IF (allocated(vectormsnew)) deallocate(vectormsnew)
                  IF (allocated(poscoefobsnew))  &
     &                 deallocate(poscoefobsnew)
                  IF (allocated(gridijkobsnew))  &
     &                 deallocate(gridijkobsnew)
!     
               ENDIF
               CASE ('jday2date')
                  OPEN(UNIT=11,file=kincfg)
                  READ(11,*) jdaytarget
                  CLOSE(11)

                  jday=-1

                  year=1950
                  DO WHILE (jday.LT.jdaytarget)
                    jday=jday+daysinyear(year)
                    year=year+1
                  ENDDO
                  year=year-1
                  jday=jday-daysinyear(year)

                  month=1
                  DO WHILE (jday.LT.jdaytarget)
                    jday=jday+daysinmonth(year,month)
                    month=month+1
                  ENDDO
                  month=month-1
                  jday=jday-daysinmonth(year,month)

                  day=jdaytarget-jday

                  PRINT '(a5,i4.4,i2.2,i2.2)','DATE:',year,month,day
               CASE DEFAULT
                  GOTO 101
               END SELECT
            CASE DEFAULT
               GOTO 101
            END SELECT
!
         CASE (1)
!     
            SELECT CASE (posit(ktextoper,'_'))
            CASE (2:bgword)        
!
               READ(ktextoper((posit(ktextoper,'_')+1): &
     &              lenv(ktextoper)),*,IOSTAT=ios) cst
               IF (ios.NE.0) GOTO 101
!     
               SELECT CASE (ktextoper(1:(posit(ktextoper,'_')-1)))
               CASE ('+')
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = vectsin(js) + cst
                     ENDIF
                  ENDDO    
               CASE ('-')
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = vectsin(js) - cst
                     ENDIF
                  ENDDO    
               CASE ('x')
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = vectsin(js) * cst
                     ENDIF
                  ENDDO            
               CASE ('/')
                  IF (cst.EQ.FREAL(0.0)) GOTO 101
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = vectsin(js) / cst
                     ENDIF
                  ENDDO    
               CASE ('min')
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = MIN( vectsin(js) , cst )
                     ENDIF
                  ENDDO
               CASE ('max')
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = MAX( vectsin(js) , cst )
                     ENDIF
                  ENDDO
               CASE ('inv')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.FREAL(0.0)) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = cst / vectsin(js)
                     ENDIF
                  ENDDO    
               CASE ('pow')
                  IF (cst.LE.FREAL(0.0)) GOTO 101
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = vectsin(js) ** cst
                     ENDIF
                  ENDDO    
               CASE ('cst')
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = cst
                     ENDIF
                  ENDDO
               CASE ('round')
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        x=vectsin(js)-INT(vectsin(js))
                        x=ANINT(x/cst)*cst
                        x=x+INT(vectsin(js))
                        vectsout(js) = x
                     ENDIF
                  ENDDO
               CASE ('mask')
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        IF (vectsin(js).EQ.0.0) THEN
                          vectsout(js) = cst
                        ELSE
                          vectsout(js) = vectsin(js)
                        ENDIF
                     ENDIF
                  ENDDO
               CASE ('zerocst')
                  DO js=1,jpssize
                     IF (vectsin(js).LT.cst) THEN
                        vectsout(js) = 0.0
                     ELSE
                        vectsout(js) = vectsin(js)
                     ENDIF
                  ENDDO
!              CASE ('masknan')
!                 DO js=1,jpssize
!                    IF (isnan(vectsin(js))) THEN
!                       vectsout(js) = cst
!                    ELSE
!                       vectsout(js) = vectsin(js)
!                    ENDIF
!                 ENDDO
               CASE ('maskcst')
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.cst) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = vectsin(js)
                     ENDIF
                  ENDDO
               CASE ('csp')
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = cst
                     ELSE
                        vectsout(js) = vectsin(js)
                     ENDIF
                  ENDDO            
               CASE ('gnoise')
                  CALL kiss_load()
                  DO js=1,jpssize
                     CALL kiss_gaussian(gran)
                     IF (vectsin(js).EQ.spvalsin) THEN
                        vectsout(js) = spvals
                     ELSE
                        IF (cst.GT.FREAL(0.0)) THEN
                           vectsout(js) = gran * cst
                        ELSE
                           vectsout(js) = gran * ABS(vectsin(js))
                        ENDIF
                     ENDIF
                  ENDDO            
                  CALL kiss_save()
               CASE ('quantiles')
                  vectsout = vectsin
                  CALL heapsort(vectsout)
                  OPEN(UNIT=11,file='quantiles.txt')
                  DO jqua=0,NINT(cst)
                    js = NINT(REAL(jqua,8)*REAL(jpssize-1,8)/cst) + 1
                    WRITE(11,*) vectsout(js)
                  ENDDO
                  CLOSE(11)
               CASE ('quantize')
                  OPEN(UNIT=11,file='quantiles.txt')
                  READ(11,*)
                  vectsout(:) = 1.
                  DO jqua=1,NINT(cst)-1
                    READ(11,*) qua
                    DO js=1,jpssize
                      IF (vectsin(js).GT.qua) THEN
                        vectsout(js) = jqua+1 !  REAL(jqua,8)/(NINT(cst)-1)
                      ENDIF
                    ENDDO
                  ENDDO
                  CLOSE(11)
               CASE DEFAULT
                  GOTO 101
               END SELECT
            CASE (0)
               SELECT CASE (ktextoper(1:lenv(ktextoper)))
               CASE ('stat')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  vectsout(:) = vectsin(:)*vectsin(:)
                  print *, 'MIN in:', MINVAL(vectsin(:))
                  print *, 'MAX in:', MAXVAL(vectsin(:))
                  print *, 'Root mean square:', &
     &             SQRT(SUM(vectsout(:))/jpssize)
               CASE ('pow2')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  vectsout(:) = vectsin(:) ** 2
               CASE ('sqrt')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  DO js=1,jpssize
                     IF (vectsin(js).LT.FREAL(0.0)) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = SQRT ( vectsin(js) )
                     ENDIF
                  ENDDO
               CASE ('abs')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  vectsout(:) = ABS ( vectsin(:) )
               CASE ('opp')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  vectsout(:) = - vectsin(:) 
               CASE ('inv')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  DO js=1,jpssize
                     IF (vectsin(js).EQ.FREAL(0.0)) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) =  FREAL(1.0) / vectsin(js)
                     ENDIF
                  ENDDO
               CASE ('log')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  DO js=1,jpssize
                     IF (vectsin(js).LE.FREAL(0.0)) THEN
                        vectsout(js) = spvals
                     ELSE
                        vectsout(js) = LOG ( vectsin(js) )
                     ENDIF
                  ENDDO
               CASE ('log10')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  DO js=1,jpssize
                     IF (vectsin(js).LE.FREAL(0.0)) THEN
                        vectsout(js) = - HUGE ( vectsin(js) )
                     ELSE
                        vectsout(js) = LOG10 ( vectsin(js) )
                     ENDIF
                  ENDDO
               CASE ('exp')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  vectsout(:) = EXP ( vectsin(:) )
               CASE ('entropy')
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  DO js=1,jpssize
                     x = vectsin(js)
                     IF ( (x.LT.FREAL(0.0)).OR.(x.GT.FREAL(1.0)) ) THEN
                        vectsout(js) = spvals
                     ELSE
                        IF ( (x.EQ.FREAL(0.0)).OR.(x.EQ.FREAL(1.0)) ) THEN
                          vectsout(js) = 0.0_kr
                        ELSE
                          vectsout(js) = - x * LOG(x) / LOG(2.0_kr)
                          x = 1.0_kr - x
                          vectsout(js) = - x * LOG(x) / LOG(2.0_kr) + vectsout(js)
                        ENDIF
                     ENDIF
                  ENDDO
               CASE ('GPRODtoGAU')
                  jpisize=10000
                  IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
                  IF (jpssize.GT.jpisize) THEN
                    allocate ( garray(0:jpisize), stat=allocok )
                    IF (allocok.NE.0) GOTO 1001
!                   ! Tabulate function
                    gpmin=minval(vectsin(:))
                    gpmax=maxval(vectsin(:))
                    DO ji=0,jpisize
                      gpval=gpmin+ji*(gpmax-gpmin)/jpisize
                      CALL gprod_to_gau(garray(ji),gpval)
                    ENDDO
!                   ! Interpolate function in table
                    DO js=1,jpssize
                      ri = jpisize*(vectsin(js)-gpmin)/(gpmax-gpmin)
                      ji = INT(ri)
                      IF (ji.EQ.jpisize) ji=ji-1
                      ri = ri - ji
                      vectsout(js) = garray(ji) + ri * ( garray(ji+1) - garray(ji) )
                    ENDDO
                    deallocate(garray)
                  ELSE
                    DO js=1,jpssize
                     CALL gprod_to_gau(vectsout(js),vectsin(js))
                    ENDDO
                  ENDIF
               CASE DEFAULT
                  GOTO 101
               END SELECT
            CASE DEFAULT
               GOTO 101
            END SELECT
         CASE (2)
            SELECT CASE (ktextoper(1:lenv(ktextoper)))
            CASE ('rstopa')
               IF (kflaganlxyo.NE.1) GOTO 104
               IF (.NOT.(validextvar(koutxyo))) GOTO 104
               jextvar = indext(koutxyo,extvartab,nbextvar)
               IF (jextvar.NE.8) GOTO 104
! writing of Xa in the restart at position now (the normally position)
               wribefore=0
               CALL writevar(koutxyo,vectsin(:),jnxyo)     
! writing of Xf in the restart at position before  (no initialization of the record           
               wribefore=1
               CALL writevar(koutxyo,vectsinref(:),jnxyo)
               wribefore=0
            CASE ('+')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = vectsin(:) + vectsinref(:)
            CASE ('-')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = vectsin(:) - vectsinref(:)
            CASE ('stat-')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = vectsin(:) - vectsinref(:)
               print *, 'MIN in:', MINVAL(vectsin(:))
               print *, 'MAX in:', MAXVAL(vectsin(:))
               print *, 'MIN inref:', MINVAL(vectsinref(:))
               print *, 'MAX inref:', MAXVAL(vectsinref(:))
               print *, 'MIN difference:', MINVAL(vectsout(:))
               print *, 'MAX difference:', MAXVAL(vectsout(:))
               print *, 'RMS difference:', &
     &             SQRT(SUM(vectsout(:)*vectsout(:))/jpssize)
            CASE ('min')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = MIN(vectsin(:),vectsinref(:))
            CASE ('max')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = MAX(vectsin(:),vectsinref(:))
            CASE ('x')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = vectsin(:) * vectsinref(:)
            CASE ('gt')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               WHERE(vectsin(:).GT.vectsinref(:)) vectsout(:) = 1.0_kr
            CASE ('dotproduct')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = vectsin(:) * vectsinref(:)
               print *, SUM(vectsout(:))
            CASE ('project')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = vectsinref(:) * vectsinref(:)
               norm = SUM(vectsout(:))
               vectsout(:) = vectsin(:) * vectsinref(:)
               print *, SUM(vectsout(:))/norm
            CASE ('/')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               DO js=1,jpssize
                  IF (vectsinref(js).EQ.FREAL(0.0)) THEN
                     vectsout(js) = spvals
                  ELSE
                     vectsout(js) = vectsin(js) / vectsinref(js)
                  ENDIF
               ENDDO
            CASE ('/sqrt')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               DO js=1,jpssize
                  IF (vectsinref(js).LE.FREAL(0.0)) THEN
                     vectsout(js) = spvals
                  ELSE
                     vectsout(js) = vectsin(js) / SQRT(vectsinref(js))
                  ENDIF
               ENDDO
            CASE ('sqrt<a2-b2>')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = vectsin(:) * vectsin(:) 
               vectsout(:) = vectsout(:) - vectsinref(:) * vectsinref(:)
               DO js=1,jpssize
                  vectsout(js) = MAX(FREAL(0.0),vectsout(js))
               ENDDO
               vectsout(:) = SQRT(vectsout(:))
            CASE ('sqrt<a2+b2>')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = vectsin(:) * vectsin(:)
               vectsout(:) = vectsout(:) + vectsinref(:) * vectsinref(:)
               vectsout(:) = SQRT(vectsout(:))
            CASE ('sqrt<a2-b2>1')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               vectsout(:) = vectsin(:) * vectsin(:)
               vectsin(:) = vectsinref(:) * vectsinref(:)
               vectsout(:) = vectsout(:) - vectsin(:)
               vectsout(:) = SQRT(MIN( FREAL(100.0)*vectsin(:), &
     &              MAX( FREAL(0.01)*vectsin(:),vectsout(:) )))
            CASE ('replace_spval')
               IF (ANY(vectsinref(:).EQ.spvalsin)) GOTO 103
               DO js=1,jpssize
                  IF (vectsin(js).NE.spvals) THEN
                     vectsout(js) = vectsin(js)
                  ELSE
                     vectsout(js) = vectsinref(js)
                  ENDIF
               ENDDO
            CASE ('GAUtoGAM')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               DO js=1,jpssize
                 nu = vectsinref(js)  ! std / mean
                 gammak = 1.0 / (nu * nu)
                 gammath = nu *nu    ! set theta = 1/k -> mean = 1
                 CALL gau_to_gam(vectsout(js),vectsin(js),gammak)  ! -> gamma(k,1)
                 vectsout(js) = vectsout(js) * gammath
               ENDDO
            CASE ('GAUtoBETA')
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               DO js=1,jpssize
                mu = vectsinref(js)
                IF (mu.EQ.0.) THEN
                  vectsout(js) = 0.
                ELSEIF (mu.EQ.1.) THEN
                  vectsout(js) = 1.
                ELSE
                  CALL gau_to_beta(vectsout(js),vectsin(js),mu,0.166_kr)
                ENDIF
               ENDDO
            CASE ('lognormal')
               obserror_type='lognormal'

               CALL kiss_load()
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsin(:).LE.0.)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               IF (ANY(vectsinref(:).LT.0.)) GOTO 103
               DO js=1,jpssize
                  CALL obserror_sample(vectsinref(js),vectsin(js),vectsout(js))
               ENDDO
               CALL kiss_save()
            CASE ('gamma')
               CALL kiss_load()
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsin(:).LE.0.)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               IF (ANY(vectsinref(:).LT.0.)) GOTO 103
               DO js=1,jpssize
                 mu = vectsinref(js) ; nu = vectsin(js) ! mu = mean, nu = std /mean
                 gammak = 1.0 / (nu * nu)
                 gammath = mu * nu *nu
                 CALL kiss_gamma(gran,gammak)
                 gran = gran * gammath
                 vectsout(js) = gran
               ENDDO
               CALL kiss_save()
            CASE ('beta')
               CALL kiss_load()
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsin(:).LE.0.)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               IF (ANY(vectsinref(:).LT.0.)) GOTO 103
               IF (ANY(vectsinref(:).GT.1.)) GOTO 103
               DO js=1,jpssize
                  mu = vectsinref(js) ; nu = vectsin(js)
                  IF (mu.EQ.0.) THEN
                    gran = 0.
                  ELSEIF (mu.EQ.1.) THEN
                    gran = 1.
                  ELSE
                    betaa = mu * nu
                    betab = (1.0 - mu) * nu
                    CALL kiss_beta(gran,betaa,betab)
                  ENDIF
                  vectsout(js) = gran
               ENDDO
               CALL kiss_save()
            CASE ('locerror')
               ! Read grid (assuming only one variable)
               jsxy=1
               IF (sxyend.NE.1) GOTO 1000

               jpisize=sxy_jpi(1)
               jpjsize=sxy_jpj(1)
!
               allocate (gridij(1:jpisize,1:jpjsize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               gridij(:,:) = type_gridij(FREAL(0.0),FREAL(0.0))

               CALL readgrd(kflaganlxyo,sxy_ord(1))

               allocate (lon(1:jpssize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001
               allocate (lat(1:jpssize), stat=allocok )
               IF (allocok.NE.0) GOTO 1001

               jk=1 ; jt=1
               CALL mk8vct(lon(1:),gridij(:,:)%longi,jk,jt, &
     &                     jsxy,nbr,kflaganlxyo)
               CALL mk8vct(lat(1:),gridij(:,:)%latj, jk,jt, &
     &                     jsxy,nbr,kflaganlxyo)
               lon = lon * deg2rad ; lat = (90. - lat) * deg2rad

               ! Compute location error
               IF (ANY(vectsin(:).EQ.spvalsin)) GOTO 103
               IF (ANY(vectsinref(:).EQ.spvalsinref)) GOTO 103
               DO js=1,jpssize
                 distmin = HUGE(nu)
                 DO js1=1,jpssize
                   IF (NINT(vectsin(js1)).EQ.NINT(vectsinref(js))) THEN
                     !dist=max(abs(lon(js1)-lon(js)),abs(lat(js1)-lat(js)))
                     dist=sph_dist(lon(js1),lon(js),lat(js1),lat(js))
                     IF (dist.LT.distmin) distmin = dist
                   ENDIF
                 ENDDO
                 vectsout(js) = distmin * earthrad
               ENDDO
            CASE DEFAULT
               GOTO 101
            END SELECT
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
         SELECT CASE (ktextoper(1:lenv(ktextoper)))
         CASE ('rstopa','grpobs','jday2date')
! --- Nothing
         CASE DEFAULT
! -3- Write out vector :
! ----------------------
!
            IF (nprint.GE.2) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'ALGOOPER : ', &
     &              'writing the out vector'
            ENDIF
            SELECT CASE (kflaganlxyo)
            CASE (1)
               CALL writevar(koutxyo,vectsout(:),jnxyo)
            CASE (2)
               CALL writedta(koutxyo,vectsout(:))
            CASE (3)
               IF (ANY(vectsout(:).EQ.spvals)) GOTO 102
               CALL writeobs(koutxyo,vectsout(:),vectorms(:), &
     &              gridijkobs(:),poscoefobs(:,:))
            CASE DEFAULT
               GOTO 1000
            END SELECT
         END SELECT
!
      ENDDO
!
! --- deallocation
      IF (allocated(nam_groupoper)) deallocate (nam_groupoper)
      IF (allocated(vectsout)) deallocate (vectsout)
      IF (allocated(vectsin)) deallocate (vectsin)
      IF (allocated(vectsinref)) deallocate (vectsinref)
      IF (allocated(vectorms)) deallocate (vectorms)
      IF (allocated(poscoefobs)) deallocate (poscoefobs)
      IF (allocated(gridijkobs)) deallocate (gridijkobs)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'algooper','algooper')
 1001 CALL printerror2(0,1001,3,'algooper','algooper')
 1003 CALL printerror2(0,1003,3,'algooper','algooper')
!
 101  WRITE (texterror,*) 'operation not allowed : ', &
     &     ktextoper(1:lenv(ktextoper))
      CALL printerror2(0,101,3,'algooper','algooper',comment=texterror)
!
 102  WRITE (texterror,*) 'spval operation not allowed for obs: ', &
     &     ktextoper(1:lenv(ktextoper))
      CALL printerror2(0,102,3,'algooper','algooper',comment=texterror)
!
 103  WRITE (texterror,*) 'spval operation not allowed : ', &
     &     ktextoper(1:lenv(ktextoper))
      CALL printerror2(0,103,3,'algooper','algooper',comment=texterror)
!
 104  WRITE (texterror,*) ktextoper(1:lenv(ktextoper)), &
     &     ' is a specific action for building restart opa'
      CALL printerror2(0,104,3,'algooper','algooper',comment=texterror)
!
 105  WRITE (texterror,*) 'operation not allowed for this file : ', &
     &     kincfg(1:lenv(kincfg))
      CALL printerror2(0,105,3,'algooper','algooper',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algooperzon(koutz, &
     &     kflagnbz,kincfg,kinz,kinrefz, &
     &     ktextoper,kconfigz)
!---------------------------------------------------------------------
!
!  Purpose : Tools for arithmetic operation on state and basis :
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
      use mod_spacexyo , only : jpyend,jpz, &
     &     pt1bubidx, pt2bubidx, pt3bubidx, pt4bubidx, &
     &     pt1dtalon, pt1dtalat, pt1dtadepth, pt1dtatime, &
     &     pt1bublon, pt1bublat, pt1bubdepth, pt1bubtime, &
     &     bubblk2, bubblk3, bubblk4, vectptbub
      use hiozon
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      INTEGER, intent(in) :: kflagnbz
      CHARACTER(len=*), intent(in) :: koutz, &
     &     kincfg,kinz,kinrefz, &
     &     ktextoper,kconfigz
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jbub,jz,jdta,jbub2,jbub3,jbub4
      INTEGER :: flagcfg,ios,jpgroup,jgroup,jpzsize
      LOGICAL :: lectinfo,lmodprint,zlectinfo
      CHARACTER(len=bgword), dimension(:), allocatable :: nam_groupoper
      BIGREAL :: cst
      INTEGER :: kjpz,kzon_jpi,kzon_jpj,kzon_jpk,kzon_jpt
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine algooperzon &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      IF (nprint.GE.1) print *, &
     &     '--> Computation of operation ...'
      lectinfo = .FALSE.
      jpbubend(:)=0
!
! -1.0- Read header zon :
! -----------------------
!
! -1.0.1- Read header config.zon :
! --------------------------------
!
      CALL evalhdrzon(kconfigz,zon_jpi,zon_jpj,zon_jpk,zon_jpt, &
     &     jpbubend(1),jpz)
!
! -1.0.2- Read header inz :
! -------------------------
!
      IF (kflagnbz.GE.1) THEN
         IF (nprint.GE.2) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOOPERZON : ', &
     &           'Reading inzon vector'
         ENDIF
!
         CALL evalhdrzon(kinz,kzon_jpi,kzon_jpj,kzon_jpk,kzon_jpt, &
     &        jpbubend(2),kjpz)
         IF (kzon_jpi.NE.zon_jpi) GOTO 102
         IF (kzon_jpj.NE.zon_jpj) GOTO 102
         IF (kzon_jpk.NE.zon_jpk) GOTO 102
         IF (kzon_jpt.NE.zon_jpt) GOTO 102
         IF (kjpz.NE.jpz) GOTO 102
!     
      ENDIF
!
! -1.0.3- Read header inrefz :
! -------------------------
!
      IF (kflagnbz.EQ.2) THEN
         IF (nprint.GE.2) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOOPERZON : ', &
     &           'Reading inzonref vector'
         ENDIF
!
         CALL evalhdrzon(kinrefz,kzon_jpi,kzon_jpj,kzon_jpk,kzon_jpt, &
     &        jpbubend(3),kjpz)
         IF (kzon_jpi.NE.zon_jpi) GOTO 103
         IF (kzon_jpj.NE.zon_jpj) GOTO 103
         IF (kzon_jpk.NE.zon_jpk) GOTO 103
         IF (kzon_jpt.NE.zon_jpt) GOTO 103
         IF (kjpz.NE.jpz) GOTO 103
!
      ENDIF
!
! -1.1- Read pointer zon :
! -----------------------
!
      jpzsize=jpz
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
! -1.1.1- Read pointer config.zon :
! ---------------------------------
!
      allocate ( pt1bubidx(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt1bubidx(:,:) = 0
!
      CALL readptzon (kconfigz,pt1bubidx,pt1dtalon,pt1dtalat, &
     &     pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth, &
     &     pt1bubtime)
!
! -1.1.2- Read pointer inz vector :
! -------------------------------
!
      IF (kflagnbz.GE.1) THEN
!     
         allocate ( pt2bubidx(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         pt2bubidx(:,:) = 0
!
         CALL readptzon (kinz,pt2bubidx,pt1dtalon,pt1dtalat, &
     &        pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth, &
     &        pt1bubtime)
!     
      ENDIF
!
! -1.1.3- Read inrefz vector :
! -----------------------------
!
      IF (kflagnbz.EQ.2) THEN
!     
         allocate ( pt3bubidx(1:jpzsize,1:dtaend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         pt3bubidx(:,:) = 0
!     
         CALL readptzon (kinrefz,pt3bubidx,pt1dtalon,pt1dtalat, &
     &        pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth, &
     &        pt1bubtime)
!     
      ENDIF
!
! -1.2- Read bub zon :
! --------------------
!
! -1.1.2- Read pointer inz vector :
! -------------------------------
!
      IF (kflagnbz.GE.1) THEN
!     
! --- allocation vectptbub
         allocate ( vectptbub(1:MAXVAL(jpbubend)), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectptbub(:) = 0
!
         allocate ( bubblk2(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &        1:jpbubend(2),1:1), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         bubblk2(:,:,:,:,:,:) = FREAL(0.0)
!     
         zlectinfo = .FALSE.
         vectptbub(1:jpbubend(2)) = (/ (jbub, jbub = 1,jpbubend(2)) /)
         CALL readnbubzon(kinz,vectptbub(1:jpbubend(2)), &
     &        bubblk2(:,:,:,:,1:jpbubend(2),1),zlectinfo)
!
      ENDIF
!
! -1.1.3- Read inrefz vector :
! -----------------------------
!
      IF (kflagnbz.EQ.2) THEN
!
         allocate ( bubblk3(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &        1:jpbubend(3),1:1), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         bubblk3(:,:,:,:,:,:) = FREAL(0.0)
!     
         zlectinfo = .FALSE.
         vectptbub(1:jpbubend(3)) = (/ (jbub, jbub = 1,jpbubend(3)) /)
         CALL readnbubzon(kinrefz,vectptbub(1:jpbubend(3)), &
     &        bubblk3(:,:,:,:,1:jpbubend(3),1),zlectinfo)
!     
      ENDIF
!
      IF (allocated(vectptbub)) deallocate (vectptbub)
!
! -1.2- Write outz vector :
! -------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ALGOOPERZON : ', &
     &        'writing the outzon vector'
      ENDIF
!     
      allocate ( pt4bubidx(1:jpzsize,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pt4bubidx(:,:) = 0
!
      SELECT CASE(kflagnbz)
      CASE (0)
         jpbubend(4)=dtaend
         DO jdta=1,dtaend
            jbub=jdta
            pt4bubidx(:,jdta)=jbub
         ENDDO
      CASE (1)
         jpbubend(4)=jpbubend(2)
         pt4bubidx(:,:)=pt2bubidx(:,:)
      CASE (2)
         IF ((jpbubend(2).NE.1).AND.jpbubend(3).NE.1) THEN
               jpbubend(4)=jpzsize*dtaend
               jbub = 0
               DO jz=1,jpz
               DO jdta=1,dtaend
                  jbub=jbub + 1
                  pt4bubidx(jz,jdta)=jbub
               ENDDO
               ENDDO
         ELSE
            IF (jpbubend(2).EQ.1) THEN 
               jpbubend(4)=jpbubend(3)
               pt4bubidx(:,:)=pt3bubidx(:,:)
            ELSE
               jpbubend(4)=jpbubend(2)
               pt4bubidx(:,:)=pt2bubidx(:,:)
            ENDIF           
         ENDIF
      CASE DEFAULT
            GOTO 1000
      END SELECt
!
! -1.2.1- Write header outz vector :
! ----------------------------------
!
      CALL writehdrzon (koutz,zon_jpi,zon_jpj,zon_jpk, &
     &     zon_jpt,jpbubend(4),jpz)
!
! -1.2.2- Write pointer outz vector :
! ----------------------------------
!
      CALL writeptzon (koutz,pt4bubidx,pt1dtalon,pt1dtalat, &
     &     pt1dtadepth,pt1dtatime,pt1bublon,pt1bublat,pt1bubdepth, &
     &     pt1bubtime)
!
      IF (allocated(pt1dtalon)) deallocate (pt1dtalon)
      IF (allocated(pt1dtalat)) deallocate (pt1dtalat)
      IF (allocated(pt1dtadepth)) deallocate (pt1dtadepth)
      IF (allocated(pt1dtatime)) deallocate (pt1dtatime)
      IF (allocated(pt1bublon)) deallocate (pt1bublon)
      IF (allocated(pt1bublat)) deallocate (pt1bublat)
      IF (allocated(pt1bubdepth)) deallocate (pt1bubdepth)
      IF (allocated(pt1bubtime)) deallocate (pt1bubtime)
!
      allocate ( bubblk4(1:zon_jpi,1:zon_jpj,1:zon_jpk,1:zon_jpt, &
     &     1:jpbubend(4),1:1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      bubblk4(:,:,:,:,:,:) = FREAL(0.0)
!     
! -2- Operation :
! -----------------
!
         SELECT CASE (kflagnbz)
         CASE (0)
            SELECT CASE (posit(ktextoper,'_'))
            CASE (2:bgword)        
               READ(ktextoper((posit(ktextoper,'_')+1): &
     &                 lenv(ktextoper)),*,IOSTAT=ios) cst
               IF (ios.NE.0) GOTO 101
               GOTO 101
            CASE (0)        
               SELECT CASE (ktextoper(1:lenv(ktextoper)))
               CASE DEFAULT
                  GOTO 101
               END SELECT
            CASE DEFAULT
               GOTO 101
            END SELECT
!
         CASE (1)
!     
            SELECT CASE (posit(ktextoper,'_'))
            CASE (2:bgword)        
!
               READ(ktextoper((posit(ktextoper,'_')+1): &
     &              lenv(ktextoper)),*,IOSTAT=ios) cst
               IF (ios.NE.0) GOTO 101
!     
               SELECT CASE (ktextoper(1:(posit(ktextoper,'_')-1)))
               CASE ('+')
                  bubblk4 (:,:,:,:,:,:) = bubblk2 (:,:,:,:,:,:) + cst
               CASE ('-')
                  bubblk4 (:,:,:,:,:,:) = bubblk2 (:,:,:,:,:,:) - cst
               CASE ('x')         
                  bubblk4 (:,:,:,:,:,:) = bubblk2 (:,:,:,:,:,:) * cst
               CASE ('/')
                  IF (cst.EQ.FREAL(0.0)) GOTO 101
                  bubblk4 (:,:,:,:,:,:) = bubblk2 (:,:,:,:,:,:) / cst
               CASE ('inv') 
                  IF (ANY(bubblk2 (:,:,:,:,:,:).EQ.FREAL(0.0))) GOTO 104
                  bubblk4 (:,:,:,:,:,:) = cst / bubblk2 (:,:,:,:,:,:) 
               CASE ('pow')
                  IF (cst.LE.FREAL(0.0)) GOTO 101
                  bubblk4 (:,:,:,:,:,:) = bubblk2 (:,:,:,:,:,:) ** cst
               CASE ('cst')   
                  bubblk4 (:,:,:,:,:,:) = cst
               CASE DEFAULT
                  GOTO 101
               END SELECT
            CASE (0)
               SELECT CASE (ktextoper(1:lenv(ktextoper)))
               CASE ('pow2')
                  bubblk4 (:,:,:,:,:,:) = bubblk2 (:,:,:,:,:,:) ** 2
               CASE ('sqrt')
                  IF (ANY(bubblk2 (:,:,:,:,:,:).LT.FREAL(0.0))) GOTO 104
                  bubblk4 (:,:,:,:,:,:) = SQRT ( bubblk2 (:,:,:,:,:,:) )
               CASE ('abs')
                  bubblk4 (:,:,:,:,:,:) = ABS ( bubblk2 (:,:,:,:,:,:) )
               CASE ('opp')
                  bubblk4 (:,:,:,:,:,:) = - bubblk2 (:,:,:,:,:,:) 
               CASE ('inv')
                  IF (ANY(bubblk2 (:,:,:,:,:,:).EQ.FREAL(0.0))) GOTO 104
                  bubblk4 (:,:,:,:,:,:) = FREAL(1.0) / bubblk2 (:,:,:,:,:,:)
               CASE DEFAULT
                  GOTO 101
               END SELECT
            CASE DEFAULT
               GOTO 101
            END SELECT
         CASE (2)
            SELECT CASE (ktextoper(1:lenv(ktextoper)))
            CASE ('+')
         IF ((jpbubend(2).NE.1).AND.jpbubend(3).NE.1) THEN
            DO jz=1,jpz
            DO jdta=1,dtaend
               jbub2=pt2bubidx(jz,jdta)
               jbub3=pt3bubidx(jz,jdta)
               jbub4=pt4bubidx(jz,jdta)
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) + &
     &              bubblk3 (:,:,:,:,jbub3,:)
            ENDDO
            ENDDO
         ELSE
            DO jbub=1,jpbubend(4)
               IF (jpbubend(2).EQ.1) THEN
                  jbub2=1
               ELSE
                  jbub2=jbub
               ENDIF
               IF (jpbubend(3).EQ.1) THEN
                  jbub3=1
               ELSE
                  jbub3=jbub
               ENDIF
               jbub4=jbub
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) + &
     &              bubblk3 (:,:,:,:,jbub3,:)
            ENDDO
         ENDIF
            CASE ('-')
         IF ((jpbubend(2).NE.1).AND.jpbubend(3).NE.1) THEN
            DO jz=1,jpz
            DO jdta=1,dtaend
               jbub2=pt2bubidx(jz,jdta)
               jbub3=pt3bubidx(jz,jdta)
               jbub4=pt4bubidx(jz,jdta)
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) - &
     &              bubblk3 (:,:,:,:,jbub3,:)
            ENDDO
            ENDDO
         ELSE
            DO jbub=1,jpbubend(4)
               IF (jpbubend(2).EQ.1) THEN
                  jbub2=1
               ELSE
                  jbub2=jbub
               ENDIF
               IF (jpbubend(3).EQ.1) THEN
                  jbub3=1
               ELSE
                  jbub3=jbub
               ENDIF
               jbub4=jbub
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) - &
     &              bubblk3 (:,:,:,:,jbub3,:)
            ENDDO
         ENDIF
            CASE ('min')
         IF ((jpbubend(2).NE.1).AND.jpbubend(3).NE.1) THEN
            DO jz=1,jpz
            DO jdta=1,dtaend
               jbub2=pt2bubidx(jz,jdta)
               jbub3=pt3bubidx(jz,jdta)
               jbub4=pt4bubidx(jz,jdta)
               bubblk4 (:,:,:,:,jbub4,:) = MIN(bubblk2 (:,:,:,:,jbub2,:) &
     &              ,bubblk3 (:,:,:,:,jbub3,:))
            ENDDO
            ENDDO
         ELSE
            DO jbub=1,jpbubend(4)
               IF (jpbubend(2).EQ.1) THEN
                  jbub2=1
               ELSE
                  jbub2=jbub
               ENDIF
               IF (jpbubend(3).EQ.1) THEN
                  jbub3=1
               ELSE
                  jbub3=jbub
               ENDIF
               jbub4=jbub
               bubblk4 (:,:,:,:,jbub4,:) = MIN(bubblk2 (:,:,:,:,jbub2,:) &
     &              ,bubblk3 (:,:,:,:,jbub3,:))
            ENDDO
         ENDIF
            CASE ('max')
         IF ((jpbubend(2).NE.1).AND.jpbubend(3).NE.1) THEN
            DO jz=1,jpz
            DO jdta=1,dtaend
               jbub2=pt2bubidx(jz,jdta)
               jbub3=pt3bubidx(jz,jdta)
               jbub4=pt4bubidx(jz,jdta)
               bubblk4 (:,:,:,:,jbub4,:) = MAX(bubblk2 (:,:,:,:,jbub2,:) &
     &              ,bubblk3 (:,:,:,:,jbub3,:))
            ENDDO
            ENDDO
         ELSE
            DO jbub=1,jpbubend(4)
               IF (jpbubend(2).EQ.1) THEN
                  jbub2=1
               ELSE
                  jbub2=jbub
               ENDIF
               IF (jpbubend(3).EQ.1) THEN
                  jbub3=1
               ELSE
                  jbub3=jbub
               ENDIF
               jbub4=jbub
               bubblk4 (:,:,:,:,jbub4,:) = MAX(bubblk2 (:,:,:,:,jbub2,:) &
     &              ,bubblk3 (:,:,:,:,jbub3,:))
            ENDDO
         ENDIF
            CASE ('x')
         IF ((jpbubend(2).NE.1).AND.jpbubend(3).NE.1) THEN
            DO jz=1,jpz
            DO jdta=1,dtaend
               jbub2=pt2bubidx(jz,jdta)
               jbub3=pt3bubidx(jz,jdta)
               jbub4=pt4bubidx(jz,jdta)
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) * &
     &              bubblk3 (:,:,:,:,jbub3,:)
            ENDDO
            ENDDO
         ELSE
            DO jbub=1,jpbubend(4)
               IF (jpbubend(2).EQ.1) THEN
                  jbub2=1
               ELSE
                  jbub2=jbub
               ENDIF
               IF (jpbubend(3).EQ.1) THEN
                  jbub3=1
               ELSE
                  jbub3=jbub
               ENDIF
               jbub4=jbub
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) *  &
     &              bubblk3 (:,:,:,:,jbub3,:)
            ENDDO
         ENDIF
            CASE ('/')
               IF (ANY(bubblk3 (:,:,:,:,:,:).EQ.FREAL(0.0))) GOTO 106
         IF ((jpbubend(2).NE.1).AND.jpbubend(3).NE.1) THEN
            DO jz=1,jpz
            DO jdta=1,dtaend
               jbub2=pt2bubidx(jz,jdta)
               jbub3=pt3bubidx(jz,jdta)
               jbub4=pt4bubidx(jz,jdta)
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) / &
     &              bubblk3 (:,:,:,:,jbub3,:)
            ENDDO
            ENDDO
         ELSE
            DO jbub=1,jpbubend(4)
               IF (jpbubend(2).EQ.1) THEN
                  jbub2=1
               ELSE
                  jbub2=jbub
               ENDIF
               IF (jpbubend(3).EQ.1) THEN
                  jbub3=1
               ELSE
                  jbub3=jbub
               ENDIF
               jbub4=jbub
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) /  &
     &              bubblk3 (:,:,:,:,jbub3,:)
            ENDDO
         ENDIF
            CASE ('/sqrt')
               IF (ANY(bubblk3 (:,:,:,:,:,:).LE.FREAL(0.0))) GOTO 106
         IF ((jpbubend(2).NE.1).AND.jpbubend(3).NE.1) THEN
            DO jz=1,jpz
            DO jdta=1,dtaend
               jbub2=pt2bubidx(jz,jdta)
               jbub3=pt3bubidx(jz,jdta)
               jbub4=pt4bubidx(jz,jdta)
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) / &
     &              SQRT(bubblk3 (:,:,:,:,jbub3,:))
            ENDDO
            ENDDO
         ELSE
            DO jbub=1,jpbubend(4)
               IF (jpbubend(2).EQ.1) THEN
                  jbub2=1
               ELSE
                  jbub2=jbub
               ENDIF
               IF (jpbubend(3).EQ.1) THEN
                  jbub3=1
               ELSE
                  jbub3=jbub
               ENDIF
               jbub4=jbub
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) /  &
     &              SQRT(bubblk3 (:,:,:,:,jbub3,:))
            ENDDO
         ENDIF
            CASE ('sqrt<a2-b2>')
         IF ((jpbubend(2).NE.1).AND.jpbubend(3).NE.1) THEN
            DO jz=1,jpz
            DO jdta=1,dtaend
               jbub2=pt2bubidx(jz,jdta)
               jbub3=pt3bubidx(jz,jdta)
               jbub4=pt4bubidx(jz,jdta)
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) *  &
     &              bubblk2 (:,:,:,:,jbub2,:)
               bubblk4 (:,:,:,:,jbub4,:) = bubblk4 (:,:,:,:,jbub4,:) -  &
     &              bubblk3 (:,:,:,:,jbub3,:)* bubblk3 (:,:,:,:,jbub3,:)
               bubblk4 (:,:,:,:,jbub4,:) = SQRT(MAX( FREAL(0.0) ,  &
     &              bubblk4 (:,:,:,:,jbub4,:) ))
            ENDDO
            ENDDO
         ELSE
            DO jbub=1,jpbubend(4)
               IF (jpbubend(2).EQ.1) THEN
                  jbub2=1
               ELSE
                  jbub2=jbub
               ENDIF
               IF (jpbubend(3).EQ.1) THEN
                  jbub3=1
               ELSE
                  jbub3=jbub
               ENDIF
               jbub4=jbub
               bubblk4 (:,:,:,:,jbub4,:) = bubblk2 (:,:,:,:,jbub2,:) *  &
     &              bubblk2 (:,:,:,:,jbub2,:)
               bubblk4 (:,:,:,:,jbub4,:) = bubblk4 (:,:,:,:,jbub4,:) -  &
     &              bubblk3 (:,:,:,:,jbub3,:)* bubblk3 (:,:,:,:,jbub3,:)
               bubblk4 (:,:,:,:,jbub4,:) = SQRT(MAX( FREAL(0.0) ,  &
     &              bubblk4 (:,:,:,:,jbub4,:) ))
            ENDDO
         ENDIF
            CASE DEFAULT
               GOTO 101
            END SELECT
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
! -3- Write out vector :
! ----------------------
!
         IF (nprint.GE.2) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'ALGOOPERZON : ', &
     &              'writing the outzon vector'
         ENDIF
         jbub=1
         CALL writenbubzon(koutz, &
     &        bubblk4(:,:,:,:,1:jpbubend(4),1),jbub)
!
! --- deallocation
      IF (allocated(nam_groupoper)) deallocate (nam_groupoper)
      IF (allocated(pt2bubidx)) deallocate (pt2bubidx)
      IF (allocated(pt3bubidx)) deallocate (pt3bubidx)
      IF (allocated(pt4bubidx)) deallocate (pt4bubidx)
      IF (allocated(bubblk2)) deallocate (bubblk2)
      IF (allocated(bubblk3)) deallocate (bubblk3)
      IF (allocated(bubblk4)) deallocate (bubblk4)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'algooper','algooperzon')
 1001 CALL printerror2(0,1001,3,'algooper','algooperzon')
!
 101  WRITE (texterror,*) 'operation not allowed : ', &
     &     ktextoper(1:lenv(ktextoper))
      CALL printerror2(0,101,3,'algooper','algooperzon', &
     &     comment=texterror)
!
 102  WRITE (texterror,*) 'file ',kinz(1:lenv(kinz)), &
     &     ' incoherent with ',kconfigz(1:lenv(kconfigz))
      CALL printerror2(0,102,3,'algooper','algooperzon', &
     &     comment=texterror)
!
 103  WRITE (texterror,*) 'file ',kinrefz(1:lenv(kinrefz)), &
     &     ' incoherent with ',kconfigz(1:lenv(kconfigz))
      CALL printerror2(0,103,3,'algooper','algooperzon', &
     &     comment=texterror)
!
 104  WRITE (texterror,*) ' bubblk2 from inzon ', &
     &     'contains null or negative values'
      CALL printerror2(0,104,3,'algooper','algooperzon', &
     &     comment=texterror)
!
 105  WRITE (texterror,*) 'operation not allowed for this file : ', &
     &     kincfg(1:lenv(kincfg))
      CALL printerror2(0,105,3,'algooper','algooperzon', &
     &     comment=texterror)
!
 106  WRITE (texterror,*) ' bubblk3 from inzonref ', &
     &     'contains null or negative values'
      CALL printerror2(0,106,3,'algooper','algooperzon', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION daysinyear(year)
!---------------------------------------------------------------------
!
!  Purpose : Get number of days in year
!  -------
!---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
      INTEGER, intent(in) :: year
      INTEGER :: daysinyear

      IF (MOD(year,4).EQ.0) THEN
        daysinyear = 366
      ELSE
        daysinyear = 365
      ENDIF

      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION daysinmonth(year,month)
!---------------------------------------------------------------------
!
!  Purpose : Get number of days in year
!  -------
!---------------------------------------------------------------------
      use mod_main
      IMPLICIT NONE
      INTEGER, intent(in) :: year,month
      INTEGER :: daysinmonth

      SELECT CASE (month)
      CASE (1,3,5,7,8,10,12)
        daysinmonth = 31
      CASE (4,6,9,11)
        daysinmonth = 30
      CASE (2)
        daysinmonth = 28
      END SELECT

      IF ((month.EQ.2).AND.(MOD(year,4).EQ.0)) THEN
        daysinmonth = 29
      ENDIF

      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algooper
