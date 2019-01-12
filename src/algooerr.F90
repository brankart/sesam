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
! ---                    ALGOOERR.F90                             ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 00-01 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE mkerrtodta
! --- SUBROUTINE mkoperstddta
! --- SUBROUTINE mkerrobastodta
! --- SUBROUTINE mkbiastodta
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algooerr
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC mkerrtodta,mkoperstddta,mkerrobastodta,mkbiastodta

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkerrtodta (kfnouterrdta,kfnincfg,kfninerrdta)
!---------------------------------------------------------------------
!
!  Purpose :
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
      use mod_spacexyo , only : jpyend,spvaldta
      use hioxyo
      use utilroa, only : mkyorms
      use utilmkh
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnouterrdta,kfnincfg
      CHARACTER(len=*), intent(in), optional :: kfninerrdta
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vecty
      BIGREAL, dimension(:), allocatable, save :: vectyrms
!
      INTEGER :: allocok,jpysize
      INTEGER :: flagcfg,flagxyo,jnxyo
      LOGICAL :: lectinfo
      INTEGER :: jdta,inddta,jy
      CHARACTER(len=bgword) kdtanam,comment
      BIGREAL :: kdtaval
      CHARACTER(len=hgword) :: text
      CHARACTER(len=1) :: textexclusion
!----------------------------------------------------------------------
!
      jpysize=jpyend
! --- allocation vectyrms
      allocate ( vectyrms(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectyrms(:) = FREAL(0.0)
! --- allocation vecty
      IF (present(kfninerrdta)) THEN
         allocate ( vecty(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vecty(:) = FREAL(0.0)
      ENDIF
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&         routine mkerrtodta               &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      lmoyect=.FALSE.
!     
! -1.- Compute the vectyrms file :
! ---------------------------------
!
      textexclusion='#'
!
      CALL openfile(numfil,kfnincfg)
!     
      DO jdta=1,dtaend
         inddta=dta_ord(jdta)
         text=''
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,ERR=102) kdtanam
!         print *, 'kdtanam',kdtanam(1:lenv(kdtanam)),
!     $        'dta_nam',dta_nam(inddta)(1:lenv(dta_nam(inddta)))
         IF (kdtanam(1:lenv(kdtanam)).NE. &
     &        (dta_nam(inddta)(1:lenv(dta_nam(inddta))))) GOTO 102
         text=''
         text=readnextline(numfil,textexclusion)
         READ(text,FMT=*,ERR=103) kdtaval
!         print *,'kdtaval',kdtaval
         IF (kdtaval.EQ.FREAL(0.0))  print *, 'Warning : dta_rms(', &
     &        dta_nam(inddta)(1:lenv(dta_nam(inddta))),')=',kdtaval
         dta_rms(inddta)=kdtaval
      ENDDO
!
      flagxyo=2
      CALL mkyorms (vectyrms,flagxyo)
!     
! -2.- Read the vecty file :
! --------------------------
!
      IF (present(kfninerrdta)) THEN
         flagxyo=2
         lectinfo = .TRUE.
         jnxyo=1
         CALL readxyo (kfninerrdta,vecty(:), &
     &        jnxyo,lectinfo,flagxyo) 
      ENDIF     
!     
! -3.- Compute the vectyrms file :
! -------------------------------
!
      IF (present(kfninerrdta)) THEN
         DO jy=1,jpysize
            IF (vecty(jy).NE.spvaldta) THEN
               vectyrms(jy)=SQRT(vectyrms(jy)**2 + vecty(jy)**2)
            ENDIF
         ENDDO
      ENDIF     
!
! -4.- write the dta file :
! -------------------------
!
      CALL writedta (kfnouterrdta,vectyrms(:))
!
! --- deallocation
      IF (allocated(vectyrms)) deallocate(vectyrms)
      IF (allocated(vecty)) deallocate(vecty)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algooerr','mkerrtodta')
 1001 CALL printerror2(0,1001,3,'algooerr','mkerrtodta')
!
 102   WRITE (texterror,*) 'dtanam not valid : ',kdtanam
       CALL printerror2(0,102,3,'algooerr','mkerrtodta', &
     &      comment=texterror)
 103   WRITE (texterror,*) 'dtaval not valid : ',kdtaval
       CALL printerror2(0,103,3,'algooerr','mkerrtodta', &
     &      comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkoperstddta (kfninstddta,kfninyo, &
     &     kfnoutdta,ktextoper,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose :
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
      use mod_spacexyo , only : jpoend,jpitpend,jpyend,spvaldta, &
     &     poscoefobs
      use hioxyo
      use utilmkh
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninstddta,kfninyo, &
     &     kfnoutdta,ktextoper
      CHARACTER(len=*), intent(in), optional :: kconfigo
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vectoin
      BIGREAL, dimension(:), allocatable, save :: vectyin
      BIGREAL, dimension(:), allocatable, save :: vectyinstd
      BIGREAL, dimension(:), allocatable, save :: vectyoutdta
!
      INTEGER :: allocok,jpysize,jposize,jpitpsize
      INTEGER :: flagcfg,flagxyo,jnxyo
      LOGICAL :: lectinfo,warningspvaldta
      INTEGER :: jdta,inddta,jy,ios
      CHARACTER(len=bgword) kdtanam,comment
      BIGREAL :: kdtaval,val1,val2,alpha,alpha1,alpha2,coeflimite
      INTEGER, dimension(:), allocatable :: vectynbpoint
!----------------------------------------------------------------------
!
      jpysize=jpyend
      jposize=jpoend
      jpitpsize=jpitpend
! --- allocation vectyinstd
      allocate ( vectyinstd(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectyinstd(:) = FREAL(0.0)
! --- allocation vectyin
      allocate ( vectyin(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectyin(:) = FREAL(0.0)
      IF (present(kconfigo)) THEN
! --- allocation vectoin
         allocate ( vectoin(1:jposize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectoin(:) = FREAL(0.0)
      ENDIF
! --- allocation vectyoutdta
      allocate ( vectyoutdta(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectyoutdta(:) = FREAL(0.0)
! --- allocation vectynbpoint
      allocate ( vectynbpoint(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectynbpoint(:) = 0
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&         routine mkoperstddta             &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      lmoyect=.FALSE.
      lectinfo = .FALSE.
      jnxyo=1
!
! -0.1- reading or define config.obs :
! ------------------------------------
!
      IF (present(kconfigo)) THEN
! --- allocation poscoefobs
         allocate ( poscoefobs(1:jposize,1:jpitpsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
         flagcfg=3
         CALL readcfgobs (kconfigo,flagcfg, &
     &        kposcoefobs=poscoefobs(:,:))
      ENDIF
!     
! -1.- Read the vectyinstd file :
! -------------------------------
!
      flagxyo=2
      CALL readxyo (kfninstddta,vectyinstd(:), &
     &     jnxyo,lectinfo,flagxyo) 
!     
! -2.- Read the vectyinerr file :
! -------------------------------
!
      IF (present(kconfigo)) THEN
         flagxyo=3
         CALL readxyo (kfninyo,vectoin(:), &
     &        jnxyo,lectinfo,flagxyo,kposcoefobs=poscoefobs(:,:))
         vectynbpoint(:)=0
         coeflimite=FREAL(1.0)
         CALL mkhotoy(vectoin(:),vectyin(:), &
     &           poscoefobs(:,:),spvaldta,kcoeflimite=coeflimite, &
     &           kvectynbpoint=vectynbpoint(:))

      ELSE
         flagxyo=2
         CALL readxyo (kfninyo,vectyin(:), &
     &        jnxyo,lectinfo,flagxyo)
      ENDIF 
!     
! -3.- Compute the operation :
! -------------------------------
!
      vectyoutdta(:) = spvaldta
      SELECT CASE (posit(ktextoper,'_'))
      CASE (2:bgword)        
         READ(ktextoper((posit(ktextoper,'_')+1): &
     &        lenv(ktextoper)),*,IOSTAT=ios) alpha
         IF (ios.NE.0) GOTO 101
         SELECT CASE (ktextoper(1:(posit(ktextoper,'_')-1)))
         CASE ('+')
            warningspvaldta=.FALSE.
            DO jy=1,jpysize
               IF (vectyinstd(jy).EQ.spvaldta) THEN
                  IF (vectyin(jy).EQ.spvaldta) THEN
                     vectyoutdta(jy) = spvaldta
                     warningspvaldta=.TRUE.
                  ELSE
                     vectyoutdta(jy) = ABS( vectyin(jy) )
                  ENDIF
               ELSE
                  IF (vectyin(jy).EQ.spvaldta) THEN
                     val1 = vectyinstd(jy)
                     val2 = FREAL(0.0)
                     alpha1=FREAL(1.)
                     alpha2=FREAL(0.)
                  ELSE
                     val1 = vectyinstd(jy)
                     val2 = vectyin(jy)
                     alpha1=FREAL(1.)-alpha
                     alpha2=alpha
                  ENDIF
                  vectyoutdta(jy) = SQRT( alpha1*val1*val1 + &
     &                 alpha2*val2*val2 )
               ENDIF
            ENDDO
            IF (warningspvaldta) PRINT *,'WARNING : existing spvaldta'
         CASE DEFAULT
            GOTO 101
         END SELECT
      CASE (0)        
         SELECT CASE (ktextoper(1:lenv(ktextoper)))
         CASE ('+')
            DO jy=1,jpysize
               IF (vectyinstd(jy).EQ.spvaldta) THEN
                  val1 = FREAL(0.0)
               ELSE
                  val1 = vectyinstd(jy)
               ENDIF
               IF (vectyin(jy).EQ.spvaldta) THEN
                  val2 = FREAL(0.0)
               ELSE
                  val2 = vectyin(jy)
               ENDIF
               vectyoutdta(jy) = SQRT( val1*val1 + &
     &              val2*val2 )
            ENDDO
         CASE ('-')
            DO jy=1,jpysize
               IF ((vectyinstd(jy).NE.spvaldta).AND. &
     &              (vectyin(jy).NE.spvaldta)) THEN
                  vectyoutdta(jy) = vectyinstd(jy)**2 - &
     &                 vectyin(jy)**2 
                  IF (vectyoutdta(jy).GE.FREAL(-0.01)) THEN
                     vectyoutdta(jy) = SQRT(ABS(vectyoutdta(jy)))
                  ELSE
                     vectyoutdta(jy) = spvaldta
                  ENDIF
               ENDIF
            ENDDO
         CASE DEFAULT 
            GOTO 101
         END SELECT
      CASE DEFAULT 
         GOTO 1000
      END SELECT
!
! -4.- write the dta file :
! -------------------------
!
      CALL writedta (kfnoutdta,vectyoutdta(:))
!
! --- deallocation
      IF (allocated(vectyinstd)) deallocate(vectyinstd)
      IF (allocated(vectyin)) deallocate(vectyin)
      IF (allocated(vectyoutdta)) deallocate(vectyoutdta)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algooerr','mkoperstddta')
 1001 CALL printerror2(0,1001,3,'algooerr','mkoperstddta')
!
 101  WRITE (texterror,*) 'operation not allowed : ', &
     &     ktextoper(1:lenv(ktextoper))
      CALL printerror2(0,101,3,'algooerr','algooerr',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkerrobastodta (kfninobas,kfnouterrdta, &
     &     kfninybas,kfnoutstddta,kfninybasref)
!---------------------------------------------------------------------
!
!  Purpose :
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
      use mod_spacexyo , only : &
     &     poscoefobs,jpoend,jpitpend,jprend,jpyend,spvaldta
      use hioxyo
      use hiobas
      use utilmkh
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobas,kfnouterrdta
      CHARACTER(len=*), intent(in), optional :: kfninybas,kfnoutstddta, &
     &     kfninybasref
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:), allocatable, save :: vectosmooth
      BIGREAL, dimension(:), allocatable, save :: vectof
      BIGREAL, dimension(:), allocatable, save :: vectoa
      BIGREAL, dimension(:), allocatable, save :: vectoweight
      BIGREAL, dimension(:,:), allocatable, save :: baseoref
      BIGREAL, dimension(:,:), allocatable, save :: baseyref
      BIGREAL, dimension(:), allocatable, save :: vectyerrohho
      BIGREAL, dimension(:), allocatable, save :: vectyerrof
      BIGREAL, dimension(:), allocatable, save :: vectyerroa
      BIGREAL, dimension(:), allocatable, save :: vecty
      INTEGER, dimension(:), allocatable :: vecty_int1
      INTEGER, dimension(:), allocatable :: vecty_int2
!
      INTEGER :: allocok,jpysize,jprsize,jprefsize,jposize,jpitpsize
      INTEGER :: numaction
      INTEGER, dimension(:), allocatable :: vectynbpoint
      INTEGER :: flagcfg,flagxyo,jnxyo
      INTEGER :: numjr,serie,jr,jy,jrbasdeb,jrbasfin,jref
      CHARACTER(len=bgword) :: fnameo,fnamey,vctnam
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbaseo,valbasey,valbaseyref
      BIGREAL :: coeflimite,coef,norm,valX,valY
!----------------------------------------------------------------------
!
      numaction=1
      IF (present(kfnoutstddta)) numaction=4
      IF ((present(kfninybas)).AND.(present(kfnoutstddta))) numaction=2
      IF ((present(kfninybas)).AND.(present(kfnoutstddta)) &
     &     .AND.(present(kfninybasref)))  numaction=3
      jpysize=jpyend
      jprsize=jprend
      numjr=1
      serie=0
      CALL fildirbas (fnameo,kfninobas,jprsize,numjr,serie)
      IF (numaction.EQ.3) THEN
         CALL fildirbas (fnamey,kfninybasref,jprefsize,numjr,serie)
      ENDIF
! --- allocation vecty
      allocate ( vecty(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecty(:) = FREAL(0.0)
      SELECT CASE (numaction)
      CASE (2)
! --- allocation baseyref
! --- XY,X,Y,E(XY),E(X),E(Y),E(X2),E(Y2)
         allocate ( baseyref(1:jpysize,1:8), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         baseyref(:,:) = FREAL(0.0)
! --- allocation vectyerrohho
         allocate ( vectyerrohho(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectyerrohho(:) = FREAL(0.0)
! --- allocation vecty_int1
         allocate ( vecty_int1(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vecty_int1(:) = 0
      CASE (1)
! --- allocation vectyerrohho
         allocate ( vectyerrohho(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectyerrohho(:) = FREAL(0.0)
! --- allocation vecty_int1
         allocate ( vecty_int1(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vecty_int1(:) = 0
      CASE (4)
! --- allocation vectyerrohho
         allocate ( vectyerrohho(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectyerrohho(:) = FREAL(0.0)
! --- allocation vectyerrof
         allocate ( vectyerrof(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectyerrof(:) = FREAL(0.0)
! --- allocation vecty_int1
         allocate ( vecty_int1(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vecty_int1(:) = 0
      CASE (3)
! --- allocation baseyref
         allocate ( baseyref(1:jpysize,jprefsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         baseyref(:,:) = FREAL(0.0)
! --- allocation vectyerroa
         allocate ( vectyerroa(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectyerroa(:) = FREAL(0.0)
      CASE DEFAULT
         GOTO 1000
      END SELECT
      SELECT CASE (numaction)
      CASE (1,4)
! --- Nothing
      CASE (2,3)
! --- allocation vectyerrof
         allocate ( vectyerrof(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectyerrof(:) = FREAL(0.0)
! --- allocation vecty_int2
         allocate ( vecty_int2(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vecty_int2(:) = 0
      CASE DEFAULT
         GOTO 1000
      END SELECT      
! --- allocation vectynbpoint
      allocate ( vectynbpoint(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectynbpoint(:) = 0
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&      routine mkerrobastodta              &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      lmoyect=.FALSE.
!
      SELECT CASE (numaction)
      CASE (1,2,4)
! --- Nothing
      CASE (3)
         CALL readinfobas(kfninybasref,valbaseyref)
         IF (valbaseyref.LT.0) GOTO 101
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      SELECT CASE (numaction)
      CASE (1,4)
! --- Nothing
      CASE (2,3)
         CALL readinfobas(kfninybas,valbasey)
         IF (valbasey.NE.0) GOTO 101
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      CALL readinfobas(kfninobas,valbaseo)
      IF (valbaseo.NE.0) GOTO 101
!
! -1.- load the ref.dta.bas :
! --------------------------
!
      SELECT CASE (numaction)
      CASE (1,2,4)
! --- Nothing
      CASE (3)
! --- load base.xy
         jnxyo=1
         flagxyo=2
         lectinfo=.FALSE.
         jrbasdeb=2
         jrbasfin=jprefsize
         CALL readbas(kfninybasref,baseyref(:,:),jnxyo, &
     &        jrbasdeb,jrbasfin,lectinfo,flagxyo)
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -1.- compute the obs.bas :
! --------------------------
!
      serie=1
      DO jr=1,jprsize
!
! -1.1- Read the header config obs file :
! ---------------------------------------
!
         lmodprint=(MOD(jr-1,(jprsize/5+1)).EQ.0)
         IF ((nprint.GE.1).AND.((lmodprint))) print *,  &
     &        'Jr number : ',jr,'/',jprsize
         numjr=jr
         CALL fildirbas (vctnam,kfninobas,jprend,numjr,serie)
         WRITE(fnameo,'("./",A,"/",A)') kfninobas(1:lenv(kfninobas)), &
     &        vctnam(1:lenv(vctnam)) 
         SELECT CASE (numaction)
         CASE (1,4)
! --- Nothing
         CASE (2,3)
            CALL fildirbas (vctnam,kfninybas,jprend,numjr,serie)
            WRITE(fnamey,'("./",A,"/",A)') kfninybas(1:lenv(kfninybas)), &
     &           vctnam(1:lenv(vctnam)) 
         CASE DEFAULT
            GOTO 1000
         END SELECT
         CALL evalhdrobs (fnameo,jpoend,jpitpend)
!
! -1.2- Allocate the config obs file :
! ---------------------------------------
!
         jposize=jpoend
         jpitpsize=jpitpend
! --- allocation vecto
         allocate ( vecto(1:jposize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vecto(:) = FREAL(0.0)
         SELECT CASE (numaction)
         CASE (1,2,3)
! --- allocation vectosmooth
            allocate ( vectosmooth(1:jposize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vectosmooth(:) = FREAL(0.0)
         CASE (4)
! --- Nothing
         CASE DEFAULT
            GOTO 1000
         END SELECT
         SELECT CASE (numaction)
         CASE (1,2,4)
! --- Nothing
         CASE (3)
! --- allocation baseoref
            allocate ( baseoref(1:jposize,jprefsize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            baseoref(:,:) = FREAL(0.0)
! --- allocation vectoa
            allocate ( vectoa(1:jposize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vectoa(:) = FREAL(0.0)
         CASE DEFAULT
            GOTO 1000
         END SELECT
         SELECT CASE (numaction)
         CASE (1,4)
! --- Nothing
         CASE (2,3)
! --- allocation vectof
            allocate ( vectof(1:jposize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vectof(:) = FREAL(0.0)
         CASE DEFAULT
            GOTO 1000
         END SELECT
! --- allocation poscoefobs
         allocate ( poscoefobs(1:jposize,1:jpitpsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
         IF (largweight.OR.largoestd) THEN
! --- allocation vectoweight
            allocate ( vectoweight(1:jposize), stat=allocok )
            IF (allocok.NE.0) GOTO 1001
            vectoweight(:) = FREAL(0.0)
         ENDIF
!
! -1.3- Read the config obs file :
! --------------------------------
!
         flagcfg=3
         CALL readcfgobs (fnameo,flagcfg,kposcoefobs=poscoefobs(:,:))
         flagxyo=3
         jnxyo=1
         CALL readxyo (fnameo,vecto(:), &
     &        jnxyo,lectinfo,flagxyo,poscoefobs(:,:))  
!
         SELECT CASE (numaction)
         CASE (2,3)
! --- Nothing
         CASE (1)
! --- make the config dta
            vecty(:)=spvaldta
            CALL mkhotoy(vecto(:),vecty(:), &
     &           poscoefobs(:,:),spvaldta)
!
! --- make the smooth config obs
            CALL mkhytoo(vecty(:),vectosmooth(:),poscoefobs(:,:))
!
! --- compute the vectosmooth 
            vectosmooth(:) = ( vecto(:) - vectosmooth(:) ) **2
!
! --- make the config dta
            vecty(:)=spvaldta
            vectynbpoint(:)=0
            coeflimite=FREAL(1.0)
            CALL mkhotoy(vectosmooth(:),vecty(:), &
     &           poscoefobs(:,:),spvaldta,kcoeflimite=coeflimite, &
     &           kvectynbpoint=vectynbpoint(:) )
!     
! --- Filter and compute variance 
            DO jy=1,jpysize
               IF (vectynbpoint(jy).GE.2) THEN
                  vectyerrohho(jy) = vectyerrohho(jy) +  &
     &                 vecty(jy)
                  vecty_int1(jy) = vecty_int1(jy) + 1
               ENDIF
            ENDDO
         CASE (4)
! --- make the config dta
            vecty(:)=spvaldta
            vectynbpoint(:)=0
            coeflimite=FREAL(1.0)
            CALL mkhotoy(vecto(:),vecty(:), &
     &           poscoefobs(:,:),spvaldta,kcoeflimite=coeflimite, &
     &           kvectynbpoint=vectynbpoint(:) )
!     
! --- Filter and compute mean and std
            DO jy=1,jpysize
               IF (vectynbpoint(jy).GE.1) THEN
                  vectyerrohho(jy) = vectyerrohho(jy) +  &
     &                 vecty(jy)
                  vectyerrof(jy) = vectyerrof(jy) +  &
     &                 vecty(jy)*vecty(jy)
                  vecty_int1(jy) = vecty_int1(jy) + 1
               ENDIF
            ENDDO
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
         SELECT CASE (numaction)
         CASE (1,4)
! --- Nothing
         CASE (2,3)
!
! --- load Of
            flagxyo=3
            jnxyo=1
            CALL readxyo (fnamey,vectof(:), &
     &           jnxyo,lectinfo,flagxyo,poscoefobs(:,:))  
!
            IF (numaction.EQ.2) THEN
! --- make XY in space obs
!               vectosmooth(:)=vectof(:)*vecto(:)
!               baseyref(:,1)=spvaldta
!               vectynbpoint(:)=0
!               coeflimite=FREAL(1.0)
!               CALL mkhotoy(vectosmooth(:),baseyref(:,1),
!     $              poscoefobs(:,:),spvaldta,kcoeflimite=coeflimite,
!     $              kvectynbpoint=vectynbpoint(:) )
! --- load X
               flagxyo=2
               jnxyo=1
               CALL readxyo (fnamey,baseyref(:,2), &
     &              jnxyo,lectinfo,flagxyo)  
! --- make Y
               baseyref(:,3)=spvaldta
               vectynbpoint(:)=0
               coeflimite=FREAL(1.0)
               CALL mkhotoy(vecto(:),baseyref(:,3), &
     &              poscoefobs(:,:),spvaldta,kcoeflimite=coeflimite, &
     &              kvectynbpoint=vectynbpoint(:) )
!     
               DO jy=1,jpysize
                  IF (vectynbpoint(jy).GE.1) THEN
! --- make XY in space obs
!                     baseyref(jy,1)=baseyref(jy,2)*baseyref(jy,3)
! --- 4=>E(XY)
                     baseyref(jy,4) = baseyref(jy,4) +  &
     &                    baseyref(jy,2)*baseyref(jy,3)
! --- 5=>E(X)
                     baseyref(jy,5) = baseyref(jy,5) +  &
     &                    baseyref(jy,2)
! --- 6=>E(Y)
                     baseyref(jy,6) = baseyref(jy,6) +  &
     &                    baseyref(jy,3)
! --- 7=>E(X2)
                     baseyref(jy,7) = baseyref(jy,7) +  &
     &                    baseyref(jy,2)*baseyref(jy,2)
! --- 8=>E(Y2)
                     baseyref(jy,8) = baseyref(jy,8) +  &
     &                    baseyref(jy,3)*baseyref(jy,3)
!
                     vecty_int1(jy) = vecty_int1(jy) + 1
                  ENDIF
               ENDDO
            ENDIF
!
            IF (largbias) THEN
               flagxyo=3
               jnxyo=1
               CALL readxyo (argbias,vectosmooth(:), &
     &              jnxyo,lectinfo,flagxyo,poscoefobs(:,:))  
               vecto(:) = vecto(:) + vectosmooth(:)
            ENDIF
!
! --- compute the vectof 
            IF (numaction.EQ.3) THEN
               vectoa(:) = ( vecto(:) - vectof(:) )
               vectof(:) = ( vectoa(:) ) **2
            ELSE
               vectof(:) = ( vecto(:) - vectof(:) ) **2
            ENDIF
!
! --- make the config dta
            vecty(:)=spvaldta
            vectynbpoint(:)=0
            coeflimite=FREAL(1.0)
            CALL mkhotoy(vectof(:),vecty(:), &
     &           poscoefobs(:,:),spvaldta,kcoeflimite=coeflimite, &
     &           kvectynbpoint=vectynbpoint(:) )
!     
! --- Filter and compute variance 
            DO jy=1,jpysize
               IF (vectynbpoint(jy).GE.1) THEN
                  vectyerrof(jy) = vectyerrof(jy) +  &
     &                 vecty(jy)
                  vecty_int2(jy) = vecty_int2(jy) + 1
               ENDIF
            ENDDO
!
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
         SELECT CASE (numaction)
         CASE (1,2,4)
! --- Nothing
         CASE (3)
!
            IF (largoestd) THEN
               flagxyo=3
               jnxyo=1
               CALL readxyo (argoestd,vecto(:), &
     &              jnxyo,lectinfo,flagxyo,poscoefobs(:,:))  
            ENDIF
            IF (largweight) THEN
                  flagxyo=3
                  jnxyo=1
                  CALL readxyo (argweight,vectoweight(:), &
     &                 jnxyo,lectinfo,flagxyo,poscoefobs(:,:)) 
                  IF (largoestd) vectoweight(:) = vectoweight(:)  &
     &                 / vecto(:)
            ELSE
                  IF (largoestd) vectoweight(:) = FREAL(1.0)  &
     &              / vecto(:)
            ENDIF
!
            vectof(:) = vectoa(:)
            IF (largweight.OR.largoestd) vectof(:) = &
     &           vectof(:)*vectoweight(:)
!
            DO jref=2,jprefsize
!
               IF (largweight.OR.largoestd) THEN
                  vecto(:) = vectosmooth(:)*vectoweight(:)
                  norm=DOT_PRODUCT(vecto(:),vecto(:))
               ELSE
                  norm=DOT_PRODUCT(vectosmooth(:),vectosmooth(:))
               ENDIF
               coef=DOT_PRODUCT(vectof(:),vectosmooth(:))
               coef=coef/norm
               vectoa(:) = vectoa(:) - coef*vectosmooth(:)
!
            ENDDO
!
            vectoa(:) = ( vectoa(:) ) **2
!
! --- make the config dta
            vecty(:)=spvaldta
            vectynbpoint(:)=0
            coeflimite=FREAL(1.0)
            CALL mkhotoy(vectoa(:),vecty(:), &
     &           poscoefobs(:,:),spvaldta,kcoeflimite=coeflimite, &
     &           kvectynbpoint=vectynbpoint(:) )
!     
! --- Filter and compute variance 
            DO jy=1,jpysize
               IF (vectynbpoint(jy).GE.6) THEN
                  vectyerroa(jy) = vectyerroa(jy) +  &
     &                 vecty(jy)
               ENDIF
            ENDDO
!
         CASE DEFAULT
            GOTO 1000
         END SELECT
!
! -1.10- deallocate :
! --------------------
!
         IF (allocated(vecto)) deallocate(vecto)
         IF (allocated(vectoweight)) deallocate(vectoweight)
         IF (allocated(vectof)) deallocate(vectof)
         IF (allocated(vectoa)) deallocate(vectoa)
         IF (allocated(vectosmooth)) deallocate(vectosmooth)
         IF (allocated(poscoefobs)) deallocate(poscoefobs)
!
      ENDDO
!
! -2.- Compute variance :
! -------------------------------
!
      SELECT CASE (numaction)
      CASE (3)
! --- Nothing
      CASE (1)
         DO jy=1,jpysize
            IF (vecty_int1(jy).GE.1) THEN
               vectyerrohho(jy)=SQRT(vectyerrohho(jy)/FREAL(vecty_int1(jy)-1))
            ELSE
               vectyerrohho(jy)=spvaldta
            ENDIF
         ENDDO
      CASE (2)
! --- XY,X,Y,E(XY),E(X),E(Y),E(X2),E(Y2)
! --- vectyerrohho=(E(XY)-E(X)E(Y))/SQRT(E(X2)-E(X)E(X))/SQRT(E(Y2)-E(Y)E(Y))
         DO jy=1,jpysize
            IF (vecty_int1(jy).GE.10) THEN
               baseyref(jy,5:6)=baseyref(jy,5:6)/FREAL(vecty_int1(jy))
               valX=baseyref(jy,7)-FREAL(vecty_int1(jy))*baseyref(jy,5)*baseyref(jy,5)
               valY=baseyref(jy,8)-FREAL(vecty_int1(jy))*baseyref(jy,6)*baseyref(jy,6)
               IF ((valX.GE.FREAL(0.0001)) &
     &              .AND.(valX.GE.FREAL(0.0001))) THEN
               
                  vectyerrohho(jy) = (baseyref(jy,4)- &
     &                 FREAL(vecty_int1(jy))* &
     &                 baseyref(jy,5)*baseyref(jy,6))/ &
     &                 (SQRT(valX)*SQRT(valY))
               ELSE
                  vectyerrohho(jy)=spvaldta
               ENDIF
            ELSE
               vectyerrohho(jy)=spvaldta
            ENDIF
         ENDDO
      CASE (4)
         DO jy=1,jpysize
            IF (vecty_int1(jy).GE.1) THEN
               vectyerrohho(jy)=vectyerrohho(jy)/FREAL(vecty_int1(jy))
            ELSE
               vectyerrohho(jy)=spvaldta
            ENDIF
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
      SELECT CASE (numaction)
      CASE (1)
! --- Nothing
      CASE (4)
         DO jy=1,jpysize
            IF (vecty_int1(jy).GE.1) THEN
               vectyerrof(jy)=SQRT(vectyerrof(jy)/FREAL(vecty_int1(jy))  &
     &              - vectyerrohho(jy)*vectyerrohho(jy))
            ELSE
               vectyerrof(jy)=spvaldta
            ENDIF
         ENDDO
      CASE (2,3)
         DO jy=1,jpysize
            IF (vecty_int2(jy).GE.6) THEN
               vectyerrof(jy)=SQRT(vectyerrof(jy)/FREAL(vecty_int2(jy)-1))
            ELSE
               vectyerrof(jy)=spvaldta
            ENDIF
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
      SELECT CASE (numaction)
      CASE (1,2,4)
! --- Nothing
      CASE (3)
         DO jy=1,jpysize
            IF (vecty_int2(jy).GE.6) THEN
               vectyerroa(jy)=SQRT(vectyerroa(jy)/FREAL(vecty_int2(jy)-1))
            ELSE
               vectyerroa(jy)=spvaldta
            ENDIF
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! -4.- write the dta file :
! -------------------------
!
      SELECT CASE (numaction)
      CASE (1)
         CALL writedta (kfnouterrdta,vectyerrohho(:))
      CASE (2)
         CALL writedta (kfnouterrdta,vectyerrohho(:))
         CALL writedta (kfnoutstddta,vectyerrof(:))
      CASE (3)
         CALL writedta (kfnoutstddta,vectyerrof(:))
         CALL writedta (kfnouterrdta,vectyerroa(:))
      CASE (4)
         CALL writedta (kfnouterrdta,vectyerrohho(:))
         CALL writedta (kfnoutstddta,vectyerrof(:))
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! --- deallocation
      IF (allocated(vecty)) deallocate(vecty)
      IF (allocated(vectyerrohho)) deallocate(vectyerrohho)
      IF (allocated(vectyerrof)) deallocate(vectyerrof)
      IF (allocated(vectyerroa)) deallocate(vectyerroa)
      IF (allocated(vectynbpoint)) deallocate(vectynbpoint)
      IF (allocated(vecty_int1)) deallocate(vecty_int1)
      IF (allocated(vecty_int2)) deallocate(vecty_int2)
      IF (allocated(baseyref)) deallocate(baseyref)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algooerr','mkerrobastodta')
 1001 CALL printerror2(0,1001,3,'algooerr','mkerrobastodta')
!
 101  WRITE (texterror,*) 'inobas not valid'
      CALL printerror2(0,101,3,'algooerr','mkerrobastodta', &
     &      comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkbiastodta (kfninobas,kfninybas,kfnoutbiasdta)
!---------------------------------------------------------------------
!
!  Purpose :
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
      use mod_spacexyo , only : &
     &     poscoefobs,jpoend,jpitpend,jprend,jpyend,spvaldta
      use hioxyo
      use hiobas
      use utilmkh
      use utilvalid
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobas,kfninybas,kfnoutbiasdta
!----------------------------------------------------------------------
! local declarations
! ==================
!
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:), allocatable, save :: vectof
      BIGREAL, dimension(:), allocatable, save :: vectyf
      BIGREAL, dimension(:), allocatable, save :: vecty
      INTEGER, dimension(:), allocatable :: vecty_int
!
      INTEGER :: allocok,jpysize,jprsize,jposize,jpitpsize
      INTEGER, dimension(:), allocatable :: vectynbpoint
      INTEGER :: flagcfg,flagxyo,jnxyo
      INTEGER :: numjr,serie,jr,jy
      CHARACTER(len=bgword) :: fnameo,fnamey,vctnam
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbaseo,valbasey
      BIGREAL :: coeflimite,coef,norm
!----------------------------------------------------------------------
!
      jpysize=jpyend
      jprsize=jprend
      numjr=1
      serie=0
      CALL fildirbas (fnameo,kfninobas,jprsize,numjr,serie)
! --- allocation vecty
      allocate ( vecty(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecty(:) = FREAL(0.0)
! --- allocation vectyf
      allocate ( vectyf(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectyf(:) = FREAL(0.0)
! --- allocation vecty_int
      allocate ( vecty_int(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecty_int(:) = 0
! --- allocation vectynbpoint
      allocate ( vectynbpoint(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectynbpoint(:) = 0
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&   routine mkbiastodta     &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      lmoyect=.FALSE.
!
      CALL readinfobas(kfninybas,valbasey)
      IF (valbasey.NE.0) GOTO 102
!
      CALL readinfobas(kfninobas,valbaseo)
      IF (valbaseo.NE.0) GOTO 101
!
! -1.- compute the obs.bas :
! --------------------------
!
      serie=1
      DO jr=1,jprsize
         lmodprint=(MOD(jr-1,(jprsize/5+1)).EQ.0)
         IF ((nprint.GE.1).AND.((lmodprint))) print *,  &
     &        'Jr number : ',jr,'/',jprsize
!
! -1.1- Read the header config obs file :
! ---------------------------------------
!
         numjr=jr
         CALL fildirbas (vctnam,kfninobas,jprend,numjr,serie)
         WRITE(fnameo,'("./",A,"/",A)') kfninobas(1:lenv(kfninobas)), &
     &        vctnam(1:lenv(vctnam)) 
!
         CALL fildirbas (vctnam,kfninybas,jprend,numjr,serie)
         WRITE(fnamey,'("./",A,"/",A)') kfninybas(1:lenv(kfninybas)), &
     &        vctnam(1:lenv(vctnam)) 
!
         CALL evalhdrobs (fnameo,jpoend,jpitpend)
!
! -1.2- Allocate the config obs file :
! ---------------------------------------
!
         jposize=jpoend
         jpitpsize=jpitpend
! --- allocation vecto
         allocate ( vecto(1:jposize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vecto(:) = FREAL(0.0)
! --- allocation vectof
         allocate ( vectof(1:jposize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectof(:) = FREAL(0.0)
! --- allocation poscoefobs
         allocate ( poscoefobs(1:jposize,1:jpitpsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
!
! -1.3- Read the config obs file :
! --------------------------------
!
         flagcfg=3
         CALL readcfgobs (fnameo,flagcfg,kposcoefobs=poscoefobs(:,:))
!
! -1.4- Read the obs file :
! -------------------------
!
         flagxyo=3
         jnxyo=1
         CALL readxyo (fnameo,vecto(:), &
     &        jnxyo,lectinfo,flagxyo,poscoefobs(:,:))  
!
! -1.5- Read the xf file :
! ------------------------
!
         flagxyo=3
         jnxyo=1
         CALL readxyo (fnamey,vectof(:), &
     &        jnxyo,lectinfo,flagxyo,poscoefobs(:,:))  
!
! -1.6- Computation of the bias :
! --------------------------------
!
         vecto(:) = vectof(:) - vecto(:)
!
         vecty(:)=spvaldta
         vectynbpoint(:)=0
         coeflimite=FREAL(1.0)
         CALL mkhotoy(vecto(:),vecty(:), &
     &        poscoefobs(:,:),spvaldta,kcoeflimite=coeflimite, &
     &           kvectynbpoint=vectynbpoint(:) )
!
! -1.7- Computation of the mean bias :
! ------------------------------------
!
         DO jy=1,jpysize
            IF (vectynbpoint(jy).GE.1) THEN
               vectyf(jy) = vectyf(jy) + vecty(jy)
               vecty_int(jy) = vecty_int(jy) + 1
            ENDIF
         ENDDO
!
! -1.10- deallocate :
! --------------------
!
         IF (allocated(vecto)) deallocate(vecto)
         IF (allocated(vectof)) deallocate(vectof)
         IF (allocated(poscoefobs)) deallocate(poscoefobs)
!
      ENDDO
!
! -2.- Computation of the mean bias :
! ------------------------------------
!
      DO jy=1,jpysize
         IF (vecty_int(jy).GE.1) THEN
            vectyf(jy) = vectyf(jy) / FREAL(vecty_int(jy))
         ENDIF
      ENDDO
!
! -3.- write the dta file :
! -------------------------
!
      CALL writedta (kfnoutbiasdta,vectyf(:))
!
! --- deallocation
      IF (allocated(vecty)) deallocate(vecty)
      IF (allocated(vectyf)) deallocate(vectyf)
      IF (allocated(vectynbpoint)) deallocate(vectynbpoint)
      IF (allocated(vecty_int)) deallocate(vecty_int)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algooerr','mkbiastodta')
 1001 CALL printerror2(0,1001,3,'algooerr','mkbiastodta')
!
 101   WRITE (texterror,*) 'inobas not valid'
       CALL printerror2(0,101,3,'algooerr','mkbiastodta', &
     &      comment=texterror)
 102   WRITE (texterror,*) 'inybas not valid'
       CALL printerror2(0,101,3,'algooerr','mkbiastodta', &
     &      comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algooerr
