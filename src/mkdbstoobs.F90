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
! ---                   MKDBSTOOBS.F90                            ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12 (C.E. Testut)                        ---
! --- modification : 99-02 (L. Parent)                          ---
! --- modification : 99-11 (J.M. Brankart)                      ---
! --- modification : 99-12 (C.E. Testut)                        ---
! --- modification : 10-04 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE dbstoobs
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkdbstoobs
      use mod_main
      use mknullobs
      use mkobsoper
      use mkinterp2dh
      use mkinterp3d
      IMPLICIT NONE
      PRIVATE

      PUBLIC dbstoobs

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE dbstoobs(kargindbs,kargoutobs,kargaffectobs)
!---------------------------------------------------------------------
!
!  Purpose : Extract observation from database
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
      use mod_spacexyo , only : jpdbsend
      use hioxyo
      use hiogrd
      use lioadbs
      use lioncdbs
      use liocrg
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kargindbs,kargoutobs,kargaffectobs
!----------------------------------------------------------------------
! local declarations
! ==================
      LOGICAL :: found
      INTEGER :: jextdbs,flagxy
      INTEGER :: inddta,jobs,indobs,inddbs,spos,jobsfound
      CHARACTER(len=bgword) :: fname
      INTEGER :: ltext,allocok,natgrd,jk
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&          routine dbstoobs                &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -1.- Identify observation index
! -------------------------------
!
      found=.FALSE.
      LOOP1 : DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         ltext=lenv(obs_nam(indobs,inddbs))
         IF ((kargaffectobs(1:ltext).EQ.obs_nam(indobs,inddbs)(1:ltext)) &
     &        .AND.((kargaffectobs((ltext+1):(ltext+1))).EQ.(' '))) THEN
            inddta=indobs
            found=.TRUE.
            EXIT LOOP1
         ENDIF            
      ENDDO LOOP1
      IF (.NOT.found) GOTO 101
!
! -2.- Allocate and read model horizontal grid
! --------------------------------------------
      flagxy=2
      IF ( (dta_jpi(inddta).GT.1) .OR. (dta_jpj(inddta).GT.1) ) THEN
!
! --- select horizontal grid type
        natgrd=dtangrd(inddta)
        IF (natgrd.LE.2) THEN
! --- allocation longi
           allocate ( longi(1:dta_jpi(inddta)) , stat = allocok )
           IF (allocok.GT.0) GOTO 1001
           longi(:) = FREAL(0.0)
! --- allocation latj
           allocate ( latj(1:dta_jpj(inddta)) , stat = allocok )
           IF (allocok.GT.0) GOTO 1001
           latj(:) = FREAL(0.0)
        ELSE
! --- allocation gridij
           allocate ( gridij(1:dta_jpi(inddta),1:dta_jpj(inddta)) ,  &
     &                stat = allocok )
           IF (allocok.GT.0) GOTO 1001
           gridij(:,:) = type_gridij(FREAL(0.0),FREAL(0.0))
        ENDIF
!
! --- read grid file
        CALL readgrd(flagxy,inddta)
!
      ENDIF
!
! -3.- Allocate and read vertical coordinate
! ------------------------------------------
!
      IF (dta_jpk(inddta).GT.1)  THEN
!
! --- allocation levk
         allocate ( levk(1:dta_jpk(inddta)) , stat = allocok )
         IF (allocok.GT.0) GOTO 1001
         levk = FREAL(0.0)
!
! --- read grid file
         IF (traditional) THEN
           DO jk=1,dta_jpk(inddta)
              levk(jk)=var_lev(jk,inddta)
           ENDDO
         ELSE
           CALL readlev(flagxy,inddta)
         ENDIF
!               
      ENDIF
!
! -3.- Allocate and read time coordinate
! --------------------------------------
!
      IF (dta_jpt(inddta).GT.1)  THEN
!
! --- allocation levk
         allocate ( time(1:dta_jpt(inddta)) , stat = allocok )
         IF (allocok.GT.0) GOTO 1001
         time = FREAL(0.0)
!
! --- read grid file
         CALL readtime(flagxy,inddta)
!               
      ENDIF
!
! -5.- Extract observations
! -------------------------
!
      CALL extractobs (kargindbs,kargoutobs,jobs)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkdbstoobs','dbstoobs')
 1001 CALL printerror2(0,1001,3,'mkdbstoobs','dbstoobs')
!
 101  WRITE (texterror,*) 'variable affectation error'
      CALL printerror2(0,101,3,'mkdbstoobs','dbstoobs', &
     &     comment=texterror)
 102  WRITE (texterror,*) 'you need an even frequency'
      CALL printerror2(0,102,3,'mkdbstoobs','dbstoobs', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE extractobs (kfnindbs,kfnoutobs,kjobs)
!---------------------------------------------------------------------
!
!  Purpose :
!  -------
!  Extract observations and define the observation operator
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
      use mod_spacexyo , only : jpoend,spvaldta,jpyend,jpdbsend, &
     &     gridijkobs,poscoefobs
      use hioxyo
      use lioadbs
      use lioncdbs
      use liocrg
      use utilroa, only : mkyorms
      use utilobs
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfnindbs,kfnoutobs
      INTEGER, intent(in) :: kjobs
!----------------------------------------------------------------------
! local declarations
! ==================
! vectdbs    : available Io vector
! vectdbsrms : associated observation error (obsolete)
! gridijdbs  : observation locations (2D case)
! gridijkdbs : observation locations (3D case)
! griddbs    : observation locations (4D case)
      BIGREAL, dimension(:), allocatable, save :: vectdbs
      BIGREAL, dimension(:), allocatable, save :: vectdbsrms
      TYPE (type_gridij), dimension(:), allocatable :: gridijdbs
      TYPE (type_gridijk), dimension(:), allocatable :: gridijkdbs
      TYPE (type_grid4d), dimension(:), allocatable :: griddbs
!
      BIGREAL, dimension(:), allocatable, save :: vecty
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      INTEGER :: allocok,jpysize,jpdbssize,jposize,jpitpsize
      CHARACTER(len=bgword) :: fname
      BIGREAL :: spvaldbs,scaledbs
      LOGICAL, parameter :: ldbsrms = .FALSE.
      INTEGER :: jextobs,jobs,indobs,jpoendloc,indice,spos
      INTEGER :: inddbs,jodebloc,jofinloc,jpoloc,jpitploc
      LOGICAL :: lectinfo
      INTEGER :: flagxyo,jnxyo,jo,jodeb,jofin,jitp,jogood
      CHARACTER(len=bgword) :: dbsobsnam
!
      INTEGER :: jpicrg,jpjcrg,jpkcrg,jextdbs,jicrg,jjcrg,ios,jj
      BIGREAL, dimension(:,:), allocatable :: matcrgbias
      BIGREAL, dimension(:), allocatable :: latcrgbias
      BIGREAL, dimension(:), allocatable :: loncrgbias
      BIGREAL :: spvaldbsbias,obsbias,ricrg,rjcrg
      INTEGER, dimension(:), allocatable :: pto
!----------------------------------------------------------------------
!
! -0.- Initialisation :
! ---------------------
!
      jpysize=jpyend
      jpdbssize=jpdbsend
      jposize=jpdbssize
!
      jpoend=jposize
!
      indobs=obs_ord(kjobs)
      inddbs=obsnord(kjobs)
      dbsobsnam=obs_nam(indobs,inddbs)
!
      IF (traditional) THEN
        IF ( (dta_dim(indobs).NE.2).AND.(dta_dim(indobs).NE.3) ) GOTO 1000
        jpitpsize = 2 ** dta_dim(indobs)
      ELSE
        jpitpsize=1
        IF (dta_jpi(indobs).GT.1) jpitpsize = jpitpsize * 2
        IF (dta_jpj(indobs).GT.1) jpitpsize = jpitpsize * 2
        IF (dta_jpk(indobs).GT.1) jpitpsize = jpitpsize * 2
        IF (dta_jpt(indobs).GT.1) jpitpsize = jpitpsize * 2
      ENDIF
!
! dynamic allocation
! ==================
! --- allocation vectdbs
      allocate ( vectdbs(1:jpdbssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectdbs(:) = FREAL(0.0)
! --- allocation vectdbsrms
      allocate ( vectdbsrms(1:1), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectdbsrms(:) = FREAL(0.0)
! 
      IF (traditional) THEN
        IF (dta_dim(indobs).EQ.2) THEN
! --- allocation gridijdbs
          allocate ( gridijdbs(1:jpdbssize), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
          gridijdbs(:) = type_gridij(FREAL(0.0),FREAL(0.0))
        ELSE
! --- allocation gridijkdbs
          allocate ( gridijkdbs(1:jpdbssize), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
          gridijkdbs(:) = type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
        ENDIF
      ELSE
! --- allocation griddbs
          allocate ( griddbs(1:jpdbssize), stat=allocok )
          IF (allocok.NE.0) GOTO 1001
          griddbs(:) = type_grid4d(FREAL(0.0),FREAL(0.0),FREAL(0.0),FREAL(0.0))
      ENDIF
! --- allocation vecto
      allocate ( vecto(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecto(:) = FREAL(0.0)
! --- allocation vectorms
      allocate ( vectorms(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectorms(:) = FREAL(0.0)
! --- allocation gridijkobs
      allocate ( gridijkobs(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
! --- allocation poscoefobs
      allocate ( poscoefobs(1:jposize,1:jpitpsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mkextractobs &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -1.- Read the mask for the extraction of observatoions (reducedta file)
! -----------------------------------------------------------------------
!		
      IF (largreducedta) THEN
! --- allocation vecty
         allocate ( vecty(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vecty(:) = FREAL(0.0)
         flagxyo=2
         jnxyo=1
         CALL readxyo (argreducedta,vecty(:), &
     &        jnxyo,lectinfo,flagxyo)  
      ENDIF
!
! -2.- Read the observation values and locations from the dbs file
! ----------------------------------------------------------------
!		
      jextdbs=indext(kfnindbs,extdbstab,nbextdbs)
      SELECT CASE (jextdbs)
      CASE (2)
! ====>    -indbs *.crg -outobs *.obs
         IF (dta_dim(indobs).EQ.2) THEN
            CALL readcrg (kfnindbs,vectdbs(:),dbsobsnam, &
     &                 spvaldbs,kgridij=gridijdbs(:))
         ELSE
            CALL readcrg (kfnindbs,vectdbs(:),dbsobsnam, &
     &                 spvaldbs,kgridijk=gridijkdbs(:))
         ENDIF
      CASE (3)
! ====>    -indbs *.ncdbs -outobs *.obs
         IF (traditional) THEN
           IF (dta_dim(indobs).EQ.2) THEN
              CALL readncdbs (kfnindbs,vectdbs(:),kjobs, &
     &                 spvaldbs,kgridij=gridijdbs(:))
           ELSE
              CALL readncdbs (kfnindbs,vectdbs(:),kjobs, &
     &                 spvaldbs,kgridijk=gridijkdbs(:))
           ENDIF
         ELSE
           CALL readncdbs (kfnindbs,vectdbs(:),kjobs, &
     &                spvaldbs,kgrid=griddbs(:))
         ENDIF
      CASE (4)
! ====>    -indbs *.adbs -outobs *.obs
         IF (traditional) THEN
           IF (dta_dim(indobs).EQ.2) THEN
              CALL readadbs (kfnindbs,vectdbs(:),kjobs, &
     &                  spvaldbs,kgridij=gridijdbs(:))
           ELSE
              CALL readadbs (kfnindbs,vectdbs(:),kjobs, &
     &                  spvaldbs,kgridijk=gridijkdbs(:))
           ENDIF
         ELSE
           CALL readadbs (kfnindbs,vectdbs(:),kjobs, &
     &                spvaldbs,kgrid=griddbs(:))
         ENDIF
      CASE DEFAULT
         GOTO 1000
      END SELECT 
!
! -3.- Compute the observation operator
! -------------------------------------
!
      indobs=obs_ord(kjobs)
!
      IF (traditional) THEN
        IF (dta_dim(indobs).EQ.2) THEN
          IF (largreducedta) THEN
              CALL interp2dh (vectdbs(:),vectdbsrms(:),gridijdbs(:), &
     &        vecto(:),vectorms(:),gridijkobs(:),poscoefobs(:,:), &
     &        jpoendloc,ldbsrms,spvaldbs,kjobs,kvectreducedta=vecty(:))
          ELSE
              CALL interp2dh (vectdbs(:),vectdbsrms(:),gridijdbs(:), &
     &        vecto(:),vectorms(:),gridijkobs(:),poscoefobs(:,:), &
     &        jpoendloc,ldbsrms,spvaldbs,kjobs)
          ENDIF
        ELSE
          IF (largreducedta) THEN
              CALL interp3d (vectdbs(:),vectdbsrms(:),gridijkdbs(:), &
     &        vecto(:),vectorms(:),gridijkobs(:),poscoefobs(:,:), &
     &        jpoendloc,ldbsrms,spvaldbs,kjobs,kvectreducedta=vecty(:))
          ELSE
              CALL interp3d (vectdbs(:),vectdbsrms(:),gridijkdbs(:), &
     &        vecto(:),vectorms(:),gridijkobs(:),poscoefobs(:,:), &
     &        jpoendloc,ldbsrms,spvaldbs,kjobs)
          ENDIF
        ENDIF
      ELSE
        IF (largreducedta) THEN
          CALL obsoper (vectdbs(:),griddbs(:), &
     &         vecto(:),gridijkobs(:),poscoefobs(:,:), &
     &         jpoendloc,spvaldbs,kjobs,kvectreducedta=vecty(:))
        ELSE
          CALL obsoper (vectdbs(:),griddbs(:), &
     &         vecto(:),gridijkobs(:),poscoefobs(:,:), &
     &         jpoendloc,spvaldbs,kjobs)
        ENDIF
      ENDIF
!
! -4.- Define some observation operator variables
! -----------------------------------------------
!
      indice=1
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         obs_ind(indobs,inddbs)=indice
         IF (jobs.EQ.kjobs) THEN
            obs_nbr(indobs,inddbs)=jpoendloc
            obs_itp(indobs,inddbs)=jpitpsize
         ELSE
            obs_nbr(indobs,inddbs)=0
            obs_itp(indobs,inddbs)=0
         ENDIF
         indice=indice+obs_nbr(indobs,inddbs)
      ENDDO
      jpoend=indice-1
!
! -5.- Read/compute observation RMS error
! ---------------------------------------
!
      IF ( jpoend .NE. 0 ) THEN
      flagxyo=3
      IF (largoestd) THEN
         lectinfo = .TRUE.
         jnxyo=1
         CALL readxyo (argoestd,vectorms(:jpoend), &
     &        jnxyo,lectinfo,flagxyo,poscoefobs(:jpoend,:))      
      ELSE
         CALL mkyorms (vectorms(:jpoend),flagxyo)
      ENDIF
      ENDIF
!
! -6.- Set the vertical level for 2D datasets
! -------------------------------------------
!
      IF ( jpoend .NE. 0 ) THEN
      IF (dta_dim(obs_ord(kjobs)).EQ.2) THEN
         gridijkobs(:jpoend)%levk = FREAL (1.0)
      ENDIF
      ENDIF
!
! -7.- Shift the observations if requested
! ----------------------------------------
!
      IF ( jpoend .NE. 0 ) THEN
      IF (largbiasdbs) THEN
         IF (dta_dim(obs_ord(kjobs)).NE.2) GOTO 103
         jextdbs=indext(argbiasdbs,extdbstab,nbextdbs)
         IF (jextdbs.NE.2) GOTO 102
!
         CALL evalhdrcrgsize(argbiasdbs,jpicrg,jpjcrg,jpkcrg)
!
! --- allocation loncrgbias
         allocate ( loncrgbias(1:jpicrg), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         loncrgbias(:) = FREAL(0.0)
! --- allocation latcrgbias
         allocate ( latcrgbias(1:jpjcrg), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         latcrgbias(:) = FREAL(0.0)
! --- allocation matcrgbias
         allocate ( matcrgbias(1:jpicrg,1:jpjcrg), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         matcrgbias(:,:) = FREAL(0.0)
! --- allocation pto
         allocate ( pto(1:jpoend), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         pto(:) = 0
!
         CALL readcrgbias(argbiasdbs,dbsobsnam,matcrgbias, &
     &                    loncrgbias,latcrgbias,spvaldbsbias)
!
         jo = 1
         jogood = 1
         DO WHILE (jo.LE.jpoend)
            CALL crglocate(jicrg,jjcrg,ricrg,rjcrg,gridijkobs(jo), &
     &                     loncrgbias,latcrgbias)
            IF ( (nprint.GE.3) .AND. (MOD(jo-1,5000).EQ.0) ) THEN
               print *, 'Adding bias :',jo,'/',jpoend
            ENDIF
            CALL crginterp(obsbias,matcrgbias,spvaldbsbias, &
     &                     jicrg,jjcrg,ricrg,rjcrg)
            IF (obsbias.NE.spvaldbsbias) THEN
              vecto(jo) = vecto(jo) + obsbias
              pto(jogood) = jo
              jogood = jogood + 1
            ENDIF
            jo = jo + 1
         ENDDO
!
         jpoend = jogood - 1
!
         vecto(1:jpoend) = vecto(pto(1:jpoend))
         vectorms(1:jpoend-1) = vectorms(pto(1:jpoend))
         gridijkobs(1:jpoend) = gridijkobs(pto(1:jpoend))
         poscoefobs(1:jpoend,:) = poscoefobs(pto(1:jpoend),:)
!
         IF (nprint.GE.2) THEN
            print *, 'Bias has been added; jpoend = ', jpoend
         ENDIF
!
      ENDIF
      ENDIF
!
! -8.- Scale the observations if requested
! ----------------------------------------
!
      IF ( jpoend .NE. 0 ) THEN
      IF (largscale) THEN
         READ(argscale,*,IOSTAT=ios) scaledbs
         IF (ios.NE.0) GOTO 101
         vecto(1:jpoend) = vecto(1:jpoend) * scaledbs
      ENDIF
      ENDIF
!
! -9.- Write observation file
! ---------------------------
!
      IF ( jpoend .NE. 0 ) THEN
      jodebloc=1
      jofinloc=jpoend
      jpitploc=jpitpsize
      CALL writesingleobs (kfnoutobs,vecto(jodebloc:jofinloc), &
     &     vectorms(jodebloc:jofinloc), &
     &     gridijkobs(jodebloc:jofinloc), &
     &     poscoefobs(jodebloc:jofinloc,:jpitploc), &
     &     kjobs)
      ELSE
         CALL nullobs(dbsobsnam,kfnoutobs)
      ENDIF
!
      IF (allocated(vecty)) deallocate(vecty)
      IF (allocated(pto)) deallocate(pto)
      IF (allocated(vectdbs)) deallocate(vectdbs)
      IF (allocated(vectdbsrms)) deallocate(vectdbsrms)
      IF (allocated(gridijdbs)) deallocate(gridijdbs)
      IF (allocated(vecto)) deallocate(vecto)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'mkextractobs','extractobs')
 1001 CALL printerror2(0,1001,3,'mkextractobs','extractobs')
!
 101  WRITE (texterror,*) 'argument scale not valid :', &
     &     argscale(1:lenv(argscale))
      CALL printerror2(0,101,3,'mkcrgtoobs','mkcrgtoobs', &
     &     comment=texterror)
 102  WRITE (texterror,*) 'argument biasdbs not valid :', &
     &     argbiasdbs(1:lenv(argbiasdbs))
      CALL printerror2(0,102,3,'mkcrgtoobs','mkcrgtoobs', &
     &     comment=texterror)
 103  WRITE (texterror,*) 'a bias can only be removed for 2D dbs :', &
     &     argbiasdbs(1:lenv(argbiasdbs))
      CALL printerror2(0,103,3,'mkcrgtoobs','mkcrgtoobs', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkdbstoobs
