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
! ---                    MKDTATOOBS.F90                           ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-05 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE mkdtatoobsbyobs
! --- SUBROUTINE mkdtatoobswithoutobs
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkdtatoobs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC mkdtatoobsbyobs,mkdtatoobswithoutobs

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkdtatoobsbyobs (kfninobs,kfnindta,kfnoutobs)
!---------------------------------------------------------------------
!
!  Purpose :
!  -------
!       build OBS from a DTA 
!       with the structure of filein.OBS
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
      use mod_spacexyo , only : gridijkobs,poscoefobs,jpoend,jpitpend
      use utilroa, only : mkyorms
      use hioxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      CHARACTER(len=*), intent(in) :: kfninobs,kfnindta,kfnoutobs
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      INTEGER :: allocok,jposize,jpitpsize
      INTEGER :: flagcfg,flagxyo,jnxyo
      LOGICAL :: lectinfo
!----------------------------------------------------------------------
!
      jposize=jpoend
      jpitpsize=jpitpend
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
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mkdtatoobsbyobs &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
! -1.- Read the config obs file :
! --------------------------------
!
      flagcfg=3
      CALL readcfgobs (argconfigobs,flagcfg,kposcoefobs=poscoefobs(:,:))
      flagcfg=2
      CALL readcfgobs (argconfigobs,flagcfg,kgridijkobs=gridijkobs(:))
!
! -2.- Compute/read rms values :
! ------------------------------
!
      flagxyo=3
      IF (largoestd) THEN
         lectinfo = .TRUE.
         jnxyo=1
         CALL readxyo (argoestd,vectorms(:), &
     &        jnxyo,lectinfo,flagxyo,poscoefobs(:,:))      
      ELSE
         CALL mkyorms (vectorms(:),flagxyo)
      ENDIF
!
! -3.- Read the vecty file :
! --------------------------
!
      flagxyo=3
      lectinfo = .TRUE.
      jnxyo=1
      CALL readxyo (kfnindta,vecto(:), &
     &     jnxyo,lectinfo,flagxyo,poscoefobs(:,:))      
!
! -4.- write the obs file :
! -------------------------
!
      CALL writeobs(kfnoutobs,vecto(:),vectorms(:),gridijkobs(:), &
     &     poscoefobs(:,:))
!
! --- deallocation
      IF (allocated(vecto)) deallocate(vecto)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkdtatoobs','mkdtatoobsbyobs')
 1001 CALL printerror2(0,1001,3,'mkdtatoobs','mkdtatoobsbyobs')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkdtatoobswithoutobs (kfnindta,kfnoutobs)
!---------------------------------------------------------------------
!
!  Purpose :
!  -------
!       build OBS from a DTA 
!  Method :
!  ------
!  Input :		: no
!  -----
!  Output :		: no
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : gridijkobs,poscoefobs, &
     &     jpoend,jpitpend,jpyend,spvaldta
      use hioxyo
      use utilroa, only : mkyorms
      use utilmkh
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      CHARACTER(len=*), intent(in) :: kfnindta,kfnoutobs
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:), allocatable, save :: vecty
      BIGREAL, dimension(:), allocatable, save :: vectreducedta
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      INTEGER :: allocok,jposize,jpitpsize,jpysize
      LOGICAL :: lectinfo,nospvaldta
      INTEGER :: jo,jy,flagcfg,flagxyo,jnxyo
      INTEGER :: jdta,inddta,inddtamsk,ji,jj,jk,jt,jobs,indobs,inddbs
      INTEGER :: jpifin,jpjfin,jpkfin,jptfin
      BIGREAL :: spvaldtamoyect
!----------------------------------------------------------------------
!
      jpitpsize=1
      jpysize=jpyend
      jposize=0
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         IF (obs_dim(indobs,inddbs).EQ.1) jposize = jposize +  &
     &        dta_jpi(indobs)
         IF (obs_dim(indobs,inddbs).EQ.2) jposize = jposize + &
     &        dta_jpi(indobs)*dta_jpj(indobs)
         IF (obs_dim(indobs,inddbs).EQ.3) jposize = jposize +  &
     &        dta_jpi(indobs)*dta_jpj(indobs)*dta_jpk(indobs)
         IF (obs_dim(indobs,inddbs).EQ.3) jposize = jposize + &
     &        dta_jpi(indobs)*dta_jpj(indobs)* &
     &        dta_jpk(indobs)*dta_jpt(indobs)
      ENDDO
! --- allocation vecto
      allocate ( vecto(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecto(:) = FREAL(0.0)
! --- allocation vecty
      allocate ( vecty(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecty(:) = FREAL(0.0)
      IF (largreducedta) THEN
! --- allocation vectreducedta
         allocate ( vectreducedta(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectreducedta(:) = FREAL(0.0)
      ENDIF
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
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mkdtatoobswithout &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      jpoend=jposize
      jpitpend=jpitpsize
!
! -1.- Read the vecty file :
! --------------------------
!
      flagxyo=2
      lectinfo = .TRUE.
      jnxyo=1
      CALL readxyo (kfnindta,vecty(:), &
     &     jnxyo,lectinfo,flagxyo)      
!
! -2.- Read the reducedta file :
! ------------------------------
!
      IF (largreducedta) THEN
         flagxyo=2
         lectinfo = .TRUE.
         jnxyo=1
         CALL readxyo (argreducedta,vectreducedta(:), &
     &        jnxyo,lectinfo,flagxyo)   
      ENDIF   
!
! -3.- make the gridijkobs :
! --------------------------
!
      jo=1
      DO jobs=1,obsend
         indobs = obs_ord(jobs)
         inddbs=obsnord(jobs)
         inddtamsk=0
         LOOP1 : DO jdta=1,dtaend
            IF (dta_ord(jdta).EQ.indobs) THEN
               inddtamsk=(varend-1)+jdta
               EXIT LOOP1
            ENDIF
         ENDDO LOOP1
         IF (inddtamsk.EQ.0) GOTO 1000
         jpifin=dta_jpi(indobs)
         jpjfin=dta_jpj(indobs)
         jpkfin=dta_jpk(indobs)
         jptfin=dta_jpt(indobs)
         IF (obs_dim(indobs,inddbs).LT.1) jpifin=1
         IF (obs_dim(indobs,inddbs).LT.2) jpjfin=1
         IF (obs_dim(indobs,inddbs).LT.3) jpkfin=1
         IF (obs_dim(indobs,inddbs).LT.4) jptfin=1
         IF ((.NOT.(lmoyect)).OR.((dta_ect(indobs).EQ.FREAL(1.0)) &
     &        .AND.(dta_moy(indobs).EQ.FREAL(0.0)))) THEN
            spvaldtamoyect=spvaldta
         ELSE
            spvaldtamoyect=spvaldta*dta_ect(indobs)-dta_moy(indobs)
         ENDIF
         obs_ind(indobs,inddbs)=jo
         jy=dta_ind(indobs)
         nospvaldta=.TRUE.
         IF (largreducedta) THEN
            DO jt=1,jptfin
            DO jk=1,jpkfin
            DO jj=1,jpjfin
            DO ji=1,jpifin
               IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
                  IF ((vecty(jy).NE.spvaldtamoyect) &
     &                 .AND.(vectreducedta(jy).NE.FREAL(0.0))) THEN
                     gridijkobs(jo)%longi=ji
                     gridijkobs(jo)%latj=jj
                     gridijkobs(jo)%levk=jk
                     poscoefobs(jo,:) = type_poscoef(jy,FREAL(1.0))
                     jo=jo+1
                  ELSE
                     nospvaldta=.FALSE.
                  ENDIF
                  jy=jy+1
               ENDIF
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ELSE
            DO jt=1,jptfin
            DO jk=1,jpkfin
            DO jj=1,jpjfin
            DO ji=1,jpifin
               IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
                  IF (vecty(jy).NE.spvaldtamoyect) THEN
                     gridijkobs(jo)%longi=ji
                     gridijkobs(jo)%latj=jj
                     gridijkobs(jo)%levk=jk
                     poscoefobs(jo,:) = type_poscoef(jy,FREAL(1.0))
                     jo=jo+1
                  ELSE
                     nospvaldta=.FALSE.
                  ENDIF
                  jy=jy+1
               ENDIF
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         obs_nbr(indobs,inddbs)=jo-obs_ind(indobs,inddbs)
         IF (obs_dim(indobs,inddbs).EQ.dta_dim(indobs)) THEN
            IF ((jy-dta_ind(indobs)).NE.dta_nbr(indobs)) GOTO 1000
            IF (nospvaldta) THEN
               IF (obs_nbr(indobs,inddbs).NE.dta_nbr(indobs)) GOTO 1000
            ELSE
! --- unverifiable
               print *,'existence of spvaldta for dtanam=', &
     &              dta_nam(indobs)(1:lenv(dta_nam(indobs))), &
     &              ' => pas de test de coherence'
            ENDIF
         ELSE
! --- unverifiable
            print *,'(obsdim(', &
     &           obs_nam(indobs,inddbs)(1:lenv(obs_nam(indobs,inddbs))), &
     &           ')=',obs_dim(indobs,inddbs),')<(dtadim(', &
     &            dta_nam(indobs)(1:lenv(dta_nam(indobs))), &
     &           ')=',dta_dim(indobs),') => no test of coherence'
         ENDIF
         obs_itp(indobs,inddbs)=1
      ENDDO
      jpoend=(jo-1)
      IF (nprint.GE.1) THEN
         WRITE (numout,*) ' '
            WRITE (numout,*) ' CONFIGURATION OBS :'
            WRITE (numout,*) ' -------------------'
            WRITE (numout,*) ' '
            WRITE (numout,30) 'Variables',' obs_ind ',' obs_nbr ', &
     &           ' jpitploc',' ratio(%)'
            DO jobs=1,obsend
               indobs=obs_ord(jobs)
               inddbs=obsnord(jobs)
               WRITE (numout,40) obs_nam(indobs,inddbs), &
     &              obs_ind(indobs,inddbs),obs_nbr(indobs,inddbs), &
     &              obs_itp(indobs,inddbs), &
     &              INT(FREAL(obs_nbr(indobs,inddbs)*100)/FREAL(jpoend))
            ENDDO
         WRITE (numout,*) ' '
         WRITE (numout,*) ' jpoend = ',jpoend
         WRITE (numout,*) ' jpitpend = ',jpitpend
         WRITE (numout,*) ' '
      ENDIF
!
! -3.- Compute/read rms values :
! -------------------------------
!
      flagxyo=3
      IF (largoestd) THEN
         lectinfo = .TRUE.
         jnxyo=1
         CALL readxyo (argoestd,vectorms(:jpoend), &
     &        jnxyo,lectinfo,flagxyo,poscoefobs(:jpoend,:))      
      ELSE
         CALL mkyorms (vectorms(:jpoend),flagxyo)
      ENDIF
!
! -4.- Read the vecty file :
! --------------------------
!
      CALL mkhytoo(vecty(:),vecto(:jpoend),poscoefobs(:jpoend,:))
!
! -5.- write the obs file :
! -------------------------
!
      CALL writeobs(kfnoutobs,vecto(:jpoend),vectorms(:jpoend) &
     &     ,gridijkobs(:jpoend),poscoefobs(:jpoend,:))
!
! --- deallocation
      IF (allocated(vecty)) deallocate(vecty)
      IF (allocated(vectreducedta)) deallocate(vectreducedta)
      IF (allocated(vecto)) deallocate(vecto)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
!
      RETURN
!
! --- format definitions
!
 10   FORMAT(2X,A9,3(1X,"|",2X,A9,1X))
! 20 FORMAT(5X,A3,4X,"|",1X,3(1PE11.5E2,3X))
 20   FORMAT(5X,A3,4X,"|",1X,3(I11,3X))
 30   FORMAT(2X,A9,4(1X,"|",2X,A9,1X))
 40   FORMAT(5X,A3,4X,"|",1X,4(I11,3X))
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkdtatoobs','mkdtatoobswithoutobs')
 1001 CALL printerror2(0,1001,3,'mkdtatoobs','mkdtatoobswithoutobs')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkdtatoobs
