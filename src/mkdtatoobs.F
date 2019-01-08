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
C ---                    MKDTATOOBS.F                           ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-05 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                        ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE mkdtatoobsbyobs
C --- SUBROUTINE mkdtatoobswithoutobs
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkdtatoobs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC mkdtatoobsbyobs,mkdtatoobswithoutobs

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkdtatoobsbyobs (kfninobs,kfnindta,kfnoutobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
CCC  -------
CCC       build OBS from a DTA 
CCC       with the structure of filein.OBS
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
      use mod_spacexyo , only : gridijkobs,poscoefobs,jpoend,jpitpend
      use utilroa, only : mkyorms
      use hioxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ==================
      CHARACTER(len=*), intent(in) :: kfninobs,kfnindta,kfnoutobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:), allocatable, save :: vectorms
C
      INTEGER :: allocok,jposize,jpitpsize
      INTEGER :: flagcfg,flagxyo,jnxyo
      LOGICAL :: lectinfo
CC----------------------------------------------------------------------
C
      jposize=jpoend
      jpitpsize=jpitpend
C --- allocation vecto
      allocate ( vecto(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecto(:) = FREAL(0.0)
C --- allocation vectorms
      allocate ( vectorms(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectorms(:) = FREAL(0.0)
C --- allocation gridijkobs
      allocate ( gridijkobs(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
C --- allocation poscoefobs
      allocate ( poscoefobs(1:jposize,1:jpitpsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
CCC---------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mkdtatoobsbyobs &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
C
C -0.- Initialisation :
C ---------------------
C
C -1.- Read the config obs file :
C --------------------------------
C
      flagcfg=3
      CALL readcfgobs (argconfigobs,flagcfg,kposcoefobs=poscoefobs(:,:))
      flagcfg=2
      CALL readcfgobs (argconfigobs,flagcfg,kgridijkobs=gridijkobs(:))
C
C -2.- Compute/read rms values :
C ------------------------------
C
      flagxyo=3
      IF (largoestd) THEN
         lectinfo = .TRUE.
         jnxyo=1
         CALL readxyo (argoestd,vectorms(:),
     $        jnxyo,lectinfo,flagxyo,poscoefobs(:,:))      
      ELSE
         CALL mkyorms (vectorms(:),flagxyo)
      ENDIF
C
C -3.- Read the vecty file :
C --------------------------
C
      flagxyo=3
      lectinfo = .TRUE.
      jnxyo=1
      CALL readxyo (kfnindta,vecto(:),
     $     jnxyo,lectinfo,flagxyo,poscoefobs(:,:))      
C
C -4.- write the obs file :
C -------------------------
C
      CALL writeobs(kfnoutobs,vecto(:),vectorms(:),gridijkobs(:),
     $     poscoefobs(:,:))
C
C --- deallocation
      IF (allocated(vecto)) deallocate(vecto)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mkdtatoobs','mkdtatoobsbyobs')
 1001 CALL printerror2(0,1001,3,'mkdtatoobs','mkdtatoobsbyobs')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkdtatoobswithoutobs (kfnindta,kfnoutobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
CCC  -------
CCC       build OBS from a DTA 
CCC  Method :
CCC  ------
CCC  Input :		: no
CCC  -----
CCC  Output :		: no
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      use mod_spacexyo , only : gridijkobs,poscoefobs,
     $     jpoend,jpitpend,jpyend,spvaldta
      use hioxyo
      use utilroa, only : mkyorms
      use utilmkh
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ==================
      CHARACTER(len=*), intent(in) :: kfnindta,kfnoutobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:), allocatable, save :: vecty
      BIGREAL, dimension(:), allocatable, save :: vectreducedta
      BIGREAL, dimension(:), allocatable, save :: vectorms
C
      INTEGER :: allocok,jposize,jpitpsize,jpysize
      LOGICAL :: lectinfo,nospvaldta
      INTEGER :: jo,jy,flagcfg,flagxyo,jnxyo
      INTEGER :: jdta,inddta,inddtamsk,ji,jj,jk,jt,jobs,indobs,inddbs
      INTEGER :: jpifin,jpjfin,jpkfin,jptfin
      BIGREAL :: spvaldtamoyect
CC----------------------------------------------------------------------
C
      jpitpsize=1
      jpysize=jpyend
      jposize=0
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         IF (obs_dim(indobs,inddbs).EQ.1) jposize = jposize + 
     $        dta_jpi(indobs)
         IF (obs_dim(indobs,inddbs).EQ.2) jposize = jposize +
     $        dta_jpi(indobs)*dta_jpj(indobs)
         IF (obs_dim(indobs,inddbs).EQ.3) jposize = jposize + 
     $        dta_jpi(indobs)*dta_jpj(indobs)*dta_jpk(indobs)
         IF (obs_dim(indobs,inddbs).EQ.3) jposize = jposize +
     $        dta_jpi(indobs)*dta_jpj(indobs)*
     $        dta_jpk(indobs)*dta_jpt(indobs)
      ENDDO
C --- allocation vecto
      allocate ( vecto(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecto(:) = FREAL(0.0)
C --- allocation vecty
      allocate ( vecty(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecty(:) = FREAL(0.0)
      IF (largreducedta) THEN
C --- allocation vectreducedta
         allocate ( vectreducedta(1:jpysize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         vectreducedta(:) = FREAL(0.0)
      ENDIF
C --- allocation vectorms
      allocate ( vectorms(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectorms(:) = FREAL(0.0)
C --- allocation gridijkobs
      allocate ( gridijkobs(1:jposize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
C --- allocation poscoefobs
      allocate ( poscoefobs(1:jposize,1:jpitpsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
CCC---------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mkdtatoobswithout &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
C
C -0.- Initialisation :
C ---------------------
C
      jpoend=jposize
      jpitpend=jpitpsize
C
C -1.- Read the vecty file :
C --------------------------
C
      flagxyo=2
      lectinfo = .TRUE.
      jnxyo=1
      CALL readxyo (kfnindta,vecty(:),
     $     jnxyo,lectinfo,flagxyo)      
C
C -2.- Read the reducedta file :
C ------------------------------
C
      IF (largreducedta) THEN
         flagxyo=2
         lectinfo = .TRUE.
         jnxyo=1
         CALL readxyo (argreducedta,vectreducedta(:),
     $        jnxyo,lectinfo,flagxyo)   
      ENDIF   
C
C -3.- make the gridijkobs :
C --------------------------
C
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
         IF ((.NOT.(lmoyect)).OR.((dta_ect(indobs).EQ.FREAL(1.0))
     $        .AND.(dta_moy(indobs).EQ.FREAL(0.0)))) THEN
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
                  IF ((vecty(jy).NE.spvaldtamoyect)
     $                 .AND.(vectreducedta(jy).NE.FREAL(0.0))) THEN
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
C --- unverifiable
               print *,'existence of spvaldta for dtanam=',
     $              dta_nam(indobs)(1:lenv(dta_nam(indobs))),
     $              ' => pas de test de coherence'
            ENDIF
         ELSE
C --- unverifiable
            print *,'(obsdim(',
     $           obs_nam(indobs,inddbs)(1:lenv(obs_nam(indobs,inddbs))),
     $           ')=',obs_dim(indobs,inddbs),')<(dtadim(',
     $            dta_nam(indobs)(1:lenv(dta_nam(indobs))),
     $           ')=',dta_dim(indobs),') => no test of coherence'
         ENDIF
         obs_itp(indobs,inddbs)=1
      ENDDO
      jpoend=(jo-1)
      IF (nprint.GE.1) THEN
         WRITE (numout,*) ' '
            WRITE (numout,*) ' CONFIGURATION OBS :'
            WRITE (numout,*) ' -------------------'
            WRITE (numout,*) ' '
            WRITE (numout,30) 'Variables',' obs_ind ',' obs_nbr ',
     $           ' jpitploc',' ratio(%)'
            DO jobs=1,obsend
               indobs=obs_ord(jobs)
               inddbs=obsnord(jobs)
               WRITE (numout,40) obs_nam(indobs,inddbs),
     $              obs_ind(indobs,inddbs),obs_nbr(indobs,inddbs),
     $              obs_itp(indobs,inddbs),
     $              INT(FREAL(obs_nbr(indobs,inddbs)*100)/FREAL(jpoend))
            ENDDO
         WRITE (numout,*) ' '
         WRITE (numout,*) ' jpoend = ',jpoend
         WRITE (numout,*) ' jpitpend = ',jpitpend
         WRITE (numout,*) ' '
      ENDIF
C
C -3.- Compute/read rms values :
C -------------------------------
C
      flagxyo=3
      IF (largoestd) THEN
         lectinfo = .TRUE.
         jnxyo=1
         CALL readxyo (argoestd,vectorms(:jpoend),
     $        jnxyo,lectinfo,flagxyo,poscoefobs(:jpoend,:))      
      ELSE
         CALL mkyorms (vectorms(:jpoend),flagxyo)
      ENDIF
C
C -4.- Read the vecty file :
C --------------------------
C
      CALL mkhytoo(vecty(:),vecto(:jpoend),poscoefobs(:jpoend,:))
C
C -5.- write the obs file :
C -------------------------
C
      CALL writeobs(kfnoutobs,vecto(:jpoend),vectorms(:jpoend)
     $     ,gridijkobs(:jpoend),poscoefobs(:jpoend,:))
C
C --- deallocation
      IF (allocated(vecty)) deallocate(vecty)
      IF (allocated(vectreducedta)) deallocate(vectreducedta)
      IF (allocated(vecto)) deallocate(vecto)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
C
      RETURN
C
C --- format definitions
C
 10   FORMAT(2X,A9,3(1X,"|",2X,A9,1X))
C 20 FORMAT(5X,A3,4X,"|",1X,3(1PE11.5E2,3X))
 20   FORMAT(5X,A3,4X,"|",1X,3(I11,3X))
 30   FORMAT(2X,A9,4(1X,"|",2X,A9,1X))
 40   FORMAT(5X,A3,4X,"|",1X,4(I11,3X))
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mkdtatoobs','mkdtatoobswithoutobs')
 1001 CALL printerror2(0,1001,3,'mkdtatoobs','mkdtatoobswithoutobs')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkdtatoobs
