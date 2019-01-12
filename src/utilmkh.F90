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
! ---                    UTILMK.F90                               ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 98-06 (C.E. Testut)                        ---
! --- revised      : 99-05 (C.E. Testut)                        ---
! --- revised      : 00-03 (J.M. Brankart)                      ---
! --- revised      : 01-06 (C.E. Testut)                        ---
! --- revised      : 03-04 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE mkhytoo      : Extract Vo vector from Vy vector
! --- SUBROUTINE mkhytoocompt : Extract Vo vector from Vy vector:
! ---                           set Vo to the product of all Vy
! ---                           connected (by Hyo) values
! --- SUBROUTINE mkhytou      : Extract Vu vector (subset of Vo vector
! ---                           inside one local data section)
! ---                           from Vy vector
! --- SUBROUTINE mkhotoy      : Reconstruct (roughly) Vy vector
! ---                           from Vo vector. This is a (very)
! ---                           crude inverse of Hyo operator
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilmkh
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC mkhytoo,mkhytoocompt,mkhytou,mkhotoy

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkhytoo(kvectyin,kvectoout,kposcoefobs)
!---------------------------------------------------------------------
!
!  Purpose : Extract Vo vector from Vy vector
!  -------
!  Method :  Use Hyo operator
!  ------
!  Input :   kvectyin : Vy vector
!  -----
!  Output :  kvectoout : Vo vector
!  ------    kposcoefobs : Hyo operator
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(in) :: kvectyin
      BIGREAL, dimension(:), intent(out) :: kvectoout
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jposize,jpitpsize
      INTEGER :: jy,jo,jodeb,jofin,jitp,jpitpfin
      INTEGER :: jobs,indobs,inddbs
      LOGICAL :: ltransfert
      BIGREAL :: sxy_moy,sxy_ect
!----------------------------------------------------------------------
!
! Check size of input arrays
      jposize = size(kvectoout,1)
      jpitpsize = size(kposcoefobs,2)
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 102
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         IF (jpitpsize.LT.obs_itp(indobs,inddbs)) GOTO 102
      ENDDO
!
! Check if Vy and Vo vectors have been centered and reduced
! using the same statistics. Decide to make correction
! if necessary.
      ltransfert=.FALSE.
      IF (lmoyect) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltransfert=(ltransfert.OR. &
     &           (obs_ect(indobs,inddbs).NE.dta_ect(indobs)))
            ltransfert=(ltransfert.OR. &
     &           (obs_moy(indobs,inddbs).NE.dta_moy(indobs)))
         ENDDO
      ENDIF
!
      kvectoout(:) = FREAL(0.0)
!
! Extract Vo vector from Vy vector using Hyo operator
! Make centering/reduction correction if necessary
      IF (ltransfert) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            sxy_moy=(obs_moy(indobs,inddbs)-dta_moy(indobs)) &
     &           /obs_ect(indobs,inddbs)
            sxy_ect=obs_ect(indobs,inddbs)/dta_ect(indobs)
            jodeb=obs_ind(indobs,inddbs)
            jofin=jodeb-1+obs_nbr(indobs,inddbs)
            jpitpfin=obs_itp(indobs,inddbs)
            DO jitp=1,jpitpfin
            DO jo=jodeb,jofin
               kvectoout(jo) = kvectoout(jo) + &
     &             (kvectyin(kposcoefobs(jo,jitp)%pos)/sxy_ect-sxy_moy) &
     &              *kposcoefobs(jo,jitp)%coef
            ENDDO
            ENDDO
         ENDDO
      ELSE
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            jodeb=obs_ind(indobs,inddbs)
            jofin=jodeb-1+obs_nbr(indobs,inddbs)
            jpitpfin=obs_itp(indobs,inddbs)
            DO jitp=1,jpitpfin
            DO jo=jodeb,jofin
               kvectoout(jo) = kvectoout(jo) + &
     &              kvectyin(kposcoefobs(jo,jitp)%pos) &
     &              *kposcoefobs(jo,jitp)%coef
            ENDDO
            ENDDO
         ENDDO
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilmkh','mkhytoo')
 1001 CALL printerror2(0,1001,3,'utilmkh','mkhytoo')
!
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,3,'utilmkh','mkhytoo',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkhytoocompt(kvectyin,kvectoout,kposcoefobs)
!---------------------------------------------------------------------
!
!  Purpose : Extract Vo vector from Vy vector (Vy contains only 0 and 1):
!  -------   set Vo to the product of all Vy connected (by Hyo) values
!
!  Method :  Use Hyo operator
!  ------
!  Input :   kvectyin : Vy vector
!  -----
!  Output :  kvectoout : Vo vector
!  ------    kposcoefobs : Hyo operator
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(in) :: kvectyin
      BIGREAL, dimension(:), intent(out) :: kvectoout
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jposize,jpitpsize
      INTEGER :: jy,jo,jodeb,jofin,jitp,jpitpfin
      INTEGER :: jobs,indobs,inddbs
      LOGICAL :: ltransfert
      BIGREAL :: sxy_moy,sxy_ect
!----------------------------------------------------------------------
!
! Compute size of input arrays
      jposize = size(kvectoout,1)
      jpitpsize = size(kposcoefobs,2)
!
! Extract Vo vector from Vy vector using Hyo operator
! Vo value is equal to the product of all Vy connected values
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         jodeb=obs_ind(indobs,inddbs)
         jofin=jodeb-1+obs_nbr(indobs,inddbs)
         jpitpfin=obs_itp(indobs,inddbs)
         kvectoout(jodeb:jofin) =  &
     &        kvectyin(kposcoefobs(jodeb:jofin,1)%pos)
         DO jitp=2,jpitpfin
            kvectoout(jodeb:jofin) =  &
     &           kvectoout(jodeb:jofin) * &
     &           kvectyin(kposcoefobs(jodeb:jofin,jitp)%pos)
            
         ENDDO
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilmkh','mkhytoocompt')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkhytou(kvectyin,kvectoout,kposcoefobs)
!---------------------------------------------------------------------
!
!  Purpose : Extract Vu vector (subset of Vo vector inside one
!  -------   local data section) from Vy vector
!
!  Method :  Use Hyu operator (subset of Hyo operator)
!  ------
!  Input :   kvectyin : Vy vector
!  -----
!  Output :  kvectoout : Vu vector (subset of Vo vector)
!  ------    kposcoefobs : Hyu operator (subset of Hyo operator)
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(in) :: kvectyin
      BIGREAL, dimension(:), intent(out) :: kvectoout
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jposize,jpitpsize
      INTEGER :: jy,jo,jodeb,jofin,jitp,jpitpfin
      INTEGER :: jobs,indobs,inddbs
      LOGICAL :: ltransfert
      BIGREAL :: sxy_moy,sxy_ect
!----------------------------------------------------------------------
!
! Check size of input arrays
      jposize = size(kvectoout,1)
      jpitpsize = size(kposcoefobs,2)
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 1000
!
! Check if Vy and Vo vectors have been centered and reduced
! using the same statistics. Decide to make correction
! if necessary.
      ltransfert=.FALSE.
      IF (lmoyect) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltransfert=(ltransfert.OR. &
     &           (obs_ect(indobs,inddbs).NE.dta_ect(indobs)))
            ltransfert=(ltransfert.OR. &
     &           (obs_moy(indobs,inddbs).NE.dta_moy(indobs)))
         ENDDO
      ENDIF
!
! Extract Vu vector from Vy vector using Hyu operator
! Centering/reduction correction is impossible in such case
      IF (.NOT.ltransfert) THEN
         DO jo=1,jposize
            kvectoout(jo) = kvectyin(MAX(1,kposcoefobs(jo,1)%pos)) &
     &           *kposcoefobs(jo,1)%coef
         ENDDO
         DO jitp=2,jpitpsize
            DO jo=1,jposize
               kvectoout(jo) = kvectoout(jo) + &
     &              kvectyin(MAX(1,kposcoefobs(jo,jitp)%pos)) &
     &              *kposcoefobs(jo,jitp)%coef
            ENDDO
         ENDDO
      ELSE
         GOTO 102
      ENDIF
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilmkh','mkhytou')
!
 102  WRITE (texterror,*) 'Cannot extract subset of Vo vector from ', &
     &     'Vy vector: centering/reduction statistics must be the same'
      CALL printerror2(0,102,3,'utilmkh','mkhytou',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkhotoy(kvectoin,kvectyout,kposcoefobs,kspvaly, &
     &     kcoeflimite,kvectynbpoint)
!---------------------------------------------------------------------
!
!  Purpose : Reconstruct (roughly) Vy vector from Vo vector
!  -------   This is a (very) crude inverse of Hyo operator
!
!  Method :  Fill output Vy vector with special values
!  ------    Count the number of observations available
!            to estimate each Vy value
!            Only observations such that interpolation coefficient
!            (in Hyo) is larger than specified value (limit coefficient)
!            are taken into account
!            Estimate Vy vector from Vo vector, by making the
!            average of available observations
!
!  Input :   kvectoin : Vo vector
!  -----     kposcoefobs : Hyo operator
!            kspvaly : special value
!            kcoeflimite : limit coefficient (inside ]0,1])
!           
!  Output :  kvectyout : Vy vector
!  ------    kvectynbpoint : number of observations used
!                            to compute each Vy value
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpyend
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(in) :: kvectoin
      BIGREAL, dimension(:), intent(out) :: kvectyout
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
      BIGREAL, intent(in) :: kspvaly
      BIGREAL, intent(in), optional :: kcoeflimite
      INTEGER, dimension(:), intent(out), optional :: kvectynbpoint      
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable :: vectycompteur
      INTEGER :: allocok,jpysize,jposize,jpitpsize
      INTEGER :: jy,jydeb,jyfin,jo,jodeb,jofin,jitp,jpitpfin
      INTEGER :: jobs,indobs,inddbs,jdta,inddta
      LOGICAL :: ltransfert
      BIGREAL :: sxy_moy,sxy_ect,epsilon,coeflimite
!----------------------------------------------------------------------
!
! Set default limit coefficient value
      epsilon=FREAL(0.001)
      IF (present(kcoeflimite)) THEN
         coeflimite=kcoeflimite
      ELSE
         coeflimite=FREAL(0.75)
      ENDIF
!
! Check size of input arrays
      jposize = size(kvectoin,1)
      jpitpsize = size(kposcoefobs,2)
      jpysize=jpyend
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 102
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         IF (jpitpsize.LT.obs_itp(indobs,inddbs)) GOTO 102
      ENDDO
!
! --- allocation vectycompteur
      allocate ( vectycompteur(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectycompteur(:) = FREAL(0.0)
!
! Check if Vy and Vo vectors have been centered and reduced
! using the same statistics. Decide to make correction
! if necessary.
      ltransfert=.FALSE.
      IF (lmoyect) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltransfert=(ltransfert.OR. &
     &           (obs_ect(indobs,inddbs).NE.dta_ect(indobs)))
            ltransfert=(ltransfert.OR. &
     &           (obs_moy(indobs,inddbs).NE.dta_moy(indobs)))
         ENDDO
      ENDIF
!
! Fill output Vy vector with special values
! Center/reduce special value if necessary (!)
      IF (lmoyect) THEN
         DO jdta=1,dtaend
            inddta=dta_ord(jdta)
            sxy_moy=dta_moy(inddta)/dta_ect(inddta)
            sxy_ect=dta_ect(inddta)
            jydeb=dta_ind(inddta)
            jyfin=dta_ind(inddta)+dta_nbr(inddta)-1
            IF ((sxy_ect.EQ.FREAL(1.0)).AND.(sxy_moy.EQ.FREAL(0.0))) THEN
               kvectyout(jydeb:jyfin) = kspvaly
            ELSE
               kvectyout(jydeb:jyfin) = kspvaly/sxy_ect-sxy_moy
            ENDIF
         ENDDO
      ELSE
         kvectyout(:) = kspvaly
      ENDIF
!
! Count the number of observations available
! to estimate each Vy value
! Only observations such that interpolation coefficient
! (in Hyo) is larger than specified value (limit coefficient)
! are taken into account
! Where there is at least one observation available,
! initialize output Vy vector to zero
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         jodeb=obs_ind(indobs,inddbs)
         jofin=jodeb-1+obs_nbr(indobs,inddbs)
         jpitpfin=obs_itp(indobs,inddbs)
         DO jitp=1,jpitpfin
         DO jo=jodeb,jofin
            vectycompteur(kposcoefobs(jo,jitp)%pos) =  &
     &        (FREAL(MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite)))+ &
     &        FREAL(1-MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite))) &
     &           *epsilon) + &
     &           vectycompteur(kposcoefobs(jo,jitp)%pos)
         ENDDO
         ENDDO
         DO jitp=1,jpitpfin
         DO jo=jodeb,jofin
            IF (vectycompteur(kposcoefobs(jo,jitp)%pos).NE.FREAL(0.0))  &
     &           kvectyout(kposcoefobs(jo,jitp)%pos)=FREAL(0.0)
         ENDDO
         ENDDO
      ENDDO
!
! Estimate Vy vector from Vo vector, by making the
! average of available observations
! Make centering/reduction correction if necessary
      IF (ltransfert) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            sxy_moy=(dta_moy(indobs)-obs_moy(indobs,inddbs)) &
     &           /dta_ect(indobs)
            sxy_ect=dta_ect(indobs)/obs_ect(indobs,inddbs)
            jodeb=obs_ind(indobs,inddbs)
            jofin=jodeb-1+obs_nbr(indobs,inddbs)
            jpitpfin=obs_itp(indobs,inddbs)
            DO jitp=1,jpitpfin
            DO jo=jodeb,jofin
               kvectyout(kposcoefobs(jo,jitp)%pos)=  &
     &              kvectyout(kposcoefobs(jo,jitp)%pos) + &
     &              ((( kvectoin(jo) * &
     &         (FREAL(MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite)))+ &
     &         FREAL(1-MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite))) &
     &              *epsilon) ) &
     &              / sxy_ect-sxy_moy) /  &
     &              vectycompteur(kposcoefobs(jo,jitp)%pos) )
            ENDDO
            ENDDO
         ENDDO
      ELSE
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            jodeb=obs_ind(indobs,inddbs)
            jofin=jodeb-1+obs_nbr(indobs,inddbs)
            jpitpfin=obs_itp(indobs,inddbs)
            DO jitp=1,jpitpfin
            DO jo=jodeb,jofin
               kvectyout(kposcoefobs(jo,jitp)%pos)=  &
     &              kvectyout(kposcoefobs(jo,jitp)%pos) + &
     &              ( kvectoin(jo) * &
     &         (FREAL(MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite)))+ &
     &         FREAL(1-MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite))) &
     &              *epsilon) )
            ENDDO
            ENDDO
         ENDDO
         DO jy=1,jpyend
            IF (vectycompteur(jy).GT.FREAL(0.0)) THEN
               kvectyout(jy) = kvectyout(jy) / vectycompteur(jy)
            ENDIF
         ENDDO
      ENDIF
      IF (present(kvectynbpoint)) THEN
         kvectynbpoint(:)=INT(vectycompteur(:))
      ENDIF
!
      IF (allocated(vectycompteur)) deallocate (vectycompteur)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilmkh','mkhotoy')
 1001 CALL printerror2(0,1001,3,'utilmkh','mkhotoy')
!
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,1,'utilmkh','mkhotoy',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilmkh
