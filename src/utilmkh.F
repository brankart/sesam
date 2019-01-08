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
C ---                    UTILMK.F                               ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 98-06 (C.E. Testut)                        ---
C --- revised      : 99-05 (C.E. Testut)                        ---
C --- revised      : 00-03 (J.M. Brankart)                      ---
C --- revised      : 01-06 (C.E. Testut)                        ---
C --- revised      : 03-04 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE mkhytoo      : Extract Vo vector from Vy vector
C --- SUBROUTINE mkhytoocompt : Extract Vo vector from Vy vector:
C ---                           set Vo to the product of all Vy
C ---                           connected (by Hyo) values
C --- SUBROUTINE mkhytou      : Extract Vu vector (subset of Vo vector
C ---                           inside one local data section)
C ---                           from Vy vector
C --- SUBROUTINE mkhotoy      : Reconstruct (roughly) Vy vector
C ---                           from Vo vector. This is a (very)
C ---                           crude inverse of Hyo operator
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilmkh
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC mkhytoo,mkhytoocompt,mkhytou,mkhotoy

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkhytoo(kvectyin,kvectoout,kposcoefobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Extract Vo vector from Vy vector
CCC  -------
CCC  Method :  Use Hyo operator
CCC  ------
CCC  Input :   kvectyin : Vy vector
CCC  -----
CCC  Output :  kvectoout : Vo vector
CCC  ------    kposcoefobs : Hyo operator
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(in) :: kvectyin
      BIGREAL, dimension(:), intent(out) :: kvectoout
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jposize,jpitpsize
      INTEGER :: jy,jo,jodeb,jofin,jitp,jpitpfin
      INTEGER :: jobs,indobs,inddbs
      LOGICAL :: ltransfert
      BIGREAL :: sxy_moy,sxy_ect
CC----------------------------------------------------------------------
C
C Check size of input arrays
      jposize = size(kvectoout,1)
      jpitpsize = size(kposcoefobs,2)
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 102
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         IF (jpitpsize.LT.obs_itp(indobs,inddbs)) GOTO 102
      ENDDO
C
C Check if Vy and Vo vectors have been centered and reduced
C using the same statistics. Decide to make correction
C if necessary.
      ltransfert=.FALSE.
      IF (lmoyect) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltransfert=(ltransfert.OR.
     $           (obs_ect(indobs,inddbs).NE.dta_ect(indobs)))
            ltransfert=(ltransfert.OR.
     $           (obs_moy(indobs,inddbs).NE.dta_moy(indobs)))
         ENDDO
      ENDIF
C
      kvectoout(:) = FREAL(0.0)
C
C Extract Vo vector from Vy vector using Hyo operator
C Make centering/reduction correction if necessary
      IF (ltransfert) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            sxy_moy=(obs_moy(indobs,inddbs)-dta_moy(indobs))
     $           /obs_ect(indobs,inddbs)
            sxy_ect=obs_ect(indobs,inddbs)/dta_ect(indobs)
            jodeb=obs_ind(indobs,inddbs)
            jofin=jodeb-1+obs_nbr(indobs,inddbs)
            jpitpfin=obs_itp(indobs,inddbs)
            DO jitp=1,jpitpfin
            DO jo=jodeb,jofin
               kvectoout(jo) = kvectoout(jo) +
     $             (kvectyin(kposcoefobs(jo,jitp)%pos)/sxy_ect-sxy_moy)
     $              *kposcoefobs(jo,jitp)%coef
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
               kvectoout(jo) = kvectoout(jo) +
     $              kvectyin(kposcoefobs(jo,jitp)%pos)
     $              *kposcoefobs(jo,jitp)%coef
            ENDDO
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilmkh','mkhytoo')
 1001 CALL printerror2(0,1001,3,'utilmkh','mkhytoo')
C
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,3,'utilmkh','mkhytoo',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkhytoocompt(kvectyin,kvectoout,kposcoefobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Extract Vo vector from Vy vector (Vy contains only 0 and 1):
CCC  -------   set Vo to the product of all Vy connected (by Hyo) values
CCC
CCC  Method :  Use Hyo operator
CCC  ------
CCC  Input :   kvectyin : Vy vector
CCC  -----
CCC  Output :  kvectoout : Vo vector
CCC  ------    kposcoefobs : Hyo operator
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(in) :: kvectyin
      BIGREAL, dimension(:), intent(out) :: kvectoout
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jposize,jpitpsize
      INTEGER :: jy,jo,jodeb,jofin,jitp,jpitpfin
      INTEGER :: jobs,indobs,inddbs
      LOGICAL :: ltransfert
      BIGREAL :: sxy_moy,sxy_ect
CC----------------------------------------------------------------------
C
C Compute size of input arrays
      jposize = size(kvectoout,1)
      jpitpsize = size(kposcoefobs,2)
C
C Extract Vo vector from Vy vector using Hyo operator
C Vo value is equal to the product of all Vy connected values
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         jodeb=obs_ind(indobs,inddbs)
         jofin=jodeb-1+obs_nbr(indobs,inddbs)
         jpitpfin=obs_itp(indobs,inddbs)
         kvectoout(jodeb:jofin) = 
     $        kvectyin(kposcoefobs(jodeb:jofin,1)%pos)
         DO jitp=2,jpitpfin
            kvectoout(jodeb:jofin) = 
     $           kvectoout(jodeb:jofin) *
     $           kvectyin(kposcoefobs(jodeb:jofin,jitp)%pos)
            
         ENDDO
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilmkh','mkhytoocompt')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkhytou(kvectyin,kvectoout,kposcoefobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Extract Vu vector (subset of Vo vector inside one
CCC  -------   local data section) from Vy vector
CCC
CCC  Method :  Use Hyu operator (subset of Hyo operator)
CCC  ------
CCC  Input :   kvectyin : Vy vector
CCC  -----
CCC  Output :  kvectoout : Vu vector (subset of Vo vector)
CCC  ------    kposcoefobs : Hyu operator (subset of Hyo operator)
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(in) :: kvectyin
      BIGREAL, dimension(:), intent(out) :: kvectoout
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jposize,jpitpsize
      INTEGER :: jy,jo,jodeb,jofin,jitp,jpitpfin
      INTEGER :: jobs,indobs,inddbs
      LOGICAL :: ltransfert
      BIGREAL :: sxy_moy,sxy_ect
CC----------------------------------------------------------------------
C
C Check size of input arrays
      jposize = size(kvectoout,1)
      jpitpsize = size(kposcoefobs,2)
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 1000
C
C Check if Vy and Vo vectors have been centered and reduced
C using the same statistics. Decide to make correction
C if necessary.
      ltransfert=.FALSE.
      IF (lmoyect) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltransfert=(ltransfert.OR.
     $           (obs_ect(indobs,inddbs).NE.dta_ect(indobs)))
            ltransfert=(ltransfert.OR.
     $           (obs_moy(indobs,inddbs).NE.dta_moy(indobs)))
         ENDDO
      ENDIF
C
C Extract Vu vector from Vy vector using Hyu operator
C Centering/reduction correction is impossible in such case
      IF (.NOT.ltransfert) THEN
         DO jo=1,jposize
            kvectoout(jo) = kvectyin(MAX(1,kposcoefobs(jo,1)%pos))
     $           *kposcoefobs(jo,1)%coef
         ENDDO
         DO jitp=2,jpitpsize
            DO jo=1,jposize
               kvectoout(jo) = kvectoout(jo) +
     $              kvectyin(MAX(1,kposcoefobs(jo,jitp)%pos))
     $              *kposcoefobs(jo,jitp)%coef
            ENDDO
         ENDDO
      ELSE
         GOTO 102
      ENDIF
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilmkh','mkhytou')
C
 102  WRITE (texterror,*) 'Cannot extract subset of Vo vector from ',
     $     'Vy vector: centering/reduction statistics must be the same'
      CALL printerror2(0,102,3,'utilmkh','mkhytou',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkhotoy(kvectoin,kvectyout,kposcoefobs,kspvaly,
     $     kcoeflimite,kvectynbpoint)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Reconstruct (roughly) Vy vector from Vo vector
CCC  -------   This is a (very) crude inverse of Hyo operator
CCC
CCC  Method :  Fill output Vy vector with special values
CCC  ------    Count the number of observations available
CCC            to estimate each Vy value
CCC            Only observations such that interpolation coefficient
CCC            (in Hyo) is larger than specified value (limit coefficient)
CCC            are taken into account
CCC            Estimate Vy vector from Vo vector, by making the
CCC            average of available observations
CCC
CCC  Input :   kvectoin : Vo vector
CCC  -----     kposcoefobs : Hyo operator
CCC            kspvaly : special value
CCC            kcoeflimite : limit coefficient (inside ]0,1])
CCC           
CCC  Output :  kvectyout : Vy vector
CCC  ------    kvectynbpoint : number of observations used
CCC                            to compute each Vy value
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpyend
      use mod_mask
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(in) :: kvectoin
      BIGREAL, dimension(:), intent(out) :: kvectyout
      TYPE (type_poscoef), dimension(:,:), intent(in) :: kposcoefobs
      BIGREAL, intent(in) :: kspvaly
      BIGREAL, intent(in), optional :: kcoeflimite
      INTEGER, dimension(:), intent(out), optional :: kvectynbpoint      
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable :: vectycompteur
      INTEGER :: allocok,jpysize,jposize,jpitpsize
      INTEGER :: jy,jydeb,jyfin,jo,jodeb,jofin,jitp,jpitpfin
      INTEGER :: jobs,indobs,inddbs,jdta,inddta
      LOGICAL :: ltransfert
      BIGREAL :: sxy_moy,sxy_ect,epsilon,coeflimite
CC----------------------------------------------------------------------
C
C Set default limit coefficient value
      epsilon=FREAL(0.001)
      IF (present(kcoeflimite)) THEN
         coeflimite=kcoeflimite
      ELSE
         coeflimite=FREAL(0.75)
      ENDIF
C
C Check size of input arrays
      jposize = size(kvectoin,1)
      jpitpsize = size(kposcoefobs,2)
      jpysize=jpyend
      IF (jposize.NE.size(kposcoefobs,1)) GOTO 102
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         IF (jpitpsize.LT.obs_itp(indobs,inddbs)) GOTO 102
      ENDDO
C
C --- allocation vectycompteur
      allocate ( vectycompteur(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vectycompteur(:) = FREAL(0.0)
C
C Check if Vy and Vo vectors have been centered and reduced
C using the same statistics. Decide to make correction
C if necessary.
      ltransfert=.FALSE.
      IF (lmoyect) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltransfert=(ltransfert.OR.
     $           (obs_ect(indobs,inddbs).NE.dta_ect(indobs)))
            ltransfert=(ltransfert.OR.
     $           (obs_moy(indobs,inddbs).NE.dta_moy(indobs)))
         ENDDO
      ENDIF
C
C Fill output Vy vector with special values
C Center/reduce special value if necessary (!)
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
C
C Count the number of observations available
C to estimate each Vy value
C Only observations such that interpolation coefficient
C (in Hyo) is larger than specified value (limit coefficient)
C are taken into account
C Where there is at least one observation available,
C initialize output Vy vector to zero
      DO jobs=1,obsend
         indobs=obs_ord(jobs)
         inddbs=obsnord(jobs)
         jodeb=obs_ind(indobs,inddbs)
         jofin=jodeb-1+obs_nbr(indobs,inddbs)
         jpitpfin=obs_itp(indobs,inddbs)
         DO jitp=1,jpitpfin
         DO jo=jodeb,jofin
            vectycompteur(kposcoefobs(jo,jitp)%pos) = 
     $        (FREAL(MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite)))+
     $        FREAL(1-MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite)))
     $           *epsilon) +
     $           vectycompteur(kposcoefobs(jo,jitp)%pos)
         ENDDO
         ENDDO
         DO jitp=1,jpitpfin
         DO jo=jodeb,jofin
            IF (vectycompteur(kposcoefobs(jo,jitp)%pos).NE.FREAL(0.0)) 
     $           kvectyout(kposcoefobs(jo,jitp)%pos)=FREAL(0.0)
         ENDDO
         ENDDO
      ENDDO
C
C Estimate Vy vector from Vo vector, by making the
C average of available observations
C Make centering/reduction correction if necessary
      IF (ltransfert) THEN
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            sxy_moy=(dta_moy(indobs)-obs_moy(indobs,inddbs))
     $           /dta_ect(indobs)
            sxy_ect=dta_ect(indobs)/obs_ect(indobs,inddbs)
            jodeb=obs_ind(indobs,inddbs)
            jofin=jodeb-1+obs_nbr(indobs,inddbs)
            jpitpfin=obs_itp(indobs,inddbs)
            DO jitp=1,jpitpfin
            DO jo=jodeb,jofin
               kvectyout(kposcoefobs(jo,jitp)%pos)= 
     $              kvectyout(kposcoefobs(jo,jitp)%pos) +
     $              ((( kvectoin(jo) *
     $         (FREAL(MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite)))+
     $         FREAL(1-MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite)))
     $              *epsilon) )
     $              / sxy_ect-sxy_moy) / 
     $              vectycompteur(kposcoefobs(jo,jitp)%pos) )
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
               kvectyout(kposcoefobs(jo,jitp)%pos)= 
     $              kvectyout(kposcoefobs(jo,jitp)%pos) +
     $              ( kvectoin(jo) *
     $         (FREAL(MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite)))+
     $         FREAL(1-MIN(1,int(kposcoefobs(jo,jitp)%coef+coeflimite)))
     $              *epsilon) )
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
C
      IF (allocated(vectycompteur)) deallocate (vectycompteur)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilmkh','mkhotoy')
 1001 CALL printerror2(0,1001,3,'utilmkh','mkhotoy')
C
 102  WRITE (texterror,*) 'Incompatible input array sizes'
      CALL printerror2(0,102,1,'utilmkh','mkhotoy',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilmkh
