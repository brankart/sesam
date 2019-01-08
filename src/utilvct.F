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
C ---                   UTILVCT.F                               ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12 (C.E. Testut)                        ---
C --- modification : 99-05 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                        ---
C --- modification : 03-04 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE mk4vct   : Paste 2D slice (single precision) of
C ---                       variable field on Vx or Vy segment
C --- SUBROUTINE mk8vct   : Paste 2D slice (double precision) of
C ---                       variable field on Vx or Vy segment
C --- SUBROUTINE unmk4vct : Extract 2D slice (single precision) of
C ---                       variable field from Vx or Vy segment
C --- SUBROUTINE unmk8vct : Extract 2D slice (double precision) of
C ---                       variable field from Vx or Vy segment
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilvct
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC mk4vct,mk8vct,unmk4vct,unmk8vct

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mk4vct(kvectsout,kptabij,jkin,jtin,
     $     kjsxy,ksompartsxynbr,kflagxyo) 
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Paste 2D slice (single precision) of variable field
CCC  -------   on Vx or Vy 1D vector segment. Exclude masked values.
CCC            Perform centering/reduction if required
CCC
CCC  Method :  Loop over 2D array dimensions
CCC  ------    Paste umasked values on 1D Vx or Vy segment
CCC
CCC  Input :   kptabij        : 2D slice of variable field
CCC  -----     jkin           : Index of 2D slice in 3rd dimension (level)
CCC            jtin           : Index of 2D slice in 4th dimension (time)
CCC            kjsxy          : Index of variable field
CCC            kflagxyo       : Vector type (1=Vx,2=Vy)
CCC
CCC  Output :  kvectsout      : Segment of 1D Vx or Vy vector
CCC  ------    ksompartsxynbr : Size of the segment
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(out) :: kvectsout
      BIGREAL4, dimension(:,:), intent(in) :: kptabij
      INTEGER, intent(in) :: jkin,jtin,kjsxy,kflagxyo
      INTEGER, intent(out) :: ksompartsxynbr
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpssize,kjpifin,kjpjfin
      INTEGER :: ji, jj, js
      INTEGER :: jpifin, jpjfin
      INTEGER :: indsxy,indmsk
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj
      BIGREAL :: sxy_moy, sxy_ect
CC----------------------------------------------------------------------
C
C Set size of input arrays
      jpssize=size(kvectsout,1)
      kjpifin=size(kptabij,1)
      kjpjfin=size(kptabij,2)
C
C Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE (1)
C --- var
         indsxy  = var_ord(kjsxy)
         sxy_dim = var_dim(indsxy)
         sxy_jpi = var_jpi(indsxy)
         sxy_jpj = var_jpj(indsxy)
         sxy_moy = var_moy(indsxy)
         sxy_ect = var_ect(indsxy)
         indmsk=kjsxy-1
      CASE (2)
C --- dta
         indsxy  = dta_ord(kjsxy)
         sxy_dim = dta_dim(indsxy)
         sxy_jpi = dta_jpi(indsxy)
         sxy_jpj = dta_jpj(indsxy)
         sxy_moy = dta_moy(indsxy)
         sxy_ect = dta_ect(indsxy)
         indmsk=kjsxy+varend-1
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Set 2D array size
      SELECT CASE (sxy_dim)
      CASE (1)
C ==>  1D
         jpifin = sxy_jpi
         jpjfin = 1
      CASE (2,3,4)
C ==>  2D
         jpifin = sxy_jpi
         jpjfin = sxy_jpj
      CASE DEFAULT
         GOTO 1000
      END SELECT
      js = 0
      ksompartsxynbr = 0
C
C Loop over 2D array dimensions
C Paste umasked values on 1D Vx or Vy segment
C Perform centering/reduction if required
      IF (lmoyect) THEN
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) THEN
               js=js+1
               kvectsout(js)=(FREAL(kptabij(ji,jj))
     $              -sxy_moy)/sxy_ect
            ENDIF
         ENDDO
         ENDDO
      ELSE
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) THEN
               js=js+1
               kvectsout(js)=FREAL(kptabij(ji,jj))
            ENDIF
         ENDDO
         ENDDO
      ENDIF
C
      IF (js.GT.jpssize) GOTO 1000
      ksompartsxynbr = js
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilvct','mk4vct')
 1001 CALL printerror2(0,1001,3,'utilvct','mk4vct')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mk8vct(kvectsout,kptabij,jkin,jtin,
     $     kjsxy,ksompartsxynbr,kflagxyo) 
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Paste 2D slice (double precision) of variable field
CCC  -------   on Vx or Vy 1D vector segment. Exclude masked values.
CCC            Perform centering/reduction if required
CCC
CCC  Method :  Loop over 2D array dimensions
CCC  ------    Paste umasked values on 1D Vx or Vy segment
CCC
CCC  Input :   kptabij        : 2D slice of variable field
CCC  -----     jkin           : Index of 2D slice in 3rd dimension (level)
CCC            jtin           : Index of 2D slice in 4th dimension (time)
CCC            kjsxy          : Index of variable field
CCC            kflagxyo       : Vector type (1=Vx,2=Vy)
CCC
CCC  Output :  kvectsout      : Segment of 1D Vx or Vy vector
CCC  ------    ksompartsxynbr : Size of the segment
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(out) :: kvectsout
      BIGREAL8, dimension(:,:), intent(in) :: kptabij
      INTEGER, intent(in) :: jkin,jtin,kjsxy,kflagxyo
      INTEGER, intent(out) :: ksompartsxynbr
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpssize,kjpifin,kjpjfin
      INTEGER :: ji, jj, js
      INTEGER :: jpifin, jpjfin
      INTEGER :: indsxy,indmsk
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj
      BIGREAL :: sxy_moy, sxy_ect
CC----------------------------------------------------------------------
C
C Set size of input arrays
      jpssize=size(kvectsout,1)
      kjpifin=size(kptabij,1)
      kjpjfin=size(kptabij,2)
C
C Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE (1)
C --- var
         indsxy  = var_ord(kjsxy)
         sxy_dim = var_dim(indsxy)
         sxy_jpi = var_jpi(indsxy)
         sxy_jpj = var_jpj(indsxy)
         sxy_moy = var_moy(indsxy)
         sxy_ect = var_ect(indsxy)
         indmsk=kjsxy-1
      CASE (2)
C --- dta
         indsxy  = dta_ord(kjsxy)
         sxy_dim = dta_dim(indsxy)
         sxy_jpi = dta_jpi(indsxy)
         sxy_jpj = dta_jpj(indsxy)
         sxy_moy = dta_moy(indsxy)
         sxy_ect = dta_ect(indsxy)
         indmsk=kjsxy+varend-1
      CASE(3)
         GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Set 2D array size
      SELECT CASE (sxy_dim)
      CASE (1)
C ==>  1D
         jpifin = sxy_jpi
         jpjfin = 1
      CASE (2,3,4)
C ==>  2D
         jpifin = sxy_jpi
         jpjfin = sxy_jpj
      CASE DEFAULT
         GOTO 1000
      END SELECT
      js = 0
      ksompartsxynbr = 0
C
C Loop over 2D array dimensions
C Paste umasked values on 1D Vx or Vy segment
C Perform centering/reduction if required
      IF (lmoyect) THEN
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0)
     $           THEN
               js=js+1
               kvectsout(js)=(FREAL(kptabij(ji,jj))
     $              -sxy_moy)/sxy_ect
            ENDIF
         ENDDO
         ENDDO
      ELSE
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0)
     $           THEN
               js=js+1
               kvectsout(js)=FREAL(kptabij(ji,jj))
            ENDIF
         ENDDO
         ENDDO
      ENDIF
C
      IF (js.GT.jpssize) GOTO 1000
      ksompartsxynbr = js
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilvct','mk8vct')
 1001 CALL printerror2(0,1001,3,'utilvct','mk8vct')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE unmk4vct(kvectsin,kptabij,jkin,jtin,
     $     kjsxy,ksompartsxynbr,kspval,kflagxyo) 
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Extract 2D slice (single precision) of variable field
CCC  -------   from Vx or Vy 1D vector segment. Flag masked values.
CCC            Perform centering/reduction if required.
CCC
CCC  Method :  Loop over 2D array dimensions
CCC  ------    Fill umasked values with 1D Vx or Vy segment
CCC
CCC  Input :   kvectsin       : Segment of 1D Vx or Vy vector
CCC  -----     jkin           : Index of 2D slice in 3rd dimension (level)
CCC            jtin           : Index of 2D slice in 4th dimension (time)
CCC            kjsxy          : Index of variable field
CCC            kflagxyo       : Vector type (1=Vx,2=Vy)
CCC
CCC  Output :  kptabij        : 2D slice of variable field
CCC  ------    ksompartsxynbr : Size of the segment
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(in) :: kvectsin
      BIGREAL4, dimension(:,:), intent(out) :: kptabij
      INTEGER, intent(in) :: jkin,jtin,kjsxy,kflagxyo
      INTEGER, intent(out) :: ksompartsxynbr
      BIGREAL4, intent(in) :: kspval
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpisize,jpjsize
      INTEGER :: ji, jj, js
      INTEGER :: jpifin, jpjfin
      INTEGER :: indsxy,indmsk
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj
      BIGREAL :: sxy_moy, sxy_ect
CC----------------------------------------------------------------------
C
C Set size of input arrays
      jpisize = size(kptabij,1)
      jpjsize = size(kptabij,2)
C
C Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE (1)
C --- var
         indsxy  = var_ord(kjsxy)
         sxy_dim = var_dim(indsxy)
         sxy_jpi = var_jpi(indsxy)
         sxy_jpj = var_jpj(indsxy)
         sxy_moy = var_moy(indsxy)
         sxy_ect = var_ect(indsxy)
         indmsk=kjsxy-1
      CASE (2)
C --- dta
         indsxy  = dta_ord(kjsxy)
         sxy_dim = dta_dim(indsxy)
         sxy_jpi = dta_jpi(indsxy)
         sxy_jpj = dta_jpj(indsxy)
         sxy_moy = dta_moy(indsxy)
         sxy_ect = dta_ect(indsxy)
         indmsk=kjsxy+varend-1
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Set 2D array size
      SELECT CASE (sxy_dim)
      CASE (1)
C ==>  1D
         jpifin = sxy_jpi
         jpjfin = 1
      CASE (2,3,4)
C ==>  2D
         jpifin = sxy_jpi
         jpjfin = sxy_jpj
      CASE DEFAULT
         GOTO 1000
      END SELECT
      js = 0
      ksompartsxynbr = 0
C
C Loop over 2D array dimensions
C Fill unmasked value with 1D Vx or Vy segment
C Flag masked values
C Perform centering/reduction if required
      IF (lmoyect) THEN
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) THEN
               js=js+1
               kptabij(ji,jj)=FREAL4(kvectsin(js))
     $              *sxy_ect+sxy_moy
            ELSE
               kptabij(ji,jj)=kspval
            ENDIF
         ENDDO
         ENDDO
      ELSE
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) THEN
               js=js+1
               kptabij(ji,jj)=FREAL4(kvectsin(js))
            ELSE
               kptabij(ji,jj)=kspval
            ENDIF
         ENDDO
         ENDDO
      ENDIF
C
      ksompartsxynbr = js
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilvct','unmk4vct')
 1001 CALL printerror2(0,1001,3,'utilvct','mk4vct')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE unmk8vct(kvectsin,kptabij,jkin,jtin,
     $     kjsxy,ksompartsxynbr,kspval,kflagxyo) 
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Extract 2D slice (double precision) of variable field
CCC  -------   from Vx or Vy 1D vector segment. Flag masked values.
CCC            Perform centering/reduction if required.
CCC
CCC  Method :  Loop over 2D array dimensions
CCC  ------    Fill umasked values with 1D Vx or Vy segment
CCC
CCC  Input :   kvectsin       : Segment of 1D Vx or Vy vector
CCC  -----     jkin           : Index of 2D slice in 3rd dimension (level)
CCC            jtin           : Index of 2D slice in 4th dimension (time)
CCC            kjsxy          : Index of variable field
CCC            kflagxyo       : Vector type (1=Vx,2=Vy)
CCC
CCC  Output :  kptabij        : 2D slice of variable field
CCC  ------    ksompartsxynbr : Size of the segment
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(in) :: kvectsin
      BIGREAL8, dimension(:,:), intent(out) :: kptabij
      INTEGER, intent(in) :: jkin,jtin,kjsxy,kflagxyo
      INTEGER, intent(out) :: ksompartsxynbr
      BIGREAL8, intent(in) :: kspval
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jpisize,jpjsize
      INTEGER :: ji, jj, js
      INTEGER :: jpifin, jpjfin
      INTEGER :: indsxy,indmsk
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj
      BIGREAL :: sxy_moy, sxy_ect
CC----------------------------------------------------------------------
C
C Set size of input arrays
      jpisize = size(kptabij,1)
      jpjsize = size(kptabij,2)
C
C Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE (1)
C --- var
         indsxy  = var_ord(kjsxy)
         sxy_dim = var_dim(indsxy)
         sxy_jpi = var_jpi(indsxy)
         sxy_jpj = var_jpj(indsxy)
         sxy_moy = var_moy(indsxy)
         sxy_ect = var_ect(indsxy)
         indmsk=kjsxy-1
      CASE (2)
C --- dta
         indsxy  = dta_ord(kjsxy)
         sxy_dim = dta_dim(indsxy)
         sxy_jpi = dta_jpi(indsxy)
         sxy_jpj = dta_jpj(indsxy)
         sxy_moy = dta_moy(indsxy)
         sxy_ect = dta_ect(indsxy)
         indmsk=kjsxy+varend-1
      CASE(3)
         GOTO 1000
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Set 2D array size
      SELECT CASE (sxy_dim)
      CASE (1)
C ==>  1D
         jpifin = sxy_jpi
         jpjfin = 1
      CASE (2,3,4)
C ==>  2D
         jpifin = sxy_jpi
         jpjfin = sxy_jpj
      CASE DEFAULT
         GOTO 1000
      END SELECT
      js = 0
      ksompartsxynbr = 0
C
C Loop over 2D array dimensions
C Fill unmasked value with 1D Vx or Vy segment
C Flag masked values
C Perform centering/reduction if required
      IF (lmoyect) THEN
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) THEN
               js=js+1
               kptabij(ji,jj)=FREAL8(kvectsin(js)
     $              *sxy_ect+sxy_moy)
            ELSE
               kptabij(ji,jj)=kspval
            ENDIF
         ENDDO
         ENDDO
      ELSE
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) THEN
               js=js+1
               kptabij(ji,jj)=FREAL8(kvectsin(js))
            ELSE
               kptabij(ji,jj)=kspval
            ENDIF
         ENDDO
         ENDDO
      ENDIF
C
      ksompartsxynbr = js
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'utilvct','unmk8vct')
 1001 CALL printerror2(0,1001,3,'utilvct','mk8vct')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilvct
