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
! ---                   UTILVCT.F90                               ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12 (C.E. Testut)                        ---
! --- modification : 99-05 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 03-04 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE mk4vct   : Paste 2D slice (single precision) of
! ---                       variable field on Vx or Vy segment
! --- SUBROUTINE mk8vct   : Paste 2D slice (double precision) of
! ---                       variable field on Vx or Vy segment
! --- SUBROUTINE unmk4vct : Extract 2D slice (single precision) of
! ---                       variable field from Vx or Vy segment
! --- SUBROUTINE unmk8vct : Extract 2D slice (double precision) of
! ---                       variable field from Vx or Vy segment
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilvct
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC mk4vct,mk8vct,unmk4vct,unmk8vct

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mk4vct(kvectsout,kptabij,jkin,jtin, &
     &     kjsxy,ksompartsxynbr,kflagxyo) 
!---------------------------------------------------------------------
!
!  Purpose : Paste 2D slice (single precision) of variable field
!  -------   on Vx or Vy 1D vector segment. Exclude masked values.
!            Perform centering/reduction if required
!
!  Method :  Loop over 2D array dimensions
!  ------    Paste umasked values on 1D Vx or Vy segment
!
!  Input :   kptabij        : 2D slice of variable field
!  -----     jkin           : Index of 2D slice in 3rd dimension (level)
!            jtin           : Index of 2D slice in 4th dimension (time)
!            kjsxy          : Index of variable field
!            kflagxyo       : Vector type (1=Vx,2=Vy)
!
!  Output :  kvectsout      : Segment of 1D Vx or Vy vector
!  ------    ksompartsxynbr : Size of the segment
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(out) :: kvectsout
      BIGREAL4, dimension(:,:), intent(in) :: kptabij
      INTEGER, intent(in) :: jkin,jtin,kjsxy,kflagxyo
      INTEGER, intent(out) :: ksompartsxynbr
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpssize,kjpifin,kjpjfin
      INTEGER :: ji, jj, js
      INTEGER :: jpifin, jpjfin
      INTEGER :: indsxy,indmsk
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj
      BIGREAL :: sxy_moy, sxy_ect
!----------------------------------------------------------------------
!
! Set size of input arrays
      jpssize=size(kvectsout,1)
      kjpifin=size(kptabij,1)
      kjpjfin=size(kptabij,2)
!
! Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE (1)
! --- var
         indsxy  = var_ord(kjsxy)
         sxy_dim = var_dim(indsxy)
         sxy_jpi = var_jpi(indsxy)
         sxy_jpj = var_jpj(indsxy)
         sxy_moy = var_moy(indsxy)
         sxy_ect = var_ect(indsxy)
         indmsk=kjsxy-1
      CASE (2)
! --- dta
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
!
! Set 2D array size
      SELECT CASE (sxy_dim)
      CASE (1)
! ==>  1D
         jpifin = sxy_jpi
         jpjfin = 1
      CASE (2,3,4)
! ==>  2D
         jpifin = sxy_jpi
         jpjfin = sxy_jpj
      CASE DEFAULT
         GOTO 1000
      END SELECT
      js = 0
      ksompartsxynbr = 0
!
! Loop over 2D array dimensions
! Paste umasked values on 1D Vx or Vy segment
! Perform centering/reduction if required
      IF (lmoyect) THEN
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) THEN
               js=js+1
               kvectsout(js)=(FREAL(kptabij(ji,jj)) &
     &              -sxy_moy)/sxy_ect
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
!
      IF (js.GT.jpssize) GOTO 1000
      ksompartsxynbr = js
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilvct','mk4vct')
 1001 CALL printerror2(0,1001,3,'utilvct','mk4vct')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mk8vct(kvectsout,kptabij,jkin,jtin, &
     &     kjsxy,ksompartsxynbr,kflagxyo) 
!---------------------------------------------------------------------
!
!  Purpose : Paste 2D slice (double precision) of variable field
!  -------   on Vx or Vy 1D vector segment. Exclude masked values.
!            Perform centering/reduction if required
!
!  Method :  Loop over 2D array dimensions
!  ------    Paste umasked values on 1D Vx or Vy segment
!
!  Input :   kptabij        : 2D slice of variable field
!  -----     jkin           : Index of 2D slice in 3rd dimension (level)
!            jtin           : Index of 2D slice in 4th dimension (time)
!            kjsxy          : Index of variable field
!            kflagxyo       : Vector type (1=Vx,2=Vy)
!
!  Output :  kvectsout      : Segment of 1D Vx or Vy vector
!  ------    ksompartsxynbr : Size of the segment
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(out) :: kvectsout
      BIGREAL8, dimension(:,:), intent(in) :: kptabij
      INTEGER, intent(in) :: jkin,jtin,kjsxy,kflagxyo
      INTEGER, intent(out) :: ksompartsxynbr
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpssize,kjpifin,kjpjfin
      INTEGER :: ji, jj, js
      INTEGER :: jpifin, jpjfin
      INTEGER :: indsxy,indmsk
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj
      BIGREAL :: sxy_moy, sxy_ect
!----------------------------------------------------------------------
!
! Set size of input arrays
      jpssize=size(kvectsout,1)
      kjpifin=size(kptabij,1)
      kjpjfin=size(kptabij,2)
!
! Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE (1)
! --- var
         indsxy  = var_ord(kjsxy)
         sxy_dim = var_dim(indsxy)
         sxy_jpi = var_jpi(indsxy)
         sxy_jpj = var_jpj(indsxy)
         sxy_moy = var_moy(indsxy)
         sxy_ect = var_ect(indsxy)
         indmsk=kjsxy-1
      CASE (2)
! --- dta
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
!
! Set 2D array size
      SELECT CASE (sxy_dim)
      CASE (1)
! ==>  1D
         jpifin = sxy_jpi
         jpjfin = 1
      CASE (2,3,4)
! ==>  2D
         jpifin = sxy_jpi
         jpjfin = sxy_jpj
      CASE DEFAULT
         GOTO 1000
      END SELECT
      js = 0
      ksompartsxynbr = 0
!
! Loop over 2D array dimensions
! Paste umasked values on 1D Vx or Vy segment
! Perform centering/reduction if required
      IF (lmoyect) THEN
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) &
     &           THEN
               js=js+1
               kvectsout(js)=(FREAL(kptabij(ji,jj)) &
     &              -sxy_moy)/sxy_ect
            ENDIF
         ENDDO
         ENDDO
      ELSE
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) &
     &           THEN
               js=js+1
               kvectsout(js)=FREAL(kptabij(ji,jj))
            ENDIF
         ENDDO
         ENDDO
      ENDIF
!
      IF (js.GT.jpssize) GOTO 1000
      ksompartsxynbr = js
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilvct','mk8vct')
 1001 CALL printerror2(0,1001,3,'utilvct','mk8vct')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE unmk4vct(kvectsin,kptabij,jkin,jtin, &
     &     kjsxy,ksompartsxynbr,kspval,kflagxyo) 
!---------------------------------------------------------------------
!
!  Purpose : Extract 2D slice (single precision) of variable field
!  -------   from Vx or Vy 1D vector segment. Flag masked values.
!            Perform centering/reduction if required.
!
!  Method :  Loop over 2D array dimensions
!  ------    Fill umasked values with 1D Vx or Vy segment
!
!  Input :   kvectsin       : Segment of 1D Vx or Vy vector
!  -----     jkin           : Index of 2D slice in 3rd dimension (level)
!            jtin           : Index of 2D slice in 4th dimension (time)
!            kjsxy          : Index of variable field
!            kflagxyo       : Vector type (1=Vx,2=Vy)
!
!  Output :  kptabij        : 2D slice of variable field
!  ------    ksompartsxynbr : Size of the segment
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(in) :: kvectsin
      BIGREAL4, dimension(:,:), intent(out) :: kptabij
      INTEGER, intent(in) :: jkin,jtin,kjsxy,kflagxyo
      INTEGER, intent(out) :: ksompartsxynbr
      BIGREAL4, intent(in) :: kspval
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpisize,jpjsize
      INTEGER :: ji, jj, js
      INTEGER :: jpifin, jpjfin
      INTEGER :: indsxy,indmsk
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj
      BIGREAL :: sxy_moy, sxy_ect
!----------------------------------------------------------------------
!
! Set size of input arrays
      jpisize = size(kptabij,1)
      jpjsize = size(kptabij,2)
!
! Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE (1)
! --- var
         indsxy  = var_ord(kjsxy)
         sxy_dim = var_dim(indsxy)
         sxy_jpi = var_jpi(indsxy)
         sxy_jpj = var_jpj(indsxy)
         sxy_moy = var_moy(indsxy)
         sxy_ect = var_ect(indsxy)
         indmsk=kjsxy-1
      CASE (2)
! --- dta
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
!
! Set 2D array size
      SELECT CASE (sxy_dim)
      CASE (1)
! ==>  1D
         jpifin = sxy_jpi
         jpjfin = 1
      CASE (2,3,4)
! ==>  2D
         jpifin = sxy_jpi
         jpjfin = sxy_jpj
      CASE DEFAULT
         GOTO 1000
      END SELECT
      js = 0
      ksompartsxynbr = 0
!
! Loop over 2D array dimensions
! Fill unmasked value with 1D Vx or Vy segment
! Flag masked values
! Perform centering/reduction if required
      IF (lmoyect) THEN
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) THEN
               js=js+1
               kptabij(ji,jj)=FREAL4(kvectsin(js)) &
     &              *sxy_ect+sxy_moy
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
!
      ksompartsxynbr = js
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilvct','unmk4vct')
 1001 CALL printerror2(0,1001,3,'utilvct','mk4vct')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE unmk8vct(kvectsin,kptabij,jkin,jtin, &
     &     kjsxy,ksompartsxynbr,kspval,kflagxyo) 
!---------------------------------------------------------------------
!
!  Purpose : Extract 2D slice (double precision) of variable field
!  -------   from Vx or Vy 1D vector segment. Flag masked values.
!            Perform centering/reduction if required.
!
!  Method :  Loop over 2D array dimensions
!  ------    Fill umasked values with 1D Vx or Vy segment
!
!  Input :   kvectsin       : Segment of 1D Vx or Vy vector
!  -----     jkin           : Index of 2D slice in 3rd dimension (level)
!            jtin           : Index of 2D slice in 4th dimension (time)
!            kjsxy          : Index of variable field
!            kflagxyo       : Vector type (1=Vx,2=Vy)
!
!  Output :  kptabij        : 2D slice of variable field
!  ------    ksompartsxynbr : Size of the segment
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_mask
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(in) :: kvectsin
      BIGREAL8, dimension(:,:), intent(out) :: kptabij
      INTEGER, intent(in) :: jkin,jtin,kjsxy,kflagxyo
      INTEGER, intent(out) :: ksompartsxynbr
      BIGREAL8, intent(in) :: kspval
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jpisize,jpjsize
      INTEGER :: ji, jj, js
      INTEGER :: jpifin, jpjfin
      INTEGER :: indsxy,indmsk
      INTEGER :: sxy_dim,sxy_jpi,sxy_jpj
      BIGREAL :: sxy_moy, sxy_ect
!----------------------------------------------------------------------
!
! Set size of input arrays
      jpisize = size(kptabij,1)
      jpjsize = size(kptabij,2)
!
! Get variable characteristics from SESAM configuration
      SELECT CASE(kflagxyo)
      CASE (1)
! --- var
         indsxy  = var_ord(kjsxy)
         sxy_dim = var_dim(indsxy)
         sxy_jpi = var_jpi(indsxy)
         sxy_jpj = var_jpj(indsxy)
         sxy_moy = var_moy(indsxy)
         sxy_ect = var_ect(indsxy)
         indmsk=kjsxy-1
      CASE (2)
! --- dta
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
!
! Set 2D array size
      SELECT CASE (sxy_dim)
      CASE (1)
! ==>  1D
         jpifin = sxy_jpi
         jpjfin = 1
      CASE (2,3,4)
! ==>  2D
         jpifin = sxy_jpi
         jpjfin = sxy_jpj
      CASE DEFAULT
         GOTO 1000
      END SELECT
      js = 0
      ksompartsxynbr = 0
!
! Loop over 2D array dimensions
! Fill unmasked value with 1D Vx or Vy segment
! Flag masked values
! Perform centering/reduction if required
      IF (lmoyect) THEN
         DO jj=1,jpjfin
         DO ji=1,jpifin
            IF (IBITS(mask(ji,jj,jkin,jtin),indmsk,1).NE.0) THEN
               js=js+1
               kptabij(ji,jj)=FREAL8(kvectsin(js) &
     &              *sxy_ect+sxy_moy)
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
!
      ksompartsxynbr = js
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'utilvct','unmk8vct')
 1001 CALL printerror2(0,1001,3,'utilvct','mk8vct')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilvct
