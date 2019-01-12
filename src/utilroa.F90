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
! ---                  UTILROA.F90                               ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 98-06 (C.E. Testut)                        ---
! --- modification : 99-05 (C.E. Testut)                        ---
! --- modification : 99-12 (J.M. Brankart)                      ---
! --- modification : 00-09 (J.M. Brankart)                      ---
! --- modification : 07-11 (J.M. Brankart)                      ---
! --- modification : 08-10 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- Set of routines used in the Reduced Order Analysis (ROA) algorithm
! --- 
! --- SUBROUTINE  algoker_u : Compute SVD decomposition of the ROA kernel matrix
! --- SUBROUTINE  algocoef_u : Compute innovation in the reduced eigenspace
! --- SUBROUTINE  algocoef_u_sf : Compute the ROA coefficients (U method)
! --- SUBROUTINE  algomatc_u : Compute matrix transforming Sf into Sa
! --- SUBROUTINE  algocoeflimit_u : Saturate the ROA coefficients (U method)
! --- SUBROUTINE  algocalcGi : Compute Gamma matrices for obs segments
! --- SUBROUTINE  algoker_Gi : SVD of the sum of Gamma matrices for obs segments
! --- SUBROUTINE  algocalcDi : Compute delta vectors for obs segments
! --- SUBROUTINE  algosum_Di : Sum up the delta vectors for obs segments
! --- SUBROUTINE  mkyorms : Fill observation error standard deviation vector
! ---                       using parameters from SESAM configuration
! --- SUBROUTINE  mkrs : Compute eigenvalues and eigenvectors
! ---                    of a symmetric real matrix
! --- SUBROUTINE  mkcoeftodta : Save local z-vector in dta file
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilroa
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC algoker_u, algocoef_u, algocoef_u_sf, algomatc_u
      PUBLIC algocoeflimit_u, algocalcGi, algoker_Gi, algocalcDi
      PUBLIC algosum_Di, mkyorms, mkrs, mkcoeftodta, algocalciRi

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoker_u(kbasesr,kvectssqrdiagRi,kmatUrr,klambda, &
     &                     kbasesq,kmatBqr)
!---------------------------------------------------------------------
!
!  Purpose : Compute SVD decomposition of the ROA kernel matrix
!  -------
!  Method : (U,lambda)  = SVD ( trsp(HSf) inv(R) (HSf) )
!  ------
!  Input :  kbasesr         : HSf
!  -----    kvectssqrdiagRi : inv(R)
!           kforgfact       : rho
!           kbasesq         : SQRT of obs. inv. correl. matrix
!
!  Output : kmatUrr         : U
!  ------   klambda         : Lambda
!           kmatBqr         : SQRT of obs. correl. matrix in reduced space
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
       BIGREAL, dimension(:,:), intent(in) :: kbasesr
       BIGREAL, dimension(:), intent(in) :: kvectssqrdiagRi
       BIGREAL, dimension(:,:), intent(out) :: kmatUrr
       BIGREAL, dimension(:), intent(out) :: klambda
       BIGREAL, dimension(:,:), optional, intent(in) :: kbasesq
       BIGREAL, dimension(:,:), optional, intent(out) :: kmatBqr
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jpssize,jprsize,jpqsize
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      INTEGER :: jqbasdeb,jqbasfin,jqbas,jqbas1,jqbas2
      INTEGER :: jqmatdeb,jqmatfin,jqmat,jqmat1,jqmat2
      INTEGER :: nvpnull
      BIGREAL, dimension(:,:), allocatable :: matsr
      BIGREAL, dimension(:,:), allocatable :: matrr
!----------------------------------------------------------------------
!
      jpssize=size(kbasesr,1)
      jprsize=size(kbasesr,2)
      IF (present(kbasesq)) THEN
         jpqsize=size(kbasesq,2)
      ENDIF
!
! --- allocation matsr
      allocate ( matsr(1:jpssize,1:jprsize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      matsr(:,:) = FREAL(0.0)
! --- allocation matrr
      allocate ( matrr(1:jprsize,1:jprsize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      matrr(:,:) = FREAL(0.0)
!
! -0.- Initialisation :
! ---------------------
!
      IF (jpssize.NE.size(kvectssqrdiagRi,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
!
      IF (present(kbasesq)) THEN
         IF (jpssize.NE.size(kbasesq,1)) GOTO 1000
         IF (.NOT.present(kmatBqr)) GOTO 1000
         IF (jpqsize.NE.size(kmatBqr,1)) GOTO 1000
         IF (jprsize.NE.size(kmatBqr,2)) GOTO 1000
      ENDIF
!
      jrbasdeb=2
      jrbasfin=jprsize
      jrmatdeb=1
      jrmatfin=jprsize-1
!
      IF (present(kbasesq)) THEN
         jqbasdeb=2
         jqbasfin=jpqsize
         jqmatdeb=1
         jqmatfin=jpqsize-1
      ENDIF
!
! -1.- Compute the ROA kernel matrix (trsp((HSf)inv(R)(HSf))
! -----------------------------------------------------------
!
      DO jrbas=jrbasdeb,jrbasfin
         jrmat = jrbas - 1
         matsr(:,jrmat) = kbasesr(:,jrbas) * kvectssqrdiagRi(:)
      ENDDO
!
      IF (.NOT.present(kbasesq)) THEN
         DO jrmat2=jrmatdeb,jrmatfin
            DO jrmat1=jrmat2,jrmatfin
               matrr(jrmat1,jrmat2) =  &
     &              DOT_PRODUCT(matsr(:,jrmat1),matsr(:,jrmat2))
            ENDDO
            matrr(jrmat2,jrmat2:jrmatfin)=matrr(jrmat2:jrmatfin,jrmat2)
         ENDDO
      ELSE
         DO jrmat2=jrmatdeb,jrmatfin
            DO jqbas1=jqbasdeb,jqbasfin
               jqmat1 = jqbas1 - 1
               kmatBqr(jqmat1,jrmat2) =  &
     &              DOT_PRODUCT(kbasesq(:,jqbas1),matsr(:,jrmat2))
            ENDDO
         ENDDO
         DO jrmat2=jrmatdeb,jrmatfin
            DO jrmat1=jrmat2,jrmatfin
               matrr(jrmat1,jrmat2) =  &
     &              DOT_PRODUCT(kmatBqr(jqmatdeb:jqmatfin,jrmat1), &
     &              kmatBqr(jqmatdeb:jqmatfin,jrmat2))
            ENDDO
            matrr(jrmat2,jrmat2:jrmatfin)=matrr(jrmat2:jrmatfin,jrmat2)
         ENDDO
      ENDIF
!
! -2.- Compute SVD decomposition of C
! -----------------------------------
!
      nvpnull=0
      CALL mkrs (matrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin), &
     &     kmatUrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin), &
     &     klambda(jrbasdeb:jrbasfin),nvpnull)
!
! --- deallocate arrays
      IF (allocated(matsr)) deallocate(matsr)
      IF (allocated(matrr)) deallocate(matrr)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'utilroa','algoker_u')
 1001 CALL printerror2(0,1001,3,'utilroa','algoker_u')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocoef_u(kbasesr,kvectsinnov,kvectssqrdiagRi, &
     &                      kmatUrr,kdelta,kbasesq,kmatBqr)
!---------------------------------------------------------------------
!
!  Purpose : Compute innovation in the reduced eigenspace
!  -------
!  Method : delta = U^ trsp(HSf) inv(R) (y-Hxf)
!  ------
!  Input :  kbasesr         : HSf
!  -----    kvectsinnov     : y-Hxf
!           kvectssqrdiagRi : inv(R)
!           kmatUrr         : U
!           kbasesq         : SQRT of inv. obs. correl. matrix
!  Output : kdelta          : delta
!  ------   kmatBqr         : SQRT of obs. correl. matrix in reduced space
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
       BIGREAL, dimension(:,:), intent(in) :: kbasesr
       BIGREAL, dimension(:), intent(in) :: kvectsinnov
       BIGREAL, dimension(:), intent(in) :: kvectssqrdiagRi
       BIGREAL, dimension(:,:), intent(in) :: kmatUrr
       BIGREAL, dimension(:), intent(out) :: kdelta
       BIGREAL, dimension(:,:), optional, intent(in) :: kbasesq
       BIGREAL, dimension(:,:), optional, intent(in) :: kmatBqr
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jpssize,jprsize,jpqsize
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      INTEGER :: jqbasdeb,jqbasfin,jqbas,jqbas1,jqbas2
      INTEGER :: jqmatdeb,jqmatfin,jqmat,jqmat1,jqmat2
      BIGREAL, dimension(:), allocatable :: kcoefr1
      BIGREAL, dimension(:), allocatable :: kcoefq1
      BIGREAL, dimension(:), allocatable :: vects
!----------------------------------------------------------------------
!
      jpssize=size(kbasesr,1)
      jprsize=size(kbasesr,2)
      IF (present(kbasesq)) THEN
         jpqsize=size(kbasesq,2)
      ENDIF
! --- allocation kcoefr1
      allocate ( kcoefr1(1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      kcoefr1(:) = FREAL(0.0)
! --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
! --- allocation kcoefq1
      IF (present(kbasesq)) THEN
         allocate ( kcoefq1(1:jpqsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         kcoefq1(:) = FREAL(0.0)
      ENDIF
!
! -0.- Initialisation
! -------------------
!
      IF (jpssize.NE.size(kvectsinnov,1)) GOTO 1000
      IF (jpssize.NE.size(kvectssqrdiagRi,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
      IF (jprsize.NE.size(kdelta,1)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
!
      IF (present(kbasesq)) THEN
         IF (jpssize.NE.size(kbasesq,1)) GOTO 1000
         IF (.NOT.present(kmatBqr)) GOTO 1000
         IF (jpqsize.NE.size(kmatBqr,1)) GOTO 1000
         IF (jprsize.NE.size(kmatBqr,2)) GOTO 1000
      ENDIF
!
      jrbasdeb=2
      jrbasfin=jprend
      jrmatdeb=1
      jrmatfin=jprend-1
!
      IF (present(kbasesq)) THEN
         jqbasdeb=2
         jqbasfin=jpqsize
         jqmatdeb=1
         jqmatfin=jpqsize-1
      ENDIF
!
! -1.- Compute innovation vector and weight by observation error
! --------------------------------------------------------------
!
      IF (.NOT.present(kbasesq)) THEN
         vects(:) = kvectsinnov(:) * &
     &           kvectssqrdiagRi(:) * kvectssqrdiagRi(:)
      ELSE
         vects(:) = kvectsinnov(:) * kvectssqrdiagRi(:)
      ENDIF
!
! -2.- Compute innovation in reduced space
! ----------------------------------------
!
      IF (.NOT.present(kbasesq)) THEN
         DO jrbas=jrbasdeb,jrbasfin
            kcoefr1(jrbas)=DOT_PRODUCT(kbasesr(:,jrbas),vects(:))
         ENDDO
      ELSE
         DO jqbas=jqbasdeb,jqbasfin
            kcoefq1(jqbas)=DOT_PRODUCT(kbasesq(:,jqbas),vects(:))
         ENDDO
         DO jrbas=jrbasdeb,jrbasfin
            jrmat = jrbas - 1
            kcoefr1(jrbas)=DOT_PRODUCT(kmatBqr(jqmatdeb:jqmatfin,jrmat), &
     &                                kcoefq1(jqbasdeb:jqbasfin))
         ENDDO
      ENDIF
!
! -3.- Transform to analysis eigenpsace
! -------------------------------------
!
      DO jrbas=jrbasdeb,jrbasfin
         jrmat = jrbas - 1
         kdelta(jrbas)=DOT_PRODUCT(kmatUrr(jrmatdeb:jrmatfin,jrmat), &
     &             kcoefr1(jrbasdeb:jrbasfin))
      ENDDO
!
! --- deallocate arrays
      IF (allocated(kcoefr1)) deallocate (kcoefr1)
      IF (allocated(kcoefq1)) deallocate (kcoefq1)
      IF (allocated(vects)) deallocate(vects)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilroa','algocoef_u')
 1001 CALL printerror2(0,1001,3,'utilroa','algocoef_u')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocoef_u_sf(kcoefr,klambda,kdelta,kmatUrr,kforgfact)
!---------------------------------------------------------------------
!   
!  Purpose : Compute the ROA coefficients (U method)
!  -------
!  Method : c = U inv[rho I + Lambda] delta
!  ------
!  Input :  klambda         : Lambda
!  -----    kdelta          : delta
!           kmatUrr         : U
!           kforgfact       : rho
!  Output : kcoefr          : ROA coefficients
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(out) :: kcoefr
      BIGREAL, dimension(:), intent(in) :: klambda
      BIGREAL, dimension(:), intent(in) :: kdelta
      BIGREAL, dimension(:,:), intent(in) :: kmatUrr
      BIGREAL, intent(in) :: kforgfact
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok, jprsize
      INTEGER :: jrbasdeb, jrbasfin, jrbas
      INTEGER :: jrmatdeb, jrmatfin, jrmat
      BIGREAL, dimension(:), allocatable :: kcoefr1
!----------------------------------------------------------------------
!
      jprsize=size(kcoefr,1)
      IF (jprsize.NE.size(klambda,1)) GOTO 1000
      IF (jprsize.NE.size(kdelta,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
! --- allocation kcoefr1
      allocate ( kcoefr1(1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      kcoefr1(:) = FREAL(0.0)
!
! -0.- Initialisation :
! ---------------------
!
      jrbasdeb=2
      jrbasfin=jprend
      jrmatdeb=1
      jrmatfin=jprend-1
!
! -1.- Compute the ROA coefficients
! ----------------------------------
!
      kcoefr1(jrbasdeb:jrbasfin) = kdelta(jrbasdeb:jrbasfin) / &
     &             ( kforgfact + klambda(jrbasdeb:jrbasfin) )
!
! -2.- Transform the ROA coefficient back to the original basis
! --------------------------------------------------------------
!
      DO jrbas=jrbasdeb,jrbasfin
         jrmat = jrbas - 1
         kcoefr(jrbas) = DOT_PRODUCT(kcoefr1(jrbasdeb:jrbasfin), &
     &            kmatUrr(jrmat,jrmatdeb:jrmatfin))
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilroa','algocoef_u_sf')
 1001 CALL printerror2(0,1001,3,'utilroa','algocoef_u_sf')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algomatc_u(klambda,kmatUrr,kforgfact)
!---------------------------------------------------------------------
!  
!  Purpose : Compute matrix transforming Sf into Sa
!  -------
!  Method : L = U inv[sqrt(rho I + Lambda)] U^
!  ------
!  Input :  klambda   : Lambda
!  -----    kmatUrr   : U
!           kforgfact : rho
!  Output : kmatUrr   : L
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(in) :: klambda
      BIGREAL, dimension(:,:), intent(inout) :: kmatUrr
      BIGREAL, intent(in) :: kforgfact
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok, jprsize
      INTEGER :: jrbasdeb, jrbasfin, jrbas, jrbas1, jrbas2
      INTEGER :: jrmatdeb, jrmatfin, jrmat, jrmat1, jrmat2
      BIGREAL, dimension(:), allocatable :: kcoefr1
      BIGREAL, dimension(:,:), allocatable :: matrr
!----------------------------------------------------------------------
!
      jprsize=size(klambda,1)
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
! --- allocation kcoefr1
      allocate ( kcoefr1(1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      kcoefr1(:) = FREAL(0.0)
! --- allocation matrr
      allocate ( matrr(1:jprsize,1:jprsize) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      matrr(:,:) = FREAL(0.0)
!
! -0.- Initialisation :
! ---------------------
!
      jrbasdeb=2
      jrbasfin=jprend
      jrmatdeb=1
      jrmatfin=jprend-1
!
! -1.- Compute diagonal reduction matrix
! --------------------------------------
!
      kcoefr1(jrbasdeb:jrbasfin) = SQRT( FREAL(1.0) / &
     &             ( kforgfact + klambda(jrbasdeb:jrbasfin) ) )
!
! -2.- Compute Sf to Sa transformation matrix
! -------------------------------------------
!
! --- L = U [ rho + lambda ]^.5
!     DO jrbas1=jrbasdeb,jrbasfin
!           jrmat1 = jrbas1 - 1
!           matrr(jrmat1,jrmatdeb:jrmatfin) = 
!    $           kcoefr1(jrbasdeb:jrbasfin)*
!    $           kmatUrr(jrmat1,jrmatdeb:jrmatfin)
!     ENDDO
! --- L = U  [ rho + lambda ]^.5  U^
      DO jrbas1=jrbasdeb,jrbasfin
         DO jrbas2=jrbasdeb,jrbasfin
            jrmat1 = jrbas1 - 1
            jrmat2 = jrbas2 - 1
            matrr(jrmat1,jrmat2) =  &
     &           DOT_PRODUCT(kcoefr1(jrbasdeb:jrbasfin) * &
     &               kmatUrr(jrmat1,jrmatdeb:jrmatfin), &
     &               kmatUrr(jrmat2,jrmatdeb:jrmatfin))
         ENDDO
      ENDDO
!
      kmatUrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin) = &
     &          matrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilroa','algomatc_u')
 1001 CALL printerror2(0,1001,3,'utilroa','algomatc_u')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocoeflimit_u(kcoefr,kmatUrr,kcoefrmax,kforgfact)
!---------------------------------------------------------------------
!
!  Purpose : Saturate the ROA coefficients (U method)
!  -------
!  Method : c = min ( c , cmax )
!  ------
!  Input :  kcoefr : c
!  -----    kmatUrr : U
!           kcoefrmax : cmax
!           kforgfact : rho
!
!  Output : kcoefr : c
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
       BIGREAL, dimension(:), intent(inout) :: kcoefr
       BIGREAL, dimension(:,:), intent(in) :: kmatUrr
       BIGREAL, intent(in) :: kcoefrmax,kforgfact
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL :: absval,maxval
      INTEGER :: jprsize, allocok
      INTEGER :: jrbas, jrbasdeb, jrbasfin
      INTEGER :: jrmat, jrmatdeb, jrmatfin
      BIGREAL, dimension(:), allocatable :: ucoefr
!----------------------------------------------------------------------
!
      jprsize=size(kcoefr,1)
      IF (jprsize.NE.jprend) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
! --- allocation ucoefr
      allocate ( ucoefr(1:jprsize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      ucoefr(:) = FREAL(0.0)
!
      jrbasdeb=2
      jrbasfin=jprsize
      jrmatdeb=1
      jrmatfin=jprsize-1
!
! -1.- Transform ROA coefficients to analysis eigenbasis
! ------------------------------------------------------
!
      DO jrbas=jrbasdeb,jrbasfin
         jrmat = jrbas - 1
         ucoefr(jrbas) = DOT_PRODUCT(kcoefr(jrbasdeb:jrbasfin), &
     &        kmatUrr(jrmatdeb:jrmatfin,jrmat))
      ENDDO
!
! -2.- Limit ROA coefficients
! ---------------------------
!
      maxval = ABS(kcoefrmax)/FREAL(SQRT(kforgfact))
      DO jrbas=jrbasdeb,jrbasfin
         absval = ABS(ucoefr(jrbas))
         absval = MIN(absval,maxval)
         ucoefr(jrbas) = SIGN(absval,ucoefr(jrbas))
      ENDDO
!
! -5.- Transform the ROA coefficients back to the original basis
! --------------------------------------------------------------
!
      DO jrbas=jrbasdeb,jrbasfin
         jrmat = jrbas - 1
         kcoefr(jrbas) = DOT_PRODUCT(ucoefr(jrbasdeb:jrbasfin), &
     &        kmatUrr(jrmat,jrmatdeb:jrmatfin))
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilroa','algocoeflimit_u')
 1001 CALL printerror2(0,1001,3,'utilroa','algocoeflimit_u')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocalcGi(kbasesr,kvectssqrdiagRi, &
     &                      kvectspart,kmatGrri)
!---------------------------------------------------------------------
!
!  Purpose : Compute Gamma matrices for obs segments
!  -------
!  Method : Gamma_i  = ( trsp(H_i Sf) inv(R_i) (H_i Sf) )
!  ------
!  Input :  kbasesr         : HSf
!  -----    kvectssqrdiagRi : inv(R)
!           kvectspart      : partition of the obs vector
!
!  Output : kmatGrri        : Gamma_i
!  ------
!---------------------------------------------------------------------
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
!----------------------------------------------------------------------
! header declarations
! ===================
       BIGREAL, dimension(:,:), intent(in) :: kbasesr
       BIGREAL, dimension(:), intent(in) :: kvectssqrdiagRi
       BIGREAL, dimension(:), intent(in) :: kvectspart
       BIGREAL, dimension(:,:,:), intent(out) :: kmatGrri
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jpssize,jprsize,jpisize,jiobs,jpu,js
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      INTEGER :: nvpnull
      BIGREAL, dimension(:,:), allocatable :: matsr
      INTEGER, dimension(:), allocatable :: ptu, partidx
!----------------------------------------------------------------------
!
      jpssize=size(kbasesr,1)
      jprsize=size(kbasesr,2)
      jpisize=size(kmatGrri,3)
      IF (jpssize.NE.size(kvectssqrdiagRi,1)) GOTO 1000
      IF (jpssize.NE.size(kvectspart,1)) GOTO 1000
      IF (jprsize.NE.size(kmatGrri,1)) GOTO 1000
      IF (jprsize.NE.size(kmatGrri,2)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
!
! --- allocation matsr
      allocate ( matsr(1:jpssize,1:jprsize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      matsr(:,:) = FREAL(0.0)
! --- allocation ptu
      allocate ( ptu(1:jpssize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
! --- allocation partidx
      allocate ( partidx(1:jpssize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      partidx(:) = NINT(kvectspart(:))
!
      IF (MAXVAL(partidx).GT.jpisize) GOTO 1000
      IF (MINVAL(partidx).LT.1) GOTO 1000
!
! -0.- Initialisation :
! ---------------------
!
      jrbasdeb=2
      jrbasfin=jprsize
      jrmatdeb=1
      jrmatfin=jprsize-1
!
! -1.- Compute the ROA kernel matrix (trsp((HSf)inv(R)(HSf))
! -----------------------------------------------------------
!
      DO jrbas=jrbasdeb,jrbasfin
        jrmat = jrbas - 1
        matsr(:,jrmat) = kbasesr(:,jrbas) * kvectssqrdiagRi(:)
      ENDDO
!
      DO jiobs=1,jpisize
!
        jpu = COUNT(partidx(1:jpssize).EQ.jiobs)
        ptu(:jpu) = PACK( (/ (js, js=1,jpssize) /) , &
     &              partidx(:) .EQ. jiobs )
!
        DO jrmat2=jrmatdeb,jrmatfin
          DO jrmat1=jrmat2,jrmatfin
             kmatGrri(jrmat1,jrmat2,jiobs) =  &
     &              DOT_PRODUCT(matsr(ptu(:jpu),jrmat1), &
     &                         matsr(ptu(:jpu),jrmat2))
          ENDDO
          kmatGrri(jrmat2,jrmat2:jrmatfin,jiobs) = &
     &         kmatGrri(jrmat2:jrmatfin,jrmat2,jiobs)
        ENDDO
!
      ENDDO
!
! --- deallocate arrays
      IF (allocated(matsr)) deallocate(matsr)
      IF (allocated(partidx)) deallocate(partidx)
      IF (allocated(ptu)) deallocate(ptu)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'utilroa','algocalGi')
 1001 CALL printerror2(0,1001,3,'utilroa','algocalGi')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoker_Gi(kmatGrri,kmatUrr,klambda,beta)
!---------------------------------------------------------------------
!
!  Purpose : SVD of the sum of Gamma matrices for obs segments
!  -------
!  Method : (U,lambda)  = SVD ( SUM_i Gamma_i / beta_i )
!  ------
!  Input :  kmatGrri   : Gamma_i
!  -----    beta       : beta_i
!
!  Output : kmatUrr    : U
!  ------   klambda    : Lambda
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
       BIGREAL, dimension(:,:,:), intent(in) :: kmatGrri
       BIGREAL, dimension(:,:), intent(out) :: kmatUrr
       BIGREAL, dimension(:), intent(out) :: klambda
       BIGREAL, dimension(:), optional, intent(in) :: beta
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jprsize,jpisize,jiobs
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      INTEGER :: nvpnull
      BIGREAL, dimension(:,:), allocatable :: matrr
!----------------------------------------------------------------------
!
      jprsize=size(kmatGrri,1)
      jpisize=size(kmatGrri,3)
      IF (jprsize.NE.size(kmatGrri,2)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
      IF (jprsize.NE.size(klambda,1)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
      IF (present(beta)) THEN
        IF (jpisize.NE.size(beta,1)) GOTO 1000
      ENDIF
!
! --- allocation matrr
      allocate ( matrr(1:jprsize,1:jprsize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
!
      jrbasdeb=2
      jrbasfin=jprsize
      jrmatdeb=1
      jrmatfin=jprsize-1
!
! -1.- Compute the sum of the Gamma_i matrices
! --------------------------------------------
!
      matrr(:,:) = FREAL(0.0) 
      IF (present(beta)) THEN
        DO jiobs=1,jpisize
          matrr(:,:) = matrr(:,:) + kmatGrri(:,:,jiobs) / beta(jiobs)
        ENDDO
      ELSE
        DO jiobs=1,jpisize
          matrr(:,:) = matrr(:,:) + kmatGrri(:,:,jiobs)
        ENDDO
      ENDIF
!
! -2.- Compute SVD decomposition of C
! -----------------------------------
!
      nvpnull=0
      CALL mkrs (matrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin), &
     &     kmatUrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin), &
     &     klambda(jrbasdeb:jrbasfin),nvpnull)
!
! --- deallocate arrays
      IF (allocated(matrr)) deallocate(matrr)
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'utilroa','algoker_Gi')
 1001 CALL printerror2(0,1001,3,'utilroa','algoker_Gi')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocalcDi(kbasesr,kvectsinnov,kvectssqrdiagRi, &
     &                      kvectspart,kdeltari)
!---------------------------------------------------------------------
!
!  Purpose : Compute delta vectors for obs segments
!  -------
!  Method : delta_i = trsp(H_i Sf) inv(R_i) (y_i-H_i xf)
!  ------
!  Input :  kbasesr         : HSf
!  -----    kvectsinnov     : y-Hxf
!           kvectssqrdiagRi : inv(R)
!           kvectspart      : partition of the obs vector
!  Output : kdeltari        : delta_i
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
       BIGREAL, dimension(:,:), intent(in) :: kbasesr
       BIGREAL, dimension(:), intent(in) :: kvectsinnov
       BIGREAL, dimension(:), intent(in) :: kvectssqrdiagRi
       BIGREAL, dimension(:), intent(in) :: kvectspart
       BIGREAL, dimension(:,:), intent(out) :: kdeltari
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jpssize,jprsize,jpisize,jiobs,jpu,js
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      BIGREAL, dimension(:), allocatable :: vects
      INTEGER, dimension(:), allocatable :: ptu, partidx
!----------------------------------------------------------------------
!
      jpssize=size(kbasesr,1)
      jprsize=size(kbasesr,2)
      jpisize=size(kdeltari,2)
      IF (jpssize.NE.size(kvectsinnov,1)) GOTO 1000
      IF (jpssize.NE.size(kvectspart,1)) GOTO 1000
      IF (jpssize.NE.size(kvectssqrdiagRi,1)) GOTO 1000
      IF (jprsize.NE.size(kdeltari,1)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
! --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
! --- allocation ptu
      allocate ( ptu(1:jpssize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
! --- allocation partidx
      allocate ( partidx(1:jpssize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      partidx(:) = NINT(kvectspart(:))
!
      IF (MAXVAL(partidx).GT.jpisize) GOTO 1000
      IF (MINVAL(partidx).LT.1) GOTO 1000
!
! -0.- Initialisation
! -------------------
!
      jrbasdeb=2
      jrbasfin=jprend
      jrmatdeb=1
      jrmatfin=jprend-1
!
! -1.- Compute innovation vector and weight by observation error
! --------------------------------------------------------------
!
      vects(:) = kvectsinnov(:) * &
     &           kvectssqrdiagRi(:) * kvectssqrdiagRi(:)
!
! -2.- Compute innovation in reduced space for all segments
! ---------------------------------------------------------
!
      DO jiobs=1,jpisize
!
        jpu = COUNT(partidx(1:jpssize).EQ.jiobs)
        ptu(:jpu) = PACK( (/ (js, js=1,jpssize) /) , &
     &              partidx(:) .EQ. jiobs )
!
        DO jrbas=jrbasdeb,jrbasfin
          kdeltari(jrbas,jiobs)=DOT_PRODUCT(kbasesr(ptu(:jpu),jrbas), &
     &                                      vects(ptu(:jpu)))
        ENDDO
!
      ENDDO
!
! --- deallocate arrays
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(partidx)) deallocate(partidx)
      IF (allocated(ptu)) deallocate(ptu)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilroa','algocoef_u')
 1001 CALL printerror2(0,1001,3,'utilroa','algocoef_u')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocalciRi(kvectsinnov,kvectssqrdiagRi, &
     &                       kiRi,ky,kvectspart)
!---------------------------------------------------------------------
!
!  Purpose : Compute iT R-1 i for each segment of the obs.
!  -------
!  Method : trsp(i) inv(R_i) (i)
!  ------
!  Input :  kvectsinnov     : y-Hxf
!  -----    kvectssqrdiagRi : inv(R)
!           kvectspart      : partition of the obs vector
!  Output : kiRi            : iT R-1 i
!  ------   ky              : number of observations
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpyend,jpoend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
       BIGREAL, dimension(:), intent(in) :: kvectsinnov
       BIGREAL, dimension(:), intent(in) :: kvectssqrdiagRi
       BIGREAL, dimension(:), intent(out) :: kiRi
       BIGREAL, dimension(:), intent(out) :: ky
       BIGREAL, dimension(:), optional, intent(in) :: kvectspart
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jpssize,jpisize,jiobs,jpu,js
      BIGREAL, dimension(:), allocatable :: vects
      INTEGER, dimension(:), allocatable :: ptu, partidx
!----------------------------------------------------------------------
!
      jpssize=size(kvectsinnov,1)
      jpisize=size(kiRi,1)
      IF (jpssize.NE.size(kvectssqrdiagRi,1)) GOTO 1000
      IF ((jpisize.GT.1).AND.(.NOT.present(kvectspart))) GOTO 1000
! --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)

      IF (jpisize.GT.1) THEN
        IF (jpssize.NE.size(kvectspart,1)) GOTO 1000
! --- allocation ptu
        allocate ( ptu(1:jpssize) , stat=allocok )
        IF (allocok.GT.0) GOTO 1001
! --- allocation partidx
        allocate ( partidx(1:jpssize) , stat=allocok )
        IF (allocok.GT.0) GOTO 1001
        partidx(:) = NINT(kvectspart(:))
!
        IF (MAXVAL(partidx).GT.jpisize) GOTO 1000
        IF (MINVAL(partidx).LT.1) GOTO 1000
      ENDIF
!
! -1.- Compute innovation vector and weight by observation error
! --------------------------------------------------------------
!
      vects(:) = kvectsinnov(:) * kvectssqrdiagRi(:)
!
! -2.- Compute innovation in reduced space for all segments
! ---------------------------------------------------------
!
      IF (jpisize.GT.1) THEN

        DO jiobs=1,jpisize
!
          jpu = COUNT(partidx(1:jpssize).EQ.jiobs)
          ptu(:jpu) = PACK( (/ (js, js=1,jpssize) /) , &
     &                partidx(:) .EQ. jiobs )
!
          kiRi(jiobs)=DOT_PRODUCT(vects(ptu(:jpu)), &
     &                            vects(ptu(:jpu)))
          ky(jiobs)=FREAL(jpu)
!
        ENDDO
!
      ELSE
!
        kiRi(1)=DOT_PRODUCT(vects(:),vects(:))
        ky(1)=FREAL(jpssize)
!
      ENDIF
!
! --- deallocate arrays
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(partidx)) deallocate(partidx)
      IF (allocated(ptu)) deallocate(ptu)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilroa','algocalciRi')
 1001 CALL printerror2(0,1001,3,'utilroa','algocalciRi')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algosum_Di(kdeltari,kdelta,beta)
!---------------------------------------------------------------------
!
!  Purpose : Sum up the delta vectors for obs segments
!  -------
!  Method : delta  =  SUM_i delta_i / beta_i
!  ------
!  Input :  kdeltari   : delta_i
!  -----    beta       : beta_i
!
!  Output : kdelta     : delta
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
       BIGREAL, dimension(:,:), intent(in) :: kdeltari
       BIGREAL, dimension(:), intent(out) :: kdelta
       BIGREAL, dimension(:), optional, intent(in) :: beta
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jprsize,jpisize,jiobs
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      INTEGER :: nvpnull
      BIGREAL, dimension(:,:), allocatable :: matrr
!----------------------------------------------------------------------
!
      jprsize=size(kdeltari,1)
      jpisize=size(kdeltari,2)
      IF (jprsize.NE.size(kdelta,1)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
      IF (present(beta)) THEN
        IF (jpisize.NE.size(beta,1)) GOTO 1000
      ENDIF
!
      kdelta(:) = FREAL(0.0)
      IF (present(beta)) THEN
        DO jiobs=1,jpisize
          kdelta(:) = kdelta(:) + kdeltari(:,jiobs) / beta(jiobs)
        ENDDO
      ELSE
        DO jiobs=1,jpisize
          kdelta(:) = kdelta(:) + kdeltari(:,jiobs)
        ENDDO
      ENDIF
!
      RETURN
!
! --- error management section
!
 1000 CALL printerror2(0,1000,1,'utilroa','algosum_Di')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkyorms (vectsrms,kflagxyo) 
!---------------------------------------------------------------------
!
!  Purpose : Fill observation error standard deviation vector
!  -------   using parameters in SESAM configuration
!
!  Method : Loop on variables and affect the corresponding
!  ------   error value to the corresponding segment of the Vy vector
!
!  Input :  kflagxyo : vector type (2=Vy, 3=Vo)
!  -----
!  Output : vectsrms : observation error standard deviation vector
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(out) :: vectsrms
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: js,jsdeb,jsfin
      INTEGER :: jdta,inddta,jobs,indobs,inddbs
      BIGREAL :: valeur
!----------------------------------------------------------------------
!
      SELECT CASE(kflagxyo)
      CASE(1)
         GOTO 1000
      CASE(2)
! Vy vector
         jsdeb=1
         DO jdta=1,dtaend
            inddta=dta_ord(jdta)
            jsfin=jsdeb+dta_nbr(inddta)-1
            IF (lmoyect) THEN
               valeur=ABS(dta_rms(inddta))/dta_ect(inddta)
            ELSE
               valeur=ABS(dta_rms(inddta))
            ENDIF
            vectsrms(jsdeb:jsfin)=valeur
            IF (nprint.GE.2) THEN
               WRITE(numout,*) 'dtanam=',dta_nam(inddta) &
     &              (1:lenv(dta_nam(inddta)))
               WRITE(numout,*) 'dtarms=',dta_rms(inddta)
            ENDIF
            jsdeb=jsfin+1
         ENDDO
      CASE(3)
! Vo vector
         jsdeb=1
         DO jobs = 1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            jsfin=jsdeb+obs_nbr(indobs,inddbs)-1
            IF (lmoyect) THEN
               valeur=ABS(obs_rms(indobs,inddbs))/obs_ect(indobs,inddbs)
            ELSE
               valeur=ABS(obs_rms(indobs,inddbs))
            ENDIF
            vectsrms(jsdeb:jsfin)=valeur
            IF (nprint.GE.2) THEN
               WRITE(numout,*) 'obsnam=',obs_nam(indobs,inddbs) &
     &              (1:lenv(obs_nam(indobs,inddbs)))
               WRITE(numout,*) 'obsrms=',obs_rms(indobs,inddbs)
            ENDIF
            jsdeb=jsfin+1
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilroa','mkyorms')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkrs (kcov,vctp_rr,valp_r,nvpnull) 
!---------------------------------------------------------------------
!
!  Purpose : Compute eigenvalues and eigenvectors
!  -------   of a symmetric real matrix
!
!  Method : Use EISPACK routine 'rs'
!  ------
!  Input :  kcov : input matrix
!  -----
!  Output : valp_r : eigenvalues
!  ------   vctp_rr : eigenvectors
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use utilmath1
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:,:), intent(in) :: kcov
      BIGREAL, dimension(:,:), intent(out) :: vctp_rr
      BIGREAL, dimension(:), intent(out) :: valp_r
      INTEGER, intent(in) :: nvpnull
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: jprsize,jr1,jrdeb,nvectindep,lrs
      BIGREAL, dimension(:), allocatable :: work1,work2
      INTEGER :: ierr,allocok
!----------------------------------------------------------------------
!
! -0.- Initialisation
! -------------------
!
      jprsize = size(kcov,1)
      IF (jprsize .NE.size(kcov,2)) GOTO 1000
      IF (jprsize .NE.size(vctp_rr,1)) GOTO 1000
      IF (jprsize .NE.size(vctp_rr,2)) GOTO 1000
      IF (jprsize .NE.size(valp_r,1)) GOTO 1000
!
      allocate(work1(jprsize),stat=allocok)
      IF (allocok.GT.0) GOTO 1001
      allocate(work2(jprsize),stat=allocok)
      IF (allocok.GT.0) GOTO 1001
!
      valp_r(:) = FREAL(0.0)
      vctp_rr(:,:) = FREAL(0.0)
      nvectindep=jprsize-nvpnull
      jrdeb=nvpnull+1
!
! -1.- Computation of the eigenvalues and eigenvectors of the kcov matrix
! -----------------------------------------------------------------------
!
! --- EISPACK : routine package to solve eigenvalue problems
! --- rs      : eispack routine computing eigenvalues and eigenvectors
! ---           of a symmetric real matrix (see utilmath1.F90)
! ---
! --- Usage of subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
! ---
! --- input
! ---
! ---       nm  : must be set to the row dimension of the two-dimensional
! ---             array parameters as declared in the calling program
! ---             dimension statement.
! ---
! ---       n   : is the order of the matrix  a.
! ---
! ---       a   : contains the real symmetric matrix.
! ---
! ---       matz: is an integer variable set equal to zero if
! ---             only eigenvalues are desired.  otherwise it is set to
! ---             any non-zero integer for both eigenvalues and eigenvectors.
! ---
! --- output
! ---
! ---       w   : contains the eigenvalues in ascending order.
! ---
! ---       z   : contains the eigenvectors if matz is not zero.
! ---
! ---       ierr: is an integer output variable set equal to an error
! ---             completion code described in the documentation for tqlrat
! ---             and tql2.  the normal completion code is zero.
! ---
! ---       fv1,fv2 :  are temporary storage arrays.
! ---
      lrs = 1
      CALL rs(jprsize,jprsize,kcov,valp_r,lrs,vctp_rr,work1,work2,ierr)
!
! -2.- Reorganisation of valp_r and vctp_rr arrays
! ------------------------------------------------
! Eigenvalues are sorted in descending order from jrdeb=2 to jprsize
! The first value (jrdeb-1=1) correspond to the smallest eigenvalue (equal to 0)
!
      work1(:)=valp_r(:)
      valp_r(jrdeb:jprsize)=work1(jprsize:jrdeb:-1)
!
      DO jr1=1,jprsize
         work1(:)=vctp_rr(jr1,:)
         vctp_rr(jr1,jrdeb:jprsize)=work1(jprsize:jrdeb:-1)
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilroa','mkrs')
 1001 CALL printerror2(0,1001,3,'utilroa','mkrs')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkcoeftodta (kfile,kvectz,kparty)
!---------------------------------------------------------------------
!
!  Purpose : Save local z-vector (constant for
!  -------   every subsystems) in dta file
!
!  Method : Fill variables in Vy vector, with the value
!  ------   of the subsystem to which they belong
!            
!  Input :  kfile  : Name of the output file
!  -----    kvectz : z-vector
!           kparty : partition of the Vy object into subsystems
!
!  Output : none
!  ------
!---------------------------------------------------------------------
! modules
! =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpyend, jpxend, jprend, jpz
      use hioxyo
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfile
      BIGREAL, dimension(:), intent(in) :: kvectz
      BIGREAL, dimension(:), intent(in) :: kparty
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: allocok,jpysize,jpzsize
      BIGREAL, dimension(:), allocatable :: diagy
!----------------------------------------------------------------------
!
      jpzsize=size(kvectz,1)
      jpysize=size(kparty,1)
      IF (jpysize.NE.jpyend) GOTO 1000
!
! --- allocation diagy
      allocate ( diagy(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      diagy(:) = FREAL(0.0)
!---------------------------------------------------------------------
!
! -1.- Fill Vy vector with the local quantity
! ---------------------------------------------
!
      diagy(1:jpysize) = kvectz(NINT(kparty(1:jpysize)))
!
! -2.- Save Vy vector in dta file
! -------------------------------
!
      CALL writedta(kfile,diagy(:))
!
! --- deallocation
      IF (allocated(diagy)) deallocate(diagy)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algoutilroa','mkcoeftodta')
 1001 CALL printerror2(0,1001,3,'algoutilroa','mkcoeftodta')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilroa
