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
C ---                  UTILROA.F                               ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 98-06 (C.E. Testut)                        ---
C --- modification : 99-05 (C.E. Testut)                        ---
C --- modification : 99-12 (J.M. Brankart)                      ---
C --- modification : 00-09 (J.M. Brankart)                      ---
C --- modification : 07-11 (J.M. Brankart)                      ---
C --- modification : 08-10 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- Set of routines used in the Reduced Order Analysis (ROA) algorithm
C --- 
C --- SUBROUTINE  algoker_u : Compute SVD decomposition of the ROA kernel matrix
C --- SUBROUTINE  algocoef_u : Compute innovation in the reduced eigenspace
C --- SUBROUTINE  algocoef_u_sf : Compute the ROA coefficients (U method)
C --- SUBROUTINE  algomatc_u : Compute matrix transforming Sf into Sa
C --- SUBROUTINE  algocoeflimit_u : Saturate the ROA coefficients (U method)
C --- SUBROUTINE  algocalcGi : Compute Gamma matrices for obs segments
C --- SUBROUTINE  algoker_Gi : SVD of the sum of Gamma matrices for obs segments
C --- SUBROUTINE  algocalcDi : Compute delta vectors for obs segments
C --- SUBROUTINE  algosum_Di : Sum up the delta vectors for obs segments
C --- SUBROUTINE  mkyorms : Fill observation error standard deviation vector
C ---                       using parameters from SESAM configuration
C --- SUBROUTINE  mkrs : Compute eigenvalues and eigenvectors
C ---                    of a symmetric real matrix
C --- SUBROUTINE  mkcoeftodta : Save local z-vector in dta file
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilroa
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC algoker_u, algocoef_u, algocoef_u_sf, algomatc_u
      PUBLIC algocoeflimit_u, algocalcGi, algoker_Gi, algocalcDi
      PUBLIC algosum_Di, mkyorms, mkrs, mkcoeftodta, algocalciRi

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoker_u(kbasesr,kvectssqrdiagRi,kmatUrr,klambda,
     $                     kbasesq,kmatBqr)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute SVD decomposition of the ROA kernel matrix
CCC  -------
CCC  Method : (U,lambda)  = SVD ( trsp(HSf) inv(R) (HSf) )
CCC  ------
CCC  Input :  kbasesr         : HSf
CCC  -----    kvectssqrdiagRi : inv(R)
CCC           kforgfact       : rho
CCC           kbasesq         : SQRT of obs. inv. correl. matrix
CCC
CCC  Output : kmatUrr         : U
CCC  ------   klambda         : Lambda
CCC           kmatBqr         : SQRT of obs. correl. matrix in reduced space
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
       BIGREAL, dimension(:,:), intent(in) :: kbasesr
       BIGREAL, dimension(:), intent(in) :: kvectssqrdiagRi
       BIGREAL, dimension(:,:), intent(out) :: kmatUrr
       BIGREAL, dimension(:), intent(out) :: klambda
       BIGREAL, dimension(:,:), optional, intent(in) :: kbasesq
       BIGREAL, dimension(:,:), optional, intent(out) :: kmatBqr
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpssize,jprsize,jpqsize
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      INTEGER :: jqbasdeb,jqbasfin,jqbas,jqbas1,jqbas2
      INTEGER :: jqmatdeb,jqmatfin,jqmat,jqmat1,jqmat2
      INTEGER :: nvpnull
      BIGREAL, dimension(:,:), allocatable :: matsr
      BIGREAL, dimension(:,:), allocatable :: matrr
CC----------------------------------------------------------------------
C
      jpssize=size(kbasesr,1)
      jprsize=size(kbasesr,2)
      IF (present(kbasesq)) THEN
         jpqsize=size(kbasesq,2)
      ENDIF
C
C --- allocation matsr
      allocate ( matsr(1:jpssize,1:jprsize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      matsr(:,:) = FREAL(0.0)
C --- allocation matrr
      allocate ( matrr(1:jprsize,1:jprsize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      matrr(:,:) = FREAL(0.0)
C
C -0.- Initialisation :
C ---------------------
C
      IF (jpssize.NE.size(kvectssqrdiagRi,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
C
      IF (present(kbasesq)) THEN
         IF (jpssize.NE.size(kbasesq,1)) GOTO 1000
         IF (.NOT.present(kmatBqr)) GOTO 1000
         IF (jpqsize.NE.size(kmatBqr,1)) GOTO 1000
         IF (jprsize.NE.size(kmatBqr,2)) GOTO 1000
      ENDIF
C
      jrbasdeb=2
      jrbasfin=jprsize
      jrmatdeb=1
      jrmatfin=jprsize-1
C
      IF (present(kbasesq)) THEN
         jqbasdeb=2
         jqbasfin=jpqsize
         jqmatdeb=1
         jqmatfin=jpqsize-1
      ENDIF
C
C -1.- Compute the ROA kernel matrix (trsp((HSf)inv(R)(HSf))
C -----------------------------------------------------------
C
      DO jrbas=jrbasdeb,jrbasfin
         jrmat = jrbas - 1
         matsr(:,jrmat) = kbasesr(:,jrbas) * kvectssqrdiagRi(:)
      ENDDO
C
      IF (.NOT.present(kbasesq)) THEN
         DO jrmat2=jrmatdeb,jrmatfin
            DO jrmat1=jrmat2,jrmatfin
               matrr(jrmat1,jrmat2) = 
     $              DOT_PRODUCT(matsr(:,jrmat1),matsr(:,jrmat2))
            ENDDO
            matrr(jrmat2,jrmat2:jrmatfin)=matrr(jrmat2:jrmatfin,jrmat2)
         ENDDO
      ELSE
         DO jrmat2=jrmatdeb,jrmatfin
            DO jqbas1=jqbasdeb,jqbasfin
               jqmat1 = jqbas1 - 1
               kmatBqr(jqmat1,jrmat2) = 
     $              DOT_PRODUCT(kbasesq(:,jqbas1),matsr(:,jrmat2))
            ENDDO
         ENDDO
         DO jrmat2=jrmatdeb,jrmatfin
            DO jrmat1=jrmat2,jrmatfin
               matrr(jrmat1,jrmat2) = 
     $              DOT_PRODUCT(kmatBqr(jqmatdeb:jqmatfin,jrmat1),
     $              kmatBqr(jqmatdeb:jqmatfin,jrmat2))
            ENDDO
            matrr(jrmat2,jrmat2:jrmatfin)=matrr(jrmat2:jrmatfin,jrmat2)
         ENDDO
      ENDIF
C
C -2.- Compute SVD decomposition of C
C -----------------------------------
C
      nvpnull=0
      CALL mkrs (matrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin),
     $     kmatUrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin),
     $     klambda(jrbasdeb:jrbasfin),nvpnull)
C
C --- deallocate arrays
      IF (allocated(matsr)) deallocate(matsr)
      IF (allocated(matrr)) deallocate(matrr)
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'utilroa','algoker_u')
 1001 CALL printerror2(0,1001,3,'utilroa','algoker_u')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocoef_u(kbasesr,kvectsinnov,kvectssqrdiagRi,
     $                      kmatUrr,kdelta,kbasesq,kmatBqr)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute innovation in the reduced eigenspace
CCC  -------
CCC  Method : delta = U^ trsp(HSf) inv(R) (y-Hxf)
CCC  ------
CCC  Input :  kbasesr         : HSf
CCC  -----    kvectsinnov     : y-Hxf
CCC           kvectssqrdiagRi : inv(R)
CCC           kmatUrr         : U
CCC           kbasesq         : SQRT of inv. obs. correl. matrix
CCC  Output : kdelta          : delta
CCC  ------   kmatBqr         : SQRT of obs. correl. matrix in reduced space
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
       BIGREAL, dimension(:,:), intent(in) :: kbasesr
       BIGREAL, dimension(:), intent(in) :: kvectsinnov
       BIGREAL, dimension(:), intent(in) :: kvectssqrdiagRi
       BIGREAL, dimension(:,:), intent(in) :: kmatUrr
       BIGREAL, dimension(:), intent(out) :: kdelta
       BIGREAL, dimension(:,:), optional, intent(in) :: kbasesq
       BIGREAL, dimension(:,:), optional, intent(in) :: kmatBqr
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpssize,jprsize,jpqsize
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      INTEGER :: jqbasdeb,jqbasfin,jqbas,jqbas1,jqbas2
      INTEGER :: jqmatdeb,jqmatfin,jqmat,jqmat1,jqmat2
      BIGREAL, dimension(:), allocatable :: kcoefr1
      BIGREAL, dimension(:), allocatable :: kcoefq1
      BIGREAL, dimension(:), allocatable :: vects
CC----------------------------------------------------------------------
C
      jpssize=size(kbasesr,1)
      jprsize=size(kbasesr,2)
      IF (present(kbasesq)) THEN
         jpqsize=size(kbasesq,2)
      ENDIF
C --- allocation kcoefr1
      allocate ( kcoefr1(1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      kcoefr1(:) = FREAL(0.0)
C --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
C --- allocation kcoefq1
      IF (present(kbasesq)) THEN
         allocate ( kcoefq1(1:jpqsize), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
         kcoefq1(:) = FREAL(0.0)
      ENDIF
C
C -0.- Initialisation
C -------------------
C
      IF (jpssize.NE.size(kvectsinnov,1)) GOTO 1000
      IF (jpssize.NE.size(kvectssqrdiagRi,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
      IF (jprsize.NE.size(kdelta,1)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
C
      IF (present(kbasesq)) THEN
         IF (jpssize.NE.size(kbasesq,1)) GOTO 1000
         IF (.NOT.present(kmatBqr)) GOTO 1000
         IF (jpqsize.NE.size(kmatBqr,1)) GOTO 1000
         IF (jprsize.NE.size(kmatBqr,2)) GOTO 1000
      ENDIF
C
      jrbasdeb=2
      jrbasfin=jprend
      jrmatdeb=1
      jrmatfin=jprend-1
C
      IF (present(kbasesq)) THEN
         jqbasdeb=2
         jqbasfin=jpqsize
         jqmatdeb=1
         jqmatfin=jpqsize-1
      ENDIF
C
C -1.- Compute innovation vector and weight by observation error
C --------------------------------------------------------------
C
      IF (.NOT.present(kbasesq)) THEN
         vects(:) = kvectsinnov(:) *
     $           kvectssqrdiagRi(:) * kvectssqrdiagRi(:)
      ELSE
         vects(:) = kvectsinnov(:) * kvectssqrdiagRi(:)
      ENDIF
C
C -2.- Compute innovation in reduced space
C ----------------------------------------
C
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
            kcoefr1(jrbas)=DOT_PRODUCT(kmatBqr(jqmatdeb:jqmatfin,jrmat),
     $                                kcoefq1(jqbasdeb:jqbasfin))
         ENDDO
      ENDIF
C
C -3.- Transform to analysis eigenpsace
C -------------------------------------
C
      DO jrbas=jrbasdeb,jrbasfin
         jrmat = jrbas - 1
         kdelta(jrbas)=DOT_PRODUCT(kmatUrr(jrmatdeb:jrmatfin,jrmat),
     $             kcoefr1(jrbasdeb:jrbasfin))
      ENDDO
C
C --- deallocate arrays
      IF (allocated(kcoefr1)) deallocate (kcoefr1)
      IF (allocated(kcoefq1)) deallocate (kcoefq1)
      IF (allocated(vects)) deallocate(vects)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilroa','algocoef_u')
 1001 CALL printerror2(0,1001,3,'utilroa','algocoef_u')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocoef_u_sf(kcoefr,klambda,kdelta,kmatUrr,kforgfact)
CCC---------------------------------------------------------------------
CCC   
CCC  Purpose : Compute the ROA coefficients (U method)
CCC  -------
CCC  Method : c = U inv[rho I + Lambda] delta
CCC  ------
CCC  Input :  klambda         : Lambda
CCC  -----    kdelta          : delta
CCC           kmatUrr         : U
CCC           kforgfact       : rho
CCC  Output : kcoefr          : ROA coefficients
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(out) :: kcoefr
      BIGREAL, dimension(:), intent(in) :: klambda
      BIGREAL, dimension(:), intent(in) :: kdelta
      BIGREAL, dimension(:,:), intent(in) :: kmatUrr
      BIGREAL, intent(in) :: kforgfact
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok, jprsize
      INTEGER :: jrbasdeb, jrbasfin, jrbas
      INTEGER :: jrmatdeb, jrmatfin, jrmat
      BIGREAL, dimension(:), allocatable :: kcoefr1
CC----------------------------------------------------------------------
C
      jprsize=size(kcoefr,1)
      IF (jprsize.NE.size(klambda,1)) GOTO 1000
      IF (jprsize.NE.size(kdelta,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
C --- allocation kcoefr1
      allocate ( kcoefr1(1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      kcoefr1(:) = FREAL(0.0)
C
C -0.- Initialisation :
C ---------------------
C
      jrbasdeb=2
      jrbasfin=jprend
      jrmatdeb=1
      jrmatfin=jprend-1
C
C -1.- Compute the ROA coefficients
C ----------------------------------
C
      kcoefr1(jrbasdeb:jrbasfin) = kdelta(jrbasdeb:jrbasfin) /
     $             ( kforgfact + klambda(jrbasdeb:jrbasfin) )
C
C -2.- Transform the ROA coefficient back to the original basis
C --------------------------------------------------------------
C
      DO jrbas=jrbasdeb,jrbasfin
         jrmat = jrbas - 1
         kcoefr(jrbas) = DOT_PRODUCT(kcoefr1(jrbasdeb:jrbasfin),
     $            kmatUrr(jrmat,jrmatdeb:jrmatfin))
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilroa','algocoef_u_sf')
 1001 CALL printerror2(0,1001,3,'utilroa','algocoef_u_sf')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algomatc_u(klambda,kmatUrr,kforgfact)
CCC---------------------------------------------------------------------
CCC  
CCC  Purpose : Compute matrix transforming Sf into Sa
CCC  -------
CCC  Method : L = U inv[sqrt(rho I + Lambda)] U^
CCC  ------
CCC  Input :  klambda   : Lambda
CCC  -----    kmatUrr   : U
CCC           kforgfact : rho
CCC  Output : kmatUrr   : L
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(in) :: klambda
      BIGREAL, dimension(:,:), intent(inout) :: kmatUrr
      BIGREAL, intent(in) :: kforgfact
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok, jprsize
      INTEGER :: jrbasdeb, jrbasfin, jrbas, jrbas1, jrbas2
      INTEGER :: jrmatdeb, jrmatfin, jrmat, jrmat1, jrmat2
      BIGREAL, dimension(:), allocatable :: kcoefr1
      BIGREAL, dimension(:,:), allocatable :: matrr
CC----------------------------------------------------------------------
C
      jprsize=size(klambda,1)
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
C --- allocation kcoefr1
      allocate ( kcoefr1(1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      kcoefr1(:) = FREAL(0.0)
C --- allocation matrr
      allocate ( matrr(1:jprsize,1:jprsize) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      matrr(:,:) = FREAL(0.0)
C
C -0.- Initialisation :
C ---------------------
C
      jrbasdeb=2
      jrbasfin=jprend
      jrmatdeb=1
      jrmatfin=jprend-1
C
C -1.- Compute diagonal reduction matrix
C --------------------------------------
C
      kcoefr1(jrbasdeb:jrbasfin) = SQRT( FREAL(1.0) /
     $             ( kforgfact + klambda(jrbasdeb:jrbasfin) ) )
C
C -2.- Compute Sf to Sa transformation matrix
C -------------------------------------------
C
C --- L = U [ rho + lambda ]^.5
C     DO jrbas1=jrbasdeb,jrbasfin
C           jrmat1 = jrbas1 - 1
C           matrr(jrmat1,jrmatdeb:jrmatfin) = 
C    $           kcoefr1(jrbasdeb:jrbasfin)*
C    $           kmatUrr(jrmat1,jrmatdeb:jrmatfin)
C     ENDDO
C --- L = U  [ rho + lambda ]^.5  U^
      DO jrbas1=jrbasdeb,jrbasfin
         DO jrbas2=jrbasdeb,jrbasfin
            jrmat1 = jrbas1 - 1
            jrmat2 = jrbas2 - 1
            matrr(jrmat1,jrmat2) = 
     $           DOT_PRODUCT(kcoefr1(jrbasdeb:jrbasfin) *
     $               kmatUrr(jrmat1,jrmatdeb:jrmatfin),
     $               kmatUrr(jrmat2,jrmatdeb:jrmatfin))
         ENDDO
      ENDDO
C
      kmatUrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin) =
     $          matrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilroa','algomatc_u')
 1001 CALL printerror2(0,1001,3,'utilroa','algomatc_u')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocoeflimit_u(kcoefr,kmatUrr,kcoefrmax,kforgfact)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Saturate the ROA coefficients (U method)
CCC  -------
CCC  Method : c = min ( c , cmax )
CCC  ------
CCC  Input :  kcoefr : c
CCC  -----    kmatUrr : U
CCC           kcoefrmax : cmax
CCC           kforgfact : rho
CCC
CCC  Output : kcoefr : c
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
       BIGREAL, dimension(:), intent(inout) :: kcoefr
       BIGREAL, dimension(:,:), intent(in) :: kmatUrr
       BIGREAL, intent(in) :: kcoefrmax,kforgfact
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL :: absval,maxval
      INTEGER :: jprsize, allocok
      INTEGER :: jrbas, jrbasdeb, jrbasfin
      INTEGER :: jrmat, jrmatdeb, jrmatfin
      BIGREAL, dimension(:), allocatable :: ucoefr
CC----------------------------------------------------------------------
C
      jprsize=size(kcoefr,1)
      IF (jprsize.NE.jprend) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,1)) GOTO 1000
      IF (jprsize.NE.size(kmatUrr,2)) GOTO 1000
C --- allocation ucoefr
      allocate ( ucoefr(1:jprsize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      ucoefr(:) = FREAL(0.0)
C
      jrbasdeb=2
      jrbasfin=jprsize
      jrmatdeb=1
      jrmatfin=jprsize-1
C
C -1.- Transform ROA coefficients to analysis eigenbasis
C ------------------------------------------------------
C
      DO jrbas=jrbasdeb,jrbasfin
         jrmat = jrbas - 1
         ucoefr(jrbas) = DOT_PRODUCT(kcoefr(jrbasdeb:jrbasfin),
     $        kmatUrr(jrmatdeb:jrmatfin,jrmat))
      ENDDO
C
C -2.- Limit ROA coefficients
C ---------------------------
C
      maxval = ABS(kcoefrmax)/FREAL(SQRT(kforgfact))
      DO jrbas=jrbasdeb,jrbasfin
         absval = ABS(ucoefr(jrbas))
         absval = MIN(absval,maxval)
         ucoefr(jrbas) = SIGN(absval,ucoefr(jrbas))
      ENDDO
C
C -5.- Transform the ROA coefficients back to the original basis
C --------------------------------------------------------------
C
      DO jrbas=jrbasdeb,jrbasfin
         jrmat = jrbas - 1
         kcoefr(jrbas) = DOT_PRODUCT(ucoefr(jrbasdeb:jrbasfin),
     $        kmatUrr(jrmat,jrmatdeb:jrmatfin))
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilroa','algocoeflimit_u')
 1001 CALL printerror2(0,1001,3,'utilroa','algocoeflimit_u')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocalcGi(kbasesr,kvectssqrdiagRi,
     $                      kvectspart,kmatGrri)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute Gamma matrices for obs segments
CCC  -------
CCC  Method : Gamma_i  = ( trsp(H_i Sf) inv(R_i) (H_i Sf) )
CCC  ------
CCC  Input :  kbasesr         : HSf
CCC  -----    kvectssqrdiagRi : inv(R)
CCC           kvectspart      : partition of the obs vector
CCC
CCC  Output : kmatGrri        : Gamma_i
CCC  ------
CCC---------------------------------------------------------------------
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
CC----------------------------------------------------------------------
CC header declarations
CC ===================
       BIGREAL, dimension(:,:), intent(in) :: kbasesr
       BIGREAL, dimension(:), intent(in) :: kvectssqrdiagRi
       BIGREAL, dimension(:), intent(in) :: kvectspart
       BIGREAL, dimension(:,:,:), intent(out) :: kmatGrri
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpssize,jprsize,jpisize,jiobs,jpu,js
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      INTEGER :: nvpnull
      BIGREAL, dimension(:,:), allocatable :: matsr
      INTEGER, dimension(:), allocatable :: ptu, partidx
CC----------------------------------------------------------------------
C
      jpssize=size(kbasesr,1)
      jprsize=size(kbasesr,2)
      jpisize=size(kmatGrri,3)
      IF (jpssize.NE.size(kvectssqrdiagRi,1)) GOTO 1000
      IF (jpssize.NE.size(kvectspart,1)) GOTO 1000
      IF (jprsize.NE.size(kmatGrri,1)) GOTO 1000
      IF (jprsize.NE.size(kmatGrri,2)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
C
C --- allocation matsr
      allocate ( matsr(1:jpssize,1:jprsize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      matsr(:,:) = FREAL(0.0)
C --- allocation ptu
      allocate ( ptu(1:jpssize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
C --- allocation partidx
      allocate ( partidx(1:jpssize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      partidx(:) = NINT(kvectspart(:))
C
      IF (MAXVAL(partidx).GT.jpisize) GOTO 1000
      IF (MINVAL(partidx).LT.1) GOTO 1000
C
C -0.- Initialisation :
C ---------------------
C
      jrbasdeb=2
      jrbasfin=jprsize
      jrmatdeb=1
      jrmatfin=jprsize-1
C
C -1.- Compute the ROA kernel matrix (trsp((HSf)inv(R)(HSf))
C -----------------------------------------------------------
C
      DO jrbas=jrbasdeb,jrbasfin
        jrmat = jrbas - 1
        matsr(:,jrmat) = kbasesr(:,jrbas) * kvectssqrdiagRi(:)
      ENDDO
C
      DO jiobs=1,jpisize
C
        jpu = COUNT(partidx(1:jpssize).EQ.jiobs)
        ptu(:jpu) = PACK( (/ (js, js=1,jpssize) /) ,
     $              partidx(:) .EQ. jiobs )
C
        DO jrmat2=jrmatdeb,jrmatfin
          DO jrmat1=jrmat2,jrmatfin
             kmatGrri(jrmat1,jrmat2,jiobs) = 
     $              DOT_PRODUCT(matsr(ptu(:jpu),jrmat1),
     $                         matsr(ptu(:jpu),jrmat2))
          ENDDO
          kmatGrri(jrmat2,jrmat2:jrmatfin,jiobs) =
     $         kmatGrri(jrmat2:jrmatfin,jrmat2,jiobs)
        ENDDO
C
      ENDDO
C
C --- deallocate arrays
      IF (allocated(matsr)) deallocate(matsr)
      IF (allocated(partidx)) deallocate(partidx)
      IF (allocated(ptu)) deallocate(ptu)
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'utilroa','algocalGi')
 1001 CALL printerror2(0,1001,3,'utilroa','algocalGi')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algoker_Gi(kmatGrri,kmatUrr,klambda,beta)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : SVD of the sum of Gamma matrices for obs segments
CCC  -------
CCC  Method : (U,lambda)  = SVD ( SUM_i Gamma_i / beta_i )
CCC  ------
CCC  Input :  kmatGrri   : Gamma_i
CCC  -----    beta       : beta_i
CCC
CCC  Output : kmatUrr    : U
CCC  ------   klambda    : Lambda
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
       BIGREAL, dimension(:,:,:), intent(in) :: kmatGrri
       BIGREAL, dimension(:,:), intent(out) :: kmatUrr
       BIGREAL, dimension(:), intent(out) :: klambda
       BIGREAL, dimension(:), optional, intent(in) :: beta
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jprsize,jpisize,jiobs
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      INTEGER :: nvpnull
      BIGREAL, dimension(:,:), allocatable :: matrr
CC----------------------------------------------------------------------
C
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
C
C --- allocation matrr
      allocate ( matrr(1:jprsize,1:jprsize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
C
      jrbasdeb=2
      jrbasfin=jprsize
      jrmatdeb=1
      jrmatfin=jprsize-1
C
C -1.- Compute the sum of the Gamma_i matrices
C --------------------------------------------
C
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
C
C -2.- Compute SVD decomposition of C
C -----------------------------------
C
      nvpnull=0
      CALL mkrs (matrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin),
     $     kmatUrr(jrmatdeb:jrmatfin,jrmatdeb:jrmatfin),
     $     klambda(jrbasdeb:jrbasfin),nvpnull)
C
C --- deallocate arrays
      IF (allocated(matrr)) deallocate(matrr)
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'utilroa','algoker_Gi')
 1001 CALL printerror2(0,1001,3,'utilroa','algoker_Gi')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocalcDi(kbasesr,kvectsinnov,kvectssqrdiagRi,
     $                      kvectspart,kdeltari)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute delta vectors for obs segments
CCC  -------
CCC  Method : delta_i = trsp(H_i Sf) inv(R_i) (y_i-H_i xf)
CCC  ------
CCC  Input :  kbasesr         : HSf
CCC  -----    kvectsinnov     : y-Hxf
CCC           kvectssqrdiagRi : inv(R)
CCC           kvectspart      : partition of the obs vector
CCC  Output : kdeltari        : delta_i
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
       BIGREAL, dimension(:,:), intent(in) :: kbasesr
       BIGREAL, dimension(:), intent(in) :: kvectsinnov
       BIGREAL, dimension(:), intent(in) :: kvectssqrdiagRi
       BIGREAL, dimension(:), intent(in) :: kvectspart
       BIGREAL, dimension(:,:), intent(out) :: kdeltari
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpssize,jprsize,jpisize,jiobs,jpu,js
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      BIGREAL, dimension(:), allocatable :: vects
      INTEGER, dimension(:), allocatable :: ptu, partidx
CC----------------------------------------------------------------------
C
      jpssize=size(kbasesr,1)
      jprsize=size(kbasesr,2)
      jpisize=size(kdeltari,2)
      IF (jpssize.NE.size(kvectsinnov,1)) GOTO 1000
      IF (jpssize.NE.size(kvectspart,1)) GOTO 1000
      IF (jpssize.NE.size(kvectssqrdiagRi,1)) GOTO 1000
      IF (jprsize.NE.size(kdeltari,1)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
C --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
C --- allocation ptu
      allocate ( ptu(1:jpssize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
C --- allocation partidx
      allocate ( partidx(1:jpssize) , stat=allocok )
      IF (allocok.GT.0) GOTO 1001
      partidx(:) = NINT(kvectspart(:))
C
      IF (MAXVAL(partidx).GT.jpisize) GOTO 1000
      IF (MINVAL(partidx).LT.1) GOTO 1000
C
C -0.- Initialisation
C -------------------
C
      jrbasdeb=2
      jrbasfin=jprend
      jrmatdeb=1
      jrmatfin=jprend-1
C
C -1.- Compute innovation vector and weight by observation error
C --------------------------------------------------------------
C
      vects(:) = kvectsinnov(:) *
     $           kvectssqrdiagRi(:) * kvectssqrdiagRi(:)
C
C -2.- Compute innovation in reduced space for all segments
C ---------------------------------------------------------
C
      DO jiobs=1,jpisize
C
        jpu = COUNT(partidx(1:jpssize).EQ.jiobs)
        ptu(:jpu) = PACK( (/ (js, js=1,jpssize) /) ,
     $              partidx(:) .EQ. jiobs )
C
        DO jrbas=jrbasdeb,jrbasfin
          kdeltari(jrbas,jiobs)=DOT_PRODUCT(kbasesr(ptu(:jpu),jrbas),
     $                                      vects(ptu(:jpu)))
        ENDDO
C
      ENDDO
C
C --- deallocate arrays
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(partidx)) deallocate(partidx)
      IF (allocated(ptu)) deallocate(ptu)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilroa','algocoef_u')
 1001 CALL printerror2(0,1001,3,'utilroa','algocoef_u')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algocalciRi(kvectsinnov,kvectssqrdiagRi,
     $                       kiRi,ky,kvectspart)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute iT R-1 i for each segment of the obs.
CCC  -------
CCC  Method : trsp(i) inv(R_i) (i)
CCC  ------
CCC  Input :  kvectsinnov     : y-Hxf
CCC  -----    kvectssqrdiagRi : inv(R)
CCC           kvectspart      : partition of the obs vector
CCC  Output : kiRi            : iT R-1 i
CCC  ------   ky              : number of observations
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpyend,jpoend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
       BIGREAL, dimension(:), intent(in) :: kvectsinnov
       BIGREAL, dimension(:), intent(in) :: kvectssqrdiagRi
       BIGREAL, dimension(:), intent(out) :: kiRi
       BIGREAL, dimension(:), intent(out) :: ky
       BIGREAL, dimension(:), optional, intent(in) :: kvectspart
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpssize,jpisize,jiobs,jpu,js
      BIGREAL, dimension(:), allocatable :: vects
      INTEGER, dimension(:), allocatable :: ptu, partidx
CC----------------------------------------------------------------------
C
      jpssize=size(kvectsinnov,1)
      jpisize=size(kiRi,1)
      IF (jpssize.NE.size(kvectssqrdiagRi,1)) GOTO 1000
      IF ((jpisize.GT.1).AND.(.NOT.present(kvectspart))) GOTO 1000
C --- allocation vects
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)

      IF (jpisize.GT.1) THEN
        IF (jpssize.NE.size(kvectspart,1)) GOTO 1000
C --- allocation ptu
        allocate ( ptu(1:jpssize) , stat=allocok )
        IF (allocok.GT.0) GOTO 1001
C --- allocation partidx
        allocate ( partidx(1:jpssize) , stat=allocok )
        IF (allocok.GT.0) GOTO 1001
        partidx(:) = NINT(kvectspart(:))
C
        IF (MAXVAL(partidx).GT.jpisize) GOTO 1000
        IF (MINVAL(partidx).LT.1) GOTO 1000
      ENDIF
C
C -1.- Compute innovation vector and weight by observation error
C --------------------------------------------------------------
C
      vects(:) = kvectsinnov(:) * kvectssqrdiagRi(:)
C
C -2.- Compute innovation in reduced space for all segments
C ---------------------------------------------------------
C
      IF (jpisize.GT.1) THEN

        DO jiobs=1,jpisize
C
          jpu = COUNT(partidx(1:jpssize).EQ.jiobs)
          ptu(:jpu) = PACK( (/ (js, js=1,jpssize) /) ,
     $                partidx(:) .EQ. jiobs )
C
          kiRi(jiobs)=DOT_PRODUCT(vects(ptu(:jpu)),
     $                            vects(ptu(:jpu)))
          ky(jiobs)=FREAL(jpu)
C
        ENDDO
C
      ELSE
C
        kiRi(1)=DOT_PRODUCT(vects(:),vects(:))
        ky(1)=FREAL(jpssize)
C
      ENDIF
C
C --- deallocate arrays
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(partidx)) deallocate(partidx)
      IF (allocated(ptu)) deallocate(ptu)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilroa','algocalciRi')
 1001 CALL printerror2(0,1001,3,'utilroa','algocalciRi')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE algosum_Di(kdeltari,kdelta,beta)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Sum up the delta vectors for obs segments
CCC  -------
CCC  Method : delta  =  SUM_i delta_i / beta_i
CCC  ------
CCC  Input :  kdeltari   : delta_i
CCC  -----    beta       : beta_i
CCC
CCC  Output : kdelta     : delta
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jprend,jpyend,jpoend
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
       BIGREAL, dimension(:,:), intent(in) :: kdeltari
       BIGREAL, dimension(:), intent(out) :: kdelta
       BIGREAL, dimension(:), optional, intent(in) :: beta
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jprsize,jpisize,jiobs
      INTEGER :: jrbasdeb,jrbasfin,jrbas,jrbas1,jrbas2
      INTEGER :: jrmatdeb,jrmatfin,jrmat,jrmat1,jrmat2
      INTEGER :: nvpnull
      BIGREAL, dimension(:,:), allocatable :: matrr
CC----------------------------------------------------------------------
C
      jprsize=size(kdeltari,1)
      jpisize=size(kdeltari,2)
      IF (jprsize.NE.size(kdelta,1)) GOTO 1000
      IF (jprsize.NE.jprend) GOTO 1000
      IF (present(beta)) THEN
        IF (jpisize.NE.size(beta,1)) GOTO 1000
      ENDIF
C
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
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'utilroa','algosum_Di')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkyorms (vectsrms,kflagxyo) 
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Fill observation error standard deviation vector
CCC  -------   using parameters in SESAM configuration
CCC
CCC  Method : Loop on variables and affect the corresponding
CCC  ------   error value to the corresponding segment of the Vy vector
CCC
CCC  Input :  kflagxyo : vector type (2=Vy, 3=Vo)
CCC  -----
CCC  Output : vectsrms : observation error standard deviation vector
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:), intent(out) :: vectsrms
      INTEGER, intent(in) :: kflagxyo
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: js,jsdeb,jsfin
      INTEGER :: jdta,inddta,jobs,indobs,inddbs
      BIGREAL :: valeur
CC----------------------------------------------------------------------
C
      SELECT CASE(kflagxyo)
      CASE(1)
         GOTO 1000
      CASE(2)
C Vy vector
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
               WRITE(numout,*) 'dtanam=',dta_nam(inddta)
     $              (1:lenv(dta_nam(inddta)))
               WRITE(numout,*) 'dtarms=',dta_rms(inddta)
            ENDIF
            jsdeb=jsfin+1
         ENDDO
      CASE(3)
C Vo vector
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
               WRITE(numout,*) 'obsnam=',obs_nam(indobs,inddbs)
     $              (1:lenv(obs_nam(indobs,inddbs)))
               WRITE(numout,*) 'obsrms=',obs_rms(indobs,inddbs)
            ENDIF
            jsdeb=jsfin+1
         ENDDO
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilroa','mkyorms')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkrs (kcov,vctp_rr,valp_r,nvpnull) 
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Compute eigenvalues and eigenvectors
CCC  -------   of a symmetric real matrix
CCC
CCC  Method : Use EISPACK routine 'rs'
CCC  ------
CCC  Input :  kcov : input matrix
CCC  -----
CCC  Output : valp_r : eigenvalues
CCC  ------   vctp_rr : eigenvectors
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use utilmath1
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      BIGREAL, dimension(:,:), intent(in) :: kcov
      BIGREAL, dimension(:,:), intent(out) :: vctp_rr
      BIGREAL, dimension(:), intent(out) :: valp_r
      INTEGER, intent(in) :: nvpnull
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jprsize,jr1,jrdeb,nvectindep,lrs
      BIGREAL, dimension(:), allocatable :: work1,work2
      INTEGER :: ierr,allocok
CC----------------------------------------------------------------------
C
C -0.- Initialisation
C -------------------
C
      jprsize = size(kcov,1)
      IF (jprsize .NE.size(kcov,2)) GOTO 1000
      IF (jprsize .NE.size(vctp_rr,1)) GOTO 1000
      IF (jprsize .NE.size(vctp_rr,2)) GOTO 1000
      IF (jprsize .NE.size(valp_r,1)) GOTO 1000
C
      allocate(work1(jprsize),stat=allocok)
      IF (allocok.GT.0) GOTO 1001
      allocate(work2(jprsize),stat=allocok)
      IF (allocok.GT.0) GOTO 1001
C
      valp_r(:) = FREAL(0.0)
      vctp_rr(:,:) = FREAL(0.0)
      nvectindep=jprsize-nvpnull
      jrdeb=nvpnull+1
C
C -1.- Computation of the eigenvalues and eigenvectors of the kcov matrix
C -----------------------------------------------------------------------
C
C --- EISPACK : routine package to solve eigenvalue problems
C --- rs      : eispack routine computing eigenvalues and eigenvectors
C ---           of a symmetric real matrix (see utilmath1.F)
C ---
C --- Usage of subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
C ---
C --- input
C ---
C ---       nm  : must be set to the row dimension of the two-dimensional
C ---             array parameters as declared in the calling program
C ---             dimension statement.
C ---
C ---       n   : is the order of the matrix  a.
C ---
C ---       a   : contains the real symmetric matrix.
C ---
C ---       matz: is an integer variable set equal to zero if
C ---             only eigenvalues are desired.  otherwise it is set to
C ---             any non-zero integer for both eigenvalues and eigenvectors.
C ---
C --- output
C ---
C ---       w   : contains the eigenvalues in ascending order.
C ---
C ---       z   : contains the eigenvectors if matz is not zero.
C ---
C ---       ierr: is an integer output variable set equal to an error
C ---             completion code described in the documentation for tqlrat
C ---             and tql2.  the normal completion code is zero.
C ---
C ---       fv1,fv2 :  are temporary storage arrays.
C ---
      lrs = 1
      CALL rs(jprsize,jprsize,kcov,valp_r,lrs,vctp_rr,work1,work2,ierr)
C
C -2.- Reorganisation of valp_r and vctp_rr arrays
C ------------------------------------------------
C Eigenvalues are sorted in descending order from jrdeb=2 to jprsize
C The first value (jrdeb-1=1) correspond to the smallest eigenvalue (equal to 0)
C
      work1(:)=valp_r(:)
      valp_r(jrdeb:jprsize)=work1(jprsize:jrdeb:-1)
C
      DO jr1=1,jprsize
         work1(:)=vctp_rr(jr1,:)
         vctp_rr(jr1,jrdeb:jprsize)=work1(jprsize:jrdeb:-1)
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilroa','mkrs')
 1001 CALL printerror2(0,1001,3,'utilroa','mkrs')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkcoeftodta (kfile,kvectz,kparty)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Save local z-vector (constant for
CCC  -------   every subsystems) in dta file
CCC
CCC  Method : Fill variables in Vy vector, with the value
CCC  ------   of the subsystem to which they belong
CCC            
CCC  Input :  kfile  : Name of the output file
CCC  -----    kvectz : z-vector
CCC           kparty : partition of the Vy object into subsystems
CCC
CCC  Output : none
CCC  ------
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_cfgxyo
      use mod_spacexyo , only : jpyend, jpxend, jprend, jpz
      use hioxyo
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfile
      BIGREAL, dimension(:), intent(in) :: kvectz
      BIGREAL, dimension(:), intent(in) :: kparty
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,jpysize,jpzsize
      BIGREAL, dimension(:), allocatable :: diagy
CC----------------------------------------------------------------------
C
      jpzsize=size(kvectz,1)
      jpysize=size(kparty,1)
      IF (jpysize.NE.jpyend) GOTO 1000
C
C --- allocation diagy
      allocate ( diagy(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      diagy(:) = FREAL(0.0)
CCC---------------------------------------------------------------------
C
C -1.- Fill Vy vector with the local quantity
C ---------------------------------------------
C
      diagy(1:jpysize) = kvectz(NINT(kparty(1:jpysize)))
C
C -2.- Save Vy vector in dta file
C -------------------------------
C
      CALL writedta(kfile,diagy(:))
C
C --- deallocation
      IF (allocated(diagy)) deallocate(diagy)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algoutilroa','mkcoeftodta')
 1001 CALL printerror2(0,1001,3,'algoutilroa','mkcoeftodta')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilroa
