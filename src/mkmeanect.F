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
C ---                    MKMEANECT.F                            ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12 (C.E. Testut)                        ---
C --- modification : 99-05 (C.E. Testut)                        ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE meanect
C --- SUBROUTINE meanstdbub
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkmeanect
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC meanect,meanstdnbub

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE meanect(kvectsmean,kvectsect,kbasesr)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
CCC  -------
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
      use mod_spacexyo , only : jpxend
      use utilclc
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
       BIGREAL, dimension(:), intent(out) :: kvectsmean,kvectsect
       BIGREAL, dimension(:,:), intent(in) :: kbasesr
CC----------------------------------------------------------------------
CC local declarations
CC ==================
        INTEGER jr,jprsize,jpssize
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '      *******************'
         WRITE(numout,*) '      * routine meanect *'
         WRITE(numout,*) '      *******************'
         WRITE(numout,*) ' '
      ENDIF
C
C -0.- Initialisation :
C ---------------------
C
      jpssize = size(kbasesr,1)
      jprsize = size(kbasesr,2)
C
C -1.- Compute the mean vector :
C ------------------------------
C
      CALL mkmean (kvectsmean(:),kbasesr(:,:))
C
C -2.- Compute the standard deviation :
C -------------------------------------
C
      kvectsect(:) = -(kvectsmean(:))**2
C
      DO jr=1,jprsize 
         kvectsect(:) = kvectsect(:) +
     $           ((kbasesr(:,jr)*kbasesr(:,jr))
     $           /FREAL(jprsize))
      ENDDO
C
C      DO jr=kjrbasdeb,kjrbasfin
C         DO js=1,kjps
C            kvectsect(js)=kvectsect(js)+
C     $           ((kbasesr(js,jr)*kbasesr(js,jr))
C     $           /FREAL(kjrbasfin-kjrbasdeb+1))
C	 ENDDO
C      ENDDO
C
      kvectsect(:)=SQRT(ABS(kvectsect(:)))
C
      RETURN
C
 10   FORMAT(2X,A9,2(1X,"|",2X,A9,1X)) 
 11   FORMAT(1X,A10,1X,"|",1X,2(1PE11.5E2,3X))
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mkmeanect','meanect')
 1001 CALL printerror2(0,1001,3,'mkmeanect','meanect')
C
 101  WRITE (texterror,*) 'error in argument of the routine :'
      CALL printerror2(0,101,3,'mkmeanect','meanect',
     $     comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE meanstdnbub(kbubmean,kbubect,kbubr)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
CCC  -------
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
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ==================
      BIGREAL, dimension(:,:,:,:,:), intent(out) :: kbubmean,kbubect
      BIGREAL, dimension(:,:,:,:,:,:), intent(in) :: kbubr
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER jr,jprsize,jpisize,jpjsize,jpksize,jptsize,jpbubsize
CC----------------------------------------------------------------------
C
      IF (nprint.GE.2) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '      **********************'
         WRITE(numout,*) '      * routine meanstdbub *'
         WRITE(numout,*) '      **********************'
         WRITE(numout,*) ' '
      ENDIF
C
C -0.- Initialisation :
C ---------------------
C
      jpisize = size(kbubr,1)
      jpjsize = size(kbubr,2)
      jpksize = size(kbubr,3)
      jptsize = size(kbubr,4)
      jpbubsize = size(kbubr,5)
      jprsize = size(kbubr,6)
C
C -1.- Compute the mean vector :
C ------------------------------
C
      kbubmean(:,:,:,:,:)=FREAL(0.0)
!CDIR NODEP
      DO jr=1,jprsize
         kbubmean(:,:,:,:,:) = kbubmean(:,:,:,:,:) + kbubr(:,:,:,:,:,jr)
      ENDDO
      kbubmean(:,:,:,:,:) = kbubmean(:,:,:,:,:) / jprsize
C
C -2.- Compute the standard deviation :
C -------------------------------------
C
      kbubmean(:,:,:,:,:) = -(kbubmean(:,:,:,:,:))**2
C
!CDIR NODEP
      DO jr=1,jprsize 
         kbubect(:,:,:,:,:) = kbubect(:,:,:,:,:) +
     $           ((kbubr(:,:,:,:,:,jr)*kbubr(:,:,:,:,:,jr))
     $           /FREAL(jprsize))
      ENDDO
C
      kbubect(:,:,:,:,:)=SQRT(ABS(kbubect(:,:,:,:,:)))
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mkmeanect','mkmeanstdbub')
 1001 CALL printerror2(0,1001,3,'mkmeanect','mkmeanstdbub')
C
 101  WRITE (texterror,*) 'error in argument of the routine :'
      CALL printerror2(0,101,3,'mkmeanect','mkmeanstdbub',
     $     comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkmeanect
