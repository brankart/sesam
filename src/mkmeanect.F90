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
! ---                    MKMEANECT.F90                            ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12 (C.E. Testut)                        ---
! --- modification : 99-05 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE meanect
! --- SUBROUTINE meanstdbub
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkmeanect
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC meanect,meanstdnbub

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE meanect(kvectsmean,kvectsect,kbasesr)
!---------------------------------------------------------------------
!
!  Purpose :
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
      use mod_spacexyo , only : jpxend
      use utilclc
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
       BIGREAL, dimension(:), intent(out) :: kvectsmean,kvectsect
       BIGREAL, dimension(:,:), intent(in) :: kbasesr
!----------------------------------------------------------------------
! local declarations
! ==================
        INTEGER jr,jprsize,jpssize
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '      *******************'
         WRITE(numout,*) '      * routine meanect *'
         WRITE(numout,*) '      *******************'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      jpssize = size(kbasesr,1)
      jprsize = size(kbasesr,2)
!
! -1.- Compute the mean vector :
! ------------------------------
!
      CALL mkmean (kvectsmean(:),kbasesr(:,:))
!
! -2.- Compute the standard deviation :
! -------------------------------------
!
      kvectsect(:) = -(kvectsmean(:))**2
!
      DO jr=1,jprsize 
         kvectsect(:) = kvectsect(:) + &
     &           ((kbasesr(:,jr)*kbasesr(:,jr)) &
     &           /FREAL(jprsize))
      ENDDO
!
!      DO jr=kjrbasdeb,kjrbasfin
!         DO js=1,kjps
!            kvectsect(js)=kvectsect(js)+
!     $           ((kbasesr(js,jr)*kbasesr(js,jr))
!     $           /FREAL(kjrbasfin-kjrbasdeb+1))
!	 ENDDO
!      ENDDO
!
      kvectsect(:)=SQRT(ABS(kvectsect(:)))
!
      RETURN
!
 10   FORMAT(2X,A9,2(1X,"|",2X,A9,1X)) 
 11   FORMAT(1X,A10,1X,"|",1X,2(1PE11.5E2,3X))
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkmeanect','meanect')
 1001 CALL printerror2(0,1001,3,'mkmeanect','meanect')
!
 101  WRITE (texterror,*) 'error in argument of the routine :'
      CALL printerror2(0,101,3,'mkmeanect','meanect', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE meanstdnbub(kbubmean,kbubect,kbubr)
!---------------------------------------------------------------------
!
!  Purpose :
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
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ==================
      BIGREAL, dimension(:,:,:,:,:), intent(out) :: kbubmean,kbubect
      BIGREAL, dimension(:,:,:,:,:,:), intent(in) :: kbubr
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER jr,jprsize,jpisize,jpjsize,jpksize,jptsize,jpbubsize
!----------------------------------------------------------------------
!
      IF (nprint.GE.2) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '      **********************'
         WRITE(numout,*) '      * routine meanstdbub *'
         WRITE(numout,*) '      **********************'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      jpisize = size(kbubr,1)
      jpjsize = size(kbubr,2)
      jpksize = size(kbubr,3)
      jptsize = size(kbubr,4)
      jpbubsize = size(kbubr,5)
      jprsize = size(kbubr,6)
!
! -1.- Compute the mean vector :
! ------------------------------
!
      kbubmean(:,:,:,:,:)=FREAL(0.0)
!CDIR NODEP
      DO jr=1,jprsize
         kbubmean(:,:,:,:,:) = kbubmean(:,:,:,:,:) + kbubr(:,:,:,:,:,jr)
      ENDDO
      kbubmean(:,:,:,:,:) = kbubmean(:,:,:,:,:) / jprsize
!
! -2.- Compute the standard deviation :
! -------------------------------------
!
      kbubmean(:,:,:,:,:) = -(kbubmean(:,:,:,:,:))**2
!
!CDIR NODEP
      DO jr=1,jprsize 
         kbubect(:,:,:,:,:) = kbubect(:,:,:,:,:) + &
     &           ((kbubr(:,:,:,:,:,jr)*kbubr(:,:,:,:,:,jr)) &
     &           /FREAL(jprsize))
      ENDDO
!
      kbubect(:,:,:,:,:)=SQRT(ABS(kbubect(:,:,:,:,:)))
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkmeanect','mkmeanstdbub')
 1001 CALL printerror2(0,1001,3,'mkmeanect','mkmeanstdbub')
!
 101  WRITE (texterror,*) 'error in argument of the routine :'
      CALL printerror2(0,101,3,'mkmeanect','mkmeanstdbub', &
     &     comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkmeanect
