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
! ---                    MKOBSTODTA.F90                           ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-05 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE obstodta
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkobstodta
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC obstodta

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE obstodta (kfninobs,kfnoutdta)
!---------------------------------------------------------------------
!
!  Purpose :
!  -------
!     build a DTA from an OBS
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
      use mod_spacexyo , only : jpyend,jpoend,jpitpend,spvaldta, &
     &     gridijkobs,poscoefobs
      use hioxyo
      use utilmkh
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kfninobs,kfnoutdta
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:), allocatable, save :: vecty
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      INTEGER :: allocok,jposize,jpitpsize,jpysize
      INTEGER :: flagcfg
!----------------------------------------------------------------------
!
      jposize=jpoend
      jpitpsize=jpitpend
      jpysize=jpyend
! --- allocation vecty
      allocate ( vecty(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecty(:) = FREAL(0.0)
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
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine obstodta   &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&'
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
! -2.- Read the obs file :
! ------------------------
!
      CALL readvalobs (kfninobs,vecto(:))
!
! -3.- make the dta file :
! ------------------------
!
      CALL mkhotoy(vecto(:),vecty(:),poscoefobs(:,:),spvaldta)
!
! -4.- write the dta file :
! -------------------------
!
      CALL writedta (kfnoutdta,vecty(:))
!
! --- deallocation
      IF (allocated(vecty)) deallocate(vecty)
      IF (allocated(vecto)) deallocate(vecto)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mkobstodta','obstodta')
 1001 CALL printerror2(0,1001,3,'mkobstodta','obstodta')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkobstodta
