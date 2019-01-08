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
C ---                    MKOBSTODTA.F                           ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-05 (C.E. Testut)                        ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE obstodta
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkobstodta
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC obstodta

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE obstodta (kfninobs,kfnoutdta)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
CCC  -------
CCC     build a DTA from an OBS
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
      use mod_spacexyo , only : jpyend,jpoend,jpitpend,spvaldta,
     $     gridijkobs,poscoefobs
      use hioxyo
      use utilmkh
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kfninobs,kfnoutdta
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable, save :: vecto
      BIGREAL, dimension(:), allocatable, save :: vecty
      BIGREAL, dimension(:), allocatable, save :: vectorms
C
      INTEGER :: allocok,jposize,jpitpsize,jpysize
      INTEGER :: flagcfg
CC----------------------------------------------------------------------
C
      jposize=jpoend
      jpitpsize=jpitpend
      jpysize=jpyend
C --- allocation vecty
      allocate ( vecty(1:jpysize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vecty(:) = FREAL(0.0)
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
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine obstodta   &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&'
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
C -2.- Read the obs file :
C ------------------------
C
      CALL readvalobs (kfninobs,vecto(:))
C
C -3.- make the dta file :
C ------------------------
C
      CALL mkhotoy(vecto(:),vecty(:),poscoefobs(:,:),spvaldta)
C
C -4.- write the dta file :
C -------------------------
C
      CALL writedta (kfnoutdta,vecty(:))
C
C --- deallocation
      IF (allocated(vecty)) deallocate(vecty)
      IF (allocated(vecto)) deallocate(vecto)
      IF (allocated(vectorms)) deallocate(vectorms)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mkobstodta','obstodta')
 1001 CALL printerror2(0,1001,3,'mkobstodta','obstodta')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkobstodta
