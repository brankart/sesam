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
! ---                   ALGOVARI.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2007-03 (E. Cosme)                         ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE calcvari
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algovari
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcvari

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcvari (kinxyobas,koutxyo,kflagxyo,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose : Compute variances from reduced order covariance matrix
!
!  Method : Read reduced order covariance matrix from input
!  ------   directory (kinxyobas), compute variances for
!           all variables, write the result in output
!           file (koutxyo).
!
!  Input : kinxyobas   : Cxyo input directory
!  -----   koutxyo     : Vxyo output file
!          kflagxyo    : Vector type (1=Vx,2=Vy,3=Vo)
!          kconfigo    : Observation operator (Vo) inputy file
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : spvalvar,spvaldta,spvalobs, &
     &     jpoend,jpitpend,jpx,jpxend,jpyend,jprend, &
     &     poscoefobs,gridijkobs,arraynx_jindxbeg
      use hioxyo
      use hiobas
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinxyobas,koutxyo,kconfigo
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:,:), allocatable, save :: basesr
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:), allocatable, save :: vectorms
!
      BIGREAL, dimension(:), allocatable :: vectr
      INTEGER :: allocok,jpssize,jpysize,jpitpsize,jprsize
      INTEGER :: jnxyo,flagcfg,js,jr
      BIGREAL :: spval
      LOGICAL :: lectinfo
      INTEGER :: valbase, jrbasdeb,jrbasfin
!----------------------------------------------------------------------
!
      jprsize=jprend
      jpysize=jpyend
      jpitpsize=1
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modvari/algovari :'
         WRITE(numout,*) '         compute variances (Vxyo) from'
         WRITE(numout,*) '         reduced order cov. matrix (Cxyo)'
      ENDIF
!
      SELECT CASE (kflagxyo)
      CASE (1)
         jpssize=jpx
         spval=spvalvar
      CASE (2)
         jpssize=jpyend
         spval=spvaldta
      CASE (3)
         jpssize=jpoend
         jpitpsize=jpitpend
         spval=spvalobs
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Allocate Vxyo arrays
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
!
! Allocate Cxyo array
      allocate ( basesr(1:jpssize,1:jprend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      basesr(:,:) = FREAL(0.0)
!
      IF (kflagxyo.EQ.3) THEN
!
! Allocate poscoefobs array
        allocate ( poscoefobs(1:jpssize,1:jpitpsize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        poscoefobs(:,:) = type_poscoef(0,FREAL(0.0))
!
! Allocate gridijkobs array
        allocate ( gridijkobs(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        gridijkobs(:)=type_gridijk(FREAL(0.0),FREAL(0.0),FREAL(0.0))
!
! Allocate vectorms array
        allocate ( vectorms(1:jpssize), stat=allocok )
        IF (allocok.NE.0) GOTO 1001
        vectorms(:) = FREAL(0.0)
!
! Read poscoefobs, vectorms and gridijkobs arrays
        flagcfg=1
        CALL readcfgobs (kconfigo,flagcfg, &
     &        kvectorms=vectorms(:))
        flagcfg=2
        CALL readcfgobs (kconfigo,flagcfg, &
     &        kgridijkobs=gridijkobs(:))
        flagcfg=3
        CALL readcfgobs (kconfigo,flagcfg, &
     &        kposcoefobs=poscoefobs(:,:))
      ENDIF
!
! Read header information from reduced order cov. matrix directory
      CALL readinfobas(kinxyobas,valbase)
      IF (valbase.LT.1) GOTO 1000
!
      DO jnxyo=1,limjpnxyo(MIN(3,kflagxyo))
!
! -1.- Read input reduced order covariance matrix
! -----------------------------------------------
!
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING the covariance matrix'
        ENDIF
!
        jrbasdeb=1
        jrbasfin=jprsize
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxyobas,basesr(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
!
! -2.- Compute variances
! ----------------------
!
        IF (limjpnxyo(MIN(3,kflagxyo)).EQ.1) THEN
          DO jr=jrbasdeb+1,jrbasfin
            vects(:) = vects(:) + basesr(:,jr)*basesr(:,jr)
          ENDDO
        ELSE
          DO jr=jrbasdeb+1,jrbasfin
            vects(:) = vects(:) + basesr(:,jr)*basesr(:,jr)
          ENDDO
        ENDIF
!     
! -3.- Write variances
! --------------------
!
        SELECT CASE (kflagxyo)
        CASE (1)
           CALL writevar (koutxyo,vects(:),jnxyo)
        CASE (2)
           CALL writedta (koutxyo,vects(:))
        CASE (3)
           CALL writeobs (koutxyo,vects(:),vectorms(:), &
     &              gridijkobs(:),poscoefobs(:,:))
        CASE DEFAULT
           GOTO 1000
        END SELECT
!
      ENDDO
!
! --- deallocation
      IF (allocated(vects)) deallocate(vects)
      IF (allocated(basesr)) deallocate(basesr)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algovari','calcvari')
 1001 CALL printerror2(0,1001,3,'algovari','calcvari')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algovari
