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
! ---                   ALGORANK.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2021-09 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE algorank
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algorank
      use mod_main
      use ensdam_score_ranks
      IMPLICIT NONE
      PRIVATE

      PUBLIC calcrank

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE calcrank (kinxyo,koutxyo,kinxbasref, &
     &                     kflagxyo,kconfigo)
!---------------------------------------------------------------------
!
!  Purpose : Compute rank of data within input ensemble
!
!  Method : Read ensemble from input
!  ------   directory (kinxbasref), read vector of data
!           from input file (kinxyo)
!           write the resultng ranks in output file (koutxyo).
!
!  Input : kinxyo      : Vxyo input file
!  -----   koutxyo     : Vxyo output file
!          kinxbasref  : Cxyo ensemble directory
!          kflagxyo    : Vector type (1=Vx,2=Vy,3=Vo)
!          kconfigo    : Observation operator (Vo) inputy file
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : &
     &     jpoend,jpitpend,jpx,jpxend,jpyend,jprend, &
     &     poscoefobs,gridijkobs,arraynx_jindxbeg
      use hioxyo
      use hiobas
      use ensdam_anatra
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kinxyo,koutxyo,kinxbasref
      CHARACTER(len=*), intent(in) :: kconfigo
      INTEGER, intent(in) :: kflagxyo
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable, save :: vects
      BIGREAL, dimension(:,:), allocatable, save :: ens
      INTEGER, dimension(:), allocatable :: ranks
      BIGREAL, dimension(:), allocatable :: vectorms
!
      INTEGER :: allocok,jpssize,jpitpsize,jprsize
      INTEGER :: jnxyo,js,jr,jp,jjproc
      LOGICAL :: lectinfo,lmodprint
      INTEGER :: valbase,jrbasdeb,jrbasfin,flagcfg
!----------------------------------------------------------------------
!
      jprsize=jprend
      jpitpsize=1
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modrank/algorank :'
         WRITE(numout,*) '         compute ranks'
      ENDIF
!
      SELECT CASE (kflagxyo)
      CASE (1)
         jpssize=jpx
      CASE (2)
         jpssize=jpyend
      CASE (3)
         jpssize=jpoend
         jpitpsize=jpitpend
      CASE DEFAULT
         GOTO 1000
      END SELECT
!
! Allocate Cxyo ensemble array
      allocate ( ens(1:jpssize,1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ens(:,:) = FREAL(0.0)
!
! Allocate Vxyo arrays for verification data
      allocate ( vects(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      vects(:) = FREAL(0.0)
!
! Allocate Vxyo arrays for ranks
      allocate ( ranks(1:jpssize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      ranks(:) = 0
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
      IF (jpproc.GT.limjpnxyo(MIN(3,kflagxyo))) GOTO 1003
      DO jnxyo=1+jproc,limjpnxyo(MIN(3,kflagxyo)),jpproc
!
! -1.- Read input ensemble
! ------------------------
!
        IF ((nprint.GE.2).AND.(jnxyo.EQ.1)) THEN
          WRITE(numout,*) '    ==> READING input ensemble'
        ENDIF
!
        jrbasdeb=1
        jrbasfin=jprsize
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readbas(kinxbasref,ens(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo)
        CASE (3)
          CALL readbas(kinxbasref,ens(:,:),jnxyo,jrbasdeb,jrbasfin, &
     &           lectinfo,kflagxyo,poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT
!     
! -4.- Read vector of verification data
! -------------------------------------
!
        lectinfo=.FALSE.
        SELECT CASE (kflagxyo)
        CASE (1,2)
          CALL readxyo (kinxyo,vects(:),jnxyo,lectinfo,kflagxyo)
        CASE (3)
          CALL readxyo (kinxyo,vects(:),jnxyo,lectinfo,kflagxyo, &
     &                  poscoefobs(:,:))
        CASE DEFAULT
          GOTO 1000
        END SELECT

!
! -3.- Compute ranks
! ------------------
!
        CALL compute_ranks( ens, vects, ranks )
        vects(:) = ranks(:)
!     
! -4.- Write output vector
! ------------------------
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
      IF (allocated(ranks)) deallocate(ranks)
      IF (allocated(ens)) deallocate(ens)
      IF (allocated(poscoefobs)) deallocate(poscoefobs)
      IF (allocated(gridijkobs)) deallocate(gridijkobs)
      IF (allocated(vectorms)) deallocate(vectorms)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algorank','algorank')
 1001 CALL printerror2(0,1001,3,'algorank','algorank')
 1003 CALL printerror2(0,1003,3,'algorank','algorank')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algorank
