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
! ---                    MKZONTODTA.F90                           ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2009-07 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! ---
! --- SUBROUTINE mkrztovar
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkrztovar
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC rztovar

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE rztovar(karginrz,karginpartvar,kargoutvar,kargincfg)
!---------------------------------------------------------------------
!
!  Purpose : Interface an rz variable
!  -------   into Vx vector
!
!  Method : Read variable from rz file, paste it
!  ------   on Vx vector (using the partition),
!           and write Vx vector in var file
!
!  Input :  karginrz      : rz input filename
!  -----    karginpartvar : partition file
!           kargoutvar    : output Vx filename
!           kargincfg     : variable name
!
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpz
      use hioxyo
      use liocrz
      use utilfiles
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: karginrz,karginpartvar, &
     &                                kargoutvar,kargincfg
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable :: vectx
      BIGREAL, dimension(:), allocatable :: vectpart
      BIGREAL, dimension(:), allocatable :: vectz
!
      INTEGER :: allocok,jpxsize,jprsize,jpzsize,jpiobs
      INTEGER :: jx,jnxyo,flagxyo
      LOGICAL :: lectinfo
      CHARACTER(len=1) :: textexclusion
      CHARACTER(len=hgword) :: line
      CHARACTER(len=bgword) :: rznam
      INTEGER :: rzdim,rzidx
      BIGREAL, dimension(:,:), allocatable :: mat2d
!----------------------------------------------------------------------
      jpxsize=jpxend
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modzone/mkrztovar :'
         WRITE(numout,*) '         interface rz file to Vx vector'
      ENDIF
!
! Get jpz size from rz file
      CALL evalhdrcrz(karginrz,jprsize,jpzsize,jpiobs)
!
! Allocate Vx output array
      allocate ( vectx(0:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
! Allocate Vx partition array
      allocate ( vectpart(0:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
! Allocate Vz input array
      allocate ( vectz(1:jpzsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
!
! Read configuration file
      CALL openfile(numfil,kargincfg)
      textexclusion='#'
      line=readnextline(numfil,textexclusion)
      READ(line,*) rznam
      line=readnextline(numfil,textexclusion)
      READ(line,*) rzdim, rzidx
      CLOSE(numfil)
!
! Read Vz vector from rz file
      SELECT CASE (rzdim)
      CASE(1)
        CALL readveccrz(karginrz,rznam,vectz(1:jpzsize))
      CASE(2)
        allocate(mat2d(1:jpzsize,1:jpiobs),stat=allocok)
        IF (allocok.NE.0) GOTO 1001
        CALL readmat2crz(karginrz,rznam,mat2d(:,:))
        vectz(1:jpzsize)=mat2d(1:jpzsize,rzidx)
        deallocate(mat2d)
      CASE DEFAULT
        GOTO 102
      END SELECT
!
! Read partition
      lectinfo=.FALSE. ; flagxyo=1 ; jnxyo=1
      CALL readvar(karginpartvar,vectpart(1:jpxsize), &
     &             jnxyo,lectinfo,flagxyo)
!
! Paste Vz vetcor into Vx vector
      IF (NINT(MINVAL(vectpart)).LT.0) GOTO 101
      IF (NINT(MAXVAL(vectpart)).GT.jpzsize) GOTO 101
!
      vectx(1:jpxsize) = vectz(NINT(vectpart(1:jpxsize)))
!
! Write Vx vector in var file
      CALL writevar (kargoutvar,vectx(1:jpxsize),jnxyo)
!
! --- deallocation
      IF (allocated(vectx)) deallocate(vectx)
      IF (allocated(vectz)) deallocate(vectz)
      IF (allocated(vectpart)) deallocate(vectpart)
!
      RETURN
!
! --- error management
!
 1001 CALL printerror2(0,1001,3,'mkrztovar','zontovar')
!
 101  WRITE (texterror,*) 'Incompatible partition file'
      CALL printerror2(0,101,3,'mkrztovar','rztovar', &
     &                 comment=texterror)
 102  WRITE (texterror,*) 'Bad dimension in cfg file'
      CALL printerror2(0,102,3,'mkrztovar','rztovar', &
     &                 comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkrztovar
