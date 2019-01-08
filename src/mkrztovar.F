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
C ---                    MKZONTODTA.F                           ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2009-07 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C ---
C --- SUBROUTINE mkrztovar
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkrztovar
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC rztovar

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE rztovar(karginrz,karginpartvar,kargoutvar,kargincfg)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Interface an rz variable
CCC  -------   into Vx vector
CCC
CCC  Method : Read variable from rz file, paste it
CCC  ------   on Vx vector (using the partition),
CCC           and write Vx vector in var file
CCC
CCC  Input :  karginrz      : rz input filename
CCC  -----    karginpartvar : partition file
CCC           kargoutvar    : output Vx filename
CCC           kargincfg     : variable name
CCC
CCC---------------------------------------------------------------------
CC modules
CC =======
      use mod_main
      use mod_cfgxyo
      use mod_spacexyo , only : jpxend,jpz
      use hioxyo
      use liocrz
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: karginrz,karginpartvar,
     $                                kargoutvar,kargincfg
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable :: vectx
      BIGREAL, dimension(:), allocatable :: vectpart
      BIGREAL, dimension(:), allocatable :: vectz
C
      INTEGER :: allocok,jpxsize,jprsize,jpzsize,jpiobs
      INTEGER :: jx,jnxyo,flagxyo
      LOGICAL :: lectinfo
      CHARACTER(len=1) :: textexclusion
      CHARACTER(len=hgword) :: line
      CHARACTER(len=bgword) :: rznam
      INTEGER :: rzdim,rzidx
      BIGREAL, dimension(:,:), allocatable :: mat2d
CC----------------------------------------------------------------------
      jpxsize=jpxend
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/modzone/mkrztovar :'
         WRITE(numout,*) '         interface rz file to Vx vector'
      ENDIF
C
C Get jpz size from rz file
      CALL evalhdrcrz(karginrz,jprsize,jpzsize,jpiobs)
C
C Allocate Vx output array
      allocate ( vectx(0:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C Allocate Vx partition array
      allocate ( vectpart(0:jpxsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C Allocate Vz input array
      allocate ( vectz(1:jpzsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
C Read configuration file
      CALL openfile(numfil,kargincfg)
      textexclusion='#'
      line=readnextline(numfil,textexclusion)
      READ(line,*) rznam
      line=readnextline(numfil,textexclusion)
      READ(line,*) rzdim, rzidx
      CLOSE(numfil)
C
C Read Vz vector from rz file
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
C
C Read partition
      lectinfo=.FALSE. ; flagxyo=1 ; jnxyo=1
      CALL readvar(karginpartvar,vectpart(1:jpxsize),
     $             jnxyo,lectinfo,flagxyo)
C
C Paste Vz vetcor into Vx vector
      IF (NINT(MINVAL(vectpart)).LT.0) GOTO 101
      IF (NINT(MAXVAL(vectpart)).GT.jpzsize) GOTO 101
C
      vectx(1:jpxsize) = vectz(NINT(vectpart(1:jpxsize)))
C
C Write Vx vector in var file
      CALL writevar (kargoutvar,vectx(1:jpxsize),jnxyo)
C
C --- deallocation
      IF (allocated(vectx)) deallocate(vectx)
      IF (allocated(vectz)) deallocate(vectz)
      IF (allocated(vectpart)) deallocate(vectpart)
C
      RETURN
C
C --- error management
C
 1001 CALL printerror2(0,1001,3,'mkrztovar','zontovar')
C
 101  WRITE (texterror,*) 'Incompatible partition file'
      CALL printerror2(0,101,3,'mkrztovar','rztovar',
     $                 comment=texterror)
 102  WRITE (texterror,*) 'Bad dimension in cfg file'
      CALL printerror2(0,102,3,'mkrztovar','rztovar',
     $                 comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkrztovar
