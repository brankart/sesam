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
! ---                  HIOCNT.F90                                 ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 01-06  (J.M. Brankart)                     ---
! --- revised      : 03-03  (J.M. Brankart)                     ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE evalhdrcnt : Read contour file header
! --- SUBROUTINE readcnt    : Read contour file
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE hiocnt
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC evalhdrcnt,readcnt

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE evalhdrcnt (kfilecnt)
!---------------------------------------------------------------------
!
!  Purpose : Read contour file header
!  ------
!  Method : Read first uncommented record of contour ASCII file
!  ------
!  Input :  kfilecnt  : filename
!  -----
!  Output : Common variables (described in module mod_cont):
!  ------   jpc, jppend, jplayend
!
!---------------------------------------------------------------------
! modules
      use mod_main
      use mod_cfgxyo
      use mod_cont
      use utilfiles
      IMPLICIT NONE
!
! header declarations
      CHARACTER(len=*), intent(in) :: kfilecnt
!
! local declarations
      CHARACTER(len=bgword) :: fline
!
! Open contour file
      CALL openfile(10,kfilecnt)
!
! Read contour file header (1st uncommented line)
!
      fline='#'
      DO WHILE (fline(1:1).EQ.'#')
         READ(10,'(a)') fline
      ENDDO
      READ(fline,*) jpc, jppend, jplayend
!
      CLOSE(10)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiocnt','evalhdrcnt')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE readcnt (kfilecnt,ktypcnt)
!---------------------------------------------------------------------
!
!  Purpose : Read contour file
!  -------
!  Method : Loop over contours and read each contour definition
!  ------
!  Input :  kfilecnt : contour filename
!  -----    ktypcnt  : contour type
!
!  Output : Common variables (in module mod_cont):
!  ------   jpp, jplay, contij, contidx, contlevmin, contlevmax
!
!---------------------------------------------------------------------
! modules
      use mod_main
      use mod_cfgxyo
      use mod_cont
      use utilfiles
      IMPLICIT NONE
!
! header declarations
      CHARACTER(len=*), intent(in) :: kfilecnt
      INTEGER, intent(in) :: ktypcnt
!---------------------------------------------------------------------
! local declarations
      INTEGER :: jc, jp, jlay
      CHARACTER(len=bgword) :: fline
!---------------------------------------------------------------------
!
! Open contour file
      CALL openfile(10,kfilecnt)
!
! Read contour file header (1st uncommented line)
!
      fline='#'
      DO WHILE (fline(1:1).EQ.'#')
         READ(10,'(a)') fline
      ENDDO
      READ(fline,*) jpc, jppend, jplayend
!
! Loop over contours
      DO jc = 1,jpc
!
! Read number of vertices defining each contour
         fline='#'
         DO WHILE (fline(1:1).EQ.'#')
            READ(10,'(a)') fline
         ENDDO
         READ(fline,*) jpp(jc)
!
! Read position of vertices defining each contour
         DO jp = 1,jpp(jc)
            fline='#'
            DO WHILE (fline(1:1).EQ.'#')
               READ(10,'(a)') fline
            ENDDO
            READ(fline,*) contij(jp,jc)%longi, contij(jp,jc)%latj
         ENDDO
!
! Read number of horizontal slices in cylinders
! defined by the polygonal contours
         fline='#'
         DO WHILE (fline(1:1).EQ.'#')
            READ(10,'(a)') fline
         ENDDO
         READ(fline,*) jplay(jc)
!
! Read minimum and maximum level defining each slice
! Read parameters of linear function (if ktypcnt=1)
! defined inside polygonal slice
         DO jlay = 1,jplay(jc)
            fline='#'
            DO WHILE (fline(1:1).EQ.'#')
               READ(10,'(a)') fline
            ENDDO
!
            SELECT CASE(ktypcnt)
            CASE(0)
               READ(fline,*) contidx(jlay,jc),  &
     &                   contlevmin(jlay,jc), contlevmax(jlay,jc)
            CASE(1)
               READ(fline,*) contidx(jlay,jc),  &
     &                   contlevmin(jlay,jc), contlevmax(jlay,jc), &
     &                   contf0(jlay,jc), contfx(jlay,jc), &
     &                   contfy(jlay,jc), contfz(jlay,jc)
            CASE DEFAULT
               GOTO 1000
            END SELECT
!
         ENDDO
      ENDDO
!
      CLOSE(10)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'hiocnt','readcnt')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE hiocnt
