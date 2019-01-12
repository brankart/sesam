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
! ---                    MKNULLOBS.F90                            ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 99-12 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  nullobs
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mknullobs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC nullobs

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE nullobs(kargnullobs,kargoutobs)
!---------------------------------------------------------------------
!
!  Purpose : Null Observation operator
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
      use lioobs
      use liocobs
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: kargnullobs,kargoutobs
!----------------------------------------------------------------------
! local declarations
! ==================
      LOGICAL :: found
      INTEGER :: jobs,indobs,inddbs,spos,jobsfound,jextobs
      CHARACTER(len=bgword) :: fname
      INTEGER :: ltext
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mknullobs &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
!
! -I.- searching indobs inddbs for creation of trap file with jpoloc=0 :
! --------------------------------------------------------------------
!
         found=.FALSE.
         LOOP3 : DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltext=lenv(obs_nam(indobs,inddbs))
            IF ((kargnullobs(1:ltext).EQ.obs_nam(indobs,inddbs)(1:ltext)) &
     &           .AND.((kargnullobs((ltext+1):(ltext+1))).EQ.(' '))) THEN
               found=.TRUE.
               jobsfound=jobs
               EXIT LOOP3
            ENDIF            
         ENDDO LOOP3
         DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            obs_ind(indobs,inddbs)=0
            obs_nbr(indobs,inddbs)=0
         ENDDO
!
         jextobs=indext(kargoutobs,extobstab,nbextobs)
         IF (extobsunit(jextobs)) GOTO 1000
!
         IF (found) THEN
            indobs=obs_ord(jobsfound)
            inddbs=obsnord(jobsfound)
            spos = posit(kargoutobs,etoile)
            WRITE (fname,'(a,a,a)') kargoutobs(1:(spos-1)), &
     &           obsonam(indobs,inddbs) &
     &           (1:lenv(obsonam(indobs,inddbs))), &
     &           kargoutobs((spos+1):lenv(kargoutobs))
            SELECT CASE (jextobs)
            CASE (1)
              CALL writehdrfileobs (fname,jobsfound)
            CASE (2)
              CALL writehdrfilecobs(fname,jobsfound)
            END SELECT
         ELSE
            DO jobs=1,obsend
               indobs=obs_ord(jobs)
               inddbs=obsnord(jobs)
               spos = posit(kargoutobs,etoile)
               WRITE (fname,'(a,a,a)') kargoutobs(1:(spos-1)), &
     &              obsonam(indobs,inddbs) &
     &              (1:lenv(obsonam(indobs,inddbs))), &
     &              kargoutobs((spos+1):lenv(kargoutobs))
               SELECT CASE (jextobs)
               CASE (1)
                 CALL writehdrfileobs (fname,jobs)
               CASE (2)
                 CALL writehdrfilecobs(fname,jobs)
               END SELECT
            ENDDO
         ENDIF
         IF (nprint.GE.1) THEN
            WRITE(numout,10) '-'
            WRITE(numout,11) 'ord','ind','obs','ndbs','dim','nmax','moy','ect'
            WRITE(numout,10) '-'
            DO jobs=1,obsend
               indobs=obs_ord(jobs)
               inddbs=obsnord(jobs)
               WRITE(numout,12) jobs,obs_ord(jobs),obsnord(jobs), &
     &              obs_nam(indobs,inddbs), &
     &              obs_dim(indobs,inddbs),obs_nbr(indobs,inddbs), &
     &              obs_moy(indobs,inddbs),obs_ect(indobs,inddbs)
            ENDDO
            WRITE(numout,10) '-'
            WRITE(numout,*) ' '
         ENDIF
!
      RETURN
!
! --- format definitions
!
 10   FORMAT (A1,64("-"))
 11   FORMAT ("|",8(2X,A4,2X,"|"))
 12   FORMAT ("|",3(I5,3X,"|"),3X,A5,"|",2X,I3,"D",2X,"|", &
     &     I7,1X,"|",2(E7.1,1X,"|"))
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'mknullobs','nullobs')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mknullobs
