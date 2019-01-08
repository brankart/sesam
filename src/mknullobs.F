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
C ---                    MKNULLOBS.F                            ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 99-12 (C.E. Testut)                        ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  nullobs
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mknullobs
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC nullobs

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE nullobs(kargnullobs,kargoutobs)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Null Observation operator
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
      use lioobs
      use liocobs
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: kargnullobs,kargoutobs
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      LOGICAL :: found
      INTEGER :: jobs,indobs,inddbs,spos,jobsfound,jextobs
      CHARACTER(len=bgword) :: fname
      INTEGER :: ltext
CC----------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine mknullobs &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
C
C -I.- searching indobs inddbs for creation of trap file with jpoloc=0 :
C --------------------------------------------------------------------
C
         found=.FALSE.
         LOOP3 : DO jobs=1,obsend
            indobs=obs_ord(jobs)
            inddbs=obsnord(jobs)
            ltext=lenv(obs_nam(indobs,inddbs))
            IF ((kargnullobs(1:ltext).EQ.obs_nam(indobs,inddbs)(1:ltext))
     $           .AND.((kargnullobs((ltext+1):(ltext+1))).EQ.(' '))) THEN
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
C
         jextobs=indext(kargoutobs,extobstab,nbextobs)
         IF (extobsunit(jextobs)) GOTO 1000
C
         IF (found) THEN
            indobs=obs_ord(jobsfound)
            inddbs=obsnord(jobsfound)
            spos = posit(kargoutobs,etoile)
            WRITE (fname,'(a,a,a)') kargoutobs(1:(spos-1)),
     $           obsonam(indobs,inddbs)
     $           (1:lenv(obsonam(indobs,inddbs))),
     $           kargoutobs((spos+1):lenv(kargoutobs))
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
               WRITE (fname,'(a,a,a)') kargoutobs(1:(spos-1)),
     $              obsonam(indobs,inddbs)
     $              (1:lenv(obsonam(indobs,inddbs))),
     $              kargoutobs((spos+1):lenv(kargoutobs))
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
               WRITE(numout,12) jobs,obs_ord(jobs),obsnord(jobs),
     $              obs_nam(indobs,inddbs),
     $              obs_dim(indobs,inddbs),obs_nbr(indobs,inddbs),
     $              obs_moy(indobs,inddbs),obs_ect(indobs,inddbs)
            ENDDO
            WRITE(numout,10) '-'
            WRITE(numout,*) ' '
         ENDIF
C
      RETURN
C
C --- format definitions
C
 10   FORMAT (A1,64("-"))
 11   FORMAT ("|",8(2X,A4,2X,"|"))
 12   FORMAT ("|",3(I5,3X,"|"),3X,A5,"|",2X,I3,"D",2X,"|",
     $     I7,1X,"|",2(E7.1,1X,"|"))
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mknullobs','nullobs')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mknullobs
