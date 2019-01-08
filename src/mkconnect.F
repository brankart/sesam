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
C ---                    MKCONNECT.F                            ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 02-02 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE connect_init
C --- SUBROUTINE connect_close
C --- SUBROUTINE mkconnecty
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE mkconnect
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC connect_init,connect_close,mkconnecty

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE connect_init
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
CCC  -------
CCC     Initialize arrays for management of grid connections
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
      use mod_mask
      use mod_coord
      use utilfiles
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: allocok,numfila
      INTEGER :: maxjpi,maxjpj,maxjpk,maxjpt,jcon,jline,jplcon,jlcon
      INTEGER :: ji,jj,jk,jt,jy,jdta,inddta,inddtamsk,jlconmax
      INTEGER :: jimin,jimax,jjmin,jjmax
      CHARACTER(len=hgword) :: line
      CHARACTER(len=1) :: textexclusion
      CHARACTER(len=3) :: key,val1,val2
CCC-------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine connect_init &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
C
      maxjpi=MAXVAL(dta_jpi(dta_ord(1:dtaend)))
      maxjpj=MAXVAL(dta_jpj(dta_ord(1:dtaend)))
      maxjpk=MAXVAL(dta_jpk(dta_ord(1:dtaend)))
      maxjpt=MAXVAL(dta_jpt(dta_ord(1:dtaend)))
C
C Read connection description file
C
      textexclusion='#'
      numfila=10
      CALL openfile(numfila,argconnect)
C
      line=readnextline(numfila,textexclusion)
      READ(line,*) jpcon
C
      IF (jpcon.GT.0) THEN
C
C --- allocation dir1con
         allocate ( dir1con(1:jpcon,2), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C --- allocation dir2con
         allocate ( dir2con(1:jpcon,2), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C --- allocation signcon
         allocate ( signcon(1:jpcon,2), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C --- allocation idxcon
         allocate ( idxcon(1:jpcon,2), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C --- allocation mincon
         allocate ( mincon(1:jpcon,2), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C --- allocation maxcon
         allocate ( maxcon(1:jpcon,2), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C
      ENDIF
C
      DO jcon=1,jpcon
         line=readnextline(numfila,textexclusion)
         READ(line,'(A,x,A,x,A)') key,val1,val2
         IF (key(1:3).NE.'dir') GOTO 102
C 
         SELECT CASE(val1(1:1))
         CASE('i')
            dir1con(jcon,1)=1
            dir2con(jcon,1)=2
         CASE('j')
            dir1con(jcon,1)=2
            dir2con(jcon,1)=1
         CASE DEFAULT
            GOTO 102
         END SELECT
C
         SELECT CASE(val1(2:2))
         CASE('+')
            signcon(jcon,1)=+1
         CASE('-')
            signcon(jcon,1)=-1
         CASE DEFAULT
            GOTO 102
         END SELECT

         SELECT CASE(val2(1:1))
         CASE('i')
            dir1con(jcon,2)=1
            dir2con(jcon,2)=2
         CASE('j')
            dir1con(jcon,2)=2
            dir2con(jcon,2)=1
         CASE DEFAULT
            GOTO 102
         END SELECT
C
         SELECT CASE(val2(2:2))
         CASE('+')
            signcon(jcon,2)=+1
         CASE('-')
            signcon(jcon,2)=-1
         CASE DEFAULT
            GOTO 102
         END SELECT
C
         line=readnextline(numfila,textexclusion)
         READ(line,'(A,x,A,x,A)') key,val1,val2
         IF (key(1:3).NE.'idx') GOTO 102
         IF (val1(1:3).EQ.'min') THEN
            idxcon(jcon,1)=1
         ELSEIF (val1(1:3).EQ.'max') THEN
            IF (dir1con(jcon,1).EQ.1)THEN
               idxcon(jcon,1)=maxjpi
            ELSE
               idxcon(jcon,1)=maxjpj
            ENDIF
         ELSE
            READ(val1,*) idxcon(jcon,1)
         ENDIF
         IF (val2(1:3).EQ.'min') THEN
            idxcon(jcon,2)=1
         ELSEIF (val2(1:3).EQ.'max') THEN
            IF (dir1con(jcon,2).EQ.1)THEN
               idxcon(jcon,2)=maxjpi
            ELSE
               idxcon(jcon,2)=maxjpj
            ENDIF
         ELSE
            READ(val2,*) idxcon(jcon,2)
         ENDIF
C
         line=readnextline(numfila,textexclusion)
         READ(line,'(A,x,A,x,A)') key,val1,val2
         IF (key(1:3).NE.'ini') GOTO 102
         IF (val1(1:3).EQ.'min') THEN
            IF (dir2con(jcon,1).EQ.1)THEN
               mincon(jcon,1)=1-zon_jpi
            ELSE
               mincon(jcon,1)=1-zon_jpj
            ENDIF
         ELSEIF (val1(1:3).EQ.'max') THEN
            IF (dir2con(jcon,1).EQ.1)THEN
               mincon(jcon,1)=maxjpi+zon_jpi-1
            ELSE
               mincon(jcon,1)=maxjpj+zon_jpj-1
            ENDIF
         ELSE
            READ(val1,*) mincon(jcon,1)
         ENDIF
         IF (val2(1:3).EQ.'min') THEN
            IF (dir2con(jcon,2).EQ.1)THEN
               mincon(jcon,2)=1-zon_jpi
            ELSE
               mincon(jcon,2)=1-zon_jpj
            ENDIF
         ELSEIF (val2(1:3).EQ.'max') THEN
            IF (dir2con(jcon,2).EQ.1)THEN
               mincon(jcon,2)=maxjpi+zon_jpi-1
            ELSE
               mincon(jcon,2)=maxjpj+zon_jpj-1
            ENDIF
         ELSE
            READ(val2,*) mincon(jcon,2)
         ENDIF
C
         line=readnextline(numfila,textexclusion)
         READ(line,'(A,x,A,x,A)') key,val1,val2
         IF (key(1:3).NE.'end') GOTO 102
         IF (val1(1:3).EQ.'min') THEN
            IF (dir2con(jcon,1).EQ.1)THEN
               maxcon(jcon,1)=1-zon_jpi
            ELSE
               maxcon(jcon,1)=1-zon_jpj
            ENDIF
         ELSEIF (val1(1:3).EQ.'max') THEN
            IF (dir2con(jcon,1).EQ.1)THEN
               maxcon(jcon,1)=maxjpi+zon_jpi-1
            ELSE
               maxcon(jcon,1)=maxjpj+zon_jpj-1
            ENDIF
         ELSE
            READ(val1,*) mincon(jcon,1)
         ENDIF
         IF (val2(1:3).EQ.'min') THEN
            IF (dir2con(jcon,2).EQ.1)THEN
               maxcon(jcon,2)=1-zon_jpi
            ELSE
               maxcon(jcon,2)=1-zon_jpj
            ENDIF
         ELSEIF (val2(1:3).EQ.'max') THEN
            IF (dir2con(jcon,2).EQ.1)THEN
               maxcon(jcon,2)=maxjpi+zon_jpi-1
            ELSE
               maxcon(jcon,2)=maxjpj+zon_jpj-1
            ENDIF
         ELSE
            READ(val2,*) mincon(jcon,2)
         ENDIF
C
      ENDDO
C
      CLOSE(unit=numfila)
C
C Test coherence of connections
C
      DO jcon=1,jpcon
         DO jline=1,2
            IF (dir1con(jcon,jline).EQ.1) THEN
               IF (idxcon(jcon,jline).LT.1) GOTO 103
               IF (idxcon(jcon,jline).GT.maxjpi) GOTO 103
            ENDIF
            IF (dir1con(jcon,jline).EQ.2) THEN
               IF (idxcon(jcon,jline).LT.1) GOTO 103
               IF (idxcon(jcon,jline).GT.maxjpj) GOTO 103
            ENDIF
         ENDDO
         IF (mincon(jcon,1).GT.maxcon(jcon,1)) GOTO 103
         IF (dir1con(jcon,1).EQ.dir1con(jcon,2)) THEN
            IF (signcon(jcon,1).EQ.signcon(jcon,2)) THEN
               IF ( (maxcon(jcon,1)-mincon(jcon,1)).NE.
     $              (mincon(jcon,2)-maxcon(jcon,2)) ) GOTO 103
            ELSE
               IF ( (maxcon(jcon,1)-mincon(jcon,1)).NE.
     $              (maxcon(jcon,2)-mincon(jcon,2)) ) GOTO 103
            ENDIF
         ELSE
            IF (signcon(jcon,1).EQ.signcon(jcon,2)) THEN
               IF ( (maxcon(jcon,1)-mincon(jcon,1)).NE.
     $              (maxcon(jcon,2)-mincon(jcon,2)) ) GOTO 103
            ELSE
               IF ( (maxcon(jcon,1)-mincon(jcon,1)).NE.
     $              (mincon(jcon,2)-maxcon(jcon,2)) ) GOTO 103
            ENDIF
         ENDIF
      ENDDO
C
C Allocate and initialize Vy pointers
C
C --- allocation pty 
      allocate ( pty(1:maxjpi,1:maxjpj,1:maxjpk,1:maxjpt,1:dtaend), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      pty(:,:,:,:,:)=0
C
      jy=0
      DO jdta=1,dtaend
         inddta=dta_ord(jdta)
         inddtamsk=jdta-1+varend
         DO jt=1,dta_jpt(inddta)
         DO jk=1,dta_jpk(inddta)
            DO jj=1,dta_jpj(inddta)
            DO ji=1,dta_jpi(inddta)
               IF (IBITS(mask(ji,jj,jk,jt),inddtamsk,1).NE.0) THEN
                  jy=jy+1
                  pty(ji,jj,jk,jt,jdta)=jy
               ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
C
C Allocate and initialize connection line pointers
C
      IF (jpcon.GT.0) THEN
C
         jplcon=MAXVAL(maxcon(1:jpcon,1)-mincon(1:jpcon,1))
         jplcon=jplcon*maxjpk*maxjpt*dtaend
C
C --- allocation ptycon
         allocate ( ptycon(1:jplcon,1:jpcon,2), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C --- allocation ptycon
         allocate ( ptyconmask(1:jplcon), stat=allocok )
         IF (allocok.NE.0) GOTO 1001
C
         DO jline=1,2
         DO jcon=1,jpcon
C
            jlcon=0
            DO jdta=1,dtaend
               inddta=dta_ord(jdta)
               DO jt=1,dta_jpt(inddta)
               DO jk=1,dta_jpk(inddta)
                  IF (dir1con(jcon,jline).EQ.1) THEN
                     ji=idxcon(jcon,jline)
                     jjmin=MAX(1,mincon(jcon,jline))
                     jjmax=MIN(maxcon(jcon,jline),maxjpj)
                     DO jj=jjmin,jjmax
                        jlcon=jlcon+1
                        ptycon(jlcon,jcon,jline)=pty(ji,jj,jk,jt,jdta)
                     ENDDO
                  ENDIF
                  IF (dir1con(jcon,jline).EQ.2) THEN
                     jj=idxcon(jcon,jline)
                     jimin=MAX(1,mincon(jcon,jline))
                     jimax=MIN(maxcon(jcon,jline),maxjpi)
                     DO ji=jimin,jimax
                        jlcon=jlcon+1
                        ptycon(jlcon,jcon,jline)=pty(ji,jj,jk,jt,jdta)
                     ENDDO
                  ENDIF
               ENDDO
               ENDDO
            ENDDO
            jlconmax=jlcon
C
            DO jlcon=jlconmax+1,jplcon
               ptycon(jlcon,jcon,jline)=ptycon(jlconmax,jcon,jline)
            ENDDO
C
         ENDDO
         ENDDO
C
      ENDIF
C
C Allocate work arrays for connection management
C
C --- allocation ptbubcon
      allocate ( ptbubcon(1:zon_jpi*zon_jpj*zon_jpk*zon_jpt,0:4), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C --- allocation ptbubmask
      allocate ( ptbubmask(1:zon_jpi*zon_jpj*zon_jpk*zon_jpt), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'mkconnect','connect_init')
 1001 CALL printerror2(0,1001,3,'mkconnect','connect_init')
C
 102  WRITE (texterror,*) 'Bad connection file'
      CALL printerror2(0,102,3,'mkconnect','connect_init',comment=texterror)
 103  WRITE (texterror,*) 'Incoherent connections'
      CALL printerror2(0,102,3,'mkconnect','connect_init',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE connect_close
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
CCC  -------
CCC     Close arrays for management of grid connexions
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
      use mod_mask
      use mod_coord
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
CCC-------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '& routine connect_close &'
         WRITE(numout,*) '&&&&&&&&&&&&&&&&&&&&&&&&&'
         WRITE(numout,*) '&'
         WRITE(numout,*) ' '
      ENDIF
C
      IF (allocated(pty)) deallocate(pty)
      IF (allocated(ptbubcon)) deallocate(ptbubcon)
      IF (allocated(ptbubmask)) deallocate(ptbubmask)
      IF (allocated(dir1con)) deallocate(dir1con)
      IF (allocated(dir2con)) deallocate(dir2con)
      IF (allocated(signcon)) deallocate(signcon)
      IF (allocated(idxcon)) deallocate(idxcon)
      IF (allocated(mincon)) deallocate(mincon)
      IF (allocated(maxcon)) deallocate(maxcon)
      IF (allocated(ptycon)) deallocate(ptycon)
      IF (allocated(ptyconmask)) deallocate(ptyconmask)
C
      RETURN
C
C --- error management section
C
 1000 CALL printerror2(0,1000,1,'mkconnect','connect_close')
 1001 CALL printerror2(0,1001,3,'mkconnect','connect_close')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkconnecty(kvecty)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
CCC  -------
CCC     Identify connected line in Vy object
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
      use mod_coord
      use mod_spacexyo , only : jpyend
      IMPLICIT NONE
C----------------------------------------------------------------------
CC header declarations
CC ==================
      BIGREAL, dimension(0:), intent(inout) :: kvecty
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER :: jcon,jpysize
CC---------------------------------------------------------------------
      jpysize=size(kvecty,1)
      IF (jpysize.NE.jpyend+1) GOTO 1000
C
      DO jcon=1,jpcon
         ptyconmask(:)=kvecty(ptycon(:,jcon,1)).EQ.FREAL(0.0)
         IF (ALL(ptyconmask)) THEN
            kvecty(ptycon(:,jcon,1))=kvecty(ptycon(:,jcon,2))
         ELSE
            kvecty(ptycon(:,jcon,2))=kvecty(ptycon(:,jcon,1))
C        ELSEIF (ALL((vecty(ptycon(:,jcon,2))).EQ.FREAL(0.0))) THEN
C           kvecty(ptycon(:,jcon,2))=kvecty(ptycon(:,jcon,1))
C        ELSE
C           GOTO 102
         ENDIF
      ENDDO
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'mkconnect','mkconnecty')
 1001 CALL printerror2(0,1001,3,'mkconnect','mkconnecty')
C
 102  WRITE (texterror,*) 'Zone overlap both connection lines'
      CALL printerror2(0,102,3,'mkconnect','mkconnecty',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE mkconnect
