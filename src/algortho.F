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
C ---                    ALGORTHO.F                             ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 97-12  ( C.E. Testut)                      ---
C --- modification : 99-05 (C.E. Testut)                        ---
C --- modification : 01-06 (C.E. Testut)                        ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- 
C --- SUBROUTINE  algortho
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algortho
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC chkortho

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE chkortho (dirnambas,kflagxyo,kbasesr,
     $        kjrbasdeb,kjrbasfin,klargxyoweight,
     $        kargxyoweight,kvectsweight)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose :
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
      use mod_spacexyo , only : jprend
      use hioxyo
      use hiobas
      use utilclc
      use utilprint
      IMPLICIT NONE
CC----------------------------------------------------------------------
CC header declarations
CC ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      BIGREAL, dimension(:,:), intent(inout) :: kbasesr
      INTEGER, intent(in) :: kjrbasdeb,kjrbasfin,kflagxyo
      LOGICAL, intent(in) :: klargxyoweight
      CHARACTER(len=*), optional, intent(in) :: kargxyoweight
      BIGREAL, dimension(:), optional, intent(inout) :: kvectsweight
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      BIGREAL, dimension(:), allocatable :: norm2,norm21
      BIGREAL, dimension(:,:), allocatable :: tabscal,tabscal1
C
      INTEGER :: allocok,jpssize,jprsize
      LOGICAL :: lectinfo,all
      INTEGER :: jr,jnxyo
      CHARACTER(len=bgword) :: textinfo
CC----------------------------------------------------------------------
C
      jpssize=size(kbasesr,1)
      jprsize=size(kbasesr,2)
C --- allocation norm2
      allocate ( norm2(1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      norm2(:) = FREAL(0.0)
C --- allocation tabscal
      allocate ( tabscal(1:jprsize,1:jprsize) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabscal(:,:) = FREAL(0.0)
C --- allocation tabscal1
      allocate ( tabscal1(1:jprsize,1:jprsize) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabscal1(:,:) = FREAL(0.0)
CCC---------------------------------------------------------------------
C
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '      *************************'
         WRITE(numout,*) '      *     routine algortho  *'
         WRITE(numout,*) '      *************************'
         WRITE(numout,*) ' '
      ENDIF
C
C -0.- Initialisation :
C ---------------------
C
      IF ((.NOT.(present(kvectsweight))).AND.(klargxyoweight)) GOTO 1000
      IF (klargxyoweight) THEN 
         IF (jpssize.NE.size(kvectsweight,1)) GOTO 101
      ENDIF
      IF (jprsize.NE.jprend) GOTO 101
      lectinfo=.FALSE.
C
C -1.- Assess the orthogonality :
C -------------------------------
C
C -1.1- Orthogonality computation :
C -----------------------_---------
C
      DO jnxyo=1,limjpnxyo(kflagxyo)
         IF (limjpnxyo(kflagxyo).NE.1) THEN
            CALL readbas(dirnambas,kbasesr(:,:),jnxyo,
     $          kjrbasdeb,kjrbasfin,lectinfo,kflagxyo)
            IF (klargxyoweight) CALL readvar(kargxyoweight,
     $           kvectsweight(:),jnxyo,lectinfo,kflagxyo)
         ENDIF
         IF (klargxyoweight) THEN
            DO jr=kjrbasdeb,kjrbasfin
               CALL mkcovweight(tabscal1(jr:kjrbasfin,jr:jr),
     $              kbasesr(:,jr:kjrbasfin),kvectsweight(:))
            ENDDO
         ELSE
            DO jr=kjrbasdeb,kjrbasfin
               CALL mkcov (tabscal1(jr:kjrbasfin,jr:jr),
     $              kbasesr(:,jr:kjrbasfin))
            ENDDO
         ENDIF
         DO jr=kjrbasdeb,kjrbasfin
            tabscal(jr:kjrbasfin,jr) = tabscal(jr:kjrbasfin,jr)
     $           + tabscal1(jr:kjrbasfin,jr)
         ENDDO
      ENDDO
C
      norm2(kjrbasdeb:kjrbasfin) = 
     $     (/ (tabscal(jr,jr),jr=kjrbasdeb,kjrbasfin) /)
C
      DO jr=kjrbasdeb,kjrbasfin
         tabscal(jr:kjrbasfin,jr)=tabscal(jr:kjrbasfin,jr)
     $           / (SQRT(norm2(jr))*SQRT(norm2(jr:kjrbasfin)))
         tabscal(jr,jr:kjrbasfin)=tabscal(jr:kjrbasfin,jr)
      ENDDO
C
C -1.2-  Print control :
C ----------------------
C
      all=.FALSE.
      WRITE (textinfo,'(a,i1.1,a)')
     $     'FROM NORM OF BASIS (METHOD ',kflagxyo,')'
      CALL printrapport(tabscal(:,kjrbasdeb:kjrbasfin)
     $     ,norm2(kjrbasdeb:kjrbasfin),textinfo,kjrbasdeb,all) 
C
C -1.3- Print control orthogonality :
C -----------------------------------
C
      WRITE (textinfo,'(a,i1.1,a)')
     $     'FROM BASIS WITH (RE)ORTHOGONALISATION (METHOD ',kflagxyo,')'
      CALL printorthog(tabscal(kjrbasdeb:kjrbasfin,kjrbasdeb:kjrbasfin),
     $     textinfo,kjrbasdeb,kjrbasdeb)
C
C --- deallocation
      IF (allocated(norm2)) deallocate (norm2)
      IF (allocated(tabscal)) deallocate (tabscal)
      IF (allocated(tabscal1)) deallocate (tabscal1)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'algortho','algortho')
 1001 CALL printerror2(0,1001,3,'algortho','algortho')
C
 101  WRITE (texterror,*) 'arguments are false'
      CALL printerror2(0,101,3,'algortho','algortho',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algortho
