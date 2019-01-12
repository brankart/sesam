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
! ---                    ALGORTHO.F90                             ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 97-12  ( C.E. Testut)                      ---
! --- modification : 99-05 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE  algortho
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE algortho
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC chkortho

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE chkortho (dirnambas,kflagxyo,kbasesr, &
     &        kjrbasdeb,kjrbasfin,klargxyoweight, &
     &        kargxyoweight,kvectsweight)
!---------------------------------------------------------------------
!
!  Purpose :
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
      use mod_spacexyo , only : jprend
      use hioxyo
      use hiobas
      use utilclc
      use utilprint
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      CHARACTER(len=*), intent(in) :: dirnambas
      BIGREAL, dimension(:,:), intent(inout) :: kbasesr
      INTEGER, intent(in) :: kjrbasdeb,kjrbasfin,kflagxyo
      LOGICAL, intent(in) :: klargxyoweight
      CHARACTER(len=*), optional, intent(in) :: kargxyoweight
      BIGREAL, dimension(:), optional, intent(inout) :: kvectsweight
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:), allocatable :: norm2,norm21
      BIGREAL, dimension(:,:), allocatable :: tabscal,tabscal1
!
      INTEGER :: allocok,jpssize,jprsize
      LOGICAL :: lectinfo,all
      INTEGER :: jr,jnxyo
      CHARACTER(len=bgword) :: textinfo
!----------------------------------------------------------------------
!
      jpssize=size(kbasesr,1)
      jprsize=size(kbasesr,2)
! --- allocation norm2
      allocate ( norm2(1:jprsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      norm2(:) = FREAL(0.0)
! --- allocation tabscal
      allocate ( tabscal(1:jprsize,1:jprsize) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabscal(:,:) = FREAL(0.0)
! --- allocation tabscal1
      allocate ( tabscal1(1:jprsize,1:jprsize) , stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      tabscal1(:,:) = FREAL(0.0)
!---------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '      *************************'
         WRITE(numout,*) '      *     routine algortho  *'
         WRITE(numout,*) '      *************************'
         WRITE(numout,*) ' '
      ENDIF
!
! -0.- Initialisation :
! ---------------------
!
      IF ((.NOT.(present(kvectsweight))).AND.(klargxyoweight)) GOTO 1000
      IF (klargxyoweight) THEN 
         IF (jpssize.NE.size(kvectsweight,1)) GOTO 101
      ENDIF
      IF (jprsize.NE.jprend) GOTO 101
      lectinfo=.FALSE.
!
! -1.- Assess the orthogonality :
! -------------------------------
!
! -1.1- Orthogonality computation :
! -----------------------_---------
!
      DO jnxyo=1,limjpnxyo(kflagxyo)
         IF (limjpnxyo(kflagxyo).NE.1) THEN
            CALL readbas(dirnambas,kbasesr(:,:),jnxyo, &
     &          kjrbasdeb,kjrbasfin,lectinfo,kflagxyo)
            IF (klargxyoweight) CALL readvar(kargxyoweight, &
     &           kvectsweight(:),jnxyo,lectinfo,kflagxyo)
         ENDIF
         IF (klargxyoweight) THEN
            DO jr=kjrbasdeb,kjrbasfin
               CALL mkcovweight(tabscal1(jr:kjrbasfin,jr:jr), &
     &              kbasesr(:,jr:kjrbasfin),kvectsweight(:))
            ENDDO
         ELSE
            DO jr=kjrbasdeb,kjrbasfin
               CALL mkcov (tabscal1(jr:kjrbasfin,jr:jr), &
     &              kbasesr(:,jr:kjrbasfin))
            ENDDO
         ENDIF
         DO jr=kjrbasdeb,kjrbasfin
            tabscal(jr:kjrbasfin,jr) = tabscal(jr:kjrbasfin,jr) &
     &           + tabscal1(jr:kjrbasfin,jr)
         ENDDO
      ENDDO
!
      norm2(kjrbasdeb:kjrbasfin) =  &
     &     (/ (tabscal(jr,jr),jr=kjrbasdeb,kjrbasfin) /)
!
      DO jr=kjrbasdeb,kjrbasfin
         tabscal(jr:kjrbasfin,jr)=tabscal(jr:kjrbasfin,jr) &
     &           / (SQRT(norm2(jr))*SQRT(norm2(jr:kjrbasfin)))
         tabscal(jr,jr:kjrbasfin)=tabscal(jr:kjrbasfin,jr)
      ENDDO
!
! -1.2-  Print control :
! ----------------------
!
      all=.FALSE.
      WRITE (textinfo,'(a,i1.1,a)') &
     &     'FROM NORM OF BASIS (METHOD ',kflagxyo,')'
      CALL printrapport(tabscal(:,kjrbasdeb:kjrbasfin) &
     &     ,norm2(kjrbasdeb:kjrbasfin),textinfo,kjrbasdeb,all) 
!
! -1.3- Print control orthogonality :
! -----------------------------------
!
      WRITE (textinfo,'(a,i1.1,a)') &
     &     'FROM BASIS WITH (RE)ORTHOGONALISATION (METHOD ',kflagxyo,')'
      CALL printorthog(tabscal(kjrbasdeb:kjrbasfin,kjrbasdeb:kjrbasfin), &
     &     textinfo,kjrbasdeb,kjrbasdeb)
!
! --- deallocation
      IF (allocated(norm2)) deallocate (norm2)
      IF (allocated(tabscal)) deallocate (tabscal)
      IF (allocated(tabscal1)) deallocate (tabscal1)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'algortho','algortho')
 1001 CALL printerror2(0,1001,3,'algortho','algortho')
!
 101  WRITE (texterror,*) 'arguments are false'
      CALL printerror2(0,101,3,'algortho','algortho',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE algortho
