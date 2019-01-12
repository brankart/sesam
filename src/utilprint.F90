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
! ---                  UTILPRINT.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 98-06 (C.E. Testut)                        ---
! --- modification : 01-06 (C.E. Testut)                        ---
! --- modification : 07-11 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE printtab_r     : Print r-dimensional array
! --- SUBROUTINE printtab_rr    : Print rr-dimensional array
! --- SUBROUTINE printrapport   : Print eigenvalues and/or eigenvectors
! --- SUBROUTINE printamplitude : Print EOFs time amplitudes
! --- SUBROUTINE printorthog    : Print orthogonality report
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilprint
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC printtab_r,printtab_rr,printrapport,printorthog
      PUBLIC printamplitude

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE printtab_r (tab_r,titre,variable,jrdeb) 
!
!     Print r-dimensional array in 'numout'
!
      use mod_cfgxyo
!
      IMPLICIT NONE
!
      BIGREAL, dimension(:), intent(in) :: tab_r
      CHARACTER(len=*), intent(in) :: titre,variable
      INTEGER, intent(in) :: jrdeb
      INTEGER :: jprsize,jr
!
      jprsize = size(tab_r,1)
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '-------------------------------------------'
         WRITE(numout,*) '==> ',titre(1:lenv(titre))
         WRITE(numout,*) 'array : ',variable(1:lenv(variable))
         WRITE(numout,*) '-------------------------------------------'
         DO jr=1,jprsize
            WRITE(numout,*) variable(1:lenv(variable)) &
     &           ,'(',jr+jrdeb-1,')=',tab_r(jr)
         ENDDO
         WRITE(numout,*) '-------------------------------------------'
      ENDIF
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilprint','printtab_r')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE printtab_rr (tab_rr,titre,variable,jr1deb,jr2deb)
!
!     Print rr-dimensional array in 'numout'
!
      use mod_cfgxyo
!
      IMPLICIT NONE
!
      BIGREAL, dimension(:,:), intent(in) :: tab_rr
      CHARACTER(len=*), intent(in) :: titre,variable
      INTEGER, intent(in) :: jr1deb,jr2deb
      INTEGER :: jpr1size,jpr2size,jr1,jr2
!
      jpr1size = size(tab_rr,1)
      jpr2size = size(tab_rr,2)
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '-------------------------------'
         WRITE(numout,*) '==> ',titre(1:lenv(titre))
         WRITE(numout,*) 'array : ',variable(1:lenv(variable))
         WRITE(numout,*) '-------------------------------'
         DO jr2=1,jpr2size
            WRITE(numout,*) variable(1:lenv(variable)) &
     &           ,'(',jr1deb,'-',jpr1size+jr1deb-1,',',jr2+jr2deb-1, &
     &           ')=',tab_rr(:,jr2)
            WRITE(numout,*)
         ENDDO
         WRITE(numout,*) '-------------------------------'
      ENDIF
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilprint','printtab_rr')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE printrapport (vctp_rr,valp_r,titre,jr2deb,all) 
!
!     Print eigenvalues and/or eigenvectors in 'numout'
!
      use mod_cfgxyo
!
      IMPLICIT NONE
!
      BIGREAL, dimension(:,:), intent(in) :: vctp_rr
      BIGREAL, dimension(:), intent(in) :: valp_r
      CHARACTER(len=*), intent(in) :: titre
      INTEGER, intent(in) :: jr2deb
      LOGICAL, intent(in) :: all
      INTEGER :: jpr1size,jpr2size,jr,jr1,jr2
      BIGREAL :: somvalp
      BIGREAL, dimension(1:size(valp_r,1)) :: ampl,somampl
!
! Get and check size of arrays, initialize computation arrays
      jpr1size = size(vctp_rr,1)
      jpr2size = size(vctp_rr,2)
      IF (jpr2size.NE.size(valp_r,1)) GOTO 1000
      somvalp = FREAL(0.0)
      ampl(:) = FREAL(0.0)
      somampl(:) = FREAL(0.0)
!
! -1.- Compute % of explained variance
! ------------------------------------
! --- compute sum of eigenvalues
      somvalp= SUM(valp_r(:),DIM=1)
! --- compute % of amplitude
      ampl(:) = valp_r(:)*(FREAL(100.0)/somvalp)
! --- compute % of explained variance
      DO jr1=1,jpr2size
         somampl(jr1) = SUM(ampl(:jr1),DIM=1)
      ENDDO
!
! -2.- Print eigenvalues
! ----------------------
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         WRITE(numout,*) '----------------------------------------'
         WRITE(numout,*) '==> eighenvalues :',titre(1:lenv(titre))
         WRITE(numout,*) '----------------------------------------'
         WRITE(numout,*) ' jr | valp | SQRT(valp) | % variance    ', &
     &        '| total of % variance'
         DO jr=1,jpr2size
            IF (valp_r(jr).GE.FREAL(0.0)) THEN
               WRITE(numout,'(i5,3e26.15,f20.13)') jr+jr2deb-1, &
     &              valp_r(jr),SQRT(valp_r(jr)),ampl(jr),somampl(jr)
            ELSE
               WRITE(numout,'(i5,e26.15,a3,e26.15,f20.13)') jr+jr2deb-1, &
     &              valp_r(jr),'nan',ampl(jr),somampl(jr)
            ENDIF
         ENDDO
         WRITE(numout,*) '----------------------------------------'
      ENDIF
!
! -3.- Print eigenvectors
! -----------------------
      IF (all) THEN
      IF ((nprint.GE.0).AND.(jproc.EQ.0)) THEN
         WRITE(numout,*) '-----------------------------------------'
         WRITE(numout,*) '==> eighenvectors :',titre(1:lenv(titre))
         WRITE(numout,*) '-----------------------------------------'
         DO jr2=1,jpr2size
            WRITE(numout,*) 'vctp(*,',jr2+jr2deb-1,')=',vctp_rr(:,jr2)
            WRITE(numout,*)
         ENDDO
         WRITE(numout,*) '----------------------------------------'
      ENDIF
      ENDIF
!
      RETURN
!
! --- error management
!
 1000  CALL printerror2(0,1000,1,'utilprint','printrapport')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE printamplitude (vctp_rr,dirname)
!
!     Print EOFs time amplitudes
!
! --- Module declaration
      use mod_cfgxyo
      use utilvalid , only : fildirbas
      IMPLICIT NONE
! --- Variable declaration 
      BIGREAL, dimension(:,:), intent(in) :: vctp_rr
      CHARACTER(len=*), intent(in) :: dirname
      INTEGER :: jpr1size,jpr2size,jr,jr2,serie,jprout
      CHARACTER(len=bgword) :: amplname
!
      jpr1size = size(vctp_rr,1)
      jpr2size = size(vctp_rr,2)
!
      serie=2
!
      DO jr=1,jpr2size
        CALL fildirbas(amplname,dirname,jprout,jr,serie)
        amplname=dirname(1:lenv(dirname))//'/'//amplname
        OPEN (unit=10,file=amplname)
        DO jr2=1,jpr2size
          WRITE(10,*) vctp_rr(jr2,jr)
        ENDDO
        CLOSE(10)
      ENDDO
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilprint','printamplitude')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE  printorthog (tabscal,titre,jr1deb,jr2deb)
!
!     Print orthogonality report
!
      use mod_cfgxyo
!
      IMPLICIT NONE
!
      BIGREAL, dimension(:,:), intent(in) :: tabscal
      CHARACTER(len=*), intent(in) :: titre
      INTEGER, intent(in) :: jr1deb,jr2deb
      INTEGER :: jpr1size,jpr2size,jr1,jr2
!
      jpr1size = size(tabscal,1)
      jpr2size = size(tabscal,2)
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '--------------------------------------------'
         WRITE(numout,*) 'ORTHOGONALITY (SCALAR PRODUCT) OF THE EOFS:', &
     &        titre(1:lenv(titre))
         WRITE(numout,*) '--------------------------------------------'
         DO jr1=1,jpr1size
            WRITE(numout,*) 'FOR jr1=',jr1+jr1deb-1
            WRITE(numout,*) '-------------'
            WRITE(numout,*)
            DO jr2=jr1,jpr2size
               WRITE(numout,*) 'EOFs jr1=',jr1+jr1deb-1, &
     &            ' AND jr2=',jr2+jr2deb-1 ,' : ',tabscal(jr1,jr2)
            ENDDO
            WRITE(numout,*)
         ENDDO
         WRITE(numout,*) '--------------------------------------------'
      ENDIF
      RETURN
!
! --- error management
!
 1000  CALL printerror2(0,1000,1,'utilprint','printorthog')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilprint
