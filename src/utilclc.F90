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
! ---                    UTILCLC.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 98-06 (C.E. Testut)                        ---
! --- revised      : 99-05 (C.E. Testut)                        ---
! --- revised      : 01-06 (C.E. Testut)                        ---
! --- revised      : 03-04 (J.M. Brankart)                      ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! ---                                                
! --- SUBROUTINE mkcov               : Compute C = transp(S) S
! --- SUBROUTINE mkcovweight         : Compute C = transp(S) W W S
! --- SUBROUTINE mkmean              : Compute Sm = mean of S columns
! --- SUBROUTINE prodmat_zr_rr_vias  : Compute S(z*r)=S(z*r)*C(r*r)
! ---                                  using TMP(s*r) as working space
! --- SUBROUTINE prodmat_zr_rd_vias  : Compute S(z*d)=S(z*r)*C(r*d)
! ---                                  using TMP(s*r) as working space
! --- SUBROUTINE prodmat_wr_rrz_vias : Compute S(w*r)=S(w*r)*C(r*r,z)
! ---                                  using TMP(s*r) as working space
! ---                                  (use C corresponding to local
! ---                                  subsystem z in scalar products)
! --- SUBROUTINE prodmat_wr_rdz_vias : Compute S(w*d)=S(w*r)*C(r*d,z)
! ---                                  using TMP(s*r) as working space
! ---                                  (use C corresponding to local
! ---                                  subsystem z in scalar products)
! --- SUBROUTINE invmat_rr           : Compute inverse of C(r*r)
! --- 
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilclc
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC mkcov, mkcovweight, mkmean
      PUBLIC prodmat_zr_rr_vias,prodmat_zr_rd_vias
      PUBLIC prodmat_wr_rrz_vias,prodmat_wr_rdz_vias,invmat_rr

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkcov(cov,matsr) 
!
! Compute C = transp(S) S
!
! Arguments:  cov   : C(r*r) (output)
!             matsr : S(s*r) (input)
!
! --- Module declaration
      IMPLICIT NONE
! --- Variable declaration
      BIGREAL, dimension(:,:), intent(out) :: cov
      BIGREAL, dimension(:,:), intent(in) :: matsr
      INTEGER :: jpr1size,jpr2size,jr1,jr2
!
      jpr1size = size(cov,1)
      jpr2size = size(cov,2)
      IF (jpr1size.GT.size(matsr,2)) GOTO 1000
      IF (jpr2size.GT.size(matsr,2)) GOTO 1000
      cov(:,:) = FREAL(0.0)
! 
      DO jr1=1,jpr1size
      DO jr2=1,jpr2size
         cov(jr1,jr2)=DOT_PRODUCT(matsr(:,jr1),matsr(:,jr2))
      ENDDO
      ENDDO
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilclc','mkcov')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkcovweight (covweight,matsr,vectsweight) 
!
!
! Compute C = transp(S) W W S
!
! Arguments:  covweight   : C(r*r) (output)
!             matsr       : S(s*r) (input)
!             vectsweight : W(s*s) diagonal matrix (input)
!
! --- Module declaration
      IMPLICIT NONE
! --- Variable declaration
      BIGREAL, dimension(:,:), intent(out) :: covweight
      BIGREAL, dimension(:,:), intent(in) :: matsr
      BIGREAL, dimension(:), intent(in) :: vectsweight
      INTEGER :: jpssize,jpr1size,jpr2size,jr1,jr2
!
      jpssize = size(matsr,1)
      jpr1size = size(covweight,1)
      jpr2size = size(covweight,2)
      IF (jpr1size.GT.size(matsr,2)) GOTO 1000
      IF (jpr2size.GT.size(matsr,2)) GOTO 1000
      IF (jpssize.NE.size(vectsweight,1)) GOTO 1000
      covweight(:,:) = FREAL(0.0)
!
      DO jr1=1,jpr1size
      DO jr2=1,jpr2size
         covweight(jr1,jr2)=DOT_PRODUCT(matsr(:,jr1)*vectsweight(:), &
     &        matsr(:,jr2)*vectsweight(:))
      ENDDO
      ENDDO
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilclc','mkcovweight')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkmean (mean_s,mat_sr)
!
! Compute Sm = mean of S columns
!
! Arguments:  covweight   : Sm(s) (output)
!             matsr       : S(s*r) (input)
!
! --- Module declaration
      IMPLICIT NONE
! --- Variable declaration
      BIGREAL, dimension(:), intent(out) :: mean_s
      BIGREAL, dimension(:,:), intent(in) :: mat_sr
      INTEGER :: jprsize,jr
!
      jprsize = size(mat_sr,2)
      IF (size(mean_s,1).NE.size(mat_sr,1)) GOTO 1000
      mean_s(:) = FREAL(0.0)
!
      DO jr=1,jprsize
         mean_s(:) = mean_s(:) + mat_sr(:,jr)
      ENDDO
      mean_s(:) = mean_s(:) / jprsize
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilclc','mkmean')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE prodmat_zr_rr_vias(mat_zr,mat_sr,mat_rr, &
     &     jrdebproj,jrfinproj)
!
! Compute S(z*r)=S(z*r)*C(r*r) using TMP(s*r) as working space
!
! Arguments:  mat_zr : S(z*r) (input and output)
!             mat_sr : TMP(s*r) (output)
!             mat_rr : C(r*r) (input)
!             jrdebproj : first index of S to compute
!             jrfinproj : last index of S to compute
!
! --- Module declaration
      IMPLICIT NONE
! --- Variable declaration   
      BIGREAL, dimension(:,:), intent(inout) :: mat_zr
      BIGREAL, dimension(:,:), intent(out) :: mat_sr
      BIGREAL, dimension(:,:), intent(in) :: mat_rr
      INTEGER, intent(in) :: jrdebproj,jrfinproj
      INTEGER :: jpzsize,jpssize,jprsize,jzdeb,jzfin,jsfin,jr1,jr2,jr
!
      jpzsize = size(mat_zr,1)
      jpssize = size(mat_sr,1)
      jprsize = size(mat_rr,1)
      IF (jprsize.NE.size(mat_zr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_sr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_rr,2)) GOTO 1000
      IF (jprsize.LT.jrfinproj) GOTO 1000
!
      mat_sr(:,:) = FREAL(0.0)
!
      DO jzdeb=1,jpzsize,jpssize
         jzfin=MIN(jpssize+jzdeb-1,jpzsize)
         jsfin=MIN(jpssize,jzfin-jzdeb+1)
         DO jr=1,jprsize
            mat_sr(1:jsfin,jr)=mat_zr(jzdeb:jzfin,jr)
         ENDDO
         mat_zr(jzdeb:jzfin,:) = FREAL(0.0)
         DO jr1=jrdebproj,jrfinproj
            DO jr2=1,jprsize
               mat_zr(jzdeb:jzfin,jr1)=mat_zr(jzdeb:jzfin,jr1) + &
     &                mat_sr(1:jsfin,jr2)*mat_rr(jr2,jr1)       
            ENDDO
         ENDDO
      ENDDO
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilclc','prodmat_zr_rr_vias')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE prodmat_zr_rd_vias(mat_zr,mat_sr,mat_rd, &
     &     jddebproj,jdfinproj)
!
! Compute S(z*d)=S(z*r)*C(r*d) using TMP(s*r) as working space
!
! Arguments:  mat_zr : S(z*r) (input and output), d <= r
!             mat_sr : TMP(s*r) (output)
!             mat_rr : C(r*d) (input)
!             jddebproj : first index of S to compute
!             jdfinproj : last index of S to compute
!
! --- Module declaration
      IMPLICIT NONE
! --- Variable declaration    
      BIGREAL, dimension(:,:), intent(inout) :: mat_zr
      BIGREAL, dimension(:,:), intent(out) :: mat_sr
      BIGREAL, dimension(:,:), intent(in) :: mat_rd
      INTEGER, intent(in) :: jddebproj,jdfinproj
      INTEGER :: jpzsize,jpssize,jprsize,jpdsize, &
     &     jzdeb,jzfin,jsfin,jd,jr
!
      jpzsize = size(mat_zr,1)
      jpssize = size(mat_sr,1)
      jprsize = size(mat_rd,1)
      jpdsize = size(mat_rd,2)
      IF (jprsize.NE.size(mat_zr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_sr,2)) GOTO 1000
      IF (jpdsize.LT.jdfinproj) GOTO 1000
      IF (jpdsize.GT.jprsize) GOTO 1000
!
      mat_sr(:,:) = FREAL(0.0)
!
      DO jzdeb=1,jpzsize,jpssize
         jzfin=MIN(jpssize+jzdeb-1,jpzsize)
         jsfin=MIN(jpssize,jzfin-jzdeb+1)
         DO jr=1,jprsize
            mat_sr(1:jsfin,jr)=mat_zr(jzdeb:jzfin,jr)
         ENDDO
         mat_zr(jzdeb:jzfin,:) = FREAL(0.0)
         DO jd=jddebproj,jdfinproj
            DO jr=1,jprsize
               mat_zr(jzdeb:jzfin,jd)=mat_zr(jzdeb:jzfin,jd) + &
     &                mat_sr(1:jsfin,jr)*mat_rd(jr,jd)       
            ENDDO
         ENDDO
      ENDDO
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilclc','prodmat_zr_rd_vias')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE invmat_rr (mat_rr)
!
! Compute inverse of a symmetric matrix C using Choleski decomposition
!
! Argument:  mat_rr : C (input and output)
!
! --- Module declaration
      use mod_cfgxyo
      IMPLICIT NONE
! --- Variable declaration     
      BIGREAL, dimension(:,:), intent(inout) :: mat_rr
      BIGREAL, dimension(:,:), allocatable :: invl
      BIGREAL, dimension(:), allocatable :: p
      BIGREAL :: sum
      INTEGER :: jprsize,jr1,jr2,jk
      INTEGER :: allocok
!
      jprsize = size(mat_rr,1)
      IF (jprsize.NE.size(mat_rr,2)) GOTO 1000
!
      allocate(invl(jprsize,jprsize),stat=allocok)
      IF (allocok.NE.0) GOTO 1001
      allocate(p(jprsize),stat=allocok)
      IF (allocok.NE.0) GOTO 1001
!
! -1.- Choleski decomposision
! ---------------------------
! Compute L such that: C = L L^T
!
      DO jr1=1,jprsize
      DO jr2=1,jprsize
        sum = mat_rr(jr1,jr2)
        DO jk=jr1-1,1,-1
          sum = sum - mat_rr(jr1,jk) * mat_rr(jr2,jk)
        ENDDO
        IF (jr1.EQ.jr2) THEN
          IF (sum.LE.0.) GOTO 101
          p(jr1) = SQRT(sum)
        ELSE
          mat_rr(jr2,jr1)=sum/p(jr1)
        ENDIF
      ENDDO
      ENDDO
!
      DO jr2=1,jprsize
         mat_rr(1:jr2-1,jr2) = FREAL(0.0)
      ENDDO
!
! -2.- Inverse computation
! ------------------------
! Inverse L
!
      DO jr1=1,jprsize
        mat_rr(jr1,jr1)=1.0/p(jr1)
        DO jr2=jr1+1,jprsize
          sum = 0.0
          DO jk=jr1,jr2-1
            sum = sum - mat_rr(jr2,jk) * mat_rr(jk,jr1)
          ENDDO
          mat_rr(jr2,jr1) = sum / p(jr2)
        ENDDO
      ENDDO
      invl = mat_rr
!
! Inverse A
!
      DO jr1=1,jprsize
      DO jr2=1,jprsize
        mat_rr(jr1,jr2) = DOT_PRODUCT(invl(:,jr1),invl(:,jr2))
      ENDDO
      ENDDO
!
      IF (allocated(invl)) deallocate(invl)
      IF (allocated(p)) deallocate(p)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilclc','invmat_rr')
 1001 CALL printerror2(0,1000,3,'utilclc','invmat_rr')
!
 101  WRITE (texterror,*) 'Error in Choleski decomposition'
      CALL printerror2(0,101,3,'utilclc','invmat_rr',comment=texterror)
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE prodmat_wr_rrz_vias(mat_wr,mat_sr,mat_rrz,vectx, &
     &     jrdebproj,jrfinproj)
!
! Compute S(w*r)=S(w*r)*C(r*r,z) using TMP(s*r) as working space
! (use C corresponding to local subsystem z in scalar products)
!
! Arguments:  mat_wr : S(w*r) (input and output)
!             mat_sr : TMP(s*r) (output)
!             mat_rr : C(r*r) (input)
!             jddebproj : first index of S to compute
!             jdfinproj : last index of S to compute
!             vectx  : partition vector containing subsystem indices
!
! --- Module declaration
      use mod_cfgxyo
      IMPLICIT NONE
! --- Variable declaration    
      BIGREAL, dimension(:,:), intent(inout) :: mat_wr
      BIGREAL, dimension(:,:), intent(out) :: mat_sr
      BIGREAL, dimension(:,:,0:), intent(in) :: mat_rrz
      BIGREAL, dimension(:), intent(in) :: vectx
      INTEGER, intent(in) :: jrdebproj,jrfinproj
      INTEGER :: jpzsize,jpssize,jprsize,jpwsize
      INTEGER :: jwdeb,jwfin,jsfin,jr1,jr2,jr,js
!
      jpwsize = size(mat_wr,1)
      jpssize = size(mat_sr,1)
      jprsize = size(mat_rrz,1)
      jpzsize = size(mat_rrz,3)
      IF (jprsize.NE.size(mat_wr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_sr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_rrz,2)) GOTO 1000
      IF (jpwsize.NE.size(vectx,1)) GOTO 1000
      IF (jprsize.LT.jrfinproj) GOTO 1000
!
      mat_sr(:,:) = FREAL(0.0)
!
      DO jwdeb=1,jpwsize,jpssize
         jwfin=MIN(jpssize+jwdeb-1,jpwsize)
         jsfin=MIN(jpssize,jwfin-jwdeb+1)
         DO jr=1,jprsize
            mat_sr(1:jsfin,jr)=mat_wr(jwdeb:jwfin,jr)
         ENDDO
         mat_wr(jwdeb:jwfin,:) = FREAL(0.0)
         DO jr1=jrdebproj,jrfinproj
            DO jr2=1,jprsize
               DO js = 1,jsfin
                  mat_wr(jwdeb+js-1,jr1)=mat_wr(jwdeb+js-1,jr1) + &
     &                mat_sr(js,jr2) * &
     &                mat_rrz(jr2,jr1,nint(vectx(jwdeb+js-1)))
               ENDDO
            ENDDO
         ENDDO
 
      ENDDO
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilclc','prodmat_wr_rrz_vias')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE prodmat_wr_rdz_vias(mat_wr,mat_sr,mat_rdz,vectx, &
     &     jddebproj,jdfinproj)
!
! Compute S(w*d)=S(w*r)*C(r*d,z) using TMP(s*r) as working space
! (use C corresponding to local subsystem z in scalar products)
!
! Arguments:  mat_wr : S(w*r) (input and output), d <= r
!             mat_sr : TMP(s*r) (output)
!             mat_rr : C(r*d) (input)
!             jddebproj : first index of S to compute
!             jdfinproj : last index of S to compute
!             vectx  : partition vector containing subsystem indices
!
! --- Module declaration
      IMPLICIT NONE
! --- Variable declaration    
      BIGREAL, dimension(:,:), intent(inout) :: mat_wr
      BIGREAL, dimension(:,:), intent(out) :: mat_sr
      BIGREAL, dimension(:,:,0:), intent(in) :: mat_rdz
      BIGREAL, dimension(:), intent(in) :: vectx
      INTEGER, intent(in) :: jddebproj,jdfinproj
      INTEGER :: jpzsize,jpssize,jprsize,jpdsize,jpwsize
      INTEGER :: jwdeb,jwfin,jsfin,jr,jd
!
      jpwsize = size(mat_wr,1)
      jpssize = size(mat_sr,1)
      jprsize = size(mat_rdz,1)
      jpdsize = size(mat_rdz,2)
      jpzsize = size(mat_rdz,3)
      IF (jprsize.NE.size(mat_wr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_sr,2)) GOTO 1000
      IF (jpwsize.NE.size(vectx,1)) GOTO 1000
      IF (jpdsize.LT.jdfinproj) GOTO 1000
      IF (jpdsize.GT.jprsize) GOTO 1000
!
      mat_sr(:,:) = FREAL(0.0)
!
      DO jwdeb=1,jpwsize,jpssize
         jwfin=MIN(jpssize+jwdeb-1,jpwsize)
         jsfin=MIN(jpssize,jwfin-jwdeb+1)
         DO jr=1,jprsize
            mat_sr(1:jsfin,jr)=mat_wr(jwdeb:jwfin,jr)
         ENDDO
         mat_wr(jwdeb:jwfin,:) = FREAL(0.0)
         DO jd=jddebproj,jdfinproj
            DO jr=1,jprsize
               mat_wr(jwdeb:jwfin,jd)=mat_wr(jwdeb:jwfin,jd) + &
     &                mat_sr(1:jsfin,jr) *  &
     &                mat_rdz(jr,jd,nint(vectx(jwdeb:jwfin)))
            ENDDO
         ENDDO
 
      ENDDO
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilclc','prodmat_wr_rdz_vias')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilclc
