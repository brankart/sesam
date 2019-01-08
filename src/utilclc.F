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
C ---                    UTILCLC.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 98-06 (C.E. Testut)                        ---
C --- revised      : 99-05 (C.E. Testut)                        ---
C --- revised      : 01-06 (C.E. Testut)                        ---
C --- revised      : 03-04 (J.M. Brankart)                      ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C ---                                                
C --- SUBROUTINE mkcov               : Compute C = transp(S) S
C --- SUBROUTINE mkcovweight         : Compute C = transp(S) W W S
C --- SUBROUTINE mkmean              : Compute Sm = mean of S columns
C --- SUBROUTINE prodmat_zr_rr_vias  : Compute S(z*r)=S(z*r)*C(r*r)
C ---                                  using TMP(s*r) as working space
C --- SUBROUTINE prodmat_zr_rd_vias  : Compute S(z*d)=S(z*r)*C(r*d)
C ---                                  using TMP(s*r) as working space
C --- SUBROUTINE prodmat_wr_rrz_vias : Compute S(w*r)=S(w*r)*C(r*r,z)
C ---                                  using TMP(s*r) as working space
C ---                                  (use C corresponding to local
C ---                                  subsystem z in scalar products)
C --- SUBROUTINE prodmat_wr_rdz_vias : Compute S(w*d)=S(w*r)*C(r*d,z)
C ---                                  using TMP(s*r) as working space
C ---                                  (use C corresponding to local
C ---                                  subsystem z in scalar products)
C --- SUBROUTINE invmat_rr           : Compute inverse of C(r*r)
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilclc
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC mkcov, mkcovweight, mkmean
      PUBLIC prodmat_zr_rr_vias,prodmat_zr_rd_vias
      PUBLIC prodmat_wr_rrz_vias,prodmat_wr_rdz_vias,invmat_rr

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkcov(cov,matsr) 
C
C Compute C = transp(S) S
C
C Arguments:  cov   : C(r*r) (output)
C             matsr : S(s*r) (input)
C
C --- Module declaration
      IMPLICIT NONE
C --- Variable declaration
      BIGREAL, dimension(:,:), intent(out) :: cov
      BIGREAL, dimension(:,:), intent(in) :: matsr
      INTEGER :: jpr1size,jpr2size,jr1,jr2
C
      jpr1size = size(cov,1)
      jpr2size = size(cov,2)
      IF (jpr1size.GT.size(matsr,2)) GOTO 1000
      IF (jpr2size.GT.size(matsr,2)) GOTO 1000
      cov(:,:) = FREAL(0.0)
C 
      DO jr1=1,jpr1size
      DO jr2=1,jpr2size
         cov(jr1,jr2)=DOT_PRODUCT(matsr(:,jr1),matsr(:,jr2))
      ENDDO
      ENDDO
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilclc','mkcov')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkcovweight (covweight,matsr,vectsweight) 
C
C
C Compute C = transp(S) W W S
C
C Arguments:  covweight   : C(r*r) (output)
C             matsr       : S(s*r) (input)
C             vectsweight : W(s*s) diagonal matrix (input)
C
C --- Module declaration
      IMPLICIT NONE
C --- Variable declaration
      BIGREAL, dimension(:,:), intent(out) :: covweight
      BIGREAL, dimension(:,:), intent(in) :: matsr
      BIGREAL, dimension(:), intent(in) :: vectsweight
      INTEGER :: jpssize,jpr1size,jpr2size,jr1,jr2
C
      jpssize = size(matsr,1)
      jpr1size = size(covweight,1)
      jpr2size = size(covweight,2)
      IF (jpr1size.GT.size(matsr,2)) GOTO 1000
      IF (jpr2size.GT.size(matsr,2)) GOTO 1000
      IF (jpssize.NE.size(vectsweight,1)) GOTO 1000
      covweight(:,:) = FREAL(0.0)
C
      DO jr1=1,jpr1size
      DO jr2=1,jpr2size
         covweight(jr1,jr2)=DOT_PRODUCT(matsr(:,jr1)*vectsweight(:),
     $        matsr(:,jr2)*vectsweight(:))
      ENDDO
      ENDDO
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilclc','mkcovweight')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE mkmean (mean_s,mat_sr)
C
C Compute Sm = mean of S columns
C
C Arguments:  covweight   : Sm(s) (output)
C             matsr       : S(s*r) (input)
C
C --- Module declaration
      IMPLICIT NONE
C --- Variable declaration
      BIGREAL, dimension(:), intent(out) :: mean_s
      BIGREAL, dimension(:,:), intent(in) :: mat_sr
      INTEGER :: jprsize,jr
C
      jprsize = size(mat_sr,2)
      IF (size(mean_s,1).NE.size(mat_sr,1)) GOTO 1000
      mean_s(:) = FREAL(0.0)
C
      DO jr=1,jprsize
         mean_s(:) = mean_s(:) + mat_sr(:,jr)
      ENDDO
      mean_s(:) = mean_s(:) / jprsize
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilclc','mkmean')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE prodmat_zr_rr_vias(mat_zr,mat_sr,mat_rr,
     $     jrdebproj,jrfinproj)
C
C Compute S(z*r)=S(z*r)*C(r*r) using TMP(s*r) as working space
C
C Arguments:  mat_zr : S(z*r) (input and output)
C             mat_sr : TMP(s*r) (output)
C             mat_rr : C(r*r) (input)
C             jrdebproj : first index of S to compute
C             jrfinproj : last index of S to compute
C
C --- Module declaration
      IMPLICIT NONE
C --- Variable declaration   
      BIGREAL, dimension(:,:), intent(inout) :: mat_zr
      BIGREAL, dimension(:,:), intent(out) :: mat_sr
      BIGREAL, dimension(:,:), intent(in) :: mat_rr
      INTEGER, intent(in) :: jrdebproj,jrfinproj
      INTEGER :: jpzsize,jpssize,jprsize,jzdeb,jzfin,jsfin,jr1,jr2,jr
C
      jpzsize = size(mat_zr,1)
      jpssize = size(mat_sr,1)
      jprsize = size(mat_rr,1)
      IF (jprsize.NE.size(mat_zr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_sr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_rr,2)) GOTO 1000
      IF (jprsize.LT.jrfinproj) GOTO 1000
C
      mat_sr(:,:) = FREAL(0.0)
C
      DO jzdeb=1,jpzsize,jpssize
         jzfin=MIN(jpssize+jzdeb-1,jpzsize)
         jsfin=MIN(jpssize,jzfin-jzdeb+1)
         DO jr=1,jprsize
            mat_sr(1:jsfin,jr)=mat_zr(jzdeb:jzfin,jr)
         ENDDO
         mat_zr(jzdeb:jzfin,:) = FREAL(0.0)
         DO jr1=jrdebproj,jrfinproj
            DO jr2=1,jprsize
               mat_zr(jzdeb:jzfin,jr1)=mat_zr(jzdeb:jzfin,jr1) +
     $                mat_sr(1:jsfin,jr2)*mat_rr(jr2,jr1)       
            ENDDO
         ENDDO
      ENDDO
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilclc','prodmat_zr_rr_vias')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE prodmat_zr_rd_vias(mat_zr,mat_sr,mat_rd,
     $     jddebproj,jdfinproj)
C
C Compute S(z*d)=S(z*r)*C(r*d) using TMP(s*r) as working space
C
C Arguments:  mat_zr : S(z*r) (input and output), d <= r
C             mat_sr : TMP(s*r) (output)
C             mat_rr : C(r*d) (input)
C             jddebproj : first index of S to compute
C             jdfinproj : last index of S to compute
C
C --- Module declaration
      IMPLICIT NONE
C --- Variable declaration    
      BIGREAL, dimension(:,:), intent(inout) :: mat_zr
      BIGREAL, dimension(:,:), intent(out) :: mat_sr
      BIGREAL, dimension(:,:), intent(in) :: mat_rd
      INTEGER, intent(in) :: jddebproj,jdfinproj
      INTEGER :: jpzsize,jpssize,jprsize,jpdsize,
     $     jzdeb,jzfin,jsfin,jd,jr
C
      jpzsize = size(mat_zr,1)
      jpssize = size(mat_sr,1)
      jprsize = size(mat_rd,1)
      jpdsize = size(mat_rd,2)
      IF (jprsize.NE.size(mat_zr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_sr,2)) GOTO 1000
      IF (jpdsize.LT.jdfinproj) GOTO 1000
      IF (jpdsize.GT.jprsize) GOTO 1000
C
      mat_sr(:,:) = FREAL(0.0)
C
      DO jzdeb=1,jpzsize,jpssize
         jzfin=MIN(jpssize+jzdeb-1,jpzsize)
         jsfin=MIN(jpssize,jzfin-jzdeb+1)
         DO jr=1,jprsize
            mat_sr(1:jsfin,jr)=mat_zr(jzdeb:jzfin,jr)
         ENDDO
         mat_zr(jzdeb:jzfin,:) = FREAL(0.0)
         DO jd=jddebproj,jdfinproj
            DO jr=1,jprsize
               mat_zr(jzdeb:jzfin,jd)=mat_zr(jzdeb:jzfin,jd) +
     $                mat_sr(1:jsfin,jr)*mat_rd(jr,jd)       
            ENDDO
         ENDDO
      ENDDO
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilclc','prodmat_zr_rd_vias')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE invmat_rr (mat_rr)
C
C Compute inverse of a symmetric matrix C using Choleski decomposition
C
C Argument:  mat_rr : C (input and output)
C
C --- Module declaration
      use mod_cfgxyo
      IMPLICIT NONE
C --- Variable declaration     
      BIGREAL, dimension(:,:), intent(inout) :: mat_rr
      BIGREAL, dimension(:,:), allocatable :: invl
      BIGREAL, dimension(:), allocatable :: p
      BIGREAL :: sum
      INTEGER :: jprsize,jr1,jr2,jk
      INTEGER :: allocok
C
      jprsize = size(mat_rr,1)
      IF (jprsize.NE.size(mat_rr,2)) GOTO 1000
C
      allocate(invl(jprsize,jprsize),stat=allocok)
      IF (allocok.NE.0) GOTO 1001
      allocate(p(jprsize),stat=allocok)
      IF (allocok.NE.0) GOTO 1001
C
C -1.- Choleski decomposision
C ---------------------------
C Compute L such that: C = L L^T
C
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
C
      DO jr2=1,jprsize
         mat_rr(1:jr2-1,jr2) = FREAL(0.0)
      ENDDO
C
C -2.- Inverse computation
C ------------------------
C Inverse L
C
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
C
C Inverse A
C
      DO jr1=1,jprsize
      DO jr2=1,jprsize
        mat_rr(jr1,jr2) = DOT_PRODUCT(invl(:,jr1),invl(:,jr2))
      ENDDO
      ENDDO
C
      IF (allocated(invl)) deallocate(invl)
      IF (allocated(p)) deallocate(p)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilclc','invmat_rr')
 1001 CALL printerror2(0,1000,3,'utilclc','invmat_rr')
C
 101  WRITE (texterror,*) 'Error in Choleski decomposition'
      CALL printerror2(0,101,3,'utilclc','invmat_rr',comment=texterror)
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE prodmat_wr_rrz_vias(mat_wr,mat_sr,mat_rrz,vectx,
     $     jrdebproj,jrfinproj)
C
C Compute S(w*r)=S(w*r)*C(r*r,z) using TMP(s*r) as working space
C (use C corresponding to local subsystem z in scalar products)
C
C Arguments:  mat_wr : S(w*r) (input and output)
C             mat_sr : TMP(s*r) (output)
C             mat_rr : C(r*r) (input)
C             jddebproj : first index of S to compute
C             jdfinproj : last index of S to compute
C             vectx  : partition vector containing subsystem indices
C
C --- Module declaration
      use mod_cfgxyo
      IMPLICIT NONE
C --- Variable declaration    
      BIGREAL, dimension(:,:), intent(inout) :: mat_wr
      BIGREAL, dimension(:,:), intent(out) :: mat_sr
      BIGREAL, dimension(:,:,0:), intent(in) :: mat_rrz
      BIGREAL, dimension(:), intent(in) :: vectx
      INTEGER, intent(in) :: jrdebproj,jrfinproj
      INTEGER :: jpzsize,jpssize,jprsize,jpwsize
      INTEGER :: jwdeb,jwfin,jsfin,jr1,jr2,jr,js
C
      jpwsize = size(mat_wr,1)
      jpssize = size(mat_sr,1)
      jprsize = size(mat_rrz,1)
      jpzsize = size(mat_rrz,3)
      IF (jprsize.NE.size(mat_wr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_sr,2)) GOTO 1000
      IF (jprsize.NE.size(mat_rrz,2)) GOTO 1000
      IF (jpwsize.NE.size(vectx,1)) GOTO 1000
      IF (jprsize.LT.jrfinproj) GOTO 1000
C
      mat_sr(:,:) = FREAL(0.0)
C
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
                  mat_wr(jwdeb+js-1,jr1)=mat_wr(jwdeb+js-1,jr1) +
     $                mat_sr(js,jr2) *
     $                mat_rrz(jr2,jr1,nint(vectx(jwdeb+js-1)))
               ENDDO
            ENDDO
         ENDDO
 
      ENDDO
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilclc','prodmat_wr_rrz_vias')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE prodmat_wr_rdz_vias(mat_wr,mat_sr,mat_rdz,vectx,
     $     jddebproj,jdfinproj)
C
C Compute S(w*d)=S(w*r)*C(r*d,z) using TMP(s*r) as working space
C (use C corresponding to local subsystem z in scalar products)
C
C Arguments:  mat_wr : S(w*r) (input and output), d <= r
C             mat_sr : TMP(s*r) (output)
C             mat_rr : C(r*d) (input)
C             jddebproj : first index of S to compute
C             jdfinproj : last index of S to compute
C             vectx  : partition vector containing subsystem indices
C
C --- Module declaration
      IMPLICIT NONE
C --- Variable declaration    
      BIGREAL, dimension(:,:), intent(inout) :: mat_wr
      BIGREAL, dimension(:,:), intent(out) :: mat_sr
      BIGREAL, dimension(:,:,0:), intent(in) :: mat_rdz
      BIGREAL, dimension(:), intent(in) :: vectx
      INTEGER, intent(in) :: jddebproj,jdfinproj
      INTEGER :: jpzsize,jpssize,jprsize,jpdsize,jpwsize
      INTEGER :: jwdeb,jwfin,jsfin,jr,jd
C
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
C
      mat_sr(:,:) = FREAL(0.0)
C
      DO jwdeb=1,jpwsize,jpssize
         jwfin=MIN(jpssize+jwdeb-1,jpwsize)
         jsfin=MIN(jpssize,jwfin-jwdeb+1)
         DO jr=1,jprsize
            mat_sr(1:jsfin,jr)=mat_wr(jwdeb:jwfin,jr)
         ENDDO
         mat_wr(jwdeb:jwfin,:) = FREAL(0.0)
         DO jd=jddebproj,jdfinproj
            DO jr=1,jprsize
               mat_wr(jwdeb:jwfin,jd)=mat_wr(jwdeb:jwfin,jd) +
     $                mat_sr(1:jsfin,jr) * 
     $                mat_rdz(jr,jd,nint(vectx(jwdeb:jwfin)))
            ENDDO
         ENDDO
 
      ENDDO
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilclc','prodmat_wr_rdz_vias')
C
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilclc
