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
C ---                   UTILZONE.F                              ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 2000-02 (J.M. Brankart)                    ---
C --- revised      : 2003-04 (J.M. Brankart)                    ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C ---
C --- FUNCTION inarea    : Check if input point is inside
C ---                      a polygonal contour
C --- FUNCTION fxyzcnt   : List of function for piecewise
C ---                      field definition
C --- FUNCTION fctcorrel : Evaluate correlation function
C ---
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilzone
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC inarea, fxyzcnt, fctcorrel

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION inarea (kgrdpt,kjc)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Check if input point is inside polygonal contour number kjc
CCC  -------
CCC  Method : Use a meridional research line including input point
CCC  ------   Count the number of crossing between contour and research line
CCC           on both side of input point
CCC           If this number is odd, the point is inside the contour
CCC           If this number is even, the point is outside the contour
CCC
CCC  Input :  kgrdpt : input point coordinates
CCC  -----    kjc    : contour index
CCC
CCC  Output : inarea : TRUE=inside (or on the border line)
CCC  ------            FALSE=outside
CCC
CCC---------------------------------------------------------------------
      use mod_main
      use mod_cfgxyo
      use mod_cont
      IMPLICIT NONE
      INTEGER, intent(in) :: kjc
      TYPE (type_gridij), intent(in) :: kgrdpt
      LOGICAL :: inarea
CCC---------------------------------------------------------------------
      LOGICAL :: nonborder
      INTEGER :: cinf, csup, allocok, jp, cshift, jk
      BIGREAL :: scont, delta
      INTEGER, dimension(:), allocatable :: cidx
CCC---------------------------------------------------------------------
C
C Check if polygonal contour is closed
C (last point is identical to first point)
      IF ( (contij(1,kjc)%longi.NE.contij(jpp(kjc),kjc)%longi)
     $  .OR. (contij(1,kjc)%latj.NE.contij(jpp(kjc),kjc)%latj) ) GOTO 102
C
C -1.- Build cidx array: ordered list of polygon vertex indices
C -------------------------------------------------------------
C The polygonal segment (between 2 succesive vertices) 
C will be checked successively following this order
C
      allocate ( cidx(1:jpp(kjc)), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      cidx(:) = 0
C
C Check if input point is not identical to first contour vertex
C (if yes, input point point is on the border line and we
C go directly to the end of the program [nonborder=FALSE])
      nonborder = .FALSE.
      nonborder = nonborder .OR. (kgrdpt%longi.NE.contij(1,kjc)%longi)
      nonborder = nonborder .OR. (kgrdpt%latj.NE.contij(1,kjc)%latj)
C
      IF (nonborder) THEN
C
C Find first vertex point with a different longitude than input point
C (i.e. which does not belong to the research line)
         jp = 1
         DO WHILE ( (kgrdpt%longi.EQ.contij(jp,kjc)%longi) 
     $                        .AND. (jp.LT.jpp(kjc)) )
            jp = jp + 1
         ENDDO
         IF (jp.EQ.jpp(kjc)) GOTO 103
C
C This is the first vertex of the list,
C the other vertices are just shifted with respect to
C the input list (so that the order remains 'area on the left')
         cshift = jp - 1
         DO jp = 1,jpp(kjc)
            cidx(jp) = jp + cshift
            IF (cidx(jp).GT.jpp(kjc)) cidx(jp) = cidx(jp) - jpp(kjc) + 1
         ENDDO
C
      ENDIF
C
C -2.- Loop over polygonal segments following
C the order specified in cidx array
C -------------------------------------------
C
      cinf = 0     
      csup = 0     
C
      jp = 1
      DO WHILE ( nonborder .AND. (jp.LT.jpp(kjc)) )
C
C Check position of the meridional research line
C with respect to current polygonal segment
C scont = 0  : the second vertex belongs to the research line
C              (the first vertex never does...)
C scont < 0  : the research line is crossing the segment
C scont > 0  : the research line is not crossing the segment
C              (-> go directly to next polygonal segment)
         scont = (contij(cidx(jp),kjc)%longi-kgrdpt%longi) *
     $              (contij(cidx(jp+1),kjc)%longi-kgrdpt%longi)
C
         IF (scont.EQ.FREAL(0.0)) THEN
C The second vertex belongs to the research line
C This is the problematic case because several vertices can
C be aligned on the research line
C
C Find next vertex not belonging to the research line
C scont < 0  : the research line is crossing the contour
C              even if it contains one or several successive segment
C scont > 0  : the research line is tangent to the contour
C              and contains one or several successive segment
C              (thos does not modified parity...)
            jk = 1
            DO WHILE (scont.EQ.FREAL(0.0)) 
               jk = jk + 1
               scont = (contij(cidx(jp),kjc)%longi-kgrdpt%longi) *
     $                 (contij(cidx(jp+jk),kjc)%longi-kgrdpt%longi)
            ENDDO
C
            IF (scont.LT.FREAL(0.0)) THEN
C
C Check if crossing segment contains the input point
               scont = (contij(cidx(jp+1),kjc)%latj-kgrdpt%latj) *
     $                   (contij(cidx(jp+jk-1),kjc)%latj-kgrdpt%latj)
               IF (scont.LE.FREAL(0.0)) THEN
C Input point is on the border line
                  nonborder = .FALSE.
               ELSE
C
C Check if tangential segment is North or South of the input point
                  IF (contij(cidx(jp+1),kjc)%latj .GT. kgrdpt%latj) THEN
C The contour is crossing the research line North of the input point
                     csup = csup + 1
                  ELSE
C The contour is crossing the research line South of the input point
                     cinf = cinf + 1
                  ENDIF
               ENDIF
            ELSE
C
C Check if tangential segment contains the input point
               scont = (contij(cidx(jp+1),kjc)%latj-kgrdpt%latj) *
     $                   (contij(cidx(jp+jk-1),kjc)%latj-kgrdpt%latj)
               IF (scont.LE.FREAL(0.0)) THEN
C Input point is on the border line
                  nonborder = .FALSE.
               ENDIF
            ENDIF
C
C Go to next vertex not belonging to the research line
            jp = jp + jk
C
         ELSE
C
            IF (scont.LT.FREAL(0.0)) THEN
C The research line is crossing the segment
C
C Check if input point lays between minimum and maximum
C latitudes of the segment
C If not (the usual case), the positionning requires
C less computational time
               scont = (contij(cidx(jp),kjc)%latj-kgrdpt%latj) *
     $                   (contij(cidx(jp+1),kjc)%latj-kgrdpt%latj)
               IF (scont.GT.FREAL(0.0)) THEN
C Input point latitude is outside segment latitude range
                  IF (contij(cidx(jp),kjc)%latj .GT. kgrdpt%latj) THEN
C The segment is crossing the research line North of the input point
                     csup = csup + 1
                  ELSE
C The segment is crossing the research line South of the input point
                     cinf = cinf + 1
                  ENDIF
               ELSE
C Input point latitude is inside segment latitude range
                  delta = (contij(cidx(jp+1),kjc)%latj -
     $                             contij(cidx(jp),kjc)%latj)
     $                  / (contij(cidx(jp+1),kjc)%longi - 
     $                             contij(cidx(jp),kjc)%longi)
                  scont = contij(cidx(jp),kjc)%latj - kgrdpt%latj 
     $                       - delta *
     $                   ( contij(cidx(jp),kjc)%longi - kgrdpt%longi )
                  IF (scont.EQ.FREAL(0.0)) THEN
C Input point is on the border line
                     nonborder = .FALSE.
                  ELSE
                     IF (scont.GT.FREAL(0.0)) THEN
C The segment is crossing the research line North of the input point
                        csup = csup + 1
                     ELSE
C The segment is crossing the research line South of the input point
                        cinf = cinf + 1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
C
C Go to next vertex
            jp = jp + 1
C
         ENDIF
C
      ENDDO
C
C -3.- Take final decision
C ------------------------
C
      IF (nonborder) THEN
C Input point is not on the border line
         IF (MOD(csup,2).EQ.0) THEN
            IF (MOD(cinf,2).EQ.0) THEN
C The number of contour segments crossing research line
C is even on both side of input point
C ==> Input point is outside the contour
               inarea = .FALSE.
            ELSE
               GOTO 1000
            ENDIF
         ELSE
            IF (MOD(cinf,2).EQ.0) THEN
               GOTO 1000
            ELSE
C The number of contour segments crossing research line
C is odd on both side of input point
C ==> Input point is inside the contour
               inarea = .TRUE.
            ENDIF
         ENDIF
      ELSE
C Input point is on the border line
         inarea = .TRUE.
      ENDIF
C
      IF (allocated(cidx)) deallocate(cidx)
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilzone','inarea')
 1001 CALL printerror2(0,1001,3,'utilzone','inarea')
C
 102  WRITE (texterror,*) 'Polygonal contour ',kjc,' is not closed!'
      CALL printerror2(0,102,3,'utilzone','inarea',comment=texterror)
 103  WRITE (texterror,*) 'All points of contour ',kjc,' are colinear!'
      CALL printerror2(0,103,3,'utilzone','inarea',comment=texterror)
C
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION fxyzcnt (kx,kdimcnt,kjlay,kjc)
CCC---------------------------------------------------------------------
CCC  
CCC  Purpose : List of function for piecewise field definition
CCC  ------- 
CCC  Method : Set function type as a function of coordinate,
CCC  ------   contour index and slice index
CCC           Compute function value
CCC
CCC  Input :  kx      : function argument (in [-1,1])
CCC  -----    kdimcnt : coordinate along which the function is defined
CCC           kjc     : contour index
CCC           kjlay   : slice index
CCC
CCC  Output : fxyzcnt : function value (=f(kx))
CCC  ------
CCC---------------------------------------------------------------------
      use mod_main
      use mod_cfgxyo
      use mod_cont
      IMPLICIT NONE
      BIGREAL, intent(in) :: kx
      INTEGER, intent(in) :: kdimcnt
      INTEGER, intent(in) :: kjlay
      INTEGER, intent(in) :: kjc
      BIGREAL :: fxyzcnt
CCC---------------------------------------------------------------------
      BIGREAL :: pival, fx
      INTEGER :: typefct
CCC---------------------------------------------------------------------
C
C Set function type as a function of coordinate,
C contour index and slice index
C (contfx, contfy, contfz are defined in "mod_cont.F")
      SELECT CASE(kdimcnt)
      CASE(1)
         typefct = contfx(kjlay,kjc)
      CASE(2)
         typefct = contfy(kjlay,kjc)
      CASE(3)
         typefct = contfz(kjlay,kjc)
      CASE DEFAULT
         GOTO 1000
      END SELECT
C
C Transform argument from [-1,1] to [0,1] interval
      fx = FREAL(0.5) * ( FREAL(1.0) + kx )
C
C Compute function value
      SELECT CASE(typefct)
      CASE(0)
C --- f(x)=1
         fxyzcnt = FREAL(1.0)
      CASE(1)
C --- f(x)=x
         fxyzcnt = fx
      CASE(2)
C --- f(x)=[1-cos(x*pi)]/2
         pival = FREAL(2.0) * ASIN(FREAL(1.0))
         fxyzcnt = FREAL(0.5) * ( FREAL(1.0) - COS(fx*pival) )
      CASE DEFAULT
         GOTO 102
      END SELECT
CCC---------------------------------------------------------------------
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilzone','fxyzcnt')
 1001 CALL printerror2(0,1001,3,'utilzone','fxyzcnt')
C
 102  WRITE (texterror,*) 'Bad function type'
      CALL printerror2(0,102,3,'utilzone','fxyzcnt',comment=texterror)
C
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION fctcorrel (jil,jjl,jkl,jtl,jic,jjc,jkc,jtc,
     $                    clgt_i,clgt_j,clgt_k,clgt_t,
     $                    ctyp_i,ctyp_j,ctyp_k,ctyp_t)
CCC---------------------------------------------------------------------
CCC
CCC  Purpose : Evaluate correlation function between 2 positions
CCC  ------- 
CCC  Method :  Use requested function type to evaluate correlation
CCC  ------    as a function of the distance between the 2 positions
CCC
CCC  Input :   jic,jjc,jkc,jtc : indices of first position in model grid
CCC  -----     jil,jjl,jkl,jtl : indices of second position in model grid
CCC            clgt_i,clgt_j,clgt_k,clgt_t : correlation scale
CCC                                          along each coordinate
CCC            ctyp_i,ctyp_j,ctyp_k,ctyp_t : correlation function type
CCC                                          along each coordinate
CCC
CCC  Output :  fctcorrel : correlation value
CCC  ------
CCC---------------------------------------------------------------------
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
CCC---------------------------------------------------------------------
      INTEGER, intent(in) :: jil,jjl,jkl,jtl,jic,jjc,jkc,jtc
      BIGREAL, intent(in) :: clgt_i,clgt_j,clgt_k,clgt_t
      INTEGER, intent(in) :: ctyp_i,ctyp_j,ctyp_k,ctyp_t
      BIGREAL :: fctcorrel
CCC---------------------------------------------------------------------
      BIGREAL :: li, lj, lk, lt, lr, alr, blr
      BIGREAL :: fij, fk, ft
      BIGREAL, parameter :: a = 2.103803_kr
      BIGREAL, parameter :: b = 3.336912_kr
CCC---------------------------------------------------------------------
C Compute relative distance between 2 positions
      li = ABS( FREAL(jic -jil) ) / clgt_i
      lj = ABS( FREAL(jjc -jjl) ) / clgt_j
      lk = ABS( FREAL(jkc -jkl) ) / clgt_k
      lt = ABS( FREAL(jtc -jtl) ) / clgt_t
      lr = SQRT( li * li + lj * lj ) 
C
C Correlation on the horizontal
      SELECT CASE(ctyp_i)
      CASE(0)
         IF (lr.EQ.FREAL(0.0)) THEN
            fij = FREAL(1.0)
         ELSE
            fij = FREAL(0.0)
         ENDIF
      CASE(1)
         fij = EXP( - lr * lr )
      CASE(2)
         alr = a * lr
         fij = ( 1.0_kr + alr + alr * alr / 3.0_kr ) * exp(-alr)
      CASE(3)
         blr = b * lr
         fij = ( 1.0_kr + blr + blr * blr / 6.0_kr 
     $                  - blr * blr * blr / 6.0_kr ) * exp(-blr)
      CASE(4)
         alr = a * lr
         fij = ( 1.0_kr + alr - alr * alr * alr / 3.0_kr ) * exp(-alr)
      CASE DEFAULT
         GOTO 101
      END SELECT
C
C Correlation on the vertical
      SELECT CASE(ctyp_k)
      CASE(0)
         IF (lr.EQ.FREAL(0.0)) THEN
            fk = FREAL(1.0)
         ELSE
            fk = FREAL(0.0)
         ENDIF
      CASE(1)
         fk = EXP( - lk * lk )
      CASE DEFAULT
         GOTO 101
      END SELECT
C
C Correlation in time
      SELECT CASE(ctyp_t)
      CASE(0)
         IF (lr.EQ.FREAL(0.0)) THEN
            ft = FREAL(1.0)
         ELSE
            ft = FREAL(0.0)
         ENDIF
      CASE(1)
         ft = EXP( - lt * lt )
      CASE DEFAULT
         GOTO 101
      END SELECT
C
      fctcorrel = fij * fk * ft
C
      RETURN
C
C --- error management
C
 1000 CALL printerror2(0,1000,1,'utilzone','fctcorrel')
C
 101  WRITE (texterror,*) 'Bad correlation type'
      CALL printerror2(0,101,3,'utilzone','fctcorrel',comment=texterror)
C
      END FUNCTION
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilzone
