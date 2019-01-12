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
! ---                   UTILZONE.F90                              ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2000-02 (J.M. Brankart)                    ---
! --- revised      : 2003-04 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! ---
! --- FUNCTION inarea    : Check if input point is inside
! ---                      a polygonal contour
! --- FUNCTION fxyzcnt   : List of function for piecewise
! ---                      field definition
! --- FUNCTION fctcorrel : Evaluate correlation function
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilzone
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC inarea, fxyzcnt, fctcorrel

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION inarea (kgrdpt,kjc)
!---------------------------------------------------------------------
!
!  Purpose : Check if input point is inside polygonal contour number kjc
!  -------
!  Method : Use a meridional research line including input point
!  ------   Count the number of crossing between contour and research line
!           on both side of input point
!           If this number is odd, the point is inside the contour
!           If this number is even, the point is outside the contour
!
!  Input :  kgrdpt : input point coordinates
!  -----    kjc    : contour index
!
!  Output : inarea : TRUE=inside (or on the border line)
!  ------            FALSE=outside
!
!---------------------------------------------------------------------
      use mod_main
      use mod_cfgxyo
      use mod_cont
      IMPLICIT NONE
      INTEGER, intent(in) :: kjc
      TYPE (type_gridij), intent(in) :: kgrdpt
      LOGICAL :: inarea
!---------------------------------------------------------------------
      LOGICAL :: nonborder
      INTEGER :: cinf, csup, allocok, jp, cshift, jk
      BIGREAL :: scont, delta
      INTEGER, dimension(:), allocatable :: cidx
!---------------------------------------------------------------------
!
! Check if polygonal contour is closed
! (last point is identical to first point)
      IF ( (contij(1,kjc)%longi.NE.contij(jpp(kjc),kjc)%longi) &
     &  .OR. (contij(1,kjc)%latj.NE.contij(jpp(kjc),kjc)%latj) ) GOTO 102
!
! -1.- Build cidx array: ordered list of polygon vertex indices
! -------------------------------------------------------------
! The polygonal segment (between 2 succesive vertices) 
! will be checked successively following this order
!
      allocate ( cidx(1:jpp(kjc)), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      cidx(:) = 0
!
! Check if input point is not identical to first contour vertex
! (if yes, input point point is on the border line and we
! go directly to the end of the program [nonborder=FALSE])
      nonborder = .FALSE.
      nonborder = nonborder .OR. (kgrdpt%longi.NE.contij(1,kjc)%longi)
      nonborder = nonborder .OR. (kgrdpt%latj.NE.contij(1,kjc)%latj)
!
      IF (nonborder) THEN
!
! Find first vertex point with a different longitude than input point
! (i.e. which does not belong to the research line)
         jp = 1
         DO WHILE ( (kgrdpt%longi.EQ.contij(jp,kjc)%longi)  &
     &                        .AND. (jp.LT.jpp(kjc)) )
            jp = jp + 1
         ENDDO
         IF (jp.EQ.jpp(kjc)) GOTO 103
!
! This is the first vertex of the list,
! the other vertices are just shifted with respect to
! the input list (so that the order remains 'area on the left')
         cshift = jp - 1
         DO jp = 1,jpp(kjc)
            cidx(jp) = jp + cshift
            IF (cidx(jp).GT.jpp(kjc)) cidx(jp) = cidx(jp) - jpp(kjc) + 1
         ENDDO
!
      ENDIF
!
! -2.- Loop over polygonal segments following
! the order specified in cidx array
! -------------------------------------------
!
      cinf = 0     
      csup = 0     
!
      jp = 1
      DO WHILE ( nonborder .AND. (jp.LT.jpp(kjc)) )
!
! Check position of the meridional research line
! with respect to current polygonal segment
! scont = 0  : the second vertex belongs to the research line
!              (the first vertex never does...)
! scont < 0  : the research line is crossing the segment
! scont > 0  : the research line is not crossing the segment
!              (-> go directly to next polygonal segment)
         scont = (contij(cidx(jp),kjc)%longi-kgrdpt%longi) * &
     &              (contij(cidx(jp+1),kjc)%longi-kgrdpt%longi)
!
         IF (scont.EQ.FREAL(0.0)) THEN
! The second vertex belongs to the research line
! This is the problematic case because several vertices can
! be aligned on the research line
!
! Find next vertex not belonging to the research line
! scont < 0  : the research line is crossing the contour
!              even if it contains one or several successive segment
! scont > 0  : the research line is tangent to the contour
!              and contains one or several successive segment
!              (thos does not modified parity...)
            jk = 1
            DO WHILE (scont.EQ.FREAL(0.0)) 
               jk = jk + 1
               scont = (contij(cidx(jp),kjc)%longi-kgrdpt%longi) * &
     &                 (contij(cidx(jp+jk),kjc)%longi-kgrdpt%longi)
            ENDDO
!
            IF (scont.LT.FREAL(0.0)) THEN
!
! Check if crossing segment contains the input point
               scont = (contij(cidx(jp+1),kjc)%latj-kgrdpt%latj) * &
     &                   (contij(cidx(jp+jk-1),kjc)%latj-kgrdpt%latj)
               IF (scont.LE.FREAL(0.0)) THEN
! Input point is on the border line
                  nonborder = .FALSE.
               ELSE
!
! Check if tangential segment is North or South of the input point
                  IF (contij(cidx(jp+1),kjc)%latj .GT. kgrdpt%latj) THEN
! The contour is crossing the research line North of the input point
                     csup = csup + 1
                  ELSE
! The contour is crossing the research line South of the input point
                     cinf = cinf + 1
                  ENDIF
               ENDIF
            ELSE
!
! Check if tangential segment contains the input point
               scont = (contij(cidx(jp+1),kjc)%latj-kgrdpt%latj) * &
     &                   (contij(cidx(jp+jk-1),kjc)%latj-kgrdpt%latj)
               IF (scont.LE.FREAL(0.0)) THEN
! Input point is on the border line
                  nonborder = .FALSE.
               ENDIF
            ENDIF
!
! Go to next vertex not belonging to the research line
            jp = jp + jk
!
         ELSE
!
            IF (scont.LT.FREAL(0.0)) THEN
! The research line is crossing the segment
!
! Check if input point lays between minimum and maximum
! latitudes of the segment
! If not (the usual case), the positionning requires
! less computational time
               scont = (contij(cidx(jp),kjc)%latj-kgrdpt%latj) * &
     &                   (contij(cidx(jp+1),kjc)%latj-kgrdpt%latj)
               IF (scont.GT.FREAL(0.0)) THEN
! Input point latitude is outside segment latitude range
                  IF (contij(cidx(jp),kjc)%latj .GT. kgrdpt%latj) THEN
! The segment is crossing the research line North of the input point
                     csup = csup + 1
                  ELSE
! The segment is crossing the research line South of the input point
                     cinf = cinf + 1
                  ENDIF
               ELSE
! Input point latitude is inside segment latitude range
                  delta = (contij(cidx(jp+1),kjc)%latj - &
     &                             contij(cidx(jp),kjc)%latj) &
     &                  / (contij(cidx(jp+1),kjc)%longi -  &
     &                             contij(cidx(jp),kjc)%longi)
                  scont = contij(cidx(jp),kjc)%latj - kgrdpt%latj  &
     &                       - delta * &
     &                   ( contij(cidx(jp),kjc)%longi - kgrdpt%longi )
                  IF (scont.EQ.FREAL(0.0)) THEN
! Input point is on the border line
                     nonborder = .FALSE.
                  ELSE
                     IF (scont.GT.FREAL(0.0)) THEN
! The segment is crossing the research line North of the input point
                        csup = csup + 1
                     ELSE
! The segment is crossing the research line South of the input point
                        cinf = cinf + 1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
!
! Go to next vertex
            jp = jp + 1
!
         ENDIF
!
      ENDDO
!
! -3.- Take final decision
! ------------------------
!
      IF (nonborder) THEN
! Input point is not on the border line
         IF (MOD(csup,2).EQ.0) THEN
            IF (MOD(cinf,2).EQ.0) THEN
! The number of contour segments crossing research line
! is even on both side of input point
! ==> Input point is outside the contour
               inarea = .FALSE.
            ELSE
               GOTO 1000
            ENDIF
         ELSE
            IF (MOD(cinf,2).EQ.0) THEN
               GOTO 1000
            ELSE
! The number of contour segments crossing research line
! is odd on both side of input point
! ==> Input point is inside the contour
               inarea = .TRUE.
            ENDIF
         ENDIF
      ELSE
! Input point is on the border line
         inarea = .TRUE.
      ENDIF
!
      IF (allocated(cidx)) deallocate(cidx)
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilzone','inarea')
 1001 CALL printerror2(0,1001,3,'utilzone','inarea')
!
 102  WRITE (texterror,*) 'Polygonal contour ',kjc,' is not closed!'
      CALL printerror2(0,102,3,'utilzone','inarea',comment=texterror)
 103  WRITE (texterror,*) 'All points of contour ',kjc,' are colinear!'
      CALL printerror2(0,103,3,'utilzone','inarea',comment=texterror)
!
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION fxyzcnt (kx,kdimcnt,kjlay,kjc)
!---------------------------------------------------------------------
!  
!  Purpose : List of function for piecewise field definition
!  ------- 
!  Method : Set function type as a function of coordinate,
!  ------   contour index and slice index
!           Compute function value
!
!  Input :  kx      : function argument (in [-1,1])
!  -----    kdimcnt : coordinate along which the function is defined
!           kjc     : contour index
!           kjlay   : slice index
!
!  Output : fxyzcnt : function value (=f(kx))
!  ------
!---------------------------------------------------------------------
      use mod_main
      use mod_cfgxyo
      use mod_cont
      IMPLICIT NONE
      BIGREAL, intent(in) :: kx
      INTEGER, intent(in) :: kdimcnt
      INTEGER, intent(in) :: kjlay
      INTEGER, intent(in) :: kjc
      BIGREAL :: fxyzcnt
!---------------------------------------------------------------------
      BIGREAL :: pival, fx
      INTEGER :: typefct
!---------------------------------------------------------------------
!
! Set function type as a function of coordinate,
! contour index and slice index
! (contfx, contfy, contfz are defined in "mod_cont.F90")
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
!
! Transform argument from [-1,1] to [0,1] interval
      fx = FREAL(0.5) * ( FREAL(1.0) + kx )
!
! Compute function value
      SELECT CASE(typefct)
      CASE(0)
! --- f(x)=1
         fxyzcnt = FREAL(1.0)
      CASE(1)
! --- f(x)=x
         fxyzcnt = fx
      CASE(2)
! --- f(x)=[1-cos(x*pi)]/2
         pival = FREAL(2.0) * ASIN(FREAL(1.0))
         fxyzcnt = FREAL(0.5) * ( FREAL(1.0) - COS(fx*pival) )
      CASE DEFAULT
         GOTO 102
      END SELECT
!---------------------------------------------------------------------
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilzone','fxyzcnt')
 1001 CALL printerror2(0,1001,3,'utilzone','fxyzcnt')
!
 102  WRITE (texterror,*) 'Bad function type'
      CALL printerror2(0,102,3,'utilzone','fxyzcnt',comment=texterror)
!
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION fctcorrel (jil,jjl,jkl,jtl,jic,jjc,jkc,jtc, &
     &                    clgt_i,clgt_j,clgt_k,clgt_t, &
     &                    ctyp_i,ctyp_j,ctyp_k,ctyp_t)
!---------------------------------------------------------------------
!
!  Purpose : Evaluate correlation function between 2 positions
!  ------- 
!  Method :  Use requested function type to evaluate correlation
!  ------    as a function of the distance between the 2 positions
!
!  Input :   jic,jjc,jkc,jtc : indices of first position in model grid
!  -----     jil,jjl,jkl,jtl : indices of second position in model grid
!            clgt_i,clgt_j,clgt_k,clgt_t : correlation scale
!                                          along each coordinate
!            ctyp_i,ctyp_j,ctyp_k,ctyp_t : correlation function type
!                                          along each coordinate
!
!  Output :  fctcorrel : correlation value
!  ------
!---------------------------------------------------------------------
      use mod_main
      use mod_cfgxyo
      IMPLICIT NONE
!---------------------------------------------------------------------
      INTEGER, intent(in) :: jil,jjl,jkl,jtl,jic,jjc,jkc,jtc
      BIGREAL, intent(in) :: clgt_i,clgt_j,clgt_k,clgt_t
      INTEGER, intent(in) :: ctyp_i,ctyp_j,ctyp_k,ctyp_t
      BIGREAL :: fctcorrel
!---------------------------------------------------------------------
      BIGREAL :: li, lj, lk, lt, lr, alr, blr
      BIGREAL :: fij, fk, ft
      BIGREAL, parameter :: a = 2.103803_kr
      BIGREAL, parameter :: b = 3.336912_kr
!---------------------------------------------------------------------
! Compute relative distance between 2 positions
      li = ABS( FREAL(jic -jil) ) / clgt_i
      lj = ABS( FREAL(jjc -jjl) ) / clgt_j
      lk = ABS( FREAL(jkc -jkl) ) / clgt_k
      lt = ABS( FREAL(jtc -jtl) ) / clgt_t
      lr = SQRT( li * li + lj * lj ) 
!
! Correlation on the horizontal
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
         fij = ( 1.0_kr + blr + blr * blr / 6.0_kr  &
     &                  - blr * blr * blr / 6.0_kr ) * exp(-blr)
      CASE(4)
         alr = a * lr
         fij = ( 1.0_kr + alr - alr * alr * alr / 3.0_kr ) * exp(-alr)
      CASE DEFAULT
         GOTO 101
      END SELECT
!
! Correlation on the vertical
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
!
! Correlation in time
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
!
      fctcorrel = fij * fk * ft
!
      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilzone','fctcorrel')
!
 101  WRITE (texterror,*) 'Bad correlation type'
      CALL printerror2(0,101,3,'utilzone','fctcorrel',comment=texterror)
!
      END FUNCTION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilzone
