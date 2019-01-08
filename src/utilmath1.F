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
C ---                    UTILMATH1.F                             ---
C ---                                                           ---
C ---                                                           ---
C --- original     : 98-04 (C.E. Testut)                        ---
C ---                                                           ---
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h"
C -----------------------------------------------------------------
C --- Subroutines extracted from the eigensystem
C --- subroutine package (eispack)
C --- 
C --- SUBROUTINE rs    : find the eigenvalues and eigenvectors
C ---                    of a real symmetric matrix
C --- FUNCTION pythag  : finds dsqrt(a**2+b**2) without overflow
C ---                    or destructive underflow
C --- SUBROUTINE tql1  : finds the eigenvalues of a symmetric
C ---                    tridiagonal matrix by the ql method
C --- SUBROUTINE tql2  : finds the eigenvalues and eigenvectors of a
C ---                    symmetric tridiagonal matrix by the ql method.
C ---                    the eigenvectors of a full symmetric matrix can also
C ---                    be found if tred2  has been used to reduce this
C ---                    full matrix to tridiagonal form
C --- SUBROUTINE tred1 : reduces a real symmetric matrix to a symmetric
C ---                    tridiagonal matrix using orthogonal similarity
C ---                    transformations
C --- SUBROUTINE tred2 : reduces a real symmetric matrix to a symmetric
C ---                    tridiagonal matrix using and accumulating
C ---                    orthogonal similarity transformations
C --- 
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilmath1
      use mod_main
      IMPLICIT NONE
      PRIVATE

      PUBLIC rs

      CONTAINS
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
      IMPLICIT NONE
C
      INTEGER, intent(in) :: n,nm,matz
      INTEGER, intent(out) :: ierr
      BIGREAL, dimension(:,:), intent(in) :: a
      BIGREAL, dimension(:), intent(out) :: w
      BIGREAL, dimension(:,:), intent(out) :: z
      BIGREAL, dimension(:), intent(inout) :: fv1,fv2
C
C      integer n,nm,ierr,matz
C
C      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
C
C      BIGREAL a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for
c           tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      IF (n .gt. nm) THEN
        ierr = 1 ; RETURN
      ENDIF
c
      SELECT CASE(matz)
      CASE (0)
c       .......... find eigenvalues only ..........
        CALL tred1(nm,n,a,w,fv1,fv2)
        CALL tql1(n,w,fv1,ierr)
      CASE DEFAULT
c       .......... find both eigenvalues and eigenvectors ..........
        CALL tred2(nm,n,a,w,fv1,z)
        CALL tql2(nm,n,w,fv1,z,ierr)
      END SELECT
C
      RETURN
      END SUBROUTINE
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function pythag(a,b)
      IMPLICIT NONE
C
      BIGREAL pythag
      BIGREAL a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      BIGREAL p,r,s,t,u
      p = FREAL(MAX(ABS(a),ABS(b)))
      if (p .eq. FREAL(0.0)) go to 20
      r = (FREAL(MIN(ABS(a),ABS(b)))/p)**2
   10 continue
         t = FREAL(4.0) + r
         if (t .eq. FREAL(4.0)) go to 20
         s = r/t
         u = FREAL(1.0) + FREAL(2.0)*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end function
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine tql1(n,d,e,ierr)
C
      IMPLICIT NONE
C
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      BIGREAL d(n),e(n)
      BIGREAL c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
c
c     this subroutine is a translation of the algol procedure tql1,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = FREAL(0.0)
      tst1 = FREAL(0.0)
      e(n) = FREAL(0.0)
c
      do 290 l = 1, n
         j = 0
         h = FREAL(ABS(d(l)) + ABS(e(l)))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + FREAL(ABS(e(m)))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (FREAL(2.0) * e(l))
         r = pythag(p,FREAL(1.0))
         d(l) = e(l) / (p + FREAL(SIGN(r,p)))
         d(l1) = e(l) * (p + FREAL(SIGN(r,p)))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = FREAL(1.0)
         c2 = c
         el1 = e(l1)
         s = FREAL(0.0)
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + FREAL(ABS(e(l)))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end subroutine
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine tql2(nm,n,d,e,z,ierr)
      IMPLICIT NONE
C
      INTEGER, intent(in) :: n,nm
      INTEGER, intent(out) :: ierr
      BIGREAL, dimension(:,:), intent(inout) :: z
      BIGREAL, dimension(:), intent(out) :: d,e

      integer i,j,k,l,m,ii,l1,l2,mml
      BIGREAL c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = FREAL(0.0)
      tst1 = FREAL(0.0)
      e(n) = FREAL(0.0)
c
      do 240 l = 1, n
         j = 0
         h = FREAL(ABS(d(l)) + ABS(e(l)))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + FREAL(ABS(e(m)))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (FREAL(2.0) * e(l))
         r = pythag(p,FREAL(1.0))
         d(l) = e(l) / (p + FREAL(SIGN(r,p)))
         d(l1) = e(l) * (p + FREAL(SIGN(r,p)))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = FREAL(1.0)
         c2 = c
         el1 = e(l1)
         s = FREAL(0.0)
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
C         tst2 = tst1 + dabs(e(l))
         tst2 = tst1 + FREAL(ABS(e(l)))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end subroutine
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine tred1(nm,n,a,d,e,e2)
      IMPLICIT NONE
C
      integer i,j,k,l,n,ii,nm,jp1
      BIGREAL a(nm,n),d(n),e(n),e2(n)
      BIGREAL f,g,h,scale
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = FREAL(0.0)
         scale = FREAL(0.0)
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + FREAL(ABS(d(k)))
c
         if (scale .ne. FREAL(0.0)) go to 140
c
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = FREAL(0.0)
  125    continue
c
  130    e(i) = FREAL(0.0)
         e2(i) = FREAL(0.0)
         go to 300
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -FREAL(SIGN(SQRT(h),f))
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = FREAL(0.0)
c
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = FREAL(0.0)
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         h = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280    continue
c
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
c
  300 continue
c
      return
      end subroutine
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -----------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine tred2(nm,n,a,d,e,z)
      IMPLICIT NONE
C
      INTEGER, intent(in) :: n,nm
      BIGREAL, dimension(:,:), intent(in) :: a
      BIGREAL, dimension(:), intent(out) :: d,e
      BIGREAL, dimension(:,:), intent(out) :: z

      integer i,j,k,l,ii,jp1
      BIGREAL f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
C         h = 0.0d0
         h = FREAL(0.0)
         scale = FREAL(0.0)
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + FREAL(ABS(d(k)))
c
         if (scale .ne. FREAL(0.0)) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = FREAL(0.0)
            z(j,i) = FREAL(0.0)
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -FREAL(SIGN(SQRT(h),f))
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = FREAL(0.0)
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = FREAL(0.0)
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = FREAL(0.0)
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = FREAL(1.0)
         h = d(i)
         if (h .eq. FREAL(0.0)) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = FREAL(0.0)
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = FREAL(0.0)
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = FREAL(0.0)
  520 continue
c
      z(n,n) = FREAL(1.0)
      e(1) = FREAL(0.0)
      return
      end subroutine
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilmath1
