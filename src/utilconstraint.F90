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
! ---                   UTILCONSTRAINT.F90                      ---
! ---                                                           ---
! ---                                                           ---
! --- original     : 2022-04 (J.M. Brankart)                    ---
! ---                                                           ---
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include "config.main.h90"
! -----------------------------------------------------------------
! --- 
! --- SUBROUTINE apply_constraint : apply dynamical constraint
! ---
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE utilconstraint
      use mod_main
      use flowsampler_grid
      use flowsampler_units
      use flowsampler_kinematics
      use flowsampler_dynamics

      IMPLICIT NONE

      PUBLIC apply_constraint

      ! Public variables
      LOGICAL, PUBLIC, save :: dyn_constraint=.FALSE. ! Apply dynamical constraint

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE apply_constraint(vectx)
!---------------------------------------------------------------------
!
!  Purpose : Apply dynamical constraint on input vector
!
!  Method : extract the required 2D slice of variables from 1D vector
!  ------   compute dependent variables from free variables
!
!  Input : vectx      : Vx vector 
!  -----
!---------------------------------------------------------------------
! modules
! =======
      use mod_main
      use mod_cfgxyo
      use mod_coord
      use mod_spacexyo , only : jpxend
      use utilvct
      IMPLICIT NONE
!----------------------------------------------------------------------
! header declarations
! ===================
      BIGREAL, dimension(:), intent(inout) :: vectx
!----------------------------------------------------------------------
! local declarations
! ==================
      BIGREAL, dimension(:,:), allocatable :: adt           ! absolute dynamic topography
      BIGREAL, dimension(:,:), allocatable :: u, u0         ! zonal velocity
      BIGREAL, dimension(:,:), allocatable :: v, v0         ! meridional velocity
      BIGREAL, dimension(:,:), allocatable :: omega, omega0 ! relative vorticity
      BIGREAL, dimension(:,:), allocatable :: xi            ! advection misfit
      INTEGER :: flagxyo,jk,jt,jvar,indvar,jx
      INTEGER :: jpisize,jpjsize,jpksize,jptsize
      INTEGER :: nbr,allocok
      BIGREAL :: dummy, dt
!----------------------------------------------------------------------
!
      IF (nprint.GE.1) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/.../aply_constraint :'
         WRITE(numout,*) '         apply dynamical constraint'
      ENDIF

      flagxyo=1
! Check dimensions : the full state vector must be recomposed
! before calling this routine
      IF (SIZE(vectx,1).NE.jpxend) GOTO 1000

! Check that grid is regular
      IF (MAXVAL(varngrd(1:varend)).GT.1) GOTO 1000

! Allocate 2D physical arrays
      jpisize=MAXVAL(var_jpi(1:varend))
      jpjsize=MAXVAL(var_jpj(1:varend))
      allocate (adt(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      allocate (u(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      allocate (v(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      allocate (omega(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      allocate (xi(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      allocate (u0(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      allocate (v0(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001
      allocate (omega0(1:jpisize,1:jpjsize), stat=allocok )
      IF (allocok.NE.0) GOTO 1001

! Initialize grid in flowsampler
      nlon = jpisize
      nlat = jpjsize
      lonmin = longi(1)
      lonmax = longi(jpisize)
      latmin = latj(1)
      latmax = latj(jpjsize)
      CALL defgrid()
      physical_units=.TRUE.
      ellipsoid_correction=.FALSE.
      normalize_residual=.TRUE.

! Loop on 2D slice in time and along the vertical
      jpksize=MAXVAL(var_jpk(1:varend))
      jptsize=MAXVAL(var_jpt(1:varend))
      DO jt=1,jptsize
      DO jk=1,jpksize

! Get free variable (adt)
        DO jvar = 1,varend
          indvar=var_ord(jvar)
          IF (var_jpk(indvar).NE.jpksize) GOTO 1000
          IF (var_jpt(indvar).NE.jptsize) GOTO 1000

          IF (var_nam(indvar).EQ.'ADT  ') THEN
            adt(:,:) = FREAL(0.0)

            jx=var_sli_idx(indvar,jk,jt)
            CALL unmk8vct(vectx(jx:),adt(:,:),jk,jt, &
     &                    jvar,nbr,spval,flagxyo)

          ENDIF

        ENDDO

! Apply constraint
        u0(:,:) = u(:,:) ; v0(:,:) = v(:,:) ; omega0(:,:) = omega(:,:)
        CALL velocity(adt,u,v)
        CALL vorticity(u,v,omega)

        IF (jt.GT.1) THEN
          dt= 86400. ; dummy=0.
          print *, 'timestep',jt
          !print *, 'u0',minval(u0),maxval(u0)
          !print *, 'v0',minval(v0),maxval(v0)
          !print *, 'omega0',minval(omega0),maxval(omega0)
          !print *, 'u',minval(u),maxval(u)
          !print *, 'v',minval(v),maxval(v)
          !print *, 'omega',minval(omega),maxval(omega)
          CALL advection(xi,omega0,u0,v0,omega,u,v,dummy,dt)
          print *, 'xi',minval(xi),maxval(xi)
          !xi = xi * dt
        ELSE
          xi(:,:) = 0.0
        ENDIF

! Add contribution to cost function (if cost)

! Save dependent variables (u,v,omega,xi) in Vx vector (if not cost)

        DO jvar = 1,varend
          indvar=var_ord(jvar)
          IF (var_jpk(indvar).NE.jpksize) GOTO 1000
          IF (var_jpt(indvar).NE.jptsize) GOTO 1000

          IF (var_nam(indvar).EQ.'U    ') THEN
            jx=var_sli_idx(indvar,jk,jt)
            CALL mk8vct(vectx(jx:),u(:,:),jk,jt, &
     &                  jvar,nbr,flagxyo)
          ENDIF

          IF (var_nam(indvar).EQ.'V    ') THEN
            jx=var_sli_idx(indvar,jk,jt)
            CALL mk8vct(vectx(jx:),v(:,:),jk,jt, &
     &                  jvar,nbr,flagxyo)
          ENDIF

          IF (var_nam(indvar).EQ.'OMEGA') THEN
            jx=var_sli_idx(indvar,jk,jt)
            CALL mk8vct(vectx(jx:),omega(:,:),jk,jt, &
     &                  jvar,nbr,flagxyo)
          ENDIF

          IF (var_nam(indvar).EQ.'XI   ') THEN
            jx=var_sli_idx(indvar,jk,jt)
            CALL mk8vct(vectx(jx:),xi(:,:),jk,jt, &
     &                  jvar,nbr,flagxyo)
          ENDIF

        ENDDO

      ENDDO
      ENDDO

! --- deallocation
      IF (allocated(adt)) deallocate(adt)
      IF (allocated(u)) deallocate(u)
      IF (allocated(v)) deallocate(v)
      IF (allocated(omega)) deallocate(omega)
      IF (allocated(xi)) deallocate(xi)
      IF (allocated(u0)) deallocate(u0)
      IF (allocated(v0)) deallocate(v0)
      IF (allocated(omega0)) deallocate(omega0)

      RETURN
!
! --- error management
!
 1000 CALL printerror2(0,1000,1,'utilconstraint','apply_constraint')
 1001 CALL printerror2(0,1001,3,'utilconstraint','apply_constraint')
!
      END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      END MODULE utilconstraint
