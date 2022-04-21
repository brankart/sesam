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
      BIGREAL, PUBLIC, save :: dyn_constraint_std=1.0 ! Constraint error std

      ! Private variables
      LOGICAL, save :: first_call=.TRUE.
      BIGREAL, dimension(:,:), allocatable, save :: adt           ! absolute dynamic topography
      BIGREAL, dimension(:,:), allocatable, save :: u, u0         ! zonal velocity
      BIGREAL, dimension(:,:), allocatable, save :: v, v0         ! meridional velocity
      BIGREAL, dimension(:,:), allocatable, save :: omega, omega0 ! relative vorticity
      BIGREAL, dimension(:,:), allocatable, save :: xi            ! advection misfit

      CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! -----------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE apply_constraint(vectx,cost_jdyn)
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
      BIGREAL, optional, intent(out) :: cost_jdyn
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER :: flagxyo,ji,jj,jk,jt,jvar,indvar,jx,idx
      INTEGER :: jpisize,jpjsize,jpksize,jptsize
      INTEGER :: nbr,allocok
      BIGREAL :: dummy, dt, variance
      LOGICAL :: cost_eval
!----------------------------------------------------------------------
!
      IF (first_call.AND.(nprint.GE.1)) THEN
         WRITE(numout,*) '*** ROUTINE : sesam/.../aply_constraint :'
         WRITE(numout,*) '         apply dynamical constraint'
      ENDIF

      flagxyo=1
! Check dimensions : the full state vector must be recomposed
! before calling this routine
      IF (SIZE(vectx,1).NE.jpxend) GOTO 1000

! Initialize evaluation of cost function
      cost_eval = PRESENT(cost_jdyn)
      IF (cost_eval) cost_jdyn = 0.0_kr

! Check that grid is regular
      IF (MAXVAL(varngrd(1:varend)).GT.1) GOTO 1000

! Allocate 2D physical arrays
      jpisize=MAXVAL(var_jpi(1:varend))
      jpjsize=MAXVAL(var_jpj(1:varend))
      jpksize=MAXVAL(var_jpk(1:varend))
      jptsize=MAXVAL(var_jpt(1:varend))

      IF (first_call) THEN
        first_call=.FALSE.
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
      ENDIF

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
      DO jt=1,jptsize
      DO jk=1,jpksize

! Get free variable (adt)
        DO jvar = 1,varend
          indvar=var_ord(jvar)
          IF (var_jpk(indvar).NE.jpksize) GOTO 1000
          IF (var_jpt(indvar).NE.jptsize) GOTO 1000

          IF (var_nam(indvar).EQ.'ADT  ') THEN
            jx=var_sli_idx(indvar,jk,jt)
            CALL unmk8vct(vectx(jx:),adt(:,:),jk,jt, &
     &                    jvar,nbr,spval,flagxyo)
          ENDIF

        ENDDO

! Apply constraint
        u0(:,:) = u(:,:) ; v0(:,:) = v(:,:) ; omega0(:,:) = omega(:,:)
        !if (jproc==0) print *, 'ok5 adt',minval(adt),maxval(adt)
        CALL velocity(adt,u,v)
        CALL vorticity(u,v,omega)

        IF (jt.GT.1) THEN
          dt= 86400. ; dummy=0.
          !if (jproc==0) print *, 'ok6 u0',minval(u0),maxval(u0)
          !if (jproc==0) print *, 'ok6 v0',minval(v0),maxval(v0)
          !if (jproc==0) print *, 'ok6 omega0',minval(omega0),maxval(omega0)
          !if (jproc==0) print *, 'ok6 u',minval(u),maxval(u)
          !if (jproc==0) print *, 'ok6 v',minval(v),maxval(v)
          !if (jproc==0) print *, 'ok6 omega',minval(omega),maxval(omega)
          !xi(:,:) = 0.0
          CALL advection(xi,omega0,u0,v0,omega,u,v,dummy,dt)
          !if (jproc==0) print *, 'xi',jt,minval(xi),maxval(xi)
        ELSE
          xi(:,:) = 0.0
        ENDIF

! Add contribution to cost function (if cost)
        IF (cost_eval) THEN
          IF (jt.GT.1) THEN
            DO idx=jproc,jpisize*jpjsize-1,jpproc
              ji = 1 + MOD(idx,jpisize)
              jj = 1 + idx/jpisize
              cost_jdyn=cost_jdyn+xi(ji,jj)*xi(ji,jj)
            ENDDO
          ENDIF
        ENDIF

! Save dependent variables (u,v,omega,xi) in Vx vector (if not cost)
        IF (.NOT.cost_eval) THEN
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
        ENDIF

      ENDDO
      ENDDO

! Finalize evaluation of cost function
      IF (cost_eval) THEN
        variance = dyn_constraint_std * dyn_constraint_std
        cost_jdyn = 0.5_kr * cost_jdyn / variance
#if defined MPI
        CALL mpi_allreduce (mpi_in_place, cost_jdyn, 1,  &
     &       mpi_double_precision,mpi_sum,mpi_comm_world,mpi_code)
#endif
      ENDIF

 900     RETURN
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
