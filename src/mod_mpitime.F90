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
!************************************************************************
!                                                                       !
!   Module : MPI_Time           (Version : 1.9)                         !
!                                                                       !
!   Goal   : Measure and print elapsed and cpu user times               !
!                                                                       !
!   Usage  : Insert a "USE MPI_TIME" instruction inside the MPI Fortran !
!            program unit to instrument, then make calls to the         !
!            MPI_TIMER subroutine as shown in the example below :       !
!                                                                       !
!           PROGRAM foo                                                 !
!             USE MPI_TIME                                              !
!             ...                                                       !
!             CALL MPI_INIT(ierr)                                       !
!                                                                       !
!             !... Set elapsed and CPU user times                       !
!             CALL MPI_TIMER(0)                                         !
!             ...                                                       !
!             ... Instruction block to instrument ...                   !
!             ...                                                       !
!             !... Measure and print elapsed and CPU user times         !
!             CALL MPI_TIMER(1)                                         !
!                                                                       !
!             CALL MPI_FINALIZE(ierr)                                   !
!           END PROGRAM foo                                             !
!                                                                       !
!   Notes   : 1) MPI_TIMER subroutine is collective over all processes  !
!                of MPI_COMM_WORLD communicator.                        !
!             2) User must be aware that, on some machines, default     !
!                CPU user time includes also MPI wait on communication  !
!                to complete.                                           !
!                                                                       !
!   Output : At normal termination, MPI process of rank 0 prints        !
!            on stdout the elapsed and cpu user times of all processes. !
!                                                                       !
!   Author : Jalel Chergui                                              !
!                                                                       !
! Permission is garanted to copy and distribute this file or modified   !
! versions of this file for no fee, provided the copyright notice and   !
! this permission notice are preserved on all copies.                   !
!                                                                       !
! (C) Mars 2002, Jalel.Chergui@idris.fr, CNRS/IDRIS, France.            !
!************************************************************************

      MODULE MOD_MPITIME
        IMPLICIT NONE
        PRIVATE

!       !... Common Timing parameters
        INTEGER, PARAMETER         :: p = SELECTED_REAL_KIND(12)  
        REAL(kind=p)               :: Eoverhead, Coverhead
        REAL(kind=p), DIMENSION(2) :: Etime, Ctime

#if defined MPI
        PUBLIC :: MPI_Timer

        CONTAINS

        SUBROUTINE MPI_Timer(flag)
          IMPLICIT NONE

!         !... Header files
          INCLUDE "mpif.h"

!         !... Input dummy parameter
          INTEGER, INTENT(IN) :: flag

!         !... Local variables
          INTEGER                                 :: rank, nb_procs, i, code
          INTEGER, ALLOCATABLE, DIMENSION(:)      :: All_Rank
          REAL(KIND=p), ALLOCATABLE, DIMENSION(:) :: All_Etime, All_Ctime
          REAL(KIND=p)                            :: Total_Etime, Max_Etime, &
     &                                               Min_Etime, Avg_Etime, &
     &                                               Total_Ctime, Max_Ctime, &
     &                                               Min_Ctime, Avg_Ctime, dummy
          CHARACTER(LEN=80), dimension(6) :: lignes
          CHARACTER(LEN=80)               :: hline, lasthline        
          CHARACTER(LEN=960)              :: FMT

          SELECT CASE(flag)
            CASE(0)
              Eoverhead = MPI_WTIME()
              Eoverhead = MPI_WTIME() - Eoverhead
              CALL CPU_TIME(dummy)
              CALL CPU_TIME(Coverhead)
              if (dummy < 0.0_p) &
     &          PRINT *,"Warning, MPI_TIMER: CPU user time is not available on this machine."
              Coverhead = Coverhead - dummy
!             !... Start elapsed and CPU time counters
              Etime(1) = MPI_WTIME()
              CALL CPU_TIME(Ctime(1))
      
      CASE(1)
!             !... Final CPU and elapsed times
              CALL CPU_TIME(Ctime(2))
              Etime(2) = MPI_WTIME() - Etime(1) - Eoverhead - Coverhead
              Ctime(2) = Ctime(2) - Ctime(1) - Coverhead
              CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)
              CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
              IF ( rank == 0) ALLOCATE(All_Etime(nb_procs), &
     &                                 All_Ctime(nb_procs), &
     &                                 All_Rank(nb_procs) )
              CALL MPI_GATHER(Etime(2), 1, MPI_DOUBLE_PRECISION, &
     &                        All_Etime, 1, MPI_DOUBLE_PRECISION, &
     &                        0, MPI_COMM_WORLD, code)
              CALL MPI_GATHER(Ctime(2), 1, MPI_DOUBLE_PRECISION, &
     &                        All_Ctime, 1, MPI_DOUBLE_PRECISION, &
     &                        0, MPI_COMM_WORLD, code)
              IF ( rank == 0) THEN
                All_Rank(:) = (/ (i,i=0,nb_procs-1) /)
      
!               !... Compute elapse user time
                Total_Etime = SUM(All_Etime(:))
                Avg_Etime   = Total_Etime/REAL(nb_procs,KIND=p)
                Max_Etime   = MAXVAL(All_Etime(:))
                Min_Etime   = MINVAL(All_Etime(:))
                IF( Min_Etime <= 0.0_p ) THEN
                  PRINT *,"Warning, MPI_TIMER: Measured elapsed user time seems to be too"
                  PRINT *,"                    short compared to the clock precision.    "
                  PRINT *,"                    Timings could be erroneous.               "
                END IF

!               !... Compute CPU user time
                Total_Ctime = SUM(All_Ctime(:))
                Avg_Ctime   = Total_Ctime/REAL(nb_procs,KIND=p)
                Max_Ctime   = MAXVAL(All_Ctime(:))
                Min_Ctime   = MINVAL(All_Ctime(:))
                IF( Min_Ctime <= 0.0_p ) THEN
                  PRINT *,"Warning, MPI_TIMER: Measured CPU user time seems to be too"
                  PRINT *,"                    short compared to the clock precision."
                  PRINT *,"                    Timings could be erroneous."
                END IF

!               !... Output Format
                hline   ='10X,18("-"),"|",19("-"),"|",16("-"),/,'
                lasthline='10X,18("-"),"|",19("-"),"|",16("-"),//)'
                lignes(1)='(//,10X,"Process Rank",6(" "),"|","Elapsed time (sec.)|","CPU Time (sec.)",/,'
                lignes(2)='    (10X,I4,14(" "),"|",F12.3,7(" "),"|",F12.3,/),'
                WRITE(lignes(2)(1:4),'(I4)') nb_procs
                lignes(3)='10X,"Total time (sec.) |",F12.3,7(" "),"|",F12.3,/,'
                lignes(4)='10X,"Min.  time (sec.) |",F12.3,7(" "),"|",F12.3,/,'
                lignes(5)='10X,"Max.  time (sec.) |",F12.3,7(" "),"|",F12.3,/,'
                lignes(6)='10X,"Avg.  time (sec.) |",F12.3,7(" "),"|",F12.3,/,'
                FMT=TRIM(lignes(1))//TRIM(hline)//TRIM(lignes(2))//TRIM(hline) &
     &            //TRIM(lignes(3))//TRIM(hline)//TRIM(lignes(4))//TRIM(hline) &
     &            //TRIM(lignes(5))//TRIM(hline)//TRIM(lignes(6))//TRIM(lasthline)
                WRITE(*, TRIM(FMT)) &
     &              (All_rank(i),All_Etime(i),All_Ctime(i),i=1, nb_procs), &
     &              Total_Etime, Total_Ctime, Min_Etime, Min_Ctime, &
     &              Max_Etime, Max_Ctime, Avg_Etime, Avg_Ctime
                DEALLOCATE(All_Etime, All_Ctime, All_rank)
              END IF

            CASE DEFAULT
              PRINT *,"Error, MPI_TIMER: Invalid input parameter"

          END SELECT
        END SUBROUTINE MPI_Timer
#endif
      END MODULE MOD_MPITIME
