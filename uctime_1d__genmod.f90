        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 27 16:11:55 2015
        MODULE UCTIME_1D__genmod
          INTERFACE 
            SUBROUTINE UCTIME_1D(Q,DELX,DELT,RTIME,DTMIN)
              USE GRID_DIMEN
              REAL(KIND=8) :: Q(3,NX-1)
              REAL(KIND=8) :: DELX(NX-1)
              REAL(KIND=8) :: DELT(NX-1)
              REAL(KIND=8) :: RTIME
              REAL(KIND=8) :: DTMIN
            END SUBROUTINE UCTIME_1D
          END INTERFACE 
        END MODULE UCTIME_1D__genmod
