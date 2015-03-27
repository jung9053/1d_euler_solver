        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 27 01:41:32 2015
        MODULE EULER_FWD_1D__genmod
          INTERFACE 
            SUBROUTINE EULER_FWD_1D(Q,RES,DELX,DELT)
              USE GRID_DIMEN
              REAL(KIND=8) :: Q(3,NX-1)
              REAL(KIND=8) :: RES(3,NX-1)
              REAL(KIND=8) :: DELX(NX-1)
              REAL(KIND=8) :: DELT(NX-1)
            END SUBROUTINE EULER_FWD_1D
          END INTERFACE 
        END MODULE EULER_FWD_1D__genmod
