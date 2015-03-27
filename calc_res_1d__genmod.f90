        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 27 01:41:32 2015
        MODULE CALC_RES_1D__genmod
          INTERFACE 
            SUBROUTINE CALC_RES_1D(Q,RES)
              USE GRID_DIMEN
              REAL(KIND=8) :: Q(3,NX-1)
              REAL(KIND=8) :: RES(3,NX-1)
            END SUBROUTINE CALC_RES_1D
          END INTERFACE 
        END MODULE CALC_RES_1D__genmod
