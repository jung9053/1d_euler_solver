        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 27 01:41:32 2015
        MODULE CALC_INV_FLUX_1D__genmod
          INTERFACE 
            SUBROUTINE CALC_INV_FLUX_1D(QF1,QF2,RES)
              USE GRID_DIMEN
              REAL(KIND=8) :: QF1(3,NX)
              REAL(KIND=8) :: QF2(3,NX)
              REAL(KIND=8) :: RES(3,NX-1)
            END SUBROUTINE CALC_INV_FLUX_1D
          END INTERFACE 
        END MODULE CALC_INV_FLUX_1D__genmod
