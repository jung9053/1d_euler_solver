        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 27 02:05:23 2015
        MODULE QFACE_1D__genmod
          INTERFACE 
            SUBROUTINE QFACE_1D(QF1,QF2,Q)
              USE GRID_DIMEN
              REAL(KIND=8) :: QF1(3,NX)
              REAL(KIND=8) :: QF2(3,NX)
              REAL(KIND=8) :: Q(3,NX-1)
            END SUBROUTINE QFACE_1D
          END INTERFACE 
        END MODULE QFACE_1D__genmod
