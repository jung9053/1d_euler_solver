        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 27 09:30:05 2015
        MODULE WENO5_RECON__genmod
          INTERFACE 
            SUBROUTINE WENO5_RECON(Q,QF1,QF2)
              USE GRID_DIMEN
              REAL(KIND=8) :: Q(3,NX-1)
              REAL(KIND=8) :: QF1(3,NX)
              REAL(KIND=8) :: QF2(3,NX)
            END SUBROUTINE WENO5_RECON
          END INTERFACE 
        END MODULE WENO5_RECON__genmod
