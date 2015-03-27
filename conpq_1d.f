      subroutine conpq_1d(n,q)
      use param
      use flow_property
      implicit none
      integer :: n
      real(8) :: q(nvar,n)
      !local
      integer :: i

      do i=1,n
         q(3,i) = q(3,i)*gm1i + 0.5d0*q(1,i)*(q(2,i)*q(2,i))
         q(2,i) = q(2,i)*q(1,i)
      end do

      return
      end subroutine conpq_1d

      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine conqp_1d(n,q)
      use param
      use flow_property
      implicit none
      integer :: n
      real(8) :: q(nvar,n)
      !local
      integer :: i

      do i=1,n
         q(2,i) = q(2,i)/q(1,i)
         q(3,i) = gm1*(q(3,i) - 0.5d0*q(1,i)*(q(2,i)*q(2,i)))
      end do

      return
      end subroutine conqp_1d

      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine conpq1_1d(q)
      use param
      use flow_property
      implicit none
      real(8) :: q(nvar)

      q(3) = q(3)*gm1i + 0.5d0*q(1)*(q(2)*q(2))
      q(2) = q(2)*q(1)

      return
      end subroutine conpq1_1d

      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine conqp1_1d(q)
      use param
      use flow_property
      implicit none
      real(8) :: q(nvar)

      q(2) = q(2)/q(1)
      q(3) = gm1*(q(3) - 0.5d0*q(1)*(q(2)*q(2)))

      return
      end subroutine conqp1_1d

