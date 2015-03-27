      subroutine muscl_recon(q,qf1,qf2)
      use param
      use grid_dimen
      implicit none
      !..global
      real(8) :: q(nvar,nx-1)
      real(8) :: qf1(nvar,nx),qf2(nvar,nx)
      !..local
      real(8) :: thm,thp,q2i,q2i1,a1,a2,q3i,q3qt,eps
      !..local array
      real(8),allocatable,dimension(:,:) :: q2
      integer :: i,j
     
      allocate(q2(nvar,nx-1))

      thm = 1.0-1.0/3.0
      thp = 1.0+1.0/3.0
      eps = 1e-10 

      do j=1,nvar

        do i=1,nx-1
          q2(j,i) = q(j,i+1) - q(j,i)
        enddo

        do i=2,nx-2
          q2i  = q2(j,i)  
          q2i1 = q2(j,i-1)
          a1   = 3.0*(q2i*q2i1+eps)
          a2   = 2.0*(q2i-q2i1)**2 + a1
          q3i  = a1/a2
          q3qt = 0.25*q3i
          qf1(j,i+1) = q(j,i) + q3qt*( thm*q2i1 + thp*q2i)
          qf2(j,i) = q(j,i) - q3qt*( thm*q2i1 + thp*q2i)
        enddo

        qf1(j,2)    = q(j,1)
        qf1(j,nx)   = q(j,nx-1)
        qf2(j,1)    = q(j,1)
        qf2(j,nx-1) = q(j,nx-1)
      enddo
      deallocate(q2)

      return
      end subroutine muscl_recon
