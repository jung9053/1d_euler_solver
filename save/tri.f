      subroutine tri(imin,imax,a,b,c,r)
      implicit none
      !..global
      integer :: imin,imax
      real(8) :: a(imax),b(imax),c(imax),r(imax)
      
      !..local
      integer :: i,j
      real(8) :: factor,z
      real(8),allocatable,dimension(:) :: x

      allocate(x(imax))

      if(b(1) == 0.0) then
         write(*,*) 'reorder the equation needed'
         stop
      endif

      do i=1,imax
        x(i) = 0.0
      enddo

      x(1) = c(1)/b(1)
      r(1) = r(1)/b(1)

      !forward sweep
      do i = imin+1,imax
         z = 1.0/(b(i)-a(i)*x(i-1))
         x(i) = c(i)*z
         r(i) = (r(i)-a(i)*r(i-1))*z         
      enddo

      !backward sweep
      do i=imin,imax-1
           j = imax-i
           r(j) = r(j)-x(j)*r(j+1)
      enddo
    
      deallocate(x)

      return
      end subroutine
