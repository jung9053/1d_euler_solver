      subroutine tri(imin,imax,a,b,c,r)
      implicit none
      !..global
      integer :: imin,imax
      real(8) :: a(imax),b(imax),c(imax),r(imax)
      
      !..local
      integer :: i,j
      real(8) :: factor,z,beta
      real(8),allocatable,dimension(:) :: u,gam

      allocate(gam(imax),u(imax))

      if(b(1) == 0.0) then
         write(*,*) 'reorder the equation needed'
         stop
      endif


      do i=1,imax
        u(i) = 0.0
        gam(i) = 0.0
      enddo


      beta  = b(1)
      u(1)  = r(1)/beta

     
      !forward sweep
      do i = 2,imax
         gam(i) = c(i-1)/beta
         beta = b(i)-a(i)*gam(i)
         if(beta==0) then
            write(*,*) 'fail in trisolver'
            stop
         endif

         u(i) = (r(i)-a(i)*u(i-1))/beta        
      enddo

      !backward sweep
      do i=1,imax-1
           j = imax-i
           u(j) = u(j)-gam(j+1)*u(j+1)
      enddo
    
      do i=1,imax
         r(i) = u(i)
      enddo
      deallocate(u,gam)

      return
      end subroutine
