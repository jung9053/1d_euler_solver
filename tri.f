      subroutine tri(imin,imax,a,b,c,r)
      implicit none
      !..global
      integer :: imin,imax
      real(8) :: a(imax),b(imax),c(imax),r(imax)
      
      !..local
      integer :: i,j
      real(8) :: factor,z

      do i=imin,imax-3
         factor = a(i+1) /b(i)
         b(i+1) = b(i+1) - factor*c(i)
         r(i+1) = r(i+1) - factor*r(i)
      enddo

      factor = a(imax-1)/b(imax-2)
      b(imax-1) = b(imax-1) -factor*c(imax-2)
      r(imax-1) = r(imax-1) -factor*r(imax-2)

      factor = a(imax)/b(imax-1)
      b(imax) = b(imax) - factor *c(imax-1)
      r(imax) = r(imax) - factor *r(imax-1)
      
      r(imax) = r(imax)/b(imax)
      do i=imax-1,imin,-1
         r(i) = (r(i) - c(i)*r(i+1))/b(i)
      enddo

      return
      end subroutine
