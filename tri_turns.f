      subroutine tri(imin,imax,a,b,c,r)
      implicit none
      !..global
      integer :: imin,imax
      real(8) :: a(imax),b(imax),c(imax),r(imax)
      
      !..local
      integer :: i
      real(8) :: factor
    

      do i = imin,imax-3
        
        write(*,*) i, b(i)
        if (b(i) .eq. 0.0) then
            write(*,*) "1:error at tridiagonal solver(zero diagonal)"
            stop
            return
        endif

        factor = a(i+1) / b(i)
        b(i+1) = b(i+1) - factor * c(i)
        r(i+1) = r(i+1) - factor * r(i)
      enddo
     
      if (b(imax-2).eq.0.0) THEN
        write(*,*) "2:error at tridiagonal solver(zero diagonal)"
        stop
        return
      endif

      factor = a(imax-1) / b(imax-2)
      b(imax-1) = b(imax-1) - factor * c(imax-2)
      r(imax-1) = r(imax-1) - factor * r(imax-2)
      
      factor = a(imax) / b(imax-1)
      b(imax) = b(imax) - factor * c(imax-1)
      r(imax) = r(imax) - factor * r(imax-1)
      
      r(imax) = r(imax) / b(imax)
      do i = imax-1, imin, -1
        r(i) = (r(i) - c(i)*r(i+1)) / b(i)
      enddo
     

      return
      end subroutine
