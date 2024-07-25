module math

implicit none
complex(8),parameter::czero=dcmplx(0.0d0,0.0d0)
complex(8),parameter::cone=dcmplx(1.0d0,0.0d0)
contains

   function gaussian(r,b,U0,sigma)
      real(8)::gaussian
      real(8),intent(in)::r,b,U0,sigma

      gaussian = U0 * exp( - (r-b)**2 / 2.0d0 / sigma**2 ) - U0/2

   end function gaussian

   ! complex identity matrix
   function zeye(n)
      complex(8)::zeye(n,n)
      integer,intent(in)::n
      integer::i
      zeye=czero
      do i=1,n
         zeye(i,i)=cone
      enddo
   end function zeye
   
   ! matrix inversion
    subroutine invert(A, nn)
        integer :: info, nn
        integer, dimension(:), allocatable :: ipiv
        complex(8), dimension(nn, nn), intent(inout) :: A
        complex(8), dimension(:), allocatable :: work
        allocate (work(nn*nn))
        allocate (ipiv(nn))
        call zgetrf(nn, nn, A, nn, ipiv, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetrf failed, info=', info
            A = czero
        else
            call zgetri(nn, A, nn, ipiv, work, nn*nn, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetri failed, info=', info
                A = czero
            end if
        end if
    end subroutine invert
    
    subroutine normalize_real_matrix(A,nn)
      real(8), dimension(nn, nn), intent(inout) :: A
      real(8)::tmp
      integer,intent(in):: nn
      integer::i
      do i = 1,nn
         tmp = sqrt(sum( abs(A(:,i))**2 ))
         A(:,i) = A(:,i) / tmp
      enddo
    end subroutine normalize_real_matrix
      
   

end module math
