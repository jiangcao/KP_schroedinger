module math

implicit none

contains

function gaussian(r,b,U0,sigma)
real(8)::gaussian
real(8),intent(in)::r,b,U0,sigma

gaussian = U0 * exp( - (r-b)**2 / 2.0d0 / sigma**2 ) - U0/2

end function gaussian

end module math
