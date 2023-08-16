PROGRAM main
use kpHam
implicit none

integer :: nbnd
real(8)::gam1,gam2,gam3,delta,kappa,dd(3),kpt(3),Ek(6),kpt1(3),kpt2(3)
complex(8),allocatable::KPcoeff(:,:,:)
integer :: nk
integer :: ik

open(unit=10,file='input',status='unknown')
 read(10,*) nbnd
 read(10,*) gam1
 read(10,*) gam2
 read(10,*) gam3
 read(10,*) kappa
 read(10,*) delta
close(10)

allocate(KPcoeff(nbnd,nbnd,10))
print *, 'nbnd=',nbnd
if (nbnd==4) then
    call generate_LK4(KPcoeff, gam1,gam2,gam3,kappa)
end if

if (nbnd==6) then
    call generate_LK6(KPcoeff, gam1,gam2,gam3,delta,kappa)
end if

call save_KP_coeff('kpcoeff.dat', nbnd, KPcoeff)


dd = (/ 1.0d-9 , 1.0d-9 , 2.5d-10 /)

call save_Ham_blocks('Ham_blocks.dat', nbnd, KPcoeff,dd)

call test_bandstructure(nbnd,KPcoeff,dd)

deallocate(KPcoeff)

END PROGRAM main
