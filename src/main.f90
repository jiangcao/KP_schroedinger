PROGRAM main
use kpHam
implicit none

integer :: nbnd
real(8)::gam1,gam2,gam3,delta,kappa,dd(3),kpt(3),Ek(6),kpt1(3),kpt2(3)
complex(8),allocatable::KPcoeff(:,:,:),H00(:,:),H10(:,:)
integer :: nk, ny, nz, direction
integer :: ik

open(unit=10,file='input',status='unknown')
 read(10,*) nbnd
 read(10,*) gam1
 read(10,*) gam2
 read(10,*) gam3
 read(10,*) kappa
 read(10,*) delta
 read(10,*) direction
 read(10,*) Ny,Nz
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

allocate(H00(Ny*Nz*nbnd,Ny*Nz*nbnd))
allocate(H10(Ny*Nz*nbnd,Ny*Nz*nbnd))

call build_wire_ham_blocks(nbnd,KPcoeff, dd, direction, Ny, Nz, H00, H10)

call save_matrix('Hii.dat', Ny*Nz*nbnd, H00)
call save_matrix('H1i.dat', Ny*Nz*nbnd, H10)

call test_transport_bandstructure(H00,H10,nbnd,ny,nz,dd(direction))

deallocate(KPcoeff)
deallocate(H10,H00)

END PROGRAM main
