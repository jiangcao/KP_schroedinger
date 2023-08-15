PROGRAM main
use kpHam
implicit none

integer :: nbnd
real(8)::gam1,gam2,gam3,delta,kappa,dd(3),kpt(3),Ek(4),kpt1(3),kpt2(3)
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

if (nbnd==4) then
    call generate_LK4(KPcoeff, gam1,gam2,gam3,kappa)
end if

if (nbnd==6) then
    call generate_LK6(KPcoeff, gam1,gam2,gam3,delta,kappa)
end if

call save_KP_coeff('kpcoeff.dat', nbnd, KPcoeff)

dd = (/ 1.0d-10 , 1.0d-10 , 1.0d-10 /)

nk = 1000

open(unit=10,file='ek0.dat',status='unknown')
open(unit=11,file='ek.dat',status='unknown')

kpt1 = (/ -1.0d9 , 0.0d0 , 0.0d0 /)
kpt2 = (/ +1.0d9 , 0.0d0 , 0.0d0 /)

do ik = 1,nk

  kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1
  
  call calc_bands_at(nbnd,KPcoeff,dd,kpt,Ek)
  
  write(11,'(5E18.6)') dble(ik)/dble(nk), Ek
  
  call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)
  
  write(10,'(5E18.6)') dble(ik)/dble(nk), Ek
  
enddo

kpt1 = (/ 0.0d9 , -1.0d9 , 0.0d0 /)
kpt2 = (/ 0.0d9 ,  1.0d9 , 0.0d0 /)

do ik = 1,nk

  kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1
  
  call calc_bands_at(nbnd,KPcoeff,dd,kpt,Ek)
  
  write(11,'(5E18.6)') dble(ik)/dble(nk)+1.0d0, Ek
  
  call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)
  
  write(10,'(5E18.6)') dble(ik)/dble(nk)+1.0d0, Ek
  
enddo

kpt1 = (/ 0.0d9 , 0.0d9 , -1.0d9 /)
kpt2 = (/ 0.0d9 , 0.0d9 ,  1.0d9 /)

do ik = 1,nk

  kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1
  
  call calc_bands_at(nbnd ,KPcoeff,dd,kpt,Ek)
  
  write(11,'(5E18.6)') dble(ik)/dble(nk)+2.0d0, Ek
  
  call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)
  
  write(10,'(5E18.6)') dble(ik)/dble(nk)+2.0d0, Ek
  
enddo

kpt1 = (/ -1.0d9 , 0.0d9 , -1.0d9 /)
kpt2 = (/  1.0d9 , 0.0d9 ,  1.0d9 /)

do ik = 1,nk

  kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1
  
  call calc_bands_at(nbnd ,KPcoeff,dd,kpt,Ek)
  
  write(11,'(5E18.6)') dble(ik)/dble(nk)+3.0d0, Ek
  
  call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)
  
  write(10,'(5E18.6)') dble(ik)/dble(nk)+3.0d0, Ek
  
enddo

kpt1 = (/ -1.0d9 , -1.0d9 , 0.0d9 /)
kpt2 = (/  1.0d9 ,  1.0d9 , 0.0d9 /)

do ik = 1,nk

  kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1
  
  call calc_bands_at(nbnd ,KPcoeff,dd,kpt,Ek)
  
  write(11,'(5E18.6)') dble(ik)/dble(nk)+4.0d0, Ek
  
  call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)
  
  write(10,'(5E18.6)') dble(ik)/dble(nk)+4.0d0, Ek
  
enddo


kpt1 = (/ -1.0d9 , -1.0d9 , -1.0d9 /)
kpt2 = (/  1.0d9 ,  1.0d9 ,  1.0d9 /)

do ik = 1,nk

  kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1
  
  call calc_bands_at(nbnd ,KPcoeff,dd,kpt,Ek)
  
  write(11,'(5E18.6)') dble(ik)/dble(nk)+5.0d0, Ek
  
  call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)
  
  write(10,'(5E18.6)') dble(ik)/dble(nk)+5.0d0, Ek
  
enddo


close(10)
close(11)

deallocate(KPcoeff)

END PROGRAM main
