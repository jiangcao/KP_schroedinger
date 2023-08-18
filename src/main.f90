! Copyright (c) 2023 Jiang Cao, ETH Zurich 
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors 
!    may be used to endorse or promote products derived from this software without 
!    specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE. 
!
PROGRAM main
use kpHam
implicit none

integer :: nbnd
real(8)::gam1,gam2,gam3,delta,kappa,dd(3),kpt(3),Ek(6),kpt1(3),kpt2(3)
complex(8),allocatable::KPcoeff(:,:,:),H00(:,:),H10(:,:)
integer :: nk, ny, nz, direction
integer :: ik , num_cpu

open(unit=10,file='input',status='unknown')
 read(10,*) nbnd
 read(10,*) gam1
 read(10,*) gam2
 read(10,*) gam3
 read(10,*) kappa
 read(10,*) delta
 read(10,*) direction
 read(10,*) Ny,Nz
 read(10,*) dd
 read(10,*) num_cpu
close(10)

call omp_set_num_threads(num_cpu)

allocate(KPcoeff(nbnd,nbnd,10))
print *, 'nbnd=',nbnd
if (nbnd==4) then
    call generate_LK4(KPcoeff, gam1,gam2,gam3,kappa)
end if

if (nbnd==6) then
    call generate_LK6(KPcoeff, gam1,gam2,gam3,delta,kappa)
end if

call save_KP_coeff('kpcoeff.dat', nbnd, KPcoeff)

call save_Ham_blocks('Ham_blocks.dat', nbnd, KPcoeff,dd)

call test_bandstructure(nbnd,KPcoeff,dd)

allocate(H00(Ny*Nz*nbnd,Ny*Nz*nbnd))
allocate(H10(Ny*Nz*nbnd,Ny*Nz*nbnd))

call build_wire_ham_blocks(nbnd,KPcoeff, dd, direction, Ny, Nz, H00, H10)

! call save_matrix('Hii.dat', Ny*Nz*nbnd, H00)
! call save_matrix('H1i.dat', Ny*Nz*nbnd, H10)

call save_matrix_csr('Hii_csr.dat', Ny*Nz*nbnd, Ny*Nz*nbnd, H00)
call save_matrix_csr('H1i_csr.dat', Ny*Nz*nbnd, Ny*Nz*nbnd, H10)

call test_transport_bandstructure(H00,H10,nbnd,ny,nz,dd(direction))

deallocate(KPcoeff)
deallocate(H10,H00)

END PROGRAM main
