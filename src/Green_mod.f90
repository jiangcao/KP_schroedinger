!
!   Green_mod.f90
!   Green Function solver module
!
!   Created by Jiang Cao on 12/7/14.
!   Copyright 2014 Jiang Cao. All rights reserved.
!
module green_mod
use matrix_mod, only : cMatrix, MUL_C, triMUL_c, size_cMatrix, trace_c, &
& invert_c, alloc2_cmatrix, free2_cmatrix, init2_cmatrix, alloc_cmatrix
use static
use mod_string, only : dbstring, string

implicit none

private

public :: rgf_backward

real(8) :: Dac, Dop_g
integer :: Nop_g, SCBA_max_iter, SCBA_tolerance

contains

!=================================================================================
!!                        Recursive Backward Green's solver                    
subroutine rgf_backward(En,mul,mur,TEMPl,TEMPr,Hii,H1i,Sii,sigma_lesser_ph,&          
                        sigma_r_ph,G_r,G_lesser,G_greater,Jdens,Gl,Gln,tr,tre)               
type(cMatrix),intent(in) :: Hii(:),H1i(:),Sii(:),sigma_lesser_ph(:),sigma_r_ph(:)
real(8),intent(in)       :: En,mul(:,:),mur(:,:),TEMPr(:,:),TEMPl(:,:)
type(cMatrix),intent(inout):: G_greater(:),G_lesser(:),G_r(:),Jdens(:),Gl(:),Gln(:)
real(8),intent(out)      :: tr,tre
!---------------------------------------------------------------------------------
integer    :: nx,M,ii,jj,kk,ll
complex(8) :: z
real(8)    :: tim
complex(8),allocatable :: sig(:,:),sig_in(:,:),sig_out(:,:),H00(:,:),H10(:,:)
complex(8),allocatable :: A(:,:),B(:,:),C(:,:),G00(:,:),GBB(:,:),sigmar(:,:),sigmal(:,:),GN0(:,:)
  nx = size(Hii)           ! <- lenght of the device
  z = dcmplx(En,0.0d0)
  !
  ! on the right contact
  ii = nx
  M  = size(Hii(ii)%m,1)
  allocate(H00(M,M))
  allocate(G00(M,M))
  allocate(GBB(M,M))
  allocate(sigmar(M,M))
  allocate(sig(M,M))
  !
  !!! H00 = H(i,i) + Sigma_ph(i) * S(i,i)
  call MUL_c(sigma_r_ph(ii)%m,Sii(ii)%m,'n','n',B)
  !
  H00 = Hii(ii)%m + B
  call sancho(M,En,Sii(ii)%m,H00,H1i(ii+1)%m,G00,GBB)
  !$omp critical
  open(unit = 10,file = 'sancho_g00.dat',position='append')
  write(10,*) En,2, -aimag(trace_c(G00))
  close(10)
  open(unit = 10,file = 'sancho_gbb.dat',position='append')
  write(10,*) En,2, -aimag(trace_c(Gbb))
  close(10)
  !$omp end critical
  !
  !!! Sigma_R = H(i,i+1) * G00 * H(i+1,i)
  !!! Gl(i) = [En*S(i,i) - H00 - Sigma_R]^-1
  call triMUL_c(H1i(ii+1)%m,G00,H1i(ii+1)%m,sigmar,'n','n','c')
  B = z*Sii(ii)%m - H00 - sigmar
  call invert_c(B)
  Gl(ii)%m = B
  !
  ! Gln(i) = Gl(i) * [Sigma_ph<(i)*S(i,i) + (-(Sigma_R - Sigma_R')*ferm(..))] * Gl(i)'
  call MUL_c(sigma_lesser_ph(ii)%m,Sii(ii)%m,'n','n',B)
  sig = - (sigmar - transpose(conjg(sigmar)))*ferm((En-mur)/(BOLTZ*TEMPr))
  !
  sig = sig + B
  call triMUL_c(Gl(ii)%m, sig, Gl(ii)%m, B, 'n', 'n', 'c')
  Gln(ii)%m = B
  deallocate(G00,GBB,sig)
  !
  allocate(A(M,M))
  ! inside device r -> l
  do ii = nx-1,2,-1
    M = size(Hii(ii)%m,1)
    if (size(H00,1).ne.M) then
    deallocate(H00,A)
    allocate(H00(M,M))
    allocate(A(M,M))
    end if
    call MUL_c(sigma_r_ph(ii)%m,Sii(ii)%m,'n','n',B)
    H00 = Hii(ii)%m + B
    !
    !!! H00 = H(i,i) + Sigma_ph(i) * S(i,i)
    !!! Gl(i) = [En*S(i,i) - H00 - H(i,i+1) * Gl(i+1) * H(i+1,i)]^-1
    call triMUL_c(H1i(ii+1)%m,Gl(ii+1)%m,H1i(ii+1)%m,B,'n','n','c')
    A = z*Sii(ii)%m - H00 - B
    call invert_c(A)
    Gl(ii)%m = A
    !
    !!! Gln(i) = Gl(i) * [Sigma_ph<(i)*S(i,i) + H(i,i+1)*Gln(i+1)*H(i+1,i)] * Gl(i)'
    call triMUL_c(H1i(ii+1)%m, Gln(ii+1)%m, H1i(ii+1)%m, B, 'n','n','c')
    call MUL_c(sigma_lesser_ph(ii)%m,Sii(ii)%m,'n','n',A)
    B = B + A
    call triMUL_c(Gl(ii)%m, B, Gl(ii)%m, A, 'n','n','c')
    Gln(ii)%m = A
  end do
  !
  ! on the left contact
  ii = 1
  M = size(Hii(ii)%m,1)
  allocate(H10(M,M))
  allocate(G00(M,M))
  allocate(GBB(M,M))
  allocate(sig(M,M))
  allocate(sigmal(M,M))
  if (size(H00,1).ne.M) then
    deallocate(H00)
    allocate(H00(M,M))
  end if
  !
  call MUL_c(sigma_r_ph(ii)%m,Sii(ii)%m,'n','n',B)
  H00 = Hii(ii)%m + B
  H10 = transpose(H1i(1)%m)
  !
  call sancho(M,En,Sii(ii)%m,H00,H10,G00,GBB)
  !
  call triMUL_c(H1i(1)%m,G00,H1i(1)%m,sigmal,'c','n','n')
  !
  !$omp critical
  open(unit = 10,file = 'sancho_g00.dat',position='append')
  write(10,*) En,1, -aimag(trace_c(G00))
  close(10)
  open(unit = 10,file = 'sancho_gbb.dat',position='append')
  write(10,*) En,1, -aimag(trace_c(Gbb))
  close(10)
  !$omp end critical
  !
  call triMUL_c(H1i(2)%m, Gl(2)%m, H1i(2)%m, B, 'n','n','c')
  A = z*Sii(1)%m - H00 - B - sigmal
  call invert_c(A)
  !
  G_r(1)%m = A
  Gl(1)%m = A
  !
  !!! Sigma^< = Sigma_11^< + Sigma_ph^< + Sigma_s^<
  call triMUL_c(H1i(2)%m, Gln(2)%m, H1i(2)%m, B ,'n','n','c')
  call MUL_c(sigma_lesser_ph(1)%m,Sii(1)%m,'n','n',A )
  sig = -(sigmal - transpose(conjg(sigmal))) * ferm((En-mul)/(BOLTZ*TEMPl))
  sig = sig + A + B
  !
  !!! G^< = G * Sigma^< * G'
  call triMUL_c(G_r(1)%m, sig, G_r(1)%m, B, 'n','n','c')
  !
  G_lesser(1)%m = B
  G_greater(1)%m = G_lesser(1)%m + (G_r(1)%m - transpose(conjg(G_r(1)%m)))
  !
  A = - (sigmal - transpose(conjg(sigmal))) * ferm((En-mul)/(BOLTZ*TEMPl))
  call MUL_c(A,G_greater(1)%m,'n','n',B)
  A = - (sigmal - transpose(conjg(sigmal))) * (ferm((En-mul)/(BOLTZ*TEMPl))-1.0d0)
  call MUL_c(A,G_lesser(1)%m,'n','n',C)
  !
  Jdens(1)%m = B - C
  !
  tim = 0.0d0
  do jj = 1,M
    tim = tim + dble(Jdens(1)%m(jj,jj))
  end do
  tre = tim
  deallocate(sigmal,sig,G00,GBB,H10)
  allocate(GN0(M,M))
  !
  ! inside device l -> r
  do ii = 2,nx
    M = size(Hii(ii)%m,1)
    !!! A = G<(i-1) * H(i-1,i) * Gl(i)' + G(i-1) * H(i-1,i) * Gln(i)
    call triMUL_c(G_lesser(ii-1)%m,H1i(ii)%m,Gl(ii)%m, A, 'n','n','c')
    call triMUL_c(G_r(ii-1)%m,H1i(ii)%m,Gln(ii)%m, B, 'n','n','n')
    A = A + B
    !!! B = H(i,i-1) * A
    !!! Jdens(i) = -2 * re(B)
    call MUL_c(H1i(ii)%m, A, 'c', 'n', B)
    Jdens(ii)%m = - 2.0d0 * dble(B(:,:))
    !
    !!! GN0 = Gl(i) * H(i,i-1) * G(i-1)
    !!! G(i) = Gl(i) + GN0 * H(i-1,i) * Gl(i)
    call MUL_c(Gl(ii)%m, H1i(ii)%m, 'n', 'c', B)
    call MUL_c(B , G_r(ii-1)%m, 'n', 'n', GN0)
    call MUL_c(GN0,H1i(ii)%m, 'n', 'n', C )
    call MUL_c(C, Gl(ii)%m, 'n', 'n', A)
    G_r(ii)%m = Gl(ii)%m + A
    !
    !!! G<(i) = Gln(i) + Gl(i) * H(i,i-1) * G<(i-1) * H(i-1,i) *Gl(i)'
    call MUL_c(Gl(ii)%m, H1i(ii)%m, 'n', 'c', B)
    call MUL_c(B, G_lesser(ii-1)%m, 'n', 'n', C)
    call MUL_c(C, H1i(ii)%m, 'n', 'n', A)
    call MUL_c(A, Gl(ii)%m, 'n','c',C)
    G_lesser(ii)%m = Gln(ii)%m + C
    !
    !!! G<(i) = G<(i) + GNO * H(i-1,i) * Gln(i)
    call MUL_c(GN0,H1i(ii)%m, 'n', 'n', B)
    call MUL_c(B, Gln(ii)%m, 'n', 'n', C)
    G_lesser(ii)%m = G_lesser(ii)%m + C
    !
    !!! G<(i) = G<(i) + Gln(i) * H(i,i-1) * GN0
    call MUL_c(Gln(ii)%m,H1i(ii)%m,'n','c',B)
    call MUL_c(B, GN0, 'n', 'c', C)
    G_lesser(ii)%m = G_lesser(ii)%m + C
    !
    !!! G>(i) = G<(i) + (G(i) - G(i)')
    G_greater(ii)%m = G_lesser(ii)%m + (G_r(ii)%m - transpose(conjg(G_r(ii)%m)))
  end do
  ii = nx
  ! on the right contact
  A = -(sigmar-transpose(conjg(sigmar)))*ferm((En-mur)/(BOLTZ*TEMPr))
  call MUL_c(A, G_greater(ii)%m,'n','n',B)
  A = -(sigmar-transpose(conjg(sigmar)))*(ferm((En-mur)/(BOLTZ*TEMPr))-1.0d0)
  call MUL_c(A,G_lesser(ii)%m, 'n', 'n',C)
  tim = 0.0d0
  do jj = 1,M
    tim = tim + dble(B(jj,jj) - C(jj,jj))
  end do
  tr = tim
  deallocate(B,A,C,GN0,sigmar)
  !
end subroutine rgf_backward






!=================================================================================
!!                      Fermi distribution function                               
elemental Function ferm(a)                                                       
Real (8),intent(in) ::  a
real (8) :: ferm
  ferm=1.0d0/(1.0d0+Exp(a))
End Function ferm


!=================================================================================
!!                           Sancho-Rubio                                
subroutine sancho(nm,E,S00,H00,H10,G00,GBB)                                      
integer i,j,k,nm,nmax
COMPLEX(8) :: z
real(8) :: E,error
REAL(8) :: TOL=1.0D-100  ! [eV]
COMPLEX(8), INTENT(IN) ::  S00(nm,nm), H00(nm,nm), H10(nm,nm)
COMPLEX(8), INTENT(OUT) :: G00(nm,nm), GBB(nm,nm)
COMPLEX(8), ALLOCATABLE :: A(:,:), B(:,:), C(:,:), tmp(:,:)
COMPLEX(8), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:), Id(:,:)
!COMPLEX(8), ALLOCATABLE :: WORK(:)
COMPLEX(8), EXTERNAL :: ZLANGE
  !
  Allocate( H_BB(nm,nm) )
  Allocate( H_SS(nm,nm) )
  Allocate( H_01(nm,nm) )
  Allocate( H_10(nm,nm) )
  Allocate( Id(nm,nm) )
  Allocate( A(nm,nm) )
  Allocate( B(nm,nm) )
  Allocate( C(nm,nm) )
  Allocate( tmp(nm,nm) )
  nmax=50
  z = dcmplx(E,1.0d-3)
  Id=0.0d0
  tmp=0.0d0
  do i=1,nm
      Id(i,i)=1.0d0
      tmp(i,i)=dcmplx(0.0d0,1.0d0)
  enddo
  H_BB = H00
  H_10 = H10
  H_01 = TRANSPOSE( CONJG( H_10 ) )
  H_SS = H00
  do i = 1, nmax
      A = z*S00 - H_BB
      call invert_c(A)
      call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_10,nm,beta,B,nm)
      call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm)
      H_SS = H_SS + C
      H_BB = H_BB + C
      call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,C,nm)
      call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_01,nm,beta,B,nm)
      call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,A,nm)
      H_10 = C
      H_BB = H_BB + A
      call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm)
      H_01 = C
      ! NORM --> inspect the diagonal of A
      error=0.0d0
      DO k=1,nm
          DO j=1,nm
              error=error+sqrt(aimag(C(k,j))**2+Dble(C(k,j))**2)
          END DO
      END DO
      !write(90,*)E,i,error
      tmp=H_SS
      IF ( abs(error) < TOL ) THEN
          !write(90,*) 'SR: Exited, abs(error)=',i,abs(error)
          EXIT
      ELSE
      END IF
      IF (i .EQ. nmax) THEN
          write(*,*) 'SEVERE warning: nmax reached in sancho!!!',error
      END IF
  enddo
  G00 = z*S00 - H_SS
  call invert_c(G00)
  !
  GBB = z*S00 - H_BB
  call invert_c(GBB)
  !
  Deallocate( tmp )
  Deallocate( A )
  Deallocate( B )
  Deallocate( C )
  Deallocate( H_BB )
  Deallocate( H_SS )
  Deallocate( H_01 )
  Deallocate( H_10 )
  Deallocate( Id )
end subroutine sancho


end module green_mod
