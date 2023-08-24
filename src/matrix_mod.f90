!
!   matrix_mod.f90
!   matrix module
!
!   Created by Jiang Cao on 12/8/14.
!   Copyright 2014 Jiang Cao. All rights reserved.
!
module matrix_mod
use mod_string
implicit none

! complex matrix type
type cMatrix
    complex(8),allocatable :: m(:,:)
end type

contains

!***************************! Matrix list operations !***************************!

!=================================================================================
!                        save trace of each matrix into file                     
subroutine print_cMatrix(Mat,filename,part,ind,x,pos,alpha)                      
type(cMatrix),intent(in)    :: mat(:)
CHARACTER(len=*),intent(in) :: filename
real(8),intent(in),optional :: alpha,ind,x(:)
CHARACTER(len=1),intent(in) :: part
CHARACTER(len=*),intent(in),optional :: pos
! ---
integer     :: ii,nx
real(8)     :: val,aa
if (present(pos)) then
    open(unit=101,file=trim(filename),position=trim(pos))
else
    open(unit=101,file=trim(filename))
endif
nx= size(Mat)
if (present(alpha)) then
    aa = alpha
else
    aa = 1.0d0
endif
do ii = 1,nx
    if (part .eq. 'R') then
        val = real(trace_c(Mat(ii)%m)) * aa
    else if (part .eq. 'A') then
        val = dble(trace_c(Mat(ii)%m)) * aa
    else
        val = aimag(trace_c(Mat(ii)%m)) * aa
    endif
    if (present(ind))then
        if (present(x)) then
            write(101,*) ind,x(ii),val
        else
            write(101,*) ind,ii,val
        end if
    else
        if (present(x)) then
            write(101,*) x(ii),val
        else
            write(101,*) ii,val
        end if
    endif
enddo
close(101)
end subroutine print_cMatrix

!=================================================================================
!                        add two lists of matrix                         
subroutine add_cMatrix(A,alpha,B,beta,C)                                         
type(cMatrix),intent(in)    :: B(:),A(:)
type(cMatrix),intent(inout) :: C(:)
real(8), intent(in)         :: alpha, beta
integer :: nx, xx
nx = size(A)
if (nx.ne.size(B)) then
    write(*,*) 'Error add_cMatrix, matrix of different sizes !'
    call abort()
endif
do xx = 1,nx
    C(xx)%m = alpha * A(xx)%m + beta * B(xx)%m
enddo
end subroutine add_cMatrix

subroutine addto_cMatrix(A,alpha,B,beta) 
implicit none
type(cMatrix),intent(in)    :: B(:)
type(cMatrix),intent(inout) :: A(:)
real(8), intent(in)         :: alpha, beta
!
integer :: nx, xx
nx = size(A)
if (nx.ne.size(B)) then
    write(*,*) 'Error add_cMatrix, matrix of different sizes !'
    call abort()
endif
do xx = 1,nx
    A(xx)%m = alpha * A(xx)%m + beta * B(xx)%m
enddo
end subroutine addto_cMatrix

!=================================================================================
!                       clone a list of matrix                                   
subroutine copy_cMatrix(fromMat,toMat)                                           
type(cMatrix),intent(in) :: fromMat(:)
type(cMatrix),intent(inout) :: toMat(:)
integer :: nx, xx
nx = size(fromMat)
do xx = 1,nx
    toMat(xx)%m = fromMat(xx)%m
enddo
end subroutine copy_cMatrix

!=================================================================================
!                       Sizes of each matrix in the list                         
Function size_cMatrix(Mat,nx) result(nm)                                         
type(cMatrix),intent(in) :: Mat(nx)
integer,intent(in)       :: nx
integer,dimension(2,nx)  :: nm
! ---
integer :: xx
do xx = 1,nx
    nm(1,xx) = size(Mat(xx)%m,1)
    nm(2,xx) = size(Mat(xx)%m,2)
enddo
return
end Function size_cMatrix

!=================================================================================
!               allocate list and each matrix by knowing the sizes               
subroutine alloc2_cMatrix(Mat,Nen,nm,template)                                        
type(cMatrix),intent(inout),allocatable :: Mat(:,:)
integer,intent(in),optional  :: nm(:,:)
integer,intent(in)           :: Nen
type(cMatrix),intent(in),optional :: template(:)
! ----
integer :: xx, nx,ee
if (present(nm)) then
nx = size(nm,2)
allocate(Mat(nx,Nen))
do ee = 1,Nen
do xx = 1,nx
allocate(Mat(xx,ee)%m(nm(1,xx),nm(2,xx)))
enddo
enddo
else if (present(template)) then
nx = size(template)
allocate(Mat(nx,Nen))
do ee = 1,Nen
do xx = 1,nx
allocate(Mat(xx,ee)%m(size(template(xx)%m,1),size(template(xx)%m,2)))
enddo
enddo
endif
end subroutine alloc2_cMatrix

!=================================================================================
!                   initiate each matrix to 0                                    
subroutine init2_cMatrix(Mat)                                                    
type(cMatrix),intent(inout) :: Mat(:,:)
integer :: xx, nx, ee, nen
nx = size(mat,1)
nen= size(mat,2)
do ee = 1,nen
do xx = 1,nx
Mat(xx,ee)%m = dcmplx(0.0d0, 0.0d0)
enddo
enddo
end subroutine init2_cMatrix

!=================================================================================
!                   free each matrix and the list                                
subroutine free2_cMatrix(mat)                                             
type(cMatrix),intent(inout),allocatable,dimension(:,:) :: mat
!
integer :: xx,nx,ee,nen
nx = size(mat,1)
nen= size(mat,2)
do ee = 1,nen
do xx = 1,nx
deallocate(mat(xx,ee)%m)
enddo
enddo
deallocate(mat)
end subroutine free2_cMatrix

!=================================================================================
!               allocate list and each matrix by knowing the sizes               
subroutine alloc_cMatrix(Mat,nm,template)                                        
type(cMatrix),intent(inout),allocatable :: Mat(:)
integer,intent(in),optional  :: nm(:,:)
type(cMatrix),intent(in),optional :: template(:)
!
integer :: xx, nx
if (present(nm)) then
    nx = size(nm,2)
    allocate(Mat(nx))
    do xx = 1,nx
        allocate(Mat(xx)%m(nm(1,xx),nm(2,xx)))
    enddo
else if (present(template)) then
    nx = size(template)
    allocate(Mat(nx))
    do xx = 1,nx
       allocate(Mat(xx)%m(size(template(xx)%m,1),size(template(xx)%m,2)))
   enddo
else
endif
end subroutine alloc_cMatrix

!=================================================================================
!                   initiate each matrix to 0                                    
subroutine init_cMatrix(Mat)                                                     
!=================================================================================
implicit none
type(cMatrix),intent(inout) :: Mat(:)
!
integer :: xx, nx
nx = size(mat)
do xx = 1,nx
    Mat(xx)%m = dcmplx(0.0d0, 0.0d0)
enddo
end subroutine init_cMatrix

!=================================================================================
!                   free each matrix and the list                                
subroutine free_cMatrix(Hii,H1i,Sii)                                             
type(cMatrix),intent(inout),allocatable,dimension(:) :: Hii
type(cMatrix),intent(inout),allocatable,dimension(:),optional :: H1i,Sii
!
integer :: xx,nx
nx = size(Hii)
do xx = 1,nx
deallocate(Hii(xx)%m)
if (present(H1i)) deallocate(H1i(xx)%m)
if (present(Sii)) deallocate(Sii(xx)%m)
enddo
deallocate(Hii)
if (present(H1i)) deallocate(H1i)
if (present(Sii)) deallocate(Sii)
end subroutine free_cMatrix


!*********************************! matrix operations !*********************************!

Function trace_c(A) result(tr)
complex(8),intent(in),dimension(:,:) :: A
complex(8) :: tr
integer :: i
if (size(A,1).ne.size(A,2))then
write(*,*) 'Error in Trace, matrix not square'
stop
end if
tr = dcmplx(0.0d0,0.0d0)
do i=1,size(A,1)
tr = tr + A(i,i)
end do
return
end Function trace_c

subroutine print_matrix_r(unit,A)
integer, intent(in) :: unit
real(8),intent(in),dimension(:,:) :: A
integer :: i,j
do i=1,size(A,1)
do j=1,size(A,2)
write(unit,*) i,j,A(i,j)
end do
end do
end subroutine print_matrix_r


subroutine show_matrix_c(A)
complex(8),intent(in),dimension(:,:) :: A
integer :: ii,jj
do ii = 1,size(A,1)
    do jj = 1,size(A,2)
!        write(*,'(A)',advance='no') ' ('//dbstring(dble(A(ii,jj)))//','//dbstring(aimag(A(ii,jj)))//') '
       write(*,'(A)') ' ('//dbstring(dble(A(ii,jj)))//','//dbstring(aimag(A(ii,jj)))//') '
    enddo
    write(*,*)
enddo
end subroutine


subroutine print_matrix_c(funit,A)
integer, intent(in) :: funit
complex(8),intent(in),dimension(:,:) :: A
integer :: i,j
do i=1,size(A,1)
    do j=1,size(A,2)
        write(funit,*) i,j,dble(A(i,j)),aimag(A(i,j))
    end do
end do
end subroutine print_matrix_c

subroutine invert_c(A)
implicit none
integer :: info,lda,lwork,n,nnz
integer,allocatable,dimension(:) :: ipiv
complex(8),intent(inout), dimension(:,:) :: A
complex(8),allocatable,dimension(:) :: work
if (size(A,1).ne.size(A,2))then
    write(*,*) 'Error in Invert! A is not square',size(A,1),size(A,2)
    call abort()
else
    n = size(A,1)
    allocate(ipiv(n))
    allocate(work(n*n))
    call zgetrf(n,n,A,n,ipiv,info)
    call zgetri(n,A,n,ipiv,work,n*n,info)
    deallocate(ipiv)
    deallocate(work)
end if
end subroutine invert_c

subroutine eye_r(Id,n)
real(8),intent(inout),dimension(:,:),Allocatable :: Id
integer,intent(in) :: n
integer :: i
if (.not.(allocated(Id)))then
    allocate(Id(n,n))
else if((size(Id,1).ne.n).or.(size(Id,2).ne.n))then
    deallocate(Id)
    allocate(Id(n,n))
end if
Id=dble(0.0d0)
do i=1,n
    Id(i,i) = dble(1.0d0)
end do
end subroutine eye_r


subroutine eye_c(Id,n)
complex(8),intent(inout),dimension(:,:),Allocatable :: Id
integer,intent(in) :: n
integer :: i
if (.not.(allocated(Id)))then
    allocate(Id(n,n))
else if((size(Id,1).ne.n).or.(size(Id,2).ne.n))then
    deallocate(Id)
    allocate(Id(n,n))
end if
Id=dcmplx(0.0d0,0.0d0)
do i=1,n
    Id(i,i) = dcmplx(1.0d0,0.0d0)
end do
end subroutine eye_c

subroutine tridiagBloc_C(H00,Hu, H)
COMPLEX(8),intent(in),dimension(:,:,:) :: H00,Hu
COMPLEX(8),intent(inout),allocatable:: H(:,:)
integer :: m,n,l
integer :: i
l = size(H00,3)
if (size(Hu,3).lt.(l-1))then
    write(*,*)"ERROR in tridiagBloc! Hu(:,:,l) l=",size(Hu,3),"H00(:,:,l) l=",size(H00,3)
    call abort()
endif
m = size(H00,1)
n = size(H00,2)
if (size(Hu,1).ne.m)then
    write(*,*)"ERROR in tridiagBloc! Hu(m,:,:) m=",size(Hu,1),"H00(m,:,:) m=",m
    call abort()
endif
if (size(Hu,2).ne.n)then
    write(*,*)"ERROR in tridiagBloc! Hu(:,n,:) n=",size(Hu,2),"H00(:,n,:) n=",n
    call abort()
endif
Allocate(H(m*l,n*l))
H = dcmplx(0.0d0,0.0d0)
do i=1,l
    H((i-1)*m+1:i*m,(i-1)*n+1:i*n) = H00(:,:,i)
    if (i>1)then
        H((i-2)*m+1:(i-1)*m,(i-1)*n+1:i*n) = Hu(:,:,i-1)
        H((i-1)*m+1:i*m,(i-2)*n+1:(i-1)*n) = conjg(TRANSPOSE(Hu(:,:,i-1)))
    end if
enddo
end subroutine tridiagBloc_C

subroutine tridiagBloc_R(H00,Hu, H)
real(8),intent(in),dimension(:,:,:) :: H00,Hu
real(8),intent(inout),Allocatable :: H(:,:)
integer :: m,n,l
integer :: i
l = size(H00,3)
if (size(Hu,3).lt.(l-1))then
    write(*,*)"ERROR in tridiagBloc! Hu(:,:,l) l=",size(Hu,3),"H00(:,:,l) l=",size(H00,3)
    call abort()
endif
m = size(H00,1)
n = size(H00,2)
if (size(Hu,1).ne.m)then
    write(*,*)"ERROR in tridiagBloc! Hu(m,:,:) m=",size(Hu,1),"H00(m,:,:) m=",m
    call abort()
endif
if (size(Hu,2).ne.n)then
    write(*,*)"ERROR in tridiagBloc! Hu(:,n,:) n=",size(Hu,2),"H00(:,n,:) n=",n
    call abort()
endif
Allocate(H(m*l,n*l))
H = 0.0d0
do i=1,l
    H((i-1)*m+1:i*m,(i-1)*n+1:i*n) = H00(:,:,i)
    if (i>1)then
        H((i-2)*m+1:(i-1)*m,(i-1)*n+1:i*n) = Hu(:,:,i-1)
        H((i-1)*m+1:i*m,(i-2)*n+1:(i-1)*n) = TRANSPOSE(Hu(:,:,i-1))
    end if
enddo
end subroutine tridiagBloc_R

subroutine triMUL_C(A,B,C,R,trA,trB,trC)
complex(8),intent(in),dimension(:,:) :: A,B,C
complex(8),intent(inout),allocatable :: R(:,:)
character,intent(in) :: trA,trB,trC
complex(8),allocatable,dimension(:,:) :: tmp
integer :: n,m,k,kb
if((trA.ne.'n').and.(trA.ne.'N').and.(trA.ne.'t').and.(trA.ne.'T').and.(trA.ne.'c').and.(trA.ne.'C'))then
write(*,*) "ERROR in triMUL_C! trA is wrong: ",trA
call abort()
endif
if((trB.ne.'n').and.(trB.ne.'N').and.(trB.ne.'t').and.(trB.ne.'T').and.(trB.ne.'c').and.(trB.ne.'C'))then
write(*,*) "ERROR in triMUL_C! trB is wrong: ",trB
call abort()
endif
if ((trA.eq.'n').or.(trA.eq.'N'))then
    k = size(A,2)
    m = size(A,1)
else
    k = size(A,1)
    m = size(A,2)
endif
if ((trB.eq.'n').or.(trB.eq.'N'))then
    kb = size(B,1)
    n = size(B,2)
else
    kb = size(B,2)
    n = size(B,1)
endif
if(k.ne.kb)then
    write(*,*) "ERROR in triMUL_C! Matrix dimension is wrong",k,kb
    call abort()
end if
call MUL_C(A,B,trA,trB,tmp)
call MUL_C(tmp,C,'n',trC,R)
deallocate(tmp)
end subroutine triMUL_C



subroutine triMUL_R(A,B,C,R,trA,trB,trC)
real(8),intent(in),dimension(:,:) :: A,B,C
real(8),intent(inout),allocatable :: R(:,:)
character,intent(in) :: trA,trB,trC
real(8),allocatable,dimension(:,:) :: tmp
integer :: n,m,k,kb
if((trA.ne.'n').and.(trA.ne.'N').and.(trA.ne.'t').and.(trA.ne.'T').and.(trA.ne.'c').and.(trA.ne.'C'))then
    write(*,*) "ERROR in triMUL_C! trA is wrong: ",trA
    call abort()
endif
if((trB.ne.'n').and.(trB.ne.'N').and.(trB.ne.'t').and.(trB.ne.'T').and.(trB.ne.'c').and.(trB.ne.'C'))then
    write(*,*) "ERROR in triMUL_C! trB is wrong: ",trB
    call abort()
endif
if ((trA.eq.'n').or.(trA.eq.'N'))then
    k = size(A,2)
    m = size(A,1)
else
    k = size(A,1)
    m = size(A,2)
endif
if ((trB.eq.'n').or.(trB.eq.'N'))then
    kb = size(B,1)
    n = size(B,2)
else
    kb = size(B,2)
    n = size(B,1)
endif
if(k.ne.kb)then
write(*,*) "ERROR in triMUL_C! Matrix dimension is wrong",k,kb
call abort()
end if
call MUL_R(A,B,trA,trB,tmp)
call MUL_R(tmp,C,'n',trC,R)
deallocate(tmp)
end subroutine triMUL_R



subroutine MUL_C(A,B,trA,trB,R)
complex(8),intent(in) :: A(:,:), B(:,:)
complex(8),intent(inout),allocatable :: R(:,:)
CHARACTER, intent(in) :: trA, trB
integer :: n,m,k,kb, lda,ldb
if((trA.ne.'n').and.(trA.ne.'N').and.(trA.ne.'t').and.(trA.ne.'T').and.(trA.ne.'c').and.(trA.ne.'C'))then
write(*,*) "ERROR in MUL_C! trA is wrong: ",trA
call abort()
endif
if((trB.ne.'n').and.(trB.ne.'N').and.(trB.ne.'t').and.(trB.ne.'T').and.(trB.ne.'c').and.(trB.ne.'C'))then
write(*,*) "ERROR in MUL_C! trB is wrong: ",trB
call abort()
endif
lda = size(A,1)
ldb = size(B,1)
if ((trA.eq.'n').or.(trA.eq.'N'))then
k = size(A,2)
m = size(A,1)
else
k = size(A,1)
m = size(A,2)
endif
if ((trB.eq.'n').or.(trB.eq.'N'))then
kb = size(B,1)
n = size(B,2)
else
kb = size(B,2)
n = size(B,1)
endif
if(k.ne.kb)then
write(*,*) "ERROR in MUL_C! Matrix dimension is wrong",k,kb
call abort()
end if
if (allocated(R)) then
    if ((size(R,1).ne.m).or.(size(R,2).ne.n)) then
        deallocate(R)
        Allocate(R(m,n))
    end if
else
    Allocate(R(m,n))
end if
R = dcmplx(0.0d0,0.0d0)
call zgemm(trA,trB,m,n,k,dcmplx(1.0d0,0.0d0),A,lda,B,ldb,dcmplx(0.0d0,0.0d0),R,m)
end subroutine MUL_C


subroutine MUL_R(A,B,trA,trB,R)
real(8),intent(in) :: A(:,:), B(:,:)
real(8),intent(inout),allocatable :: R(:,:)
CHARACTER, intent(in) :: trA, trB
integer :: n,m,k,kb, lda,ldb
if((trA.ne.'n').and.(trA.ne.'N').and.(trA.ne.'t').and.(trA.ne.'T').and.(trA.ne.'c').and.(trA.ne.'C'))then
write(*,*) "ERROR in MUL_C! trA is wrong: ",trA
call abort()
endif
if((trB.ne.'n').and.(trB.ne.'N').and.(trB.ne.'t').and.(trB.ne.'T').and.(trB.ne.'c').and.(trB.ne.'C'))then
write(*,*) "ERROR in MUL_C! trB is wrong: ",trB
call abort()
endif
lda = size(A,1)
ldb = size(B,1)
if ((trA.eq.'n').or.(trA.eq.'N'))then
k = size(A,2)
m = size(A,1)
else
k = size(A,1)
m = size(A,2)
endif
if ((trB.eq.'n').or.(trB.eq.'N'))then
kb = size(B,1)
n = size(B,2)
else
kb = size(B,2)
n = size(B,1)
endif
if(k.ne.kb)then
    write(*,*) "ERROR in MUL_C! Matrix dimension is wrong",k,kb
    call abort()
end if
if (allocated(R)) then
    if ((size(R,1).ne.m).or.(size(R,2).ne.n)) then
        deallocate(R)
        Allocate(R(m,n))
    end if
else
    Allocate(R(m,n))
end if
R = 0.0d0
call dgemm(trA,trB,m,n,k,1.0d0,A,lda,B,ldb,0.0d0,R,m)
end subroutine MUL_R

end module matrix_mod
