!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 08/2023

module kpHam

    implicit none 

    private
    
    public :: read_KP_coeff, save_KP_coeff
    public :: save_Ham_blocks
    public :: calc_bands_at, calc_KP_block, calc_kpbands_at
    public :: generate_LK4,generate_LK6
    public :: test_bandstructure
    

    complex(8), parameter :: cone = dcmplx(1.0d0,0.0d0)
    complex(8), parameter :: czero  = dcmplx(0.0d0,0.0d0)
    complex(8), parameter :: c1i  = dcmplx(0.0d0,1.0d0)
    REAL(8), PARAMETER :: pi = 3.14159265359d0
    REAL(8), PARAMETER :: twopi = 3.14159265359d0*2.0d0
    real(8), parameter :: hbar=1.0546d-34 ! m^2 kg / s
    real(8), parameter :: m0=9.109d-31 ! kg
    real(8), parameter :: eps0=8.854d-12 ! C/V/m 
    real(8), parameter :: c0=2.998d8 ! m/s
    real(8), parameter :: e0=1.6022d-19 ! C
    real(8), parameter :: hb2m=hbar**2/2.0d0/m0/e0  ! hbar^2/(2m0) eV*m^2
        
    integer, parameter :: num_k = 10 
    ! index for different types of k-terms available in the KP model
    integer, parameter :: const = 1
    integer, parameter :: dx    = 2
    integer, parameter :: dy    = 3
    integer, parameter :: dz    = 4
    integer, parameter :: dxdy  = 5
    integer, parameter :: dxdz  = 6
    integer, parameter :: dydz  = 7
    integer, parameter :: dx2   = 8
    integer, parameter :: dy2   = 9
    integer, parameter :: dz2   = 10

CONTAINS

    ! generate Luttinger 6-band KP model coefficients
    !   PRB 98, 155319 (2018) Eq (D1) (D4)
    subroutine generate_LK6(KPcoeff, gam1,gam2,gam3,delta,kap)        
        real(8), intent(in) :: gam1,gam2,gam3,kap,delta
        complex(8), intent(out) :: KPcoeff(6,6,num_k) ! KP coeff. table        
        ! ----
        complex(8),dimension(num_k) :: P,Q,S,R,D
        integer::i,m,n
        real(8),parameter :: sq2 = sqrt(2.0d0)
        real(8),parameter :: sq3o2 = sqrt(3.0d0/2.0d0)
        ! construct the coeff. table in another (better) way
        P=czero
        Q=czero
        S=czero
        R=czero
        D=czero
        !
        P(dx2) = hb2m * gam1
        P(dy2) = hb2m * gam1
        P(dz2) = hb2m * gam1
        !
        Q(dx2) = hb2m * gam2
        Q(dy2) = hb2m * gam2
        Q(dz2) = -2.0d0 * hb2m * gam2
        !
        R(dx2) = - hb2m * sqrt(3.0d0) * gam2
        R(dy2) = + hb2m * sqrt(3.0d0) * gam2
        R(dxdy) = hb2m * sqrt(3.0d0) * 2.0d0 * c1i * gam3
        !
        S(dxdz) = hb2m * 2.0d0 * sqrt(3.0d0) * gam3 
        S(dydz) = hb2m * 2.0d0 * sqrt(3.0d0) * gam3 * (-c1i) 
        !
        D(const) = delta
        !
        KPcoeff=czero
        !
        do i = 1,num_k
          KPcoeff(:,:,i) = reshape( &
                (/ P(i)+Q(i),      -S(i),      R(i) ,    czero ,         S(i)/sq2 ,        -sq2*R(i) ,&
                     czero  ,  P(i)-Q(i),     czero ,     R(i) ,         sq2*Q(i) ,      -sq3o2*S(i) ,&
                     czero  ,    czero  , P(i)-Q(i) ,     S(i) ,-sq3o2*conjg(S(i)),      -sq2*Q(i)   ,&
                     czero  ,    czero  ,   czero   , P(i)+Q(i), sq2*conjg( R(i) ), conjg( S(i) )/sq2,&
                     czero  ,    czero  ,   czero   ,   czero  ,         P(i)+D(i),         czero    ,&
                     czero  ,    czero  ,   czero   ,   czero  ,          czero   ,      P(i)+D(i) /),&
                (/6, 6/) )
          do n=1,6
            do m=1,n-1
                KPcoeff(m,n,i) = conjg( KPcoeff(n,m,i) )               
            enddo
          enddo
        enddo                    
    end subroutine generate_LK6



    ! generate Luttinger 4-band KP model coefficients
    !   PRB 98, 155319 (2018) Eq (D1) (D4)
    subroutine generate_LK4(KPcoeff, gam1,gam2,gam3,kap)        
        real(8), intent(in) :: gam1,gam2,gam3,kap
        complex(8), intent(out) :: KPcoeff(4,4,num_k) ! KP coeff. table        
        ! ----
        complex(8)::LK6(6,6,num_k)
        !
        KPcoeff=czero
        call generate_LK6(LK6,gam1,gam2,gam3,0.0d0,kap)
        KPcoeff(:,:,:)=LK6(1:4,1:4,:)
    end subroutine generate_LK4


    ! read the KP coefficients from a file
    subroutine read_KP_coeff(filename, nbnd, KPcoeff)
        character(len=*),intent(in)::filename ! file name
        integer,intent(in)::nbnd ! number of bands
        complex(8), intent(out) :: KPcoeff(nbnd,nbnd,num_k) ! KP coeff. table
        ! ----
        integer::i,j,k
        real(8)::re,im
        open(unit=11, file=filename, status='unknown')
        do i=1,nbnd
          do j=1,nbnd
            do k=1,num_k
              read(11,*) re,im
              KPcoeff(i,j,k) = dcmplx(re,im)
            enddo
          enddo
        enddo
        close(11)
    end subroutine read_KP_coeff


    ! save the KP coefficients to a file
    subroutine save_KP_coeff(filename, nbnd, KPcoeff)
        character(len=*),intent(in)::filename ! file name
        integer,intent(in)::nbnd ! number of bands
        complex(8), intent(in) :: KPcoeff(nbnd,nbnd,num_k) ! KP coeff. table
        ! ----
        integer::i,j,k        
        open(unit=11, file=filename, status='unknown')
        do i=1,nbnd
          do j=1,nbnd
            do k=1,num_k
              write(11,'(2E18.6)') dble( KPcoeff(i,j,k) ), aimag( KPcoeff(i,j,k) )              
            enddo
          enddo
        enddo
        close(11)
    end subroutine save_KP_coeff

    ! save the KP Hamiltonian blocks to a file
    subroutine save_Ham_blocks(filename, nbnd, KPcoeff,dd)
        character(len=*),intent(in)::filename ! file name
        integer,intent(in)::nbnd ! number of bands
        complex(8), intent(in) :: KPcoeff(nbnd,nbnd,num_k) ! KP coeff. table
        real(8), intent(in) :: dd(3) ! discretization step size in x-y-z
        ! ----
        integer::i(3),j(3),ix,iy,iz,m       
        complex(8)::Hij(nbnd,nbnd)
        open(unit=11, file=filename, status='unknown')                
        i = (/0,0,0/)
        !
        do ix = -1,1
          do iy = -1,1
            do iz = -1,1
              j = (/ix,iy,iz/)
              write(11,'(3I10)') j
              call calc_KP_block(i,j,dd,nbnd,KPcoeff,Hij)
              do m=1,nbnd
                write(11,'(100E18.6)') Hij(:,m)
              enddo
              write(11,*)
            enddo
          enddo
        enddo
        close(11)
    end subroutine save_Ham_blocks


    ! compute the bands at a k-point before discretization
    subroutine calc_kpbands_at(nbnd,KPcoeff,kpt,Ek,veck)
    integer, intent(in) :: nbnd ! number of bands        
        complex(8), intent(in) :: KPcoeff(nbnd,nbnd,num_k) ! KP coeff. table
        real(8), intent(in) :: kpt(3) ! k-point
        real(8), intent(out) :: Ek(nbnd) ! energies
        complex(8), intent(out), optional :: veck(nbnd,nbnd) ! eigen-vectors
        ! ---
        real(8)::kx,ky,kz                
        complex(8)::Ham(nbnd,nbnd)
        kx=kpt(1)
        ky=kpt(2)
        kz=kpt(3)
        Ham=czero        
        Ham=Ham + KPcoeff(:,:,const)
        Ham=Ham + KPcoeff(:,:,dx)*kx
        Ham=Ham + KPcoeff(:,:,dy)*ky
        Ham=Ham + KPcoeff(:,:,dz)*kz
        Ham=Ham + KPcoeff(:,:,dxdz)*kx*kz
        Ham=Ham + KPcoeff(:,:,dydz)*ky*kz
        Ham=Ham + KPcoeff(:,:,dxdy)*kx*ky
        Ham=Ham + KPcoeff(:,:,dx2)*kx*kx
        Ham=Ham + KPcoeff(:,:,dy2)*ky*ky
        Ham=Ham + KPcoeff(:,:,dz2)*kz*kz
        !
        if (present(veck)) then
            Ek = eigv(nbnd,Ham)
            veck = Ham
        else
            Ek = eig(nbnd,Ham)
        endif               
    end subroutine calc_kpbands_at


    ! compute the bands at a k-point using the Ham blocks
    subroutine calc_bands_at(nbnd,KPcoeff,dd,kpt,Ek,veck)
        integer, intent(in) :: nbnd ! number of bands
        real(8), intent(in) :: dd(3) ! discretization step size in x-y-z
        complex(8), intent(in) :: KPcoeff(nbnd,nbnd,num_k) ! KP coeff. table
        real(8), intent(in) :: kpt(3) ! k-point
        real(8), intent(out) :: Ek(nbnd) ! energies
        complex(8), intent(out), optional :: veck(nbnd,nbnd) ! eigen-vectors
        ! ---        
        integer::i(3),j(3),ix,iy,iz
        complex(8)::Ham(nbnd,nbnd),Hij(nbnd,nbnd)
        Ham=czero        
        !
        i = (/0,0,0/)
        !
        do ix = -1,1
          do iy = -1,1
            do iz = -1,1
              j = (/ix,iy,iz/)
              call calc_KP_block(i,j,dd,nbnd,KPcoeff,Hij)
              Ham = Ham + Hij * exp(-c1i * dot_product( kpt, dd*dble(j) ) )
            enddo
          enddo
        enddo
        !
        if (present(veck)) then
            Ek = eigv(nbnd,Ham)
            veck = Ham
        else
            Ek = eig(nbnd,Ham)
        endif        
    end subroutine calc_bands_at


    ! compute the KP Hamiltonian block between point i and j, based on 
    !   the KP coeff. table
    subroutine calc_KP_block(i,j,dd,nbnd,KPcoeff,Hblock)
        integer, intent(in) :: i(3),j(3)
        real(8), intent(in) :: dd(3) ! discretization step size in x-y-z
        integer, intent(in) :: nbnd ! number of bands
        complex(8), intent(in) :: KPcoeff(nbnd,nbnd,num_k) ! KP coeff. table
        complex(8), intent(out) :: Hblock(nbnd,nbnd) ! Hamiltonian block
        ! -----
        integer::ik
        complex(8)::FDcoeff(num_k)
        call calc_FD_coeff(i,j,dd,FDcoeff)        
        Hblock = czero
        do ik=1,num_k        
          Hblock = Hblock + KPcoeff(:,:,ik) * FDcoeff(ik) 
        enddo
    end subroutine calc_KP_block    


    ! compute the finite-difference coefficients between point i and j
    !   for k-terms
    subroutine calc_FD_coeff(i,j,dd,coeff)
        integer, intent(in) :: i(3), j(3) ! point x-y-z index 
        real(8), intent(in) :: dd(3) ! discretization step size in x-y-z
        complex(8), intent(out):: coeff(num_k) ! coefficients for k-terms
        ! -----
        integer::ij(3)
        real(8)::ddx,ddy,ddz
        ! k_v -> -i d/d_v 
        ! use the following rules
        !   d/dx = (psi_i+1 - psi_i-1)/2/dx
        !   d^2/dx^2 = (psi_i+1 + psi_i-1 - 2psi_i)/dx^2
        !   d^2/dx/dy = (psi_i+1,j+1 + psi_i-1,j-1 - psi_i+1,j-1 - psi_i-1,j+1) / (4 dx dy)        
        ij = i-j        
        ddx=dd(1)
        ddy=dd(2)
        ddz=dd(3)
        coeff(:) = czero
        ! onsite 
        if (all(ij == (/0,0,0/))) then
          coeff(const) = 1.0d0
          coeff(dx2) = -2.0d0/ddx/ddx * (-c1i)**2
          coeff(dy2) = -2.0d0/ddy/ddy * (-c1i)**2
          coeff(dz2) = -2.0d0/ddz/ddz * (-c1i)**2
        endif
        ! +x
        if (all(ij == (/1,0,0/))) then
          coeff(dx) = 1.0d0/2.0d0/ddx * (-c1i)
          coeff(dx2) = 1.0d0/ddx**2   * (-c1i)**2          
        endif
        ! -x
        if (all(ij == (/-1,0,0/))) then
          coeff(dx) = -1.0d0/2.0d0/ddx * (-c1i)
          coeff(dx2) = 1.0d0/ddx**2    * (-c1i)**2
        endif
        ! +y
        if (all(ij == (/0,1,0/))) then
          coeff(dy) = 1.0d0/2.0d0/ddy  * (-c1i)
          coeff(dy2) = 1.0d0/ddy**2    * (-c1i)**2
        endif
        ! -y
        if (all(ij == (/0,-1,0/))) then
          coeff(dy) = -1.0d0/2.0d0/ddy  * (-c1i)
          coeff(dy2) = 1.0d0/ddy**2     * (-c1i)**2
        endif
        ! +z
        if (all(ij == (/0,0,1/))) then
          coeff(dz) = 1.0d0/2.0d0/ddz * (-c1i)
          coeff(dz2) = 1.0d0/ddz**2   * (-c1i)**2
        endif
        ! -z
        if (all(ij == (/0,0,-1/))) then
          coeff(dz) = -1.0d0/2.0d0/ddz  * (-c1i)
          coeff(dz2) = 1.0d0/ddz**2     * (-c1i)**2
        endif
        ! +x +y
        if (all(ij == (/1,1,0/))) then
          coeff(dxdy) = 1.0d0/4.0d0/ddx/ddy * (-c1i)**2
        endif
        ! +x -y
        if (all(ij == (/1,-1,0/))) then
          coeff(dxdy) = -1.0d0/4.0d0/ddx/ddy * (-c1i)**2
        endif
        ! -x +y
        if (all(ij == (/-1,1,0/))) then
          coeff(dxdy) = -1.0d0/4.0d0/ddx/ddy * (-c1i)**2
        endif
        ! -x -y
        if (all(ij == (/-1,-1,0/))) then
          coeff(dxdy) = 1.0d0/4.0d0/ddx/ddy * (-c1i)**2
        endif        
        ! +x +z
        if (all(ij == (/1,0,1/))) then
          coeff(dxdz) = 1.0d0/4.0d0/ddx/ddz * (-c1i)**2
        endif
        ! +x -z
        if (all(ij == (/1,0,-1/))) then
          coeff(dxdz) = -1.0d0/4.0d0/ddx/ddz * (-c1i)**2
        endif 
        ! -x -z
        if (all(ij == (/-1,0,-1/))) then
          coeff(dxdz) = 1.0d0/4.0d0/ddx/ddz * (-c1i)**2
        endif 
        ! -x +z
        if (all(ij == (/-1,0,1/))) then
          coeff(dxdz) = -1.0d0/4.0d0/ddx/ddz * (-c1i)**2
        endif 
        ! +y +z
        if (all(ij == (/0,1,1/))) then
          coeff(dydz) = 1.0d0/4.0d0/ddy/ddz * (-c1i)**2
        endif 
        ! -y +z
        if (all(ij == (/0,-1,1/))) then
          coeff(dydz) = -1.0d0/4.0d0/ddy/ddz * (-c1i)**2
        endif 
        ! +y -z
        if (all(ij == (/0,1,-1/))) then
          coeff(dydz) = -1.0d0/4.0d0/ddy/ddz * (-c1i)**2
        endif 
        ! -y -z
        if (all(ij == (/0,-1,-1/))) then
          coeff(dydz) = 1.0d0/4.0d0/ddy/ddz * (-c1i)**2
        endif         
    end subroutine calc_FD_coeff


    ! vector cross-product
    FUNCTION cross(a, b)    
        REAL(8), DIMENSION(3) :: cross
        REAL(8), DIMENSION(3), INTENT(IN) :: a, b
        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)
    END FUNCTION cross
    
    
    ! calculate eigen-values of a Hermitian matrix A
    FUNCTION eig(NN, A)    
        INTEGER, INTENT(IN) :: NN
        COMPLEX(8), INTENT(INOUT), DIMENSION(:,:) :: A
        ! -----
        REAL(8) :: eig(NN)
        real(8) :: W(1:NN)
        integer :: INFO,LWORK,liwork, lrwork
        complex(8), allocatable :: work(:)
        real(8), allocatable :: RWORK(:)
        !integer, allocatable :: iwork(:) 
        lwork= max(1,2*NN-1)
        lrwork= max(1,3*NN-2)
        allocate(work(lwork))
        allocate(rwork(lrwork))
        !
        CALL zheev( 'N','U', NN, A, NN, W, WORK, LWORK, RWORK, INFO )
        !
        deallocate(work,rwork)
        if (INFO.ne.0)then
        write(*,*)'SEVERE WARNING: ZHEEV HAS FAILED. INFO=',INFO
        call abort
        endif
        eig(:)=W(:)
    END FUNCTION eig
    
    
    ! calculate eigen-values and eigen-vectors of a Hermitian matrix A
    !   upon return A will be modified and contains the eigen-vectors
    FUNCTION eigv(NN, A)    
        INTEGER, INTENT(IN) :: NN
        COMPLEX(8), INTENT(INOUT), DIMENSION(:,:) :: A
        ! -----
        REAL(8) :: eigv(NN)
        real(8) :: W(1:NN)
        integer :: INFO,LWORK,liwork, lrwork
        complex(8), allocatable :: work(:)
        real(8), allocatable :: RWORK(:)
        !integer, allocatable :: iwork(:) 
        lwork= max(1,2*NN-1)
        lrwork= max(1,3*NN-2)
        allocate(work(lwork))
        allocate(rwork(lrwork))
        !
        CALL zheev( 'V','U', NN, A, NN, W, WORK, LWORK, RWORK, INFO )
        !
        deallocate(work,rwork)
        if (INFO.ne.0)then
        write(*,*)'SEVERE WARNING: ZHEEV HAS FAILED. INFO=',INFO
        call abort
        endif
        eigv(:)=W(:)
    END FUNCTION eigv

!!!!!!!!!!!!!!!!!!!!!!!!   test  functions   !!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine test_bandstructure(nbnd,KPcoeff,dd)
        implicit none
        integer, intent(in)::nbnd
        real(8), intent(in) :: dd(3) ! discretization step size in x-y-z
        complex(8),intent(in)::KPcoeff(nbnd,nbnd,num_k)
        ! ----
        real(8)::kpt(3),kpt1(3),kpt2(3)
        real(8),allocatable::Ek(:)
        integer :: nk,ik            
        nk = 1000
        allocate(Ek(nbnd))
        !
        open(unit=10,file='ek0.dat',status='unknown')
        open(unit=11,file='ek.dat',status='unknown')
        !
        kpt1 = (/ -1.0d9 , 0.0d0 , 0.0d0 /)
        kpt2 = (/ +1.0d9 , 0.0d0 , 0.0d0 /)
        do ik = 1,nk
          kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1          
          call calc_bands_at(nbnd,KPcoeff,dd,kpt,Ek)          
          write(11,'(7E18.6)') dble(ik)/dble(nk), Ek          
          call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)          
          write(10,'(7E18.6)') dble(ik)/dble(nk), Ek          
        enddo
        !
        kpt1 = (/ 0.0d9 , -1.0d9 , 0.0d0 /)
        kpt2 = (/ 0.0d9 ,  1.0d9 , 0.0d0 /)
        do ik = 1,nk
          kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1          
          call calc_bands_at(nbnd,KPcoeff,dd,kpt,Ek)          
          write(11,'(7E18.6)') dble(ik)/dble(nk)+1.0d0, Ek          
          call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)          
          write(10,'(7E18.6)') dble(ik)/dble(nk)+1.0d0, Ek          
        enddo
        !
        kpt1 = (/ 0.0d9 , 0.0d9 , -1.0d9 /)
        kpt2 = (/ 0.0d9 , 0.0d9 ,  1.0d9 /)
        do ik = 1,nk
          kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1                    
          call calc_bands_at(nbnd ,KPcoeff,dd,kpt,Ek)          
          write(11,'(7E18.6)') dble(ik)/dble(nk)+2.0d0, Ek          
          call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)          
          write(10,'(7E18.6)') dble(ik)/dble(nk)+2.0d0, Ek          
        enddo
        !
        kpt1 = (/ -1.0d9 , 0.0d9 , -1.0d9 /)
        kpt2 = (/  1.0d9 , 0.0d9 ,  1.0d9 /)
        do ik = 1,nk
          kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1    
          call calc_bands_at(nbnd ,KPcoeff,dd,kpt,Ek)          
          write(11,'(7E18.6)') dble(ik)/dble(nk)+3.0d0, Ek          
          call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)          
          write(10,'(7E18.6)') dble(ik)/dble(nk)+3.0d0, Ek          
        enddo
        !
        kpt1 = (/ -1.0d9 , -1.0d9 , 0.0d9 /)
        kpt2 = (/  1.0d9 ,  1.0d9 , 0.0d9 /)
        do ik = 1,nk
          kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1          
          call calc_bands_at(nbnd ,KPcoeff,dd,kpt,Ek)          
          write(11,'(7E18.6)') dble(ik)/dble(nk)+4.0d0, Ek          
          call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)          
          write(10,'(7E18.6)') dble(ik)/dble(nk)+4.0d0, Ek          
        enddo
        !
        kpt1 = (/ -1.0d9 , -1.0d9 , -1.0d9 /)
        kpt2 = (/  1.0d9 ,  1.0d9 ,  1.0d9 /)
        do ik = 1,nk
          kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1          
          call calc_bands_at(nbnd ,KPcoeff,dd,kpt,Ek)          
          write(11,'(7E18.6)') dble(ik)/dble(nk)+5.0d0, Ek          
          call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)          
          write(10,'(7E18.6)') dble(ik)/dble(nk)+5.0d0, Ek          
        enddo
        !
        close(10)
        close(11)
        deallocate(Ek)
    end subroutine test_bandstructure
    
    
end module kpHam
