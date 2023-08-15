!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 08/2023

module kpHam

    implicit none 

    private
    
    public :: read_KP_coeff, save_KP_coeff
    public :: calc_bands_at, calc_KP_block, calc_kpbands_at
    public :: generate_LK4,generate_LK6

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
    real(8), parameter :: hb2m=7.62058d-20  ! hbar^2/(2m0) eV*m^2
        
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
        complex(8)::P,Q,S,R1,R2,R1p,R2p,Sp,Qp,Sp2
        !
        P = hb2m * gam1
        Q = hb2m * gam2
        R1 = hb2m * sqrt(3.0d0) * gam3
        R2 = hb2m * sqrt(3.0d0) * 2.0d0 * c1i * gam2
        S  = hb2m * 2.0d0 * sqrt(3.0d0) * gam3 
        R1p=sqrt(2.0d0)*R1
        R2p=sqrt(2.0d0)*R2
        Sp=S/sqrt(2.0d0)
        Sp2=-S*sqrt(3.0d0/2.0d0)
        Qp=sqrt(2.0d0)*Q
        !
        KPcoeff=czero
        !
        KPcoeff(:,:,const) = czero
        !
        KPcoeff(:,:,dx) = czero
        KPcoeff(:,:,dy) = czero
        KPcoeff(:,:,dz) = czero
        !
        KPcoeff(:,:,dx2) = reshape( &
                                  (/ P+Q   ,    czero ,  -R1 ,  czero , czero , +R1p  ,&
                                     czero ,    P-Q   , czero,  -R1   , Qp    , czero ,&
                                -conjg(R1) ,    czero ,  P-Q ,  czero ,  &
                                     czero ,-conjg(R1), czero,  P+Q     /),  (/6, 6/) )
        !
        KPcoeff(:,:,dy2) = reshape( &
                                  (/ P+Q   ,    czero ,   R1 ,  czero , czero , -R1p  ,&
                                     czero ,    P-Q   , czero,   R1   , Qp    , czero ,&
                                 conjg(R1) ,    czero ,  P-Q ,  czero , &
                                     czero , conjg(R1), czero,  P+Q     /),  (/6, 6/) )
                                             
        KPcoeff(:,:,dz2) = reshape( &
                                  (/ P-2.0*Q ,    czero , czero   ,  czero , czero , czero, &
                                     czero   ,  P+2.0*Q , czero   ,  czero ,-2.0*Qp, czero, &
                                     czero   ,    czero , P+2.0*Q ,  czero ,  &
                                     czero   ,    czero , czero   ,  P-2.0*Q  /),  (/6, 6/) )                                     
        KPcoeff(:,:,dxdy) = reshape( &
                                  (/ czero ,    czero ,  R2 ,  czero , czero , -R2p , &
                                     czero ,    czero , czero,  R2   , czero , czero, &
                                 conjg(R2) ,    czero , czero, czero ,  &
                                     czero , conjg(R2), czero, czero    /),  (/6, 6/) )
                                             
        KPcoeff(:,:,dxdz) = reshape( &
                                  (/ czero ,    -S ,  czero  ,  czero ,   Sp , czero , &
                                -conjg(S)  , czero ,  czero  ,  czero , czero, Sp2   , &
                                     czero , czero ,  czero  ,   S    , &
                                     czero , czero , conjg(S),  czero   /),  (/6, 6/) )
                                                                                          
        KPcoeff(:,:,dydz) = reshape( &
                                  (/ czero , c1i*S ,  czero     ,  czero , -c1i*Sp, czero    ,&
                              conjg(c1i*S) , czero ,  czero     ,  czero , czero  , -c1i*Sp2 ,&
                                     czero , czero ,  czero     ,  c1i*S , &
                                     czero , czero ,conjg(c1i*S),  czero   /),  (/6, 6/) )
                    
    end subroutine generate_LK6



    ! generate Luttinger 4-band KP model coefficients
    !   PRB 98, 155319 (2018) Eq (D1) (D4)
    subroutine generate_LK4(KPcoeff, gam1,gam2,gam3,kap)        
        real(8), intent(in) :: gam1,gam2,gam3,kap
        complex(8), intent(out) :: KPcoeff(4,4,num_k) ! KP coeff. table        
        ! ----
        complex(8)::P,Q,S,R1,R2
        !
        P = hb2m * gam1
        Q = hb2m * gam2
        R1 = hb2m * sqrt(3.0d0) * gam3
        R2 = hb2m * sqrt(3.0d0) * 2.0d0 * c1i * gam2
        S  = hb2m * 2.0d0 * sqrt(3.0d0) * gam3 
        KPcoeff=czero
        !
        KPcoeff(:,:,const) = czero
        !
        KPcoeff(:,:,dx) = czero
        KPcoeff(:,:,dy) = czero
        KPcoeff(:,:,dz) = czero
        !
        KPcoeff(:,:,dx2) = reshape( &
                                  (/ P+Q   ,    czero ,  -R1 ,  czero , &
                                     czero ,    P-Q   , czero,  -R1   , &
                                -conjg(R1) ,    czero ,  P-Q ,  czero , &
                                     czero ,-conjg(R1), czero,  P+Q     /),  (/4, 4/) )
        !
        KPcoeff(:,:,dy2) = reshape( &
                                  (/ P+Q   ,    czero ,   R1 ,  czero , &
                                     czero ,    P-Q   , czero,   R1   , &
                                 conjg(R1) ,    czero ,  P-Q ,  czero , &
                                     czero , conjg(R1), czero,  P+Q     /),  (/4, 4/) )
                                             
        KPcoeff(:,:,dz2) = reshape( &
                                  (/ P-2.0*Q ,    czero , czero   ,  czero ,  &
                                     czero   ,  P+2.0*Q , czero   ,  czero ,  &
                                     czero   ,    czero , P+2.0*Q ,  czero ,  &
                                     czero   ,    czero , czero   ,  P-2.0*Q  /),  (/4, 4/) )                                     
        KPcoeff(:,:,dxdy) = reshape( &
                                  (/ czero ,    czero ,  R2 ,  czero ,  &
                                     czero ,    czero , czero,  R2   ,  &
                                 conjg(R2) ,    czero , czero, czero ,  &
                                     czero , conjg(R2), czero, czero    /),  (/4, 4/) )
                                             
        KPcoeff(:,:,dxdz) = reshape( &
                                  (/ czero ,    -S ,  czero  ,  czero , &
                                -conjg(S)  , czero ,  czero  ,  czero , &
                                     czero , czero ,  czero  ,   S    , &
                                     czero , czero , conjg(S),  czero   /),  (/4, 4/) )
                                                                                          
        KPcoeff(:,:,dydz) = reshape( &
                                  (/ czero , c1i*S ,  czero     ,  czero , &
                              conjg(c1i*S) , czero ,  czero     ,  czero , &
                                     czero , czero ,  czero     ,  c1i*S , &
                                     czero , czero ,conjg(c1i*S),  czero   /),  (/4, 4/) )
                    
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
        real(8)::dx,dy,dz
        ! k_v -> -i d/d_v 
        ! use the following rules
        !   d/dx = (psi_i+1 - psi_i-1)/2/dx
        !   d^2/dx^2 = (psi_i+1 + psi_i-1 - 2psi_i)/dx^2
        !   d^2/dx/dy = (psi_i+1,j+1 + psi_i-1,j-1 - psi_i+1,j-1 - psi_i-1,j+1) / (4 dx dy)        
        ij = i-j        
        dx=dd(1)
        dy=dd(2)
        dz=dd(3)
        coeff(:) = czero
        ! onsite 
        if (all(ij == (/0,0,0/))) then
          coeff(const) = 1.0d0
          coeff(dx2) = -2.0d0/dx/dx * (-c1i)**2
          coeff(dy2) = -2.0d0/dy/dy * (-c1i)**2
          coeff(dz2) = -2.0d0/dz/dz * (-c1i)**2
        endif
        ! +x
        if (all(ij == (/1,0,0/))) then
          coeff(dx) = 1.0d0/2.0d0/dx * (-c1i)
          coeff(dx2) = 1.0d0/dx/dx   * (-c1i)**2          
        endif
        ! -x
        if (all(ij == (/-1,0,0/))) then
          coeff(dx) = -1.0d0/2.0d0/dx * (-c1i)
          coeff(dx2) = 1.0d0/dx/dx    * (-c1i)**2
        endif
        ! +y
        if (all(ij == (/0,1,0/))) then
          coeff(dy) = 1.0d0/2.0d0/dy  * (-c1i)
          coeff(dy2) = 1.0d0/dy/dy    * (-c1i)**2
        endif
        ! -y
        if (all(ij == (/0,-1,0/))) then
          coeff(dy) = -1.0d0/2.0d0/dy  * (-c1i)
          coeff(dy2) = 1.0d0/dy/dy     * (-c1i)**2
        endif
        ! +z
        if (all(ij == (/0,0,1/))) then
          coeff(dz) = 1.0d0/2.0d0/dz * (-c1i)
          coeff(dz2) = 1.0d0/dz/dz   * (-c1i)**2
        endif
        ! -z
        if (all(ij == (/0,0,-1/))) then
          coeff(dz) = -1.0d0/2.0d0/dz  * (-c1i)
          coeff(dz2) = 1.0d0/dz/dz     * (-c1i)**2
        endif
        ! +x +y
        if (all(ij == (/1,1,0/))) then
          coeff(dxdy) = 1.0d0/4.0d0/dx/dy * (-c1i)**2
        endif
        ! +x -y
        if (all(ij == (/1,-1,0/))) then
          coeff(dxdy) = -1.0d0/4.0d0/dx/dy * (-c1i)**2
        endif
        ! -x +y
        if (all(ij == (/-1,1,0/))) then
          coeff(dxdy) = -1.0d0/4.0d0/dx/dy * (-c1i)**2
        endif
        ! -x -y
        if (all(ij == (/-1,-1,0/))) then
          coeff(dxdy) = 1.0d0/4.0d0/dx/dy * (-c1i)**2
        endif        
        ! +x +z
        if (all(ij == (/1,0,1/))) then
          coeff(dxdz) = 1.0d0/4.0d0/dx/dz * (-c1i)**2
        endif
        ! +x -z
        if (all(ij == (/1,0,-1/))) then
          coeff(dxdz) = -1.0d0/4.0d0/dx/dz * (-c1i)**2
        endif 
        ! -x -z
        if (all(ij == (/-1,0,-1/))) then
          coeff(dxdz) = 1.0d0/4.0d0/dx/dz * (-c1i)**2
        endif 
        ! -x +z
        if (all(ij == (/-1,0,1/))) then
          coeff(dxdz) = -1.0d0/4.0d0/dx/dz * (-c1i)**2
        endif 
        ! +y +z
        if (all(ij == (/0,1,1/))) then
          coeff(dydz) = 1.0d0/4.0d0/dy/dz * (-c1i)**2
        endif 
        ! -y +z
        if (all(ij == (/0,-1,1/))) then
          coeff(dydz) = -1.0d0/4.0d0/dy/dz * (-c1i)**2
        endif 
        ! +y -z
        if (all(ij == (/0,1,-1/))) then
          coeff(dydz) = -1.0d0/4.0d0/dy/dz * (-c1i)**2
        endif 
        ! -y -z
        if (all(ij == (/0,-1,-1/))) then
          coeff(dydz) = 1.0d0/4.0d0/dy/dz * (-c1i)**2
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


end module kpHam
