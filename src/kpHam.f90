! Copyright (c) 2023 Jiang Cao, ETH Zurich 
! All rights reserved.
!
module kpHam

    implicit none 

    private
    
    public :: read_KP_coeff, save_KP_coeff
    public :: save_Ham_blocks, save_matrix, save_matrix_csr, load_matrix_csr
    public :: calc_bands_at, calc_KP_block, calc_kpbands_at
    public :: generate_LK4,generate_LK6
    public :: build_wire_ham_blocks, build_dot_ham
    public :: test_bandstructure, test_transport_bandstructure, test_schroedinger
    public :: eigv_feast
    public :: generate_BLG8
    public :: generate_LS3
    public :: rotate_basis, rotate_k_vector
    

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
    ! lookup table for the index of linear k-terms
    integer, parameter,dimension(3) :: lookup_linear_k  =  [dx, dy, dz]
    ! lookup table for the index of quadratic k-terms
    integer, parameter,dimension(3,3) :: lookup_quadratic_k  = reshape( [dx2,  dxdy, dxdz,&
                                                                         dxdy,  dy2, dydz,&
                                                                         dxdz, dydz,  dz2] , shape=[3,3] )

CONTAINS

    ! build the Hamiltonian blocks for a wire structure
    !   use KPcoeff. table, wire direction can be picked from x-y-z
    subroutine build_wire_ham_blocks(nbnd,KPcoeff, dd, wire_direction, Ny, Nz, H00, H10, vector_field)
        integer,intent(in)::nbnd ! number of bands
        complex(8), intent(in) :: KPcoeff(nbnd,nbnd,num_k) ! KP coeff. table        
        real(8), intent(in) :: dd(3) ! discretization step size in x-y-z
        integer, intent(in) :: wire_direction ! wire direction 1-3
        integer, intent(in) :: Ny,Nz ! number of points in the cross-section
        complex(8), intent(out), dimension(Ny*Nz*nbnd,Ny*Nz*nbnd) :: H00,H10 ! Ham blocks of the wire
        real(8), intent(in), optional :: vector_field(3)
        ! ----
        complex(8)::Hij(nbnd,nbnd),phi
        integer:: x(3),y(3),z(3),i(3),j(3), m,n,p,q, row,col
        real(8):: Avec(3)
        if(present(vector_field)) then
          Avec = vector_field
        else
          Avec = 0.0d0
        endif
        select case( wire_direction )
            case default  ! x
                x = (/1,0,0/)
                y = (/0,1,0/)
                z = (/0,0,1/)
            case(2) ! along y
                x = (/0,1,0/) 
                y = (/0,0,1/)
                z = (/1,0,0/)
            case(3) ! along z
                x = (/0,0,1/) 
                y = (/1,0,0/)
                z = (/0,1,0/)
        end select
        H00 = czero
        H10 = czero
        do m = 1,Ny
          do n = 1,Nz
            row = (m-1)*Nz + n
            do p = 1,Ny
              do q = 1,Nz
                col = (p-1)*Nz + q
                i = (m-1) * y + (n-1) * z
                j = (p-1) * y + (q-1) * z
                phi = dot_product(dble(j-i) * dd , Avec)
                call calc_KP_block(i,j,dd,nbnd,KPcoeff,Hij) 
                H00( (row-1)*nbnd+1:row*nbnd , (col-1)*nbnd+1:col*nbnd ) = Hij(:,:) * exp(-c1i*e0/hbar/twopi * phi)
                !
                j = (p-1) * y + (q-1) * z + x
                phi = dot_product(dble(j-i) * dd , Avec)
                call calc_KP_block(i,j,dd,nbnd,KPcoeff,Hij) 
                H10( (row-1)*nbnd+1:row*nbnd , (col-1)*nbnd+1:col*nbnd ) = Hij(:,:) * exp(c1i*e0/hbar/twopi * phi)
              enddo
            enddo
          enddo
        enddo
    end subroutine build_wire_ham_blocks

    
    ! build the dot Hamiltonian 
    subroutine build_dot_ham(nbnd,KPcoeff,dd,Nx,Ny,Nz,Ham,vector_field)
        integer,intent(in)::nbnd ! number of bands
        complex(8), intent(in) :: KPcoeff(nbnd,nbnd,num_k) ! KP coeff. table        
        real(8), intent(in) :: dd(3) ! discretization step size in x-y-z        
        integer, intent(in) :: Nx,Ny,Nz ! number of points
        complex(8), intent(out), dimension(Nx*Ny*Nz*nbnd,Nx*Ny*Nz*nbnd) :: Ham ! dot Hamiltonian        
        real(8), intent(in), optional :: vector_field(3)
        ! ----
        complex(8)::Hij(nbnd,nbnd), phi
        integer:: x(3),y(3),z(3),i(3),j(3), m,n,l,p,q,o, row,col
        real(8):: Avec(3)
        if(present(vector_field)) then
          Avec = vector_field
        else
          Avec = 0.0d0
        endif
        x = (/1,0,0/)
        y = (/0,1,0/)
        z = (/0,0,1/)
        do l = 1,Nx
          do m = 1,Ny
            do n = 1,Nz            
              row = (l-1)*Ny*Nz + (m-1)*Nz + n
              do o = 1,Nx
                do p = 1,Ny
                  do q = 1,Nz
                    col = (o-1)*Ny*Nz + (p-1)*Nz + q
                    i = (m-1) * y + (n-1) * z + (l-1) * x
                    j = (p-1) * y + (q-1) * z + (o-1) * x
                    phi = dot_product(dble(j-i) * dd , Avec)
                    call calc_KP_block(i,j,dd,nbnd,KPcoeff,Hij) 
                    Ham( (row-1)*nbnd+1:row*nbnd , (col-1)*nbnd+1:col*nbnd ) = Hij(:,:) * exp(c1i*e0/hbar/twopi * phi)                   
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
    end subroutine build_dot_ham 


    ! generate BiLayer Graphene 8-band KP model coefficients
    !   PRL 125 196402 (2020) , Supp Mat
    subroutine generate_BLG8(KPcoeff, tau,gam0,gam1,gam3,gam4,V,delta,lIA,lIB,lattice,correct_fdp,magfield)
      real(8),intent(in) :: tau,gam0,gam1,gam3,gam4,V,delta,lIA,lIB,lattice
      complex(8),intent(out) :: KPcoeff(8,8,num_k) ! KP coeff. table      
      real(8),intent(in) :: correct_fdp
      real(8),intent(in),optional :: magfield(3) ! magnetic field vector
      ! ----
      complex(8),dimension(2,2) :: sx,sy,sz,s0 ! spin operators S_x,S_y,S_z,S_+,S_-,I
      real(8) :: lIA1,lIA2,lIB1,lIB2,lexA1,lexA2,lexB1,lexB2,l0,lR
      complex(8) :: HTB(4,4,num_k) 
      complex(8) :: HSO(4,4) , HZ(2,2)
      complex(8) :: s(2,2,4,4) 
      complex(8) :: DD(num_k), VV(num_k), ff(num_k), G1(num_k),AA(num_k)
      integer :: i,m,n      
      real(8) :: Bvec(3) 
      real(8),parameter::g0 = 2.0d0
      real(8),parameter::muB = hb2m*(e0)/hbar
      if (present(magfield)) then          
          Bvec=magfield
      else
          Bvec=0.0d0          
      endif
      !
      lIA1=0.0d0
      lIA2=lIA
      lIB1=0.0d0
      lIB2=lIB
      !
      lexA1=0.0d0
      lexA2=0.0d0
      lexB1=0.0d0
      lexB2=0.0d0
      !
      l0=0.0d0
      lR=0.0d0
      !
      sx = 0.5d0 * reshape( (/0.0d0,1.0d0,1.0d0,0.0d0/) , (/2,2/) )
      sy = 0.5d0/c1i * reshape( (/0.0d0,1.0d0,-1.0d0,0.0d0/) , (/2,2/) )
      sz = 0.5d0 * reshape( (/1.0d0,0.0d0,0.0d0,-1.0d0/) , (/2,2/) )      
      s0 = reshape( (/1.0d0,0.0d0,0.0d0,1.0d0/) , (/2,2/) )
      !
      KPcoeff = czero
      DD = czero
      VV = czero
      ff = czero
      G1 = czero
      AA = czero
      !
      DD(const) = delta
      VV(const) = V      
      G1(const) = gam1
      !
      ff(dx) = - sqrt(3.0d0) / 2.0d0 * lattice * tau
      ff(dy) = + sqrt(3.0d0) / 2.0d0 * lattice * c1i       
      !
      ! To avoid Fermion doubling problem, refer to L. Susskind, Phys. Rev. D
      ! 16, 3031 (1977) and the discussion around Eq (1.7)
      if (correct_fdp .ne. 0) then
        AA(dx2) = -abs(ff(dx)) * correct_fdp
        AA(dy2) = -abs(ff(dx)) * correct_fdp        
      endif
      !
      ! H_TB      
      do i=1,num_k
        HTB(:,:,i) = reshape( (/ DD(i)+VV(i)+AA(i)  ,  gam0*ff(i)  ,   gam4*conjg(ff(i))  , G1(i)              ,  &
                                 gam0*conjg(ff(i))  ,  VV(i)-AA(i) ,   gam3*ff(i)         , gam4*conjg(ff(i))  ,  &
                                 gam4*ff(i)         ,  gam3*conjg(ff(i))   ,  -VV(i)+AA(i), gam0*ff(i)         ,  &
                                 G1(i)              ,  gam4*ff(i)  ,  gam0*conjg(ff(i))   , DD(i)-VV(i) -AA(i)   /), (/4,4/) )
      enddo
      do i=1,num_k
        do m=1,4
          do n=1,4
            KPcoeff( (m-1)*2+1:m*2 , (n-1)*2+1:n*2 , i ) = KPcoeff( (m-1)*2+1:m*2 , (n-1)*2+1:n*2 , i ) + HTB(m,n,i) * s0
          enddo
        enddo      
      enddo
      ! + H_SO 
      s = czero
      s(:,:,1,1) = sz
      s(:,:,1,2) = 0.5d0*(sx-c1i*tau*sy)
      s(:,:,2,1) = 0.5d0*(sx+c1i*tau*sy)
      s(:,:,2,2) = sz
      s(:,:,3,3) = sz
      s(:,:,3,4) = 0.5d0*(sx-c1i*tau*sy)
      s(:,:,4,3) = 0.5d0*(sx+c1i*tau*sy)
      s(:,:,4,4) = sz
      !
      HSO = reshape( (/  (tau*lIA1 - lexA1)*cone, c1i*(l0+2.0d0*lR)          , czero                , czero               , &
                        -c1i*(l0+2.0d0*lR)      , (-tau*lIB1-lexB1)*cone     , czero                , czero               , &
                        czero                   , czero                      , (tau*lIA2-lexA2)*cone,  -c1i*(l0-2.0d0*lR) , &
                        czero                   , czero                      , c1i*(l0-2.0d0*lR)    , (-tau*lIB2-lexB2)*cone /), (/4,4/) )
      !              
      do m=1,4
        do n=1,4
          KPcoeff( (m-1)*2+1:m*2 , (n-1)*2+1:n*2 , const ) = KPcoeff( (m-1)*2+1:m*2 , (n-1)*2+1:n*2 , const ) + HSO(m,n) * s(:,:,m,n)
        enddo
      enddo
      ! + H_Zeeman
      HZ = g0*muB*( Bvec(1)*sx + Bvec(2)*sy + Bvec(3)*sz )  
      do m=1,4        
        KPcoeff( (m-1)*2+1:m*2 , (m-1)*2+1:m*2 , const ) = KPcoeff( (m-1)*2+1:m*2 , (m-1)*2+1:m*2 , const ) + HZ(:,:)        
      enddo
    end subroutine generate_BLG8

    ! rotate the basis functions of the Hamiltonian
    !    given $\psi'_j = U_{ji} \psi_i$ 
    !    D' = U D U^T
    subroutine rotate_basis(nb,in_KPcoeff,out_KPcoeff,U)
      integer,intent(in)::nb ! number of bands
      complex(8),intent(in)::in_KPcoeff(nb,nb,num_k)
      complex(8),intent(out)::out_KPcoeff(nb,nb,num_k)
      real(8),intent(in)::U(nb,nb) ! rotation matrix (unitary)
      ! ----
      integer::i,j,m,n,mp,np,k
      !
      out_KPcoeff = czero      
      ! rotation for quadratic k terms  
      do k=5,10
        do mp=1,3
          do np=1,3
            do m=1,3
              do n=1,3                
                out_KPcoeff(mp,np,k) = out_KPcoeff(mp,np,k) + U(mp,m)*U(np,n)*in_KPcoeff(m,n,k)                              
              enddo
            enddo              
          enddo
        enddo            
      enddo
    end subroutine rotate_basis
    
    ! replace the original wave vector (unprimed one) by the one in the rotated frame (prime one)
    !    given $k'_j = U_{ji} k_i$ 
    subroutine rotate_k_vector(nb,in_KPcoeff,out_KPcoeff,U)
      integer,intent(in)::nb ! number of bands
      complex(8),intent(in)::in_KPcoeff(nb,nb,num_k)
      complex(8),intent(out)::out_KPcoeff(nb,nb,num_k)
      real(8),intent(in)::U(3,3) ! rotation matrix (unitary)
      ! ----
      integer::i,j,ip,jp,mp,np,k,kp
      !
      out_KPcoeff = czero
      !out_KPcoeff(:,:,const)=in_KPcoeff(:,:,const)
      do ip=1,3        
        do i=1,3
          ! rotation for linear k terms
          kp=lookup_linear_k(ip)
          k =lookup_linear_k(i)
          do mp=1,nb
            do np=1,nb
              out_KPcoeff(mp,np,kp)=out_KPcoeff(mp,np,kp) + in_KPcoeff(mp,np,k)*U(ip,i)
            enddo
          enddo
          ! rotation for quadratic k terms  
          do jp=1,3                            
            do j=1,3             
              kp=lookup_quadratic_k(ip,jp)
              k =lookup_quadratic_k(i,j)
              do mp=1,nb
                do np=1,nb
                  if (i /= j) then 
                    out_KPcoeff(mp,np,kp)=out_KPcoeff(mp,np,kp) + in_KPcoeff(mp,np,k)/2.0d0*U(ip,i)*U(jp,j)
                  else 
                    out_KPcoeff(mp,np,kp)=out_KPcoeff(mp,np,kp) + in_KPcoeff(mp,np,k)*U(ip,i)*U(jp,j)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    end subroutine rotate_k_vector


    ! generate LS 3-band KP model coefficients
    !   
    subroutine generate_LS3(KPcoeff,gam1,gam2,gam3,delta,kap_,qB_,Bvec_)        
        real(8), intent(in) :: gam1,gam2,gam3,delta        
        complex(8), intent(out) :: KPcoeff(3,3,num_k) ! KP coeff. table        
        real(8), intent(in),optional ::kap_,qB_
        real(8), intent(in),optional :: Bvec_(3)
        ! ----
        real(8)::Bvec(3),kap,qB
        integer::i,m,n  
        real(8):: Bx,By,Bz
        real(8) :: A,B,C
        A = hb2m * (gam1 + 4.0* gam2)
        B = hb2m * (gam1 - 2.0* gam2)
        C = hb2m * 6.0 * gam3
        Bvec= merge(Bvec_,0.0d0,present(Bvec_))
        kap= merge(kap_,0.0d0,present(kap_))
        qB= merge(qB_,0.0d0,present(qB_))
        print *
        print *,'B=',Bvec
        Bx=Bvec(1)
        By=Bvec(2)
        Bz=Bvec(3)
        ! construct the coeff. table   
        KPcoeff=czero
        !
        KPcoeff(1,1,dx2)=A
        KPcoeff(1,1,dy2)=B
        KPcoeff(1,1,dz2)=B
        !
        KPcoeff(2,2,dx2)=B
        KPcoeff(2,2,dy2)=A
        KPcoeff(2,2,dz2)=B
        !
        KPcoeff(3,3,dx2)=B
        KPcoeff(3,3,dy2)=B
        KPcoeff(3,3,dz2)=A
        !
        KPcoeff(2,1,dxdy)=C
        !
        KPcoeff(3,1,dxdz)=C
        !
        KPcoeff(3,2,dydz)=C
        do n=1,3
          do m=1,n-1
            do i=1,num_k
              KPcoeff(m,n,i) = conjg( KPcoeff(n,m,i) )               
            enddo
          enddo
        enddo
    end subroutine generate_LS3


    ! generate Luttinger 6-band KP model coefficients
    !   PRB 98, 155319 (2018) Eq (D1) (D4), correct an error in the paper in the R term
    subroutine generate_LK6(KPcoeff, gam1,gam2,gam3,delta,kap_,qB_,Bvec_)        
        real(8), intent(in) :: gam1,gam2,gam3,delta
        real(8), intent(in),optional ::kap_,qB_
        complex(8), intent(out) :: KPcoeff(6,6,num_k) ! KP coeff. table        
        real(8), intent(in),optional :: Bvec_(3)
        ! ----
        complex(8),dimension(num_k) :: P,Q,S,R,D
        complex(8) :: kkx,kky,kkz,qx,qy,qz
        real(8):: kapp,kappp
        complex(8),dimension(6,6)::Jx,Jy,Jz
        complex(8),dimension(6,6)::Kx,Ky,Kz
        complex(8),dimension(6,6)::Jx3,Jy3,Jz3 
        real(8)::Bvec(3),kap,qB
        integer::i,m,n  
        real(8):: Bx,By,Bz
        real(8),parameter :: sq2 = sqrt(2.0d0)
        real(8),parameter :: sq3o2 = sqrt(3.0d0/2.0d0)
        real(8),parameter :: sq3b2 = sqrt(3.0d0)/2.0d0
        Bvec= merge(Bvec_,0.0d0,present(Bvec_))
        kap= merge(kap_,0.0d0,present(kap_))
        qB= merge(qB_,0.0d0,present(qB_))
        print *
        print *,'B=',Bvec
        Bx=Bvec(1)
        By=Bvec(2)
        Bz=Bvec(3)
        ! construct the coeff. table   
        P=czero
        Q=czero
        S=czero
        R=czero
        D=czero
        !
        Jx3=czero
        Jy3=czero
        Jz3=czero
        ! 
        kapp=1.0d0+kap
        kappp=1.0d0+2.0d0*kap
        !
        Kx(:,:) = reshape( &
               (/ 0.0d0 , -sq3b2*kap,       0.0d0,      0.0d0,    sq3o2/2.0d0*kapp,                  0.0d0, &
                  0.0d0 ,     0.0d0 ,  -1.0d0*kap,      0.0d0,               0.0d0,   1.0d0/2.0d0/sq2*kapp, &
                  0.0d0 ,     0.0d0 ,      0.0d0 , -sq3b2*kap, -1.0/2.0d0/sq2*kapp,                  0.0d0, &
                  0.0d0 ,     0.0d0 ,      0.0d0 ,     0.0d0 ,              0.0d0 ,      -sq3o2/2.0d0*kapp, &
                  0.0d0 ,     0.0d0 ,      0.0d0 ,     0.0d0 ,              0.0d0 ,           -0.5d0*kappp, &
                  0.0d0 ,     0.0d0 ,      0.0d0 ,     0.0d0 ,              0.0d0 ,                0.0d0/), &
                (/6, 6/) )
        !
        Ky(:,:) = -c1i * reshape( &
               (/ 0.0d0, -sq3b2*kap ,      0.0d0 ,      0.0d0,  sq3o2/2.0d0*kapp,               0.0d0, &
                  0.0d0,      0.0d0 , -1.0d0*kap ,      0.0d0,             0.0d0,    1/2.0d0/sq2*kapp, &
                  0.0d0,      0.0d0 ,      0.0d0 , -sq3b2*kap,  1/2.0d0/sq2*kapp,               0.0d0, &
                  0.0d0,      0.0d0 ,      0.0d0 ,     0.0d0 ,            0.0d0 ,    sq3o2/2.0d0*kapp, &
                  0.0d0,      0.0d0 ,      0.0d0 ,     0.0d0 ,            0.0d0 ,        -0.5d0*kappp, &
                  0.0d0,      0.0d0 ,      0.0d0 ,     0.0d0 ,            0.0d0 ,             0.0d0/), &
                (/6, 6/) )
        !
        Kz(:,:) = reshape( &
               (/ -sq3o2**2*kap,             0.0d0,          0.0d0 ,         0.0d0,           0.0d0,            0.0d0, &
                          0.0d0, -1.0d0/2.0d0*kap ,          0.0d0 ,         0.0d0, -1.0d0/sq2*kapp,            0.0d0, &
                          0.0d0,            0.0d0 , 1.0d0/2.0d0*kap,         0.0d0,          0.0d0 ,  -1.0d0/sq2*kapp, &
                          0.0d0,            0.0d0 ,          0.0d0 ,  sq3o2**2*kap,          0.0d0 ,            0.0d0, &
                          0.0d0,            0.0d0 ,          0.0d0 ,         0.0d0,    -0.5d0*kappp,            0.0d0, &
                          0.0d0,            0.0d0 ,          0.0d0 ,         0.0d0,          0.0d0 ,      0.5d0*kappp/),&
                (/6, 6/) )
        !
        do n=1,6
          do m=1,n-1
                Kx(m,n) = conjg( Kx(n,m) )
                Ky(m,n) = conjg( Ky(n,m) )
                Kz(m,n) = conjg( Kz(n,m) )               
          enddo
        enddo
        !
        Jx(1:4,1:4)=Kx(1:4,1:4)/kap
        Jy(1:4,1:4)=Ky(1:4,1:4)/kap
        Jz(1:4,1:4)=Kz(1:4,1:4)/kap
        !
        Jx3(1:4,1:4)=matmul(matmul(Jx(1:4,1:4),Jx(1:4,1:4)),Jx(1:4,1:4))
        Jy3(1:4,1:4)=matmul(matmul(Jy(1:4,1:4),Jy(1:4,1:4)),Jy(1:4,1:4))
        Jz3(1:4,1:4)=matmul(matmul(Jz(1:4,1:4),Jz(1:4,1:4)),Jz(1:4,1:4))
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
        kkx = hb2m*2.0d0*(e0)/hbar*Bx        
        qx = hb2m*2.0d0*(e0)/hbar*qB*Bx
        !
        kky = hb2m*2.0d0*(e0)/hbar*By
        qy = hb2m*2.0d0*(e0)/hbar*qB*By
        !
        kkz = hb2m*2.0d0*(e0)/hbar*Bz
        qz = hb2m*2.0d0*(e0)/hbar*qB*Bz
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
        !
        KPcoeff(:,:,const)= KPcoeff(:,:,const) &
                           + Kx*kkx + Jx3*qx  &
                           + Ky*kky + Jy3*qy  &
                           + Kz*kkz + Jz3*qz 
    end subroutine generate_LK6



    ! generate Luttinger 4-band KP model coefficients
    !   PRB 98, 155319 (2018) Eq (D1) (D4)
    subroutine generate_LK4(KPcoeff, gam1,gam2,gam3,kap,qB,Bvec)        
        real(8), intent(in) :: gam1,gam2,gam3
        real(8), intent(in),optional::kap,qB
        real(8), intent(in),optional::Bvec(3)
        complex(8), intent(out) :: KPcoeff(4,4,num_k) ! KP coeff. table        
        ! ----
        complex(8)::LK6(6,6,num_k)
        !
        KPcoeff=czero
        call generate_LK6(LK6,gam1,gam2,gam3,0.0d0,kap,qB,Bvec)
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
    
    
    ! save a complex matrix to a file in row-column-value format
    subroutine save_matrix(filename, n, m, Mat)
        character(len=*),intent(in)::filename ! file name
        integer,intent(in)::n,m ! size
        complex(8), intent(in) :: Mat(n,m) ! matrix        
        ! ----
        integer::i,j        
        open(unit=11, file=filename, status='unknown')                
        do i=1,n
            do j=1,m
              if ( abs(Mat(i,j)) > 0.0d0 ) then
                write(11,'(2I10,2E18.6)') i,j,dble(Mat(i,j)),aimag(Mat(i,j))
              endif
            enddo            
        enddo        
        close(11)
    end subroutine save_matrix
    
    ! save a complex sparse matrix to a file in CSR format
    subroutine save_matrix_csr(filename, n, m, Mat)
        character(len=*),intent(in)::filename ! file name
        integer,intent(in)::n,m ! size
        complex(8), intent(in) :: Mat(n,m) ! matrix           
        ! ----
        integer::i,j,nnz,ind(n),k        
        nnz = count(Mat /= czero) ! number of nonzero elements        
        open(unit=11, file=filename, status='unknown')   
        write(11,*) nnz          
        k=1   
        do i=1,n
            ind(i)=k      
            do j=1,m
                if (mat(i,j) /= czero) then
                  write(11,'(1I12,2E18.6)')j, dble(Mat(i,j)),aimag(Mat(i,j))
                  k=k+1
                endif                
            enddo                  
        enddo    
        write(11,'(I12)') ind(:)
        close(11)
    end subroutine save_matrix_csr
    
    
    ! load a complex sparse matrix from a file in CSR format
    subroutine load_matrix_csr(filename, n, m, Mat)
        character(len=*),intent(in)::filename ! file name
        integer,intent(in)::n,m ! size
        complex(8), intent(out) :: Mat(n,m) ! matrix           
        ! ----
        integer::i,j,nnz,ind(n),k           
        integer,allocatable::col(:)
        complex(8),allocatable::csr(:)
        real(8)::re,im     
        open(unit=11, file=filename, status='unknown')   
        read(11,*) nnz          
        allocate(csr(nnz))
        allocate(col(nnz))
        do i=1,nnz
          read(11,*)j, re, im
          csr(i)=dcmplx(re,im)
          col(i)=j
        enddo
        do i=1,n
            read(11,'(I12)') ind(i)                
        enddo
        close(11)
        Mat=czero
        do i=1,n-1                
            do j=ind(i),ind(i+1)-1
                Mat(i,col(j))=csr(j)
            enddo
        enddo    
        do j=ind(n),nnz
            Mat(n,col(j))=csr(j)
        enddo
        deallocate(col,csr)
    end subroutine load_matrix_csr

    
    ! save a wavefunction 
    subroutine save_wavefunc(filename, n, m, Mat, label)
        character(len=*),intent(in)::filename ! file name
        integer,intent(in)::n,m ! size
        complex(8), intent(in) :: Mat(n,m) ! matrix 
        integer,intent(in),optional::label
        ! ----
        integer::i,j        
        open(unit=11, file=filename, status='unknown')                
        do i=1,n
            do j=1,m
                if (present(label)) then
                    write(11,'(3I10,1E18.6)') i,j,label,abs(Mat(i,j))
                else
                    write(11,'(2I10,1E18.6)') i,j,abs(Mat(i,j))
                endif
            enddo
            write(11,*)
        enddo        
        close(11)
    end subroutine save_wavefunc


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
          coeff(dx2) = -2.0d0/ddx**2 * (-c1i)**2
          coeff(dy2) = -2.0d0/ddy**2 * (-c1i)**2
          coeff(dz2) = -2.0d0/ddz**2 * (-c1i)**2
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
        if (INFO/=0)then
        write(*,*)'SEVERE WARNING: ZHEEV HAS FAILED. INFO=',INFO
        call abort
        endif
        eig(:)=W(:)
    END FUNCTION eig
    
    
    
    ! calculate all eigen-values and eigen-vectors of a Hermitian matrix A 
    !   within a given search interval, a wrapper to the FEAST function in MKL https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-1/feast-syev-feast-heev.html   
    !   upon return A(:,1:m) will be modified and contains the eigen-vectors
    FUNCTION eigv_feast(NN, A, emin, emax, m, m_init)    
        include 'mkl.fi'
        INTEGER, INTENT(IN) :: NN
        COMPLEX(8), INTENT(INOUT), DIMENSION(:,:) :: A
        REAL(8), INTENT(IN) :: emin, emax ! lower and upper bounds of the interval to be searched for eigenvalues
        REAL(8) :: eigv_feast(NN)
        integer,intent(out) :: m ! total number of eigenvalues found
        integer,intent(in),optional :: m_init
        ! -----        
        real(8) :: epsout
        integer :: fpm(128), m0, loop, info
        complex(8),allocatable :: x(:,:)
        real(8), allocatable :: w(:), res(:)
        if (present(m_init)) then
            m0=m_init
        else
            m0=nn/2
        endif
        allocate(x(nn,m0))
        allocate(w(m0))
        allocate(res(m0))
        !
        call feastinit (fpm)
        fpm(1)=1 ! print runtime status to the screen
        !
        call zfeast_heev('U',nn,A,nn,fpm,epsout,loop,emin,emax,m0,W,x,m,res, info)        
        !
        if (INFO/=0)then
        write(*,*)'SEVERE WARNING: zfeast_heev HAS FAILED. INFO=',INFO
        call abort
        endif
        eigv_feast(1:m)=W(1:m)
        A(:,1:m) = x(:,1:m)
        deallocate(x,w,res)
    END FUNCTION eigv_feast
    
    
    
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
        if (INFO/=0)then
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
        real(8)::kpt(3),kpt1(3),kpt2(3),kmin,kmax
        real(8)::k_path(6,7)
        real(8),allocatable::Ek(:)
        integer :: nk,ik,ipath            
        nk = 1000
        kmin= 2.0d9
        kmax= 2.0d9
        k_path=reshape( &
               [ -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, &
                  0.0,-1.0, 0.0, 0.0, 1.0, 0.0, &
                  0.0, 0.0,-1.0, 0.0, 0.0, 1.0, &
                 -1.0,-1.0, 0.0, 1.0, 1.0, 0.0, &
                 -1.0, 0.0,-1.0, 1.0, 0.0, 1.0, &
                  0.0,-1.0,-1.0, 0.0, 1.0, 1.0, &
                 -1.0,-1.0,-1.0, 1.0, 1.0, 1.0] , shape=[6,7] )
        allocate(Ek(nbnd))
        !
        print *
        print *, ' Computing bulk bandstructure ...'
        !
        open(unit=10,file='ek0.dat',status='unknown')
        open(unit=11,file='ek.dat',status='unknown')
        !
        do ipath=1,7        
          !
          kpt1(:) = k_path(1:3,ipath) / sqrt( sum(k_path(1:3,ipath)**2) ) * kmin
          
          print *, kpt1
          kpt2(:) = k_path(4:6,ipath) / sqrt( sum(k_path(4:6,ipath)**2) ) * kmax
          
          print *, kpt2
          do ik = 1,nk
            kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1          
            call calc_bands_at(nbnd,KPcoeff,dd,kpt,Ek)          
            write(11,'(20E18.6)') dble(ik)/dble(nk)+ipath-1, Ek          
            call calc_kpbands_at(nbnd, KPcoeff,kpt,Ek)          
            write(10,'(20E18.6)') dble(ik)/dble(nk)+ipath-1, Ek          
          enddo
        enddo
        !
        close(10)
        close(11)
        deallocate(Ek)
    end subroutine test_bandstructure
    
    
    subroutine test_transport_bandstructure(H00,H10,nbnd,ny,nz,dx,emin,emax)
        implicit none
        complex(8),intent(in),dimension(ny*nz*nbnd,ny*nz*nbnd) :: H00, H10
        integer,intent(in) :: ny,nz,nbnd
        real(8),intent(in) :: dx,emin,emax
        ! ----
        complex(8)::Ham(ny*nz*nbnd,ny*nz*nbnd),V(nbnd,nz,ny)
        real(8)::kx,Ek(ny*nz*nbnd),kpt1,kpt2,kpt
        integer::nk,ik,ib,m,im,i,j,fu
        character(len=100)::filename
        !
        nk = 150
        print *
        print *, ' Computing wire bandstructure ...'
        open(unit=10,file='Boundary_ek.dat',status='unknown')        
        !
        !kpt1 = -2.0d9
        !kpt2 = +2.0d9 
        kpt1 = -3.1415d0/dx
        kpt2 = +3.1415d0/dx
        do ik = 1,nk
          print *
          print *, ik, '/', nk
          kpt= (kpt2-kpt1) * dble(ik)/dble(nk) + kpt1          
          Ham(:,:) = H00(:,:) + exp( - c1i * kpt * dx )*H10(:,:) &
                              + exp( + c1i * kpt * dx )*transpose(conjg(H10(:,:)))
          !
          if ( ny*nz*nbnd < 3000 ) then
            Ek(:) = eig(ny*nz*nbnd,Ham)
            m = ny*nz*nbnd
          else
            Ek(:) = eigv_feast(ny*nz*nbnd,Ham, emin=emin, emax=emax, m=m)  
          endif
          !
          do ib=1,m
            write(10,'(2E25.14)') dble(ik)/dble(nk), Ek(ib)
          enddo
        enddo            
        close(10)
        kpt=0.0d0
        Ham(:,:) = H00(:,:) + H10(:,:) + transpose(conjg(H10(:,:)))        
        if ( ny*nz*nbnd < 3000 ) then
          Ek(:) = eigv(ny*nz*nbnd,Ham)
          m = ny*nz*nbnd
        else
          Ek(:) = eigv_feast(ny*nz*nbnd,Ham, emin=emin, emax=emax, m=m)
        endif
        ! save eigen vectors
        do im=1,m
            if ((Ek(im)>emin) .and. (Ek(im)<emax)) then
                V = reshape( Ham(:,im) ,(/nbnd, nz, ny/) )
                filename = 'psi_'//string(im)//'.dat'
                open(newunit=fu, file=filename, status='unknown')                
                write(fu,*) '# energy', Ek(im)
                write(fu,*) '# Y - Z - Component - Abs(psi)'
                do ib=1,nbnd
                    write(fu,*) '# component', ib
                    do i=1,ny
                        do j=1,nz
                            write(fu,'(3I10,1E18.6)') i,j,ib,abs(V(ib,j,i))
                        enddo
                    enddo
                enddo        
                close(fu)
            endif
        enddo        
    end subroutine test_transport_bandstructure
    
    
    subroutine test_schroedinger(H,nbnd,nx,ny,nz,emin,emax,num_modes)
        implicit none
        complex(8),intent(inout),dimension(nx*ny*nz*nbnd,nx*ny*nz*nbnd) :: H
        integer,intent(in) :: nx,ny,nz,nbnd,num_modes        
        real(8),intent(in) :: emin, emax
        ! ----
        complex(8)::V(nbnd,nz,ny,nx)
        real(8)::Ek(nx*ny*nz*nbnd)
        integer::m,ib,im,m0
        !        
        m0=nbnd*nx*ny*nz/5
        print *
        print *, ' Computing Schoedinger ...'
        !        
        Ek(:) = eigv_feast(nx*ny*nz*nbnd,H, emin=emin, emax=emax, m=m, m_init=m0)
        open(unit=10,file='dot_en.dat',status='unknown')  
        write(10,'(1E25.18)') Ek(1:m)      
        close(10)
        do im=1,min(m,num_modes)
          ! map eigen vector to real-space
          V = reshape( H(:,im) ,(/nbnd, nz, ny, nx/) )
          do ib=1,nbnd        
              !call save_wavefunc( 'dot_'//string(im)//'_'//string(ib)//'_yz.dat', nz,ny, V(ib,:,:,max(nx/2,1)) )
          enddo        
          !        
          do ib=1,nbnd
              call save_wavefunc( 'dot_'//string(im)//'_'//string(ib)//'_xy.dat', nx,ny, V(ib,max(nz/2,1),:,:) )
          enddo        
          !
        enddo
        print *, Ek(1:m)
    end subroutine test_schroedinger
    
!!!!!!!!!!!!!!!!!!!!!!!!  utility functions  !!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    FUNCTION STRING(inn)
      IMPLICIT NONE
      INTEGER, PARAMETER :: POS= 4
      INTEGER, INTENT(IN) :: inn
      CHARACTER(LEN=POS) :: STRING
      !............................................................
      INTEGER :: cifra, np, mm, num  
      IF (inn > (10**POS)-1) stop "ERRORE: (inn > (10**3)-1)  in STRING"
      num= inn
      DO np= 1, POS
         mm= pos-np
         cifra= num/(10**mm)            
         STRING(np:np)= ACHAR(48+cifra)
         num= num - cifra*(10**mm)
      END DO
    END FUNCTION STRING
    
end module kpHam
