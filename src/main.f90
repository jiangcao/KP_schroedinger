! Copyright (c) 2023 Jiang Cao, ETH Zurich 
! All rights reserved.
!
PROGRAM main
    use kpHam
    use math
    implicit none

    integer :: nbnd
    real(8)::gam0,gam1,gam2,gam3,gam4,delta,kappa,V,tau,lattice,lIA,lIB,dd(3),kpt(3),Ek(6),kpt1(3),kpt2(3),qB,Bvec(3),emin,emax,U0
    complex(8),allocatable::KPcoeff(:,:,:),rot_KPcoeff(:,:,:),H00(:,:),H10(:,:)
    integer :: nk, nx, ny, nz, direction
    integer :: num_cpu, num_modes
    integer :: iy, ix, iz, ii
    integer :: fu,rc
    character(len=10) :: kp_type
    
    real(8) :: rot_mat(3,3), rot_theta, rot_phi

    namelist /input/ rot_theta,rot_phi, rot_mat,kp_type,nbnd, gam0, gam1, gam2, gam3, gam4, kappa, delta, V, lIA, lIB, lattice, tau, direction,Nx,Ny,Nz,dd, num_cpu, qB, Bvec, emin, emax,U0,num_modes

    tau=1.0d0
    emin=0.0d0
    emax=0.25d0
    U0=0.0d0
    rot_mat=zeye(3)
    rot_theta=0.0d0
    rot_phi=0.0d0

    open (action='read', file='input', iostat=rc, newunit=fu)
    read (nml=input, iostat=rc, unit=fu)
    close(fu)
    
    call normalize_real_matrix(rot_mat, 3)

    call omp_set_num_threads(num_cpu)

    print *, 'nbnd=',nbnd
    
    call alloc_KPcoeff(nbnd, KPcoeff)
    
    select case( trim(kp_type) )
    
        case( 'LS' )
            print *, rot_mat
            call generate_LS3(KPcoeff, gam1=gam1,gam2=gam2,gam3=gam3,delta=delta,kap_=kappa,qB_=qB,Bvec_=Bvec)  
            call save_KP_coeff('kpcoeff_unrotated.dat', nbnd, KPcoeff)            
            !
            call alloc_KPcoeff(nbnd, rot_KPcoeff)
            call rotate_basis(nb=3,in_KPcoeff=KPcoeff,out_KPcoeff=rot_KPcoeff,U=rot_mat)
            call save_KP_coeff('kpcoeff_rot_stage1.dat', nbnd, rot_KPcoeff)
            !            
            call rotate_k_vector(nb=nbnd,in_KPcoeff=rot_KPcoeff,out_KPcoeff=KPcoeff,U=rot_mat)
            call save_KP_coeff('kpcoeff_rot_stage2.dat', nbnd, KPcoeff)            
            
        case( 'LK4' )
            call generate_LK4(KPcoeff, gam1,gam2,gam3,kappa,qB,Bvec)

        case( 'LK6' )    
            call generate_LK6(KPcoeff, gam1,gam2,gam3,delta,kappa,qB,Bvec)
        
        case( 'BLG8' )    
            call generate_BLG8(KPcoeff, tau,gam0,gam1,gam3,gam4,V,delta,lIA,lIB,lattice,correct_fdp=lattice/20.0,magfield=Bvec)
            
        case default
            call generate_LK4(KPcoeff, gam1,gam2,gam3,kappa,qB,Bvec)
    
    end select

    call save_KP_coeff('kpcoeff.dat', nbnd, KPcoeff)

    call save_Ham_blocks('Ham_blocks.dat', nbnd, KPcoeff,dd)

    call test_bandstructure(nbnd,KPcoeff,dd)

    allocate(H00(Ny*Nz*nbnd,Ny*Nz*nbnd))
    allocate(H10(Ny*Nz*nbnd,Ny*Nz*nbnd))

    call build_wire_ham_blocks(nbnd,KPcoeff, dd, direction, Ny, Nz, H00, H10,vector_field=[-Bvec(2),Bvec(1),0.0d0])

    ! add a Gaussian potential
    do iy = 1,Ny
      do iz = 1,Nz
        do ii = ((iy-1) * Nz + iz - 1) * nbnd + 1 , ((iy-1) * Nz + iz) * nbnd
            H00(ii,ii) = H00(ii,ii) + gaussian(r=(iy-Ny/2)*dd(2),b=0.0d0,U0=U0,sigma=Ny*dd(2)/8)
        enddo
      enddo
    enddo

    ! call save_matrix('Hii.dat', Ny*Nz*nbnd, H00)
    ! call save_matrix('H1i.dat', Ny*Nz*nbnd, H10)

    call save_matrix_csr('Hii_csr.dat', Ny*Nz*nbnd, Ny*Nz*nbnd, H00)
    call save_matrix_csr('H1i_csr.dat', Ny*Nz*nbnd, Ny*Nz*nbnd, H10)

    call test_transport_bandstructure(H00,H10,nbnd,ny,nz,dd(direction),emin,emax)

    deallocate(H10,H00)

    allocate(H00(Nx*Ny*Nz*nbnd,Nx*Ny*Nz*nbnd))

    call build_dot_ham(nbnd,KPcoeff,dd,Nx,Ny,Nz,H00,vector_field=[-Bvec(2),Bvec(1),0.0d0])

    ! add a Gaussian potential
    open(unit=11, file='potential.dat', status='unknown') 
    do ix = 1,Nx
        do iy = 1,Ny
          write(11,*) ix, iy, gaussian(r=sqrt(((iy-Ny/2)*dd(2))**2 + ((ix-Nx/2)*dd(1))**2),b=0.0d0,U0=U0,sigma=Ny*dd(2)/8)
          do iz = 1,Nz
            do ii = ((ix-1)*Ny*Nz + (iy-1)*Nz + iz - 1) * nbnd + 1 , ((ix-1)*Ny*Nz + (iy-1)*Nz + iz) * nbnd
                H00(ii,ii) = H00(ii,ii) + gaussian(r=sqrt(((iy-Ny/2)*dd(2))**2 + ((ix-Nx/2)*dd(1))**2),b=0.0d0,U0=U0,sigma=Ny*dd(2)/8)
            enddo
          enddo
        enddo
    enddo
    close(11)

    call save_matrix_csr('Ham_csr.dat', Nx*Ny*Nz*nbnd, Nx*Ny*Nz*nbnd, H00)

    ! call test_schroedinger(H00,nbnd,nx,ny,nz,emin,emax,num_modes)

    deallocate(H00)

    deallocate(KPcoeff)

END PROGRAM main
