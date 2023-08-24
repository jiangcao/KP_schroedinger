MODULE static   

IMPLICIT NONE

REAL(8), PARAMETER  :: DIEL_0=8.85418d-14 ! F/cm
!!$REAL(8), PARAMETER  :: DIEL_OX=3.9D0
!!$REAL(8), PARAMETER  :: DIEL_OX=6.7D0
REAL(8), PARAMETER  :: DIEL_METAL=0.0D0
!REAL(8), PARAMETER  :: Hk=4.0d0
REAL(8), PARAMETER  :: m0=5.6856D-16 !eV s^2 / cm^2   rest mass of electron
real(8), PARAMETER  :: m0ISO=9.10938D-28   ! gram   rest mass of electron
real(8), PARAMETER  :: gram2eV=m0/m0ISO    ! eV s^2 / cm^2 / gram
REAL(8), PARAMETER  :: hbar=6.58211899E-16 !eV s
!REAL(8), PARAMETER  :: TEMP=300.0d0 ! K
REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
REAL(8), PARAMETER  :: PII=3.14159265d0
REAL(8), PARAMETER  :: ELCH=1.60217653d-19   ! C
real(8), PARAMETER  :: v_light=2.998d10      ! cm/s
REAL(8), PARAMETER  :: x_mass=0.228d0
REAL(8), PARAMETER  :: y_mass=0.228d0
REAL(8), PARAMETER  :: z_mass=0.228d0
REAL(8), PARAMETER  :: x_mass_si=0.191d0
REAL(8), PARAMETER  :: y_mass_si=0.916d0
REAL(8), PARAMETER  :: z_mass_si=0.191d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(8), parameter :: Egap_WTe2 = 0.131d0  !! [wrt. CB of MoS_2]
real(8), parameter :: L_thick = 0.6D-7  !! [cm]

!real(8), parameter :: m_mos2_no_strain = 0.476D0
!real(8), parameter :: m_wte2_no_strain = 0.453D0
real(8), parameter :: m_mos2_strain = 0.391d0 ! 
real(8), parameter :: m_wte2_strain = 0.543d0 !
real(8), parameter :: E_delta1 = 0.0d0 ! -238.79D-3
real(8), parameter :: E_delta2 = 0.0d0 ! -45.66D-3

real(8), parameter :: m_top = m_mos2_strain  !! Used in Poisson
real(8), parameter :: m_bot = m_wte2_strain  !! Used in Poisson


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)


REAL(8), PARAMETER :: htqm=(hbar)**2/m0		 !/ELCH
real(8), parameter :: hb2m=7.6305d-16     	 ! eV*cm^2

real(8), parameter :: SCBA_error_seuil = 1.0D-10
integer, parameter :: N_max_fail = 3
integer, parameter :: N_max_LS = 3
real(8), parameter :: LineSearch_coef = 0.5d0

logical, parameter :: GMSH_DOUBLE_CHECK = .true.
logical, parameter :: FORCE_MOB_CHARG = .false.
logical, parameter :: RELAX_POISSON = .false.
logical, PARAMETER :: MASK_GREEN = .false.

END MODULE static
