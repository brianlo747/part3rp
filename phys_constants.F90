!-*- F90 -*- so emacs thinks this is an f90 file
module phys_constants
  !--------------------------------------------------------------------!
  ! author:  Michael Herzog                                            !
  ! email:   Michael.Herzog@noaa.gov                                   !
  ! date:    May 2005                                                  !
  ! version: 0.1                                                       !
  !                                                                    !
  ! purpose: define physical constants                                 !
  !                                                                    !
  !--------------------------------------------------------------------!
  use precision
  implicit none
  public
  !--------------------------------------------------------------------!
  ! physical constants                                                 !
  !                                                                    !
  ! avno     Avagadro constant            [molec/mol]                  !
  ! rgas     gas constant                 [J/mol/K]                    !
  ! boltz    Boltzmann constant           [J/molec/K]                  !
  ! radair   mean radius of air molecule  [m]                          !
  !                                                                    !
  ! grav     gravitational constant       [m2/s]                       !
  ! rearth   radius of earth              [m]                          !
  ! omega    angular velocity of earth    [s-1]                        !
  ! ps0      reference pressure for z=0   [Pa]                         !
  !--------------------------------------------------------------------!
  real(kreal), parameter ::                 &
       pi     =3.141592653589793_kreal,     &
       epsilon=1.e-30_kreal,                &
       avno   =6.02252e23_kreal,            &
       rgas   =8.3144621_kreal,             &
       boltz  =1.380485e-23_kreal,          &
       radair =1.82525e-10_kreal
#ifdef MARS
  real(kreal), parameter ::                 &
       grav   =3.711_kreal,                 &
       rearth =3.3895e6_kreal,              &
       omega  =7.088271e-5_kreal,           &
       ps0    =351.199_kreal
!!$       ps0    =661.07_kreal
#else
  real(kreal), parameter ::                 &
       grav   =9.80616_kreal,               &
       rearth =6.371229e6_kreal,            &
       omega  =7.29211e-5_kreal

  ! XXX: quick fix here to make it possible to modify reference pressure
  real(kreal) :: ps0    =101300._kreal
#endif
  !--------------------------------------------------------------------!
  ! thermodynamic properties of dry air and h2o                        !
  !                                                                    !
  ! rgasair  gas constant of dry air                         [J/kg/K]  !
  ! cpair    specific heat of dry air at constant pressure   [J/kg/K]  !
  ! cvair    specific heat of dry air at constant volume     [J/kg/K]  !
  ! rgashum  gas constant of water vapor                     [J/kg/K]  !
  ! cphum    specific heat of water vapor at const pressure  [J/kg/K]  !
  ! cvhum    specific heat of water vapor at const volume    [J/kg/K]  !
  ! heatlat  latent heat (gas to liquid)                     [J/kg]    !
  ! heatfrz  latent heat for freezing                        [J/kg]    !
  ! heatsub  latent heat for sublimation                     [J/kg]    !
  !--------------------------------------------------------------------!
#ifdef MARS
  real(kreal), parameter ::                 &
!!$       cpair  =1004.64_kreal,               &   ! Air = dry air
!!$       cvair  = 717.60_kreal,               &
       cpair  =844._kreal,                  &   ! Air = CO2
       cvair  =655._kreal,                  &
       rgasair=cpair-cvair,                 & 
       cpco2  =844._kreal,                  &   ! CO2
       cvco2  =655._kreal,                  &
       rgasco2=cpco2-cvco2,                 &
       cphum  =1864._kreal,                 &   ! H2O
       cvhum  =1402.5_kreal,                &
       rgashum=cphum-cvhum,                 &
       cpso2  =655.6_kreal,                 &   ! SO2
       cvso2  =512.0_kreal,                 &
       rgasso2=cpso2-cvso2
#else
  real(kreal), parameter ::                 &
       cpair  =1004.64_kreal,               &   ! Air = dry air
       cvair  = 717.60_kreal,               &
       rgasair=cpair-cvair,                 & 
       cpco2  =844._kreal,                  &   ! CO2
       cvco2  =655._kreal,                  &
       rgasco2=cpco2-cvco2,                 &
       cphum  =1864._kreal,                 &   ! H2O
       cvhum  =1402.5_kreal,                &
       rgashum=cphum-cvhum,                 &
       cpso2  =655.6_kreal,                 &   ! SO2
       cvso2  =512.0_kreal,                 &
       rgasso2=cpso2-cvso2
#endif
  real(kreal), parameter ::                 &
       heatlat=2500800._kreal,              &
       heatfrz= 335800._kreal,              &
       heatsub=2836600._kreal
  !--------------------------------------------------------------------!
  ! parameter for turbulence                                           !
  !                                                                    !
  !  hormin   minimal horizontal turbulent energy    [(m/s)**2]        !
  !  vermin   minimal vertical turbulent energy      [(m/s)**2]        !
  !  scamin   minimal turbulent length scale         [m]               !
  !--------------------------------------------------------------------!
  real(kreal), parameter ::                 &
       hormin=2.0e-3_kreal,                 &
       vermin=1.0e-3_kreal,                 &
       scamin=1.0_kreal
  !--------------------------------------------------------------------!
  ! limits                                                             !
  !--------------------------------------------------------------------!
  real(kreal), parameter ::                 &
       epsmach=1e-7_kreal,                  &
       epsmin =1.e-20_kreal,                &
       gasmin =1e-7_kreal,                  &
       tempmin=150._kreal
  !--------------------------------------------------------------------!
  ! real constants:                                                    !
  !--------------------------------------------------------------------!
  real(kreal), parameter ::                  &
       r0=0._kreal,                          &
       r1=1._kreal,                          &
       r2=2._kreal,                          &
       r1h=0.5_kreal,                        &
       r3h=1.5_kreal,                        &
       r1q=0.25_kreal,                       &
       r3q=0.75_kreal
end module phys_constants
