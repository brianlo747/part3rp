!> Simple microphysics implementation which only supports the creation of cloud
!> water from water vapour and rain-droplets through auto-conversion and
!> accretion, no ice-phases are included.
!TODO: Graupel Density!!

module mphys_with_ice
   !use microphysics_register, only: register_variable
   use microphysics_constants, only: kreal
   !use microphysics_register, only: n_variables

   implicit none

   !public init, rho_f, dqr_dt__autoconversion

   contains
   ! subroutine init(p_disable_rain)
   !    logical, optional :: p_disable_rain
   !
   !    call register_variable('cloud_water', 1)
   !    call register_variable('rain', 1)
   !    call register_variable('cloud_ice', 1)
   !    call register_variable('idx_graupel', 1)
   !
   ! end subroutine

   ! pure function dydt(t, y, c_m)
   !    use microphysics_register, only: idx_temp, idx_pressure
   !    use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain, idx_cice, idx_graupel
   !    use microphysics_constants, only: L_v => L_cond
   !
   !    real(kreal), dimension(:), intent(in) :: y
   !    real(kreal), intent(in) :: t, c_m
   !    real(kreal) :: dydt(size(y))
   !
   !    real(kreal) :: ql, qv, qg, qr, qd, qh
   !    real(kreal) :: rho, rho_g
   !    real(kreal) :: dqrdt_autoconv, dqrdt_accre, dqrdt_condevap, dqldt_condevap
   !    real(kreal) :: temp, pressure
   !
   !    ! OBS: it's important to make sure that return variable is initiated to
   !    ! zero
   !    dydt = 0.0_kreal
   !
   !    temp = y(idx_temp)
   !    pressure = y(idx_pressure)
   !
   !    ! pick out specific concentrations from state vectors
   !    ql = y(idx_cwater)
   !    qr = y(idx_rain)
   !    qv = y(idx_water_vapour)
   !    qi = y(idx_cice)
   !    qh = y(idx_graupel)
   !    qd = 1.0_kreal - ql - qr - qv - qi - qh
   !    qg = qv + qd
   !
   !    ! compute gas and mixture density using equation of state
   !    rho = rho_f(qd, qv, ql, qr, pressure, temp)
   !    rho_g = rho_f(qd, qv, 0.0_kreal, 0.0_kreal, pressure, temp)
   !
   !    ! compute time derivatives for each process TODO: Include other processes!
   !    dqrdt_autoconv = dqr_dt__autoconversion(ql, qg, rho_g)
   !    dqrdt_accre    = dqr_dt__accretion(ql, qg, rho_g, qr)
   !    dqldt_condevap = dqr_dt__condensation_evaporation(qv=qv, qr=qr, rho=rho, T=temp, p=pressure)
   !    dqldt_condevap = dql_dt__condensation_evaporation(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=temp, p=pressure)
   !
   !    ! combine to create time derivatives for species TODO: Include many more terms!
   !    dydt(idx_water_vapour) = -dqldt_condevap                                - dqrdt_condevap
   !    dydt(idx_cwater)       =  dqldt_condevap - dqrdt_autoconv - dqrdt_accre
   !    dydt(idx_rain)         =                   dqrdt_autoconv + dqrdt_accre + dqrdt_condevap
   !
   !    dydt(idx_temp) = L_v/c_m*dqldt_condevap
   !
   ! end function

   pure function rho_f(qd, qv, ql, qr, qi, qh, p, temp) result(rho)
      use microphysics_constants, only: R_v, R_d, rho_i, rho_l => rho_w

      real(kreal), intent(in) :: qd, qv, ql, qr, qi, qh, p, temp
      real(kreal) :: rho, rho_inv
      real(kreal), parameter :: rho_h = 470.0_kreal

      rho_inv = (qd*R_d + qv*R_v)*temp/p + (ql+qr)/rho_l + qi/rho_i + qh/rho_h

      rho = 1.0_kreal/rho_inv
   end function rho_f

   !> Condesation/evaporation of cloud-water droplets
   pure function dql_dt__condensation_evaporation(rho, rho_g, qv, ql, T, p)
      use microphysics_common, only: pv_sat_f => saturation_vapour_pressure
      use microphysics_common, only: qv_sat_f => saturation_vapour_concentration
      use microphysics_common, only: Ka_f => thermal_conductivity
      use microphysics_common, only: Dv_f => water_vapour_diffusivity
      use microphysics_constants, only: Lv => L_cond, R_v
      use microphysics_constants, only: pi, rho_l => rho_w

      real(kreal), intent(in) :: rho, rho_g, qv, ql, T, p
      real(kreal) :: dql_dt__condensation_evaporation

      real(kreal), parameter :: r0 = 0.1e-6_kreal  ! initial cloud droplet radius
      real(kreal), parameter :: N0 = 200*1.0e6_kreal  ! initial cloud droplet number

      real(kreal) :: r_c, Nc
      real(kreal) :: pv_sat, qv_sat, Sw
      real(kreal) :: Dv, Fd
      real(kreal) :: Ka, Fk

      real(kreal), parameter :: r4_3 = 4.0_kreal/3.0_kreal
      real(kreal), parameter :: r1_3 = 1.0_kreal/3.0_kreal

      ! calculate saturation concentration of water vapour
      pv_sat = pv_sat_f(T)
      qv_sat = qv_sat_f(T, p)
      Sw = qv/qv_sat

      ! calculate radius of cloud droplet
      r_c = (ql*rho/(r4_3*pi*N0*rho_l))**r1_3

      ! don't allow evaporation if the droplets are smaller than the
      ! aerosol, there's nothing to evaporate then(!)
      if (Sw < 1.0_kreal) then
         if (r_c < r0) then
            r_c = 0.0_kreal
         else
         endif
      else
         r_c = max(r_c, r0)
      endif

      ! Setting cloud droplet concentraion to constant
      Nc = N0

      ! Calculation of diffusion F_d and heat F_k terms

      Ka = Ka_f(T)
      Fk = (Lv/(R_v*T) - 1._kreal) * (Lv/(Ka*T))

      Dv = Dv_f(T, p)
      Fd = R_v*T/(pv_sat*Dv)

      ! Compute rate of change of condensate from diffusion
      dql_dt__condensation_evaporation = 4.*pi*Nc/rho*r_c*(Sw - 1.0)/(Fk + Fd)

   end function dql_dt__condensation_evaporation

   !> Condensation and evaporation of rain. Similar to cloud-water droplet
   !> condensation/evaporation but includes corrections for "ventilation" and
   !> and droplet-size is assumed to follow a Marshall-Palmer distribution:
   !>
   !>     N(r)dr = N0 exp(-l*r) dr
   pure function dqr_dt__condensation_evaporation(qv, qr, rho, T, p)
      use microphysics_common, only: pv_sat_f => saturation_vapour_pressure
      use microphysics_common, only: qv_sat_f => saturation_vapour_concentration
      use microphysics_common, only: Ka_f => thermal_conductivity
      use microphysics_common, only: Dv_f => water_vapour_diffusivity
      use microphysics_common, only: dyn_visc_f => dynamic_viscosity
      use microphysics_constants, only: Lv => L_cond, R_v
      use microphysics_constants, only: pi, rho_l => rho_w

      real(kreal), intent(in) :: qv, qr, rho, T, p
      real(kreal) :: dqr_dt__condensation_evaporation

      real(kreal), parameter :: G2p75 = 1.608359421985546_kreal ! = Gamma(2.75)

      ! droplet-size distribution constant
      real(kreal), parameter :: N0 = 1.0e7_kreal ! [m^-4]

      ! fall-speed coefficient taken from the r > 0.5mm expression for
      ! fall-speed from Herzog '98
      real(kreal), parameter :: a_r = 201._kreal
      ! reference density
      real(kreal), parameter :: rho0 = 1.12_kreal

      real(kreal) :: pv_sat, qv_sat, Sw, l_r
      real(kreal) :: Dv, Fd
      real(kreal) :: Ka, Fk
      real(kreal) :: nu, f

      ! can't do cond/evap without any rain-droplets present
      if (qr == 0.0) then
         dqr_dt__condensation_evaporation = 0.0_kreal
      else
         ! calculate super/sub-saturation
         qv_sat = qv_sat_f(T, p)
         Sw = qv/qv_sat

         ! size-distribtion length-scale
         l_r = (8.0_kreal*rho_l*pi*N0/(qr*rho))**0.25_kreal

         ! air condutivity and diffusion effects
         Ka = Ka_f(T)
         Fk = (Lv/(R_v*T) - 1.0_kreal)*Lv/(Ka*T)

         pv_sat = pv_sat_f(T)
         Dv = Dv_f(T, p)
         Fd = R_v*T/(pv_sat*Dv)

         ! compute the ventilation coefficient `f`
         ! dynamic viscosity
         nu = dyn_visc_f(T=T)

         ! Calcualte factor from ventilation coefficient
         f = 1.0_kreal + 0.22_kreal*(2.0_kreal*a_r*rho/nu)**0.5_kreal * &
         (rho0/rho)**0.25_kreal*G2p75/(l_r**0.75_kreal)

         ! compute rate of change of condensate from diffusion
         dqr_dt__condensation_evaporation = 4.0_kreal*pi/rho * &
         N0/l_r**2.0_kreal*(Sw - 1.0_kreal)/(Fk + Fd)*f
      endif
   end function

   pure function dqh_dt__condensation_evaporation(qv, qh, rho, T, p)
      use microphysics_common, only: pv_sat_f => saturation_vapour_pressure
      use microphysics_common, only: qv_sat_f => saturation_vapour_concentration
      use microphysics_common, only: Ka_f => thermal_conductivity
      use microphysics_common, only: Dv_f => water_vapour_diffusivity
      use microphysics_common, only: dyn_visc_f => dynamic_viscosity
      use microphysics_constants, only: Lv => L_cond, R_v
      use microphysics_constants, only: pi, rho_l => rho_w

      real(kreal), intent(in) :: qv, qh, rho, T, p
      real(kreal) :: dqh_dt__condensation_evaporation

      real(kreal), parameter :: G2p75 = 1.608359421985546_kreal ! = Gamma(2.75)

      ! droplet-size distribution constant
      real(kreal), parameter :: N_0h = 1.21e4_kreal ! [m^-4]

      ! fall-speed coefficient taken from the r > 0.5mm expression for
      ! fall-speed from Herzog '98
      real(kreal), parameter :: a_h = 174.7_kreal
      ! reference density
      real(kreal), parameter :: rho0 = 1.12_kreal
      real(kreal), parameter :: rho_h = 470.0_kreal

      real(kreal) :: pv_sat, qv_sat, Sw, l_h
      real(kreal) :: Dv, Fd
      real(kreal) :: Ka, Fk
      real(kreal) :: nu, f

      ! can't do cond/evap without any rain-droplets present
      if (qh == 0.0) then
         dqh_dt__condensation_evaporation = 0.0_kreal
      else
         ! calculate super/sub-saturation
         qv_sat = qv_sat_f(T, p)
         Sw = qv/qv_sat

         ! size-distribtion length-scale
         l_h = (8.0_kreal*rho_h*N_0h/(qh*rho))**0.25_kreal

         ! air condutivity and diffusion effects
         Ka = Ka_f(T)
         Fk = (Lv/(R_v*T) - 1.0_kreal)*Lv/(Ka*T)

         pv_sat = pv_sat_f(T)
         Dv = Dv_f(T, p)
         Fd = R_v*T/(pv_sat*Dv)

         ! dynamic viscosity
         nu = dyn_visc_f(T=T)

         ! Calcualte factor from ventilation coefficient
         f = 1.6_kreal + 0.3_kreal*(2.0_kreal*a_h*rho/nu)**0.5_kreal * &
         (rho0/rho)**0.25_kreal*G2p75/(l_h**0.75_kreal)

         ! compute rate of change of condensate from diffusion

         dqh_dt__condensation_evaporation = 4.0_kreal*pi/rho * &
         N_0h/l_h**2.0_kreal*(Sw - 1.0_kreal)/(Fk + Fd)*f
      endif
   end function

   pure function dqi_dt__sublimation_deposition(qi, rho, T, p) ! TODO: Check sat pressure!!
     use microphysics_constants, only: pi, rho_i, T0
     use microphysics_constants, only: Ls => L_subl, R_v
     use microphysics_common, only: pv_sat_f => saturation_vapour_pressure
     use microphysics_common, only: qv_sat_f => saturation_vapour_concentration
     use microphysics_common, only: Ka_f => thermal_conductivity
     use microphysics_common, only: Dv_f => water_vapour_diffusivity


     real(kreal) :: dqi_dt__sublimation_deposition
     real(kreal), intent(in) :: qi, rho, T, p
     real(kreal), parameter :: m_i = 1.0e-12_kreal
     real(kreal), parameter :: N_0f = 1.0e-2_kreal
     real(kreal), parameter :: beta = 0.6_kreal
     !real(kreal), parameter :: T_0 = 273.15_kreal
     real(kreal), parameter :: r4_3 = 4.0_kreal/3.0_kreal
     real(kreal), parameter :: r1_3 = 1.0_kreal/3.0_kreal

     real(kreal) :: N_f, Dv, Fd_prime, Fk_prime, Ka, qi_sat, pv_sat, Si
     real(kreal) :: r_i

     if (qi == 0.0 .OR. T > T0) then
        dqi_dt__sublimation_deposition = 0.0_kreal
     else

       N_f = N_0f * EXP(beta * (T0 - T)) ! Cooper parameterisation

       r_i = (qi*rho/(r4_3*pi*N_f*rho_i))**r1_3

       ! calculate super/sub-saturation
       qi_sat = qv_sat_f(T, p)
       Si = qi/qi_sat

       ! air condutivity and diffusion effects
       Ka = Ka_f(T)
       Fk_prime = (Ls/(R_v*T) - 1.0_kreal)*Ls/(Ka*T)

       pv_sat = pv_sat_f(T)
       Dv = Dv_f(T, p)
       Fd_prime = R_v*T/(pv_sat*Dv)

       dqi_dt__sublimation_deposition = 1.6_kreal*4.0_kreal*pi*N_f/rho*r_i * &
       (Si - 1.0)/(Fk_prime + Fd_prime)

       dqi_dt__sublimation_deposition = max(0.0_kreal, &
       dqi_dt__sublimation_deposition)
     endif
   end function

   pure function dqh_dt__sublimation_evaporation(qg, qv, qh, rho, T, p) ! TODO: Check sat pressure!!
     use microphysics_common, only: pv_sat_f => saturation_vapour_pressure
     use microphysics_common, only: qv_sat_f => saturation_vapour_concentration
     use microphysics_common, only: Ka_f => thermal_conductivity
     use microphysics_common, only: Dv_f => water_vapour_diffusivity
     use microphysics_common, only: dyn_visc_f => dynamic_viscosity
     use microphysics_constants, only: Ls => L_subl, R_v
     use microphysics_constants, only: pi, rho_l => rho_w

     real(kreal), intent(in) :: qg, qv, qh, rho, T, p
     real(kreal) :: dqh_dt__sublimation_evaporation

     real(kreal), parameter :: G2p75 = 1.608359421985546_kreal ! = Gamma(2.75)

     ! droplet-size distribution constant
     real(kreal), parameter :: N_0h = 1.21e4_kreal  ! [m^-4]

     ! fall-speed coefficient taken from the r > 0.5mm expression for
     ! fall-speed from Herzog '98
     real(kreal), parameter :: a_h = 174.7_kreal  ! [m^.5 s^-1]
     ! reference density
     real(kreal), parameter :: rho0 = 1.12_kreal
     real(kreal), parameter :: rho_h = 470.0_kreal

     real(kreal) :: pv_sat, qh_sat, Sh, l_h
     real(kreal) :: Dv, Fd_prime
     real(kreal) :: Ka, Fk_prime
     real(kreal) :: nu
     real(kreal) :: f

     ! can't do cond/evap without any graupel present
     if (qh == 0.0) then
        dqh_dt__sublimation_evaporation = 0.0_kreal
     else
        ! computer super/sub-saturation
        qh_sat = qv_sat_f(T, p)
        Sh = qh/qh_sat

        ! size-distribtion length-scale
        l_h = (8.0_kreal*rho_h*N_0h/(qh*rho))**0.25_kreal

        ! air condutivity and diffusion effects
        Ka = Ka_f(T)
        Fk_prime = (Ls/(R_v*T) - 1.0_kreal)*Ls/(Ka*T)

        pv_sat = pv_sat_f(T)
        Dv = Dv_f(T, p)
        Fd_prime = R_v*T/(pv_sat*Dv)

        ! dynamic viscosity
        nu = dyn_visc_f(T=T) / rho

        ! Calcualte factor from ventilation coefficient
        f = 1.6_kreal + 0.3_kreal*(2.0_kreal*a_h*rho/nu)**0.5_kreal * &
        (rho0/rho)**0.25_kreal*G2p75/(l_h**0.75_kreal)

        dqh_dt__sublimation_evaporation = 4.0_kreal * pi / rho * N_0h * &
        (Sh - 1)/(Fk_prime + Fd_prime) / l_h**2.0_kreal * f
     endif
   end function

   pure function dqr_dt__autoconversion(ql, qg, rho_g, T)
      real(kreal), intent(in) :: ql, qg, rho_g, T
      real(kreal) :: dqr_dt__autoconversion

      real(kreal), parameter :: k_c = 1.0e-3_kreal, a_c = 5.0e-4_kreal

      ! TODO: what happens if ql < qg ?
      dqr_dt__autoconversion = k_c*(ql - qg/rho_g*a_c)
      dqr_dt__autoconversion = max(0.0_kreal, dqr_dt__autoconversion)
   end function

   pure function dqh_dt__autoconversion_ice_graupel(qi, qg, rho_g, T)
     use microphysics_constants, only: T0

     real(kreal), intent(in) :: qi, qg, rho_g, T
     real(kreal) :: dqh_dt__autoconversion_ice_graupel

     real(kreal), parameter :: a_i = 1.0e-3_kreal

     real(kreal) :: k_i

     if (qi == 0.0 .OR. T > T0) then
       dqh_dt__autoconversion_ice_graupel  = 0.0_kreal

     else
       k_i = 1.0e-3_kreal * EXP(0.025_kreal * (T - T0)) !TODO: Check T dependence

       ! TODO: what happens if ql < qg ?
       dqh_dt__autoconversion_ice_graupel = k_i*(qi - qg/rho_g*a_i)
       dqh_dt__autoconversion_ice_graupel = max(0.0_kreal, &
       dqh_dt__autoconversion_ice_graupel)
     endif
   end function

   pure function dqr_dt__accretion_cloud_rain(ql, rho_g, qr)
      use microphysics_constants, only: pi, rho_l => rho_w

      real(kreal), intent(in) :: ql, rho_g, qr
      real(kreal) :: dqr_dt__accretion_cloud_rain

      real(kreal), parameter :: G3p5 = 3.32399614155_kreal  ! = Gamma(3.5)
      real(kreal), parameter :: N_0r = 1.e7_kreal  ! [m^-4]
      real(kreal), parameter :: a_r = 201.0_kreal  ! [m^.5 s^-1]
      real(kreal), parameter :: rho0 = 1.12_kreal

      real(kreal) :: lambda_r

      ! If there is no rain available to perform accretion there is no need to calculate the accretion rate (also avoids
      ! divide-by-zero, see https://github.com/leifdenby/unified-microphysics/issues/5)

      if (qr .eq. 0.0_kreal) then
         dqr_dt__accretion_cloud_rain = 0.0_kreal
      else
         lambda_r = (pi*(rho_l)/(qr*rho_g)*N_0r)**(0.25_kreal)

         dqr_dt__accretion_cloud_rain = pi*N_0r*a_r*(rho0/rho_g)**0.5_kreal * &
         G3p5*lambda_r**(-3.5_kreal)*ql

         dqr_dt__accretion_cloud_rain = max(0.0_kreal, &
         dqr_dt__accretion_cloud_rain)
      endif
   end function

   pure function dqh_dt__accretion_ice_graupel(qi, rho_g, qh, T) !TODO Rate too large!
     ! Equation 19
     use microphysics_constants, only: pi, rho_l => rho_w

     real(kreal), intent(in) :: qi, rho_g, qh, T
     real(kreal) :: dqh_dt__accretion_ice_graupel

     real(kreal), parameter :: G3p5 = 3.32399614155_kreal  ! = Gamma(3.5)
     real(kreal), parameter :: N_0h = 1.21e4_kreal  ! [m^-4]
     real(kreal), parameter :: a_h = 174.7_kreal  ! [m^.5 s^-1]
     real(kreal), parameter :: rho0 = 1.12_kreal
     real(kreal), parameter :: rho_h = 470.0_kreal
     real(kreal), parameter :: T_0 = 273.15_kreal

     real(kreal) :: l_h, E_hi

     ! If there is no graupel available to perform accretion there is no need to calculate the accretion rate (also avoids
     ! divide-by-zero, see https://github.com/leifdenby/unified-microphysics/issues/5)

     if (qh .eq. 0.0_kreal) then
        dqh_dt__accretion_ice_graupel = 0.0_kreal
     else
        l_h = (8.0_kreal * pi * rho_h/(qh*rho_g)*N_0h)**(0.25_kreal)

        E_hi = min(1.0_kreal, EXP(0.05_kreal * (T - T_0)))

        dqh_dt__accretion_ice_graupel = pi*E_hi*N_0h*a_h*rho_h * &
        (rho0/rho_g)**0.5_kreal*G3p5/l_h**3.5_kreal * qi

        dqh_dt__accretion_ice_graupel = max(0.0_kreal, &
        dqh_dt__accretion_ice_graupel)
     endif
   end function

   pure function dqh_dt__accretion_cloud_graupel_rain(ql, rho_g, rho, qh) !TODO: Rate too big as well?
     !Equation 20
     use microphysics_constants, only: pi, rho_l => rho_w

     real(kreal), intent(in) :: ql, rho_g, rho, qh
     real(kreal) :: dqh_dt__accretion_cloud_graupel_rain

     real(kreal), parameter :: G3p5 = 3.32399614155_kreal  ! = Gamma(3.5)
     real(kreal), parameter :: N_0h = 1.21e4_kreal  ! [m^-4]
     real(kreal), parameter :: a_h = 174.7_kreal  ! [m^.5 s^-1]
     real(kreal), parameter :: rho0 = 1.12_kreal
     real(kreal), parameter :: rho_h = 470.0_kreal

     real(kreal) :: l_h

     ! If there is no graupel available to perform accretion there is no need to calculate the accretion rate (also avoids
     ! divide-by-zero, see https://github.com/leifdenby/unified-microphysics/issues/5)

     if (qh .eq. 0.0_kreal) then
        dqh_dt__accretion_cloud_graupel_rain = 0.0_kreal
     else
        l_h = (8.0_kreal*pi*rho_h*N_0h/(qh*rho))**0.25_kreal

        dqh_dt__accretion_cloud_graupel_rain = pi*N_0h*a_h*rho_h * &
        (rho0/rho_g)**0.25_kreal*G3p5*l_h**(-3.5_kreal)*ql

        dqh_dt__accretion_cloud_graupel_rain = max(0.0_kreal, &
        dqh_dt__accretion_cloud_graupel_rain)
     endif
   end function

   pure function dqr_dt__accretion_ice_rain_graupel_i(qi, rho_g, rho, qr) !TODO: Rate too big as well?
     !Equation 21
     use microphysics_constants, only: pi, rho_l => rho_w

     real(kreal), intent(in) :: qi, rho_g, rho, qr
     real(kreal) :: dqr_dt__accretion_ice_rain_graupel_i

     real(kreal), parameter :: G3p5 = 3.32399614155_kreal  ! = Gamma(3.5)
     real(kreal), parameter :: N_0r = 1.e7_kreal  ! [m^-4]
     real(kreal), parameter :: a_r = 201.0_kreal  ! [m^.5 s^-1]
     real(kreal), parameter :: rho0 = 1.12_kreal

     real(kreal) :: lambda_r

     ! If there is no rain available to perform accretion there is no need to calculate the accretion rate (also avoids
     ! divide-by-zero, see https://github.com/leifdenby/unified-microphysics/issues/5)

     if (qr .eq. 0.0_kreal) then
        dqr_dt__accretion_ice_rain_graupel_i = 0.0_kreal
     else
        lambda_r = (8.0_kreal*pi*rho_l/(qr*rho)*N_0r)**(0.25_kreal)

        dqr_dt__accretion_ice_rain_graupel_i = pi*N_0r*a_r*rho_l * &
        (rho0/rho_g)**0.5_kreal*G3p5*lambda_r**(-3.5_kreal)*qi

        dqr_dt__accretion_ice_rain_graupel_i = max(0.0_kreal, &
        dqr_dt__accretion_ice_rain_graupel_i)
     endif
   end function

   pure function dqr_dt__accretion_ice_rain_graupel_r(qi, ql, rho, rho_g, qr) !TODO: Rate too big as well, defaulting to min?
     !Equation 21
     use microphysics_constants, only: pi, rho_l => rho_w

     real(kreal), intent(in) :: qi, ql, rho, rho_g, qr
     real(kreal) :: dqr_dt__accretion_ice_rain_graupel_r

     real(kreal), parameter :: G3p5 = 3.32399614155_kreal  ! = Gamma(3.5)
     real(kreal), parameter :: N_0r = 1.e7_kreal  ! [m^-4]
     real(kreal), parameter :: Ni = 200*1.0e6_kreal
     real(kreal), parameter :: a_r = 201.0_kreal  ! [m^.5 s^-1]
     real(kreal), parameter :: rho0 = 1.12_kreal
     real(kreal), parameter :: r4_3 = 4.0_kreal/3.0_kreal
     real(kreal), parameter :: r1_3 = 1.0_kreal/3.0_kreal

     real(kreal) :: lambda_r, r_i

     ! Radius of cloud droplet
     r_i = (ql*rho/(r4_3*pi*Ni*rho_l))**r1_3

     ! If there is no rain available to perform accretion there is no need to calculate the accretion rate (also avoids
     ! divide-by-zero, see https://github.com/leifdenby/unified-microphysics/issues/5)

     if (qr .eq. 0.0_kreal) then
        dqr_dt__accretion_ice_rain_graupel_r = 0.0_kreal
     else
        lambda_r = (8.0_kreal*pi*rho_l/(qr*rho)*N_0r)**(0.25_kreal)

        dqr_dt__accretion_ice_rain_graupel_r = lambda_r/N_0r * &
        (3.0_kreal/(4.0_kreal*pi*r_i**3.0_kreal)) * &
        dqr_dt__accretion_ice_rain_graupel_i(qi, rho_g, rho, qr)
        !dqr_dt__accretion_ice_rain_graupel_r = min(1.0_kreal, &
        !dqr_dt__accretion_ice_rain_graupel_r)
     endif
   end function

   pure function dqr_dt__accretion_graupel(qg, rho_g, qv, qh, rho, T, p, qr) !TODO: Confirm no T dep?

     use microphysics_constants, only: pi, rho_l => rho_w

     real(kreal), intent(in) :: qg, rho_g, qv, qh, rho, T, p, qr
     real(kreal) :: dqr_dt__accretion_graupel
     real(kreal), parameter :: N_0r = 1.e7_kreal  ! [m^-4]
     real(kreal), parameter :: N_0h = 1.21e4_kreal  ! [m^-4]
     real(kreal), parameter :: rho0 = 1.12_kreal
     real(kreal), parameter :: rho_h = 470.0_kreal
     real(kreal), parameter :: g = 9.80665_kreal
     real(kreal), parameter :: r4_3 = 4.0_kreal / 3.0_kreal
     real(kreal), parameter :: r1_3 = 1.0_kreal / 3.0_kreal
     real(kreal), parameter :: r1_2 = 1.0_kreal / 2.0_kreal
     real(kreal), parameter :: C_Dr = 0.54_kreal
     real(kreal), parameter :: C_Dh = 0.6_kreal
     real(kreal), parameter :: k_2 = 8.0_kreal

     real(kreal) :: lambda_r, lambda_h, wr1, wr2, wr3, wh1, wh2, wh3, f, &
     r_r, r_h, wr, wh

     r_r = (qr*rho/(r4_3*pi*N_0r*rho_l))**r1_3
     r_h = (qh*rho/(r4_3*pi*N_0h*rho_h))**r1_3
     lambda_h = (8.0_kreal*rho_h*N_0h/(qh*rho))**0.25_kreal
     lambda_r = (8.0_kreal*pi*(rho_l)/(qr*rho_g)*N_0r)**(0.25_kreal)
     wr1 = 2.0_kreal * rho_l * g * r_r**2
     wr2 = k_2 * rho_l * (rho/rho_g)**r1_2 * r_r
     wr3 = (8.0_kreal*rho_g*g/(3*C_Dr*rho_g))**r1_2
     wr = min(wr1,min(wr2,wr3))
     wh1 = 2.0_kreal * rho_h * g * r_h**2
     wh2 = k_2 * rho_l * (rho/rho_g)**r1_2 * r_r
     wh3 = (8.0_kreal*rho_g*g/(3*C_Dh*rho_g))**r1_2
     wh = min(wh1,min(wh2,wh3))
     f = ((5.0_kreal/(lambda_r**6 * lambda_h)) + &
     (2.0_kreal/(lambda_r**5 * lambda_h**2)) + &
     (0.5_kreal/(lambda_r**4 * lambda_h**3)))
     dqr_dt__accretion_graupel = pi**2 * (rho_l / rho) * &
     N_0r * N_0h * ABS(wr - wh)*f
   end function

   pure function dqr_dt__melting_graupel(qg, rho_g, qv, qh, rho, T, p, qr, ql) !TODO: Check negative sign
     use microphysics_common, only: pv_sat_f => saturation_vapour_pressure
     use microphysics_common, only: qv_sat_f => saturation_vapour_concentration
     use microphysics_common, only: Ka_f => thermal_conductivity
     use microphysics_common, only: Dv_f => water_vapour_diffusivity
     use microphysics_common, only: dyn_visc_f => dynamic_viscosity
     use microphysics_constants, only: Lv => L_cond, Lf => L_fusi, R_v
     use microphysics_constants, only: cp_l, pi, rho_l => rho_w

     real(kreal), intent(in) :: qg, rho_g, qv, qh, rho, T, p, qr, ql
     real(kreal) :: dqr_dt__melting_graupel

     real(kreal), parameter :: G2p75 = 1.608359421985546_kreal ! = Gamma(2.75)

     ! droplet-size distribution constant
     real(kreal), parameter :: N_0r = 1.e7_kreal  ! [m^-4]
     real(kreal), parameter :: N_0h = 1.21e4_kreal  ! [m^-4]

     ! fall-speed coefficient taken from the r > 0.5mm expression for
     ! fall-speed from Herzog '98
     real(kreal), parameter :: a_r = 201._kreal
     real(kreal), parameter :: a_h = 174.7_kreal  ! [m^.5 s^-1]
     ! reference density
     real(kreal), parameter :: rho0 = 1.12_kreal
     real(kreal), parameter :: rho_h = 470.0_kreal
     real(kreal), parameter :: T_0 = 273.15_kreal

     real(kreal) :: pv_sat, qv_sat, Sw, l_h, l_r
     real(kreal) :: Dv, Fd
     real(kreal) :: Ka, Fk
     real(kreal) :: nu, f1, f2, f3
     !real(kreal) :: dqh_dt__accretion_cloud_graupel, dqr_dt__accretion_graupel

     ! can't do cond/evap without any rain-droplets present
     if (qr == 0.0) then
        dqr_dt__melting_graupel = 0.0_kreal
     else
        ! computer super/sub-saturation
        qv_sat = qv_sat_f(T, p)
        Sw = qv/qv_sat

        ! size-distribtion length-scale
        l_h = (8.0_kreal*rho_h*N_0h/(qh*rho))**0.25_kreal
        l_r = (8.0_kreal*pi*(rho_l)/(qr*rho_g)*N_0r)**(0.25_kreal)

        ! air condutivity and diffusion effects
        Ka = Ka_f(T)
        Fk = (Lv/(R_v*T) - 1.0_kreal)*Lv/(Ka*T)*rho_l

        pv_sat = pv_sat_f(T)
        Dv = Dv_f(T, p)
        Fd = R_v*T/(pv_sat*Dv)*rho_l

        ! compute the ventilation coefficient `f`
        ! dynamic viscosity
        nu = dyn_visc_f(T=T) / rho

        f1 = Ka * (T - T_0) + Lv * Dv * (qv - qv_sat)
        f2 = 1.6_kreal / l_h**2.0_kreal + 0.3_kreal * G2p75 * &
        (2.0_kreal * a_h / nu)**0.5_kreal * (rho0 / rho_g)**0.25_kreal / &
        l_h**2.75_kreal
        f3 = (cp_l / Lf) * (T - T_0) * &
        (dqh_dt__accretion_cloud_graupel_rain(ql, rho_g, rho, qh) + &
         dqr_dt__accretion_graupel(qg, rho_g, qv, qh, rho, T, p, qr))

        ! compute rate of change of condensate from diffusion
        dqr_dt__melting_graupel = (qg/rho_g) * N_0h * (2 * pi / Lf) &
        * f1 * f2 - f3
     endif
   end function

   pure function dqr_dt__melting_ice(qg, rho_g, qv, qh, rho, T, p, qr, ql) !TODO: Fix from graupel to ice!
     use microphysics_common, only: pv_sat_f => saturation_vapour_pressure
     use microphysics_common, only: qv_sat_f => saturation_vapour_concentration
     use microphysics_common, only: Ka_f => thermal_conductivity
     use microphysics_common, only: Dv_f => water_vapour_diffusivity
     use microphysics_common, only: dyn_visc_f => dynamic_viscosity
     use microphysics_constants, only: Lv => L_cond, Lf => L_fusi, R_v
     use microphysics_constants, only: cp_l, pi, rho_l => rho_w

     real(kreal), intent(in) :: qg, rho_g, qv, qh, rho, T, p, qr, ql
     real(kreal) :: dqr_dt__melting_ice

     real(kreal), parameter :: G2p75 = 1.608359421985546_kreal ! = Gamma(2.75)

     ! droplet-size distribution constant
     real(kreal), parameter :: N_0r = 1.e7_kreal  ! [m^-4]
     real(kreal), parameter :: Ni = 200*1.0e6_kreal

     ! fall-speed coefficient taken from the r > 0.5mm expression for
     ! fall-speed from Herzog '98
     real(kreal), parameter :: a_r = 201._kreal
     real(kreal), parameter :: a_h = 174.7_kreal  ! [m^.5 s^-1]
     ! reference density
     real(kreal), parameter :: rho0 = 1.12_kreal
     real(kreal), parameter :: rho_h = 470.0_kreal
     real(kreal), parameter :: T_0 = 273.15_kreal

     real(kreal) :: pv_sat, qv_sat, Sw, l_h, l_r
     real(kreal) :: Dv, Fd
     real(kreal) :: Ka, Fk
     real(kreal) :: nu, f1, f2, f3
     !real(kreal) :: dqh_dt__accretion_cloud_graupel, dqr_dt__accretion_graupel

     ! can't do cond/evap without any rain-droplets present
     if (qr == 0.0) then
        dqr_dt__melting_ice = 0.0_kreal
     else
        ! computer super/sub-saturation
        qv_sat = qv_sat_f(T, p)
        Sw = qv/qv_sat

        ! size-distribtion length-scale
        l_r = (8.0_kreal*pi*(rho_l)/(qr*rho_g)*N_0r)**(0.25_kreal)

        ! air condutivity and diffusion effects
        Ka = Ka_f(T)
        Fk = (Lv/(R_v*T) - 1.0_kreal)*Lv/(Ka*T)*rho_l

        pv_sat = pv_sat_f(T)
        Dv = Dv_f(T, p)
        Fd = R_v*T/(pv_sat*Dv)*rho_l

        ! compute the ventilation coefficient `f`
        ! dynamic viscosity
        nu = dyn_visc_f(T=T) / rho

        f1 = Ka * (T - T_0) + Lv * Dv * (qv - qv_sat)
        f2 = 1.0_kreal / l_r**2.0_kreal
        !f3 = (cp_l / Lf) * (T - T_0) * &
        !(dqh_dt__accretion_cloud_graupel_rain(ql, rho_g, rho, qh) + &
        ! dqr_dt__accretion_graupel(qg, rho_g, qv, qh, rho, T, p, qr))

        ! compute rate of change of condensate from diffusion
        dqr_dt__melting_ice = (qg/rho_g) * Ni * (2 * pi / Lf) &
        * f1 * f2
     endif
   end function

   pure function dqi_dt__freezing_graupel(qh, rho, T) !TODO: Check negative sign
     use microphysics_constants, only: pi, T0, rho_l => rho_w

     real(kreal), intent(in) :: qh, rho, T
     real(kreal) :: dqi_dt__freezing_graupel
     real(kreal), parameter :: rho_h = 470.0_kreal
     real(kreal), parameter :: N_0r = 1.e7_kreal  ! [m^-4]
     real(kreal), parameter :: N_0h = 1.21e4_kreal  ! [m^-4]
     real(kreal), parameter :: A_prime = 0.66_kreal  ! [m^-4]
     real(kreal), parameter :: B_prime = 100.0_kreal  ! [m^-4]

     real(kreal) :: lambda_h
     lambda_h = (8.0_kreal*pi*rho_h*N_0h/(qh*rho))**0.25_kreal

     dqi_dt__freezing_graupel = -1280.0_kreal/lambda_h**7.0_kreal * pi**2 * &
     B_prime * N_0r * (EXP(A_prime*(T0-T))-1.0_kreal) * rho_l/rho
   end function

   pure function dqh_dt__freezing_ice(ql, rho, T) !TODO: Check negative sign
     use microphysics_constants, only: pi, T0, rho_l => rho_w

     real(kreal), intent(in) :: ql, rho, T
     real(kreal) :: dqh_dt__freezing_ice
     real(kreal), parameter :: N_c = 200*1.0e6_kreal  ! [m^-4]
     real(kreal), parameter :: A_prime = 0.66_kreal  ! [m^-4]
     real(kreal), parameter :: B_prime = 100.0_kreal  ! [m^-4]
     real(kreal), parameter :: r4_3 = 4.0_kreal / 3.0_kreal
     real(kreal), parameter :: r1_3 = 1.0_kreal / 3.0_kreal

     real(kreal) :: r_c
     r_c = (ql*rho/(r4_3*pi*N_c*rho_l))**r1_3

     dqh_dt__freezing_ice = -16.0_kreal/9.0_kreal * pi**2 * r_c**6 * &
     B_prime * N_c * (EXP(A_prime*(T0-T))-1.0_kreal) * rho_l/rho
   end function
end module mphys_with_ice
