module integrator_tests_precip
  use microphysics_constants, only: kreal
  use integrator_config, only: updraft_vel
  implicit none
  contains

  pure function rain_after_precip(rho_g, rho, temp, p, qr, del_t, del_z)

    use microphysics_constants, only: pi, rho_l => rho_w, g_earth

    use microphysics_common, only: dyn_visc_f => dynamic_viscosity

    real(kreal), intent(in) :: rho_g, rho, temp, p, qr, del_t, del_z
    real(kreal) :: rain_after_precip
    real(kreal), parameter :: N_0r = 1.e7_kreal  ! [m^-4]
    real(kreal), parameter :: r4_3 = 4.0_kreal / 3.0_kreal
    real(kreal), parameter :: r1_3 = 1.0_kreal / 3.0_kreal
    real(kreal), parameter :: r1_2 = 1.0_kreal / 2.0_kreal
    real(kreal), parameter :: C_Dr = 0.54_kreal
    real(kreal), parameter :: k_2 = 8.0_kreal

    real(kreal) :: lambda_r, wr1, wr2, wr3, r_r, wr

    if (qr .eq. 0.0_kreal) then
      rain_after_precip = 0.0_kreal
    else
      lambda_r = (8.0_kreal*pi*rho_l/(qr*rho)*N_0r)**0.25_kreal
      r_r = (qr*rho*lambda_r/(r4_3*pi*N_0r*rho_l))**r1_3

      wr1 = 2.0_kreal * rho_l * g_earth * r_r**2 / (9.0_kreal * dyn_visc_f(temp))
      wr2 = k_2 * rho_l * (rho/rho_g)**r1_2 * r_r
      wr3 = (8.0_kreal*rho_l*g_earth*r_r/(3.0_kreal*C_Dr*rho_g))**r1_2
      wr = min(wr1,min(wr2,wr3))

       ! if (wr > updraft_vel) then
       !   rain_after_precip = 0.0_kreal
       ! else
        rain_after_precip = qr * 1/(1+del_t*wr/del_z)
       ! endif
    endif
  end function

  pure function graupel_after_precip(rho_g, rho, temp, p, qh, del_t, del_z)

    use microphysics_constants, only: pi, rho_l => rho_w, g_earth, rho_h

    use microphysics_common, only: dyn_visc_f => dynamic_viscosity

    real(kreal), intent(in) :: rho_g, rho, temp, p, qh, del_t, del_z
    real(kreal) :: graupel_after_precip
    real(kreal), parameter :: N_0h = 1.21e4_kreal  ! [m^-4]
    real(kreal), parameter :: r4_3 = 4.0_kreal / 3.0_kreal
    real(kreal), parameter :: r1_3 = 1.0_kreal / 3.0_kreal
    real(kreal), parameter :: r1_2 = 1.0_kreal / 2.0_kreal
    real(kreal), parameter :: C_Dh = 0.6_kreal
    real(kreal), parameter :: k_2 = 8.0_kreal

    real(kreal) :: lambda_h, wh1, wh2, wh3, r_h, wh, eq_mass

    if (qh .eq. 0.0_kreal) then
      graupel_after_precip = 0.0_kreal
    else
      lambda_h = (8.0_kreal*pi*rho_h/(qh*rho)*N_0h)**0.25_kreal
      r_h = (qh*rho*lambda_h/(r4_3*pi*N_0h*rho_h))**r1_3

      wh1 = 2.0_kreal * rho_h * g_earth * r_h**2 / (9.0_kreal * dyn_visc_f(temp))
      wh2 = k_2 * rho_h * (rho/rho_g)**r1_2 * r_h
      wh3 = (8.0_kreal*rho_h*g_earth*r_h/(3.0_kreal*C_Dh*rho_g))**r1_2
      wh = min(wh1,min(wh2,wh3))

       ! if (wh > updraft_vel) then
       !   graupel_after_precip = 0.0_kreal
       ! else
        graupel_after_precip = qh * 1/(1+del_t*wh/del_z)
       ! endif
    endif
  end function

end module integrator_tests_precip
