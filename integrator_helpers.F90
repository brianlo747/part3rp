module isometric_integration_helpers
   use species_rate_of_change, only: dydt_mphys => dydt_mphys_with_ice
   !use microphysics_register, only: dydt_mphys => dydt_registered
   use integrator_rkf, only: integrate_with_message

   implicit none

   public integrate_isometric

   contains
      subroutine integrate_isometric(y, t0, t_end, msg_out, m_total, n)
         integer, intent(in) :: n
         real(8), intent(inout), dimension(n) :: y
         real(8), intent(in) :: t_end, t0
         integer, intent(inout) :: m_total
         character(len=100), intent(inout), optional :: msg_out

         call integrate_with_message(dydt_isometric, y, increment_state_isometric, t0, t_end, msg_out, m_total)
      end subroutine integrate_isometric

      pure function dydt_isometric(t, y) result(dydt)
         use microphysics_common, only: cv_mixture

         real(8), intent(in) :: t
         real(8), intent(in), dimension(:) :: y
         real(8) :: c_m, dydt(size(y))

         c_m = cv_mixture(y)

         dydt = dydt_mphys(t, y, c_m)
      end function dydt_isometric

      pure function fix_y_isometric(y_old, dy) result(y_new)
         real(8), intent(in), dimension(:) :: y_old, dy
         real(8), dimension(size(y_old)) :: y_new

         y_new = y_old + dy
      end function fix_y_isometric

      !> Return a new state changed by increment `dy` keeping constant volume,
      !method is intended to guarantee self-consistency of new state `y_new`
      ! TODO: Implement mass-scaling
      function increment_state_isometric(y, dy) result(y_new)
         !use microphysics_register, only: idx_pressure

         real(8), intent(in), dimension(:) :: y, dy
         real(8) :: y_new(size(y)), rho_old

         rho_old = density_eos(y)
         y_new = y + dy

         ! need to update pressure, since we assume constant volume density is
         ! unchanged, use old density to calculate update pressure
         y_new(2) = pressure_eos(y_new, rho_old)

      end function increment_state_isometric

      !> Calculate density of mixture from equation of state
      pure function density_eos(y) result(rho)
         !use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain
         !use microphysics_register, only: idx_temp, idx_pressure
         use microphysics_constants, only: R_v, R_d, rho_l => rho_w

         real(8), intent(in), dimension(:) :: y

         real(8) :: temp, p, qd, qv, ql, qr, qi, qh
         real(8) :: rho

         temp = y(1)
         p = y(2)

         ql = y(3)
         !if (idx_rain /= 0) then
            qr = y(4)
         !else
         !  qr = 0.
         !endif
         qv = y(5)
         qi = y(6)
         qh = y(7)
         qd = 1.0 - ql - qr - qv - qi - qh

         rho = 1.0/( (qd*R_d + qv*R_v)*temp/p + (ql+qr)/rho_l + qi/rho_i + qh/rho_h )
      end function density_eos

      !> Compute pressure from equation of state.
      ! Because state representation doesn't contain density (it is diagnosed)
      ! we need to supply a density so that pressure can be calculated from the
      ! other state variables
      pure function pressure_eos(y, rho) result(p)
         !use microphysics_register, only: idx_temp
         !use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain
         use microphysics_constants, only: R_v, R_d, rho_l => rho_w

         real(8), intent(in) :: rho
         real(8), intent(in), dimension(:) :: y

         real(8) :: temp, qd, qv, ql, qr
         real(8) :: p

         temp = y(1)

         ql = y(3)
         !if (idx_rain /= 0) then
            qr = y(4)
         !else
         !   qr = 0.
         !endif
         qv = y(5)
         qi = y(6)
         qh = y(7)
         qd = 1.0 - ql - qr - qv - qi - qh

         p = temp*(qd*R_d + qv*R_v)/(1./rho - (ql + qr)/rho_l - qi/rho_i - qh/rho_h)

      end function pressure_eos
end module

module species_rate_of_change
  use mphys_with_ice

  implicit none

  contains
    pure function dydt_mphys_with_ice(t, y, c_m)
       !use microphysics_register, only: idx_temp, idx_pressure
       !use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain, idx_cice, idx_graupel
       use microphysics_constants, only: L_v => L_cond
       use mphys_with_ice

       real(kreal), dimension(:), intent(in) :: y
       real(kreal), intent(in) :: t, c_m
       real(kreal) :: dydt(size(y))

       real(kreal) :: ql, qv, qg, qr, qd, qi, qh
       real(kreal) :: rho, rho_g
       real(kreal) :: dqrdt_autoconv, dqrdt_accre, dqrdt_condevap, dqldt_condevap
       real(kreal) :: temp, pressure

       ! OBS: it's important to make sure that return variable is initiated to
       ! zero
       dydt = 0.0_kreal

       temp = y(1)
       pressure = y(2)

       ! pick out specific concentrations from state vectors
       ql = y(3)           !idx_cwater
       qr = y(4)             !idx_rain
       qv = y(5)     !idx_water_vapour
       qi = y(6)             !idx_cice
       qh = y(7)          !idx_graupel
       qd = 1.0_kreal - ql - qr - qv - qi - qh
       qg = qv + qd

       ! compute gas and mixture density using equation of state
       rho = rho_f(qd, qv, ql, qr, qi, qh, p, temp)
       rho_g = rho_f(qd, qv, 0.0_kreal, 0.0_kreal, 0.0_kreal, 0.0_kreal, pressure, temp)

       ! compute time derivatives for each process TODO: Functions for each function

       dqldt_condevap = dql_dt__condensation_evaporation(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=temp, p=pressure)
       dqrdt_condevap = dqr_dt__condensation_evaporation(qv=qv, qr=qr, rho=rho, T=temp, p=pressure)
       dqhdt_condevap = dqh_dt__condensation_evaporation(qv, qh, rho, T, p)
       dqidt_sublidep = dqi_dt__sublimation_deposition(qi, rho, T, p)
       dqhdt_sublidep = dqh_dt__sublimation_evaporation(qg, qv, qh, rho, T, p)
       dqrdt_autoconv = dqr_dt__autoconversion(ql, qg, rho_g)
       dqhdt_autoconv = dqh_dt__autoconversion_ice_graupel(qi, qg, rho_g, T)
       dqrdt_accre_rc = dqr_dt__accretion_cloud_rain(ql, rho_g, qr)
       dqhdt_accre_hi = dqh_dt__accretion_ice_graupel(qi, rho_g, qh, T)
       dqldt_accre_chr= dqh_dt__accretion_cloud_graupel_rain(ql, rho_g, rho, qh)
       dqhdt_accre_hiri= dqr_dt__accretion_ice_rain_graupel_i(qi, rho_g, rho, qr)
       dqhdt_accre_hirr= dqr_dt__accretion_ice_rain_graupel_r(qi, ql, rho, rho_g, qr)
       dqhdt_accre_hr = dqr_dt__accretion_graupel(qg, rho_g, qv, qh, rho, T, p, qr)
       dqrdt_melt_rh  = dqr_dt__melting_graupel(qg, rho_g, qv, qh, rho, T, p, qr, ql)
       dqldt_melt_lc  = dqr_dt__melting_ice(qg, rho_g, qv, qh, rho, T, p, qr, ql)
       dqidt_freeze   = dqi_dt__freezing_graupel(qh, rho, T)
       dqhdt_freeze   = dqh_dt__freezing_ice(ql, rho, T)



       ! combine to create time derivatives for species TODO: Include many more terms!

       dydt(3) =  dqldt_condevap - dqrdt_autoconv - dqrdt_accre_rc - dqidt_freeze
       dydt(4) =                   dqrdt_autoconv + dqrdt_accre_rc + dqrdt_condevap
       dydt(5) = -dqldt_condevap                                - dqrdt_condevap
       dydt(6) =
       dydt(7) =

       dydt(1) = L_v/c_m*dqldt_condevap

    end function

end module
