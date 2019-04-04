module integrator_helpers
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
         use integrator_species_rate, only: dydt_mphys_with_ice

         real(8), intent(in) :: t
         real(8), intent(in), dimension(:) :: y
         real(8) :: c_m, dydt(size(y))

         c_m = cv_mixture(y)

         dydt = dydt_mphys_with_ice(t, y, c_m)
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
         use microphysics_constants, only: R_v, R_d, rho_h, rho_i, rho_l => rho_w

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
         use microphysics_constants, only: R_v, R_d, rho_h, rho_i, rho_l => rho_w

         real(8), intent(in) :: rho
         real(8), intent(in), dimension(:) :: y

         real(8) :: temp, qd, qv, ql, qr, qi, qh
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
