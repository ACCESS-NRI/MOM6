!> Functions and routines involved in calculating numerical mixing
module MOM_numerical_mixing

implicit none ; private

public numerical_mixing

!< Calculate the numerical mixing of tracer C due to advection
real function numerical_mixing(C, C_adxy, h, h_tendency, dt, C_adx, umo, C_ady, vmo, V, scale_constant, rho_ref)

  real, intent(in) :: C              !< Tracer to calculate numerical mixing for
  real, intent(in) :: C_adxy         !< Explicit horizontal advection of tracer C
  real, intent(in) :: h              !< Thickness
  real, intent(in) :: h_tendency     !< Thickness tendency
  real, intent(in) :: dt             !< Model timestep
  real, intent(in) :: C_adx          !< Explicit zonal advection of tracer C
  real, intent(in) :: umo            !< Total zonal mass transport
  real, intent(in) :: C_ady          !< Explicit meridional advection of tracer C
  real, intent(in) :: vmo            !< Total meridional mass transport
  real, intent(in) :: V              !< Volume
  real, intent(in) :: scale_constant !< Scaling for tracer e.g. Specific heat capacity for temperature
  real, intent(in) :: rho_ref        !< Reference density

  real :: nm      !< Numerical mixing due to advection of tracer C
  real :: A       !< Area of grid cells, computed from V / h
  real :: xupwind !< Upwind index in the zonal direction
  real :: yupwind !< Upwind index in the meridional direction

  nm = call thickness_weighted_variance_change(C, C_adxy, h, h_tendency, dt)
  nm = nm + call zonal_upwind_fluxes(C, C_adx, umo, A, xupwind)
  nm = nm + call meridional_upwind_fluxes(C, C_ady, vmo, A, yupwind)

end function numerical_mixing

!< Subroutine to calculate the thickness weighted variance change over a timestep
subroutine thickness_weighted_variance_change(C, C_adxy, h, h_tendency, dt)

  implicit none
  real, intent(in) :: C              !< Tracer to calculate numerical mixing for
  real, intent(in) :: C_adxy         !< Explicit horizontal advection of tracer C
  real, intent(in) :: h              !< Thickness
  real, intent(in) :: h_tendency     !< Thickness tendency

end subroutine thickness_weighted_variance_change

!< Subroutine to calculate the zonal upwind fluxes
subroutine zonal_upwind_fluxes(C, C_adx, umo, A)

  implicit none
  real, intent(in) :: C       !< Tracer to calculate numerical mixing for
  real, intent(in) :: C_adx   !< Explicit zonal advection of tracer C
  real, intent(in) :: umo     !< Total zonal mass transport
  real, intent(in) :: A       !< Area of grid cells, computed from V / h
  real, intent(in) :: xupwind !< Upwind index in the zonal direction

end subroutine zonal_upwind_fluxes

!< Subroutine to calculate the meriodional upwind flues
subroutine meridional_upwind_fluxes(C, C_ady, vmo, A)

  implicit none
  real, intent(in) :: C       !< Tracer to calculate numerical mixing for
  real, intent(in) :: C_ady   !< Explicit meridional advection of tracer C
  real, intent(in) :: vmo     !< Total meridional mass transport
  real, intent(in) :: A       !< Area of grid cells, computed from V / h
  real, intent(in) :: yupwind !< Upwind index in the meridional direction

end subroutine meridional_upwind_fluxes

!< Subroutine to caluate upind tracer values for finite difference
subroutine calculate_upwind(Cl, Cr, u)

end subroutine calculate_upwind

end module MOM_numerical_mixing