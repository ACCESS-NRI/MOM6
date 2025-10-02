!> Functions and routines involved in calculating numerical mixing
module MOM_numerical_mixing

use MOM_diag_mediator, only : diag_ctrl
use MOM_grid,          only : ocean_grid_type
use MOM_variables,     only : thermo_var_ptrs, ocean_internal_state
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

public numerical_mixing

contains

!< Calculate the suprious ``numerical'' mixing of tracer C due to advection.
subroutine numerical_mixing(G, GV, C, C_adxy, h, h_tendency, dt, C_adx, umo, C_ady, vmo, scale_constant, rho_ref, nm)

  implicit none
  type(ocean_grid_type),    intent(inout) :: G                 !< ocean grid structure
  type(verticalGrid_type),  intent(in)    :: GV                !< ocean vertical grid structure
  real,                     intent(in) :: C(:, :, :)           !< Tracer to calculate numerical mixing for
  real,                     intent(inout) :: C_adxy(:, :, :)   !< Explicit horizontal advection of tracer C
  real,                     intent(in) :: h(:, :, :)           !< Thickness
  real,                     intent(in) :: h_tendency(:, :, :)  !< Thickness tendency
  real,                     intent(in) :: dt                   !< Model timestep
  real,                     intent(inout) :: C_adx(:, :, :)    !< Explicit zonal advection of tracer C
  real,                     intent(inout) :: umo(:, :, :)      !< Total zonal mass transport
  real,                     intent(inout) :: C_ady(:, :, :)    !< Explicit meridional advection of tracer C
  real,                     intent(inout) :: vmo(:, :, :)      !< Total meridional mass transport
  real,                     intent(in) :: scale_constant       !< Scaling for tracer e.g. Specific heat capacity for T
  real,                     intent(in) :: rho_ref              !< Reference density
  real,                     intent(inout) :: nm(:, :, :)       !< temporary numerical mixing

  integer :: is, ie, js, je, nz                                !< Grid cell centre indexes
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  !< adjust for correct dimensions
  C_adxy = C_adxy / (scale_constant * rho_ref) !< units: [C]ms⁻¹
  C_adx = C_adx / (scale_constant * rho_ref)   !< units: [C]m⁻²s⁻¹
  C_ady = C_ady / (scale_constant * rho_ref)   !< units: [C]m⁻²s⁻¹
  umo = umo / rho_ref                          !< units: m³s⁻¹
  vmo = vmo / rho_ref                          !< units: m³s⁻¹

  call thickness_weighted_variance_change(C, C_adxy, h, h_tendency, dt, is, ie, js, je, nz, nm)
  call zonal_upwind_fluxes(C, C_adx, umo, G%IareaT, is, ie, js, je, nz, nm)
  call meridional_upwind_fluxes(C, C_ady, vmo, G%IareaT, is, ie, js, je, nz, nm)

end subroutine numerical_mixing

!< Subroutine to calculate the thickness weighted variance change over a timestep
subroutine thickness_weighted_variance_change(C, C_adxy, h, h_tendency, dt, is, ie, js, je, nz, nm)

  implicit none
  real, intent(in) :: C(:, :, :)               !< Tracer to calculate numerical mixing for
  real, intent(in) :: C_adxy(:, :, :)          !< Explicit horizontal advection of tracer C
  real, intent(in) :: h(:, :, :)               !< Thickness
  real, intent(in) :: h_tendency(:, :, :)      !< Thickness tendency
  real, intent(in) :: dt                       !< Model timestep
  integer, intent(in) :: is, ie, js, je, nz    !< Grid cell centre indexes
  real, intent(inout) :: nm(:, :, :)           !< Numerical mixing diagnostic to update

  integer :: i, j, k
  real :: h1, C1, hadv, Cadv

  do i = is+1, ie-2 ; do j = js+1, je-2 ; do k = 1, nz
    h1 = h(i, j, k)
    hadv = h1 + dt * h_tendency(i, j, k)
    C1 = C(i, j, k)
    Cadv = h1 * C1 + dt * C_adxy(i, j, k)
    nm(i-1, j-1, k) = (Cadv**2 / hadv - h1 * C1**2) / dt
  enddo ; enddo ; enddo

end subroutine thickness_weighted_variance_change

!< Subroutine to calculate the zonal upwind fluxes
subroutine zonal_upwind_fluxes(C, C_adx, uh, A, is, ie, js, je, nz, nm)

  implicit none
  real, intent(in) :: C(:, :, :)            !< Tracer to calculate numerical mixing for
  real, intent(in) :: C_adx(:, :, :)        !< Explicit zonal advection of tracer C
  real, intent(in) :: uh(:, :, :)           !< Zonal transport
  real, intent(in) :: A(:, :, :)            !< Area of grid cells

  integer, intent(in) :: is, ie, js, je, nz !< Grid cell centre indexes
  real, intent(inout) :: nm(:, :, :)        !< Numerical mixing diagnostic to update

  integer :: i, j, k
  real :: Cupwind(ie-2, je-3, nz)           !< Empty variable for the upwind values of C
  real :: east, west                        !< east and west locations for zonal derivative

  call zonal_upwind_values(uh, C, Cupwind)

  do i = is+1, ie-2 ; do j = js+1, je-2 ; do k = 1, nz
    east = 2 * C_adx(i, j, k) * Cupwind(i, j, k) - uh(i, j, k) * Cupwind(i, j, k)**2
    west = 2 * C_adx(i-1, j, k) * Cupwind(i-1, j, k) - uh(i-1, j, k) * Cupwind(i-1, j, k)**2
    nm(i-1, j-1, k) = nm(i-1, j-1, k) + ((east - west) / A(i, j, k))
  enddo ; enddo ; enddo

end subroutine zonal_upwind_fluxes

!< Subroutine to calculate upwind values in zonal direction
subroutine zonal_upwind_values(u, C, Cupwind, is, ie, js, je, nz)

  implicit none
  real, intent(in) :: u(:, :, :)            !< Zonal transport
  real, intent(in) :: C(:, :, :)            !< Tracer
  integer, intent(in) :: is, ie, js, je, nz !< Grid cell centre indexes
  real, intent(inout) :: Cupwind(:, :, :)   !< Zonal upwind values of C calculated using u

  integer :: i, j, k

  do i = is, ie-2 ; do j = js, je-3 ; do k = 1, nz
    if (u(i, j, k) >= 0) then
      Cupwind(i, j, k) = C(i, j, k)
    elseif (u(i, j, k) < 0) then
      Cupwind(i, j, k) = C(i+1, j, k)
    endif
  enddo ; enddo ; enddo

end subroutine zonal_upwind_values

!< Subroutine to calculate the meriodional upwind flues
subroutine meridional_upwind_fluxes(C, C_ady, vh, A, is, ie, js, je, nz, nm )

  implicit none
  real, intent(in) :: C(:, :, :)        !< Tracer to calculate numerical mixing for
  real, intent(in) :: C_ady(:, :, :)    !< Explicit zonal advection of tracer C
  real, intent(in) :: vh(:, :, :)       !< Meridional transport
  real, intent(in) :: A(:, :, :)        !< Area of grid cells

  integer, intent(in) :: is, ie, js, je, nz !< Grid cell centre indexes
  real, intent(inout) :: nm(:, :, :)

  integer :: i, j, k
  real :: Cupwind(ie-3, je-2, nz)
  real :: north, south

  call meridional_upwind_values(vh, C, Cupwind)

  do i = is+1, ie-2 ; do j = js+1, je-2 ; do k = 1, nz
    north = 2 * C_ady(i, j, k) * Cupwind(i, j, k) - vh(i, j, k) * Cupwind(i, j, k)**2
    south = 2 * C_ady(i, j-1, k) * Cupwind(i, j-1, k) - vh(i, j-1, k) * Cupwind(i, j-1, k)**2
    nm(i-1, j-1, k) = nm(i-1, j-1, k) + ((north - south) / A(i, j, k))
  enddo ; enddo ; enddo

end subroutine meridional_upwind_fluxes

subroutine meridional_upwind_values(v, C, Cupwind, is, ie, js, je, nz)

  implicit none
  real, intent(in) :: v(:, :, :)            !< Meridional transport
  real, intent(in) :: C(:, :, :)            !< Tracer
  integer, intent(in) :: is, ie, js, je, nz !< Grid cell centre indexes
  real, intent(inout) :: Cupwind(:, :, :)   !< Meridional upwind values of C calculated using v

  integer :: i, j, k

  do i = is, ie-3 ; do j = js, js-2 ; do k = 1, nz
    if (v(i, j, k) >= 0) then
      Cupwind(i, j, k) = C(i, j, k)
    elseif (v(i, j, k) < 0) then
      Cupwind(i, j, k) = C(i, j+1, k)
    endif
  enddo ; enddo ; enddo

end subroutine meridional_upwind_values

end module MOM_numerical_mixing