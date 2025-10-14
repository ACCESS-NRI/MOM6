!> Functions and routines involved in calculating numerical mixing of tracers due to advection
module MOM_numerical_mixing

use MOM_diag_mediator, only : diag_ctrl
use MOM_grid,          only : ocean_grid_type
use MOM_tracer_types,  only : tracer_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

public numerical_mixing

contains

!< Calculate the suprious ``numerical'' mixing of tracer C due to advection.
subroutine numerical_mixing(G, GV, Tr, h, h_tendency, dt, umo, vmo, scale_constant, rho_ref, nm)

  implicit none
  type(ocean_grid_type),    intent(in) :: G                    !< Ocean grid structure
  type(verticalGrid_type),  intent(in) :: GV                   !< Ocean vertical grid structure
  type(tracer_type),        intent(in) :: Tr                   !< Tracer 
  real,                     intent(in) :: h(:, :, :)           !< Thickness
  real,                     intent(in) :: h_tendency(:, :, :)  !< Thickness tendency
  real,                     intent(in) :: dt                   !< Model timestep
  real,                     intent(inout) :: umo(:, :, :)      !< Total zonal mass transport
  real,                     intent(inout) :: vmo(:, :, :)      !< Total meridional mass transport
  real,                     intent(in) :: scale_constant       !< Scaling for tracer e.g. Specific heat capacity for T
  real,                     intent(in) :: rho_ref              !< Reference density
  real,                     intent(inout) :: nm(:, :, :)       !< Numerical mixing diagnostic

  integer :: nz         !< Grid cell layer indexes
  real :: Tr_adv_scale  !< Scaling required for advection terms to ensure dimensions are correct
                        ! e.g. for temperature need to divide by specific heat capacity * rho_ref
  nz = GV%ke
  Tr_adv_scale = scale_constant * rho_ref

  !< adjust for correct dimensions
  umo = umo / rho_ref  !< units: m³s⁻¹
  vmo = vmo / rho_ref  !< units: m³s⁻¹

  call thickness_weighted_variance_change(Tr, Tr_adv_scale, h, h_tendency, dt, G, nz, nm)
  call zonal_upwind_fluxes(Tr, Tr_adv_scale, umo, G, nz, nm)
  call meridional_upwind_fluxes(Tr, Tr_adv_scale, vmo, G, nz, nm)

end subroutine numerical_mixing

!< Subroutine to calculate the thickness weighted variance change over a timestep
subroutine thickness_weighted_variance_change(Tr, Tr_adv_scale, h, h_tendency, dt, G, nz, nm)

  implicit none
  type(tracer_type),     intent(in) :: Tr                   !< Tracer
  real,                  intent(in) :: Tr_adv_scale         !< Scaling for tracer advection
  real,                  intent(in) :: h(:, :, :)           !< Thickness
  real,                  intent(in) :: h_tendency(:, :, :)  !< Thickness tendency
  real,                  intent(in) :: dt                   !< Model timestep
  type(ocean_grid_type), intent(in) :: G                    !< Ocean grid structure
  integer,               intent(in) :: nz                   !< Grid cell layer indexes
  real,                  intent(inout) :: nm(:, :, :)       !< Numerical mixing diagnostic to update

  integer :: is, ie, js, je   !< Grid cell centre indexes
  integer :: i, j, k          !< Counters
  real :: h1, C1, hadv, Cadv  !< Temporary grid cell variables

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do i = is, ie ; do j = js, je ; do k = 1, nz
    h1 = h(i, j, k)
    hadv = h1 + dt * h_tendency(i, j, k)
    C1 = Tr%t(i, j, k)
    Cadv = h1 * C1 + dt * Tr%advection_xy(i, j, k) / Tr_adv_scale
    nm(i, j, k) = (Cadv**2 / hadv - h1 * C1**2) / dt
  enddo ; enddo ; enddo

end subroutine thickness_weighted_variance_change

!< Subroutine to calculate the zonal upwind fluxes
subroutine zonal_upwind_fluxes(Tr, Tr_adv_scale, umo, G, nz, nm)

  implicit none
  type(tracer_type),     intent(in) :: Tr              !< Tracer
  real,                  intent(in) :: Tr_adv_scale    !< Scaling for tracer advection
  real,                  intent(in) :: umo(:, :, :)    !< Zonal mass transport
  type(ocean_grid_type), intent(in) :: G               !< Ocean grid structure for inverse area
  integer,               intent(in) :: nz              !< Grid cell layer indexes
  real,                  intent(inout) :: nm(:, :, :)  !< Numerical mixing diagnostic to update

  integer :: is, ie, js, je         !< Grid cell centre indexes
  integer :: i, j, k                !< Counters
  real :: Cupwind(G%iec, G%jec, nz) !< Empty variable for the upwind values of C
  real :: east, west                !< East and West positions for zonal derivative

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  call zonal_upwind_values(umo, Tr, Cupwind, G, nz)

  do i = is, ie-1 ; do j = js, je ; do k = 1, nz
    east = 2 * (Tr%ad_x(i, j, k) / Tr_adv_scale) * Cupwind(i, j, k) - umo(i, j, k) * Cupwind(i, j, k)**2
    west = 2 * (Tr%ad_x(i-1, j, k) / Tr_adv_scale) * Cupwind(i-1, j, k) - umo(i-1, j, k) * Cupwind(i-1, j, k)**2
    nm(i, j, k) = nm(i, j, k) + ((east - west) * G%IareaT(i, j)) 
  enddo ; enddo ; enddo

end subroutine zonal_upwind_fluxes

!< Subroutine to calculate upwind values in zonal direction
subroutine zonal_upwind_values(u, Tr, Cupwind, G, nz)

  implicit none
  real,                  intent(in) :: u(:, :, :)           !< Zonal transport
  type(tracer_type),     intent(in) :: Tr                   !< Tracer
  real,                  intent(inout) :: Cupwind(:, :, :)  !< Zonal upwind values of C calculated using v
  type(ocean_grid_type), intent(in) :: G                    !< Ocean grid structure for inverse area
  integer,               intent(in) :: nz                   !< Grid cell layer indexes

  integer :: is, ie, js, je  !< Grid cell centre indexes
  integer :: i, j, k         !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do i = is, ie ; do j = js, je ; do k = 1, nz
    if (u(i, j, k) >= 0) then
      Cupwind(i, j, k) = Tr%t(i, j, k)
    elseif (u(i, j, k) < 0) then
      Cupwind(i, j, k) = Tr%t(i+1, j, k)
    endif
  enddo ; enddo ; enddo

end subroutine zonal_upwind_values

!< Subroutine to calculate the meriodional upwind flues
subroutine meridional_upwind_fluxes(Tr, Tr_adv_scale, vmo, G, nz, nm)

  implicit none
  type(tracer_type),     intent(in) :: Tr            !< Tracer
  real,                  intent(in) :: Tr_adv_scale  !< Scaling for tracer advection
  real,                  intent(in) :: vmo(:, :, :)  !< Meridional mass transport
  type(ocean_grid_type), intent(in) :: G             !< Ocean grid structure for inverse area
  integer, intent(in) :: nz                          !< Grid cell layer indexes
  real, intent(inout) :: nm(:, :, :)                 !< Numerical mixing diagnostic to update

  integer :: is, ie, js, je          !< Grid cell centre indexes
  integer :: i, j, k                 !< Counters
  real :: Cupwind(G%iec, G%jec, nz)  !< Empty variable for the meridional upwind tracer values
  real :: north, south               !< North and South positions for meridional derivative
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  call meridional_upwind_values(vmo, Tr, Cupwind, G, nz)

  do i = is, ie ; do j = js, je ; do k = 1, nz
    north = 2 * (Tr%ad_y(i, j, k) / Tr_adv_scale)* Cupwind(i, j, k) - vmo(i, j, k) * Cupwind(i, j, k)**2
    south = 2 * (Tr%ad_y(i, j-1, k) / Tr_adv_scale)* Cupwind(i, j-1, k) - vmo(i, j-1, k) * Cupwind(i, j-1, k)**2
    nm(i, j, k) = nm(i, j, k) + ((north - south) * G%IareaT(i, j))
  enddo ; enddo ; enddo

end subroutine meridional_upwind_fluxes

subroutine meridional_upwind_values(v, Tr, Cupwind, G, nz)

  implicit none
  real,                  intent(in) :: v(:, :, :)           !< Meridional transport
  type(tracer_type),     intent(in) :: Tr                   !< Tracer
  real,                  intent(inout) :: Cupwind(:, :, :)  !< Meridional upwind values of C calculated using v
  type(ocean_grid_type), intent(in) :: G                    !< Ocean grid structure for inverse area
  integer,               intent(in) :: nz                   !< Grid cell layer indexes

  integer :: is, ie, js, je  !< Grid cell centre indexes
  integer :: i, j, k         !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do i = is, ie ; do j = js, js ; do k = 1, nz
    if (v(i, j, k) >= 0) then
      Cupwind(i, j, k) = Tr%t(i, j, k)
    elseif (v(i, j, k) < 0) then
      Cupwind(i, j, k) = Tr%t(i, j+1, k)
    endif
  enddo ; enddo ; enddo

end subroutine meridional_upwind_values

end module MOM_numerical_mixing
