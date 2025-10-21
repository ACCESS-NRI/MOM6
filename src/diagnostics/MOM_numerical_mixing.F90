!> Functions and routines involved in calculating numerical mixing of tracers due to advection
module MOM_numerical_mixing

use MOM_diag_mediator, only : diag_ctrl
use MOM_grid,          only : ocean_grid_type
use MOM_tracer_types,  only : tracer_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public numerical_mixing

contains

!< Calculate the suprious ``numerical'' mixing of tracer C due to advection.
subroutine numerical_mixing(G, GV, Tr, h, h_tendency, dt, Idt, uhtr, vhtr, scale_constant, nm)

  implicit none
  type(ocean_grid_type),   intent(in) :: G                    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV                   !< Ocean vertical grid structure
  type(tracer_type),       intent(in) :: Tr                   !< Tracer
  real,                    intent(in) :: h(:, :, :)           !< Thickness
  real,                    intent(in) :: h_tendency(:, :, :)  !< Thickness tendency
  real,                    intent(in) :: dt                   !< Model timestep
  real,                    intent(in) :: Idt                  !< Inverse model timestep
  real,                    intent(in) :: uhtr(:, :, :)        !< Total zonal mass transport
  real,                    intent(in) :: vhtr(:, :, :)        !< Total meridional mass transport
  real,                    intent(in) :: scale_constant       !< Scaling for tracer e.g. Specific heat capacity for T
  real,                 intent(inout) :: nm(:, :, :)          !< Numerical mixing diagnostic

  !< Local variables
  real :: Tr_adv_scale  !< Scaling required for advection terms to ensure dimensions are correct
                        ! e.g. for temperature need to divide by specific heat capacity * rho_ref
  Tr_adv_scale = scale_constant * 1035


  call thickness_weighted_variance_change(Tr, Tr_adv_scale, h, h_tendency, dt, Idt, G, GV, nm)
!  call zonal_upwind_fluxes(Tr, Tr_adv_scale, uhtr, Idt, G, GV, nm)
!  call meridional_upwind_fluxes(Tr, Tr_adv_scale, vhtr, Idt, G, GV, nm)

end subroutine numerical_mixing

!< Subroutine to calculate the thickness weighted variance change over a timestep
subroutine thickness_weighted_variance_change(Tr, Tr_adv_scale, h, h_tendency, dt, Idt, G, GV, nm)

  implicit none
  type(tracer_type),       intent(in) :: Tr                   !< Tracer
  real,                    intent(in) :: Tr_adv_scale         !< Scaling for tracer advection
  real,                    intent(in) :: h(:, :, :)           !< Thickness
  real,                    intent(in) :: h_tendency(:, :, :)  !< Thickness tendency
  real,                    intent(in) :: dt                   !< Model timestep
  real,                    intent(in) :: Idt                  !< Inverse model timestep
  type(ocean_grid_type),   intent(in) :: G                    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV                   !< Ocean vertical grid structure
  real,                 intent(inout) :: nm(:, :, :)          !< Numerical mixing diagnostic to update

  !< Local variables
  integer :: is, ie, js, je, nz  !< Grid cell centre and layer indexes
  integer :: i, j, k             !< Counters
  real :: h1, C1, hadv, Cadv     !< Temporary variables for thickness and tracer at current timestep
                                 !< and the changes in thickness and tracer due to advection.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do k = 1, nz
    do j = js, je ; do i = is, ie
      h1 = h(i, j, k)
      hadv = h1 + dt * h_tendency(i, j, k)
      C1 = Tr%t(i, j, k)
      Cadv = h1 * C1 + dt * Tr%advection_xy(i, j, k) / Tr_adv_scale
      nm(i, j, k) = (Cadv**2 / hadv - h1 * C1**2) * Idt
    enddo ; enddo
  enddo

end subroutine thickness_weighted_variance_change

!< Subroutine to calculate the zonal upwind fluxes
subroutine zonal_upwind_fluxes(Tr, Tr_adv_scale, uhtr, Idt, G, GV, nm)

  implicit none
  type(tracer_type),       intent(in) :: Tr             !< Tracer
  real,                    intent(in) :: Tr_adv_scale   !< Scaling for tracer advection
  real,                    intent(in) :: uhtr(:, :, :)  !< Accumulated zonal fluxes
  real,                    intent(in) :: Idt            !< Inverse model timestep
  type(ocean_grid_type),   intent(in) :: G              !< Ocean grid structure for inverse area
  type(verticalGrid_type), intent(in) :: GV             !< Ocean vertical grid structure
  real,                 intent(inout) :: nm(:, :, :)    !< Numerical mixing diagnostic to update

  !< Local variables
  integer :: is, ie, js, je, nz                          !< Grid cell centre and layer indexes
  integer :: i, j, k                                     !< Counters
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: Cupwind  !< Empty variable for the upwind values of C
  real :: east, west                                     !< East and West positions for zonal derivative

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call zonal_upwind_values(uhtr, Tr, G, nz, Cupwind)

  do k =1, nz
    do j = js, je ; do i = is, ie
      east = 2 * (Tr%ad_x(I, j, k)   / Tr_adv_scale) * Cupwind(I, j, k)   - (uhtr(I, j, k)   * Idt) * Cupwind(I, j, k)**2
      west = 2 * (Tr%ad_x(I-1, j, k) / Tr_adv_scale) * Cupwind(I-1, j, k) - (uhtr(I-1, j, k) * Idt) * Cupwind(I-1, j, k)**2
      nm(i, j, k) = nm(i, j, k) + ((east - west) * G%IareaT(i, j))
    enddo ; enddo
  enddo

end subroutine zonal_upwind_fluxes

!< Subroutine to calculate upwind values in zonal direction
subroutine zonal_upwind_values(uhtr, Tr, G, nz, Cupwind)

  implicit none
  real,                     intent(in) :: uhtr(:, :, :)     !< Accumulates zonal transport
  type(tracer_type),        intent(in) :: Tr                !< Tracer
  type(ocean_grid_type),    intent(in) :: G                 !< Ocean grid structure for inverse area
  integer,                  intent(in) :: nz                !< Grid cell layer indexes
  real,                  intent(inout) :: Cupwind(:, :, :)  !< Zonal upwind values of C calculated using uhtr

  !< Local variables
  integer :: is, ie, js, je  !< Grid cell centre indexes
  integer :: i, j, k         !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do k = 1, nz
    do j = js, je ; do I = is-1, ie
      if (uhtr(I, j, k) >= 0) then
        Cupwind(I, j, k) = Tr%t(i, j, k)
      elseif (uhtr(I, j, k) < 0) then
        Cupwind(I, j, k) = Tr%t(i+1, j, k)
      endif
    enddo ; enddo
  enddo

end subroutine zonal_upwind_values

!< Subroutine to calculate the meriodional upwind flues
subroutine meridional_upwind_fluxes(Tr, Tr_adv_scale, vhtr, Idt, G, GV, nm)

  implicit none
  type(tracer_type),       intent(in) :: Tr             !< Tracer
  real,                    intent(in) :: Tr_adv_scale   !< Scaling for tracer advection
  real,                    intent(in) :: vhtr(:, :, :)  !< Meridional mass transport
  real,                    intent(in) :: Idt            !< Inverse model timestep
  type(ocean_grid_type),   intent(in) :: G              !< Ocean grid structure for inverse area
  type(verticalGrid_type), intent(in) :: GV             !< Ocean vertical grid structure
  real,                 intent(inout) :: nm(:, :, :)    !< Numerical mixing diagnostic to update

  !< Local variables
  integer :: is, ie, js, je, nz                          !< Grid cell centre and layer indexes
  integer :: i, j, k                                     !< Counters
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: Cupwind  !< Empty variable for the meridional upwind tracer values
  real :: north, south                                   !< North and South positions for meridional derivative

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call meridional_upwind_values(vhtr, Tr, G, nz, Cupwind)

  do k = 1, nz
    do j = js, je ; do i = is, ie
      north = 2 * (Tr%ad_y(i, J, k)   / Tr_adv_scale) * Cupwind(i, J, k)   - (vhtr(i, J, k)   * Idt) * Cupwind(i, J, k)**2
      south = 2 * (Tr%ad_y(i, J-1, k) / Tr_adv_scale) * Cupwind(i, J-1, k) - (vhtr(i, J-1, k) * Idt) * Cupwind(i, J-1, k)**2
      nm(i, j, k) = nm(i, j, k) + ((north - south) * G%IareaT(i, j))
    enddo ; enddo
  enddo

end subroutine meridional_upwind_fluxes

subroutine meridional_upwind_values(vhtr, Tr, G, nz, Cupwind)

  implicit none
  real,                  intent(in) :: vhtr(:, :, :)     !< Accumulated meridional transport
  type(tracer_type),     intent(in) :: Tr                !< Tracer
  type(ocean_grid_type), intent(in) :: G                 !< Ocean grid structure for inverse area
  integer,               intent(in) :: nz                !< Grid cell layer indexes
  real,               intent(inout) :: Cupwind(:, :, :)  !< Meridional upwind values of C calculated using vhtr

  !< Local variables
  integer :: is, ie, js, je  !< Grid cell centre indexes
  integer :: i, j, k         !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do k = 1, nz
    do j = js, je+1 ; do i = is, ie
      if (vhtr(i, J-1, k) >= 0) then
        Cupwind(i, J-1, k) = Tr%t(i, j-1, k)
      elseif (vhtr(i, J-1, k) < 0) then
        Cupwind(i, J-1, k) = Tr%t(i, j, k)
      endif
    enddo ; enddo
  enddo

end subroutine meridional_upwind_values

end module MOM_numerical_mixing
