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
subroutine numerical_mixing(G, GV, Tr, h, h_tendency, dt, Idt, uhtr, vhtr, scale_constant, x_upwind, y_upwind, nm)

  implicit none
  type(ocean_grid_type),   intent(in) :: G                    !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV                   !< Ocean vertical grid structure
  type(tracer_type),       intent(in) :: Tr                   !< Tracer
  real,                    intent(in) :: h(:, :, :)           !< Thickness
  real,                    intent(in) :: h_tendency(:, :, :)  !< Thickness tendency
  real,                    intent(in) :: dt                   !< Model timestep
  real,                    intent(in) :: Idt                  !< Inverse model timestep
  real,                    intent(in) :: uhtr(:, :, :)        !< Accumulated zonal transport
  real,                    intent(in) :: vhtr(:, :, :)        !< Accumulated meridional transport
  real,                    intent(in) :: scale_constant       !< Scaling for tracer e.g. Specific heat capacity for T
  real,                 intent(inout) :: x_upwind(:, :, :)    !< Zonal upwind values for tracer
  real,                 intent(inout) :: y_upwind(:, :, :)    !< Meridional upwind values for tracer
  real,                 intent(inout) :: nm(:, :, :)          !< Numerical mixing diagnostic

  !< Local variables
  real :: Tr_adv_scale  !< Scaling required for advection terms to ensure dimensions are correct
                        ! e.g. for temperature need to divide by specific heat capacity * rho_ref
  real :: mass_transport_scale !< Scaling required for transforming accumulated fluxes into m3 s-1.
  Tr_adv_scale = scale_constant * GV%Rho0
  mass_transport_scale =(Idt * GV%H_to_RZ) / GV%Rho0

  ! call thickness_weighted_variance_change(Tr, Tr_adv_scale, h, h_tendency, dt, Idt, G, GV, nm)
  call zonal_upwind_fluxes(Tr, Tr_adv_scale, uhtr, mass_transport_scale, G, GV, x_upwind, nm)
  call meridional_upwind_fluxes(Tr, Tr_adv_scale, vhtr, mass_transport_scale, G, GV, y_upwind, nm)

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
      Cadv = h1 * C1 +  dt * (Tr%advection_xy(i, j, k)  / Tr_adv_scale)
      nm(i, j, k) = ( (Cadv**2 / hadv) - (h1 * C1**2) ) * Idt
    enddo ; enddo
  enddo

end subroutine thickness_weighted_variance_change

!< Subroutine to calculate the zonal upwind fluxes
subroutine zonal_upwind_fluxes(Tr, Tr_adv_scale, uhtr, mass_transport_scale, G, GV, x_upwind, nm)

  implicit none
  type(tracer_type),       intent(in) :: Tr                    !< Tracer
  real,                    intent(in) :: Tr_adv_scale          !< Scaling for tracer advection
  real,                    intent(in) :: uhtr(:, :, :)         !< Accumulated zonal transport
  real,                    intent(in) :: mass_transport_scale  !< Inverse model timestep
  type(ocean_grid_type),   intent(in) :: G                     !< Ocean grid structure for inverse area
  type(verticalGrid_type), intent(in) :: GV                    !< Ocean vertical grid structure
  real,                 intent(inout) :: x_upwind(:, :, :)     !< Zonal upwind value for tracer
  real,                 intent(inout) :: nm(:, :, :)           !< Numerical mixing diagnostic to update

  !< Local variables
  integer :: is, ie, js, je, nz                          !< Grid cell centre and layer indexes
  integer :: i, j, k                                     !< Counters
  real :: east, west                                     !< East and West positions for zonal derivative
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: u_trans  !< Zonal transport

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  
  u_trans(:, :, :) = 0.
  do k=1,nz ; do j=js,je ; do I=is-1,ie
    u_trans(I,j,k) = uhtr(I,j,k) * mass_transport_scale
  enddo ; enddo ; enddo

  call zonal_upwind_values(Tr, G, nz, u_trans, x_upwind)

  do k =1, nz
    do j = js, je ; do i = is, ie
      ! east = 2 * (Tr%ad_x(I, j, k)   / Tr_adv_scale) * x_upwind(I, j, k)   - u_trans(I, j, k)   * x_upwind(I, j, k)**2
      ! west = 2 * (Tr%ad_x(I-1, j, k) / Tr_adv_scale) * x_upwind(I-1, j, k) - u_trans(I-1, j, k) * x_upwind(I-1, j, k)**2
      east = 2 * (Tr%ad_x(I, j, k)   / Tr_adv_scale) - u_trans(I, j, k)
      west = 2 * (Tr%ad_x(I-1, j, k) / Tr_adv_scale) - u_trans(I-1, j, k)
      nm(i, j, k) = nm(i, j, k) + ((east - west) * G%IareaT(i, j))
    enddo ; enddo
  enddo

end subroutine zonal_upwind_fluxes

!< Subroutine to calculate upwind values in zonal direction
subroutine zonal_upwind_values(Tr, G, nz, u_trans, x_upwind)

  implicit none
  type(tracer_type),        intent(in) :: Tr                 !< Tracer
  type(ocean_grid_type),    intent(in) :: G                  !< Ocean grid structure for inverse area
  integer,                  intent(in) :: nz                 !< Grid cell layer indexes
  real,                     intent(in) :: u_trans(:, :, :)   !< Zonal transport
  real,                  intent(inout) :: x_upwind(:, :, :)  !< Zonal upwind values of C calculated using u_trans

  !< Local variables
  integer :: is, ie, js, je  !< Grid cell centre indexes
  integer :: i, j, k         !< Counters

  Is = G%IscB ; Ie = G%IecB ; js = G%jsc ; je = G%jec

  do k = 1, nz
    do j = js, je ; do I = Is, Ie
      if (u_trans(I, j, k) >= 0) then
        x_upwind(I, j, k) = Tr%t(i, j, k)
      elseif (u_trans(I, j, k) < 0) then
        x_upwind(I, j, k) = Tr%t(i+1, j, k)
      endif
    enddo ; enddo
  enddo

end subroutine zonal_upwind_values

!< Subroutine to calculate the meriodional upwind flues
subroutine meridional_upwind_fluxes(Tr, Tr_adv_scale, vhtr, mass_transport_scale, G, GV, y_upwind, nm)

  implicit none
  type(tracer_type),       intent(in) :: Tr                    !< Tracer
  real,                    intent(in) :: Tr_adv_scale          !< Scaling for tracer advection
  real,                    intent(in) :: vhtr(:, :, :)         !< Accumulated meridional transport
  real,                    intent(in) :: mass_transport_scale  !< Scaling for mass transport
  type(ocean_grid_type),   intent(in) :: G                     !< Ocean grid structure for inverse area
  type(verticalGrid_type), intent(in) :: GV                    !< Ocean vertical grid structure
  real,                 intent(inout) :: y_upwind(:, :, :)     !< Meridional upwind tracer values
  real,                 intent(inout) :: nm(:, :, :)           !< Numerical mixing diagnostic to update

  !< Local variables
  integer :: is, ie, js, je, nz                          !< Grid cell centre and layer indexes
  integer :: i, j, k                                     !< Counters
  real :: north, south                                   !< North and South positions for meridional derivative
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: v_trans  !< Meridional transport
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  v_trans(:, :, :) = 0.
  do k=1,nz ; do J=js-1,je ; do i=is,ie
    v_trans(i,J,k) = vhtr(i,J,k) * mass_transport_scale
  enddo ; enddo ; enddo

  call meridional_upwind_values(Tr, G, nz, v_trans, y_upwind)


  do k = 1, nz
    do j = js, je ; do i = is, ie
      ! north = 2 * (Tr%ad_y(i, J, k)   / Tr_adv_scale) * y_upwind(i, J, k)   - v_trans(i, J, k)   * y_upwind(i, J, k)**2
      ! south = 2 * (Tr%ad_y(i, J-1, k) / Tr_adv_scale) * y_upwind(i, J-1, k) - v_trans(i, J-1, k) * y_upwind(i, J-1, k)**2
      north = 2 * (Tr%ad_y(i, J, k)   / Tr_adv_scale) - v_trans(i, J, k)  
      south = 2 * (Tr%ad_y(i, J-1, k) / Tr_adv_scale) - v_trans(i, J-1, k)
      nm(i, j, k) = nm(i, j, k) + ((north - south) * G%IareaT(i, j))
    enddo ; enddo
  enddo

end subroutine meridional_upwind_fluxes

subroutine meridional_upwind_values(Tr, G, nz, v_trans, y_upwind)

  implicit none
  type(tracer_type),     intent(in) :: Tr                 !< Tracer
  type(ocean_grid_type), intent(in) :: G                  !< Ocean grid structure for inverse area
  integer,               intent(in) :: nz                 !< Grid cell layer indexes
  real,                  intent(in) :: v_trans(:, :, :)   !< Meridional transport
  real,               intent(inout) :: y_upwind(:, :, :)  !< Meridional upwind values of C calculated using v_trans

  !< Local variables
  integer :: is, ie, js, je  !< Grid cell centre indexes
  integer :: i, j, k         !< Counters

  is = G%isc ; ie = G%iec ; Js = G%JscB ; Je = G%JecB

  do k = 1, nz
    do J = Js, Je ; do i = is, ie
      if (v_trans(i, J, k) >= 0) then
        y_upwind(i, J, k) = Tr%t(i, j, k)
      elseif (v_trans(i, J, k) < 0) then
        y_upwind(i, J, k) = Tr%t(i, j+1, k)
      endif
    enddo ; enddo
  enddo

end subroutine meridional_upwind_values

end module MOM_numerical_mixing
