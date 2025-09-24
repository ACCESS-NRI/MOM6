!> Functions and routines involved in calculating numerical mixing
module MOM_numerical_mixing

use MOM_grid,              only : ocean_grid_type
use MOM_variables,         only : thermo_var_ptrs, ocean_internal_state
use MOM_verticalGrid,      only : verticalGrid_type
implicit none ; private

public numerical_mixing

contains

!< Calculate numerical mixing from saved output.
subroutine numerical_mixing(MIS, tv, G, GV, nm)
  type(ocean_internal_state), intent(in) :: MIS          !< For "MOM Internal State" a set of pointers to
                                                         !! the fields and accelerations making up ocean
                                                         !! internal physical state.
  type(thermo_var_ptrs),      intent(in) :: tv           !< A structure pointing to various
                                                         !! thermodynamic variables.
  type(ocean_grid_type),      intent(in) :: G            !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV           !< The ocean's vertical grid structure.
  real, intent(inout)                    :: nm(:, :, :)  !< numerical mxiing output

  !< adjust for correct dimensions
  C_adxy = C_adxy / (tv%C_p * rho_ref)  !< units: [C]ms⁻¹
  C_adx = C_adx / (tv%C_p  * rho_ref)   !< units: [C]m⁻²s⁻¹
  C_ady = C_ady / (tv%C_p  * rho_ref)   !< units: [C]m⁻²s⁻¹

  call thickness_weighted_variance_change(MIS%T, C_adxy, MIS%h, h_tendency, dt, G%iec, G%jec, GV%ke, nm)
  call zonal_upwind_fluxes(MIS%T, C_adx, MIS%uh, G%IareaT, G%iec, G%jec, GV%ke, nm)
  call meridional_upwind_fluxes(MIS%T, C_ady, MIS%vh, G%IareaT, G%iec, G%jec, GV%ke, nm)

end subroutine numerical_mixing

!< Subroutine to calculate the thickness weighted variance change over a timestep
subroutine thickness_weighted_variance_change(C, C_adxy, h, h_tendency, dt, Nx, Ny, Nz, nm)

  implicit none
  real, intent(in) :: C(:, :, :, :)               !< Tracer to calculate numerical mixing for
  real, intent(in) :: C_adxy(:, :, :, :)          !< Explicit horizontal advection of tracer C
  real, intent(in) :: h(:, :, :, :)               !< Thickness
  real, intent(in) :: h_tendency(:, :, :, :)      !< Thickness tendency
  real, intent(in) :: dt                          !< Model timestep
  integer, intent(in) :: Nx, Ny, Nz               !< Dimension lengths
  real, intent(inout) :: nm(:, :, :)

  integer :: i, j, k
  real :: h1, C1, hadv, Cadv

  do i = 2, Nx-2
    do j = 2, Ny-2
      do k = 1, Nz
        h1 = h(i, j, k, 1)
        hadv = h1 + dt * h_tendency(i, j, k, 1)
        C1 = C(i, j, k, 1)
        Cadv = (h1 * C1 + dt * C_adxy(i, j, k, 1)) / hadv
        nm(i-1, j-1, k) = (hadv * Cadv**2 - h1 * C1**2) / dt
      enddo
    enddo
  enddo

end subroutine thickness_weighted_variance_change

!< Subroutine to calculate the zonal upwind fluxes
subroutine zonal_upwind_fluxes(C, C_adx, uh, V, Nx, Ny, Nz, nm)

  implicit none
  real, intent(in) :: C(:, :, :, :)        !< Tracer to calculate numerical mixing for
  real, intent(in) :: C_adx(:, :, :, :)    !< Explicit zonal advection of tracer C
  real, intent(in) :: uh(:, :, :, :)      !< Zonal mass transport
  real, intent(in) :: V(:, :, :, :)        !< Area of grid cells, computed from V / h

  integer, intent(in) :: Nx, Ny, Nz        !< Dimension lengths
  real, intent(inout) :: nm(:, :, :)

  integer :: i, j, k
  real :: Cupwind(Nx-2, Ny-3, Nz)
  real :: east, west

  call zonal_upwind_values(uh, C, Cupwind)

  do i = 2, Nx-2
    do j = 2, Ny-2
      do k = 1, Nz
        east = 2 * C_adx(i, j, k, 1) * Cupwind(i, j, k) - uh(i, j, k, 1) * Cupwind(i, j, k)**2
        west = 2 * C_adx(i-1, j, k, 1) * Cupwind(i-1, j, k) - uh(i-1, j, k, 1) * Cupwind(i-1, j, k)**2
        nm(i-1, j-1, k) = nm(i-1, j-1, k) + ((east - west) / V(i, j, k, 1))
      enddo
    enddo
  enddo

end subroutine zonal_upwind_fluxes

!< Subroutine to calculate upwind values in zonal direction
subroutine zonal_upwind_values(u, C, Cupwind)

  implicit none
  real, intent(in) :: u(:, :, :, :)
  real, intent(in) :: C(:, :, :, :)
  real, intent(inout) :: Cupwind(:, :, :)

  integer :: i, j, k
  integer :: xdim, ydim, zdim

  xdim = size(Cupwind, dim = 1)
  ydim = size(Cupwind, dim = 2)
  zdim = size(Cupwind, dim = 3)

  do i = 1, xdim
    do j = 1, ydim
      do k = 1, zdim
        if (u(i, j, k, 1) >= 0) then
          Cupwind(i, j, k) = C(i, j, k, 1)
        else if (u(i, j, k, 1) < 0) then
          Cupwind(i, j, k) = C(i+1, j, k, 1)
        end if
      end do
    end do
  end do

end subroutine zonal_upwind_values

!< Subroutine to calculate the meriodional upwind flues
subroutine meridional_upwind_fluxes(C, C_ady, vh, V, Nx, Ny, Nz, nm )

  implicit none
  real, intent(in) :: C(:, :, :, :)        !< Tracer to calculate numerical mixing for
  real, intent(in) :: C_ady(:, :, :, :)    !< Explicit zonal advection of tracer C
  real, intent(in) :: vh(:, :, :, :)       !< Meridional transport
  real, intent(in) :: V(:, :, :, :)        !< Area of grid cells, computed from V / h

  integer, intent(in) :: Nx, Ny, Nz        !< Dimension lengths
  real, intent(inout) :: nm(:, :, :)

  integer :: i, j, k
  real :: Cupwind(Nx-3, Ny-2, Nz)
  real :: north, south

  call meridional_upwind_values(vh, C, Cupwind)

  do i = 2, Nx-2
    do j = 2, Ny-2
      do k = 1, Nz
        north = 2 * C_ady(i, j, k, 1) * Cupwind(i, j, k) - vh(i, j, k, 1) * Cupwind(i, j, k)**2
        south = 2 * C_ady(i, j-1, k, 1) * Cupwind(i, j-1, k) - vh(i, j-1, k, 1) * Cupwind(i, j-1, k)**2
        nm(i-1, j-1, k) = nm(i-1, j-1, k) + ((north - south) / V(i, j, k, 1))
      enddo
    enddo
  enddo

end subroutine meridional_upwind_fluxes

subroutine meridional_upwind_values(v, C, Cupwind)

  implicit none
  real, intent(in) :: v(:, :, :, :)
  real, intent(in) :: C(:, :, :, :)
  real, intent(inout) :: Cupwind(:, :, :)

  integer :: i, j, k
  integer :: xdim, ydim, zdim

  xdim = size(Cupwind, dim = 1)
  ydim = size(Cupwind, dim = 2)
  zdim = size(Cupwind, dim = 3)

  do i = 1, xdim
    do j = 1, ydim
      do k = 1, zdim
        if (v(i, j, k, 1) >= 0) then
          Cupwind(i, j, k) = C(i, j, k, 1)
        else if (v(i, j, k, 1) < 0) then
          Cupwind(i, j, k) = C(i, j+1, k, 1)
        end if
      end do
    end do
  end do

end subroutine meridional_upwind_values

end module MOM_numerical_mixing