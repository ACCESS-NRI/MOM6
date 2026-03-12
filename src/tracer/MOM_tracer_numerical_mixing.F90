!> Functions and routines involved in calculating numerical mixing of tracers due to advection
module MOM_tracer_numerical_mixing

use MOM_diag_mediator, only : diag_ctrl, diag_grid_storage
use MOM_grid,          only : ocean_grid_type
use MOM_tracer_types,  only : tracer_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public numerical_mixing, variance_advection, variance_flux, east_west_upoints, Tr_east_west_upoints

contains

!< Calculate the spurious ``numerical'' mixing of tracer due to advection.
subroutine numerical_mixing(G, GV, Tr, h, diag_pre_dyn, dt_trans, Idt, uhtr, vhtr, nm)

  type(ocean_grid_type),   intent(in) :: G                !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV               !< Ocean vertical grid structure
  type(tracer_type),       intent(in) :: Tr               !< Pointer to the tracer regsitry
  real,                    intent(in) :: h(:,:,:)         !< The updated layer thicknesses [H ~> m or kg m-2]
  type(diag_grid_storage), intent(in) :: diag_pre_dyn     !< Stored grids from before dynamics
  real,                    intent(in) :: dt_trans         !< The transport time interval [T ~> s]
  real,                    intent(in) :: Idt              !< The inverse of the time interval [T-1 ~> s-1]
  real,                    intent(in) :: uhtr(:,:,:)      !< Accumulated zonal thickness fluxes
                                                          !! used to advect tracers [H L2 ~> m3 or kg]
  real,                    intent(in) :: vhtr(:,:,:)      !< Accumulated meridional thickness fluxes
                                                          !! used to advect tracers [H L2 ~> m3 or kg]
  real,                 intent(inout) :: nm(:,:,:)        !< Numerical mixing diagnostic [CU2 H T-1 ~> conc2 m s-1]

  ! Upwind variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: x_upwind ! zonal upwind values for tracer [CU ~> conc]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: y_upwind ! meridional upwind values for tracer [CU ~> conc]

  x_upwind(:,:,:) = 0.
  y_upwind(:,:,:) = 0.

  call thickness_weighted_variance_advection(Tr, h, diag_pre_dyn, dt_trans, Idt, G, GV, nm)
  call thickness_weighted_zonal_variance_flux(Tr, uhtr, G, GV, Idt, x_upwind, nm)
  call thickness_weighted_meridional_variance_flux(Tr, vhtr, G, GV, Idt, y_upwind, nm)

end subroutine numerical_mixing

!< Subroutine for the variance advection, likely will remove once numerical mixing is sorted out
subroutine variance_advection(G, GV, Tr, h, diag_pre_dyn, dt_trans, Idt, va)

  type(ocean_grid_type),   intent(in) :: G             !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV            !< Ocean vertical grid structure
  type(tracer_type),       intent(in) :: Tr            !< Pointer to the tracer regsitry
  real,                    intent(in) :: h(:,:,:)      !< The updated layer thicknesses [H ~> m or kg m-2]
  type(diag_grid_storage), intent(in) :: diag_pre_dyn  !< Stored grids from before dynamics
  real,                    intent(in) :: dt_trans      !< The transport time interval [T ~> s]
  real,                    intent(in) :: Idt           !< The inverse of the time interval [T-1 ~> s-1]
  real,                 intent(inout) :: va(:,:,:)     !< Thickness weighted variance advection
                                                       !! [CU2 H T-1 ~> conc2 m s-1]

  call thickness_weighted_variance_advection(Tr, h, diag_pre_dyn, dt_trans, Idt, G, GV, va)

end subroutine variance_advection

!< Subroutine for the horizontal variance flux, likely will remove once numerical mixing is sorted out
subroutine variance_flux(G, GV, Tr, Idt, uhtr, vhtr, vf)

  type(ocean_grid_type),   intent(in) :: G                !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV               !< Ocean vertical grid structure
  type(tracer_type),       intent(in) :: Tr               !< Pointer to the tracer regsitry
  real,                    intent(in) :: Idt              !< The inverse of the time interval [T-1 ~> s-1]
  real,                    intent(in) :: uhtr(:,:,:)      !< Accumulated zonal thickness fluxes
                                                          !! used to advect tracers [H L2 ~> m3 or kg]
  real,                    intent(in) :: vhtr(:,:,:)      !< Accumulated meridional thickness fluxes
                                                          !! used to advect tracers [H L2 ~> m3 or kg]
  real,                 intent(inout) :: vf(:,:,:)        !< Horizontal thickness weighted variance flux
                                                          !! [CU2 H T-1 ~> conc2 m s-1]
  ! Try with local variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: x_upwind ! zonal upwind values for tracer [CU ~> conc]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: y_upwind ! meridional upwind values for tracer [CU ~> conc]

  x_upwind(:,:,:) = 0.
  y_upwind(:,:,:) = 0.

  call thickness_weighted_zonal_variance_flux(Tr, uhtr, G, GV, Idt, x_upwind, vf)
  call thickness_weighted_meridional_variance_flux(Tr, vhtr, G, GV, Idt, y_upwind, vf)

end subroutine variance_flux

!< Subroutine to calculate the thickness weighted variance advection over the transport timestep.
subroutine thickness_weighted_variance_advection(Tr, h, diag_pre_dyn, dt, Idt, G, GV, res)

  type(tracer_type),       intent(in) :: Tr            !< Pointer to the tracer registry
  real,                    intent(in) :: h(:,:,:)      !< The updated layer thicknesses [H ~> m or kg m-2]
  type(diag_grid_storage), intent(in) :: diag_pre_dyn  !< Stored grids from before dynamics
  real,                    intent(in) :: dt            !< The transport time interval [T ~> s]
  real,                    intent(in) :: Idt           !< The inverse of the time interval [T-1 ~> s-1]
  type(ocean_grid_type),   intent(in) :: G             !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV            !< Ocean vertical grid structure
  real,                 intent(inout) :: res(:,:,:)    !< Array to store result in [CU2 H T-1 ~> conc2 m s-1]

  !< Local variables
  integer :: is, ie, js, je, nz        !< Grid cell centre and layer indexes
  integer :: i, j, k                   !< Counters
  real :: h_prev, C_prev, Ihadv, Cadv  !< Thickness and tracer at previous timestep, inverse
                                       !< updated thickness and non-invers tracer after advection.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do k=1,nz ; do j=js,je ; do i=is,ie
    h_prev = diag_pre_dyn%h_state(i,j,k)
    Ihadv = 1 / h(i,j,k)
    C_prev = Tr%t_prev(i,j,k)
    Cadv = h_prev * C_prev + dt * Tr%advection_xy(i,j,k)
    res(i,j,k) = ( (Ihadv * Cadv**2) - (h_prev * C_prev**2) ) * Idt
  enddo ; enddo ; enddo

end subroutine thickness_weighted_variance_advection

!< Subroutine to calculate the thickness weigthed zonal variance flux. The spatial derivatives are calucated
!! from upwind values.
subroutine thickness_weighted_zonal_variance_flux(Tr, uhtr, G, GV, Idt, x_upwind, res)

  type(tracer_type),       intent(in) :: Tr               !< Pointer to the tracer registry
  real,                    intent(in) :: uhtr(:,:,:)      !< Accumulated zonal thickness fluxes
                                                          !! used to advect tracers [H L2 ~> m3 or kg]
  type(ocean_grid_type),   intent(in) :: G                !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV               !< Ocean vertical grid structure
  real,                    intent(in) :: Idt              !< Inverse transport time intervale [T-1 ~> s-1]
  real,                 intent(inout) :: x_upwind(:,:,:)  !< Zonal upwind tracer value [CU ~> conc]
  real,                 intent(inout) :: res(:,:,:)       !< Array to store the result in [CU2 H T-1 ~> conc2 m s-1]

  !< Local variables
  integer :: is, ie, js, je, nz           !< Grid cell centre and layer indexes
  integer :: i, j, k                      !< Counters
  real :: east, west                      !< East and West for zonal derivative [CU2 H T-1 ~> conc2 m s-1]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call zonal_upwind_values(Tr, G, nz, uhtr, x_upwind)

  do k=1,nz ;  do j=js,je ; do i=is,ie
    east = (2 * (Tr%ad_x(I,j,k)  *x_upwind(I,j,k)))   - ((Idt*uhtr(I+1,j,k)) * (x_upwind(I,j,k)  *x_upwind(I,j,k)))
    west = (2 * (Tr%ad_x(I-1,j,k)*x_upwind(I-1,j,k))) - ((Idt*uhtr(I,j,k))   * (x_upwind(I-1,j,k)*x_upwind(I-1,j,k)))
    res(i,j,k) = res(i,j,k) + ((east - west) * G%IareaT(i,j))
    ! This code passes the thickness dimensional test but does not accurately calculate numerical mixing
    ! east = (2 * Tr%ad_x(I,j,k)   * x_upwind(I,j,k))   - (Idt * uhtr(I,j,k) * x_upwind(I,j,k)**2)
    ! west = (2 * Tr%ad_x(I-1,j,k) * x_upwind(I-1,j,k)) - (Idt * uhtr(I-1,j,k)   * x_upwind(I-1,j,k)**2)
  enddo ; enddo; enddo

end subroutine thickness_weighted_zonal_variance_flux

!< Subroutine to calculate upwind values in zonal direction
subroutine zonal_upwind_values(Tr, G, nz, uhtr, x_upwind)

  type(tracer_type),        intent(in) :: Tr               !< Pointer to the tracer registry
  type(ocean_grid_type),    intent(in) :: G                !< Ocean grid structure
  integer,                  intent(in) :: nz               !< Vertical extent of domain
  real,                     intent(in) :: uhtr(:,:,:)      !< Accumulated zonal transport
                                                           !! used to advect tracers [H L2 ~> m3 or kg]
  real,                  intent(inout) :: x_upwind(:,:,:)  !< Zonal upwind values [CU ~> conc]

  !< Local variables
  integer :: is, ie, js, je  !< Grid cell centre indexes
  integer :: i, j, k         !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do k=1,nz ;  do j=js,je ; do I=is-1,ie
    if (uhtr(I+1,j,k) >= 0) then
      x_upwind(I,j,k) = Tr%t_prev(i,j,k)
    elseif (uhtr(I+1,j,k) < 0) then
      x_upwind(I,j,k) = Tr%t_prev(i+1,j,k)
    endif
  enddo ; enddo ; enddo

end subroutine zonal_upwind_values

!< Subroutine to calculate the thickness weighted meriodional variance flux. The spatial derivative is calculated
!! from upwind values.
subroutine thickness_weighted_meridional_variance_flux(Tr, vhtr, G, GV, Idt, y_upwind, res)

  type(tracer_type),       intent(in) :: Tr               !< Tracer
  real,                    intent(in) :: vhtr(:,:,:)      !< Accumulated meridional thickness fluxes
                                                          !! used to advect tracers [H L2 ~> m3 or kg]
  type(ocean_grid_type),   intent(in) :: G                !< Ocean grid structure
  type(verticalGrid_type), intent(in) :: GV               !< Ocean vertical grid structure
  real,                    intent(in) :: Idt              !< Inverse model timestep
  real,                 intent(inout) :: y_upwind(:,:,:)  !< Meridional upwind tracer values [CU ~> conc]
  real,                 intent(inout) :: res(:,:,:)       !< Array to store the result in [CU2 H T-1 ~> conc2 m s-1]

  !< Local variables
  integer :: is, ie, js, je, nz           !< Grid cell centre and layer indexes
  integer :: i, j, k                      !< Counters
  real :: north, south                    !< North and South positions for meridional derivative

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call meridional_upwind_values(Tr, G, nz, vhtr, y_upwind)

  do k=1,nz ; do j=js,je ; do i=is,ie
    north = (2 * (Tr%ad_y(i,J,k)  * y_upwind(i,J,k)))   - ((Idt*vhtr(i,J+1,k)) * (y_upwind(i,J,k)  *y_upwind(i,J,k)))
    south = (2 * (Tr%ad_y(i,J-1,k)* y_upwind(i,J-1,k))) - ((Idt*vhtr(i,J,k))   * (y_upwind(i,J-1,k)*y_upwind(i,J-1,k)))
    res(i,j,k) = res(i,j,k) + ((north - south) * G%IareaT(i,j))
    ! This code passes the thickness dimensional test but is not correct for the numerical mixing diagnostic
    ! north = (2 * Tr%ad_y(i,J,k)   * y_upwind(i,J,k))   - (Idt * vhtr(i,J,k) * y_upwind(i,J,k)**2)
    ! south = (2 * Tr%ad_y(i,J-1,k) * y_upwind(i,J-1,k)) - (Idt * vhtr(i,J-1,k)   * y_upwind(i,J-1,k)**2)
  enddo ; enddo ; enddo

end subroutine thickness_weighted_meridional_variance_flux

!< Subroutine to calculate upwind value in the meridional direction
subroutine meridional_upwind_values(Tr, G, nz, vhtr, y_upwind)

  type(tracer_type),     intent(in) :: Tr               !< Tracer
  type(ocean_grid_type), intent(in) :: G                !< Ocean grid structure for inverse area
  integer,               intent(in) :: nz               !< Grid cell layer indexes
  real,                  intent(in) :: vhtr(:,:,:)      !< Accumulated meridional thickness fluxes
                                                        !! used to advect tracers [H L2 ~> m3 or kg]
  real,               intent(inout) :: y_upwind(:,:,:)  !< Meridional upwind values of C [CU ~> conc]

  !< Local variables
  integer :: is, ie, js, je  !< Grid cell centre indexes
  integer :: i, j, k         !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do k=1,nz ; do J=js-1,je ; do i=is,ie
    if (vhtr(i,J+1,k) >= 0) then
      y_upwind(i,J,k) = Tr%t_prev(i,j,k)
    elseif (vhtr(i,J+1,k) < 0) then
      y_upwind(i,J,k) = Tr%t_prev(i,j+1,k)
    endif
  enddo ; enddo ; enddo

end subroutine meridional_upwind_values

! Subroutine to construct a MWE that demonstrates the hack that I have above where I have the `uhtr` variable
! indexed one value higher than `TR%ad_x`. I should then be able to compare this output to the equivalent
! saved variables and demonstrate that the index online and offline does not seem to match. Or in doing this I will
! find the error in my ways!
subroutine east_west_upoints(uhtr, G, GV, uhtr_eu, uhtr_wu)

  real,  dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(in) :: uhtr !< the variable to save the east faces of
  type(ocean_grid_type), intent(in) :: G          !< Ocean grid structure for indexes
  type(verticalGrid_type), intent(in)    :: GV     !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr_eu !< The east u point value of `var`
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhtr_wu !< The west u point value of `var`

  !< Local variables
  integer :: is, ie, js, je  !< Grid cell centre indexes
  integer :: i, j, k         !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do k=1,nz ;  do j=js,je ; do i=is,ie
    uhtr_wu(i,j,k) = uhtr(I-1,j,k)
    uhtr_eu(i,j,k) = uhtr(I,j,k)
  enddo ; enddo; enddo

end subroutine east_west_upoints

subroutine Tr_east_west_upoints(Tr, G, nz, adx_eu, adx_wu)

  type(tracer_type),     intent(in) :: Tr               !< Tracer
  type(ocean_grid_type), intent(in) :: G          !< Ocean grid structure for indexes
  integer,               intent(in) :: nz         !< number of vertical levels

  real,                  intent(inout) :: adx_eu(:,:,:) !< The east u point value of `var`
  real,                  intent(inout) :: adx_wu(:,:,:) !< The west u point value of `var`

  !< Local variables
  integer :: is, ie, js, je  !< Grid cell centre indexes
  integer :: i, j, k         !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  do k=1,nz ;  do j=js,je ; do i=is,ie
    adx_wu(i,j,k) = Tr%ad_x(I-1,j,k)
    adx_eu(i,j,k) = Tr%ad_x(I,j,k)
  enddo ; enddo; enddo

end subroutine Tr_east_west_upoints

end module MOM_tracer_numerical_mixing
