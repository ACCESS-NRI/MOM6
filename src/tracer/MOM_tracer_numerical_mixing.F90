!> Functions and routines involved in calculating numerical mixing of tracers due to advection schemes.
module MOM_tracer_numerical_mixing

use MOM_grid,          only : ocean_grid_type
use MOM_tracer_types,  only : tracer_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public numerical_mixing

contains

!< Calculate the spurious ``numerical'' mixing of tracer due to advection schemes.
subroutine numerical_mixing(G, GV, Tr, h, h_diag, dt_trans, Idt, uhtr, vhtr, nm)

  type(ocean_grid_type),                        intent(in) :: G         !< Ocean grid structure
  type(verticalGrid_type),                      intent(in) :: GV        !< Ocean vertical grid structure
  type(tracer_type),                            intent(in) :: Tr        !< Pointer to the tracer regsitry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h         !< Updated layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h_diag    !< Layer thicknesses prior to
                                                                        !! dynamics [H ~> m or kg m-2]
  real,                                         intent(in) :: dt_trans  !< The transport time interval [T ~> s]
  real,                                         intent(in) :: Idt       !< Inverse time interval [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   intent(in) :: uhtr      !< Accumulated zonal thickness fluxes
                                                                        !! used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   intent(in) :: vhtr      !< Accumulated meridional thickness fluxes
                                                                        !! used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: nm        !< Numerical mixing diagnostic
                                                                        !! [CU2 H T-1 ~> conc2 m s-1]

  ! Upwind variables
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: x_upwind ! zonal upwind values for tracer [CU ~> conc]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: y_upwind ! meridional upwind values for tracer [CU ~> conc]

  x_upwind(:,:,:) = 0.
  y_upwind(:,:,:) = 0.

  call thickness_weighted_variance_advection(G, GV, Tr, h, h_diag, dt_trans, Idt, nm)
  call thickness_weighted_zonal_variance_flux(G, GV, Tr, uhtr, Idt, x_upwind, nm)
  call thickness_weighted_meridional_variance_flux(G, GV, Tr, vhtr, Idt, y_upwind, nm)

end subroutine numerical_mixing

!< Subroutine to calculate the thickness weighted variance advection over the transport timestep.
subroutine thickness_weighted_variance_advection(G, GV, Tr, h, h_diag, dt, Idt, res)

  type(ocean_grid_type),                        intent(in) :: G       !< Ocean grid structure
  type(verticalGrid_type),                      intent(in) :: GV      !< Ocean vertical grid structure
  type(tracer_type),                            intent(in) :: Tr      !< Pointer to the tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h       !< Updated thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h_diag  !< Layer thicknesses prior to
                                                                      !! dynamics [H ~> m or kg m-2]
  real,                                         intent(in) :: dt      !< Transport time interval [T ~> s]
  real,                                         intent(in) :: Idt     !< Inverse time interval [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: res     !< Array to store result in
                                                                      !! [CU2 H T-1 ~> conc2 m s-1]

  !< Local variables
  integer :: is, ie, js, je, nz        !< Grid cell centre and layer indexes
  integer :: i, j, k                   !< Counters
  real :: h_prev, C_prev, Ihadv, Cadv  !< Thickness and tracer at previous timestep, inverse
                                       !< updated thickness and non-invers tracer after advection.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do k=1,nz ; do j=js,je ; do i=is,ie
    h_prev = h_diag(i,j,k)
    Ihadv = 1 / h(i,j,k)
    C_prev = Tr%t_prev(i,j,k)
    Cadv = h_prev * C_prev + dt * Tr%advection_xy(i,j,k)
    res(i,j,k) = ( (Ihadv * Cadv**2) - (h_prev * C_prev**2) ) * Idt
  enddo ; enddo ; enddo

end subroutine thickness_weighted_variance_advection

!< Subroutine to calculate the thickness weigthed zonal variance flux. The spatial derivatives are calucated
!! from upwind values.
subroutine thickness_weighted_zonal_variance_flux(G, GV, Tr, uhtr, Idt, x_upwind, res)

  type(ocean_grid_type),                         intent(in) :: G         !< Ocean grid structure
  type(verticalGrid_type),                       intent(in) :: GV        !< Ocean vertical grid structure
  type(tracer_type),                             intent(in) :: Tr        !< Pointer to the tracer registry
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),    intent(in) :: uhtr      !< Accumulated zonal thickness fluxes
                                                                         !! used to advect tracers [H L2 ~> m3 or kg]
  real,                                          intent(in) :: Idt       !< Inverse time interval [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: x_upwind  !< Zonal upwind tracer [CU ~> conc]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: res       !< Array to store result in
                                                                         !![CU2 H T-1 ~> conc2 m s-1]

  !< Local variables
  integer :: is, ie, js, je, nz  !< Grid cell centre and layer indexes
  integer :: i, j, k             !< Counters
  real :: east, west             !< East and West for zonal derivative [CU2 H T-1 ~> conc2 m s-1]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call zonal_upwind_values(G, GV, Tr, uhtr, x_upwind)

  do k=1,nz ;  do j=js,je ; do i=is,ie
    east = (2 * (Tr%ad_x(I,j,k)  *x_upwind(I,j,k)))   - ((Idt*uhtr(I,j,k))   * (x_upwind(I,j,k)  *x_upwind(I,j,k)))
    west = (2 * (Tr%ad_x(I-1,j,k)*x_upwind(I-1,j,k))) - ((Idt*uhtr(I-1,j,k)) * (x_upwind(I-1,j,k)*x_upwind(I-1,j,k)))
    res(i,j,k) = res(i,j,k) + ((east - west) * G%IareaT(i,j))
  enddo ; enddo; enddo

end subroutine thickness_weighted_zonal_variance_flux

!< Subroutine to calculate upwind values in zonal direction
subroutine zonal_upwind_values(G, GV, Tr, uhtr, x_upwind)

  type(ocean_grid_type),                         intent(in) :: G         !< Ocean grid structure
  type(verticalGrid_type),                       intent(in) :: GV        !< Ocean vertical grid structure
  type(tracer_type),                             intent(in) :: Tr        !< Pointer to the tracer registry
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),    intent(in) :: uhtr      !< Accumulated zonal thickness fluxes
                                                                         !! used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: x_upwind  !< zonal upwind tracer [CU ~> conc]

  !< Local variables
  integer :: is, ie, js, je, nz  !< Grid cell centre indexes
  integer :: i, j, k             !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do k=1,nz ;  do j=js,je ; do I=is-1,ie
    if (uhtr(I,j,k) >= 0) then
      x_upwind(I,j,k) = Tr%t_prev(i,j,k)
    elseif (uhtr(I,j,k) < 0) then
      x_upwind(I,j,k) = Tr%t_prev(i+1,j,k)
    endif
  enddo ; enddo ; enddo

end subroutine zonal_upwind_values

!< Subroutine to calculate the thickness weighted meriodional variance flux. The spatial derivative is calculated
!! from upwind values.
subroutine thickness_weighted_meridional_variance_flux(G, GV, Tr, vhtr, Idt, y_upwind, res)

  type(ocean_grid_type),                         intent(in) :: G         !< Ocean grid structure
  type(verticalGrid_type),                       intent(in) :: GV        !< Ocean vertical grid structure
  type(tracer_type),                             intent(in) :: Tr        !< Pointer to tracer registry
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),    intent(in) :: vhtr      !< Accumulated meridional thickness fluxes used
                                                                        !! to advect tracers [H L2 ~> m3 or kg]
  real,                                          intent(in) :: Idt       !< Inverse timestep [T-1 ~> S-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: y_upwind  !< Meridional upwind tracer values [CU ~> conc]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: res       !< Array to store result in
                                                                        !![CU2 H T-1 ~> conc2 m s-1]

  !< Local variables
  integer :: is, ie, js, je, nz  !< Grid cell centre and layer indexes
  integer :: i, j, k             !< Counters
  real :: north, south           !< North and South positions for meridional derivative

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call meridional_upwind_values(G, GV, Tr, vhtr, y_upwind)

  do k=1,nz ; do j=js,je ; do i=is,ie
    north = (2 * (Tr%ad_y(i,J,k)  *y_upwind(i,J,k)))   - ((Idt*vhtr(i,J,k))   * (y_upwind(i,J,k)  *y_upwind(i,J,k)))
    south = (2 * (Tr%ad_y(i,J-1,k)*y_upwind(i,J-1,k))) - ((Idt*vhtr(i,J-1,k)) * (y_upwind(i,J-1,k)*y_upwind(i,J-1,k)))
    res(i,j,k) = res(i,j,k) + ((north - south) * G%IareaT(i,j))
  enddo ; enddo ; enddo

end subroutine thickness_weighted_meridional_variance_flux

!< Subroutine to calculate upwind value in the meridional direction
subroutine meridional_upwind_values(G, GV, Tr, vhtr, y_upwind)

  type(ocean_grid_type),                        intent(in) :: G         !< Ocean grid structure
  type(verticalGrid_type),                      intent(in) :: GV        !< Ocean vertical grid structure
  type(tracer_type),                            intent(in) :: Tr        !< Pointer to tracer registry
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   intent(in) :: vhtr      !< Accumulated meridional thickness fluxes used
                                                                        !! to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),intent(inout) :: y_upwind  !< Meridional upwind tracer values [CU ~> conc]

  !< Local variables
  integer :: is, ie, js, je, nz  !< Grid cell centre indexes
  integer :: i, j, k             !< Counters

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do k=1,nz ; do J=js-1,je ; do i=is,ie
    if (vhtr(i,J,k) >= 0) then
      y_upwind(i,J,k) = Tr%t_prev(i,j,k)
    elseif (vhtr(i,J,k) < 0) then
      y_upwind(i,J,k) = Tr%t_prev(i,j+1,k)
    endif
  enddo ; enddo ; enddo

end subroutine meridional_upwind_values

end module MOM_tracer_numerical_mixing
