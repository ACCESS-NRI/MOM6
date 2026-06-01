!> Functions and routines involved in calculating numerical mixing of tracers due to advection schemes.
module MOM_tracer_numerical_mixing

use MOM_grid,          only : ocean_grid_type
use MOM_tracer_types,  only : tracer_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public advection_scheme_variance_production

contains

!< Calculate the spurious variance production of tracer `Tr` due to the advection schemes.
subroutine advection_scheme_variance_production(G, GV, Tr, h_diag, h, dt_trans, Idt, uhtr, vhtr, asvp)

  type(ocean_grid_type),                        intent(in) :: G         !< Ocean grid structure
  type(verticalGrid_type),                      intent(in) :: GV        !< Ocean vertical grid structure
  type(tracer_type),                            intent(in) :: Tr        !< Pointer to the tracer regsitry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h_diag    !< Previous thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h         !< Updated layer thicknesses [H ~> m or kg m-2]
  real,                                         intent(in) :: dt_trans  !< The transport time interval [T ~> s]
  real,                                         intent(in) :: Idt       !< Inverse time interval [T-1 ~> s-1]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),   intent(in) :: uhtr      !< Accumulated zonal thickness fluxes
                                                                        !! used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),   intent(in) :: vhtr      !< Accumulated meridional thickness fluxes
                                                                        !! used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: asvp      !< Advection scheme varianve production
                                                                        !! diagnostic [CU2 H T-1 ~> conc2 m s-1]

  ! call thickness_weighted_variance_advection(G, GV, Tr, h_diag, h, dt_trans, Idt, asvp)
  call thickness_weighted_variance_flux_divergence(G, Gv, Tr, uhtr, vhtr, Idt, asvp)

end subroutine advection_scheme_variance_production

!< Subroutine to calculate the thickness weighted variance advection over the transport timestep.
subroutine thickness_weighted_variance_advection(G, GV, Tr, h_diag, h, dt, Idt, asvp)

  type(ocean_grid_type),                        intent(in) :: G       !< Ocean grid structure
  type(verticalGrid_type),                      intent(in) :: GV      !< Ocean vertical grid structure
  type(tracer_type),                            intent(in) :: Tr      !< Pointer to the tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h_diag  !< Previous thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),    intent(in) :: h       !< Updated thicknesses [H ~> m or kg m-2]
  real,                                         intent(in) :: dt      !< Transport time interval [T ~> s]
  real,                                         intent(in) :: Idt     !< Inverse time interval [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: asvp    !< Array to store asvp in
                                                                      !! [CU2 H T-1 ~> conc2 m s-1]

  !< Local variables
  integer :: is, ie, js, je, nz  !< Grid cell centre and layer indexes
  integer :: i, j, k             !< Counters
  real :: Ih                     !< Inverse updated thickness [H-1 ~> m-1]
  real :: ht_prev                !< Thickness weighted tracer prior to dynamics [CU H ~> conc m]
  real :: ht_adv                 !< Thickness weighted tracer after lateral advection [CU H ~> conc m]

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  do k=1,nz ; do j=js,je ; do i=is,ie
    Ih = 1 / h(i,j,k)
    ht_prev = h_diag(i,j,k) * Tr%t_prev(i,j,k)
    ht_adv = ht_prev + dt * Tr%advection_xy(i,j,k)
    asvp(i,j,k) = ( (Ih * (ht_adv*ht_adv)) - (ht_prev * Tr%t_prev(i,j,k)) ) * Idt
  enddo ; enddo ; enddo

end subroutine thickness_weighted_variance_advection

!< Subroutine to calculate the divergence of the variance flux due to the advection scheme
subroutine thickness_weighted_variance_flux_divergence(G, Gv, Tr, uhtr, vhtr, Idt, asvp)

  type(ocean_grid_type),                         intent(in) :: G     !< Ocean grid structure
  type(verticalGrid_type),                       intent(in) :: GV    !< Ocean vertical grid structure
  type(tracer_type),                             intent(in) :: Tr    !< Pointer to the tracer registry
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),    intent(in) :: uhtr  !< Accumulated zonal thickness fluxes
                                                                     !! used to advect tracers [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),    intent(in) :: vhtr  !< Accumulated meridional thickness fluxes
                                                                     !! used to advect tracers [H L2 ~> m3 or kg]
  real,                                          intent(in) :: Idt   !< Inverse time interval [T-1 ~> s-1]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: asvp  !< Array to store asvp in
                                                                     !! [CU2 H T-1 ~> conc2 m s-1]

  !< Local variables
  integer :: is, ie, js, je, nz                           !< Grid cell centre and layer indexes
  integer :: i, j, k                                      !< Counters
  real :: east, west, north, south                        !< east, west, north and south faces of h point
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: x_upwind  !< Zonal upwind values for tracer [CU ~> conc]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: y_upwind  !< Meridional upwind values for tracer [CU ~> conc]

  x_upwind(:,:,:) = 0.
  y_upwind(:,:,:) = 0.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call zonal_upwind_values(G, GV, Tr, uhtr, x_upwind)
  call meridional_upwind_values(G, GV, Tr, vhtr, y_upwind)

  do k=1,nz ;  do j=js,je ; do i=is,ie
     east = (2 * (Tr%ad_x(I,j,k)  *x_upwind(I,j,k)))   - ((Idt*uhtr(I,j,k))   * (x_upwind(I,j,k)  *x_upwind(I,j,k)))
     west = (2 * (Tr%ad_x(I-1,j,k)*x_upwind(I-1,j,k))) - ((Idt*uhtr(I-1,j,k)) * (x_upwind(I-1,j,k)*x_upwind(I-1,j,k)))
    north = (2 * (Tr%ad_y(i,J,k)  *y_upwind(i,J,k)))   - ((Idt*vhtr(i,J,k))   * (y_upwind(i,J,k)  *y_upwind(i,J,k)))
    south = (2 * (Tr%ad_y(i,J-1,k)*y_upwind(i,J-1,k))) - ((Idt*vhtr(i,J-1,k)) * (y_upwind(i,J-1,k)*y_upwind(i,J-1,k)))
    asvp(i,j,k) = asvp(i,j,k) + (((east - west) + (north - south)) * G%IareaT(i,j))
  enddo ; enddo; enddo

end subroutine thickness_weighted_variance_flux_divergence

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

!< Subroutine to calculate upwind values in the meridional direction
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
