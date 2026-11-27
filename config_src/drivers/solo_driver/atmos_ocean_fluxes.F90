! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A dummy version of atmos_ocean_fluxes_mod module for
!! use when the vastly larger FMS package is not needed.
module atmos_ocean_fluxes_mod

use coupler_types_mod, only : coupler_1d_bc_type

implicit none ; private

public :: aof_set_coupler_flux
public :: atmos_ocean_type_fluxes_init, atmos_ocean_fluxes_init

contains

!> This subroutine duplicates an interface used by the FMS coupler, but only
!! returns a value of -1.  None of the arguments are used for anything.
function aof_set_coupler_flux(name, flux_type, implementation, atm_tr_index,     &
                              param, flag, mol_wt, ice_restart_file, ocean_restart_file, &
                              units, caller, verbosity)  result (coupler_index)

  character(len=*),                intent(in) :: name !< An unused argument
  character(len=*),                intent(in) :: flux_type !< An unused argument
  character(len=*),                intent(in) :: implementation !< An unused argument
  integer,               optional, intent(in) :: atm_tr_index !< An unused argument
  real,    dimension(:), optional, intent(in) :: param !< An unused argument that would be used to
                                                   !! pass parameters for flux parameterizations
                                                   !! in other contexts [various]
  logical, dimension(:), optional, intent(in) :: flag !< An unused argument
  real,                  optional, intent(in) :: mol_wt !< An unused argument that would usually be
                                                   !! the tracer's molecular weight [g mol-1]
  character(len=*),      optional, intent(in) :: ice_restart_file !< An unused argument
  character(len=*),      optional, intent(in) :: ocean_restart_file !< An unused argument
  character(len=*),      optional, intent(in) :: units !< An unused argument
  character(len=*),      optional, intent(in) :: caller !< An unused argument
  integer,               optional, intent(in) :: verbosity !< An unused argument

  ! None of these arguments are used for anything.

  integer :: coupler_index
  coupler_index = -1

end function aof_set_coupler_flux

!> This subroutine duplicates an interface used by the FMS coupler, but does nothing.
!! None of the arguments are used for anything.
subroutine atmos_ocean_type_fluxes_init(verbosity)
  integer, optional, intent(in) :: verbosity !< An unused argument
end subroutine atmos_ocean_type_fluxes_init

!> This subroutine duplicates an interface used by the FMS coupler, but does nothing.
!! None of the arguments are used for anything.
subroutine atmos_ocean_fluxes_init(gas_fluxes, gas_fields_atm, gas_fields_ice, verbosity)
  type(coupler_1d_bc_type), intent(inout) :: gas_fluxes !< An unused argument
  type(coupler_1d_bc_type), intent(inout) :: gas_fields_atm !< An unused argument
  type(coupler_1d_bc_type), intent(inout) :: gas_fields_ice !< An unused argument
  integer,        optional, intent(in)    :: verbosity !< An unused argument
end subroutine atmos_ocean_fluxes_init

end module atmos_ocean_fluxes_mod
