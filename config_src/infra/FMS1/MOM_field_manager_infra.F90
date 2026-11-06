! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module wraps the FMS field manager module
module MOM_field_manager_infra

use field_manager_mod, only : FMS_field_manager_init => field_manager_init

implicit none ; private

public :: field_manager_init

contains

!> Read the field table and initialize the FMS field manager.
subroutine field_manager_init()
  call FMS_field_manager_init()
end subroutine field_manager_init

end module MOM_field_manager_infra
