!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

!> a small module for getting a free energy as engrad call

module approxg_module
  use crest_parameters
  use calc_type
  use modelhessian_core
  use thermochem_module
  use strucrd
  implicit none
  private

  public :: modh_engrad

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine modh_engrad(mol,calc,dg,dggrad)
    type(coord),intent(in) :: mol
    type(calculation_settings),intent(inout) :: calc
    real(wp),intent(out) :: dg
    real(wp),intent(out) :: dggrad(:)
    integer :: n3

    type(mhparam) :: mhset

    dg = 0.0_wp
    dggrad(:) = 0.0_wp

    !> setup
    n3 = mol%nat*3
    if (calc%approxg_dim .ne. mol%nat) then
      !$omp critical
      if (allocated(calc%approxg_hess)) deallocate (calc%approxg_hess)
      allocate (calc%approxg_hess(n3,n3),source=0.0_wp)

      if (allocated(calc%approxg_h)) deallocate (calc%approxg_h) 
      allocate (calc%approxg_h(n3*(n3+1)/2),source=0.0_wp)


      calc%approxg_dim = mol%nat
      !$omp end critical
    else 
      calc%approxg_hess(:,:) = 0.0_wp
      calc%approxg_h(:) = 0.0_wp
    end if

    call ddvopt(mol%xyz,mol%nat,calc%approxg_h,mol%at,mhset)

  end subroutine modh_engrad

!========================================================================================!
!========================================================================================!
end module approxg_module
