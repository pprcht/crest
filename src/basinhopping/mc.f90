!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2024 Philipp Pracht, David Wales
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

module bh_mc_module
  use crest_parameters
  use strucrd,only:coord
  use crest_calculator
  use bh_class_module
  use bh_step_module
  implicit none
  private

  logical,parameter :: debug = .true.
!  logical,parameter :: debug = .false.

  public :: mc

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine mc(calc,mol,bh)
    implicit none
    !> IN/OUTPUT
    type(calcdata),intent(inout) :: calc  !> potential settings
    type(coord),intent(inout)    :: mol   !> molecular system
    type(bh_class),intent(inout) :: bh    !> BH settings
    !> LOCAL
    type(coord) :: tmpmol    !> copy to take steps
    integer :: iter

    do iter = 1,bh%maxsteps

!>--- Take the step 
      call takestep(mol,bh,tmpmol)

!>--- Quench it


!>--- Accept/reject


!>--- Update structures

    end do

  end subroutine mc

!=========================================================================================!

  function mcaccept(mol,bh) result(accept)
!**************************************
!* The regular MC acceptance condition
!**************************************
    implicit none
    logical :: accept
    type(coord),intent(in) :: mol
    type(bh_calss),intent(in) :: bh
    real(wp) :: eold,enew,temp
    real(wp) :: random,fact
    accept = .false.
    eold = bh%emin
    enew = mol%energy
    temp = bh%temp*kB  !> Kelvin to a.u.


    if (enew .lt. eold) then
      accept = .true.
    else
      call random_number(random)
      fact = exp(-(enew-eold)/temp)
      if (fact .gt. random) accept = .true.
    end if

  end function mcaccept

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!
end module bh_mc_module
