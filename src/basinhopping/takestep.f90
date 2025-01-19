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

module bh_step_module
  use crest_parameters
  use strucrd,only:coord
  use crest_calculator
  use bh_class_module
  implicit none
  private

  logical,parameter :: debug = .true.
!  logical,parameter :: debug = .false.

  public :: takestep,steptypestr

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  function steptypestr(steptype) result(str)
    implicit none
    integer,intent(in) :: steptype
    character(len=:),allocatable :: str
    select case (steptype)
    case default !> Cartesian
      str = 'Cartesian'
    case (1) !> natural internals
      str = 'internal '
    case (2) !> dihedral only
      str = 'dihedral '
    case (3) !> intermolecular (CMA,tilt) only
      str = 'intermol '
    end select
  end function steptypestr

!=========================================================================================!

  subroutine takestep(mol,calc,bh,newmol)
    implicit none
    !> IN/OUTPUT
    type(coord),intent(in)       :: mol   !> molecular system
    type(calcdata),intent(inout) :: calc
    type(bh_class),intent(inout) :: bh    !> BH settings
    type(coord),intent(out)      :: newmol
    !> LOCAL

    select case (bh%steptype)
    case(0) !> Cartesian
      newmol = mol
      call takestep_cart(newmol,bh%stepsize(1),calc)
    case default
      error stop 'Steptype not implemented yet'
    end select

  end subroutine takestep

!=========================================================================================!

  subroutine takestep_cart(newmol,stepsize,calc)
    implicit none
    type(coord),intent(inout) :: newmol
    real(wp),intent(in) :: stepsize
    type(calcdata),intent(inout) :: calc
    real(wp) :: r(3)
    integer :: i
    do i = 1,newmol%nat
      if (calc%nfreeze > 0) then
        if (calc%freezelist(i)) cycle
      end if
      call random_number(r)
      r(:) = (r(:)-0.5_wp)*2.0_wp
      newmol%xyz(:,i) = newmol%xyz(:,i)+r(:)*stepsize
    end do
  end subroutine takestep_cart

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!
end module bh_step_module
