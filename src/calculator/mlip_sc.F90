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

!> module mlip_sc
!> A module containing routines for calling MLIPs though persistent python instances
!> enabled through the fortbridge submodule

!=========================================================================================!
module mlip_sc
  use iso_fortran_env,only:wp => real64
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove,command
#ifdef WITH_FORTBRIDGE
  use fortbridge_client
#endif
  implicit none
  !>--- private module variables and parameters
  private


  public :: mlip_engrad

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine mlip_engrad(mol,energy,gradient,iostatus)
    type(coord),intent(in) :: mol
    real(wp),intent(out)   :: energy
    real(wp),intent(out)   :: gradient(3,mol%nat)
    integer,intent(out)    :: iostatus

    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    iostatus = 1

  end subroutine mlip_engrad

!========================================================================================!
end module mlip_sc
