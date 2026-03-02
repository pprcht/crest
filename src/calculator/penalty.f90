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

!> a small module for getting a penalty contribution, e.g. from the RMSD potential used in metadynamics

module penalty_module
  use crest_parameters
  use strucrd
  use irmsd_module
  implicit none
  private

  public :: penalty_params
  type :: penalty_params
    type(coord),pointer  :: biaslist(:)
    real(wp) :: alpha = 1.0_wp
    real(wp) :: kpush = 0.01_wp
    real(wp),allocatable :: ramp(:)
    real(wp),allocatable :: gradtmp(:,:)
    type(rmsd_core_cache) :: ccache
  end type penalty_params

  public :: rmsd_engrad

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine rmsd_engrad(mol,ppars,energy,grad,iostatus)
    type(coord),intent(in) :: mol
    type(penalty_params),intent(inout) :: ppars
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: grad(:,:)
    integer,intent(out) :: iostatus
    integer :: nall,io,ii
    real(wp) :: etmp,rmsdval,dEdr

    iostatus = 0
    energy = 0.0_wp
    grad(:,:) = 0.0_wp
    rmsdval = 0.0_wp
    nall = size(ppars%biaslist,1)

    do ii = 1,nall

      rmsdval = rmsd(mol,ppars%biaslist(ii),gradient=ppars%gradtmp,ccache=ppars%ccache)

      etmp = ppars%kpush*(-ppars%alpha*rmsdval**2)
      dEdr = -2.0_wp*ppars%alpha*etmp*rmsdval

      energy = energy+etmp
      grad(:,:) = grad(:,:)+dEdr*ppars%gradtmp(:,:)
    end do

  end subroutine rmsd_engrad

!========================================================================================!
!========================================================================================!
end module penalty_module
