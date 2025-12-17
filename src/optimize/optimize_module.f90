!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

!> This module wrapps the different optimization algorithms,
!> i.e., this is what can be called for geometry opt.

module optimize_module
  use iso_fortran_env,only:wp => real64
  use crest_parameters
  use crest_calculator
  use strucrd
  use ancopt_module
  use gradientdescent_module
  use rfo_module
  use optimize_utils
  use thermochem_module
  use hessian_reconstruct
  use hessian_tools
  implicit none
  private

  public :: optimize_geometry
  public :: print_opt_data

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine optimize_geometry(mol,molnew,calc,etot,grd,pr,wr,iostatus)
    implicit none
    !> Input
    type(coord)    :: mol
    type(calcdata) :: calc
    logical,intent(in)        :: pr
    logical,intent(in)        :: wr
    !> Output
    type(coord)   :: molnew
    integer,intent(out)       :: iostatus
    real(wp),intent(inout)    :: etot
    real(wp),intent(inout)    :: grd(3,mol%nat)
    real(wp),allocatable :: H_inv(:,:), freq(:)
    integer :: nat3
    integer :: io

    iostatus = -1
    !> do NOT overwrite original geometry
    !$omp critical
    molnew%at = mol%at
    molnew%xyz = mol%xyz
    molnew%nat = mol%nat
    !$omp end critical
    nat3 = 3*mol%nat

    !> Check for optimization-individual calculation setup
    if (calc%optnewinit) then
      !$omp critical
      call calc%dealloc_params()
      !$omp end critical
    end if

    !> Check if Hessian Reconstruct is called
    if (calc%do_HU) then
      allocate (calc%chess)
      call calc%chess%alloc(mol%nat,calc%hu_steps,calc%hguess)
    end if

    !> initial singlepoint
    call engrad(molnew,calc,etot,grd,iostatus)

    !> optimization
    select case (calc%opt_engine)
    case (0)
      call ancopt(molnew,calc,etot,grd,pr,wr,iostatus)
    case (1)
      !> l-bfgs goes here
      write (stdout,'(a)') 'L-BFGS currently not implemented'
      stop
    case (2)
      !> rfo goes here
      call rfopt(molnew,calc,etot,grd,pr,wr,iostatus)
    case (-1)
      call gradientdescent(molnew,calc,etot,grd,pr,wr,iostatus)
    case default
      write (stdout,'(a)') 'Unknown optimization engine!'
      stop
    end select
    molnew%energy = etot

    if (calc%do_HU) then !> Hessian construction and post-processing happen here
      !print*, "Energies", calc%chess%energy
      !print*, "Gradients", calc%chess%gradient
      !print*, "Coords", calc%chess%coords
      !print*, "Order", calc%chess%order

      call calc%chess%construct_hessian_bfgs()

      !allocate(H_inv(size(calc%chess%B,1),size(calc%chess%B,2)))
      !H_inv(:,:) = invert_matrix(calc%chess%B)

      print*
      print*,"THERMO FROM MY OWN SHITTY HESSIAN"
      print*

      call calc_thermo_from_hess(molnew,calc%chess%B,pr, &
      & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
      & calc%ht,calc%gt,calc%stot)

      print*
      print*,"THERMO FROM BFGS"
      print*

      !call calc_thermo_from_hess(molnew,calc%chess%H,pr, &
      !& calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
      !& calc%ht,calc%gt,calc%stot)

      call mass_weight_hess(molnew%nat,molnew%at,nat3,calc%chess%H(:,:))

      allocate(freq(nat3))

      call frequencies(molnew%nat,molnew%at,molnew%xyz,nat3,calc%chess%H(:,:),freq,io)

      call calcthermo(molnew%nat,molnew%at,mol%xyz,freq,pr,calc%ithr,calc%fscal,calc%sthr, &
          & calc%nt,calc%temperatures, &
          &      calc%et,calc%ht,calc%gt,calc%stot)

      !write(stdout,*) "et:", calc%et
      !write(stdout,*) "ht:", calc%ht
      !write(stdout,*) "gt:", calc%gt
      !write(stdout,*) "stot:", calc%stot

      call calc%chess%dealloc()
      deallocate (calc%chess)
    end if

    return
  end subroutine optimize_geometry

!========================================================================================!

  subroutine print_opt_data(calc,ich)
    implicit none
    type(calcdata) :: calc
    integer,intent(in) :: ich
    integer :: tight
    real(wp) :: ethr,gthr

    write (ich,'(1x,a)',advance='no') 'Optimization engine: '
    select case (calc%opt_engine)
    case (0)
      write (ich,'(a)') 'ANCOPT'
    case (1)
      write (ich,'(a)') 'L-BFGS'
    case (2)
      write (ich,'(a)') 'RFO'
    case (-1)
      write (ich,'(a)') 'Gradient Descent'
    case default
      write (ich,'(a)') 'Unknown'
    end select
    if (calc%opt_engine >= 0) then
      write (ich,'(1x,a)',advance='no') 'Hessian update type: '
      select case (calc%iupdat)
      case (0)
        write (ich,'(a)') 'BFGS'
      case (1)
        write (ich,'(a)') 'Powell'
      case (2)
        write (ich,'(a)') 'SR1'
      case (3)
        write (ich,'(a)') 'Bofill'
      case (4)
        write (ich,'(a)') 'Farkas-Schlegel'
      end select
    end if

    tight = calc%optlev
    call get_optthr(0,tight,calc,ethr,gthr)
    write (ich,'(1x,a,e10.3,a,e10.3,a)') 'E/G convergence criteria: ',&
    & ethr,' Eh,',gthr,' Eh/a0'

    write (ich,'(1x,a,i0)') 'maximum optimization steps: ',calc%maxcycle

  end subroutine print_opt_data

!========================================================================================!
!========================================================================================!
end module optimize_module
