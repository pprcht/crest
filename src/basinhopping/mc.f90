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
  use optimize_module
  use axis_module
  use bh_class_module
  use bh_step_module
  use omp_lib
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
    type(coord) :: optmol    !> quenched structure
    integer :: iter,iostatus
    real(wp) :: etot
    real(wp),allocatable :: grd(:,:)
    logical :: accept

   
!>--- Add input energy to Markov chain
    bh%emin = mol%energy    


    allocate(grd(3,mol%nat), source=0.0_wp)

!=======================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
    MonteCarlo : do iter = 1,bh%maxsteps
       bh%iteration = iter
     
!>--- Take the step (mol --> tmpmol)
      call takestep(mol,calc,bh,tmpmol)

!>--- Quench it (tmpmol --> optmol)
      call optimize_geometry(tmpmol,optmol,calc,etot,grd, &
      &                      .false.,.false.,iostatus)

!>--- Accept/reject
    if(iostatus == 0)then  !> successfull optimization

      if(debug)then
        write(*,*) 'Quench energy',etot,' Markov chain energy',bh%emin 
      endif 

      accept = mcaccept(optmol,bh)
      if( accept )then
        call axis(optmol%nat,optmol%at,optmol%xyz)

        if(debug) write(*,*) "accepted quench"
      !> check duplicates here

      else
        if(debug) write(*,*) 'Quench does not fulfill MC criterion'
        cycle MonteCarlo
      endif
    else
       if(debug) write(*,*) "Quench failed"
       cycle MonteCarlo
    endif 

!>--- Update structures
      mol = optmol
      call bh%add(mol)

    end do MonteCarlo
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=======================================================================================!

!>--- Stats


    deallocate(grd)
  end subroutine mc

!=========================================================================================!

  function mcaccept(mol,bh) result(accept)
!**************************************
!* The regular MC acceptance condition
!**************************************
    implicit none
    logical :: accept
    type(coord),intent(in) :: mol
    type(bh_class),intent(in) :: bh
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
