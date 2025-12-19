!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> Implementation of whatever, for testing implementations
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> Input/Output:
!>  env  -  crest's systemdata object
!>  tim  -  timer object
!>-----------------------------------------------
subroutine crest_playground(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc
  real(wp) :: accuracy,etemp

  integer :: V,maxgen
  integer,allocatable :: A(:,:)
  logical,allocatable :: rings(:,:)
  integer,allocatable :: tmp(:)
  logical :: connected,fail

  real(wp) :: energy
  real(wp),allocatable :: grad(:,:),geo(:,:),csv(:,:),q(:)

!========================================================================================!
  call tim%start(14,'Test implementation')
!========================================================================================!
  !call system('figlet welcome')
  write (*,*) "              _                          "
  write (*,*) "__      _____| | ___ ___  _ __ ___   ___ "
  write (*,*) "\ \ /\ / / _ \ |/ __/ _ \| '_ ` _ \ / _ \"
  write (*,*) " \ V  V /  __/ | (_| (_) | | | | | |  __/"
  write (*,*) "  \_/\_/ \___|_|\___\___/|_| |_| |_|\___|"
  write (*,*)
!========================================================================================!
!  call env%ref%to(mol)
!  write (*,*)
!  write (*,*) 'Input structure:'
!  call mol%append(stdout)
!  write (*,*)
!!========================================================================================!
!
!  allocate (grad(3,mol%nat),source=0.0_wp)
!  call env2calc(env,calc,mol)
!  calc%calcs(1)%rdwbo = .true.
!  call calc%info(stdout)
!
!  call engrad(mol,calc,energy,grad,io)
!  call calculation_summary(calc,mol,energy,grad)
!========================================================================================!
  block
    use construct_mod
    type(coord) :: base,side,new,newnew
    type(coord),allocatable :: splitlist(:)
    integer,allocatable :: alignmap(:,:)

    !call base%open("base.xyz")
    !call side%open("side.xyz")

    open (newunit=ich,file='molbuilder.xyz')
    !call base%append(ich)
    !call side%append(ich)

    !allocate (alignmap(3,2),source=0)

    !alignmap(1:3,1) = [9,7,8]
    !alignmap(1:3,2) = [3,1,2]
    !call attach(base,side,alignmap,new)
    !call new%append(ich)


    call new%open("struc.xyz")
    !call split(new, [8,9],base,side)
    call split(new, [6,7,8],splitlist,alignmap)

    do i=1,size(splitlist,1)
      call splitlist(i)%append(ich)
    enddo

    call attach(splitlist(1), splitlist(2), alignmap,newnew)
 
    call newnew%append(ich)
    !call base%append(ich)
    !call side%append(ich)
    close (ich) 
    
  end block

!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_playground
