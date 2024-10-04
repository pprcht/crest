!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2024 Philipp Pracht
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
!> Implementation for standalone sorting
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
!> Input/Output:
!>  env  -  crest's systemdata object
!>  tim  -  timer object
!>-----------------------------------------------
subroutine crest_sort(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd 
  use cregen_interface
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich 
  logical :: pr,wr
!========================================================================================!
  integer :: nall
  type(coord),allocatable :: structures(:)
  integer,allocatable :: groups(:) 


!========================================================================================!


  write(stdout,'(a,a,a)',advance='no') '> Read ensemble ',trim(env%ensemblename),' ... '
  flush(stdout)
  call rdensemble(env%ensemblename,nall,structures)
  allocate(groups(nall), source=0)
  write(stdout,'(i0,a)') nall,' structures!' 
  write(stdout,*)

!========================================================================================!
  call tim%start(11,'Sorting') 

  select case(env%sortmode)

  case('isort')
!>--- Assigning structures to conformers based on RTHR,with canonical atom IDs
    call underline('Assigning conformers based on iRMSD and RTHR')
    call cregen_irmsd_sort(env,nall,structures,groups,allcanon=.true.,printlvl=2)    


  case('isort_noid')
!>--- Assigning structures to conformers based on RTHR, WITHOUT canonical atom IDs
    call underline('Assigning conformers based on iRMSD and RTHR')
    call cregen_irmsd_sort(env,nall,structures,groups,allcanon=.false.,printlvl=2)    


  case('all','allpair')
!>--- all unique pairs of the ensemble (only suitable for small ensembles)
    call underline('Running all unique pair RMSDs incl. atom permutation')
    call cregen_irmsd_all(nall,structures,2)

  case default
!>--- all unique pairs of the ensemble (only suitable for small ensembles)
    call cregen_irmsd_all(nall,structures,2)
  end select

!========================================================================================!
  call tim%stop(11)
  return
end subroutine crest_sort
