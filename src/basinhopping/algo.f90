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

subroutine crest_basinhopping(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use cregen_interface,only:unionizeEnsembles
  use optimize_module
  use bh_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn,mciter
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc
  type(bh_class) :: bh

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  integer :: nall
  type(coord),allocatable :: structuredump(:)

  character(len=80) :: atmp
  character(len=*),parameter :: trjf = 'crest_quenched.xyz'
!========================================================================================!
  write (stdout,*)
  write (stdout,*) " ____            _       _   _                   _             "
  write (stdout,*) "| __ )  __ _ ___(_)_ __ | | | | ___  _ __  _ __ (_)_ __   __ _ "
  write (stdout,*) "|  _ \ / _` / __| | '_ \| |_| |/ _ \| '_ \| '_ \| | '_ \ / _` |"
  write (stdout,*) "| |_) | (_| \__ \ | | | |  _  | (_) | |_) | |_) | | | | | (_| |"
  write (stdout,*) "|____/ \__,_|___/_|_| |_|_| |_|\___/| .__/| .__/|_|_| |_|\__, |"
  write (stdout,*) "                                    |_|   |_|            |___/ "
  write (stdout,*) ""
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()

  calc = env%calc
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!
!>--- print calculation info
  call calc%info(stdout)
  write (stdout,'(a)') '> Geometry optimization settings:'
  call print_opt_data(calc,stdout,natoms=mol%nat,tag=' : ')
  write (stdout,*)

!>--- singlepoint of input structure
  allocate (grad(3,mol%nat),source=0.0_wp)
  call engrad(mol,calc,energy,grad,io)
  mol%energy = energy  !> we need this to start the Markov-chain

!>--- actual basin hopping
  if (allocated(env%bh_ref)) then
    bh = env%bh_ref
    call bh%init()
  else
    call bh%init(300.0_wp,200,20)
    bh%stepsize(1) = 1.0_wp
  end if

  call tim%start(14,'Basin-Hopping (BH)')
  nall = 0
  do mciter = 1,bh%maxiter
    if (bh%maxiter > 1) call printiter2(mciter)
    call bh%newiter()
    call mc(calc,mol,bh,verbosity=2)

    write (stdout,'(a)') 'New structures will be appended to memory ...'
    call unionizeEnsembles(nall,structuredump,bh%saved,bh%structures, &
    &                      ethr=bh%ethr,rthr=bh%rthr)
    write (stdout,'(a,i0,a)') 'Currently ',nall,' structures saved!'
  end do

!>--- dump saved minima
  open (newunit=ich,file=trjf)
  call wrensemble(ich,nall,structuredump)
  close (ich)

  if (io == 0) then
    write (stdout,*)
    write (stdout,*) 'BH run completed successfully'
    write (stdout,*) 'Successfull quenches written to ',trjf
  else
    write (stdout,*) 'WARNING: BH run terminated ABNORMALLY'
    env%iostatus_meta = status_failed
  end if
!========================================================================================!
  call tim%stop(14)

  if (allocated(structuredump)) deallocate (structuredump)
  return
end subroutine crest_basinhopping
