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

subroutine crest_basinhopping(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use bh_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc
  type(bh_class) :: bh

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)

  character(len=80) :: atmp
  character(len=*),parameter :: trjf='crest_quenched.xyz'
!========================================================================================!
  write(stdout,*)
  !call system('figlet dynamics')
  write(stdout,*) " ____            _       _   _                   _             "
  write(stdout,*) "| __ )  __ _ ___(_)_ __ | | | | ___  _ __  _ __ (_)_ __   __ _ "
  write(stdout,*) "|  _ \ / _` / __| | '_ \| |_| |/ _ \| '_ \| '_ \| | '_ \ / _` |"
  write(stdout,*) "| |_) | (_| \__ \ | | | |  _  | (_) | |_) | |_) | | | | | (_| |"
  write(stdout,*) "|____/ \__,_|___/_|_| |_|_| |_|\___/| .__/| .__/|_|_| |_|\__, |"
  write(stdout,*) "                                    |_|   |_|            |___/ "
  write(stdout,*) ""
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()

  calc = env%calc 
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!
  pr = .true.
!>--- print calculation info
  call calc%info( stdout )

!>--- singlepoint of input structure
  allocate(grad(3,mol%nat), source=0.0_wp)
  call engrad(mol,calc,energy,grad,io)
  mol%energy = energy  !> we need this to start the Markov-chain

!>--- actual basin hopping
  call bh%init(300.0_wp,50,20)
  bh%stepsize(1) = 0.75_wp

  call tim%start(14,'Basin-Hopping (BH)')

  call mc(calc,mol,bh)


  open(newunit=ich,file=trjf)
    do i=1,bh%saved
      call bh%structures(i)%append(ich) 
    enddo
  close(ich)

  if (io == 0) then
    write (stdout,*) 'BH run completed successfully'
    write (stdout,*) 'Successfull quenches written to ',trjf
  else
    write (stdout,*) 'WARNING: BH run terminated ABNORMALLY'
    env%iostatus_meta = status_failed
  end if
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_basinhopping
