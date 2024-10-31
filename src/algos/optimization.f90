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

subroutine crest_optimization(env,tim)
!***********************************************
!* subroutine crest_optimization
!* This routine implements a standalone runtype
!* to perform geometry optimization for the 
!* specified input file (read from env%ref)
!***********************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)

  character(len=80) :: atmp
  character(len=*),parameter :: partial = '∂E/∂'
!========================================================================================!
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()
  call tim%start(14,'Geometry optimization')
!========================================================================================!
  call env%ref%to(mol)
  write (stdout,*)
  call smallhead('Input structure:')
  call mol%append(stdout)
  write (stdout,*)

!========================================================================================!

  allocate (grad(3,mol%nat),source=0.0_wp)
  calc = env%calc
!>--- check if we have any calculation settings allocated
  if (calc%ncalculations < 1) then
    write (stdout,*) 'no calculations allocated'
    return
  else
    call calc%info( stdout )
  end if
  write(stdout,'(a)') repeat('-',80)

!>-- geometry optimization
  pr = .true. !> stdout printout
  wr = .true. !> write crestopt.log
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  if (io == 0) then
!>--- print out the results
    write (stdout,*) 'SUCCESS!'
    write (stdout,*) 'geometry successfully optimized!'
    write (stdout,*)
    write(stdout,'(a)') repeat('-',80)

    write (stdout,*)
    write (stdout,'(1x,a)') 'FINAL CALCULATION SUMMARY'
    write (stdout,'(1x,a)') '========================='
    call calculation_summary(calc,mol,energy,grad,molnew)

    write (stdout,*)
    write (stdout,'(a)') '> Optimized geometry written to crestopt.xyz'
    gnorm = norm2(grad)
    write (atmp,'(1x,"Etot=",f16.10,1x,"g norm=",f12.8)') energy,gnorm
    molnew%comment = trim(atmp)

    open (newunit=ich,file='crestopt.xyz')
    call molnew%append(ich)
    close (ich)
  else
    write (stdout,*) 'geometry optimization FAILED!'
    env%iostatus_meta = status_failed
  endif

  write (stdout,*)
  write (stdout,'(a)') repeat('=',42)
  write (stdout,'(1x,a,f20.10,a)') 'TOTAL ENERGY   ',energy,' Eh'
  write (stdout,'(1x,a,f20.10,a)') 'GRADIENT NORM  ',norm2(grad),' Eh/a0'
  write (stdout,'(a)') repeat('=',42)

   if(io /= 0)then
    write (stdout,*) 'WARNING: geometry optimization FAILED!'
   endif

  deallocate (grad)
!========================================================================================!
  call tim%stop(14)

!========================================================================================!
!>--- append numerical hessian calculation
  if( io == 0 .and. env%crest_ohess )then
    call env%ref%load(molnew)      !> load the optimized geometry
    call crest_numhess(env,tim) !> run the numerical hessian
  endif

!========================================================================================!
  return
end subroutine crest_optimization

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
subroutine crest_ensemble_optimization(env,tim)
!***********************************************
!* subroutine crest_ensemble_optimization
!* This routine implements a standalone runtype
!* to perform geometry optimizations along an
!* ensemble or trajectory file.
!***********************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use parallel_interface
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,c,T,Tn
  logical :: pr,wr,ex
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)

  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  character(len=10),allocatable :: comments(:)
  character(len=80) :: atmp
  real(wp) :: percent
  character(len=52) :: bar
!========================================================================================!
  write (*,*)
!>--- check for the ensemble file
  inquire (file=env%ensemblename,exist=ex)
  if (ex) then
    ensnam = env%ensemblename
  else
    write (stdout,*) '**ERROR** no ensemble file provided.'
    env%iostatus_meta = status_input
    return
  end if

!>--- start the timer
  call tim%start(14,'Ensemble optimization')

!>---- read the input ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1)then
    write (stdout,*) '**ERROR** empty ensemble file.'
    env%iostatus_meta = status_input
    return
  endif
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_oloop requires coordinates in Bohrs
  xyz = xyz / bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

!>--- set OMP parallelization
  call new_ompautoset(env,'auto',nall,T,Tn)

!========================================================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,14x,"│")') "ENSEMBLE OPTIMIZATION"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)
  write (stdout,'(1x,a,i0,a,1x,a)') 'Optimizing all ',nall,' structures of file',trim(ensnam)
!>--- call the loop
  call crest_oloop(env,nat,nall,at,xyz,eread,.true.)

!========================================================================================!
!>--- output
  write(stdout,'(/,a,a,a)') 'Rewriting ',ensemblefile,' in the correct order'// &
  & ' (failed optimizations are assigned an energy of +1.0)'
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Back to Angstroem
  xyz = xyz * bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
  allocate(comments(nall))
  do i=1,nall
     comments(i) = ''
     if(eread(i) > 0.0_wp) comments(i) = '!failed'
  enddo
  call wrensemble(ensemblefile,nat,nall,at,xyz,eread,comments)

  deallocate (eread,at,xyz)
  write(stdout,'(/,a,a,a)') 'Optimized ensemble written to <',ensemblefile,'>'

!========================================================================================!
!>--- (optional) refinement step
  if (allocated(env%refine_queue)) then
    write(stdout,*)
    call crest_refine(env,ensemblefile,ensemblefile//'.refine')
    write(stdout,'(/,a,a,a)') 'Refined ensemble written to <',ensemblefile,'.refine>'
  endif 

!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_ensemble_optimization

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
subroutine crest_ensemble_screening(env,tim)
!****************************************************
!* subroutine crest_ensemble_screening
!* This routine implements a standalone runtype
!* to perform geometry optimizations along an
!* ensemble in a multilevel step and sort in between
!****************************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use iomod 
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,c,T,Tn
  logical :: pr,wr,ex
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)

  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  character(len=80) :: atmp
  real(wp) :: percent
  character(len=52) :: bar
  logical :: multilevel(6)
!========================================================================================!
  write (*,*)
!>--- check for the ensemble file
  inquire (file=env%ensemblename,exist=ex)
  if (ex) then
    ensnam = env%ensemblename
  else
    write (stdout,*) '**ERROR** no ensemble file provided.'
    env%iostatus_meta = status_input
    return
  end if

!>--- start the timer
  call tim%start(14,'Ensemble screening')

!>---- read the input ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1)then
    write (stdout,*) '**ERROR** empty ensemble file.'
    env%iostatus_meta = status_input
    return
  endif
!>--- set OMP parallelization
  call new_ompautoset(env,'auto',nall,T,Tn)

!========================================================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",48("━"),"┑")')
  write (stdout,'(10x,"│",15x,a,15x,"│")') "ENSEMBLE SCREENING"
  write (stdout,'(10x,"┕",48("━"),"┙")')
  write (stdout,*)
  write (stdout,'(1x,''Multilevel optimization and structure screening.'')')
  write (stdout,*)
  write (stdout,'(1x,a,a)') 'Input file: ','<'//trim(ensnam)//'>'
  write (stdout,'(1x,a,i0,a)') 'Containing ',nall,' structures.'

!>--- call the loop
  call rmrfw('crest_rotamers_')
  call optlev_to_multilev(3.0d0,multilevel)
  call crest_multilevel_oloop(env,ensnam,multilevel)
  if(env%iostatus_meta .ne. 0 ) return

!>--- printout
  call catdel('cregen.out.tmp')
  write (stdout,'(/,1x,a,1x,a)') 'Final ensemble on file','<'//trim(ensemblefile)//'>'

  call rename(conformerfile,trim(ensemblefile))

!>--- clean up
  call screen_cleanup


!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_ensemble_screening

