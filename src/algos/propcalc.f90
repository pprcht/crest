!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!> Skeleton stub for property calculations on an ensemble.
!> Each mode needs to be implemented via modern algo routines.
!> The legacy implementation (system-call/I/O based) has been removed.
!> See git history for reference.

subroutine propcalc(iname,imode,env,tim)
  use crest_parameters
  use crest_data
  use cregen_interface
  implicit none
  character(len=*),intent(in) :: iname
  integer,intent(in) :: imode
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout) :: tim

  select case (imode)
  case (p_prop_hess)
    !> TODO: Hessian calculations for all conformers (was: xtb --hess)
  case (p_prop_autoir)
    !> TODO: IR spectrum averaging over populated conformers (was: autoir + xtb --ohess)
  case (p_prop_ohess)
    !> TODO: Optimization + Hessian for all conformers (was: xtb --ohess)
  case (p_prop_gsolv)
    !> TODO: Free energy in solvation, 2-step (was: xtb --sp + xtb --ohess)
  case (p_prop_reopt)
    !> TODO: Vtight reoptimization for all conformers (was: xtb --opt vtight)
  case (p_prop_multilevel:p_prop_multilevel+9)
    !> Post-search re-optimization of the conformer ensemble at the higher level.
    !> Input iname is typically crest_rotamers.xyz; output is crest_reopt.xyz.
    call crest_multilevel_reopt(iname,env,tim)
  case (p_prop_dipole)
    !> TODO: Singlepoint + dipole extraction (was: xtb --sp, grep molecular dipole)
  case (p_prop_rerank)
    !> TODO: Singlepoint + reranking (was: xtb --sp + newcregen)
  case default
    write (stdout,'(a,i0,a)') 'propcalc: mode ',imode,' not yet implemented'
  end select

end subroutine propcalc

!========================================================================================!

subroutine crest_multilevel_reopt(iname,env,tim)
!*******************************************************************
!* Read the ensemble iname, optimize all structures using the
!* calculator tagged with refine_lvl == refine%post_opt (set by the
!* A@B hybrid keyword), sort via CREGEN, and write crest_reopt.xyz.
!*
!* The refine_stage mechanism in calculator.F90 is used to activate
!* only the post-search calculator during crest_oloop.
!*
!* Input:
!*   iname  - path to input ensemble (e.g. crest_rotamers.xyz)
!* Output:
!*   crest_reopt.xyz (sorted conformer ensemble at the higher level)
!*******************************************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use parallel_interface
  use cregen_interface
  use iomod,only:drawbox,catdel
  implicit none
  character(len=*),intent(in) :: iname
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout) :: tim
  integer :: nat,nall,T,Tn,old_stage
  real(wp),allocatable :: xyz(:,:,:),eread(:)
  integer,allocatable  :: at(:)
  character(len=*),parameter :: outname = 'crest_reopt.xyz'
  logical :: ex

  inquire(file=iname,exist=ex)
  if (.not.ex) then
    write(stdout,'(a,a,a)') '**WARNING** ',trim(iname),' not found, skipping multilevel reopt'
    return
  end if

  call tim%start(16,'Multilevel reopt')

  call rdensembleparam(iname,nat,nall)
  if (nall < 1) then
    write(stdout,*) '**WARNING** empty ensemble, skipping multilevel reopt'
    call tim%stop(16)
    return
  end if
  allocate(xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(iname,nat,nall,at,xyz,eread)
! ── crest_oloop requires Bohr ────────────────────────────────────
  xyz = xyz/bohr

  call new_ompautoset(env,'auto',nall,T,Tn)

  write(stdout,*)
  call drawbox(stdout,'MULTILEVEL ENSEMBLE REOPT',charset=7,width=51,ltab=10)
  write(stdout,'(1x,a,i0,a,1x,a)') &
    & 'Re-optimizing ',nall,' structures of file ',trim(iname)

! ── activate only the post-search calculator ─────────────────────
  old_stage = env%calc%refine_stage
  env%calc%refine_stage = refine%post_opt

  call crest_oloop(env,nat,nall,at,xyz,eread,.true.)

  env%calc%refine_stage = old_stage

! ── back to Angstrom, write output ───────────────────────────────
  xyz = xyz*bohr
  call wrensemble(outname,nat,nall,at,xyz,eread)

  write(stdout,'(/,a,a,a)') 'Re-optimized ensemble written to <',outname,'>'

! ── sort via CREGEN ──────────────────────────────────────────────
  call newcregen(env,0,outname)
  call catdel('cregen.out.tmp')

  deallocate(xyz,at,eread)
  call tim%stop(16)
end subroutine crest_multilevel_reopt
