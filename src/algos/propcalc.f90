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
    !> TODO: Multilevel/hybrid reoptimization of entire CRE, e.g. GFN2@GFF
    !>       (was: xtb --opt vtight with gfnver2)
  case (p_prop_dipole)
    !> TODO: Singlepoint + dipole extraction (was: xtb --sp, grep molecular dipole)
  case (p_prop_rerank)
    !> TODO: Singlepoint + reranking (was: xtb --sp + newcregen)
  case default
    write (stdout,'(a,i0,a)') 'propcalc: mode ',imode,' not yet implemented'
  end select

end subroutine propcalc
