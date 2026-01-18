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

subroutine crest_setup_alkylize(env)
  use crest_parameters
  use crest_data
  use strucrd
  use molbuilder_classify
  implicit none
  type(systemdata),intent(inout) :: env
  type(coord_classify) :: molc
  type(coord) :: mol

  call env%ref%to(mol)

  call underline("Analyzing Input Structure")

  call setup_classify(mol,molc)
  call functional_group_classify(molc)
  if(molc%nfuncs == 0)then
    write(stdout,'(a)') 'no relevant substructures found'
  else 
    write(stdout,'(a)') 'Found the following substructure parts'
    call molc%print_funcgroups(stdout)
  endif

  !> TODO the actual checkup

  write(stdout,*)
end subroutine crest_setup_alkylize
