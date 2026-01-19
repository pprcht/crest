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

  integer :: ii,jj,kk
  integer :: splt(3)

  call env%ref%to(mol)

  call underline("Analyzing Input Structure")

  call setup_classify(mol,molc)
  call functional_group_classify(molc)
  if (molc%nfuncs == 0) then
    write (stdout,'(a)') 'no relevant substructures found'
    return
  else
    write (stdout,'(a)') 'Found the following substructure parts'
    call molc%print_funcgroups(stdout)
  end if

  do ii = 1,molc%nfuncs
    associate (func => molc%funcgroups(ii))
      if (trim(func%name) == 'alkyl') then

        !> only for propane or longer
        if (func%natms > 6) then
          write (stdout,'(a)') 'selected alkyl group for fragment dispatching'
          splt(:) = 0
          splt(1) = func%attached_to
          kk = 1
          do while (kk < 3)
            do jj = 1,molc%nat
              if (molc%Ah(jj,splt(kk)) == 1.and. &
                & .not.any(splt(:) .eq. jj).and. &
                & any(func%ids(:) .eq. jj)) then
                kk = kk+1
                splt(kk) = jj
              end if
            end do
          end do
          call env%addsplitqueue(splt)
        end if

      end if
    end associate
  end do

  write (stdout,*)
end subroutine crest_setup_alkylize

subroutine crest_proxy_nalkane(env,doreturn)
  use crest_parameters
  use crest_data
  use strucrd
  use molbuilder_classify
  use INTERNALS_mod
  implicit none
  type(systemdata),intent(inout) :: env
  logical,intent(out) :: doreturn
  type(coord) :: mol,newmol
  type(coord_classify) :: molc
  integer :: ii,jj
  integer,allocatable  :: na(:),nb(:),nc(:)
  real(wp),allocatable :: zmat(:,:)

  doreturn = .false.

  if (env%alkylize) then
    call env%ref%to(mol)
    call setup_classify(mol,molc)
    call functional_group_classify(molc)

    do ii = 1,molc%nfuncs
      if (molc%funcgroups(ii)%name == 'alkane') then
        write (stdout,'(a)') '> This substructure contains an n-alkane.'
        write (stdout,'(a)') '> SKIPPING sampling and writing linear structure.'
        doreturn = .true.

        !> ZMAT construction to make the molecule linear
        allocate (na(mol%nat),nb(mol%nat),nc(mol%nat),source=0)
        allocate (zmat(3,mol%nat),source=0.0_wp)
        call BETTER_XYZINT(mol%nat,mol%xyz,molc%A,na,nb,nc,zmat)

        !> setting internal CC dihedrals to trans-config
        do jj=1,mol%nat
           if(mol%at(jj) == 6 .and. mol%at(na(jj)) == 6 .and. &
             mol%at(nb(jj)) == 6 .and. mol%at(nc(jj)) == 6)then
              zmat(3,jj) = -pi
           endif
        enddo
        call smallhead('Internal coordinates:')
        call print_zmat(stdout,mol%nat,mol%at,zmat,na,nb,nc,.true.)
        call reconstruct_zmat_to_mol(mol%nat,mol%at,zmat,na,nb,nc,newmol)
        call newmol%write(conformerfile)

        exit
      end if
    end do

  end if

end subroutine crest_proxy_nalkane

