!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021-2025 Christoph Plett, Philipp Pracht
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

module qcg_utils
  use crest_parameters,only:stdout,wp
  use crest_data
  use iomod
  implicit none
  public

contains

!========================================================================================!
!> Convert given QCG coordinate files into (TM format)
!> Write "solute" and "solvent" coordinate files
!========================================================================================!
  subroutine inputcoords_qcg(env,solute,solvent)
    use crest_parameters
    use crest_data
    use strucrd
    use zdata
    use iomod
    implicit none

    type(systemdata),intent(inout) :: env
    type(zmolecule),intent(out) :: solute,solvent

    logical :: ex11,ex21,solu,solv
    type(coord) :: mol
    type(zmolecule) :: zmol,zmol1
    integer :: i

!--------------------Checking for input-------------!

    !Solute
    inquire (file=env%solu_file,exist=ex11)
    inquire (file='solute',exist=solu)
    if (solu) call copy('solute','solute.old') !Backup solute file
    if ((.not.ex11).and.(.not.solu)) then
      error stop 'No (valid) solute file! exit.'
    else if ((.not.ex11).and.(solu)) then
      env%solu_file = 'solute'
    end if

    !Solvent
    inquire (file=env%solv_file,exist=ex21)
    inquire (file='solvent',exist=solv)
    if (solu) call copy('solvent','solvent.old') !Backup solvent file
    if ((.not.ex21).and.(.not.solv)) then
      error stop 'No (valid) solvent file! exit.'
    else if ((.not.ex11).and.(solu)) then
      env%solu_file = 'solvent'
    end if

!---------------Handling solute---------------------!
    call mol%open(env%solu_file)
    call mol%write('solute')
    solute%nat = mol%nat
    solute%at = mol%at
    solute%xyz = mol%xyz

    !--- if the input was a SDF file, special handling
    env%sdfformat = .false.
    call checkcoordtype(env%solu_file,i)
    if (i == 31.or.i == 32) then
      !Add sdf stuff here, if somebody needs it
    end if
    !--- Add as ref structure in env
    call env%ref%load(mol)
    call mol%deallocate() 
!---------------Handling solvent---------------------!

    call mol%open(env%solv_file)
    call mol%write('solvent')
    solvent%nat = mol%nat
    solvent%at = mol%at
    solvent%xyz = mol%xyz
    call mol%deallocate()

    !--- if the input was a SDF file, special handling
    env%sdfformat = .false.
    call checkcoordtype(env%solv_file,i)
    if (i == 31.or.i == 32) then
      !Add sdf stuff here, if somebody needs it
    end if

    return
  end subroutine inputcoords_qcg

!==============================================================================!

  subroutine write_reference(env,solu,clus)
    use crest_data
    use zdata,only:zmolecule
    use iomod
    use strucrd
    implicit none
    type(systemdata):: env    ! MAIN STORAGE OS SYSTEM DATA
    type(zmolecule)            :: solu,clus
    type(zmolecule)            :: ref_mol,ref_clus
    ref_mol = solu
    call rdcoord(env%solu_file,ref_mol%nat,ref_mol%at,ref_mol%xyz) !original solute coordinates
    call remove(env%fixfile)
    ref_clus = clus
    ref_clus%xyz(1:3,1:solu%nat) = solu%xyz
    call wrc0(env%fixfile,ref_clus%nat,ref_clus%at,ref_clus%xyz)
  end subroutine write_reference

!=============================================================================!

  subroutine aver(pr,env,runs,e_tot,S,H,G,sasa,a_present,a_tot)
    use crest_parameters
    use crest_data

    implicit none
!---- Dummy
    type(systemdata),intent(in)   :: env
    integer,intent(in)            :: runs
    real(wp),intent(inout)        :: e_tot
    real(wp),intent(in),optional  :: a_tot
    real(wp),intent(out)          :: S
    real(wp),intent(out)          :: H
    real(wp),intent(out)          :: G
    real(wp),intent(out)          :: sasa
!---- Stack
    logical,intent(in)            :: pr,a_present
    integer                       :: j,jmin
    real(wp)                      :: A
    real(wp)                      :: e0
    real(wp),allocatable          :: de(:)
    real(wp),allocatable          :: p(:)
    real(wp)                      :: pmax
    real(wp)                      :: eav
    real(wp)                      :: area
    real(wp)                      :: beta
    real(wp)                      :: temp
    integer                       :: ich48
    dimension e_tot(runs)
    dimension a_tot(runs)

    temp = env%tboltz
    allocate (de(runs),source=0.0d0)
    allocate (p(runs),source=0.0d0)

    beta = 1./(temp*8.314510/4.184/1000.+1.d-14)
    e0 = e_tot(1)
    de(1:runs) = (e_tot(1:runs)-e0)
    call qcg_boltz(env,runs,de,p)

    A = 0
    eav = 0
    pmax = 0
    area = 0
    do j = 1,runs
      A = A+p(j)*log(p(j)+1.d-12)
      eav = eav+p(j)*e_tot(j)
      if (p(j) .gt. pmax) then
        pmax = p(j)
        jmin = j
      end if
      if (a_present) area = area+p(j)*a_tot(j)
    end do
    sasa = area
    S = (1./beta)*A
    H = eav
    G = eav+S
    if (pr) then
      open (newunit=ich48,file='population.dat')
      write (ich48,'(2x, ''cluster'',2x,''E_norm [Eh]'',2x, ''De [kcal]'', 4x, ''p'')')
      do j = 1,runs
        if (j .lt. 10) then
          write (ich48,'(5x,i0,3x,f11.6,5x,f6.4,3x,f6.4)') j,e_tot(j)/autokcal,de(j),p(j)
        else
          write (ich48,'(5x,i0,2x,f11.6,5x,f6.4,3x,f6.4)') j,e_tot(j)/autokcal,de(j),p(j)
        end if
      end do
      write (ich48,*)
      write (ich48,'(''Ensemble free energy [Eh]:'', f20.10)') G/autokcal
      close (ich48)
    end if

    deallocate (de,p)

  end subroutine aver

 !=============================================================================! 
end module qcg_utils
