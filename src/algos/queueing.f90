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

subroutine crest_queue_setup(env,iterate)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use construct_list
  use construct_mod
  implicit none
  type(systemdata),intent(inout) :: env
  logical,intent(out) :: iterate

  integer :: splitlayers
  integer :: ii,jj,nn
  type(coord),pointer :: reference_mol
  type(coord),target :: mol
  integer,allocatable :: splitatms(:)
  integer :: parentlayer,parentnode
  character(len=1024) :: thispath

  iterate = .false.

  if (allocated(env%splitqueue)) then

    !> check for incompatible runtypes (or rather, whitelist a few)
    if (.not.any(env%crestver == [crest_imtd,crest_imtd2,crest_sp, &
      & crest_optimize,crest_moldyn,crest_rigcon,crest_trialopt,crest_bh,crest_test])) then
      write (stdout,'(a)') '** ERROR ** Selected CREST runtype incompatible with substructure builder'
      call creststop(status_config)
    end if
    if (allocated(env%ONIOM_input).or.allocated(env%ONIOM_toml)) then
      write (stdout,'(a)') '** ERROR ** ONIOM incompatible with substructure builder'
      call creststop(status_config)
    end if

    !> if the program sees no problem, set the global boolean
    env%substructure_queue = .true.
    splitlayers = size(env%splitqueue,1)

    !> start constructing the splitheap
    env%splitheap%nlayer = splitlayers
    allocate (env%splitheap%layer(splitlayers))
    associate (heap => env%splitheap,layer => env%splitheap%layer)

      do ii = 1,heap%nlayer
        layer%id = ii

        nn = env%splitqueue(ii)%natms
        allocate (splitatms(nn))
        splitatms(:) = env%splitqueue(ii)%atms(:)

        if (ii == 1) then
          call env%ref%to(mol)
          reference_mol => mol
        else

          call pick_parent(heap,ii,splitatms,parentlayer,parentnode)
          if (parentlayer == 0) then
            call env%ref%to(mol)
            reference_mol => mol
          else
            mol = heap%layer(parentlayer)%node(parentnode)
            reference_mol => mol
          end if
        end if
        call split(reference_mol,splitatms,layer(ii)%node,layer(ii)%alignmap, &
          & ncap=layer(ii)%ncapped,position_mapping=layer(ii)%position_mapping)
        deallocate (splitatms)
        layer(ii)%nnodes = size(layer(ii)%node,1)
        call heap%map_origins_for_layer(ii)
      end do

      call heap%setup_queue()
      !write (*,*) 'endpoints selected:'
      !do ii = 1,heap%nqueue
      !  write (*,*) 'layer and node',heap%queue(ii)%layer,heap%queue(ii)%node
      !end do
      call getcwd(thispath)
      heap%origindir = trim(thispath)
    end associate
    iterate = .true.
  end if

  return

contains
  subroutine pick_parent(heap,current_layer,splitatms,parentlayer,parentnode)
    use construct_list
    implicit none
    type(construct_heap),intent(inout) :: heap
    integer,intent(inout) :: splitatms(:)
    integer,intent(in) :: current_layer
    integer,intent(out) :: parentlayer,parentnode
    integer :: ii,jj,kk,prev_layer
    logical :: matching

    parentlayer = 0
    parentnode = 0
    if (current_layer .eq. 1) return

    !> iterate through the previous layer and check which node
    !> contains all the split atoms
    prev_layer = current_layer-1
    LAYITER: do while (prev_layer >= 1)
      do ii = 1,heap%layer(prev_layer)%nnodes
        matching = .true.
        do jj = 1,size(splitatms,1)
          matching = matching.and.any(heap%layer(prev_layer)%origin(ii)%map(:) .eq. splitatms(jj))
        end do
        if (matching) then
          parentlayer = prev_layer
          parentnode = ii
          !> on the first match, exit
          exit LAYITER
        end if
      end do
      !> if no matching parent node was found, try again in one layer further up
      if (parentnode == 0) prev_layer = prev_layer-1
    end do LAYITER

    !> IMPORTANT; we need to update the splitatms with the correctly mapped indices
    if (parentnode .ne. 0) then
      do ii = 1,size(splitatms,1)
        jj = splitatms(ii)
        call heap%find_current_position(jj,parentlayer,parentnode,kk)
        splitatms(ii) = kk
      end do
      !> we also map the current node as a child node of the selected parent
      if (.not.allocated(heap%layer(parentlayer)%childlayer)) then
        ii = heap%layer(parentlayer)%nnodes
        allocate (heap%layer(parentlayer)%childlayer(ii),source=0)
      end if
      heap%layer(parentlayer)%childlayer(parentnode) = current_layer
    end if

  end subroutine pick_parent
end subroutine crest_queue_setup

subroutine crest_queue_iter(env,iterate)
  use crest_parameters
  use crest_data
  use strucrd
  use iomod
  implicit none
  type(systemdata),intent(inout) :: env
  logical,intent(out) :: iterate
  integer :: ii,jj,kk,io
  type(coord) :: mol
  character(len=10) :: atmp
  character(len=*),parameter :: dirname = 'crest_queue_'

  iterate = .false.

  if (allocated(env%splitqueue).and.env%splitheap%nqueue > 0) then
    call chdir(env%splitheap%origindir)
    ii = env%queue_iter+1
    env%queue_iter = ii

    associate (queue => env%splitheap%queue(ii))

      jj = queue%layer
      kk = queue%node

      !> create a dedicated work directory
      write(atmp,'(i0)') ii
      queue%workdir = dirname//trim(atmp)
      io = makedir(queue%workdir)
      call chdir(queue%workdir)

      mol = env%splitheap%layer(jj)%node(kk)
      call env%ref%load(mol)

    end associate
    if (ii < env%splitheap%nqueue) then
      iterate = .true.
    end if
  end if
end subroutine crest_queue_iter
