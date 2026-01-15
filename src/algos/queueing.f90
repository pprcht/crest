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

  iterate = .true.

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
      call getcwd(thispath)
      call env%ref%to(heap%originmol)
      heap%origindir = trim(thispath)
      heap%origincalc => env%calc
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
          matching = matching.and. &
            & any(heap%layer(prev_layer)%origin(ii)%map(:) .eq. splitatms(jj))
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

!=============================================================================!
!#############################################################################!
!=============================================================================!

subroutine crest_queue_iter(env,iterate)
  use crest_parameters
  use crest_data
  use strucrd
  use iomod
  use crest_calculator
  implicit none
  type(systemdata),intent(inout),target :: env
  logical,intent(out) :: iterate
  integer :: ii,jj,kk,io,nn,ll,lll,ati,atj
  type(coord) :: mol
  character(len=10) :: atmp
  character(len=*),parameter :: dirname = 'crest_queue_'

  iterate = .false.

  if (allocated(env%splitqueue).and.env%splitheap%nqueue > 0) then
!>--- important restoring to initial calc/dir
    env%calc => env%splitheap%origincalc
    call chdir(env%splitheap%origindir)

    !> next iter
    ii = env%queue_iter+1
    env%queue_iter = ii

    write (stdout,'(/,70("§"))')
    write (stdout,'(a,i0)') "§§§   QUEUE ITERATION ",ii
    write (stdout,'(70("§"))')

    jj = env%splitheap%queue(ii)%layer
    kk = env%splitheap%queue(ii)%node
    associate (heap => env%splitheap,queue => env%splitheap%queue(ii))

      !> create a dedicated work directory
      write (atmp,'(i0)') ii
      queue%workdir = dirname//trim(atmp)
      io = makedir(queue%workdir)
      call chdir(queue%workdir)
      write (stdout,'(a,t28,a,t30,a)') 'Queue work (sub-)directory',':', &
        & trim(queue%workdir)

      !> selecting output file depending on runtype
      select case (env%crestver)
      case (crest_imtd,crest_imtd2)
        queue%file = 'crest_conformers.xyz'
      case (crest_optimize)
        queue%file = 'crestopt.xyz'
      case (crest_moldyn)
        queue%file = 'crest_dynamics.trj.xyz'
      case (crest_bh)
        queue%file = 'crest_quenched.xyz'
      case default
        queue%file = 'struc.xyz'
      end select
      write (stdout,'(a,t28,a,t30,a)') 'Selected output file',':',queue%file

!>--- new calculator setup section and env update
      call queue%calc%copy(env%calc,ignore_constraints=.true.)
      !> for constraints we must be careful and map them to the new order
      call update_constraints_queue(heap,jj,kk,env%calc,queue%calc)

      call queue%calc%info(stdout)

      mol = env%splitheap%layer(jj)%node(kk)
      call env%ref%load(mol)
      call mol%write('coord')

      if (allocated(env%ref%wbo)) deallocate (env%ref%wbo)
      env%nat = mol%nat
      env%rednat = mol%nat

      env%calc => queue%calc

    end associate
    if (ii < env%splitheap%nqueue) then
      iterate = .true.
    end if

    write (stdout,*)
  end if

contains
  subroutine update_constraints_queue(heap,layer,node,refcalc,newcalc)
    use construct_list
    implicit none
    type(construct_heap) :: heap
    integer :: layer,node
    type(calcdata),intent(in) :: refcalc
    type(calcdata),intent(inout) :: newcalc
    integer :: nn,ll,lll,ati,atj,nn2
    type(constraint),allocatable :: cons(:)
    if (refcalc%nconstraints > 0) then
      nn = refcalc%nconstraints
      allocate (cons(nn))
      do ll = 1,nn
        call cons(ll)%copy(refcalc%cons(ll))
        do lll = 1,cons(ll)%n
          ati = cons(ll)%atms(lll)
          call heap%find_current_position(ati,layer,node,atj)
          cons(ll)%atms(lll) = atj !> overwrite with the current position
        end do
        if (any(cons(ll)%atms(:) .eq. 0)) then
          cons(ll)%active = .false.
        end if
      end do
      !> clean (active) constraints
      nn2 = 0
      do ll = 1,nn
        if (cons(ll)%active) nn2 = nn2+1
      end do
      if (nn2 > 0) then
        newcalc%nconstraints = nn2
        allocate (newcalc%cons(nn2))
        lll = 0
        do ll = 1,nn
          if (cons(ll)%active) then
            lll = lll+1
            call newcalc%cons(lll)%copy(cons(ll))
          end if
        end do
      end if
    end if
  end subroutine update_constraints_queue
end subroutine crest_queue_iter

!=============================================================================!
!#############################################################################!
!=============================================================================!

subroutine crest_queue_reconstruct(env,tim)
  use crest_parameters
  use crest_data
  use construct_list
  use strucrd
  implicit none
  type(systemdata),intent(inout) ::  env
  type(timer),intent(inout) :: tim
  type(coord) :: mol

  if (.not. (allocated(env%splitqueue).and.env%splitheap%nqueue > 0)) then
    return
  end if

  write (stdout,'(/,80("#"))')
  write (stdout,'(3("#"),t25,a,t78,3("#"))') 'QUEUE STRUCTURE RECONSTRUCTION'
  write (stdout,'(80("#"),/)')

  !> reset
  mol = env%splitheap%originmol
  call env%ref%load(mol)
  env%nat = mol%nat
  env%rednat = mol%nat
  env%calc => env%splitheap%origincalc
  call chdir(env%splitheap%origindir)

  call recusrive_construct(env,env%splitheap,1)

contains
  recursive subroutine recusrive_construct(env,heap,targetlayer)
    implicit none
    type(systemdata),intent(inout) :: env
    type(construct_heap),intent(inout) :: heap
    integer,intent(in) :: targetlayer

    integer :: ii,jj
    character(len=:),allocatable :: basefile,sidefile
    type(coord),allocatable :: structures_b(:)
    type(coord),allocatable :: structures_s(:)
    integer :: nall_b,nall_s,id_b,id_s
    logical :: ex

    character(len=*),parameter :: subdir_tmp = 'crest_queue_'
    character(len=:),allocatable :: subdirfile
    character(len=10) :: atmp

    associate (layer => heap%layer(targetlayer))
      if (layer%nnodes > 2) then
        write (stdout,'(a)') 'currently unhandled edge-case in layer reconstruction:'
        write (stdout,'(a,i0,a)') 'layer ',targetlayer,' was split in more than 2 structures'
        stop
      end if

      if (.not.allocated(layer%childlayer).or. &
        & all(layer%childlayer(:) .eq. 0)) then
        !> saveguard to not reconstruct layer multiple times
        if (layer%nmols > 0) return

        !> pick base an side-group files
        do ii = 1,heap%nqueue
          if (heap%queue(ii)%layer == targetlayer.and.heap%queue(ii)%node == 1) then
            basefile = heap%queue(ii)%file
            id_b = ii
          else if (heap%queue(ii)%layer == targetlayer.and.heap%queue(ii)%node == 2) then
            sidefile = heap%queue(ii)%file
            id_s = ii
          end if
        end do

        write (atmp,'(i0)') id_b
        subdirfile = subdir_tmp//trim(atmp)//'/'//basefile
        inquire (exist=ex,file=subdirfile)
        if (ex) then
          write (stdout,'(a,t28,a,t30,a)') 'Reading fragment(s) from',':',subdirfile
          call rdensemble(subdirfile,nall_b,structures_b)
        end if

        write (atmp,'(i0)') id_s
        subdirfile = subdir_tmp//trim(atmp)//'/'//sidefile
        inquire (exist=ex,file=subdirfile)
        if (ex) then
          write (stdout,'(a,t28,a,t30,a)') 'Reading fragment(s) from',':',subdirfile
          call rdensemble(subdirfile,nall_s,structures_s)
        end if

      else
        do ii = 1,layer%nnodes
          jj = layer%childlayer(ii)
          if (jj == 0) cycle
          call recusrive_construct(env,heap,jj)
        end do
      end if

    end associate
  end subroutine recusrive_construct

end subroutine crest_queue_reconstruct

