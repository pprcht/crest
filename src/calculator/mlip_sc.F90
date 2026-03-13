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

!> module mlip_sc
!> A module containing routines for calling MLIPs though persistent python instances
!> enabled through the fmlip_relay submodule

!=========================================================================================!
module mlip_sc
  use crest_parameters
  use strucrd
  use iomod
#ifdef WITH_FMLIP_RELAY
  use fmlip_relay_client
#endif
  implicit none
  !>--- private module variables and parameters
  private

  character(len=*),parameter :: basebin = 'fmlip-relay-server'

  public :: mlip_params
  type :: mlip_params
    integer :: BASE_PORT = 54320
    integer :: TIMEOUT_SEC = 120
    character(len=:),allocatable :: backend
    character(len=:),allocatable :: modelpath
    character(len=:),allocatable :: modelsize
    integer :: iid = 0
  end type mlip_params

  public :: mlip_engrad_core,fmlip_relay_init,mlips_shutdown

  integer,parameter  :: nopbc(3) = (/0,0,0/)
  integer,parameter  :: allpbc(3) = (/1,1,1/)
  real(wp),parameter :: bigcell(3,3) = reshape( &
    & (/10000.0_wp,0.0_wp,0.0_wp, &
    &   0.0_wp,10000.0_wp,0.0_wp, &
    &   0.0_wp,0.0_wp,10000.0_wp/), [3,3])

  external creststop
!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine fmlip_relay_init(MPAR,iid)
    type(mlip_params),intent(inout) :: MPAR
    integer,intent(in) :: iid
    integer :: io,tmpport
    character(len=256) :: cmd,cmd_0,cmd_1
#ifdef WITH_FMLIP_RELAY
    if (.not.allocated(MPAR%backend)) then
      write (stdout,*)
      write (stdout,*) '** ERROR ** No model backend selected for MLIP'
      write (stdout,*)
      call creststop(20)
    end if

    call checkprog_silent(basebin,verbose=.false.,iostat=io)
    if (io .ne. 0) then
      write (stdout,*)
      write (stdout,*) '** ERROR ** can not find socket server for MLIPs '//basebin
      write (stdout,*) ' Make sure you install it from the fmlip_relay subproject via pip'
      write (stdout,*)
      call creststop(20)
    end if

    !> check for already running instances that may need reinitialization
    !> or rather, shutdown first
    if (MPAR%iid .ne. 0) then
      call mlip_finalize(MPAR%iid,io)
    end if

    !> options prepping
    tmpport = MPAR%BASE_PORT+iid
    write(cmd_1,'("--dtype float64")')

    select case (MPAR%backend)
    case ('mace_off','mace_mp')

      if (allocated(MPAR%modelpath)) then
        if (.not.file_exists(MPAR%modelpath)) then
          write (stdout,*)
          write (stdout,*) '** ERROR ** model path allocated but can not find '//trim(MPAR%modelpath)
          write (stdout,*) 
          call creststop(20)
        end if
        write (cmd,'(a,1x,a,1x,i0,2(1x,a,1x,a),1x,a)') basebin,'--port',tmpport,'--backend', &
          & 'mace','--model',trim(MPAR%modelpath),trim(cmd_1)
      else
        cmd_0 = ''
        if (allocated(MPAR%modelsize)) write (cmd_0,'(a,1x,a)') '--mace_model',trim(MPAR%modelsize)
        write (cmd,'(a,1x,a,1x,i0,2(1x,a,1x,a),1x,a)') basebin,'--port',tmpport,'--backend', &
        & trim(MPAR%backend),trim(cmd_0),'',trim(cmd_1)
      end if

    case default

      if (allocated(MPAR%modelpath)) then
        if (.not.file_exists(MPAR%modelpath)) then
          write (stdout,*)
          write (stdout,*) '** ERROR ** model path allocated but can not find '//trim(MPAR%modelpath)
          write (stdout,*)
          call creststop(20)
        end if
      end if

      write (cmd,'(a,1x,a,1x,i0,2(1x,a,1x,a),1x,a)') basebin,'--port',tmpport,'--backend', &
        & trim(MPAR%backend),'--model',trim(MPAR%modelpath),trim(cmd_1)
    end select

    call mlip_init(iid,tmpport,trim(cmd)//' 2>/dev/null',MPAR%TIMEOUT_SEC,io)
    if (io /= MLIP_OK) then
      write (stdout,*)
      write (stdout,*) '** ERROR ** failed to initialize MLIP server'
      write (stdout,*)
      call creststop(1)
    end if
    !> Test it
    call mlip_ping(iid,io)
    if (io /= MLIP_OK) then
      write (stdout,*)
      write (stdout,*) '** ERROR ** failed to ping MLIP server'
      call creststop(1)
    end if

    MPAR%iid = iid

#else /* WITH_FMLIP_RELAY */
    write (stdout,*) 'Error: Compiled without fmlip-relay support!'
    write (stdout,*) 'Use -DWITH_FMLIP_RELAY=true in the setup to enable this function'
    write (stdout,*)
    call creststop(20)
#endif
  end subroutine fmlip_relay_init

  subroutine mlip_engrad_core(mol,MPAR,energy,gradient,iostatus)
    type(coord),intent(in) :: mol
    type(mlip_params)    :: MPAR
    real(wp),intent(out)   :: energy
    real(wp),intent(out)   :: gradient(3,mol%nat)
    integer,intent(out)    :: iostatus

    real(wp) :: stress(3,3)

    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    iostatus = 1

#ifdef WITH_FMLIP_RELAY
    if (allocated(mol%lat)) then
      call mlip_compute(MPAR%iid,mol%nat,mol%at,mol%xyz*autoaa,mol%lat,allpbc,0, &
      &                 energy,gradient,stress,iostatus)
    else
      call mlip_compute(MPAR%iid,mol%nat,mol%at,mol%xyz*autoaa,bigcell,nopbc,0, &
      &                 energy,gradient,stress,iostatus)
    end if

    !> CREST always works with atomic units, convert from eV and Angstroem:
    energy = energy/autoev
    gradient(:,:) = -gradient(:,:)*(1.0_wp/(autoev*aatoau))
#endif
  end subroutine mlip_engrad_core

!========================================================================================!

  subroutine mlips_shutdown()
    integer :: io
#ifdef WITH_FMLIP_RELAY
    call mlip_finalize_all(io)
#endif
  end subroutine mlips_shutdown

!========================================================================================!
end module mlip_sc
