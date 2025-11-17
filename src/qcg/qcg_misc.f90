!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 Christoph Plett, Sebastian Spicher, Philipp Pracht
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

subroutine xtb_sp_qcg(env,fname,success,eout)
!********************************************************
!*  xtb_sp_qcg
!*  A quick single point xtb calculation without wbo
!********************************************************
  use crest_parameters
  use iomod
  use crest_data
  use crest_calculator
  use strucrd
  implicit none
  type(systemdata) :: env
  character(len=*),intent(in) :: fname
  logical,intent(out) :: success
  real(wp),intent(out) :: eout

  character(len=512) :: jobcall
  character(len=*),parameter :: pipe = ' > xtb.out 2> /dev/null'
  logical,parameter :: debug = .false.
  integer :: io,T,Tn

  success = .false.
  eout = 0.0_wp

  if (env%legacy) then
!>----------------------------------------------
!> The original implementation with systemcall
    call remove('gfnff_topo')
    call remove('energy')
    call remove('charges')
    call remove('xtbrestart')

!---- setting threads
    call new_ompautoset(env,'auto',1,T,Tn)

!---- jobcall
    write (jobcall,'(a,1x,a,1x,a,'' --sp '',a,1x,a)') &
    &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)

    if (debug) write (stdout,*) trim(jobcall)
    call command(trim(jobcall),io)
    call grepval('xtb.out','| TOTAL ENERGY',success,eout)
!---- cleanup
    call remove('energy')
    call remove('charges')
    call remove('xtbrestart')
    call remove('xtbtopo.mol')
    call remove('gfnff_topo')
  else
!>---------------------------------------------
!> New implementation with calculator and api
    block
      type(calcdata) :: calc
      type(coord) :: mol
      real(wp),allocatable :: gradtmp(:,:)

      call mol%open(fname)
      allocate (gradtmp(3,mol%nat))
      call env2calc(env,calc,mol)
      if (debug) call calc%info(stdout)
      call engrad(mol,calc,eout,gradtmp,io)
      success = (io == 0)
    end block
  end if
end subroutine xtb_sp_qcg

!--------------------------------------------------------------------------------------------
! A quick single xtb optimization gets mol and overwrites it with optimized stuff
!--------------------------------------------------------------------------------------------
subroutine xtb_opt_qcg(env,mol,constrain)
  use crest_parameters
  use iomod
  use crest_data
  use qcg_coord_type
  use strucrd

  implicit none
  type(systemdata),intent(in) :: env
  type(coord_qcg),intent(inout) :: mol

  character(:),allocatable :: fname
  character(len=512) :: jobcall
  logical :: constrain,const
  real(wp) :: energy
  integer :: io,T,Tn
  character(stdout),parameter :: pipe = ' > xtb_opt.out 2> /dev/null'
  logical,parameter :: debug = .false.

  if (env%legacy) then
    !> LEGACY version with syscall

    !--- Write coordinated
    fname = 'coord'
    call mol%write(fname)

    !---- setting threads
    call new_ompautoset(env,'auto',1,T,Tn)

    !---- jobcall & Handling constraints
    if (constrain.AND.env%cts%used) then
      call write_constraint(env,fname,'xcontrol')
      call mol%write('coord.ref')
      write (jobcall,'(a,1x,a,1x,a,'' --opt --input xcontrol '',a,1x,a)') &
      &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)
    else
      write (jobcall,'(a,1x,a,1x,a,'' --opt '',a,1x,a)') &
      &     trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)
    end if

    call command(trim(jobcall),io)
    !---- cleanup
    call rdcoord('xtbopt.coord',mol%nat,mol%at,mol%xyz)
    call remove('energy')
    call remove('charges')
    call remove('xtbrestart')
    call remove('xtbtopo.mol')
    call remove('gfnff_topo')

  else

    !> NEW version with calculator
    call new_ompautoset(env,'max',1,T,Tn)
    block
      use crest_calculator
      use optimize_module
      type(calcdata) :: calc
      type(coord) :: molin,molout
      real(wp),allocatable :: gradtmp(:,:)

      allocate (gradtmp(3,mol%nat))
      molin = mol%as_coord()
      call env2calc(env,calc,molin)
      if (debug) call calc%info(stdout)

      call optimize_geometry(molin,molout,calc,energy,gradtmp,debug,.false.,io)

      deallocate (gradtmp)
      if (io == 0) then
        call mol%from_coord(molout)
      else
        write (stdout,*) 'FAILURE in QCG optimization!'
        write (stdout,*) 'Stopping run to avoid unecessary compuations'
        call creststop(status_safety)
      end if
    end block

  end if
end subroutine xtb_opt_qcg

!___________________________________________________________________________________
!
! An xTB docking on all available threads
!___________________________________________________________________________________

subroutine xtb_dock(env,fnameA,fnameB,solu,clus)
  use crest_parameters
  use iomod
  use crest_data
  use qcg_coord_type
  implicit none

  type(systemdata)               :: env
  type(coord_qcg),intent(in)     :: solu,clus
  character(len=*),intent(in)    :: fnameA,fnameB
  character(len=80)              :: pipe
  character(len=512)             :: jobcall
  integer                        :: i,ich,T,Tn

  call remove('xtb_dock.out')
  call remove('xcontrol')

  pipe = ' 2>/dev/null'

!---- writing wall pot in xcontrol
  call write_wall(env,solu%nat,solu%ell_abc,clus%ell_abc,'xcontrol')

!---- Write directed stuff, if requested
  if (allocated(env%directed_file)) then
    do i = 1,size(env%directed_number)
      if &
      & ((i == 1.and.env%directed_number(i) >= clus%nmol).OR. &
      & (env%directed_number(i) >= clus%nmol.and.env%directed_number(i-1) < clus%nmol)) &
      & then
        open (newunit=ich,file='xcontrol',status='old',position='append',action='write')
        write (ich,'("$directed")')
        write (ich,'(a,1x,a)') 'atoms:',trim(env%directed_list(i,1))
        write (ich,'("$end")')
      end if
    end do
  end if

!--- Setting threads
  call new_ompautoset(env,'auto',1,T,Tn)

!--- Jobcall docking
  write (jobcall,'(a,1x,''dock'',1x,a,1x,a,1x,a,1x,f4.2,1x,''--nfrag1'',1x,i0,1x,a,1x, &
          & ''--input xcontrol > xtb_dock.out'',a)') &
          &     trim(env%ProgName),trim(fnameA),trim(fnameB),trim(env%gfnver),&
          &     env%optlev,solu%nat,trim(env%docking_qcg_flag),trim(pipe)
  call command(trim(jobcall))

! cleanup
  call remove('wbo')
  call remove('charges')
  call remove('xtbrestart')

end subroutine xtb_dock

!___________________________________________________________________________________
!
! An xTB optimization on all available threads
!___________________________________________________________________________________

subroutine opt_cluster(env,solu,clus,fname,without_pot)
  use crest_parameters
  use iomod
  use crest_data
  use qcg_coord_type

  implicit none

  type(systemdata)              :: env
  type(coord_qcg),intent(in)    :: solu,clus
  character(len=*),intent(in)   :: fname
  logical,optional,intent(in)   :: without_pot
  character(len=*),parameter    :: pipe = ' 2>/dev/null'
  character(len=:),allocatable  :: jobcall
  integer :: T,Tn

  call remove('xtb.out')

!---- writing wall pot in xcontrol
  if (.not.without_pot) then
    call write_wall(env,solu%nat,solu%ell_abc,clus%ell_abc,'xcontrol')
  end if

!--- Setting threads
  call new_ompautoset(env,'subprocess',1,T,Tn)

!--- Jobcall optimization
  jobcall = trim(env%ProgName)//' '//trim(fname)//' --opt '//optlevflag(env%optlev)
  jobcall = trim(jobcall)//' '//trim(env%gfnver)
  if (without_pot) then
    jobcall = trim(jobcall)//' '//trim(env%solv)
  end if
  jobcall = trim(jobcall)//' > xtb_opt.out 2>/dev/null'
  call command(trim(jobcall))

! cleanup
  call remove('wbo')
  call remove('charges')
  call remove('xtbrestart')

!--- Jobcall SP for gbsa model
  if (.not.without_pot) then
    jobcall = trim(env%ProgName)//' xtbopt.coord --sp '//trim(env%gfnver)
    jobcall = trim(jobcall)//' '//trim(env%solv)
    jobcall = trim(jobcall)//' > xtb_sp.out 2>/dev/null'
    call command(trim(jobcall)) 
  end if

! cleanup
  call remove('wbo')
  call remove('charges')
  call remove('xtbrestart')

end subroutine opt_cluster

!___________________________________________________________________________________
!
! xTB docking calculation performed in parallel
!___________________________________________________________________________________

subroutine ensemble_dock(env,outer_ell_abc,nfrag1,frag1_file,frag2_file,n_shell&
        &,n_solvent,NTMP,TMPdir,conv)
  use crest_parameters
  use iomod
  use crest_data
  use qcg_coord_type

  implicit none
  type(systemdata)                :: env

  character(len=*),intent(in)     :: TMPdir     !directory name
  integer,intent(in)              :: NTMP       !number of structures to be optimized
  integer,intent(in)              :: nfrag1     !#atoms of larger fragment
  integer,intent(in)              :: conv(env%nqcgclust+1)
  real(wp),intent(in)             :: outer_ell_abc(env%nqcgclust,3)
  integer,intent(in)              :: n_shell,n_solvent

  integer                         :: i,k
  integer                         :: vz,T,Tn
  character(len=20)               :: pipe
  character(len=1024)             :: jobcall
  character(len=512)              :: thispath,tmppath
  character(len=*),intent(in)    :: frag1_file
  character(len=*),intent(in)    :: frag2_file
  character(len=64)               :: frag1
  character(len=64)               :: frag2
  real(wp)                        :: percent
  character(len=2)                :: flag
  integer                         :: funit

! some options
  pipe = '2>/dev/null'
  frag1 = 'solvent_cluster.coord'
  frag2 = 'solvent'
  call getcwd(thispath)

! setting the threads for correct parallelization
  call new_ompautoset(env,'auto',NTMP,T,Tn)

  write (jobcall,'(a,1x,''dock'',1x,a,1x,a,1x,a,1x,f4.2,1x,''--nfrag1'',1x,i0,1x,&
          & ''--input xcontrol --fast > xtb_dock.out '',a)') &
          & trim(env%ProgName),trim(frag1_file),trim(frag2_file),&
          & trim(env%gfnver),env%optlev,nfrag1,trim(pipe)

  flag = '$'
  do i = 1,NTMP
    vz = i
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    call chdirdbug(trim(tmppath))
    open (newunit=funit,file='xcontrol')
    write (funit,'(a,"fix")') trim(flag)
    write (funit,'(3x,"atoms: 1-",i0)') n_shell !Initial number of atoms (starting solvent shell)
    write (funit,'(a,"wall")') trim(flag)
    write (funit,'(3x,"potential=polynomial")')
    write (funit,'(3x,"ellipsoid:",1x,3(g0,",",1x),i0,"-",i0)') outer_ell_abc(conv(vz),:), &
            & n_shell+1,n_shell+n_solvent !Initial number of atoms (starting solvent shell)
    close (funit)
    call chdirdbug(trim(thispath))
  end do

  k = 0 !counting the finished jobs

!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,NTMP,percent,k,TMPdir,conv,n_shell,n_solvent,jobcall )
!$omp single
  do i = 1,NTMP
    vz = i
    !$omp task firstprivate( vz ) private( tmppath )
    write (tmppath,'(a,i0)') trim(TMPdir),conv(vz)
    call command('cd '//trim(tmppath)//' && '//trim(jobcall))
    !$omp critical
    k = k+1
    percent = float(k)/float(NTMP)*100
    !$omp end critical
    !$omp end task
  end do

!$omp taskwait
!$omp end single
!$omp end parallel

!___________________________________________________________________________________
  call chdirdbug(trim(thispath))

end subroutine ensemble_dock

!___________________________________________________________________________________
!
! xTB CFF optimization performed in parallel
!___________________________________________________________________________________

subroutine cff_opt(postopt,env,fname,n12,NTMP,TMPdir,conv,nothing_added)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd
  implicit none

  type(systemdata)                :: env
  character(len=*),intent(in)     :: fname      !file base name
  character(len=*),intent(in)     :: TMPdir     !directory name
  integer,intent(inout)           :: NTMP       !number of structures to be optimized
  integer,intent(inout)           :: conv(env%nqcgclust+1)
  logical,intent(in)              :: postopt
  logical,intent(in)              :: nothing_added(env%nqcgclust)
  integer                         :: i,k,n12
  integer                         :: vz,T,Tn
  integer                         :: funit
  character(len=20)               :: pipe
  character(len=512)              :: thispath,tmppath
  character(len=1024)             :: jobcall
  character(len=2)                :: flag
  real(wp)                        :: percent

! setting the threads for correct parallelization
  call new_ompautoset(env,'auto',NTMP,T,Tn)

  if (postopt) then
    write (stdout,'(2x,''Starting optimizations + SP  of structures'')')
    write (stdout,'(2x,i0,'' jobs to do.'')') NTMP
  end if

! postopt eq true => post opt run, which has to be performed in every directory !!!
  if (postopt) then
    k = 0
    NTMP = env%nqcgclust
    do i = 1,env%nqcgclust
      k = k+1
      conv(k) = i
      conv(env%nqcgclust+1) = k
    end do
  end if
  pipe = '2>/dev/null'

  call getcwd(thispath)
  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    call chdirdbug(trim(tmppath))
    open (newunit=funit,file='xcontrol')
    if (n12 .ne. 0) then
      flag = '$'
      write (funit,'(a,"fix")') trim(flag)
      write (funit,'(3x,"atoms: 1-",i0)') n12 !Initial number of atoms (starting solvent shell)
    end if
    close (funit)
    if (postopt.and.nothing_added(i)) call remove('xcontrol')
    call chdirdbug(trim(thispath))
  end do

!--- Jobcall
  write (jobcall,'(a,1x,a,1x,a,'' --input xcontrol --opt '',i0,1x,a,'' >xtb.out'')') &
  &    trim(env%ProgName),trim(fname),trim(env%gfnver),nint(env%optlev),trim(pipe)

  if (NTMP .lt. 1) then
    write (stdout,'(2x,"No structures to be optimized")')
    return
  end if

  k = 0 !counting the finished jobs
  if (postopt) call printprogbar(0.0_wp)
!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,jobcall,NTMP,percent,k,TMPdir,conv )
!$omp single
  do i = 1,NTMP
    vz = i
    !$omp task firstprivate( vz ) private( tmppath )
    write (tmppath,'(a,i0)') trim(TMPdir),conv(vz)
    call command('cd '//trim(tmppath)//' && '//trim(jobcall))
    !$omp critical
    k = k+1
    percent = float(k)/float(NTMP)*100
    if (postopt) then
      call printprogbar(percent)
    end if
    !$omp end critical
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

!__________________________________________________________________________________

  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    call chdirdbug(trim(tmppath))
    call remove('xtbrestart')
    call chdirdbug(trim(thispath))
  end do

  !create the system call for sp (needed for gbsa model)
  write (jobcall,'(a,1x,a,1x,a,'' --sp '',a,1x,a,'' >xtb_sp.out'')') &
  &    trim(env%ProgName),'xtbopt.coord',trim(env%gfnver),trim(env%solv),trim(pipe)

  if (NTMP .lt. 1) then
    write (stdout,'(2x,"Nothing to do")')
    return
  end if

  k = 0 !counting the finished jobs
  if (postopt) call printprogbar(0.0_wp)
!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,jobcall,NTMP,percent,k,TMPdir,conv )
!$omp single
  do i = 1,NTMP
    vz = i
    !$omp task firstprivate( vz ) private( tmppath )
    write (tmppath,'(a,i0)') trim(TMPdir),conv(vz)
    call command('cd '//trim(tmppath)//' && '//trim(jobcall))
    !$omp critical
    k = k+1
    percent = float(k)/float(NTMP)*100
    if (postopt) then
      call printprogbar(percent)
    end if
    !$omp end critical
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

!___________________________________________________________________________________

  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),conv(i)
    call chdirdbug(trim(tmppath))
    call remove('xtbrestart')
    !call remove('xcontrol')
    call chdirdbug(trim(thispath))
  end do

  if (postopt) then
    write (stdout,*) ''
    write (stdout,'(2x,"done.")')
  end if

end subroutine cff_opt

!___________________________________________________________________________________
!
! xTB SP performed in parallel
!___________________________________________________________________________________

subroutine ens_sp(env,fname,NTMP,TMPdir)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd
  implicit none

  type(systemdata)                :: env
  character(len=*),intent(in)     :: fname      !file base name
  character(len=*),intent(in)     :: TMPdir     !directory name
  integer,intent(inout)           :: NTMP       !number of structures to be optimized

  integer                         :: i,k
  integer                         :: vz,T,Tn
  character(len=20)               :: pipe
  character(len=512)              :: thispath,tmppath
  character(len=1024)             :: jobcall
  real(wp)                        :: percent

! setting the threads for correct parallelization
  call new_ompautoset(env,'auto',NTMP,T,Tn)

  write (stdout,'(2x,''Single point computation with GBSA model'')')
  write (stdout,'(2x,i0,'' jobs to do.'')') NTMP

  pipe = '2>/dev/null'

  call getcwd(thispath)

  if (NTMP .lt. 1) then
    write (stdout,'(2x,"No structures to be optimized")')
    return
  end if

!--- Jobcall
  write (jobcall,'(a,1x,a,1x,a,'' --sp '',a,1x,a,'' > xtb_sp.out'')') &
  &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(env%solv),trim(pipe)

  k = 0 !counting the finished jobs
  call printprogbar(0.0_wp)
!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,NTMP,percent,k,TMPdir,jobcall )
!$omp single
  do i = 1,NTMP
    vz = i
    !$omp task firstprivate( vz ) private( tmppath )
    call initsignal()
    !$omp critical
    write (tmppath,'(a,i0)') trim(TMPdir),vz
    !$omp end critical
    call command('cd '//trim(tmppath)//' && '//trim(jobcall))
    !$omp critical
    k = k+1
    percent = float(k)/float(NTMP)*100
    call printprogbar(percent)
    !$omp end critical
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

!__________________________________________________________________________________

  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),i
    call chdirdbug(trim(tmppath))
    call remove('xtbrestart')
    call chdirdbug(trim(thispath))
  end do
  write (stdout,*) ''
  write (stdout,'(2x,"done.")')

end subroutine ens_sp

!___________________________________________________________________________________
!
! xTB Freq compuatation performed in parallel
!___________________________________________________________________________________

subroutine ens_freq(env,fname,NTMP,TMPdir,opt)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd
  implicit none

  type(systemdata)                :: env
  character(len=*),intent(in)     :: fname      !file base name
  character(len=*),intent(in)     :: TMPdir     !directory name
  integer,intent(inout)           :: NTMP       !number of structures to be optimized

  integer                         :: i,k
  integer                         :: vz,T,Tn
  character(len=20)               :: pipe
  character(len=512)              :: thispath,tmppath
  character(len=1024)             :: jobcall
  real(wp)                        :: percent
  logical                         :: opt

! setting the threads for correct parallelization
  call new_ompautoset(env,'auto',NTMP,T,Tn)

  write (stdout,'(2x,''Starting reoptimizations + Frequency computation of ensemble'')')
  write (stdout,'(2x,i0,'' jobs to do.'')') NTMP

  pipe = '2>/dev/null'

  call getcwd(thispath)

  if (NTMP .lt. 1) then
    write (stdout,'(2x,"No structures to be optimized")')
    return
  end if

  k = 0 !counting the finished jobs
  call printprogbar(0.0_wp)

!--- Jobcall
  if (.not.opt) then
    write (jobcall,'(a,1x,a,1x,a,'' --hess '',a,'' >xtb_freq.out'')') &
     &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(pipe)
  else
    write (jobcall,'(a,1x,a,1x,a,'' --ohess '',a,'' >xtb_freq.out'')') &
    &    trim(env%ProgName),trim(fname),trim(env%gfnver),trim(pipe)
  end if

!___________________________________________________________________________________

!$omp parallel &
!$omp shared( vz,NTMP,percent,k,TMPdir,jobcall )
!$omp single
  do i = 1,NTMP
    vz = i
    !$omp task firstprivate( vz ) private( tmppath )
    write (tmppath,'(a,i0)') trim(TMPdir),i
    call command('cd '//trim(tmppath)//' && '//trim(jobcall))
    !$omp critical
    k = k+1
    percent = float(k)/float(NTMP)*100
    call printprogbar(percent)
    !$omp end critical
    !$omp end task
  end do
!$omp taskwait
!$omp end single
!$omp end parallel

!__________________________________________________________________________________

  do i = 1,NTMP
    write (tmppath,'(a,i0)') trim(TMPdir),i
    call chdirdbug(trim(tmppath))
    call remove('xtbrestart')
    call chdirdbug(trim(thispath))
  end do
  write (stdout,*) ''
  write (stdout,'(2x,"done.")')

end subroutine ens_freq

!============================================================!
! subroutine wr_cluster_cut
! Cuts a cluster file and and writes the parts
!
! On Input: fname          - name of the coord file
!           n1             - number of atoms fragment1
!           n2             - number of atmos fragment2
!           iter           - number of solvent molecules
!           fname_solu_cut - name of outputfile fragment1
!           fname_solv_cut - name of outputfile fragment2
!
!============================================================!

subroutine wr_cluster_cut(fname_cluster,n1,n2,iter,fname_solu_cut,fname_solv_cut)
  use crest_parameters
  use strucrd

  implicit none
  integer,intent(in)         :: n1,n2,iter
  real(wp)                    :: xyz1(3,n1)
  real(wp)                    :: xyz2(3,n2*iter)
  integer                     :: at1(n1),at2(n2*iter)
  character(len=*),intent(in) :: fname_cluster,fname_solu_cut,fname_solv_cut
  character(len=256)         :: atmp
  character(len=2)           :: a2
  integer                     :: ich,i,k,stat,io,io2

  ich = 142
  open (unit=ich,file=fname_cluster,iostat=stat)
  read (ich,'(a)') atmp
  k = 1
  do i = 1,n1
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit
    atmp = adjustl(atmp)
    call coordline(atmp,a2,xyz1(1:3,k),io2)
    at1(k) = e2i(a2)
    k = k+1
  end do
  k = 1
  do i = 1,n2*iter
    read (ich,'(a)',iostat=io) atmp
    if (io < 0) exit
    atmp = adjustl(atmp)
    call coordline(atmp,a2,xyz2(1:3,k),io2)
    at2(k) = e2i(a2)
    k = k+1
  end do

  call wrc0(fname_solu_cut,n1,at1,xyz1)
  call wrc0(fname_solv_cut,n2*iter,at2,xyz2)
  close (ich)

end subroutine wr_cluster_cut

!----------------------------------------------------------------------------
! write a wall potential in a file used as xtb input

subroutine write_wall(env,n1,rabc1,rabc12,fname)
  use crest_parameters
  use crest_data

  implicit none

  type(systemdata)     :: env
  integer,intent(in)  :: n1
  real(wp),intent(in)  :: rabc1(3),rabc12(3)
  character(len=8)    :: flag
  character(len=*)    :: fname
  integer :: funit

  open (newunit=funit,file=fname)
  flag = '$'
  write (funit,'(a,"wall")') trim(flag)
  write (funit,'(3x,"potential=polynomial")')
  write (funit,'(3x,"ellipsoid:",1x,3(g0,",",1x),"all")') rabc12
  write (funit,'(3x,"ellipsoid:",1x,3(g0,",",1x),"1-",i0)') rabc1,n1
  if (env%constrain_solu) then
    write (funit,'("$fix")')
    write (funit,'(3x,"atoms: 1-",i0)') n1
  end if
  call write_cts(funit,env%cts)
  call write_cts_biasext(funit,env%cts)
  if (env%cts%used) then !Only, if user set constrians is an $end written
    write (funit,'(a)') '$end'
  end if
  close (funit)

end subroutine write_wall

subroutine check_dock(neg_E)
  use crest_parameters
  use crest_data
  use iomod,only:minigrep,grepval

  implicit none
  real(wp)             :: int_E
  logical,intent(out) :: neg_E
  logical :: ex
  character(len=*),parameter :: filename = 'xtbscreen.xyz'

  neg_E = .false.
  int_E = 0.0_wp

  call minigrep('xtb_dock.out','  Lowest Interaction Energy: ********** kcal/mol',ex)
  if (ex) return

  call grepval('xtb_dock.out','Lowest Interaction Energy:',ex,int_E)

  if (ex.and.int_E < 0.0_wp) neg_E = .true.

end subroutine check_dock

subroutine write_constraint(env,coord_name,fname)
  use crest_parameters
  use crest_data
  use iomod

  implicit none

  type(systemdata)     :: env
  character(len=*),intent(in)    :: fname,coord_name
  integer :: funit

  call copysub(coord_name,'coord.ref')
  open (newunit=funit,file=fname)
  call write_cts(funit,env%cts)
  call write_cts_biasext(funit,env%cts)
  if (env%cts%used) then !Only, if user set constrians is an $end written
    write (funit,'(a)') '$end'
  end if
  close (funit)

end subroutine write_constraint

!==============================================================================!

subroutine get_interaction_E(env,solu,solv,clus,iter,E_inter)
  use iso_fortran_env,wp => real64
  use crest_data
  use iomod
  use qcg_coord_type
  use strucrd
  implicit none

  type(systemdata)           :: env
  type(coord_qcg),intent(in) :: solu,solv,clus
  real(wp)                   :: e_cluster,e_solute,e_solvent
  real(wp)                   :: E_inter(env%nsolv)           ! interaction energy
  integer                    :: iter
  logical                    :: e_there

  call remove('cluster.coord')

!--- Prepare input coordinate files
  call clus%write('cluster.coord')
  call wr_cluster_cut('cluster.coord',solu%nat,solv%nat,iter,'solute_cut.coord','solvent_cut.coord')

!--- Perform single point calculations and recieve energies
  call xtb_sp_qcg(env,'solute_cut.coord',e_there,e_solute)
  if (.not.e_there) write (stdout,*) 'Solute energy not found'

  call xtb_sp_qcg(env,'solvent_cut.coord',e_there,e_solvent)
  if (.not.e_there) write (stdout,*) 'Solvent energy not found'

  call xtb_sp_qcg(env,'cluster.coord',e_there,e_cluster)
  if (.not.e_there) write (stdout,*) 'Cluster energy not found'

  E_inter(iter) = e_cluster-e_solute-e_solvent

end subroutine get_interaction_E

!===============================================================================!
subroutine chdirdbug(path)
  implicit none
  character(len=*),intent(in) :: path
  logical,parameter :: debug = .true.
  character(len=500) :: debugpath
  call chdir(path)
  if (debug) then
    call getcwd(debugpath)
    write (*,'(a,a)') '>>>>>>> NOW IN ',trim(debugpath)
  end if
end subroutine chdirdbug

