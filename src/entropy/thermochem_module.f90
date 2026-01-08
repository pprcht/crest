module thermochem_module
  use crest_parameters
  use getsymmetry
  use hessian_tools
  use atmasses,only:molweight
  use iomod,only:to_lower
  use axis_module
  implicit none
  private

  public calcthermo,calc_thermo_from_hess

contains

  subroutine prepthermo(nat,at,xyz,pr,molmass,rabc,avmom,symnum,symchar,iunit)
!***********************************************************************
!* Prepare the calculation of thermodynamic properties of a structure
!* In particular, determine rotational constants and check the symmetry
!***********************************************************************
    implicit none
    integer,intent(in)     :: nat
    integer,intent(in)     :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat) !> in Angstroem
    logical,intent(in)     :: pr
    real(wp),intent(out)   :: molmass
    real(wp),intent(inout) :: rabc(3)
    real(wp),intent(out)   :: avmom
    real(wp),intent(out)   :: symnum
    integer,intent(in)     :: iunit

    real(wp) :: a,b,c
    character(len=4) :: sfsym
    character(len=3) :: sym,symchar
    real(wp),parameter :: desy = 0.1_wp
    integer,parameter  :: maxat = 200

    !>--- molecular mass in amu
    molmass = molweight(nat,at)

    if (pr) then
      write (iunit,'(1x,a,f15.2)') 'Mol. weight /amu  : ',molmass
    end if

    !>--- rotational constants in cm-1
    rabc = 0.0d0
    call axis(nat,at,xyz,rabc(1:3),avmom)
    a = rabc(3)
    b = rabc(2)
    c = rabc(1)
    rabc(1) = a
    rabc(3) = c
    if (pr) then
      write (iunit,'(1x,a,3f15.2)') 'Rot. const. /MHz  : ',rabc(1:3)
    end if
    rabc = rabc/2.99792458d+4   ! MHz to cm-1
    if (pr) then
      write (iunit,'(1x,a,3f15.2)') 'Rot. const. /cm-1 : ',rabc(1:3)
    end if

    !>--- symmetry number from rotational symmetry
    xyz = xyz/bohr
    call getsymmetry2(.false.,6,nat,at,xyz,desy,maxat,sfsym)
    xyz = xyz*bohr
    sym = sfsym(1:3)
    symchar = sym
    symnum = 1.0d0
    if (a .lt. 1.d-9.or.b .lt. 1.d-9.or.c .lt. 1.d-9) then
      if (index(sym,'d') .ne. 0) symnum = 2.0d0
    else
      call to_lower(sym)
      if (index(sym,'c2') .ne. 0) symnum = 2.0d0
      if (index(sym,'s4') .ne. 0) symnum = 2.0d0
      if (index(sym,'c3') .ne. 0) symnum = 3.0d0
      if (index(sym,'s6') .ne. 0) symnum = 3.0d0
      if (index(sym,'c4') .ne. 0) symnum = 4.0d0
      if (index(sym,'s8') .ne. 0) symnum = 4.0d0
      if (index(sym,'c5') .ne. 0) symnum = 5.0d0
      if (index(sym,'c6') .ne. 0) symnum = 6.0d0
      if (index(sym,'c7') .ne. 0) symnum = 7.0d0
      if (index(sym,'c8') .ne. 0) symnum = 8.0d0
      if (index(sym,'c9') .ne. 0) symnum = 9.0d0
      if (index(sym,'d2') .ne. 0) symnum = 4.0d0
      if (index(sym,'d3') .ne. 0) symnum = 6.0d0
      if (index(sym,'d4') .ne. 0) symnum = 8.0d0
      if (index(sym,'d5') .ne. 0) symnum = 10.0d0
      if (index(sym,'d6') .ne. 0) symnum = 12.0d0
      if (index(sym,'d7') .ne. 0) symnum = 14.0d0
      if (index(sym,'d8') .ne. 0) symnum = 16.0d0
      if (index(sym,'d9') .ne. 0) symnum = 18.0d0
      if (index(sym,'t') .ne. 0) symnum = 12.0d0
      if (index(sym,'td') .ne. 0) symnum = 12.0d0
      if (index(sym,'th') .ne. 0) symnum = 12.0d0
      if (index(sym,'o') .ne. 0) symnum = 24.0d0
      if (index(sym,'oh') .ne. 0) symnum = 24.0d0
      if (index(sym,'ih') .ne. 0) symnum = 60.0d0
    end if

    if (pr) then
      write (iunit,'(1x,a,4x,a)') 'Symmetry:',sym
    end if
    return
  end subroutine prepthermo

  subroutine calcthermo(nat,at,xyz,freq,pr,ithr,fscal,sthr,nt,temps, &
      &      et,ht,gt,stot,iunit_in)
!**************************************************************
!* Calculate thermodynamic contributions for a given structure
!* from it's frequencies (from second derivatives/the Hessian)
!* Based on xtb's "print_thermo" routine
!**************************************************************
    use crest_parameters,only:wp,bohr,stdout
    use crest_thermo
    use atmasses,only:molweight
    use iomod,only:to_lower
    implicit none
    integer,intent(in)     :: nat
    integer,intent(in)     :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)  !in Bohr
    real(wp),intent(inout) :: freq(3*nat) !in cm-1
    logical,intent(in)     :: pr
    real(wp),intent(in) :: ithr     !imag. inv. in cm-1
    real(wp),intent(in) :: fscal    !freq scaling
    real(wp),intent(in) :: sthr     !rotor cut
    integer,intent(in)  :: nt
    real(wp),intent(in) :: temps(nt)
    integer,intent(in),optional     :: iunit_in
    real(wp) :: et(nt)          !< enthalpy in Eh
    real(wp) :: ht(nt)          !< enthalpy in Eh
    real(wp) :: gt(nt)          !< free energy in Eh
    real(wp) :: stot(nt)        !< entropy in cal/molK
    real(wp) :: ts(nt)          !< entropy*T in Eh
    real(wp) :: rabc(3),a,b,c
    real(wp) :: avmom
    real(wp) :: molmass
    real(wp) :: sym
    real(wp) :: zp
    character(len=3) :: symchar
    logical :: pr2
    logical :: linear = .false.
    logical :: atom = .false.
    integer :: nvib_theo
    integer :: nvib,nimag
    real(wp) :: vibthr
    real(wp),allocatable :: vibs(:)

    integer :: i,j,iunit
    integer :: n3,rt
    real(wp) :: adum(nt)
    character(len=64) :: atmp

    character(len=*),parameter :: outfmt = &
    &  '(9x,"::",1x,a,f24.12,1x,a,1x,"::")'
    character(len=*),parameter :: dblfmt = &
    &  '(10x,":",2x,a,f24.7,1x,a,1x,":")'
    character(len=*),parameter :: intfmt = &
    &  '(10x,":",2x,a,i24,       6x,":")'
    character(len=*),parameter :: chrfmt = &
    &  '(10x,":",2x,a,a24,       6x,":")'

    real(wp),parameter :: autorcm = 219474.63067_wp
    real(wp),parameter :: rcmtoau = 1.0_wp/autorcm
    real(wp),parameter :: autocal = 627.50947428_wp*1000.0_wp

    xyz = xyz*autoaa

    if (present(iunit_in)) then
      iunit = iunit_in
    else
      iunit = stdout
    end if

    call prepthermo(nat,at,xyz,pr,molmass,rabc,avmom,sym,symchar,iunit)

    n3 = 3*nat
    allocate (vibs(n3))
    vibthr = 1.0
    a = rabc(1)
    b = rabc(2)
    c = rabc(3)

    nvib_theo = 3*nat-6
    if (c .lt. 1.d-10.or.(symchar == 'din')) linear = .true.
    if (linear) nvib_theo = 3*nat-5

    if (a+b+c .lt. 1.d-6) then
      atom = .true.
      nvib = 0
      nvib_theo = 0
    end if

    nvib = 0
    vibs = 0.0
    do i = 1,n3
      if (abs(freq(i)) .gt. vibthr) then
        nvib = nvib+1
        vibs(nvib) = freq(i)
      end if
    end do
    !> scale
    vibs(1:nvib) = vibs(1:nvib)*fscal

    !> invert imaginary modes
    nimag = 0
    do i = 1,nvib
      if (vibs(i) .lt. 0.and.vibs(i) .gt. ithr) then
        vibs(i) = -vibs(i)
        if (pr) write (stdout,'(a,i5," :",f10.2)') 'Inverting frequency',i,vibs(i)
      end if
      if (vibs(i) < 0.0) then
        nimag = nimag+1
      end if
    end do

    if (pr) then
      write (iunit,'(a)')
      write (iunit,'(10x,51("."))')
      write (iunit,'(10x,":",22x,a,22x,":")') "SETUP"
      write (iunit,'(10x,":",49("."),":")')
      write (iunit,intfmt) "# frequencies    ",nvib
      write (iunit,intfmt) "# imaginary freq.",nimag
      write (atmp,*) linear
      write (iunit,chrfmt) "linear?          ",trim(atmp)
      write (iunit,chrfmt) "symmetry         ",adjustr(symchar)
      write (iunit,intfmt) "rotational number",nint(sym)
      write (iunit,dblfmt) "scaling factor   ",fscal,"    "
      write (iunit,dblfmt) "rotor cutoff     ",sthr,"cm⁻¹"
      write (iunit,dblfmt) "imag. cutoff     ",ithr,"cm⁻¹"
      write (iunit,'(10x,":",49("."),":")')
    end if

    vibs = vibs*rcmtoau   ! thermodyn needs vibs and zp in Eh

    zp = 0.5_wp*sum(vibs(1:nvib))
    adum = abs(temps-298.15d0)
    rt = minloc(adum,1)  !temperature closest to 298.15 is the ref.
    do j = 1,nt
      if ((j == rt).and.pr) then
        pr2 = .true.
      else
        pr2 = .false.
      end if
      if (pr2) then
        call print_thermo_sthr_ts(iunit,nvib,vibs,avmom,sthr,temps(j))
      end if
      call thermodyn(iunit,a,b,c,avmom,linear,atom,sym,molmass,vibs,nvib, &
      & temps(j),sthr,et(j),ht(j),gt(j),ts(j),zp,pr2)
      stot(j) = (ts(j)/temps(j))*autocal
    end do

    if (pr) then
      write (iunit,'(a)')
      write (iunit,'(a10)',advance='no') "T/K"
      write (iunit,'(a16)',advance='no') "H(0)-H(T)+PV"
      write (iunit,'(a16)',advance='no') "H(T)/Eh"
      write (iunit,'(a16)',advance='no') "T*S/Eh"
      write (iunit,'(a16)',advance='no') "G(T)/Eh"
      write (iunit,'(a)')
      write (iunit,'(3x,72("-"))')
      do i = 1,nt
        write (iunit,'(3f10.2)',advance='no') temps(i)
        write (iunit,'(3e16.6)',advance='no') ht(i)
        write (iunit,'(3e16.6)',advance='no') et(i)
        write (iunit,'(3e16.6)',advance='no') ts(i)
        write (iunit,'(3e16.6)',advance='no') gt(i)
        if (i == rt.and.nt > 1) then
          write (iunit,'(1x,"(used)")')
        else
          write (iunit,'(a)')
        end if
      end do
      write (iunit,'(3x,72("-"))')
    end if

    xyz = xyz*aatoau

    deallocate (vibs)
    return
  end subroutine calcthermo

  subroutine calc_thermo_from_hess(mol,hess,pr,nt,temps,ithr,&
  & fscal,sthr,et,ht,gt,stot, etot)
    type(coord),intent(inout) :: mol
    integer :: nat3
    integer :: io,iunit
    logical :: pr
    real(wp) :: ithr,fscal,sthr
    real(wp),intent(in) :: temps(nt)
    integer,intent(in) :: nt
    real(wp),allocatable,intent(out) :: et(:),ht(:),gt(:),stot(:)
    real(wp),intent(inout) :: hess(:,:)
    real(wp),allocatable :: freq(:)
    real(wp), intent(in) :: etot
    real(wp) :: zpve
    integer :: nrt
    real(wp),allocatable :: int_temps(:)
    character(len=*),parameter :: outfmt = &
    &  '(10x,"::",1x,a,f24.12,1x,a,1x,"::")'

    nat3 = 3*mol%nat
    allocate (freq(nat3))
    allocate (et(nt))
    allocate (ht(nt))
    allocate (gt(nt))
    allocate (stot(nt))
    allocate (int_temps(nt))

    int_temps = abs(temps-298.15_wp)
    nrt = minloc(int_temps(:),1)

    call prj_mw_hess(mol%nat,mol%at,nat3,mol%xyz,hess)

    call frequencies(mol%nat,mol%at,mol%xyz,nat3,hess,freq,io)

    call calcthermo(mol%nat,mol%at,mol%xyz,freq,pr,ithr,fscal,sthr,nt,temps, &
        &      et,ht,gt,stot)

    zpve = et(nrt)-ht(nrt)
    write (stdout,*)
    write (stdout,'(10x,a)') repeat(':',50)
    write (stdout,'(10x,"::",7x,a,f12.2,1x,a,8x,"::")') "THERMODYNAMICS at",temps(nrt),'K'
    write (stdout,'(10x,a)') repeat(':',50)
    write (stdout,outfmt) 'TOTAL FREE ENERGY',etot+gt(nrt),'Eh'
    write (stdout,'(10x,a)') '::'//repeat('-',46)//'::'
    write (stdout,outfmt) 'total energy     ',etot,'Eh'
    write (stdout,outfmt) 'ZPVE             ',zpve,'Eh'
    write (stdout,outfmt) 'G(RRHO) w/o ZPVE ',gt(nrt)-zpve,'Eh'
    write (stdout,outfmt) 'G(RRHO) total    ',gt(nrt),'Eh'
    write (stdout,'(10x,a)') repeat(':',50)

  end subroutine calc_thermo_from_hess


end module thermochem_module
