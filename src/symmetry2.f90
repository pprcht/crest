! This file is part of xtb, modified for crest
!
! Copyright (C) 2017-2020 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

Module getsymmetry

  private
  public :: getsymmetry2

! ══════════════════════════════════════════════════════════════════════════════
contains
! ══════════════════════════════════════════════════════════════════════════════

  subroutine getsymmetry2(pr,iunit,n,iat,xyz,symthr,maxatdesy,sfsym)
    use symmetry_i,only:schoenflies
    use iso_fortran_env,only:wp => real64
    implicit none
    integer,intent(in) :: iunit
    integer :: n,iat(n),maxatdesy
    real(wp) :: xyz(3,n)
    real(wp) :: symthr
    Character(len=*) :: sfsym
    logical :: pr
    character(len=8) :: atmp
    Real(wp) :: paramar(11)  !parameter array for schoenflies

    if (n .gt. maxatdesy) then
      if (pr) write (iunit,*) 'symmetry recognition skipped because # atoms >',maxatdesy
      sfsym = 'none'
      return
    end if

    if (pr) write (iunit,'(a)')
    !parameters for symmetry recognition:
    paramar(1) = -1          ! verbose, increase for more detailed output (to stdout)
    paramar(2) = 10          ! MaxAxisOrder
    paramar(3) = 100         ! MaxOptCycles
    paramar(4) = 0.001d0     ! ToleranceSame
    paramar(5) = 0.5d0       ! TolerancePrimary
    paramar(6) = symthr      ! ToleranceFinal, THIS IS THE IMPORTANT VALUE
    paramar(7) = 0.5d0       ! MaxOptStep
    paramar(8) = 1.0D-7      ! MinOptStep
    paramar(9) = 1.0D-7      ! GradientStep
    paramar(10) = 1.0D-8     ! OptChangeThreshold
    paramar(11) = 5          ! OptChangeHits

    atmp = '        '
    call schoenflies(n,iat,xyz,atmp,paramar)

    !TM stuff (trafo table)
    sfsym(1:3) = atmp(1:3)
    if (sfsym(1:1) .eq. 'D') sfsym(1:1) = 'd'
    if (sfsym(1:1) .eq. 'C') sfsym(1:1) = 'c'
    if (sfsym(1:1) .eq. 'T') sfsym(1:1) = 't'
    if (sfsym(1:1) .eq. 'O') sfsym(1:1) = 'o'
    if (sfsym(1:1) .eq. 'I') sfsym(1:1) = 'i'
    if (sfsym(1:1) .eq. 'S') sfsym(1:1) = 's'
    if (sfsym .eq. 'dih') sfsym = 'd6h'
    if (sfsym .eq. 'civ') sfsym = 'c6v'
    if (sfsym(3:3) .gt. 'v'.or.sfsym(3:3) .lt. 'a') sfsym(3:3) = ' '

    if (pr) then
      write (iunit,'(a3,'' symmetry found (for desy threshold: '',e9.2,'')'')') sfsym,symthr
    end if
  End subroutine getsymmetry2

! ══════════════════════════════════════════════════════════════════════════════
end module getsymmetry
