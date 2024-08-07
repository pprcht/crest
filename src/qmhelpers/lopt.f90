! Copyright (C) 2017-2020 Stefan Grimme, Sebastian Ehlert (xtb)
! Copyright (C) 2024 Philipp Pracht
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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

module lopt_mod
  use iso_fortran_env,only:wp => real64
  public

  logical,private,parameter :: debug = .false.

!================================================================================!
!================================================================================!
contains  !> MODULE PROCEDURES START HERE
!================================================================================!
!================================================================================!

  subroutine lopt(init,n,no,accr,op,d)
!********************************************************************
!*  Original subroutine for Boys localization by R.Ahlrichs 2/85
!*
!* to determine vectors |i>, i=1,n
!* which maximize
!*   f = Σ(i,j=1,n)Σ(k=1,no) (<i|op(k)|i>-<j|op(k)|j>)**2
!* equivalent to maximization of
!*   g = Σ(i=1,n)Σ(k=1,no) <i|op(k)|i>**2 + <j|op(k)|j>**2
!* equivalent to minimization of
!*   h = Σ(i<j)Σ(k) <i|op(k)|j>**2
!* stationarity condition
!*   Σ(k) op(ij,k)*(op(ii,k)-op(jj,k)) = 0 all i.ne.j
!*
!* input
!*       n = dimension of lin. space
!*       no = # of operators in f or g
!*       op(ij,k) matrix representation of k'th operator
!*            symmetric and stored triangular,nn1=n(n+1)/2
!*       thi = rotation threshold (1.-6 to 1.-12 sngl or dbl prc)
!*       nboys = max. nb. of sweeps
!* result
!*        d  - contains the solution vectors as columns
!*        op - the transformed operators
!*
!*******************************************************************
    implicit real(wp) (a-h,o-z)
    !> INPUT
    logical,intent(in) :: init
    integer,intent(in) :: n,no
    real(wp),intent(in) :: accr
    !> OUTPUT
    real(wp),intent(inout) :: op(n*(n+1)/2,no)
    real(wp),intent(out) :: d(n,n)
    !> LOCAL
    integer :: auxi
    real(wp) :: auxvec(n)

    a0 = 0.d0
    a1 = 1.d0
    a2 = 2.d0
    a10 = 10.d0
    aqu = 0.25d0
    d(:,:) = 0.0_wp
!>--- accr=1.d-6 ! originally 10-14 but this setting saves a factor of 3 !
    thi = 1.d-10
    th = dmax1(thi,accr)

    first = 0

!>--- set d to unit matrix
    if (init) then
      if (debug) write (*,*) 'initialization of trafo matrix to unity'
      do i = 1,n
        do j = 1,n
10        d(j,i) = a0-accr*i+accr*i*j
        end do
20      d(i,i) = a1
      end do
    end if

    if (n .eq. 1) return

    nboys = 20000
    nn1 = n*(n+1)/2

!>---  opn = operator norm needed for cut off threshold
    opn = a0
    do k = 1,no
      do i = 1,nn1
21      opn = opn+op(i,k)**2
      end do
22    continue
    end do
    if (opn .eq. a0) return
    opn = dsqrt(opn/dble(nn1*no))*accr

    ths = aqu
    it = 0
    do isw = 1,nboys
!>--- loop over sweeps
      thsw = a0
      ii = 1
      ij = 2
      do i = 2,n
        ii = ii+i
        jj = 0
        im1 = i-1
        ip1 = i+1
        li = ii-i
        do j = 1,im1
!>--- loop over i,j rotations
!>    maximize sum(k) <i|op(k)|i>**2 + <j|op(k)|j>**2
!>    like in jacobi
          jj = jj+j
!>--- get ax , b determining the tan of the rotation angle
          ax = a0
          b = a0
          do k = 1,no
            aux1 = op(jj,k)-op(ii,k)
            aux2 = op(ij,k)
            ax = ax+aux1*aux2
30          b = b+aqu*aux1**2-aux2**2
          end do
!>--- rotation angle x0
          if (dabs(ax)+dabs(b) .lt. opn) go to 100
          x0 = aqu*datan2(ax,b)
          thsw = dmax1(thsw,dabs(x0))
          if (dabs(x0) .lt. ths) go to 100
!>--- rotation
          it = it+1
          s = dsin(x0)
          s2 = s*s
          c2 = a1-s2
          c = dsqrt(c2)
          do l = 1,n
            aux1 = d(l,j)
            d(l,j) = c*d(l,j)+s*d(l,i)
            d(l,i) = c*d(l,i)-s*aux1
40          continue
          end do
          jm1 = j-1
          cs = c*s
          c2ms2 = c2-s2
          lj = jj-j
          jp1 = j+1

          ! !$OMP  PARALLEL PRIVATE (k,l,aux1,aux2,il,jl) &
          ! !$OMP& SHARED(op) &
          ! !$OMP& DEFAULT(SHARED)
          ! !$OMP DO
          do k = 1,no
            aux1 = op(ii,k)
            aux2 = op(jj,k)
            op(ii,k) = aux1*c2+aux2*s2-a2*op(ij,k)*cs
            op(jj,k) = aux1*s2+aux2*c2+a2*op(ij,k)*cs
            op(ij,k) = cs*(aux1-aux2)+op(ij,k)*c2ms2
            if (jm1 .le. 0) go to 60

            !                   do l=1,jm1
            !                       aux1=c*op(li+l,k)-s*op(lj+l,k)
            !                       op(lj+l,k)=s*op(li+l,k)+c*op(lj+l,k)
            !                       op(li+l,k)=aux1
            !                   end do
            auxvec(1:jm1) = c*op(li+1:li+jm1,k)-s*op(lj+1:lj+jm1,k)
            op(lj+1:lj+jm1,k) = op(lj+1:lj+jm1,k)*c
            op(lj+1:lj+jm1,k) = op(lj+1:lj+jm1,k)+s*op(li+1:li+jm1,k)
            op(li+1:li+jm1,k) = auxvec(1:jm1)

60          if (im1 .eq. j) go to 80
            jl = jj+j

            do l = jp1,im1
              aux1 = c*op(li+l,k)-s*op(jl,k)
              op(jl,k) = c*op(jl,k)+s*op(li+l,k)
              op(li+l,k) = aux1
70            jl = jl+l
            end do

80          if (i .eq. n) go to 95
            jl = ii+j
            il = ii+i

            do l = ip1,n
              aux1 = c*op(jl,k)+s*op(il,k)
              op(il,k) = c*op(il,k)-s*op(jl,k)
              op(jl,k) = aux1
              il = il+l
90            jl = jl+l
            end do

95          continue
          end do
          ! !$OMP END DO
          ! !$OMP END PARALLEL

100       ij = ij+1
        end do
110     ij = ij+1
      end do

!>--- check convergence
      if (thsw .le. th) go to 210
      ths = dmin1(thsw**2,ths*aqu)
      ths = dmax1(th,ths)
      if (mod(isw,50) .eq. 0.AND.debug)&
      & write (*,'(''iteration'',i7,2x,''convergence'',d16.8)') isw,thsw
200   continue
    end do
    if (debug) write (*,300) isw,thsw
    return
210 if (debug) write (*,310) isw,thsw

    return
300 format(/' not converged in',i7,' iterations, threshold :',d16.8)
310 format(/' converged in',i7,' iterations, threshold : ',d16.8)
  end subroutine lopt

!========================================================================================!
!========================================================================================!
end module lopt_mod
