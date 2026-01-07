!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2025 Philipp Pracht
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

!> This module implemnts a simple L-BFGS with different coordinate choices

module newton_raphson_module
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use crest_calculator
  use strucrd
  use optimize_type  !> This module provides the 'optimizer' type.
  use optimize_utils
  use coordinate_transform_module
  implicit none
  private

  public :: newton_raphson

!========================================================================================!
!========================================================================================!
contains !>  MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  function lbfgs_direction(nvar,g,k,OPT) result(d)
    !*******************************************************************************
    !* Two-loop recursion routine to compute the search direction.
    !*
    !* This function uses the stored correction pairs (S and Y) and
    !* corresponding scaling factors (rho) to approximate the product
    !* H*g, where H is the inverse Hessian approximation.
    !*
    !* The algorithm proceeds in two loops:
    !* 1. The first (backward) loop computes the coefficients "alpha" and
    !*    subtracts corrections from the gradient.
    !* 2. The second (forward) loop applies the corrections in the reverse order.
    !* Finally, the result is negated to obtain the descent direction.
    !* Note, this routine is NOT called for the very first iteration (k == 0)
    !*
    !* @param nvar   Dimension of the variable space.
    !* @param k      Number of stored correction pairs (k ≤ m).
    !* @param g      Current gradient vector.
    !* @param OPT    optimizer type that stores the following variables
    !*     => S      Matrix containing the s-vectors (x_{k+1} - x_k),
    !*               of size(nvar, m).
    !*     => Y      Matrix containing the y-vectors (g_{k+1} - g_k),
    !*               of size(nvar, m).
    !*     => rho    Array of stored values 1/(y^T*s) for each correction.
    !*     => alpha  coefficients (get computed in this function)
    !*     => q      temporary workspace
    !* @param gamm  Scaling factor for the initial Hessian approximation.
    !*
    !* @return d     Computed search direction (negative approximate inverse
    !*               Hessian times g).
    !********************************************************************************
    !> INPUT
    integer,intent(in) :: nvar,k
    type(optimizer),intent(inout) :: OPT
    real(wp),intent(in) :: g(nvar)
    !> OUTPUT
    real(wp) :: d(nvar)
    !> LOCAL
    integer :: i
    external :: dppsv
    external :: dspsv
    real(wp) :: c(nvar)
    integer :: IPIV(nvar)
    integer :: info2,info3
    real(wp), allocatable :: int_hess(:)

    if (.not. allocated(int_hess)) then
      allocate(int_hess(size(OPT%hess)))
    endif
    int_hess = OPT%hess
    c = -g !> This will become the update step after system is solved
    call dppsv('L', nvar, 1, int_hess, c, nvar, info2)
    if (info2 /= 0) then
        int_hess = OPT%hess
        c = -g
        call dspsv('L', nvar, 1, int_hess, IPIV, c, nvar, info3)
        d = c  
    else
        d = c
    endif

  end function lbfgs_direction

!========================================================================================!

  subroutine newton_raphson(mol,calc,etot,grd,pr,io)
    !**************************************************************************
    !* L-BFGS Optimization Routine
    !*
    !* Performs optimization using the Limited-memory BFGS (L-BFGS) algorithm.
    !* The routine updates the coordinate vector x to approach a local minimum of the
    !* objective function. It integrates with an optimizer type (OPT) to manage the
    !* correction pairs (s and y) and related internal data using associate constructs.
    !*
    !* The main steps include:
    !*  1. Evaluating the objective function and gradient at the current x.
    !*  2. Computing the search direction via the two-loop recursion (lbfgs_direction).
    !*  3. Updating x using a fixed step (with the option to incorporate a line search).
    !*  4. Updating the correction pairs: s = x_new - x and y = g_new - g, while managing
    !*     the history using a shifting strategy when full.
    !*
    !* @param io       Integer. Output status variable (0 indicates success).
    !**************************************************************************
    implicit none
    !> INPUT
    type(coord),intent(inout) :: mol
    type(calcdata),intent(in) :: calc
    real(wp),intent(inout) :: etot
    real(wp),intent(inout) :: grd(3,mol%nat)
    logical,intent(in)     :: pr
    !> OUTPUT
    integer,intent(out)    :: io
    !> LOCAL
    type(optimizer) :: OPT
    integer :: iter,k,nvar,m
    integer :: tight,max_iter
    real(wp) :: gnorm,deltaE,energy
    real(wp) :: ethr,gthr,maxerise
    logical :: econverged,gconverged,converged,Erise
    real(wp),allocatable :: x(:),g(:),d(:),g_new(:),x_new(:),gtmp(:,:)
    real(wp) :: f,f_new,gamm,step
    integer :: ilog

    !> Prepare settings
    io = 0
    nvar = compute_nvar(mol)
    m = calc%lbfgs_histsize
    gnorm = norm2(grd)
    deltaE = huge(deltaE)
    tight = calc%optlev
    call get_optthr(mol%nat,tight,calc,ethr,gthr)
    max_iter = calc%maxcycle !> automatic setting in get_optthr or by user
    maxerise = calc%maxerise
    econverged = .false.
    gconverged = .false.
    converged = .false.

    open (newunit=ilog,file='crestopt.log.xyz')
    call mol%appendlog(ilog,etot)

    !$omp critical
    !> Allocate the vectors for position, gradient, and search direction.
    allocate (x(nvar),g(nvar),d(nvar),g_new(nvar),x_new(nvar),gtmp(3,mol%nat))
    !> Allocate matrices to store up to m correction pairs (columns correspond to each stored pair).
    call OPT%allocate2(mol%nat)
    !$omp end critical

      !> First trafo
      call transform_mol('cart2v',mol,nvar,x)
      call transform_grd('cart2v',mol,grd,nvar,g)

      iter = 0
      if (pr) then
        call print_optiter(iter)
        write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') etot
        write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') 0.0_wp
        write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")') gnorm
      end if

      LBFGS_iter: do while (.not.converged.and.iter < max_iter)
        iter = iter+1
        if (pr) call print_optiter(iter)

        if (iter == 1) then
          !> First iteration: use the steepest descent direction.
          d = -g
        else
          !> Compute the search direction using the two-loop recursion.
          d = lbfgs_direction(nvar,g,k,OPT)
        end if

        !---------------------------------------------------------
        !> A fixed step size could be used here for simplicity.
        !  In a full implementation, a line search could be used.
        !  If the energy rises, we reduce the stepsize iteratively
        !---------------------------------------------------------
        step = 0.0008_wp

        Erise = .true.
        do while (Erise)
          !> Update the position: x_new = x + step * d.
          x_new = x+step*d

          !====================================================================!
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
          !> Evaluate the objective and gradient at the new position.
          call transform_mol('v2cart',mol,nvar,x_new)
          grd = 0.0_wp
          call engrad(mol,calc,energy,gtmp,io)
          call mol%appendlog(ilog,energy)
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
          !====================================================================!
          if (io /= 0) exit

          call transform_grd('cart2v',mol,gtmp,nvar,g_new)
          gnorm = norm2(gtmp)
          deltaE = energy-etot

          Erise = (deltaE > maxerise)
          if (Erise) then
            step = step*0.25_wp
            if (pr) then
              write (*,'(" * energy rise detected, decreasing stepsize")')
            end if
          end if
        end do
        econverged = abs(deltaE) .lt. ethr
        gconverged = gnorm .lt. gthr

        call bfgs(OPT%nvar,gnorm,g_new,g,d*step,OPT%hess)

        !> Update the current position, gradient, and function value.
        x = x_new
        g = g_new
        etot = energy

        !> Optional: print iteration information.
        if (pr) then
          write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') etot
          write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') deltaE
          write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no') gnorm
          call print_convd(econverged,gconverged)
        end if
        converged = econverged.and.gconverged
      end do LBFGS_iter

    !> Final trafo
    call transform_mol('v2cart',mol,nvar,x_new)
    call transform_grd('v2cart',mol,grd,nvar,g_new)

    !> Deallocate all temporary arrays.
    deallocate (x_new,g_new,d,g,x)
  end subroutine newton_raphson

!========================================================================================!
!========================================================================================!
end module newton_raphson_module

