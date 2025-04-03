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

module lbfgs_module
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use crest_calculator
  use strucrd
  use optimize_type  !> This module provides the 'optimizer' type.
  implicit none
  private

  public :: lbfgs_optimize

!========================================================================================!
!========================================================================================!
contains !>  MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  function lbfgs_direction(nvar,g,k,OPT,gamma) result(d)
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
    !* @param gamma  Scaling factor for the initial Hessian approximation.
    !*
    !* @return d     Computed search direction (negative approximate inverse
    !*               Hessian times g).
    !********************************************************************************
    !> INPUT
    integer,intent(in) :: nvar,k
    type(optimizer),intent(inout) :: OPT
    real(wp),intent(in) :: g(nvar)
    real(wp),intent(in) :: gamma
    !> OUTPUT
    real(wp) :: d(nvar)
    !> LOCAL
    integer :: i

    associate (S => OPT%S,Y => OPT%Y,alpha => OPT%alpha,rho => OPT%rho,q => OPT%q)

      !> Initialize q with the current gradient.
      q = g

      !---------------------------------------------------------
      !> First loop (backward pass): for i = k downto 1,
      !> compute the coefficient alpha(i) and update q.
      !---------------------------------------------------------
      do i = k,1,-1
        alpha(i) = rho(i)*dot_product(S(1:nvar,i),q)
        q = q-alpha(i)*Y(1:nvar,i)
      end do

      !---------------------------------------------------------
      !> Apply the initial Hessian approximation.
      !> We use a scaled identity matrix H0 = gamma * I.
      !---------------------------------------------------------
      d = gamma*q

      !---------------------------------------------------------
      !> Second loop (forward pass): for i = 1 to k,
      !> compute the correction and update d.
      !---------------------------------------------------------
      do i = 1,k
        d = d+S(1:nvar,i)*(alpha(i)-rho(i)*dot_product(Y(1:nvar,i),d))
      end do

      !> The final search direction is the negative of d.
      d = -d

    end associate
  end function lbfgs_direction

  subroutine obj_func(x,f,g)
    !************************************************************************
    !* Dummy objective function (quadratic).
    !*
    !* This subroutine computes the value and gradient of the objective
    !* function defined as f(x) = 1/2 * x^T * x. Its gradient is simply x.
    !* Replace this with your actual function evaluations.
    !*
    !* @param x   Input variable vector.
    !* @param f   Output function value.
    !* @param g   Output gradient vector.
    !***********************************************************************
    real(wp),intent(in)  :: x(:)
    real(wp),intent(out) :: f
    real(wp),intent(out) :: g(size(x))
    f = 0.5_wp*dot_product(x,x)
    g = x
  end subroutine obj_func

  subroutine lbfgs_optimize(mol,calc,etot,grd,nvar,x,max_iter,m,tol,pr,io)
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
    !* @param nvar     Integer. Dimension of the variable space.
    !* @param x        Real(wp) array. Input coordinate vector; updated with optimized values.
    !* @param max_iter Integer. Maximum number of iterations allowed.
    !* @param m        Integer. Maximum number of stored corrections (history length).
    !* @param tol      Real(wp). Convergence tolerance (stopping criteria based on f).
    !* @param io       Integer. Output status variable (0 indicates success).
    !**************************************************************************
    implicit none
    !> INPUT
    type(coord),intent(inout) :: mol
    type(calcdata),intent(in) :: calc
    real(wp),intent(inout) :: etot
    real(wp),intent(inout) :: grd(3,mol%nat)
    integer,intent(in)     :: nvar
    real(wp),intent(inout) :: x(nvar)
    integer,intent(in)     :: max_iter,m
    real(wp),intent(in)    :: tol
    logical,intent(in)     :: pr
    !> OUTPUT
    integer,intent(out)    :: io
    !> LOCAL
    type(optimizer) :: OPT
    integer :: iter,k
    real(wp) :: gnorm,deltaE
    real(wp),allocatable :: g(:),d(:),g_new(:),x_new(:)
    real(wp) :: f,f_new,gamma,step

    !> Prepare settings
    io = 0


    !> Allocate the vectors for position, gradient, and search direction.
    allocate (g(nvar),d(nvar),g_new(nvar),x_new(nvar))
    !> Allocate matrices to store up to m correction pairs (columns correspond to each stored pair).
    call OPT%allocatelbfgs(nvar,m)
    associate (S => OPT%S,Y => OPT%Y,rho => OPT%rho,alpha => OPT%alpha)
      S = 0.0_wp
      Y = 0.0_wp
      rho = 0.0_wp
      k = 0  ! Initially, no correction pairs are stored.

      !> Evaluate the objective function and its gradient at the starting point.
      call obj_func(x,f,g)

      iter = 0
      if (pr) call print_optiter(iter)
      gnorm = norm2(grd)
      if (pr) then
        write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') etot
        write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') 0.0_wp
        write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")') gnorm
      end if

      do while (f > tol.and.iter < max_iter)
        iter = iter+1
        if (pr) call print_optiter(iter) 

        if (iter == 1) then
          !> First iteration: use the steepest descent direction.
          d = -g
        else
          !---------------------------------------------------------
          !> Determine the scaling factor gamma for the initial Hessian.
          !  Here we use the most recent correction pair.
          !---------------------------------------------------------
          if (k > 0) then
            gamma = dot_product(S(1:nvar,k),Y(1:nvar,k))/ &
                    dot_product(Y(1:nvar,k),Y(1:nvar,k))
          else
            gamma = 1.0_wp
          end if

          !> Compute the search direction using the two-loop recursion.
          d = lbfgs_direction(nvar,g,k,OPT,gamma)
        end if

        !---------------------------------------------------------
        !> A fixed step size could be used here for simplicity.
        !  In a full implementation, a line search could be used.
        !---------------------------------------------------------
        step = 1.0_wp

        !> Update the position: x_new = x + step * d.
        x_new = x+step*d

        !====================================================================!
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
        !> Evaluate the objective and gradient at the new position.
        call obj_func(x_new,f_new,g_new)
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
        !====================================================================!
        if (io /= 0) exit
        gnorm = norm2(grd) 

        !---------------------------------------------------------
        !> Compute the correction pair:
        !   s = x_new - x and y = g_new - g.
        !---------------------------------------------------------
        if (k < m) then
          !> If there is still room in the history, simply add the new pair.
          k = k+1
          S(1:nvar,k) = x_new-x
          Y(1:nvar,k) = g_new-g
        else
          !> When the history is full, shift the stored corrections and
          !> insert the new pair at the end.
          S(1:nvar,1:m-1) = S(1:nvar,2:m)
          Y(1:nvar,1:m-1) = Y(1:nvar,2:m)
          S(1:nvar,m) = x_new-x
          Y(1:nvar,m) = g_new-g
        end if

        !> Update the scaling factor for the new correction pair.
        rho(k) = 1.0_wp/dot_product(Y(1:nvar,k),S(1:nvar,k))

        !> Update the current position, gradient, and function value.
        x = x_new
        g = g_new
        f = f_new

        !> Optional: print iteration information.
        if (pr) then
          write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') etot
          write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') deltaE
          write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no') gnorm
        end if
      end do

      !> stop associating
    end associate

    !> Deallocate all temporary arrays.
    deallocate (g,d,g_new,x_new)
  end subroutine lbfgs_optimize

end module lbfgs_module

