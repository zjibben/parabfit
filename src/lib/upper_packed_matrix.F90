!!
!! UPPER_PACKED_MATRIX
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, March 2014.  Original code 10 Jan 2006.
!!

#include "f90_assert.fpp"

module upper_packed_matrix

  use kinds, only: r8
  implicit none
  private

  public :: upm_factor, upm_invert, upm_solve, upm_matmul, upm_col_sum

  interface upm_matmul
    procedure upm_matmul_vec, upm_matmul_mat
  end interface

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  UPM_INVERT
 !!
 !!  This subroutine computes the inverse of the symmetric, positive-definite
 !!  matrix A stored in upper packed storage format.  A is overwritten with its
 !! (symmetric) inverse, also in upper packed format.
 !!

  pure subroutine upm_invert (a)

    real(r8), intent(inout) :: a(:)

    integer :: l, m, p, q
    real(r8) :: s

    if (size(a) == 0) return

    !! Compute the upper Choleski factor.
    call upm_factor (a)

    !! Compute the inverse of the Choleski factor.
    a(1) = 1.0_r8 / a(1)
    p = 2
    do while (p <= size(a))
      l = 1
      q = p
      do while (l < p)
        s = a(q)
        m = p
        do while (m < q)
          a(m) = a(m) + s*a(l)
          l = l + 1
          m = m + 1
        end do
        a(m) = s*a(l)
        l = l + 1
        q = q + 1
      end do
      s = 1.0_r8 / a(q)
      a(q) = s
      s = -s
      m = p
      do while (m < q)
        a(m) = s*a(m)
        m = m + 1
      end do
      p = q + 1
    end do

    !! Compute the product of the inverse Choleski factors.
    a(1) = a(1)*a(1)
    p = 2
    do while (p <= size(a))
      l = 1
      q = p
      do while (l < p)
        s = a(q)
        m = p
        do while (m <= q)
          a(l) = a(l) + s*a(m)
          l = l + 1
          m = m + 1
        end do
        q = q + 1
      end do
      s = a(q)
      do while (l <= q)
        a(l) = s*a(l)
        l = l + 1
      end do
      p = q + 1
    end do

  end subroutine upm_invert

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! UPM_FACTOR
 !!
 !! Compute the Choleski factorization of a symmetric positive definite matrix.
 !! The matrix is provided in upper packed storage format by the argument A,
 !! and the upper Choleski factor (also in upper packed storage) is returned
 !! in A, overwriting the original matrix.
 !!

  pure subroutine upm_factor (a)

    real(r8), intent(inout) :: a(:)

    integer :: l, m, p, q
    real(r8) :: s, t

    if (size(a) == 0) return

    !! Compute the upper Choleski factor.
    a(1) = sqrt(a(1))
    p = 2
    do while (p <= size(a))
      l = 1
      q = p
      t = 0.0_r8
      do while (l < p)
        s = a(q)
        m = p
        do while (m < q)
          s = s - a(l)*a(m)
          l = l + 1
          m = m + 1
        end do
        s = s / a(l)
        a(q) = s
        t = t + s**2
        l = l + 1
        q = q + 1
        if (q > size(a)) s = sqrt(real(size(a) - q))  ! trigger an exception
      end do
      s = a(q) - t
      a(q) = sqrt(s)
      p = q + 1
    end do

  end subroutine upm_factor

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! UPM_SOLVE
 !!
 !! Solve the linear system A*x = b, where A is a symmetric matrix.  The
 !! argument A holds the upper Choleski factor (as computed by FACTOR_UPM)
 !! in upper packed storage format.  The rhs vector is provided in X, and
 !! is overwritten with the solution.
 !!

  pure subroutine upm_solve (a, x)

    real(r8), intent(in) :: a(:)
    real(r8), intent(inout) :: x(:)

    integer :: i, j, l
    real(r8) :: s

    if (size(x) == 0) return

    !! Forward substitution.
    x(1) = x(1)/a(1)
    l = 2
    do i = 2, size(x)
      s = x(i)
      do j = 1, i-1
        s = s - a(l)*x(j)
        l = l + 1
      end do
      x(i) = s/a(l)
      l = l + 1
    end do

    !! Backward substitution.
    l = l - 1
    do i = size(x), 2, -1
      x(i) = x(i)/a(l)
      l = l - 1
      do j = i-1, 1, -1
        x(j) = x(j) - a(l)*x(i)
        l = l - 1
      end do
    end do
    x(1) = x(1)/a(1)

  end subroutine upm_solve

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedures for UPM_MATMUL
 !!
 !! These functions compute a matrix-vector and matrix-matrix product when
 !! the first argument (and left factor) is a symmetric matrix stored in
 !! upper-packed storage format.  The result is either a rank-1 (vector) or
 !! rank-2 (matrix) array as appropriate.
 !!

  pure function upm_matmul_vec (a, b) result (c)

    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: c(size(b))

    integer :: i, j, l
    real(r8) :: bi, ci

    l = 1
    do i = 1, size(b)
      bi = b(i)
      ci = 0.0_r8
      do j = 1, i-1
        ci = ci + a(l)*b(j)
        c(j) = c(j) + a(l) * bi
        l = l + 1
      end do
      c(i) = ci + a(l)*bi
      l = l + 1
    end do

  end function upm_matmul_vec

  pure function upm_matmul_mat (a, b) result (c)

    real(r8), intent(in) :: a(:), b(:,:)
    real(r8) :: c(size(b,1),size(b,2))

    integer :: i, j, k, l
    real(r8) :: bi, ci

    do k = 1, size(b,2)
      l = 1
      do i = 1, size(b,1)
        bi = b(i,k)
        ci = 0.0_r8
        do j = 1, i-1
          ci = ci + a(l)*b(j,k)
          c(j,k) = c(j,k) + a(l) * bi
          l = l + 1
        end do
        c(i,k) = ci + a(l)*bi
        l = l + 1
      end do
    end do

  end function upm_matmul_mat

  subroutine upm_col_sum (a, csum)

    real(r8), intent(in)  :: a(:)
    real(r8), intent(out) :: csum(:)

    integer :: i, j, l
    real(r8) :: s

    ASSERT(size(a) == (size(csum)*(size(csum)+1))/2)

    l = 1
    do i = 1, size(csum)
      s = 0.0_r8
      do j = 1, i-1
        s = s + a(l)
        csum(j) = csum(j) + a(l)
        l = l + 1
      end do
      csum(i) = s + a(l)
      l = l + 1
    end do

  end subroutine upm_col_sum

end module upper_packed_matrix
