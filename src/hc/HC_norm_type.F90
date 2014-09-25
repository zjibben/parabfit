!!
!! HC_NORM_TYPE
!!
!! This module defines a derived type that encapsulates the norm used for
!! solution increments of the heat conduction model.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2003, July 2014
!!

#include "f90_assert.fpp"

module HC_norm_type

  use kinds, only: r8
  use HC_model_type
  implicit none
  private

  type, public :: HC_norm
    private
    type(HC_model), pointer :: model => null()  ! reference only -- do not own
    real(r8) :: abs_T_tol   ! absolute temperature tolerance
    real(r8) :: rel_T_tol   ! relative temperature tolerance
    real(r8) :: abs_H_tol   ! absolute enthalpy tolerance
    real(r8) :: rel_H_tol   ! relative enthalpy tolerance
  contains
    procedure :: init
    procedure :: compute
  end type HC_norm

contains

  subroutine init (this, model, params)

    use parameter_list_type
    use logging_services

    class(HC_norm), intent(out) :: this
    type(HC_model), intent(in), target :: model
    type(parameter_list) :: params
    
    integer :: stat
    character(:), allocatable :: context, errmsg

    this%model => model

    context = 'processing ' // params%name() // ': '
    call params%get ('temp-abs-tol', this%abs_T_tol, default=0.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    call params%get ('temp-rel-tol', this%rel_T_tol, default=0.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (this%abs_T_tol < 0.0_r8) call LS_fatal (context//'"temp-abs-tol" must be >= 0.0')
    if (this%rel_T_tol < 0.0_r8) call LS_fatal (context//'"temp-rel-tol" must be >= 0.0')
    if (this%abs_T_tol == 0.0_r8 .and. this%rel_T_tol == 0.0_r8) &
        call LS_fatal (context//'"temp-abs-tol" and "temp-rel-tol" cannot both be 0.0')

    call params%get ('enth-abs-tol', this%abs_H_tol, default=0.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    call params%get ('enth-rel-tol', this%rel_H_tol, default=0.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (this%abs_H_tol < 0.0_r8) call LS_fatal (context//'"enth-abs-tol" must be >= 0.0')
    if (this%rel_H_tol < 0.0_r8) call LS_fatal (context//'"enth-rel-tol" must be >= 0.0')
    if (this%abs_H_tol == 0.0_r8 .and. this%rel_H_tol == 0.0_r8) &
        call LS_fatal (context//'"enth-abs-tol" and "enth-rel-tol" cannot both be 0.0')

  end subroutine init


  subroutine compute (this, u, du, du_norm)

    class(HC_norm), intent(in) :: this
    real(r8), intent(in), target :: u(:), du(:)
    real(r8), intent(out) :: du_norm

    real(r8), pointer :: useg(:), duseg(:)

    ASSERT(size(u) == size(du))
    ASSERT(size(u) == this%model%num_dof())

    du_norm = 0.0_r8

    !! Cell temperature delta norm
    call this%model%get_cell_temp_view (u, useg)
    call this%model%get_cell_temp_view (du, duseg)
    du_norm = max(du_norm, maxerr(useg, duseg, this%abs_T_tol, this%rel_T_tol))

    !! Face temperature delta norm
    call this%model%get_face_temp_view (u, useg)
    call this%model%get_face_temp_view (du, duseg)
    du_norm = max(du_norm, maxerr(useg, duseg, this%abs_T_tol, this%rel_T_tol))

    !! Cell enthalpy delta norm
    call this%model%get_cell_heat_view (u, useg)
    call this%model%get_cell_heat_view (du, duseg)
    du_norm = max(du_norm, maxerr(useg, duseg, this%abs_T_tol, this%rel_T_tol))

  contains

    real(r8) function maxerr (u, du, atol, rtol)
      real(r8), intent(in) :: u(:), du(:), atol, rtol
      maxerr = maxval(abs(du)/(atol + rtol*abs(u)))
    end function

  end subroutine compute

end module HC_norm_type
