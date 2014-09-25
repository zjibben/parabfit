!!
!! HC_MODEL_TYPE
!!
!! This module defines a class that encapsulates the ODE system arising from
!! the mimetic finite difference discretization of the heat conduction equation
!! over an unstructured 3D mesh.  This is a basic version of the linear heat
!! equation with constant coefficients, and fully general Dirichlet and flux
!! boundary conditions.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, April 2014
!!

#include "f90_assert.fpp"

module HC_model_type

  use kinds, only: r8
  use unstr_mesh_type
  use mfd_disc_type
  use data_layout_type
  use bndry_func_class
  use timer_tree_type
  implicit none
  private

  type, public :: HC_model
    type(mfd_disc),   pointer :: disc => null() ! reference only -- do not own
    type(unstr_mesh), pointer :: mesh => null() ! reference only -- do not own
    !! Variable layout
    type(data_layout) :: layout
    integer :: cell_heat_segid, cell_temp_segid, face_temp_segid
    !! Boundary condition data
    class(bndry_func), allocatable :: temp_bc, flux_bc
    !! Constant model parameters (simple first step)
    real(r8) :: rho, cp, kappa, q
  contains
    procedure :: init
    procedure :: num_dof
    procedure :: H_of_T
    procedure :: T_of_H
    procedure :: dHdT
    procedure :: conductivity
    procedure :: source
    procedure :: get_cell_heat_view
    procedure :: get_cell_temp_view
    procedure :: get_face_temp_view
    procedure :: get_cell_heat_copy
    procedure :: get_cell_temp_copy
    procedure :: get_face_temp_copy
    procedure :: set_cell_heat
    procedure :: set_cell_temp
    procedure :: set_face_temp
    procedure :: residual
  end type HC_model

contains

  subroutine init (this, disc, params)

    use parameter_list_type
    use bc_factory_type
    use logging_services

    class(HC_model), intent(out) :: this
    type(mfd_disc), intent(in), target :: disc
    type(parameter_list) :: params

    integer :: stat
    character(:), allocatable :: context, errmsg
    type(bc_factory) :: bcfact

    this%disc => disc
    this%mesh => disc%mesh

    !! Create the packed layout of the model variables.
    this%cell_heat_segid = this%layout%alloc_segment(this%mesh%ncell)
    this%cell_temp_segid = this%layout%alloc_segment(this%mesh%ncell)
    this%face_temp_segid = this%layout%alloc_segment(this%mesh%nface)
    call this%layout%alloc_complete

    !! Model parameters
    context = 'processing ' // params%name() // ': '
    call params%get ('density', this%rho, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (this%rho <= 0.0_r8) call LS_fatal (context//'"density" must be > 0.0')
    call params%get ('specific-heat', this%cp, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (this%cp <= 0.0_r8) call LS_fatal (context//'"specific-heat" must be > 0.0')
    call params%get ('conductivity', this%kappa, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (this%kappa <= 0.0_r8) call LS_fatal (context//'"conductivity" must be > 0.0')
    call params%get ('source', this%q, default=0.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)

    !! Initialize the boundary condition components
    if (params%is_sublist('bc')) then
      call bcfact%init (this%mesh, params%sublist('bc'))
      call bcfact%alloc_bc ('dirichlet', this%temp_bc)
      call bcfact%alloc_bc ('flux', this%flux_bc)
    else
      call LS_fatal (context//'missing "bc" sublist parameter')
    end if

  end subroutine init

  integer function num_dof (this)
    class(HC_model), intent(in) :: this
    num_dof = this%layout%layout_size()
  end function num_dof

  subroutine H_of_T (this, temp, value)
    class(HC_model), intent(in) :: this
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: value(:)
    value = (this%rho * this%cp) * temp
  end subroutine H_of_T

  subroutine T_of_H (this, enth, value)
    class(HC_model), intent(in) :: this
    real(r8), intent(in)  :: enth(:)
    real(r8), intent(out) :: value(:)
    value = enth / (this%rho * this%cp)
  end subroutine T_of_H

  subroutine dHdT (this, temp, value)
    class(HC_model), intent(in) :: this
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: value(:)
    value = (this%rho * this%cp)
  end subroutine dHdT

  subroutine conductivity (this, temp, value)
    class(HC_model), intent(in) :: this
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: value(:)
    value = this%kappa
  end subroutine conductivity

  subroutine source (this, t, value)
    class(HC_model), intent(in) :: this
    real(r8), intent(in)  :: t
    real(r8), intent(out) :: value(:)
    value = this%q
  end subroutine source


  subroutine get_cell_heat_view (this, array, view)
    class(HC_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    call this%layout%segment_view (array, this%cell_heat_segid, view)
  end subroutine get_cell_heat_view

  subroutine get_cell_temp_view (this, array, view)
    class(HC_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    call this%layout%segment_view (array, this%cell_temp_segid, view)
  end subroutine get_cell_temp_view

  subroutine get_face_temp_view (this, array, view)
    class(HC_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    call this%layout%segment_view (array, this%face_temp_segid, view)
  end subroutine get_face_temp_view

  subroutine get_cell_heat_copy (this, array, copy)
    class(HC_model), intent(in) :: this
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    call this%layout%segment_copy (array, this%cell_heat_segid, copy)
  end subroutine get_cell_heat_copy

  subroutine get_cell_temp_copy (this, array, copy)
    class(HC_model), intent(in) :: this
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    call this%layout%segment_copy (array, this%cell_temp_segid, copy)
  end subroutine get_cell_temp_copy

  subroutine get_face_temp_copy (this, array, copy)
    class(HC_model), intent(in) :: this
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    call this%layout%segment_copy (array, this%face_temp_segid, copy)
  end subroutine get_face_temp_copy

  subroutine set_cell_heat (this, source, array)
    class(HC_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout), target :: array(:)
    real(r8), pointer :: view(:)
    call this%layout%segment_view (array, this%cell_heat_segid, view)
    view = source(:size(view))
  end subroutine set_cell_heat

  subroutine set_cell_temp (this, source, array)
    class(HC_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout), target :: array(:)
    real(r8), pointer :: view(:)
    call this%layout%segment_view (array, this%cell_temp_segid, view)
    view = source(:size(view))
  end subroutine set_cell_temp

  subroutine set_face_temp (this, source, array)
    class(HC_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout), target :: array(:)
    real(r8), pointer :: view(:)
    call this%layout%segment_view (array, this%face_temp_segid, view)
    view = source(:size(view))
  end subroutine set_face_temp


  subroutine residual (this, t, u, udot, f)

    class(HC_model), intent(inout) :: this
    real(r8), intent(in)  :: t, u(:), udot(:)
    real(r8), intent(out) :: f(:)
    target :: u, udot, f

    real(r8), pointer :: Hcell(:), Tcell(:), Fcell(:), Fface(:), Hdot(:)
    real(r8) :: Tface(this%mesh%nface), cval(this%mesh%ncell)
    real(r8), allocatable :: Fdir(:)

    call start_timer ('hc-residual')

  !!!! RESIDUAL OF THE ALGEBRAIC ENTHALPY-TEMPERATURE RELATION !!!!!!!!!!!!!!!!!

    call this%get_cell_heat_view (u, Hcell)
    call this%get_cell_temp_view (u, Tcell)
    call this%get_cell_heat_view (f, Fcell)
    call this%H_of_T (Tcell, cval)
    Fcell = Hcell - cval

  !!!! RESIDUAL OF THE HEAT EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Pre-compute the Dirichlet condition residual and impose the Dirichlet data.
    call this%get_face_temp_copy (u, Tface) ! NB: copy, will modify boundary
    call this%temp_bc%compute (t)
    associate (index => this%temp_bc%index, value => this%temp_bc%value)
      Fdir = Tface(index) - value
      Tface(index) = value
    end associate

    !! Compute the generic heat equation residual.
    call this%get_cell_temp_view (f, Fcell)
    call this%get_face_temp_view (f, Fface)
    call this%conductivity (Tcell, cval)
    call this%disc%apply_diff (cval, Tcell, Tface, Fcell, Fface)
    call this%source (t, cval)
    call this%get_cell_heat_view (udot, Hdot)
    Fcell = Fcell + this%mesh%volume*(Hdot - cval)

    !! Dirichlet condition residuals.
    associate (faces => this%temp_bc%index)
      Fface(faces) = Fdir ! overwrite with pre-computed values
    end associate

    !! Simple flux BC contribution.
    call this%flux_bc%compute (t)
    associate (index => this%flux_bc%index, value => this%flux_bc%value)
      Fface(index) = Fface(index) + this%mesh%area(index) * value
    end associate

    call stop_timer ('hc-residual')

  end subroutine residual

end module HC_model_type
