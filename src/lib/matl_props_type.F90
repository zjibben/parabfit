!!
!! MATERIAL_PROPERTIES_TYPE
!!
!! This module defines a type that encapsulates material data, such as density, sound speed,
!! immobility, etc, for all materials needed by the flow solver. For now it is a simple
!! collection of arrays giving the properties for each material, which are read from a
!! parameter list with the quantities essentially written out verbatim, but in the future
!! it could be expanded to include things like..
!!
!! TODO:
!!   * temperature dependence
!!   * read in data from a material library, given a list of materials needed
!!
!! Zechariah J Jibben <zjibben@lanl.gov>
!! February 2016
!! 

module matl_props_type

  use kinds, only: r8
  implicit none
  private

  type, public :: matl_props
    integer,  public :: nmat
    real(r8), public, allocatable :: density(:), sound_speed(:)
    logical,  public, allocatable :: is_immobile(:), is_void(:)
  contains
    procedure :: init
  end type matl_props

contains

  subroutine init (this, params)

    use parameter_list_type
    use logging_services

    class(matl_props), intent(out) :: this
    type(parameter_list) :: params

    character(:), allocatable :: context,errmsg
    integer, allocatable :: immobile_materials(:)
    integer :: stat

    context = 'processing ' // params%name() // ': '

    !! get material densities
    call params%get ('densities', this%density, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)

    this%nmat = size(this%density)

    !! identify void materials
    allocate(this%is_void(this%nmat))
    this%is_void = this%density==0.0_r8

    !! identify immobile materials
    allocate(this%is_immobile(this%nmat))
    this%is_immobile = .false.
    if (params%is_vector('immobile-materials')) then
      call params%get ('immobile-materials', immobile_materials, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      this%is_immobile(immobile_materials) = .true.
      deallocate(immobile_materials)
    end if

    !! get sound speeds
    allocate(this%sound_speed(this%nmat))
    this%sound_speed = 0.0_r8 ! default for non-void materials
    if (params%is_vector('sound-speeds')) then
      call params%get ('sound-speeds', this%sound_speed, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      ! TODO: some sort of check to ensure only void materials have a sound speed
    end if

  end subroutine init

end module matl_props_type
