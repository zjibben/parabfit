!!
!! MULTIMAT_CELL_TYPE
!!
!! This module provides a cell type that
!! describes internal material geometries
!! as child arbitrary polyhedra
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! Last revised 4 Nov 2012.
!!

module multimat_cell_type
  use kinds, only: r8
  use polyhedron_type
  implicit none
  private

  public :: multimat_cell_unit_test_suite

  ! a multimat_cell is a polyhedron itself, describing
  ! the cell geometry, and also contains an array
  ! of polyhedra each describing the geometry of a
  ! particular material
  type, extends(polyhedron), public :: multimat_cell
    integer                       :: nmat,m ! number of materials actually present in cell
    !integer,          allocatable :: mat_id(:)
    type(polyhedron), allocatable :: mat_poly(:)
  contains
    procedure :: partition
    procedure :: volumes_behind_plane
    procedure :: outward_volflux
    procedure :: interface_polygon
  end type multimat_cell
  
contains

  subroutine multimat_cell_unit_test_suite ()
    use hex_types, only: cube_v, hex_f, hex_e

    type(multimat_cell) :: cube
    logical             :: success
    real(r8)            :: vof(2), intnorm(3,2), posXflow(6), posxyzn(3), posxyn(3), posxn(3), &
        posyn(3), tmp

    call cube%init (cube_v, hex_f, hex_e)
    
    ! define face velocities [+y, -y, -x, +x, -z, +z]
    posXflow = [ 0.0_r8, 0.0_r8, -1.0_r8, 1.0_r8, 0.0_r8, 0.0_r8 ]

    ! define normals
    posxyzn = [1.0_r8, 1.0_r8, 1.0_r8] / sqrt(3.0_r8) ! positive x-y-z
    posxn   = [1.0_r8, 0.0_r8, 0.0_r8]                ! positive x-direction
    posyn   = [0.0_r8, 1.0_r8, 0.0_r8]                ! positive y-direction
    posxyn  = (posxn + posyn) / sqrt(2.0_r8) ! positive x-y plane
    
    write(*,*)
    write(*,*) 'MULTIMAT CELL TYPE'
    write(*,*) '===================================================='

    ! partitioning
    write(*,*) 'CELL PARTITIONING'
    ! TODO
    ! This is really already taken care of for two materials via
    ! the locate_plane_nd unit tests, but still needs to be done
    ! for more than two materials. Fluxing more than 2 materials
    ! also needs to be done.
    write(*,*) 'need to write tests here'

    ! fluxing
    write(*,*) 'FLUXING'

    ! cube 8/10ths filled in x direction
    vof = [0.8_r8, 0.2_r8]
    intnorm(:,1) = posxn; intnorm(:,2) = -intnorm(:,1)
    call cube%partition (vof, intnorm)

    success = fluxing_unit_test (cube, 0.25_r8*posXflow, reshape([&
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.05_r8, 0.2_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 8/10ths filled in x, fluxed in +x?  ', success

    success = fluxing_unit_test (cube, -0.25_r8*posXflow, reshape([&
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.25_r8, 0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 8/10ths filled in x, fluxed in -x?  ', success

    ! cube 8/10ths filled in -x direction
    vof = [0.8_r8, 0.2_r8]
    intnorm(:,1) = -posxn; intnorm(:,2) = -intnorm(:,1)
    call cube%partition (vof, intnorm)
    
    success = fluxing_unit_test (cube, 0.25_r8*posXflow, reshape([&
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.25_r8, 0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 8/10ths filled in -x, fluxed in +x? ', success
    
    success = fluxing_unit_test (cube, -0.25_r8*posXflow, reshape([&
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.05_r8, 0.2_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 8/10ths filled in -x, fluxed in -x? ', success

    ! cube 1/8th filled in y direction
    vof = [0.125_r8, 0.875_r8]
    intnorm(:,1) = posyn; intnorm(:,2) = -intnorm(:,1)
    call cube%partition (vof, intnorm)

    success = fluxing_unit_test (cube, 0.5_r8*posXflow, reshape([&
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.0625_r8,  0.4375_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 1/8th filled in y, fluxed in +x?    ', success

    ! cube 1/8th filled along xy diagonal
    vof = [0.125_r8, 0.875_r8]
    intnorm(:,1) = posxyn; intnorm(:,2) = -intnorm(:,1)
    call cube%partition (vof, intnorm)

    success = fluxing_unit_test (cube, 0.5_r8*posXflow, reshape([&
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.5_r8, &
         0.0_r8,  0.0_r8, &
         0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 1/8th filled in xy, fluxed in +x?   ', success

    success = fluxing_unit_test (cube, -0.5_r8*posXflow, reshape([&
         0.0_r8,   0.0_r8,   &
         0.0_r8,   0.0_r8,   &
         0.125_r8, 0.375_r8, &
         0.0_r8,   0.0_r8,   &
         0.0_r8,   0.0_r8,   &
         0.0_r8,   0.0_r8], [2,6]) )
    write(*,*) 'passed cube 1/8th filled in xy, fluxed in +x?   ', success

    ! cube 1/8th filled along xy diagonal
    vof = [0.125_r8, 0.875_r8]
    intnorm(:,1) = posxyzn; intnorm(:,2) = -intnorm(:,1)
    call cube%partition (vof, intnorm)

    tmp = (0.25_r8 - (1.0_r8 - (6.0_r8*(0.125_r8))**(1.0_r8/3.0_r8)))**3.0_r8/6.0_r8
    success = fluxing_unit_test (cube, 0.25_r8*posXflow, reshape([&
         0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, &
         tmp,    0.25_r8 - tmp, &
         0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8], [2,6]) )
    write(*,*) 'passed cube 1/8th filled in xyz, fluxed in +x?  ', success

    write(*,*) '===================================================='
    write(*,*)

  end subroutine multimat_cell_unit_test_suite

  logical function fluxing_unit_test (cell, fluxing_velocity, volflux_ex)

    use consts,      only: cutvof,nfc
    use array_utils, only: isZero

    type(multimat_cell), intent(inout) :: cell
    real(r8),            intent(in)    :: fluxing_velocity(:),volflux_ex(:,:)

    real(r8) :: outward_volflux(size(volflux_ex,dim=1),nfc)
    integer  :: nmat

    nmat = size(volflux_ex,dim=1)

    outward_volflux = cell%outward_volflux (1.0_r8, fluxing_velocity)
    fluxing_unit_test = all(isZero (outward_volflux-volflux_ex,cutvof))
    
  end function fluxing_unit_test

  ! given a set of VoFs, normals, and an order,
  ! create child polyhedra for each material
  subroutine partition (this, vof, norm)

    use consts, only: cutvof
    use locate_plane_nd_module
    use logging_services

    class(multimat_cell), intent(inout) :: this
    real(r8),             intent(in)    :: vof(:), norm(:,:)

    type(polyhedron) :: tmp(2),remainder
    integer          :: m,nm

    if (allocated(this%mat_poly)) deallocate(this%mat_poly)
    allocate(this%mat_poly(size(vof)))
    this%mat_poly(:)%nVerts = 0

    call remainder%init (this)
    
    this%nmat = count(vof > cutvof)
    nm = 0

    do m = 1,size(vof)
      if (vof(m) < cutvof) cycle
      nm = nm+1 ! update the counter of how many materials we've seen thus far
      
      ! reconstruct the plane from the remaining free space
      ! use the plane to generate the polyhedron for this material,
      ! and update the free-space polyhedron

      if (nm==this%nmat .or. (1.0_r8-cutvof)*remainder%volume () < vof(m)*this%volume ()) then
        ! if this is the final material in the cell,
        ! it gets the entire remainder of the polyhedron
        this%mat_poly(m) = remainder
        if (this%nmat==1) this%m = m ! if this is the only material, store its ID
        exit
      else
        ! if this is not the final material in the cell, split the cell
        tmp = remainder%split (locate_plane_nd (remainder, norm(:,m), vof(m)*this%volume ()))
        remainder = tmp(1)
        this%mat_poly(m) = tmp(2)
      end if
    end do
    
  end subroutine partition

  ! given a plane, find the volumes of each
  ! material behind that plane (flux volumes)
  function volumes_behind_plane (this, P) result(vol)
    use plane_type
    
    class(multimat_cell), intent(in) :: this
    class(plane),         intent(in) :: P
    real(r8)                         :: vol(size(this%mat_poly))

    integer                          :: m
    
    do m = 1,size(this%mat_poly)
      vol(m) = this%mat_poly(m)%volume_behind_plane (P)
    end do

  end function volumes_behind_plane

  function outward_volflux (this, adv_dt, fluxing_velocity, face_area)

    use consts, only: nfc,cutvof
    use plane_type
    use logging_services
    
    class(multimat_cell), intent(inout) :: this !inout because of call to volume
    real(r8),             intent(in)    :: adv_dt, fluxing_velocity(:)
    real(r8), optional,   intent(in)    :: face_area(:)
    real(r8)                            :: outward_volflux(size(this%mat_poly),nfc)

    type(plane) :: flux_plane
    real(r8)    :: xf(3)
    integer     :: f,nV,m
    
    do f = 1,nfc
      if (fluxing_velocity(f) < cutvof*this%volume ()) then
        outward_volflux(:,f) = 0.0_r8
      else
        ! find the plane equation for the back end of the flux volume
        ! WARNING: in general, this could be non-planar, just like cell faces
        flux_plane%normal = -this%face_normal(:,f)

        nV = count(this%face_vid(:,f) /= 0) ! number of vertices on this face
        xf = sum(this%x(:,this%face_vid(1:nV,f)),dim=2) / real(nV,r8) ! face center

        flux_plane%rho  = sum(xf*flux_plane%normal) + adv_dt * fluxing_velocity(f)
        
        ! find the volume of the materials behind the flux plane
        if (this%nmat==1 .and. present(face_area)) then
          ! pure hex cells are easy, don't cut up polyhedrons to calculate the flux volume
          outward_volflux(:,f) = 0.0_r8
          outward_volflux(this%m,f) = adv_dt * fluxing_velocity(f) * face_area(f)
        else
          outward_volflux(:,f) = this%volumes_behind_plane (flux_plane)
        end if

        ! make sure we have a valid outward volume flux
        if (any(outward_volflux(:,f) < 0.0_r8)) then
          write(*,*)
          write(*,'(a,i6,4es14.4)') 'f,volflux: ',f,outward_volflux(:,f)
          write(*,'(a,es10.4)') 'correct tot volflux: ', adv_dt * fluxing_velocity(f) * face_area(f)

          write(*,'(a,4es20.10)') 'flux plane n,p: ',flux_plane%normal, flux_plane%rho
          write(*,*)

          call this%print_data ()

          call LS_fatal ("in nested dissection outward_volflux: negative fluxes calculated!")
        end if
      end if
    end do


  end function outward_volflux

  ! returns a polygon for a desired interface
  type(polygon) function interface_polygon (this,m)

    use polygon_type

    class(multimat_cell), intent(in) :: this
    integer,              intent(in) :: m

    integer :: interface_face_id,nVerts

    interface_polygon%nVerts = 0
    ! if this polyhedron doesn't exist, or the polyhedron describes a pure cell,
    ! there is no interface to find
    if (m > size(this%mat_poly) .or. this%nmat<2 .or. this%mat_poly(m)%nVerts<4) return

    ! by the convention set in polyhedron_type%polyhedron_on_side_of_plane,
    ! the face corresponding to the phase interface is the last face in the polyhedron
    interface_face_id = this%mat_poly(m)%nFaces
    
    ! count how many real vertices are listed for this face (0s represent non-existent vertices)
    nVerts = count(this%mat_poly(m)%face_vid(:,interface_face_id)/=0)

    ! initialize the polyhedron with the vertices used by the interface face
    ! and the corresponding normal vector
    call interface_polygon%init (this%mat_poly(m)%x(:,this%mat_poly(m)%face_vid(1:nVerts,interface_face_id))) !, &
         !this%mat_poly(m)%face_normal(:,interface_face_id))

  end function interface_polygon

end module multimat_cell_type
