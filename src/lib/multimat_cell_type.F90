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
  use logging_services
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
    procedure, private :: volumes_behind_plane
    procedure :: outward_volflux
    procedure :: interface_polygon
  end type multimat_cell
  
contains

  subroutine multimat_cell_unit_test_suite ()
    use hex_types, only: cube_v, hex_f, hex_e

    type(multimat_cell) :: cube
    logical             :: success
    integer             :: ierr
    real(r8)            :: posXflow(6), posZflow(6), posxyzn(3), posxyn(3), posxn(3), &
        posyn(3), tmp
    real(r8), allocatable :: vof(:), intnorm(:,:)

    
    write(*,*)
    write(*,*) 'MULTIMAT CELL TYPE'
    write(*,*) '===================================================='

    ! partitioning
    write(*,*) 'CELL PARTITIONING'
    ! cube with 3-matls, one very small
    call cube%init (ierr, reshape([&
        2.8125E-01_r8,    5.6250E-01_r8,    2.8125E-01_r8,&
        3.1250E-01_r8,    5.6250E-01_r8,    2.8125E-01_r8,&
        3.1250E-01_r8,    5.9375E-01_r8,    2.8125E-01_r8,&
        2.8125E-01_r8,    5.9375E-01_r8,    2.8125E-01_r8,&
        2.8125E-01_r8,    5.6250E-01_r8,    3.1250E-01_r8,&
        3.1250E-01_r8,    5.6250E-01_r8,    3.1250E-01_r8,&
        3.1250E-01_r8,    5.9375E-01_r8,    3.1250E-01_r8,&
        2.8125E-01_r8,    5.9375E-01_r8,    3.1250E-01_r8],&
        [3,8]), hex_f, hex_e)
    vof = [8.7485114822E-01_r8,    1.2514879878E-01_r8,    5.2994441145E-08_r8]
    intnorm = reshape([&
        9.8094153325E-01_r8,   -1.1866987389E-01_r8,   -1.5385437718E-01_r8,&
        -9.8094153821E-01_r8,    1.1866966066E-01_r8,    1.5385451003E-01_r8,&
        -4.6914971675E-01_r8,    7.9578595328E-01_r8,   -3.8291416771E-01_r8], [3,3])

    !call cube%partition (vof, intnorm)
    success = partition_unit_test (cube, vof, intnorm)
    ! success = fluxing_unit_test (cube, 2.692080062668622E-003*2.2845404866E-01_r8*posZflow, &
    !     reshape([&
    !      0.0_r8, 0.0_r8, 0.0_r8, &
    !      0.0_r8, 0.0_r8, 0.0_r8, &
    !      0.0_r8, 0.0_r8, 0.0_r8, &
    !      0.0_r8, 0.0_r8, 0.0_r8, &
    !      0.0_r8, 0.0_r8, 0.0_r8, &
    !      0.0_r8, 0.0_r8, 0.0_r8], [3,6]) )

    write(*,*) '3-matl cube? ', success

    ! 3-matl cube2
    call cube%init (ierr, reshape([&
        2.96875E-01_r8,    3.75000E-01_r8,    5.31250E-01_r8,&
        3.12500E-01_r8,    3.75000E-01_r8,    5.31250E-01_r8,&
        3.12500E-01_r8,    3.90625E-01_r8,    5.31250E-01_r8,&
        2.96875E-01_r8,    3.90625E-01_r8,    5.31250E-01_r8,&
        2.96875E-01_r8,    3.75000E-01_r8,    5.46875E-01_r8,&
        3.12500E-01_r8,    3.75000E-01_r8,    5.46875E-01_r8,&
        3.12500E-01_r8,    3.90625E-01_r8,    5.46875E-01_r8,&
        2.96875E-01_r8,    3.90625E-01_r8,    5.46875E-01_r8], [3,8]), hex_f, hex_e)
    vof = [0.811404648085707_r8,  0.188595062059073_r8,   2.898552200505922E-007_r8]
    intnorm = reshape([&
        3.5592178718E-01_r8,   -9.3450685990E-01_r8,    4.0755626314E-03_r8,&
        -3.5592493730E-01_r8,    9.3450564851E-01_r8,   -4.0782237116E-03_r8,&
        5.7303886156E-01_r8,    5.7634688730E-01_r8,    5.8262400280E-01_r8], [3,3])
    success = partition_unit_test (cube, vof, intnorm)
    write(*,*) '3-matl cube2? ', success
    
    
    ! fluxing
    write(*,*) 'FLUXING'
    call cube%init (ierr, cube_v, hex_f, hex_e)

    ! define face velocities [+y, -y, -x, +x, -z, +z]
    posXflow = [ 0.0_r8, 0.0_r8, -1.0_r8, 1.0_r8,  0.0_r8, 0.0_r8 ]
    posZflow = [ 0.0_r8, 0.0_r8,  0.0_r8, 0.0_r8, -1.0_r8, 1.0_r8 ]

    ! define normals
    posxyzn = [1.0_r8, 1.0_r8, 1.0_r8] / sqrt(3.0_r8) ! positive x-y-z
    posxn   = [1.0_r8, 0.0_r8, 0.0_r8]                ! positive x-direction
    posyn   = [0.0_r8, 1.0_r8, 0.0_r8]                ! positive y-direction
    posxyn  = (posxn + posyn) / sqrt(2.0_r8) ! positive x-y plane
    
    ! cube 8/10ths filled in x direction
    deallocate(intnorm)
    allocate(intnorm(3,2))
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

    ! cube 1/8th filled along xyz diagonal
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
    integer  :: nmat, ierr
    
    nmat = size(volflux_ex,dim=1)

    write(*,*) 'WARNING: hardcoding face areas'
    outward_volflux = cell%outward_volflux (1.0_r8, fluxing_velocity, &
        [1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8],ierr)
    if (ierr /= 0) then
      call LS_fatal ("outward volflux failed")
    end if
    fluxing_unit_test = all(isZero (outward_volflux-volflux_ex,cutvof))

    ! if (.not.fluxing_unit_test) then
    !   write(*,*) 'failed: '
    !   write(*,'(6es20.10)') outward_volflux
    !   write(*,'(6es20.10)') volflux_ex
    ! end if
    
  end function fluxing_unit_test

  logical function partition_unit_test (cell, vof, intnorm)

    use consts, only: cutvof
    use array_utils, only: isZero

    type(multimat_cell), intent(inout) :: cell
    real(r8),            intent(in)    :: vof(:), intnorm(:,:)

    real(r8) :: vof_result
    integer :: m

    call cell%partition (vof, intnorm)

    partition_unit_test = .true.
    do m = 1,size(vof)
      vof_result = cell%mat_poly(m)%volume() / cell%volume()
      write(*,*) 'mv', m, vof_result, vof(m)
      if (.not.isZero(vof(m) - vof_result, cutvof)) then
        write(*,*) 'partition test failed on material: ',m
        write(*,'(a,4es20.10)') 'vof, vof_result, error: ',vof(m), vof_result, abs(vof(m)-vof_result)
        write(*,'(a,4es20.10)') 'remaining vof: ',1.0_r8 - sum(vof(:m))
        write(*,*) 'failure may simply be due to insufficient number of iterations'
        partition_unit_test = .false.
      end if
    end do

  end function partition_unit_test

  ! given a set of VoFs, normals, and an order,
  ! create child polyhedra for each material
  subroutine partition (this, vof, norm)

    use consts, only: cutvof
    use array_utils, only: isZero
    use plane_type
    use locate_plane_nd_module

    class(multimat_cell), intent(inout) :: this
    real(r8),             intent(in)    :: vof(:), norm(:,:)

    type(plane)      :: interface_plane
    type(polyhedron) :: tmp(2),remainder
    integer          :: m,nm,ierr

    ierr = 0
    if (allocated(this%mat_poly)) deallocate(this%mat_poly)
    allocate(this%mat_poly(size(vof)))
    this%mat_poly(:)%nVerts = 0
    this%m = 0
    
    call remainder%init (this)
    
    this%nmat = count(vof > cutvof)
    nm = 0
    
    do m = 1,size(vof)
      if (vof(m) < cutvof) cycle
      nm = nm+1 ! update the counter of how many materials we've seen thus far
      
      ! reconstruct the plane from the remaining free space
      ! use the plane to generate the polyhedron for this material,
      ! and update the free-space polyhedron
      if (nm==this%nmat .or. (1.0_r8-cutvof)*remainder%volume() < vof(m)*this%volume() .or. &
          isZero(remainder%volume())) then
        ! if this is the final material in the cell,
        ! it gets the entire remainder of the polyhedron
        this%mat_poly(m) = remainder
        if (this%nmat==1) this%m = m ! if this is the only material, store its ID
        exit
      else
        ! if this is not the final material in the cell, split the cell
        interface_plane = locate_plane_nd (remainder, norm(:,m), vof(m)*this%volume(), this%volume())
        !tmp = remainder%split (interface_plane)
        call remainder%split (interface_plane,tmp,ierr)

        ! this check ensures the partitions give their vofs within the requested cutvof
        ! it will fail if the Brent's iterations did not converge within the maximum
        ! number of allowed iterations
        ! if (.not.isZero(tmp(1)%volume() / this%volume() - (1.0_r8 - sum(vof(:m))), cutvof)) ierr=1

        if (ierr/=0) call partitionError ()

        remainder = tmp(1)
        this%mat_poly(m) = tmp(2)
      end if
    end do
    
  contains

    subroutine partitionError ()

      integer :: i
      real(r8) :: rem

      write(*,*)
      write(*,*) 'partition error!'
      
      write(*,*) 'cell:'
      call this%print_data ()

      write(*,*) 'vof: ',vof

      write(*,*) 'norms:'
      do i = 1,size(norm, dim=2)
        write(*,'(3es20.10)') norm(:,i)
      end do

      write(*,*) 'previous vofs:'
      rem = 1.0_r8
      do i = 1,m-1
        write(*,'(i3,a,2es20.10)') i,': ',this%mat_poly(i)%volume() / this%volume()
        rem = rem - this%mat_poly(i)%volume() / this%volume()
      end do
      write(*,*) 'remaining vof: ',rem
      rem = rem - tmp(2)%volume() / this%volume()
      write(*,*) 'remaining after last cutout vof: ',rem
      
      write(*,*)
      write(*,*) 'previous remainder polyhedron data: '
      call remainder%print_data ()
      
      write(*,*)
      write(*,*) 'material ',m,' cutout: '
      call tmp(2)%print_data ()

      write(*,*)
      write(*,*) 'remainder polyhedron data: '
      call tmp(1)%print_data ()

      write(*,*) 'interface plane:'
      call interface_plane%print_data ()

      write(*,*)
      ! write(*,*) 'remainder polyhedron volume does not match remaining volume from vof'
      ! write(*,'(a,es20.10)') 'remainder polyhedron vof: ',tmp(1)%volume() / this%volume()
      ! write(*,'(a,es20.10)') 'exact remaining vof:      ',1.0_r8 - sum(vof(:m))
      ! write(*,'(a,es20.10)') 'error:                    ',&
      !     abs(tmp(1)%volume() / this%volume() - (1.0_r8 - sum(vof(:m))))

      call LS_fatal ("partition error")
    end subroutine partitionError

  end subroutine partition

  ! given a plane, find the volumes of each
  ! material behind that plane (flux volumes)
  function volumes_behind_plane (this, P, ierr) result(vol)

    use plane_type
    
    class(multimat_cell), intent(in) :: this
    class(plane),         intent(in) :: P
    integer,              intent(out) :: ierr
    real(r8)                         :: vol(size(this%mat_poly))

    integer                          :: m

    ierr = 0
    do m = 1,size(this%mat_poly)
      vol(m) = this%mat_poly(m)%volume_behind_plane (P,ierr)
      if (ierr /= 0) then
        write(*,*) 'volumes_behind_plane failed on material ',m
        return
      end if
    end do

  end function volumes_behind_plane

  function outward_volflux (this, adv_dt, fluxing_velocity, face_area, ierr)

    use consts, only: nfc,cutvof
    use plane_type
    
    class(multimat_cell), intent(inout) :: this !inout because of call to volume
    real(r8),             intent(in)    :: adv_dt, fluxing_velocity(:), face_area(:)
    !real(r8), optional,   intent(in)    :: face_area(:) ! WARNING: this isn't really optional, and maybe should be copied into the object
    integer, intent(out) :: ierr
    real(r8)                            :: outward_volflux(size(this%mat_poly),nfc)

    type(plane) :: flux_plane
    real(r8)    :: xf(3)
    integer     :: f,nV,m

    ierr = 0
    do f = 1,nfc
      if (fluxing_velocity(f)*adv_dt*face_area(f) < cutvof*this%volume()) then
        ! if we would be fluxing very very little, don't flux anything
        outward_volflux(:,f) = 0.0_r8
      else
        ! find the plane equation for the back end of the flux volume
        ! WARNING: in general, this could be non-planar, just like cell faces
        flux_plane%normal = -this%face_normal(:,f)

        nV = count(this%face_vid(:,f) /= 0) ! number of vertices on this face
        xf = sum(this%x(:,this%face_vid(1:nV,f)),dim=2) / real(nV,r8) ! face center

        flux_plane%rho  = sum(xf*flux_plane%normal) + adv_dt * fluxing_velocity(f)
        
        ! find the volume of the materials behind the flux plane
        if (this%nmat==1) then
          ! pure hex cells are easy, don't cut up polyhedrons to calculate the flux volume
          outward_volflux(:,f) = 0.0_r8
          outward_volflux(this%m,f) = adv_dt * fluxing_velocity(f) * face_area(f)
        else
          outward_volflux(:,f) = this%volumes_behind_plane (flux_plane,ierr)
          if (ierr /= 0) then
            write(*,*) 'outward_volflux failed on face ',f
            return
          end if
        end if

        ! make sure we have a valid outward volume flux
        if (any(outward_volflux(:,f) < 0.0_r8)) then
          write(*,*)
          write(*,'(a,i6,4es14.4)') 'f,volflux: ',f,outward_volflux(:,f)
          write(*,'(a,es10.4)') 'correct tot volflux: ', adv_dt * fluxing_velocity(f) * face_area(f)

          write(*,'(a,4es20.10)') 'flux plane n,p: ',flux_plane%normal, flux_plane%rho
          write(*,*) 'nmat ',this%nmat, this%m

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
