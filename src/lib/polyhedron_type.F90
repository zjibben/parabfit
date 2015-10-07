!!
!! polyhedron_type
!!
!! This module defines an arbitrary polyhedron type, along with routines for
!! calculating volume, splitting polyhedra, locating intersections, etc.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!
!! References:
!!     1. Mirtich. Fast and Accurate Computation of Polehedral Mass Properties. Journal of Graphics Tools, 1996.
!! 

module polyhedron_type
  use kinds, only: r8
  use logging_services
  use polygon_type
  implicit none
  private

  public :: polyhedron_unit_test

  type, public :: polyhedron
    !private
    integer               :: nVerts, nEdges, nFaces       ! number of vertices, edges, and faces
    real(r8), allocatable :: x(:,:)                       ! vertex positions
    integer,  allocatable :: face_vid(:,:), edge_vid(:,:) ! face and edge IDs
  contains
    procedure :: init => init_polyhedron
    procedure :: volume
    procedure :: intersection_verts
  end type polyhedron

  real(r8), parameter :: alittle = epsilon(1.0_r8)

contains

  subroutine polyhedron_unit_test ()
    use plane_type
    use hex_types, only: hex_f, hex_e

    type(polyhedron) :: cube,pyramid
    real(r8)         :: cube_v(3,8) = [ &
         [ 0.0_r8, 0.0_r8, 0.0_r8 ], & ! vertex positions
         [ 1.0_r8, 0.0_r8, 0.0_r8 ], &
         [ 1.0_r8, 1.0_r8, 0.0_r8 ], &
         [ 0.0_r8, 1.0_r8, 0.0_r8 ], &
         [ 0.0_r8, 0.0_r8, 1.0_r8 ], &
         [ 1.0_r8, 0.0_r8, 1.0_r8 ], &
         [ 1.0_r8, 1.0_r8, 1.0_r8 ], &
         [ 0.0_r8, 1.0_r8, 1.0_r8 ]  &
         ]
    real(r8)         :: pyr_v(3,5) = [ &
         [ 0.0_r8, 0.0_r8, 0.0_r8 ], & ! vertex positions
         [ 1.0_r8, 0.0_r8, 0.0_r8 ], &
         [ 1.0_r8, 1.0_r8, 0.0_r8 ], &
         [ 0.0_r8, 1.0_r8, 0.0_r8 ], &
         [ 0.5_r8, 0.5_r8, 0.5_r8 ]  &
         ]
    integer          :: pyr_f(4,5) = [ &
         [ 1,2,3,4 ], & ! face vertices
         [ 1,2,5,0 ], &
         [ 2,3,5,0 ], &
         [ 3,4,5,0 ], &
         [ 4,1,5,0 ]  &
         ]
    integer          :: pyr_e(2,8) = [ &
         [ 1,2 ], & ! edge vertices
         [ 2,3 ], &
         [ 3,4 ], &
         [ 4,1 ], &
         [ 1,5 ], &
         [ 2,5 ], &
         [ 3,5 ], &
         [ 4,5 ]  &
         ]
    real(r8)         :: volume
    type(polygon)    :: intpoly
    type(plane)      :: P
    integer          :: i

    write(*,*)

    ! calculate the volume of a unit cube (1)
    call cube%init (cube_v, hex_f, hex_e)
    volume = cube%volume ()
    write(*,*) 'cube volume = ', volume

    ! calculate the volume of a pyramid (1/6)
    call pyramid%init (pyr_v, pyr_f, pyr_e)
    volume = pyramid%volume ()
    write(*,*) 'pyramid volume = ', volume

    ! create a plane, and return coordinates it intersects with polyhedron edges
    P%normal = [ 1.0_r8, 0.0_r8, 0.0_r8 ]
    P%rho    = 0.5_r8

    intpoly = cube%intersection_verts (P)
    
    write(*,*) 'intersection points'
    do i = 1,intpoly%nVerts
      write(*,*) intpoly%x(:,i)
    end do
    
  end subroutine polyhedron_unit_test

  subroutine init_polyhedron (this, x, face_v, edge_v)
    class(polyhedron), intent(out) :: this
    real(r8),          intent(in)  :: x(:,:)
    integer,           intent(in)  :: face_v(:,:), edge_v(:,:)
    
    integer :: f,nV

    this%nVerts = size(x,     dim=2)
    this%nEdges = size(edge_v,dim=2)
    this%nFaces = size(face_v,dim=2)
    allocate( this%x(3,this%nVerts), & !this%face(this%nFaces), &
         this%face_vid(size(face_v,dim=1),this%nVerts),&
         this%edge_vid(size(edge_v,dim=1),this%nVerts) )

    this%x = x
    this%edge_vid = edge_v
    this%face_vid = face_v
    
    ! do f = 1,this%nFaces
    !   nV = count(face_v(:,f) /= 0) ! number of vertices on this face
    !   call this%face(f)%init (x(:,face_v(1:nV,f)))
    ! end do
    
    ! note that by taking the cross product between edges described in a
    ! counter-clockwise manner, we guarantee the normal to be outward facing
    
  end subroutine init_polyhedron

  !
  ! This function calculates the volume of a polehedron following the
  ! algorithm presented by [1]. It calculates the sum of the surface
  ! integral of x over all faces, which is equal to the volume by the
  ! divergence theorem.
  !
  real(r8) function volume (this)
    class(polyhedron), intent(in) :: this
    
    integer :: f,nV
    type(polygon) :: face

    ! sum up the integral of n_x*x over all faces
    ! (could just as easily be any other direction)
    volume = 0.0_r8
    do f = 1,this%nFaces
      ! generate a polygon from this face
      ! WARNING: not sure if face%init is leaking memory with the x allocation here
      nV = count(this%face_vid(:,f) /= 0) ! number of vertices on this face
      call face%init (this%x(:,this%face_vid(1:nV,f)))

      ! calculate this face's contribution
      volume = volume + face%norm(1) * face%intXdA (1)
    end do

  end function volume

  !
  ! Given an equation of a plane and a polyhedron, return a polygon from the
  ! points where the plane intersects the polyhedron edges
  !
  type(polygon) function intersection_verts (this,P,vof)
    use plane_type

    class(polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P
    real(r8),          intent(in), optional :: vof

    integer       :: e,Nintersections
    real(r8)      :: x(3,this%nEdges),intx(3)

    ! loop through all edges
    Nintersections = 0
    do e = 1,this%nEdges
      ! check if the P intersects this edge
      if (P%intersects(this%x(:,this%edge_vid(:,e)))) then
        ! if it does, find the point where they intersect
        intx = P%intersection_point (this%x(:,this%edge_vid(:,e)))
        ! if this point wasn't already found, store it
        if (.not.any( &
             abs(intx(1)-x(1,1:Nintersections))<1e4_r8*alittle .and. &
             abs(intx(2)-x(2,1:Nintersections))<1e4_r8*alittle .and. &
             abs(intx(3)-x(3,1:Nintersections))<1e4_r8*alittle ) ) then
          
          Nintersections = Nintersections + 1
          x(:,Nintersections) = intx
        end if
      end if
    end do
    
    ! pass the intersection points to the polygon constructor
    if (Nintersections>2) then
      call intersection_verts%init (x(:,1:Nintersections))
      call intersection_verts%order () ! this probably doesn't need to be called every time this function is used
    else
      ! do e = 1,this%nVerts
      !   write(*,*) this%x(:,e)
      ! end do
      ! write(*,'(a,3f)') 'plane n: ',P%normal
      ! write(*,'(a, f)') 'plane m: ',P%rho
      
      ! do e = 1,Nintersections
      !   write(*,'(a,3f)') 'int x: ', x(:,e)
      ! end do
      ! write(*,*) 'vof',vof

      ! call LS_fatal ('no intersection -- need to handle this case')
    end if
    
  end function intersection_verts

  ! these two utility functions would be best pulled aside to a separate module
  pure logical function isZero (x)
    real(r8), intent(in) :: x
    isZero = abs(x)<1e4_r8*alittle
  end function isZero

  pure logical function eq (a,b)
    real(r8), intent(in) :: a,b
    eq = isZero (a-b)
  end function eq
  
end module polyhedron_type
