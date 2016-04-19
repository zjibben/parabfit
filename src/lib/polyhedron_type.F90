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
!!     1. Hopcroft and Kahn. A Paradigm for Robust Geometric Algorithms. Algorithmica, 1992
!! 

module polyhedron_type
  use kinds,  only: r8
  use logging_services
  use polygon_type
  implicit none
  private

  public :: polyhedron_unit_test

  type, public :: polyhedron
    !private
    integer               :: nVerts, nEdges, nFaces       ! number of vertices, edges, and faces
    real(r8), allocatable :: x(:,:),face_normal(:,:)      ! vertex positions and face outward normals
    integer,  allocatable :: face_vid(:,:), edge_vid(:,:) ! face and edge IDs
    real(r8)              :: vol                          ! volume
    ! ideally vol should be private, but Fortran then also hides it from child types
  contains
    procedure, private :: init_polyhedron
    procedure, private :: init_polyhedron_copy
    generic            :: init => init_polyhedron, init_polyhedron_copy
    procedure          :: volume
    procedure          :: intersection_verts
    procedure          :: split
    procedure          :: volume_behind_plane
    procedure          :: print_data
    !procedure, private :: edge_containing_vertices
    procedure, private :: polyhedron_on_side_of_plane
  end type polyhedron

contains

  subroutine polyhedron_unit_test ()
    use plane_type
    use hex_types,   only: hex_f, hex_e, cube_v
    use array_utils, only: isZero

    type(polyhedron) :: cube,pyramid,pyramid3,cutcube,tmp(2)
    type(polygon)    :: intpoly
    type(plane)      :: P
    real(r8)         :: volume,tmpr1,tmpr2
    logical          :: success
    real(r8)         :: pyr3_v(3,4) = reshape([ &
         0.0_r8, 0.0_r8, 0.0_r8, & ! vertex positions
         0.0_r8, 1.0_r8, 0.0_r8, &
         1.0_r8, 0.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 1.0_r8],&
         shape(pyr3_v))
    integer          :: pyr3_f(3,4) = reshape([ &
         1,2,3, & ! face vertices
         1,4,2, &
         1,3,4, &
         4,3,2],&
         shape(pyr3_f))
    integer          :: pyr3_e(2,6) = reshape([ &
         1,2, & ! edge vertices
         2,3, &
         3,4, &
         1,3, &
         4,1, &
         4,2],&
         shape(pyr3_e))
    real(r8)         :: pyr_v(3,5) = reshape([ & ! pyramid with square bottom
         0.0_r8, 0.0_r8, 0.0_r8, & ! vertex positions
         1.0_r8, 0.0_r8, 0.0_r8, &
         1.0_r8, 1.0_r8, 0.0_r8, &
         0.0_r8, 1.0_r8, 0.0_r8, &
         0.5_r8, 0.5_r8, 0.5_r8],&
         shape(pyr_v))
    integer          :: pyr_f(4,5) = reshape([ &
         1,2,3,4, & ! face vertices
         1,2,5,0, &
         2,3,5,0, &
         3,4,5,0, &
         4,1,5,0],&
         shape(pyr_f))
    integer          :: pyr_e(2,8) = reshape([ &
         1,2, & ! edge vertices
         2,3, &
         3,4, &
         4,1, &
         1,5, &
         2,5, &
         3,5, &
         4,5],&
         shape(pyr_e))
    real(r8)         :: cutcube_v(3,10) = reshape([ &
         0.9_r8, 0.0_r8, 0.0_r8, & ! cube cut with triangle face
         0.0_r8, 0.9_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 0.9_r8, &
         1.0_r8, 0.0_r8, 0.0_r8, &
         1.0_r8, 1.0_r8, 0.0_r8, &
         0.0_r8, 1.0_r8, 0.0_r8, &
         0.0_r8, 0.0_r8, 1.0_r8, &
         1.0_r8, 0.0_r8, 1.0_r8, &
         1.0_r8, 1.0_r8, 1.0_r8, &
         0.0_r8, 1.0_r8, 1.0_r8],&
         shape(cutcube_v))
    integer :: cutcube_f(5,7) = reshape([ &
          3,2,1,0,0, &
          7,3,1,4,8, &
          4,5,9,8,0, &
         10,9,5,6,0, &
         10,6,2,3,7, &
          1,2,6,5,4, &
         10,7,8,9,0],&
         shape(cutcube_f))
    integer :: cutcube_e(2,15) = reshape([ &
          1,2, &
          2,3, &
          3,1, &
          1,4, &
          4,5, &
          5,6, &
          6,2, &
          3,7, &
          7,8, &
          8,9, &
         10,9, &
         10,7, &
         10,6, &
          9,5, &
          8,4],  &
         shape(cutcube_e))

    !     real(r8)         :: pyr3_v(3,4) = [ &
    !      [ 0.0_r8, 0.0_r8, 0.0_r8 ], & ! vertex positions
    !      [ 0.0_r8, 1.0_r8, 0.0_r8 ], &
    !      [ 1.0_r8, 0.0_r8, 0.0_r8 ], &
    !      [ 0.0_r8, 0.0_r8, 1.0_r8 ]  &
    !      ]
    ! integer          :: pyr3_f(3,4) = [ &
    !      [ 1,2,3 ], & ! face vertices
    !      [ 1,4,2 ], &
    !      [ 1,3,4 ], &
    !      [ 4,3,2 ]  &
    !      ]
    ! integer          :: pyr3_e(2,6) = [ &
    !      [ 1,2 ], & ! edge vertices
    !      [ 2,3 ], &
    !      [ 3,4 ], &
    !      [ 1,3 ], &
    !      [ 4,1 ], &
    !      [ 4,2 ]  &
    !      ]
    ! real(r8)         :: pyr_v(3,5) = [ & ! pyramid with square bottom
    !      [ 0.0_r8, 0.0_r8, 0.0_r8 ], & ! vertex positions
    !      [ 1.0_r8, 0.0_r8, 0.0_r8 ], &
    !      [ 1.0_r8, 1.0_r8, 0.0_r8 ], &
    !      [ 0.0_r8, 1.0_r8, 0.0_r8 ], &
    !      [ 0.5_r8, 0.5_r8, 0.5_r8 ]  &
    !      ]
    ! integer          :: pyr_f(4,5) = [ &
    !      [ 1,2,3,4 ], & ! face vertices
    !      [ 1,2,5,0 ], &
    !      [ 2,3,5,0 ], &
    !      [ 3,4,5,0 ], &
    !      [ 4,1,5,0 ]  &
    !      ]
    ! integer          :: pyr_e(2,8) = [ &
    !      [ 1,2 ], & ! edge vertices
    !      [ 2,3 ], &
    !      [ 3,4 ], &
    !      [ 4,1 ], &
    !      [ 1,5 ], &
    !      [ 2,5 ], &
    !      [ 3,5 ], &
    !      [ 4,5 ]  &
    !      ]
    ! real(r8)         :: cutcube_v(3,10) = [ &
    !      [ 0.9_r8, 0.0_r8, 0.0_r8 ], & ! cube cut with triangle face
    !      [ 0.0_r8, 0.9_r8, 0.0_r8 ], &
    !      [ 0.0_r8, 0.0_r8, 0.9_r8 ], &
    !      [ 1.0_r8, 0.0_r8, 0.0_r8 ], &
    !      [ 1.0_r8, 1.0_r8, 0.0_r8 ], &
    !      [ 0.0_r8, 1.0_r8, 0.0_r8 ], &
    !      [ 0.0_r8, 0.0_r8, 1.0_r8 ], &
    !      [ 1.0_r8, 0.0_r8, 1.0_r8 ], &
    !      [ 1.0_r8, 1.0_r8, 1.0_r8 ], &
    !      [ 0.0_r8, 1.0_r8, 1.0_r8 ]  &
    !      ]
    ! integer :: cutcube_f(5,7) = [ &
    !      [  3,2,1,0,0 ], &
    !      [  7,3,1,4,8 ], &
    !      [  4,5,9,8,0 ], &
    !      [ 10,9,5,6,0 ], &
    !      [ 10,6,2,3,7 ], &
    !      [  1,2,6,5,4 ], &
    !      [ 10,7,8,9,0 ]  &
    !      ]
    ! integer :: cutcube_e(2,15) = [ &
    !      [  1,2 ], &
    !      [  2,3 ], &
    !      [  3,1 ], &
    !      [  1,4 ], &
    !      [  4,5 ], &
    !      [  5,6 ], &
    !      [  6,2 ], &
    !      [  3,7 ], &
    !      [  7,8 ], &
    !      [  8,9 ], &
    !      [ 10,9 ], &
    !      [ 10,7 ], &
    !      [ 10,6 ], &
    !      [  9,5 ], &
    !      [  8,4 ]  &
    !      ]


    write(*,*)
    write(*,*) 'POLYHEDRON'
    write(*,*) '===================================================='


    ! calculate the volume of a unit cube (1)
    write(*,*) 'SHAPE VOLUMES'
    call cube%init (cube_v, hex_f, hex_e)
    volume = cube%volume ()
    write(*,*) 'cube volume?                     ', isZero (volume-1.0_r8)

    ! calculate the volume of a pyramid (1/6)
    call pyramid%init (pyr_v, pyr_f, pyr_e)
    volume = pyramid%volume ()
    write(*,*) 'pyramid volume?                  ', isZero (volume-1.0_r8/6.0_r8)
    
    ! calculate the volume of a pyramid (1/6)
    call pyramid3%init (pyr3_v, pyr3_f, pyr3_e)
    volume = pyramid3%volume ()
    write(*,*) 'pyramid3 volume?                 ', isZero (volume-1.0_r8/6.0_r8)

    ! calculate the volume of a "cutcube" (1-0.9**3/6)
    call cutcube%init (cutcube_v, cutcube_f, cutcube_e)
    volume = cutcube%volume ()
    write(*,*) 'cutcube volume?                  ', isZero (volume-(1.0_r8-0.9_r8**3/6.0_r8))

    ! create a plane, and return coordinates it intersects with polyhedron edges
    write(*,*) 'SHAPE SPLITTING'
    P%normal = [ 1.0_r8, 0.0_r8, 0.0_r8 ]
    P%rho    = 0.5_r8

    ! intpoly = cube%intersection_verts (P)
    
    ! write(*,*) 'intersection points'
    ! do i = 1,intpoly%nVerts
    !   write(*,*) intpoly%x(:,i)
    ! end do

    ! split the cube vertically down the center
    !write(*,*) 'cube split volumes'
    tmp = cube%split (P)
    success = isZero (tmp(1)%volume ()-0.5_r8) .and. isZero (tmp(2)%volume ()-0.5_r8)
    write(*,*) 'vertical cut?                    ',success

    ! split the cube at an angle
    P%normal = [ 1.0_r8, 1.0_r8, 0.0_r8 ] / sqrt(2.0_r8)
    P%rho    = 1.5_r8 / sqrt(2.0_r8)
    tmp = cube%split (P)
    success = isZero (tmp(1)%volume ()-0.125_r8) .and. isZero (tmp(2)%volume ()-0.875_r8)
    write(*,*) 'xy-angle off-center cut?         ',success

    ! split the cube at an angle through the center
    P%normal = [ 1.0_r8, 1.0_r8, 0.0_r8 ] / sqrt(2.0_r8)
    P%rho    = 1.0_r8 / sqrt(2.0_r8)
    tmp = cube%split (P)
    success = isZero (tmp(1)%volume ()-0.5_r8) .and. isZero (tmp(2)%volume ()-0.5_r8)
    write(*,*) 'center xy-angle cut?             ',success

    ! split the cube at an angle through the center
    P%normal = [ 1.0_r8, 1.0_r8, 1.0_r8 ] / sqrt(3.0_r8)
    P%rho    = 1.5_r8 / sqrt(3.0_r8)
    tmp = cube%split (P)
    success = isZero (tmp(1)%volume ()-0.5_r8) .and. isZero (tmp(2)%volume ()-0.5_r8)
    write(*,*) 'center xyz-angle cut?            ',success

    ! split the cube at an angle through an offset
    P%normal = [ 1.0_r8, 1.0_r8, 1.0_r8 ] / sqrt(3.0_r8)
    P%rho    = 0.5_r8/sqrt(3.0_r8)
    tmp = cube%split (P)
    success = isZero (tmp(1)%volume ()-47.0_r8/48.0_r8) .and. isZero (tmp(2)%volume ()-1.0_r8/48.0_r8)
    write(*,*) 'off-center xyz-angle cut?        ',success

    ! split the pyramid in the x direction
    P%normal = -[ 1.0_r8, 0.0_r8, 0.0_r8 ]
    P%rho    = -0.8_r8
    tmp = pyramid3%split (P)
    success = isZero (tmp(1)%volume ()-0.992_r8/6.0_r8) .and. isZero (tmp(2)%volume ()-4e-3_r8/3.0_r8)
    write(*,*) 'pyramid3 cut?                    ',success

    ! split the cutcube
    P%normal = -[ 1.0_r8, 0.0_r8, 0.0_r8 ]
    P%rho    = -0.8_r8
    tmp = cutcube%split (P)
    tmpr1 = 1.0_r8-0.9_r8**3/6.0_r8
    tmpr2 = 0.2_r8 - 0.1_r8**3/6.0_r8
    success = isZero (tmp(1)%volume ()-(tmpr1-tmpr2)) .and. isZero (tmp(2)%volume ()-tmpr2)
    write(*,*) 'cutcube cut?                     ',success

    write(*,*) '===================================================='
    write(*,*)
    
  end subroutine polyhedron_unit_test

  subroutine init_polyhedron (this, x, face_v, edge_v, vol, face_normal)

    class(polyhedron),  intent(out) :: this
    real(r8),           intent(in)  :: x(:,:)
    integer,            intent(in)  :: face_v(:,:), edge_v(:,:)
    real(r8), optional, intent(in)  :: vol, face_normal(:,:)
    
    integer :: f,nV
    
    this%nVerts = size(x,     dim=2)
    this%nEdges = size(edge_v,dim=2)
    this%nFaces = size(face_v,dim=2)

    if (allocated(this%x))           deallocate(this%x)
    if (allocated(this%face_vid))    deallocate(this%face_vid)
    if (allocated(this%edge_vid))    deallocate(this%edge_vid)
    if (allocated(this%face_normal)) deallocate(this%face_normal)
    
    allocate(this%x(3,this%nVerts),& !this%face(this%nFaces),&
         this%face_vid(size(face_v,dim=1),this%nFaces),&
         this%edge_vid(size(edge_v,dim=1),this%nEdges),&
         this%face_normal(3,this%nFaces))
    
    this%x = x
    this%edge_vid = edge_v
    this%face_vid = face_v
    
    this%vol = merge(vol, 0.0_r8, present(vol))
    
    if (present(face_normal)) then
      this%face_normal = face_normal
    else
      this%face_normal = 0.0_r8
    end if
    
    ! if the faces are of type polygon
    ! do f = 1,this%nFaces
    !   nV = count(face_v(:,f) /= 0) ! number of vertices on this face
    !   call this%face(f)%init (x(:,face_v(1:nV,f)))
    ! end do
    
    ! note that by taking the cross product between edges described in a
    ! counter-clockwise manner, we guarantee the normal to be outward facing
    
  end subroutine init_polyhedron

  subroutine init_polyhedron_copy (this,poly)
    class(polyhedron), intent(out) :: this
    class(polyhedron), intent(in)  :: poly

    this%nVerts = poly%nVerts
    this%nEdges = poly%nEdges
    this%nFaces = poly%nFaces
    this%vol    = poly%vol

    if (allocated(this%x))           deallocate(this%x)
    if (allocated(this%face_vid))    deallocate(this%face_vid)
    if (allocated(this%edge_vid))    deallocate(this%edge_vid)
    if (allocated(this%face_normal)) deallocate(this%face_normal)

    allocate( this%x(3,this%nVerts), & !this%face(this%nFaces), &
         this%face_vid(size(poly%face_vid,dim=1),this%nFaces),&
         this%edge_vid(size(poly%edge_vid,dim=1),this%nEdges),&
         this%face_normal(3,this%nFaces) )

    this%x = poly%x
    this%edge_vid = poly%edge_vid
    this%face_vid = poly%face_vid
    this%face_normal = poly%face_normal

  end subroutine init_polyhedron_copy

  !
  ! This function calculates the volume of a polehedron following the
  ! algorithm presented by [1]. It calculates the sum of the surface
  ! integral of x over all faces, which is equal to the volume by the
  ! divergence theorem.
  !
  real(r8) function volume (this)
    
    use array_utils,     only: isZero
    use ieee_arithmetic, only: ieee_is_nan

    class(polyhedron), intent(inout) :: this
    
    integer :: f,nV
    type(polygon) :: face
    
    ! sum up the integral of n_x*x over all faces
    ! (could just as easily be any other direction)
    volume = merge(this%vol, 0.0_r8, .not.ieee_is_nan(this%vol))
    
    if (.not.allocated(this%x) .or. volume /= 0.0_r8) return
    
    do f = 1,this%nFaces
      ! generate a polygon from this face
      nV = count(this%face_vid(:,f) /= 0) ! number of vertices on this face
      call face%init (this%x(:,this%face_vid(1:nV,f)), this%face_normal(:,f))

      ! calculate this face's contribution
      if (.not.isZero (face%norm(1))) volume = volume + face%norm(1) * face%intXdA (1)
    end do
    this%vol = volume

  end function volume

  !
  ! Given an equation of a plane and a polyhedron, return a polygon from the
  ! points where the plane intersects the polyhedron edges
  !
  type(polygon) function intersection_verts (this,P,v_assoc_pe)
    use plane_type
    use array_utils, only: reverse,first_true_loc,isZero

    class(polyhedron), intent(in)  :: this
    class(plane),      intent(in)  :: P
    integer, optional, intent(out) :: v_assoc_pe(:)

    logical       :: pteq(this%nEdges)
    integer       :: e,Nintersections
    real(r8)      :: x(3,this%nEdges),intx(3)

    if (present(v_assoc_pe)) v_assoc_pe = 0

    ! loop through all edges
    Nintersections = 0
    do e = 1,this%nEdges
      ! check if the P intersects this edge
      if (P%intersects(this%x(:,this%edge_vid(:,e)))) then
        ! if it does, find the point where they intersect
        intx = P%intersection_point (this%x(:,this%edge_vid(:,e)))

        ! check if the point is already in the list
        pteq(1:Nintersections) = isZero (intx(1)-x(1,1:Nintersections)) &
             .and.               isZero (intx(2)-x(2,1:Nintersections)) &
             .and.               isZero (intx(3)-x(3,1:Nintersections))
        
        ! if this point wasn't already found, store it
        if (.not.any(pteq(1:Nintersections))) then
          Nintersections = Nintersections + 1
          x(:,Nintersections) = intx
          if (present(v_assoc_pe)) v_assoc_pe(e) = Nintersections
        else if (present(v_assoc_pe)) then
          ! if the point was already found,
          ! note that this edge intersects with that found point
          ! this particularly comes into effect when we intersect a vertex
          v_assoc_pe(e) = first_true_loc (pteq(1:Nintersections))
        end if
      end if
    end do
    
    ! pass the intersection points to the polygon constructor
    if (Nintersections>2) then
      call intersection_verts%init (x(:,1:Nintersections))

      ! this probably doesn't need to be called every time this function is used
      if (present(v_assoc_pe)) then
        call intersection_verts%order (v_assoc_pe)
      else
        call intersection_verts%order ()
      end if

      ! make sure the vertices are ordered counter-clockwise
      if (sum(intersection_verts%norm*P%normal)<0.0_r8) then
        ! reverse both x and v_assoc_pe
        intersection_verts%x = reverse (intersection_verts%x)
        if (present(v_assoc_pe)) then
          do e = 1,size(v_assoc_pe)
            if (v_assoc_pe(e)>0) v_assoc_pe(e) = Nintersections - v_assoc_pe(e)+1
          end do
        end if
        call intersection_verts%update_plane_normal ()
      end if
    end if
    
  end function intersection_verts

  ! return a list of the edge ids for edges intersected by the plane
  function intersected_edges (this,P) result(inte)
    use plane_type

    class(polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P

    integer :: e
    logical :: inte(this%nEdges)
    
    do e = 1,this%nEdges
      inte(e) = P%intersects(this%x(:,this%edge_vid(:,e)))
    end do

  end function intersected_edges

  ! return two polyhedrons produced by dividing a given polyhedron with a plane
  ! the first element returned is in front of the plane
  ! the second element returned is behind the plane
  function split (this,P)
    use consts, only: alpha
    use plane_type
    
    class(polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P
    type(polyhedron)              :: split(2)

    integer       :: v, v_assoc_pe(this%nEdges),side(this%nVerts)
    type(polygon) :: intpoly
    real(r8)      :: dist

    ! check which side of the plane each vertex lies
    ! vertices within distance alpha of the plane are considered to lie on it
    !  dist  >  alpha => side =  1
    !  dist  < -alpha => side = -1
    ! |dist| <  alpha => side =  0
    do v = 1,this%nVerts
      dist = P%signed_distance (this%x(:,v)) ! calculate the signed distance from the plane
      side(v) = merge(int(sign(1.0_r8, dist)), 0, abs(dist) > alpha) ! decide where it lies
    end do

    if (all(side>=0)) then
      split(1) = this
    else if (all(side<=0)) then
      split(2) = this
    else
      intpoly = this%intersection_verts (P,v_assoc_pe)

      split(1) = this%polyhedron_on_side_of_plane ( 1, side, intpoly, v_assoc_pe)
      split(2) = this%polyhedron_on_side_of_plane (-1, side, intpoly, v_assoc_pe)
    end if

    ! if any of the polyhedrons have a face with less than 3 vertices, throw an error
    if (any(count(split(1)%face_vid /= 0, dim=1) < 3) .or. &
        any(count(split(2)%face_vid /= 0, dim=1) < 3)) then
      call split(1)%print_data (normalized=.true.)
      call split(2)%print_data (normalized=.true.)
      call LS_fatal ("polyhedron split failed--one of the children has an invalid face")
    end if

  end function split

  ! return the volume behind (opposite normal) a plane and inside the polyhedron
  real(r8) function volume_behind_plane (this,P)

    use consts, only: alpha
    use ieee_arithmetic, only: ieee_is_nan
    use plane_type

    class(polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P

    real(r8)         :: dist
    integer          :: v, side(this%nVerts), v_assoc_pe(this%nEdges)
    type(polyhedron) :: behind, split(2)
    type(polygon)    :: intpoly

    ! check which side of the plane each vertex lies
    ! vertices within distance alpha of the plane are considered to lie on it
    !  dist  >  alpha => side =  1
    !  dist  < -alpha => side = -1
    ! |dist| <  alpha => side =  0
    do v = 1,this%nVerts
      dist = P%signed_distance (this%x(:,v)) ! calculate the signed distance from the plane
      side(v) = merge(int(sign(1.0_r8, dist)), 0, abs(dist) > alpha) ! decide where it lies
    end do
    !write(*,*) 'side',side
    if (all(side<=0)) then
      behind = this
    else if (any(side>0) .and. any(side<0)) then
      intpoly = this%intersection_verts (P,v_assoc_pe)
      behind = this%polyhedron_on_side_of_plane (-1, side, intpoly, v_assoc_pe)
    else
      volume_behind_plane = 0.0_r8
      return
    end if

    ! if any of the polyhedrons have a face with less than 3 vertices, throw an error
    if (any(count(behind%face_vid /= 0, dim=1) < 3)) then
      write(*,*) 'parent polyhedron'
      call this%print_data () !(normalized=.true.)
      write(*,*) side

      write(*,'(a,4es20.10)') 'plane rho, n: ',P%rho, P%normal
      write(*,*)

      write(*,*) 'intersecting polygon'
      call intpoly%print_data ()

      write(*,*) 'new polyhedron'
      call behind%print_data () !(normalized=.true.)

      call LS_fatal ("polyhedron split failed: invalid face")
    end if

    ! calculate the volume of the polyhedron behind the plane
    volume_behind_plane = behind%volume ()

    if (ieee_is_nan(volume_behind_plane)) then
      write(*,*)
      write(*,*) volume_behind_plane
      write(*,'(9i3)') side
      call this%print_data (normalized=.true.)
      if (allocated(intpoly%x)) then
        do v = 1,size(intpoly%x,dim=2)
          write(*,*) intpoly%x(:,v)
        end do
      end if
      call behind%print_data (normalized=.true.)
      call LS_fatal ('polyhedron volume is NaN')
    end if

  end function volume_behind_plane

  ! WARNING: need to update this routine to pass face normals down to the child polyhedron
  ! Reference [1]
  type(polyhedron) function polyhedron_on_side_of_plane (this,valid_side,side,intpoly,v_assoc_pe)
    use array_utils, only: xrange,first_true_loc

    class(polyhedron), intent(in) :: this
    integer,           intent(in) :: side(:)       ! lists which vertices are on the chosen side of the plane
    integer,           intent(in) :: valid_side    ! +/- 1, indicating the side of the plane we want
    type(polygon),     intent(in) :: intpoly       ! polygon of intersection coordinates
    integer,           intent(in) :: v_assoc_pe(:) ! intersection polygon vertex id for a given parent edge id

    integer       :: e,v,f, nV, nVerts, nParVerts, nEdges, nFaces, tmp, invalid_side, new_v, &
         p2c_vid(this%nVerts), &                   ! parent to child vertex id (given parent vertex id, give the child's id)
         edge_vid(2,this%nEdges+intpoly%nVerts), & ! could have intpoly%nVerts more edges than parent polyhedron
         face_vid(size(this%face_vid,dim=1)+2,this%nFaces+1),& ! could have 1 more face and faces could be 2 longer than parent
         edge_cont_verts(this%nEdges,this%nEdges)
    ! note: an updated planar face can only include 1 more node than the original, but I'm not sure if there is a limit
    !       to how many nodes the entirely new face can have. For cubes the number is 2.
    real(r8)      :: x(3,this%nVerts+intpoly%nVerts)

    invalid_side = -valid_side

    ! make a list of vertices for the new polyhedron
    call generate_new_verts (x,p2c_vid,nVerts,nParVerts, this,side,valid_side)
    
    ! find new edges, knowing the original structure and the interface polygon
    call find_edges (edge_vid,nEdges, this,side,valid_side,p2c_vid,nParVerts,nVerts,v_assoc_pe)
    
    ! construct a set of faces from the edge information
    edge_cont_verts = 0
    do e = 1,this%nEdges
      edge_cont_verts(this%edge_vid(1,e),this%edge_vid(2,e)) = e
      edge_cont_verts(this%edge_vid(2,e),this%edge_vid(1,e)) = e
    end do

    ! note these will not be in a particular order, like pececillo would expect for hexes
    call find_faces (face_vid,nFaces, this,side,valid_side,p2c_vid,nParVerts,nVerts,&
        v_assoc_pe,edge_cont_verts)

    ! initialize final polyhedron
    tmp = maxval(count(face_vid(:,:) /= 0,dim=1)) ! the maximum number of vertices on a face
    if (nVerts < 3) then
      call LS_fatal ("not enough vertices for a polygon!")
    end if
    call polyhedron_on_side_of_plane%init (x(:,1:nVerts), face_vid(1:tmp,1:nFaces), &
        edge_vid(:,1:nEdges))

    !write(*,*) 'poly', nVerts, tmp, nFaces, nEdges

  contains
    
    subroutine generate_new_verts (x,p2c_vid,nVerts,nParVerts, this,side,valid_side)

      real(r8),         intent(out) :: x(:,:)
      integer,          intent(out) :: p2c_vid(:), nVerts, nParVerts
      type(polyhedron), intent(in)  :: this
      integer,          intent(in)  :: side(:), valid_side

      integer :: v

      nVerts = 0; p2c_vid = 0
      do v = 1,this%nVerts
        if (side(v)==valid_side) then
          nVerts = nVerts+1
          x(:,nVerts) = this%x(:,v)
          p2c_vid(v) = nVerts ! parent to child vertex id
        end if
      end do
      nParVerts = nVerts                                    ! number of vertices acquired from parent
      x(:,nParVerts+1:nParVerts+intpoly%nVerts) = intpoly%x ! add vertices from plane-polyhedron intersections
      nVerts = nVerts + intpoly%nVerts

    end subroutine generate_new_verts

    subroutine find_edges (edge_vid,nEdges, this,side,valid_side,p2c_vid,nParVerts,nVerts,v_assoc_pe)

      use array_utils, only: containsPair

      integer,          intent(out) :: edge_vid(:,:), nEdges
      type(polyhedron), intent(in)  :: this
      integer,          intent(in)  :: side(:), valid_side, p2c_vid(:), nParVerts, nVerts, &
          v_assoc_pe(:)

      integer :: e, v, n_on_valid_side, new_edge(2)

      nEdges = 0
      do e = 1,this%nEdges
        ! how many of the edge vertices are on the valid side of the plane?
        n_on_valid_side = count(side(this%edge_vid(:,e))==valid_side)
        
        if (n_on_valid_side==2) then      ! this entire edge is on the valid side of the plane
          nEdges = nEdges + 1
          edge_vid(:,nEdges) = p2c_vid(this%edge_vid(:,e))
        else if (n_on_valid_side==1) then ! this edge was divided (or one of the vertices lies on the plane)
          ! the edge consists of the vertex on this side of the plane and the intersection point
          if (side(this%edge_vid(1,e))==valid_side) then
            new_edge = [p2c_vid(this%edge_vid(1,e)),nParVerts+v_assoc_pe(e)]
          else
            new_edge = [p2c_vid(this%edge_vid(2,e)),nParVerts+v_assoc_pe(e)]
          end if

          ! if this edge isn't already listed, add it
          ! the edge might already be listed in cases where we have very close vertices (O(alpha)),
          ! which are then combined in the new polyhedron.
          if (.not.containsPair(new_edge, edge_vid(:,1:nEdges))) then
            nEdges = nEdges + 1
            edge_vid(:,nEdges) = new_edge
          end if
        end if
      end do

      ! points on the plane make up edges with each other
      do v = nParVerts+1,nVerts-1
        nEdges = nEdges + 1
        edge_vid(:,nEdges) = [v,v+1]
      end do
      nEdges = nEdges + 1
      edge_vid(:,nEdges) = [nVerts,nParVerts+1] ! complete the loop

    end subroutine find_edges

    subroutine find_faces (face_vid,nFaces, this,side,valid_side,p2c_vid,nParVerts,nVerts,&
        v_assoc_pe,edge_cont_verts)
      
      integer,          intent(out) :: face_vid(:,:),nFaces
      type(polyhedron), intent(in)  :: this
      integer,          intent(in)  :: side(:), valid_side, p2c_vid(:), nParVerts, nVerts, &
          v_assoc_pe(:),edge_cont_verts(:,:)

      integer :: nV,f,tmp,new_v

      nFaces = 0; face_vid = 0
      do f = 1,this%nFaces
        ! if any vertices for this original face are on the valid side of the
        ! plane, then the face structure can be preserved. Otherwise, it is
        ! thrown out
        nV = count(this%face_vid(:,f) /= 0) ! number of vertices on this face
        if (any(side(this%face_vid(1:nV,f))==valid_side)) then
          nFaces = nFaces + 1
          if (all(side(this%face_vid(1:nV,f))==valid_side)) then
            ! if all vertices are on the valid side of the face,
            ! this face is preserved exactly
            face_vid(1:nV,nFaces) = p2c_vid(this%face_vid(1:nV,f))
          else
            ! if some vertices are on the valid side of the face,
            ! this face is modified

            ! start with the first valid vertex
            v = first_true_loc (side(this%face_vid(1:nV,f))==valid_side)

            ! add all valid vertices in sequence, stop when we hit one that is invalid
            ! this stopping point indicates an edge that was cut by the plane
            tmp = 1
            do while (v<=nV)
              if (side(this%face_vid(v,f))/=valid_side) exit
              face_vid(tmp,nFaces) = p2c_vid(this%face_vid(v,f))
              v = v+1; tmp = tmp+1
            end do

            ! add the new vertices associated with this face
            ! the first is the one that cuts the edge we just stopped on
            ! new_v = nParVerts + &
            !     v_assoc_pe(this%edge_containing_vertices ([this%face_vid(mod(v-2,nV)+1,f),this%face_vid(mod(v-1,nV)+1,f)]))
            new_v = nParVerts + v_assoc_pe(&
                edge_cont_verts(this%face_vid(mod(v-2,nV)+1,f),this%face_vid(mod(v-1,nV)+1,f)))
            if (.not.any(face_vid(1:tmp-1,nFaces)==new_v)) then
              face_vid(tmp  ,nFaces) = new_v
              tmp = tmp+1
            end if

            ! find the next valid vertex
            do while (side(this%face_vid(mod(v-1,nV)+1,f))/=valid_side)
              v = v+1
            end do

            ! the second new vertex is on this edge
            ! new_v = nParVerts + &
            !     v_assoc_pe(this%edge_containing_vertices ([this%face_vid(mod(v-2,nV)+1,f),this%face_vid(mod(v-1,nV)+1,f)]))
            new_v = nParVerts + v_assoc_pe(&
                edge_cont_verts(this%face_vid(mod(v-2,nV)+1,f),this%face_vid(mod(v-1,nV)+1,f)))
            if (.not.any(face_vid(1:tmp-1,nFaces)==new_v)) then
              face_vid(tmp  ,nFaces) = new_v
              tmp = tmp+1
            end if

            ! continue adding vertices that were part of the original face
            do while (v<=nV)
              face_vid(tmp,nFaces) = p2c_vid(this%face_vid(v,f))
              v = v+1; tmp = tmp+1
            end do

          end if
        end if
      end do

      ! add the new face from points lying on the plane
      nFaces = nFaces + 1
      if (valid_side>0) then ! needs to be opposite the input order
        face_vid(1:nVerts-nParVerts,nFaces) = xrange (nVerts,nParVerts+1)
      else
        face_vid(1:nVerts-nParVerts,nFaces) = xrange (nParVerts+1,nVerts)
      end if

    end subroutine find_faces

  end function polyhedron_on_side_of_plane

  ! ! return the edge associated with a pair of vertices
  ! integer function edge_containing_vertices (this,v)
  !   use array_utils, only: reverse

  !   class(polyhedron), intent(in) :: this
  !   integer,           intent(in) :: v(2)
    
  !   do edge_containing_vertices = 1,this%nEdges
  !     if ( all(this%edge_vid(:,edge_containing_vertices)==v) .or. &
  !          all(this%edge_vid(:,edge_containing_vertices)==reverse (v))) &
  !          return
  !   end do
  !   edge_containing_vertices = 0 ! result for edge not found

  ! end function edge_containing_vertices

  subroutine print_data (this,normalized)

    class(polyhedron), intent(in) :: this
    logical, optional, intent(in) :: normalized

    integer :: v,e,f
    real(r8) :: x0(3),xl(3)
    logical :: normalizedh

    normalizedh = merge(normalized, .false., present(normalized))
    
    write(*,*) 'POLYHEDRON DATA:'
    
    if (allocated(this%x)) then
      if (normalizedh) then
        x0 = minval(this%x,dim=2)
        xl = maxval(this%x,dim=2) - x0
        do v = 1,this%nVerts
          write(*,'(a,i3,a,3es20.10)') 'x ',v,':  ',(this%x(:,v)-x0)/xl
        end do
      else
        do v = 1,this%nVerts
          write(*,'(a,i3,a,3es20.10)') 'x ',v,':  ',this%x(:,v)
        end do
      end if
      write(*,*)
    end if
    
    if (allocated(this%edge_vid)) then
      do e = 1,this%nEdges
        write(*,'(a,i3,a,2i4)') 'edge ',e,':  ',this%edge_vid(:,e)
      end do
      write(*,*)
    end if
    
    if (allocated(this%face_vid)) then
      do f = 1,this%nFaces
        write(*,'(a,i3,a,10i4)') 'face ',f,':  ',this%face_vid(:,f)
      end do
      write(*,*)
    end if

    write(*,'(a,es20.10)') 'volume ',this%vol
    write(*,*)

  end subroutine print_data

end module polyhedron_type
