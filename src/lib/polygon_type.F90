!!
!! polygon_type
!!
!! This module defines an arbitrary polygon type, along with routines for
!! calculating area, splitting, locating intersections, etc.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!
!! References:
!!     1. Mirtich. Fast and Accurate Computation of Polehedral Mass Properties.
!!        Journal of Graphics Tools, 1996.
!!

#include "f90_assert.fpp"

module polygon_type

  use consts, only: ndim
  use kinds, only: r8
  use logging_services
  implicit none
  private

  type, public :: polygon
    integer :: nVerts
    real(r8), allocatable :: x(:,:), norm(:)
    real(r8) :: rot(3,3)
  contains
    procedure :: init => init_polygon
    procedure :: centroid
    procedure :: centroid2
    procedure :: area
    procedure :: area2
    procedure :: order
    procedure :: sort_order
    procedure :: basis
    procedure :: update_plane_normal
    procedure :: rotate_offset
    procedure :: print_data
  end type polygon

  type, public :: polygon_box
    integer :: n_elements
    type(polygon), allocatable :: elements(:)
  end type polygon_box

  public :: flat_polygon_box

contains

  subroutine init_polygon (this, x, norm)

    use array_utils, only: rotationMatrix

    class(polygon),     intent(out)   :: this
    real(r8),           intent(in)    :: x(:,:)
    real(r8), optional, intent(inout) :: norm(:)

    this%nVerts = size(x, dim=2)
    this%x = x

    if (present(norm)) then
      call this%update_plane_normal (norm)
    else
      call this%update_plane_normal ()
    end if

    this%rot = rotationMatrix(this%norm)

  end subroutine init_polygon

  subroutine rotate_offset (this, normal, xcen)

    use array_utils, only: rotationMatrix

    class(polygon), intent(inout) :: this
    real(r8), intent(in) :: normal(:), xcen(:)

    real(r8) :: R(3,3)
    integer :: i

    ! set up rotation matrix
    R = rotationMatrix(normal)

    ! rotate and offset vertices
    do i = 1,this%nVerts
      this%x(:,i) = matmul(R, this%x(:,i) - xcen)
    end do

    ! rotate normal vector
    this%norm = matmul(R, this%norm)

    ! update in-plane rotation matrix
    this%rot = rotationMatrix(this%norm)

  end subroutine rotate_offset

  ! calculate the normal via cross product from vectors defined by 3 vertices
  subroutine update_plane_normal (this,norm)

    use consts,        only: ndim
    use cell_geometry, only: cross_product
    use array_utils,   only: isZero,normalize, magnitude

    class(polygon),     intent(inout) :: this
    real(r8), optional, intent(inout) :: norm(:) ! return the newly calculated norm if it wasn't known

    integer :: i,j

    ! the direction of the normal is assumed from the node ordering (and assuming convex)

    ! if the polygon has >3 verteces, it could be non-planar and we should subdivide
    if (present(norm)) then
      ASSERT(size(norm)==ndim)
      this%norm = norm
    else
      if (allocated(this%norm)) deallocate(this%norm)
      allocate(this%norm(ndim))
      this%norm = 0
    end if

    if (all(isZero (this%norm))) then
      i = 3
      ! make sure we pick 3 vertices that don't all lie in the same plane
      do while (all(isZero (this%norm)) .and. i<=this%nVerts)
        this%norm = normalize(cross_product (this%x(:,2) - this%x(:,1), this%x(:,i) - this%x(:,1)))
        i = i + 1
      end do
      if (i>this%nVerts .and. all(isZero(this%norm))) then
        call this%print_data ()
        call LS_fatal ("polygon only consists of a line")
      end if
      if (present(norm)) norm = this%norm
    end if

  end subroutine update_plane_normal

  function centroid (this)
    use consts, only: ndim
    class(polygon), intent(in) :: this
    real(r8) :: centroid(ndim)
    centroid = sum(this%x, dim=2) / this%nVerts
  end function centroid

  function centroid2 (this) result(centroid)

    use consts, only: ndim
    use cell_geometry, only: cross_product

    class(polygon), intent(in) :: this

    integer :: i, j
    real(r8) :: centroid(ndim), c(2), xr(ndim,this%nVerts), a, tmp, area, &
        xc(ndim), xi(ndim), xj(ndim)

    ! ! get vertex coordinates in the plane of the polygon
    ! do i = 1,this%nVerts
    !   xr(:,i) = matmul(this%rot, this%x(:,i) - this%x(:,1))
    ! end do

    ! ! calculate the centroid in the plane of the polygon
    ! c = 0; a = 0
    ! do i = 1,this%nVerts
    !   j = modulo(i,this%nVerts) + 1
    !   tmp = xr(1,i)*xr(2,j) - xr(1,j)*xr(2,i)
    !   c = c + (xr(:2,i) + xr(:2,j)) * tmp
    !   a = a + tmp
    ! end do
    ! a = a / 2
    ! c = c / (6*a)

    ! ! rotate the centroid back into real space
    ! centroid = matmul(transpose(this%rot), [c(1), c(2), 0.0_r8]) + this%x(:,1)

    xc = this%centroid() ! translate to fake centroid to reduce floating point errors
    centroid = 0; area = 0
    do i = 1,this%nVerts
      j = modulo(i,this%nVerts) + 1
      xi = this%x(:,i) - xc
      xj = this%x(:,j) - xc

      tmp = dot_product(cross_product(xi, xj), this%norm)
      centroid = centroid + dot_product(xi + xj, this%norm) * tmp
      area = area + tmp
    end do
    area = area / 2
    centroid = centroid / (6*area) + xc

  end function centroid2

  ! calculate the area of a convex polygon
  ! assumes polygon vertices are ordered
  real(r8) function area (this)

    class(polygon), intent(in) :: this

    real(r8) :: xc(ndim), a, b, c, s
    integer :: v, w

    xc = this%centroid()

    ! polygon area is the sum of the areas associated
    ! with triangles connecting edges and the centroid
    area = 0
    do v = 1,this%nVerts
      ! next vertex
      w = modulo(v,this%nVerts) + 1

      ! triangle lengths
      a = norm2(xc - this%x(:,v))
      b = norm2(xc - this%x(:,w))
      c = norm2(this%x(:,w) - this%x(:,v))

      ! calculate the area using Heron's formula
      s = (a + b + c) / 2
      area = area + sqrt(s*(s-a)*(s-b)*(s-c))
    end do

  end function area

  function area2 (this) result(area)

    use cell_geometry, only: cross_product

    class(polygon), intent(in) :: this
    real(r8) :: area

    integer :: i, j
    real(r8) :: xi(ndim), xj(ndim), xc(ndim)

    xc = this%centroid()

    area = 0
    do i = 1,this%nVerts
      j = modulo(i,this%nVerts) + 1
      xi = this%x(:,i) - xc
      xj = this%x(:,j) - xc
      area = area + dot_product(cross_product(xi, xj), this%norm)
    end do
    area = area / 2

  end function area2

  ! order the vertices of a convex polygon
  !
  ! this is done by calculating the vector between each vertex and the polygon centroid,
  ! then the angle of that vector with respect to the x-axis in that space.
  ! this angle is used to sort the vertices
  subroutine order (this,array,index_sort)

    use array_utils, only: insertion_sort,xrange,invert,normalize,projectOnto,magnitude,isZero,&
        orthonormalBasis

    class(polygon),    intent(inout) :: this
    integer, optional, intent(inout) :: array(:) ! an array that gets sorted along with the polygon
    integer, optional, allocatable, intent(out) :: index_sort(:)

    real(r8), allocatable :: q(:,:)
    real(r8) :: xc(ndim), t(this%nVerts), t2(this%nVerts), prjx(2), xl(ndim,this%nVerts), tmp
    integer  :: i,ind(size(this%x,dim=2))

    ! calculate the location of the centroid, and the vertex coordinates with respect to the centroid
    ! normalize to avoid floating point cutoffs, since the polygon may be very tiny
    xc = this%centroid ()
    do i = 1,this%nVerts
      xl(:,i) = normalize(this%x(:,i) - xc(:))
    end do

    ! the projection coordinate directions
    ! WARNING: problems will occur here if the vertices are slightly non-planar
    q = orthonormalBasis(xl)
    if (.not.all(shape(q) >= [3,2])) then
      write(*,*)
      write(*,*) shape(q)
      do i = 1,size(q, dim=2)
        write(*,'(a,i3,a,3es20.10)') 'q  ',i,': ',q(:,i)
      end do
      write(*,*)
      do i = 1,this%nVerts
        write(*,'(a,i3,a,3es20.10)') 'xl ',i,': ',xl(:,i)
      end do
      write(*,*)
      call this%print_data ()
      call LS_fatal ("polygon ordering failed: unable to calculate polygon-plane coordinates")
    end if

    ! calculate the rotation angle
    !t(1) = 0.0_r8
    do i = 1,this%nVerts
      ! get coordinates for the vertex in the 2D plane defined by the polygon
      prjx(1) = dot_product(xl(:,i),q(:,1))
      prjx(2) = dot_product(xl(:,i),q(:,2))

      ! find the angle made by this vertex with respect to the first vertex
      t(i) = atan2(prjx(2),prjx(1))
      !write(*,'(i6, 3es20.10)') i, prjx, t(i)
    end do

    ! sort based on angle
    if (present(array) .or. present(index_sort)) then
      t2 = t
      ind = xrange (1,size(this%x,dim=2))
      call insertion_sort (ind,t2)
      if (present(index_sort)) index_sort = ind
      ind = invert (ind)
      do i = 1,size(array)
        if (array(i)>0) array(i) = ind(array(i))
      end do
    end if

    call insertion_sort (this%x,t)

    call this%update_plane_normal ()

  end subroutine order

  ! order the vertices of a convex polygon
  !
  ! this is done by calculating the vector between each vertex and the polygon centroid,
  ! then the angle of that vector with respect to the x-axis in that space.
  ! this angle is used to sort the vertices
  function sort_order(this)

    use array_utils, only: insertion_sort,xrange,normalize,orthonormalBasis

    class(polygon), intent(in) :: this
    integer, allocatable :: sort_order(:)

    real(r8), allocatable :: q(:,:)
    real(r8) :: xc(ndim), t(this%nVerts), prjx(2), xl(ndim,this%nVerts)
    integer  :: i

    ! calculate the location of the centroid, and the vertex coordinates with respect to the centroid
    ! normalize to avoid floating point cutoffs, since the polygon may be very tiny
    xc = this%centroid ()
    do i = 1,this%nVerts
      xl(:,i) = normalize(this%x(:,i) - xc(:))
    end do

    ! the projection coordinate directions
    ! WARNING: problems will occur here if the vertices are slightly non-planar
    q = orthonormalBasis(xl)
    if (.not.all(shape(q) >= [3,2])) then
      write(*,*)
      write(*,*) shape(q)
      do i = 1,size(q, dim=2)
        write(*,'(a,i3,a,3es20.10)') 'q  ',i,': ',q(:,i)
      end do
      write(*,*)
      do i = 1,this%nVerts
        write(*,'(a,i3,a,3es20.10)') 'xl ',i,': ',xl(:,i)
      end do
      write(*,*)
      call this%print_data ()
      call LS_fatal ("polygon ordering failed: unable to calculate polygon-plane coordinates")
    end if

    ! calculate the rotation angle
    !t(1) = 0.0_r8
    do i = 1,this%nVerts
      ! get coordinates for the vertex in the 2D plane defined by the polygon
      prjx(1) = dot_product(xl(:,i),q(:,1))
      prjx(2) = dot_product(xl(:,i),q(:,2))

      ! find the angle made by this vertex with respect to the first vertex
      t(i) = atan2(prjx(2),prjx(1))
      !write(*,'(i6, 3es20.10)') i, prjx, t(i)
    end do

    sort_order = xrange (1,size(this%x,dim=2))
    call insertion_sort (sort_order,t)

  end function sort_order

  ! generate an orthogonal basis for the polygon, approximately scaled to the size of the polygon
  ! the first element is the shortest dimension, the second the longest
  function basis (this)

    use array_utils, only: projectOnto, magnitude

    class(polygon), intent(in) :: this
    real(r8) :: basis(ndim,2)

    integer :: i, iN
    real(r8) :: xc(ndim), xl(ndim), minlen, maxlen, length

    ! find the minimum and maximum directions
    maxlen = 0.0_r8
    minlen = huge(1.0_r8)
    xc = this%centroid ()
    ! do i = 1,this%nVerts
    !   xl = xc - this%x(:,i)
    !   length = magnitude(xl)

    !   if (length < minlen) then
    !     minlen = length
    !     basis(:,1) = xl
    !   else if (length > maxlen) then
    !     maxlen = length
    !     basis(:,2) = xl
    !   end if
    ! end do

    ! find closest and farthest points to the centroid on the polygon edge
    do i = 1,this%nVerts
      iN = modulo(i,this%nVerts) + 1
      xl = this%x(:,i) + projectOnto(xc - this%x(:,i), this%x(:,iN) - this%x(:,i)) - xc
      length = magnitude(xl)

      if (length < minlen) then
        minlen = length
        basis(:,1) = xl
      else if (length > maxlen) then
        maxlen = length
        basis(:,2) = xl
      end if
    end do

    ! align basis(:,2) such that the two vectors are orthogonal
    ! it should be almost orthogonal as is
    basis(:,2) = basis(:,2) - projectOnto(basis(:,2),basis(:,1))

    ! call this%print_data ()
    ! print '(a,3es15.5)', 'b: ',basis(:,1)
    ! print '(a,3es15.5)', 'b: ',basis(:,2)
    ! print *

  end function basis

  subroutine print_data (this)

    class(polygon), intent(in) :: this

    integer :: v

    print *, 'POLYGON DATA:'
    if (allocated(this%x)) then
      do v = 1,this%nVerts
        print '(a,i3,a,3es15.5)', 'x ',v,':  ',this%x(:,v)
      end do
      write(*,*)
    end if

    print '(a,3es15.5)', 'norm ',this%norm

  end subroutine print_data

  type(polygon_box) function flat_polygon_box(polygon_boxes)

    type(polygon_box), intent(in) :: polygon_boxes(:)

    integer :: i, j

    flat_polygon_box%n_elements = sum(polygon_boxes%n_elements)
    allocate(flat_polygon_box%elements(flat_polygon_box%n_elements))

    j = 0
    do i = 1,size(polygon_boxes)
      if (polygon_boxes(i)%n_elements > 0) then
        flat_polygon_box%elements(j+1:j+polygon_boxes(i)%n_elements) = polygon_boxes(i)%elements
        j = j + polygon_boxes(i)%n_elements
      end if
    end do

  end function flat_polygon_box

end module polygon_type
