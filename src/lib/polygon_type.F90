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
!!     1. Mirtich. Fast and Accurate Computation of Polehedral Mass Properties. Journal of Graphics Tools, 1996.
!! 

module polygon_type
  use kinds, only: r8
  use logging_services
  implicit none
  private

  type, public :: polygon
    integer               :: nVerts
    real(r8), allocatable :: x(:,:)
    real(r8)              :: norm(3)
  contains
    procedure          :: init => init_polygon
    procedure          :: intXdA
    procedure          :: centroid
    procedure          :: order
    procedure          :: update_plane_normal
    procedure          :: print_data
  end type polygon
  
contains

  subroutine init_polygon (this, x, norm)
    class(polygon),     intent(out)   :: this
    real(r8),           intent(in)    :: x(:,:)
    real(r8), optional, intent(inout) :: norm(:)

    this%nVerts = size(x, dim=2)
    allocate(this%x(3,this%nVerts))
    this%x = x

    if (present(norm)) then
      call this%update_plane_normal (norm)
    else
      call this%update_plane_normal ()
    end if
    
  end subroutine init_polygon

  ! calculate the normal via cross product from vectors defined by 3 vertices
  subroutine update_plane_normal (this,norm)
    use cell_geometry, only: cross_product
    use array_utils,   only: isZero

    class(polygon),     intent(inout) :: this
    real(r8), optional, intent(inout) :: norm(:) ! return the newly calculated norm if it wasn't known
    
    integer :: i,j

    ! the direction of the normal is assumed from the node ordering (and assuming convex)

    ! if the polygon has >3 verteces, it could be non-planar and we should subdivide
    if (present(norm)) then
      this%norm = norm
    else
      this%norm = 0.0_r8
    end if

    if (all(isZero (this%norm))) then
      i = 3
      ! make sure we pick 3 vertices that don't all lie in the same plane
      do while (all(isZero (this%norm)) .and. i<=this%nVerts)
        this%norm = cross_product ( this%x(:,2) - this%x(:,1), this%x(:,i) - this%x(:,1))
        this%norm = this%norm / sqrt(sum(this%norm**2)) ! normalize
        i = i + 1
      end do
      if (i>this%nVerts .and. all(isZero(this%norm))) then
        write(*,*) 'nVerts',this%nVerts
        do j = 1,this%nVerts
          write(*,'(a,i3,3es14.4)') 'vertex ',j,this%x(:,j)
        end do
        call LS_fatal ("polygon only consists of a line")
      end if
      if (present(norm)) norm = this%norm
    end if

  end subroutine update_plane_normal
  
  !
  ! This function calculates the integral of X, Y, or Z (input dir)
  ! over a polygon by reducing it to the sum of the integral of
  ! dir^2 over all edges, which can be calculated algebraically.
  ! Follows the algorithm proposed by [1].
  !
  real(r8) function intXdA (this,dir)
    class(polygon), intent(in) :: this
    integer,        intent(in) :: dir

    integer :: e,eN,nEdges,A,B,C

    nEdges = size(this%x, dim=2)

    ! project into the plane that maximizes the area
    C = maxloc(abs(this%norm), dim=1)
    A = mod(C,3)+1
    B = mod(A,3)+1
    
    ! sum up the integral of dir^2 over all edges
    intXdA = 0.0_r8
    do e = 1,nEdges
      eN = mod(e,nEdges)+1 ! pick the next edge, looping back to the beginning if we are at the end
      
      if (A==dir) then
        intXdA = intXdA + (this%x(B,eN)-this%x(B,e)) * (this%x(A,eN)**2 + this%x(A,eN)*this%x(A,e) + this%x(A,e)**2) &
             / (6.0_r8*this%norm(C))
      else if (B==dir) then
        intXdA = intXdA - (this%x(A,eN)-this%x(A,e)) * (this%x(B,eN)**2 + this%x(B,eN)*this%x(B,e) + this%x(B,e)**2) &
             / (6.0_r8*this%norm(C))
      else if (C==dir) then
        intXdA = intXdA - ( &
             + this%norm(A) * (this%x(B,eN)-this%x(B,e)) * (this%x(A,eN)**2 + this%x(A,eN)*this%x(A,e) + this%x(A,e)**2) / 6.0_r8 &
             - this%norm(B) * (this%x(A,eN)-this%x(A,e)) * (this%x(B,eN)**2 + this%x(B,eN)*this%x(B,e) + this%x(B,e)**2) / 6.0_r8 &
             - sum(this%norm*this%x(:,e)) * (this%x(B,eN)-this%x(B,e)) * (this%x(A,eN) + this%x(A,e)) / 2.0_r8 &
             ) / this%norm(C)**2
      end if

    end do

  end function intXdA

  function centroid (this)
    class(polygon), intent(in) :: this
    real(r8)                   :: centroid(3)
    centroid = sum(this%x, dim=2) / real(this%nVerts,r8)
  end function centroid

  ! order the vertices of a convex polygon
  !
  ! this is done by calculating the vector between each vertex and the polygon centroid,
  ! then the angle of that vector with respect to the x-axis in that space.
  ! this angle is used to sort the vertices
  subroutine order (this,array)
    use array_utils, only: insertion_sort,xrange,invert,magnitude,projectOnto

    class(polygon),    intent(inout) :: this
    integer, optional, intent(inout) :: array(:)
    
    real(r8) :: xc(3), t(this%nVerts), t2(this%nVerts), prjx(2), xl(3,this%nVerts), magxl1, q(3,2)
    integer  :: i,ind(size(this%x,dim=2))

    ! calculate the location of the centroid, and the vertex coordinates
    ! with respect to the centroid
    xc = this%centroid ()
    do i = 1,this%nVerts
      xl(:,i) = this%x(:,i) - xc(:)
    end do

    ! the projection coordinate directions
    ! WARNING: need to ensure the second vector here is not parallel to q1
    q(:,1) = xl(:,1) / magnitude (xl(:,1))
    q(:,2) = xl(:,2) - projectOnto (xl(:,2),q(:,1))
    q(:,2) = q(:,2) / magnitude (q(:,2))
    
    ! calculate the rotation angle
    t(1) = 0.0_r8
    do i = 2,this%nVerts
      ! get coordinates for the vertex in the 2D plane defined by the polygon
      prjx(1) = sum(xl(:,i)*q(:,1))
      prjx(2) = sum(xl(:,i)*q(:,2))
      
      ! find the angle made by this vertex with respect to the first vertex
      t(i) = atan2(prjx(2),prjx(1))
    end do
    
    ! sort based on angle
    if (present(array)) then
      t2 = t
      ind = xrange (1,size(this%x,dim=2))
      call insertion_sort (ind,t2)
      ind = invert (ind)
      do i = 1,size(array)
        if (array(i)>0) array(i) = ind(array(i))
      end do
    end if

    call insertion_sort (this%x,t)

    call this%update_plane_normal ()

  end subroutine order

  subroutine print_data (this)
    class(polygon), intent(in) :: this

    integer :: v

    write(*,*) 'POLYGON DATA'
    do v = 1,this%nVerts
      write(*,'(a,i3,a,3es20.10)') 'x ',v,':  ',this%x(:,v)
    end do
    write(*,*)

    write(*,'(a,3f14.4)') 'norm ',this%norm
    write(*,*)

  end subroutine print_data

end module polygon_type
