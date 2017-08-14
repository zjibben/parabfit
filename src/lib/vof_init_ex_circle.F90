!! exactly initialize a vof field for a circle on a mesh of rectangles

module vof_init_ex_circle

  use kinds, only: r8
  use logging_services
  implicit none
  private

  public :: vof_init_circle

contains

  subroutine vof_init_circle(mesh, R, vof)

    use unstr_mesh_type

    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in) :: R
    real(r8), intent(out) :: vof(:,:)

    integer :: i

    do i = 1,mesh%ncell
      vof(2,i) = cell_vof(mesh%x(:,mesh%cnode(:,i)), R)
      vof(1,i) = 1 - vof(2,i)
    end do

  end subroutine vof_init_circle

  real(r8) function cell_vof(x,R)

    real(r8), intent(in) :: x(:,:), R

    integer :: ninside
    real(r8) :: xb(2,2), volume, xp
    logical :: is_inside(4)

    ! get the bounds & volume
    xb(:,1) = [minval(x(1,:)), minval(x(2,:))]
    xb(:,2) = [maxval(x(1,:)), maxval(x(2,:))]
    volume = (xb(1,2) - xb(1,1)) * (xb(2,2) - xb(2,1))

    ! switch to quadrant 1 if the cell doesn't cross quadrants
    if (count(xb(1,:) < 0) == 2) then
      xp = xb(1,1)
      xb(1,1) = -xb(1,2)
      xb(1,2) = -xp
    end if
    if (count(xb(2,:) < 0) == 2) then
      xp = xb(2,1)
      xb(2,1) = -xb(2,2)
      xb(2,2) = -xp
    end if

    ! check which points are inside the circle
    is_inside(1) = isInside(xb(:,1), R)
    is_inside(2) = isInside([xb(1,2), xb(2,1)], R)
    is_inside(3) = isInside(xb(:,2), R)
    is_inside(4) = isInside([xb(1,1), xb(2,2)], R)
    ninside = count(is_inside)

    select case (ninside)
    case (0)
      cell_vof = 0

      if (xb(2,1) < R .and. xb(1,1) < 0 .and. xb(1,2) > 0) then
        xp = sqrt(R**2 - xb(2,1)**2)
        cell_vof = antiderivative(xp,R) - antiderivative(-xp,R)
        cell_vof = cell_vof - 2 * xp * xb(2,1) ! subtract off the area outside the cell
      else if (xb(1,1) < R .and. xb(2,1) < 0 .and. xb(2,2) > 0) then
        xp = sqrt(R**2 - xb(1,1)**2)
        cell_vof = antiderivative(xp,R) - antiderivative(-xp,R)
        cell_vof = cell_vof - 2 * xp * xb(1,1) ! subtract off the area outside the cell
      end if
    case (1)
      xp = sqrt(R**2 - xb(2,1)**2)
      cell_vof = antiderivative(xp,R) - antiderivative(xb(1,1),R)
      cell_vof = cell_vof - xb(2,1) * (xp - xb(1,1)) ! subtract off the area outside the cell
    case (2)
      if (is_inside(2)) then
        ! integrate by x
        cell_vof = antiderivative(xb(1,2),R) - antiderivative(xb(1,1),R)
        cell_vof = cell_vof - xb(2,1) * (xb(1,2) - xb(1,1)) ! subtract off the area outside the cell
      else
        ! integrate by y
        cell_vof = antiderivative(xb(2,2),R) - antiderivative(xb(2,1),R)
        cell_vof = cell_vof - xb(1,1) * (xb(2,2) - xb(2,1)) ! subtract off the area outside the cell
      end if
    case (3)
      xp = sqrt(R**2 - xb(2,2)**2)
      cell_vof = antiderivative(xb(1,2),R) - antiderivative(xp,R)
      cell_vof = cell_vof - xb(2,1) * (xb(1,2) - xp) ! subtract off the area outside the cell
      cell_vof = cell_vof + (xp - xb(1,1)) * (xb(2,2) - xb(2,1)) ! add in the left rectangle
    case (4)
      cell_vof = volume
    end select
    cell_vof = cell_vof / volume

    if (cell_vof > 1 .or. cell_vof < 0) then

      print *, cell_vof
      print *, ninside, is_inside
      print *, xb(:,1)
      print *, xb(:,2)
      print *, xp

      call LS_Fatal ("vof init error")
    end if

    ! print *, is_inside
    ! print *, cell_vof
    ! stop

  end function cell_vof

  pure real(r8) function antiderivative(x, R)
    real(r8), intent(in) :: x, R
    antiderivative = (x * sqrt(R**2 - x**2) + R**2 * atan(x/sqrt(R**2 - x**2))) / 2
  end function antiderivative

  pure logical function isInside(x,R)
    real(r8), intent(in) :: x(:), R
    isInside = norm2(x) < R
  end function isInside

end module vof_init_ex_circle
