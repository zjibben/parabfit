module vof_io

  use kinds, only: r8
  implicit none

contains

  subroutine store_vof_field(filename, vof)

    character(*), intent(in) :: filename
    real(r8), intent(in) :: vof(:,:)

    integer :: fh

    open(newunit=fh, file=filename, action='write', access='stream')
    write(fh) vof
    close(fh)

  end subroutine store_vof_field

  subroutine read_vof_field(filename, vof)

    character(*), intent(in) :: filename
    real(r8), intent(out) :: vof(:,:)

    integer :: fh

    open(newunit=fh, file=filename, action='read', access='stream')
    read(fh) vof
    close(fh)

  end subroutine read_vof_field

  logical function file_exists(filename)
    character(*), intent(in) :: filename
    inquire(file=filename, exist=file_exists)
  end function file_exists

  subroutine store_normals(filename, normals)

    character(*), intent(in) :: filename
    real(r8), intent(in) :: normals(:,:,:)

    integer :: fh

    open(newunit=fh, file=filename, action='write', access='stream')
    write(fh) normals
    close(fh)

  end subroutine store_normals

  subroutine read_normals(filename, normals)

    character(*), intent(in) :: filename
    real(r8), intent(out) :: normals(:,:,:)

    integer :: fh

    open(newunit=fh, file=filename, action='read', access='stream')
    read(fh) normals
    close(fh)

  end subroutine read_normals

end module vof_io
