#include "f90_assert.fpp"

module unstr_mesh_partition_type

  use kinds, only: r8
  use logging_services
  use unstr_mesh_type, only : unstr_mesh
  use timer_tree_type
  implicit none
  private
  
  type, public :: unstr_mesh_partition
     private

     !! For the Fortran experts: I'd like to make this protected, that
     !! is I'd like to expose it as a read-only field within this
     !! type.  I could do that by creating a function with the same
     !! arguments as the array, however most compilers cannot see
     !! enough to inline and vectorize.
     integer,allocatable,public :: cface(:,:,:)

     !! Private internal state
     type(unstr_mesh),pointer :: mesh => null() 
     integer :: npart
     integer :: psize
     integer, allocatable :: exp_face_offset(:)

   contains

     !! Provide a means to initialize the mesh from a parameter list
     !! of an explicit size, the explicit size is used in our unit
     !! tests.
     procedure :: ump_init_size
     procedure :: ump_init_params
     generic :: init => ump_init_size, ump_init_params

     !! Scatter and gather-sum face data.
     procedure :: scatter_face => ump_scatter_face
     procedure :: gather_sum_face => ump_gather_sum_face

     !! The number of partitions
     procedure :: num_part => ump_num_part

     !! The maximum size of the partiiton
     procedure :: part_size => ump_part_size

     !! The bounds for a particular partition
     procedure :: part_start => ump_part_start
     procedure :: part_end => ump_part_end

     !! Implementation
     procedure, private :: cell_to_part => ump_cell_to_part
     procedure, private :: find_shared_faces => ump_find_shared_faces
     procedure, private :: create_cface => ump_create_cface
  end type unstr_mesh_partition
  
contains

  subroutine ump_init_size(this, mesh, psize)

    class(unstr_mesh_partition), intent(out) :: this
    type(unstr_mesh), target, intent(in) :: mesh
    integer, intent(in) :: psize

    integer, allocatable :: shared_faces(:,:)
    integer :: ncface,c,f

    call start_timer('partition init')
    
    !! Grab a reference to the mesh
    this%mesh => mesh
    
    !! Record the partition_size
    this%psize = psize

    !! Record the partition size and compute the number of partitions.
    INSIST(this%psize.gt.0)
    this%npart = mesh%ncell/this%psize
    if (0.ne.mod(mesh%ncell,this%psize)) this%npart = this%npart + 1

    !! Find faces that join cells in different partitions.
    ncface = size(this%mesh%cface,dim=1)
    allocate(shared_faces(ncface,mesh%ncell))
    call this%find_shared_faces(shared_faces)
    
    !! Compute the offsets for the expanded face array. 
    allocate(this%exp_face_offset(1+mesh%nface))
    this%exp_face_offset = 1
    do c=1,mesh%ncell
       do f=1,ncface
          if (0.ne.shared_faces(f,c)) then
             this%exp_face_offset(shared_faces(f,c)) = &
                  this%exp_face_offset(shared_faces(f,c)) + 1
          end if
       end do
    end do

    call in_place_incl_scan(this%exp_face_offset)
    
    !! this%exp_face_offset is now the mapping frome a global face
    !! array to a partitioned face array, where faces are duplicated
    !! if they are shared between partitions.  The list of entries in
    !! a partitioned face array, for a face whose global index is f
    !! is given by this%exp_face_offset(f):this%exp_face_offset(f+1)-1

    !! Create the mapping from partition, cell-in-partition,
    !! face-in-cell to the index of that face in the partitioned face
    !! array.
    call this%create_cface()

    write(*,*) 'Fraction increase in faces for partitioning', &
         1.0d0*(this%exp_face_offset(mesh%nface+1)-1)/this%mesh%nface
    write(*,*) 'Number of partitions', this%num_part()

    call stop_timer('partition init')
    
  end subroutine ump_init_size
  
  subroutine ump_init_params(this, mesh, params) 

    use parameter_list_type
    
    class(unstr_mesh_partition), intent(out) :: this
    type(unstr_mesh), target, intent(in) :: mesh
    type(parameter_list), intent(inout) :: params

    integer :: stat
    character(:), allocatable :: errmsg,context
    integer :: psize

    context = 'processing ' // params%name() // ': '

    if (params%is_parameter('partition-size')) then 
       call params%get('partition-size', psize, stat=stat, errmsg=errmsg)
       if (stat /= 0) call LS_fatal (context//errmsg)
    else
       call LS_fatal ('could not find "partition-size" in ' // params%name() )
    end if

    call ump_init_size(this, mesh, psize)

  end subroutine ump_init_params

  subroutine ump_find_shared_faces(this, shared_face)

    use unstr_mesh_tools, only: get_cell_neighbor_array

    class(unstr_mesh_partition), intent(in) :: this
    integer, intent(out) :: shared_face(:,:)

    integer :: ncface,stat,c,f,p1,p2
    integer, allocatable :: cnhbr(:,:)

    !! Generate cell neighor list
    ncface = size(this%mesh%cface,dim=1)
    allocate(cnhbr(ncface,this%mesh%ncell))
    call get_cell_neighbor_array (this%mesh%cnode, cnhbr, stat)
    if (stat /= 0) call LS_fatal ('bad mesh topology detected')

    !! Use the cell neighbor list to find the shared faces.
    shared_face = 0
    do c=1,this%mesh%ncell       
       p1 = this%cell_to_part(c)
       do f=1,ncface
          if (0.ne.cnhbr(f,c)) then
             p2 = this%cell_to_part(cnhbr(f,c))
             if (p2.lt.p1) then 
                shared_face(f,c) = find_shared_face(c,cnhbr(f,c))
             end if
          end if
       end do
    end do

  contains 
    !! Find the index of a shared face between cells ci and cj
    function find_shared_face(ci,cj) result(f)
      integer :: f
      integer,intent(in) :: ci,cj
      integer :: i,j,ncface
      f = 0
      ncface = size(this%mesh%cface,dim=1)
      do i=1,ncface
         do j=1,ncface
            if (this%mesh%cface(i,ci).eq.this%mesh%cface(j,cj) & 
                 .and.this%mesh%cface(i,ci).ne.0) f = this%mesh%cface(i,ci)
         end do
      end do
      ASSERT(f.ne.0)
    end function find_shared_face
  end subroutine ump_find_shared_faces

  !! Perform an in-place exclusive scan
  subroutine in_place_incl_scan(a)
    integer, intent(inout) :: a(:)
    integer :: aa, j, sum
    aa = a(1)
    sum = aa + 1
    do j=2,size(a)
       aa = a(j)
       a(j) = sum
       sum = sum + aa
    end do
  end subroutine in_place_incl_scan

  subroutine ump_create_cface(this)
    class(unstr_mesh_partition), intent(inout) :: this

    integer :: ncface,p,st,en,c,f,fid
    logical, allocatable :: tag(:)

    !! Allocate and initialize the tag array. 
    allocate(tag(this%exp_face_offset(this%mesh%nface+1)-1))
    tag = .false.

    !! Allocate the array that maps from partiton, cell and face into
    !! a partitioned array of face data.
    ncface = size(this%mesh%cface,dim=1)
    allocate(this%cface(ncface,this%psize,this%npart))

    !! Populate this array.  Because this array will be accessed from
    !! within loops over meshes, we would like to initialize it in a
    !! way that respects 'first-touch'.  It isn't clear how much of an
    !! effect first-touch will have, so this is more about erring on
    !! the side of safety.
    do p=1,this%num_part()
       st = this%part_start(p)
       en = this%part_end(p)
       do c=st,en
          do f=1,ncface
            fid = this%mesh%cface(f,c)
             !! Is this face shared between partitions
             if (this%exp_face_offset(fid).ne.this%exp_face_offset(fid+1)-1) then
                !! Yes it is
                if (tag(this%exp_face_offset(fid))) then
                   ASSERT(.not. tag(this%exp_face_offset(fid+1)-1))
                   this%cface(f,c-st+1,p) = this%exp_face_offset(fid+1)-1
                   tag(this%exp_face_offset(fid+1)-1) = .true.
                else
                   this%cface(f,c-st+1,p) = this%exp_face_offset(fid)
                   tag(this%exp_face_offset(fid)) = .true.
                end if
             else
               !! No, this face is not shared 
               this%cface(f,c-st+1,p) = this%exp_face_offset(fid)
             end if
          end do
       end do
    end do

  end subroutine ump_create_cface

  subroutine ump_scatter_face(this, gfd, pfd) 
    class(unstr_mesh_partition), intent(in) :: this
    real(r8), intent(in) :: gfd(:)
    real(r8), allocatable, intent(out) :: pfd(:)

    integer :: exp_size,p,j,k,m,ncface,fid

    call start_timer('partition scatter')
    INSIST(size(gfd).eq.this%mesh%nface)

    !! Make sure the partitioned face array is allocated with sufficent
    !! size.
    exp_size = this%exp_face_offset(this%mesh%nface+1)-1
    if (allocated(pfd)) then 
       if (size(pfd).le.exp_size) deallocate(pfd)
    end if
    if (.not.allocated(pfd)) allocate(pfd(exp_size))
    
    !! Scatter the data from the global face array to the partitioned
    !! face array.  So that first-touch effects place data in the
    !! desired NUMA domain, we duplicate the iteration pattern that
    !! will be used by the client code.

    ncface = size(this%mesh%cface,dim=1)
    !$OMP parallel private(p,j,k,fid)
    !$OMP do 
    do p=1,this%num_part() ! partitions
       do j = this%part_start(p),this%part_end(p) ! cells in the partition p
          do k = 1,ncface ! faces of cell j
             fid = this%mesh%cface(k,j)
             pfd(this%exp_face_offset(fid):this%exp_face_offset(fid+1)-1) = gfd(fid)
          end do
       end do
    end do
    !$OMP end do
    !$OMP end parallel
    call stop_timer('partition scatter')
  end subroutine ump_scatter_face

  subroutine ump_gather_sum_face(this, pfd, gfd)
    class(unstr_mesh_partition), intent(in) :: this
    real(r8), intent(in) :: pfd(:)
    real(r8), intent(out) :: gfd(:)

    integer :: exp_size,p,j,k,m,ncface,fid

    call start_timer('partition gather')
    exp_size = this%exp_face_offset(this%mesh%nface+1)-1
    INSIST(size(pfd).eq.exp_size)
    INSIST(size(gfd).eq.this%mesh%nface)

    !$OMP parallel
    !$OMP do
    do fid=1,this%mesh%nface
       gfd(fid) = sum(pfd(this%exp_face_offset(fid):this%exp_face_offset(fid+1)-1))
    end do
    !$OMP end do
    !$OMP end parallel

    call stop_timer('partition gather')
  end subroutine ump_gather_sum_face

  function ump_num_part(this) result(n)
    integer :: n
    class(unstr_mesh_partition), intent(in) :: this
    n = this%npart
  end function ump_num_part

  function ump_part_size(this) result(n)
    integer :: n
    class(unstr_mesh_partition), intent(in) :: this
    n = this%psize
  end function ump_part_size

  function ump_part_start(this, p) result(s)
    integer :: s
    class(unstr_mesh_partition), intent(in) :: this
    integer, intent(in) :: p
    INSIST(p.ge.1.and.p.le.this%npart)
    s = 1+(p-1)*this%psize
    ASSERT(s.gt.0)
    ASSERT(s.le.this%mesh%ncell)
  end function ump_part_start

  function ump_part_end(this, p) result(e)
    integer :: e
    class(unstr_mesh_partition), intent(in) :: this
    integer, intent(in) :: p
    INSIST(p.gt.0.and.p.le.this%npart)
    e = min(p*this%psize,this%mesh%ncell)
    ASSERT(e.gt.0)
    ASSERT(e.le.this%mesh%ncell)
  end function ump_part_end

  !! Given a cell index, get the partition index for the partition to
  !! which the cell belongs.
  function ump_cell_to_part(this, c) result(p)
    integer p
    class(unstr_mesh_partition), intent(in) :: this
    integer, intent(in) :: c
    INSIST(c.gt.0.and.c.le.this%mesh%ncell)
    p = c / this%psize
    if (0.ne.mod(c,this%psize)) p = p + 1 
    ASSERT(p.le.this%npart)
  end function ump_cell_to_part

end module unstr_mesh_partition_type
