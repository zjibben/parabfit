!!
!! BITFIELD_TYPE: For when you need more bits than fit into an integer.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 26 May 2010
!!

module bitfield_type

  implicit none
  private

  !! A derived type that stores a large, but fixed, number of bits.
  !! As for integers, bits are indexed 0 through BIT_SIZE(<bitfield>)-1.
  !! Internally, the bits are stored in a fixed-size integer array.
  public :: bitfield

  !! Intrinsic bit manipulation procedures extended to the BITFIELD type.
  !! Procedures are elemental where the instrinsics are.
  public :: bit_size, btest, ibset, ibclr, iand, ieor, ior, not

  !! Equality/inequality operators overloaded for BITFIELD operands.
  public :: operator(==), operator(/=)

  !! Constant bitfield value with all bits set to 0.
  public :: ZERO_BITFIELD

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, parameter :: NUM_CHUNK = 4 ! increase this if you need more bits
  integer, parameter :: CHUNK_BITS = bit_size(1)  ! bits in a default integer
  integer, parameter :: BITFIELD_BITS = NUM_CHUNK * CHUNK_BITS

  type :: bitfield
    private
    integer :: chunk(0:NUM_CHUNK-1)
  end type bitfield

  type(bitfield), parameter :: ZERO_BITFIELD = bitfield(0)

  interface bit_size
    module procedure bit_size_bitfield, bit_size_bitfield1
  end interface

  interface ibset
    module procedure ibset_bitfield
  end interface

  interface ibclr
    module procedure ibclr_bitfield
  end interface

  interface btest
    module procedure btest_bitfield
  end interface

  interface iand
    module procedure iand_bitfield
  end interface

  interface ieor
    module procedure ieor_bitfield
  end interface

  interface ior
    module procedure ior_bitfield
  end interface

  interface not
    module procedure not_bitfield
  end interface

  interface operator(==)
    module procedure eq_bitfield
  end interface

  interface operator(/=)
    module procedure ne_bitfield
  end interface

contains

  integer function bit_size_bitfield (bf)
    type(bitfield), intent(in) :: bf
    bit_size_bitfield = BITFIELD_BITS
  end function bit_size_bitfield

  integer function bit_size_bitfield1 (bf)
    type(bitfield), intent(in) :: bf(:)
    bit_size_bitfield1 = BITFIELD_BITS
  end function bit_size_bitfield1

  elemental logical function btest_bitfield (bf, pos)
    type(bitfield), intent(in) :: bf
    integer, intent(in) :: pos
    integer :: n, p
    n = pos / CHUNK_BITS
    p = modulo(pos, CHUNK_BITS)
    btest_bitfield = btest(bf%chunk(n),p)
  end function btest_bitfield

  elemental function ibset_bitfield (bf, pos) result (bf_out)
    type(bitfield), intent(in) :: bf
    integer, intent(in) :: pos
    type(bitfield) :: bf_out
    integer :: n, p
    bf_out = bf
    n = pos / CHUNK_BITS
    p = modulo(pos, CHUNK_BITS)
    bf_out%chunk(n) = ibset(bf%chunk(n),p)
  end function ibset_bitfield

  elemental function ibclr_bitfield (bf, pos) result (bf_out)
    type(bitfield), intent(in) :: bf
    integer, intent(in) :: pos
    type(bitfield) :: bf_out
    integer :: n, p
    bf_out = bf
    n = pos / CHUNK_BITS
    p = modulo(pos, CHUNK_BITS)
    bf_out%chunk(n) = ibclr(bf%chunk(n),p)
  end function ibclr_bitfield

  elemental function iand_bitfield (bf1, bf2) result (bf_out)
    type(bitfield), intent(in) :: bf1, bf2
    type(bitfield) :: bf_out
    bf_out%chunk = iand(bf1%chunk, bf2%chunk)
  end function iand_bitfield

  elemental function ior_bitfield (bf1, bf2) result (bf_out)
    type(bitfield), intent(in) :: bf1, bf2
    type(bitfield) :: bf_out
    bf_out%chunk = ior(bf1%chunk, bf2%chunk)
  end function ior_bitfield

  elemental function ieor_bitfield (bf1, bf2) result (bf_out)
    type(bitfield), intent(in) :: bf1, bf2
    type(bitfield) :: bf_out
    bf_out%chunk = ieor(bf1%chunk, bf2%chunk)
  end function ieor_bitfield

  elemental function not_bitfield (bf) result (bf_out)
    type(bitfield), intent(in) :: bf
    type(bitfield) :: bf_out
    bf_out%chunk = not(bf%chunk)
  end function not_bitfield

  elemental logical function eq_bitfield (bf1, bf2)
    type(bitfield), intent(in) :: bf1, bf2
    eq_bitfield = all(bf1%chunk == bf2%chunk)
  end function eq_bitfield

  elemental logical function ne_bitfield (bf1, bf2)
    type(bitfield), intent(in) :: bf1, bf2
    ne_bitfield = any(bf1%chunk /= bf2%chunk)
  end function ne_bitfield

end module bitfield_type
