!=======================================================================
! Purpose(s):
!
!   Define procedures necessary to locate the planar interfaces
!   in all interface cells.
!
!   Public Interface:
!
!     * call LOCATE_PLANE ()
!
!       Compute the constant Rho in the interface plane equation
!       X*Normal - Rho = 0 (given the plane Normal) subject to a
!       volume conservation constraint. The interface constant
!       Rho effectively locates the plane in the cell.
!
! Contains: LOCATE_PLANE
!           RHO_BRACKET
!           RHO_BRENT
!
! Author(s): Stewart J. Mosso, LANL X-HM (sjm@lanl.gov)
!            Douglas B. Kothe, LANL (dbk@lanl.gov)
!
!=======================================================================


module locate_plane_module
  use kinds, only: r8
  use hex_types, only: reconstruction_hex
  use logging_services
  implicit none
  private

  integer, parameter :: ndim = 3
  integer, parameter :: nfc = 6 ! number of faces per cell
  integer, parameter :: nvc = 8 ! number of vertices per cell
  
  type, extends(reconstruction_hex), public :: locate_plane_hex
  contains
    procedure :: unit_test
    procedure :: init => init_locate_plane_hex
    procedure :: locate_plane
    procedure, private :: rho_bracket
    procedure, private :: rho_brent
  end type locate_plane_hex

contains

  logical function unit_test (this, normal, vof, rhoex)
    class(locate_plane_hex), intent(out) :: this
    real(r8), intent(in) :: normal(3), vof, rhoex

    integer  :: iter
    real(r8) :: flux_vol_coord(3,8)

    ! set the nodes
    this%node(:,1) = [0.0_r8, 0.0_r8, 0.0_r8]
    this%node(:,2) = [1.0_r8, 0.0_r8, 0.0_r8]
    this%node(:,3) = [1.0_r8, 1.0_r8, 0.0_r8]
    this%node(:,4) = [0.0_r8, 1.0_r8, 0.0_r8]
    this%node(:,5) = [0.0_r8, 0.0_r8, 1.0_r8]
    this%node(:,6) = [1.0_r8, 0.0_r8, 1.0_r8]
    this%node(:,7) = [1.0_r8, 1.0_r8, 1.0_r8]
    this%node(:,8) = [0.0_r8, 1.0_r8, 1.0_r8]

    this%volume = this%calc_volume () ! set the volume
    this%normal = normal              ! set the normal
    this%vof = vof                    ! set the vof

    ! calculate the plane constant
    call this%locate_plane (iter)

    ! put a check here rather than a print
    write(*,'(2(a,f14.10),L)') 'plane constant = ',this%rho,',     correct plane constant = ',rhoex
    
    unit_test = (this%rho==rhoex)
    
  end function unit_test

  subroutine init_locate_plane_hex (this, norm, vof, volume, node)
    class(locate_plane_hex), intent(out) :: this
    real(r8),                intent(in)  :: norm(3), vof, volume, node(3,8)

    this%normal = norm
    this%vof    = vof
    this%volume = volume
    this%node   = node

  end subroutine init_locate_plane_hex

  !=======================================================================
  ! Purpose(s):
  !
  !   Given the equation for interface planes: X*Normal - Rho = 0,
  !   compute the parameter Rho (given the plane Normal) subject
  !   to a volume conservation constraint.  This is done in
  !   an iterative fashion, using either Brent's
  !   method.  The iteration is first bracketed by finding the
  !   interval [Rho_Min,Rho_Max] bracketing Rho.  Given Rho and the
  !   volume flux geometry, the plane equation parameters
  !   (Rho, Normal) are then used to compute the truncated
  !   flux volume (Vp).
  !
  !=======================================================================
  subroutine locate_plane (this, iter)
    use timer_tree_type
    use truncate_volume_module, only: truncvol_data

    class(locate_plane_hex), intent(inout) :: this
    integer,                 intent(out)   :: iter

    real(r8) :: Rho_Min, Rho_Max, V_Min, V_Max
    type(truncvol_data) :: trunc_vol(nfc)
    
    ! Start the locate plane timer.
    ! WARNING: Need to figure out how to make this run in parallel with OpenMP.
    !          Currently, the timer is a global variable.
    !call start_timer ("Locate Plane")

    ! Bracket the correct value of Rho [Rho_Min,Rho_Max] to insure
    ! the subsequent iteration will converge efficiently
    call rho_bracket (this, Rho_Min, Rho_Max, V_Min, V_Max, trunc_vol)

    ! Compute Rho, which parameterizes the plane
    ! characterized by the equation Normal*X - Rho = 0
    ! using Brents method iteration
    call rho_brent (this, iter, Rho_Min, Rho_Max, V_Min, V_Max, trunc_vol)

    ! Stop the locate plane timer.
    !call stop_timer("Locate Plane")

  end subroutine locate_plane

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  !=======================================================================
  ! PURPOSE -
  !   Compute the plane parameter Rho in the plane equation:
  !                    X*Normal - Rho = 0
  !   iteratively using Brents method as documented in Numerical
  !   Recipes (Press, Flannery, Teukolsky and Vetterling; Cambridge
  !   University Press, 1986), p. 253
  !=======================================================================
  subroutine rho_bracket (this, Rho_Min, Rho_Max, V_Min, V_Max, trunc_vol)
    use truncate_volume_module, only: truncate_face,face_param,truncvol_data

    class(locate_plane_hex), intent(inout) :: this
    real(r8), intent(out) :: Rho_Min, Rho_Max, V_min, V_max
    type(truncvol_data), intent(inout) :: trunc_vol(nfc)

    integer  :: f, n, v
    real(r8) :: V_v(nvc), Volume, Rho

    ! Initialize the bracketed values of truncation volume [V_Min,V_Max]
    ! and the associated plane constant [Rho_Min,Rho_Max]. 
    V_v = 0.0_r8

    ! Compute volumes truncated by the interface passing through each of the vertices.
    FACE_LOOP: do f = 1,nfc

      ! Compute and store face parameters.
      trunc_vol(f) = face_param (this, 'full_cell', f)

      VERTEX_LOOP: do v = 1,nvc

        ! Get the value of Rho for this vertex.
        This%Rho = sum(This%node(:,v) * This%Normal(:))

        ! Get this face's contribution to the truncation volume.
        Volume = truncate_face (this, trunc_vol(f))

        V_v(v) = V_v(v) + Volume

      end do VERTEX_LOOP

    end do FACE_LOOP

    ! Set up the limits on the volume bracketing. Make sure that
    ! V_v is bounded zero = V_Min <= V_v <= V_Max = Cell Volume.
    ! If V_v > Cell Volume (because of roundoff), set V_v = Cell Volume

    V_Min = 0.0_r8
    V_Max = This%Volume
    Rho_Max = -HUGE(0.0_r8)
    Rho_Min =  HUGE(0.0_r8)

    ! If any V_v is outside of allowed bounds, force V_min <= V_v <= V_max
    do v = 1,nvc
      if (V_v(v) > V_Max )  V_v(v) = V_Max
      if (V_v(v) < 0.0_r8)  V_v(v) = V_Min
      Rho_Max = MAX(Rho_Max,V_v(v))
      Rho_Min = MIN(Rho_Min,V_v(v))
    end do

    ! Make sure that V_v = V_max for at least one vertex.
    ! If not, then set the maximum V_v to V_max.
    if (Rho_Max < V_Max) then
      do v = 1,nvc
        if (V_v(v) == Rho_Max) V_v(v) = V_Max
      end do
    end if

    ! Make sure that V_v = V_min for at least one vertex.
    ! If not, then set the minimum V_v to V_min.
    if (Rho_Min > V_Min) then
      do v = 1,nvc
        if (V_v(v) == Rho_Min) V_v(v) = V_Min
      end do
    end if

    ! Now bracket Rho.
    Rho_Max  =  HUGE(0.0_r8)
    Rho_Min  = -HUGE(0.0_r8)

    do v = 1,nvc

      ! Get the value of Rho for this vertex.
      This%Rho = sum(This%node(:,v) * This%Normal(:))

      if ((V_v(v) <= This%Vof*This%Volume .and. V_v(v) > V_Min) &
           .or. V_v(v) == V_Min) then
        Rho_Min = This%Rho
        V_Min = V_v(v)
      end if

      if ((V_v(v) > This%Vof*This%Volume .and. V_v(v) < V_Max) &
           .or. V_v(v) == V_Max) then
        Rho_Max = This%Rho
        V_Max = V_v(v)
      end if

    end do

    ! At this point, [V_Min,V_Max] brackets V_v, and [Rho_Min,Rho_Max]
    ! are those plane constants corresponding to [V_Min,V_Max]. The
    ! bracketed volumes are monotonic, but the Rho's aren't necessarily
    ! monotonic, so make sure that Rho_Min < R_Max. Be sure to interchange
    ! V_Min and V_Max if R_Min and R_Max have to be interchanged.

    if (Rho_Max < Rho_Min) then
      Rho = Rho_Max
      Volume = V_Max
      Rho_Max = Rho_Min
      V_Max = V_Min
      Rho_Min = Rho
      V_Min = Volume
    end if

    ! If we haven't properly set Rho everywhere, punt.
    if (Rho_Min == -HUGE(0.0_r8) .or. Rho_Max == HUGE(0.0_r8)) &
         call LS_fatal ('RHO_BRACKET: unable to bracket plane constant Rho')

  end subroutine rho_bracket

  subroutine rho_brent (this, iter, Rho_Min, Rho_Max, V_Min, V_Max, trunc_vol)
    !=======================================================================
    ! PURPOSE -
    !   Compute the plane parameter Rho in the plane equation:
    !                    X*Normal - Rho = 0
    !   iteratively using Brents method as documented in Numerical
    !   Recipes (Press, Flannery, Teukolsky and Vetterling; Cambridge
    !   University Press, 1986), p. 253
    !=======================================================================
    use truncate_volume_module, only: truncate_volume, truncvol_data

    ! Arguments
    class(locate_plane_hex), intent(inout)  :: this
    integer                , intent(out  )  :: iter
    real(r8)           , intent(in)  :: Rho_Min, Rho_Max, V_Min, V_Max
    type(truncvol_data), intent(in) :: trunc_vol(:)

    ! Local Variables
    integer, parameter :: volume_track_iter_max = 10
    real(r8), parameter :: volume_track_iter_tol = 1.0d-8
    real(r8), parameter :: alittle = epsilon(1.0_r8)

    integer :: i
    logical  :: Quad_intrp
    real(r8) :: Rho_a, Rho_b, Rho_c, Rho_d,   &
         Rho_e, Rho_mid, Rho_tol, V_a, &
         V_b, V_c, P, Q, R, S

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize values before the iteration loop
    i = 0
    iter = 0
    Rho_a = Rho_Min
    Rho_b = Rho_Max
    Rho_c = 0.0_r8
    Rho_d = 0.0_r8
    Rho_e = 0.0_r8
    P = 0.0_r8
    Q = 0.0_r8
    R = 0.0_r8
    S = 0.0_r8
    V_a = V_Min - this%Vof*this%Volume
    V_b = V_Max - this%Vof*this%Volume
    V_c = V_b
    Quad_intrp = .false.

    ! Main iteration loop
    do while (abs(V_b/max(this%Volume,alittle)) > volume_track_iter_tol &
         .and. i < volume_track_iter_Max)

      i = i + 1

      if (V_b*V_c > 0.0_r8) then
        Rho_c = Rho_a
        V_c = V_a
        Rho_d = Rho_b - Rho_a
        Rho_e = Rho_d
      end if

      if (ABS(V_c) < ABS(V_b)) then
        Rho_a = Rho_b
        Rho_b = Rho_c
        Rho_c = Rho_a
        V_a = V_b
        V_b = V_c
        V_c = V_a
      end if

      Rho_tol = 2.0_r8*alittle*ABS(Rho_b) + 0.5_r8*volume_track_iter_tol
      Rho_mid = 0.5_r8*(Rho_c - Rho_b)

      if (ABS(Rho_mid) > Rho_tol .and. V_b /= 0.0_r8) Iter = i

      Quad_intrp = ABS(Rho_e) >= Rho_tol .and. ABS(V_a) > ABS(V_b)

      if (Quad_intrp) S = V_b/V_a

      if (Quad_intrp .and. Rho_a == Rho_c) then
        P = 2.0_r8*Rho_mid*S
        Q = 1.0_r8 - S
      end if

      if (Quad_intrp .and. Rho_a /= Rho_c) then
        Q = V_a/V_c
        R = V_b/V_c
        P = S*(2.0_r8*Rho_mid*Q*(Q - R) - (Rho_b - Rho_a)*(R - 1.0_r8))
        Q = (Q - 1.0_r8)*(R - 1.0_r8)*(S - 1.0_r8)
      end if

      if (Quad_intrp .and. P > 0.0_r8) Q = -Q
      if (Quad_intrp) P = ABS(P)

      if (Quad_intrp) R = MIN(3.0_r8*Rho_mid*Q - ABS(Rho_tol*Q), ABS(Rho_e*Q))

      if (Quad_intrp .and. 2.0_r8*P < R) then
        Rho_e = Rho_d
        Rho_d = P/Q
      end if

      if (Quad_intrp .and. 2.0_r8*P >= R) then
        Rho_d = Rho_mid
        Rho_e = Rho_d
      end if

      if (.not. Quad_intrp) then
        Rho_d = Rho_mid
        Rho_e = Rho_d
      end if

      Rho_a = Rho_b ! Move last best guess to A
      V_a  = V_b

      if (ABS(Rho_d) > Rho_tol) then ! Evaluate new trial root
        Rho_b = Rho_b + Rho_d
      else
        Rho_b = Rho_b + SIGN(Rho_tol,Rho_mid)
      end if
      this%Rho = Rho_b

      V_b = truncate_volume (this, trunc_vol) - this%Vof*this%Volume
    end do

    this%Rho = Rho_b
  end subroutine rho_brent

end module locate_plane_module
