module truncate_volume_module
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures necessary to compute hexahedral volumes
  !   truncated by a planar interface
  !
  !   Public Interface:
  !
  !     * call FACE_PARAM (option, face, K, Lambda, MUa, MUi, MUp, Nu, V1234, X)
  !
  !         Compute and store various face parameters needed for a volume
  !         truncation calculation.
  !
  !     * call TRUNCATE_VOLUME (K, Lambda, MUa, MUi, MUp, Nu, Vf, Vol, V1234, X)
  !
  !         Compute the volume truncated by a plane.
  !
  !     * call TRUNCATE_FACE (f, K, Lambda, MUa, MUi, MUp, Nu, Vf, V1234, X)
  !
  !         Compute the volume truncated at the current hex face by a plane.
  !
  ! Contains: FACE_PARAM
  !           TRUNCATE_VOLUME
  !           TRUNCATE_FACE
  !           TRUNCATE_FACE_2
  !           TRUNCATE_FACE_2
  !           TRUNCATE_FACE_4
  !           TRUNCATE_FACE_N
  !           Y_FUNCTION
  !
  ! Author(s): Douglas B. Kothe, LANL (dbk@lanl.gov)
  !            S. Jay Mosso, LANL (sjm@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  implicit none
  private

  public :: truncate_volume, truncate_face, face_param, truncvol_data

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  integer, parameter :: nvf = 4 ! number of vertices per cell face
  integer, parameter :: nfc = 6 ! number of faces
  integer, parameter :: alittle = epsilon(1.0_r8)
  real(r8), parameter :: eps(4) = (/ 1.0_r8, -1.0_r8, 1.0_r8, -1.0_r8 /)

  type truncvol_data
    real(r8) :: K(3)
    real(r8) :: Lambda
    real(r8) :: MUa(nvf)
    real(r8) :: MUi(nvf)
    real(r8) :: Nu
    real(r8) :: V1234
    real(r8) :: X(nvf,3)
    integer  :: MUp(nvf)
  end type truncvol_data

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

contains

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  subroutine face_param (cell, option, face, trunc_vol, int_flux_vol_coord)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute and store various face parameters needed
    !   for the volume truncation calculation
    !
    !=======================================================================
    use hex_types, only: reconstruction_hex
    use logging_services

    ! Arguments
    class(reconstruction_hex), intent(in) :: cell
    character(9), intent(IN) :: option
    integer, intent(IN) :: face
    type(truncvol_data), dimension(nfc), intent(inout) :: trunc_vol
    real(r8), dimension(3,8), intent(in), optional :: int_flux_vol_coord

    ! Local Variables
    integer :: i, j, n, v1, v2, v3, v4
    real(r8), dimension(3) :: Tmp1, Tmp2

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Store face vertices in an order that is
    !    counterclockwise relative to a direction
    !    that looks from outside the face into
    !    the cell.
    ! For option = 'full_cell', use the donor cell vertices (Cell_Coord)
    ! For option = 'flux_cell', use the vertices of the flux volume (Xv_flux)
    select case (face)
    case(1)
      v1 = 4; v2 = 8; v3 = 7; v4 = 3      ! Left face
    case(2)
      v1 = 5; v2 = 1; v3 = 2; v4 = 6      ! Right face
    case(3)
      v1 = 5; v2 = 8; v3 = 4; v4 = 1      ! Front face
    case(4)
      v1 = 6; v2 = 2; v3 = 3; v4 = 7      ! Back face
    case(5)
      v1 = 3; v2 = 2; v3 = 1; v4 = 4      ! Bottom face
    case(6)
      v1 = 7; v2 = 8; v3 = 5; v4 = 6      ! Top face
    end select

    if (option == 'full_cell') then
      Trunc_Vol(face)%X(1,:) = cell%node(:,v1)
      Trunc_Vol(face)%X(2,:) = cell%node(:,v2)
      Trunc_Vol(face)%X(3,:) = cell%node(:,v3)
      Trunc_Vol(face)%X(4,:) = cell%node(:,v4)
    else if (option == 'flux_cell') then
      if (.not.present(int_flux_vol_coord)) call LS_fatal('when calculating flux_cell face_param, int_flux_vol_coord is required')
      Trunc_Vol(face)%X(1,:) = Int_Flux_Vol_Coord(:,v1)
      Trunc_Vol(face)%X(2,:) = Int_Flux_Vol_Coord(:,v2)
      Trunc_Vol(face)%X(3,:) = Int_Flux_Vol_Coord(:,v3)
      Trunc_Vol(face)%X(4,:) = Int_Flux_Vol_Coord(:,v4)
    endif

    ! Compute K   (the Area Vector of the ruled surface)
    do i = 1,3
      Tmp1(i) = Trunc_Vol(face)%X(3,i) - Trunc_Vol(face)%X(1,i)
      Tmp2(i) = Trunc_Vol(face)%X(4,i) - Trunc_Vol(face)%X(2,i)
    end do

    Trunc_Vol(face)%K(1) = Tmp1(2)*Tmp2(3) - Tmp1(3)*Tmp2(2)
    Trunc_Vol(face)%K(2) = Tmp1(3)*Tmp2(1) - Tmp1(1)*Tmp2(3)
    Trunc_Vol(face)%K(3) = Tmp1(1)*Tmp2(2) - Tmp1(2)*Tmp2(1)

    ! Compute V1234
    Trunc_Vol(face)%V1234 = 0.0_r8
    do i = 1,3
      Tmp1(i) = Trunc_Vol(face)%X(1,i) - Trunc_Vol(face)%X(2,i) + &
           Trunc_Vol(face)%X(3,i) - Trunc_Vol(face)%X(4,i)
      Trunc_Vol(face)%V1234 = Trunc_Vol(face)%V1234 +       &
           Tmp1(i)*Trunc_Vol(face)%K(i)
    end do
    Trunc_Vol(face)%V1234 = 0.5_r8 * Trunc_Vol(face)%V1234

    ! Compute the Mu-i-s.  This is the normal of the interface dotted
    ! with the coordinates of each faces vertex.  Mu-p is the vertex number
    ! of the vertex for each face (1 <= Mu-p <= nvf).  When the Mu-i-s are
    ! in ascending order, it denotes which vertex the interface will
    ! pass through first, second, etc.  The variable Mu-p is the vertex
    ! number of the reordered distances.
    do j = 1, nvf
      Trunc_Vol(face)%MUi(j) = 0.0_r8
      do i = 1,3
        Trunc_Vol(face)%MUi(j) = Trunc_Vol(face)%MUi(j) +     &
             cell%Normal(i)*Trunc_Vol(face)%X(j,i)
      end do
      Trunc_Vol(face)%MUp(j) = j
      Trunc_Vol(face)%MUa(j) = Trunc_Vol(face)%MUi(j)
    end do

    ! Here Nu and Lambda are temporaries used to facilitate ordering the
    ! Mu-i into the Mu-a.  Put the minimum distance (the first vertex
    ! that the interface will pass through) into Mu-p(1) and put its
    ! facial vertex number into Mu-p(1).
    Trunc_Vol(face)%Nu = 0.0_r8
    Trunc_Vol(face)%Lambda = MIN(Trunc_Vol(face)%MUi(1),Trunc_Vol(face)%MUi(2), &
         Trunc_Vol(face)%MUi(3),Trunc_Vol(face)%MUi(4))
    do j = 1, nvf
      if (Trunc_Vol(face)%MUi(j) == Trunc_Vol(face)%Lambda .and.   &
           Trunc_Vol(face)%Nu == 0.0_r8) then
        Trunc_Vol(face)%MUp(j) = 1
        Trunc_Vol(face)%MUp(1) = j
        Trunc_Vol(face)%MUa(j) = Trunc_Vol(face)%MUa(1)
        Trunc_Vol(face)%MUa(1) = Trunc_Vol(face)%MUi(j)
        Trunc_Vol(face)%Nu = 1.0_r8
      end if
    end do

    ! Now that the minimum distance is in element 1, order the
    ! other vertices in ascending order using a bubble sort.
    if (Trunc_Vol(face)%MUa(3) > Trunc_Vol(face)%MUa(4)) then
      Trunc_Vol(face)%Lambda   = Trunc_Vol(face)%MUp(4)
      Trunc_Vol(face)%MUp   (4) = Trunc_Vol(face)%MUp(3)
      Trunc_Vol(face)%MUp   (3) = Trunc_Vol(face)%Lambda

      Trunc_Vol(face)%Lambda   = Trunc_Vol(face)%MUa(4)
      Trunc_Vol(face)%MUa   (4) = Trunc_Vol(face)%MUa(3)
      Trunc_Vol(face)%MUa   (3) = Trunc_Vol(face)%Lambda
    end if

    if (Trunc_Vol(face)%MUa(2) > Trunc_Vol(face)%MUa(3)) then
      Trunc_Vol(face)%Lambda   = Trunc_Vol(face)%MUp(3)
      Trunc_Vol(face)%MUp   (3) = Trunc_Vol(face)%MUp(2)
      Trunc_Vol(face)%MUp   (2) = Trunc_Vol(face)%Lambda

      Trunc_Vol(face)%Lambda   = Trunc_Vol(face)%MUa(3)
      Trunc_Vol(face)%MUa   (3) = Trunc_Vol(face)%MUa(2)
      Trunc_Vol(face)%MUa   (2) = Trunc_Vol(face)%Lambda
    end if

    if (Trunc_Vol(face)%MUa(3) > Trunc_Vol(face)%MUa(4)) then
      Trunc_Vol(face)%Lambda   = Trunc_Vol(face)%MUp   (4)
      Trunc_Vol(face)%MUp   (4) = Trunc_Vol(face)%MUp   (3)
      Trunc_Vol(face)%MUp   (3) = Trunc_Vol(face)%Lambda

      Trunc_Vol(face)%Lambda   = Trunc_Vol(face)%MUa   (4)
      Trunc_Vol(face)%MUa   (4) = Trunc_Vol(face)%MUa   (3)
      Trunc_Vol(face)%MUa   (3) = Trunc_Vol(face)%Lambda
    end if

    ! The definition of Lambda is given in Eqn. 11.5 of Zemach-s notes.
    Trunc_Vol(face)%Lambda = Trunc_Vol(face)%MUi(2)*Trunc_Vol(face)%MUi(4) &
         - Trunc_Vol(face)%MUi(1)*Trunc_Vol(face)%MUi(3)

    ! Nu is the face deviation vector dotted with the interface normal.
    ! If a face is a parallelogram, Nu (and B) will be zero.
    ! If a face is not a parallelogram, the magnitude of B measures
    ! the deviation from the parallelogram.
    Trunc_Vol(face)%Nu = Trunc_Vol(face)%MUi(1) - Trunc_Vol(face)%MUi(2) + &
         Trunc_Vol(face)%MUi(3) - Trunc_Vol(face)%MUi(4)

  end subroutine face_param

  real(r8) function truncate_volume (cell, trunc_vol)
    !=======================================================================
    ! Purpose(s):
    !  
    !   Compute the truncation volume.
    !
    !=======================================================================
    !use vof_data_module,  only: Cases, count_cases

    use hex_types, only: reconstruction_hex

    ! Arguments
    class(reconstruction_hex), intent(in) :: cell
    type(truncvol_data), dimension(nfc), intent(in) :: trunc_vol

    ! Local Variables
    integer :: f
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant quantities
    !if (count_cases) Cases = 0

    ! Loop over faces, accumulating the truncated volume
    Truncate_Volume     = 0.0_r8
    do f = 1, nfc
      Truncate_Volume = Truncate_Volume + truncate_face (cell, trunc_vol(f))
    end do

  end function truncate_volume

  !=======================================================================
  ! Purpose(s):
  !
  !   Compute the volume truncated at the current
  !   hex face by the plane given by X*Normal - Ro = 0
  !
  !=======================================================================
  function truncate_face (cell, trunc_vol_face) result(Vf)
    use hex_types, only: reconstruction_hex

    class(reconstruction_hex), intent(in) :: cell
    type(truncvol_data),      intent(in) :: trunc_vol_face
    real(r8)                             :: Vf

    real(r8) :: Y(nvf)

    Y = Y_function (cell, trunc_vol_face)

    ! get intersection case
    if ( trunc_vol_face%MUa(1) < cell%rho .and. cell%Rho <= trunc_vol_face%MUa(2) .and. &
         trunc_vol_face%MUa(1) /= trunc_vol_face%MUa(2)) then                                ! case 1
      Vf = truncate_face_N (1, cell%rho, Y, trunc_vol_face)
    else if ( trunc_vol_face%MUa(2) < cell%Rho .and. cell%Rho <= trunc_vol_face%MUa(3)) then ! case 2 & 5
      if (trunc_vol_face%MUp(2) == (1 + mod(trunc_vol_face%MUp(1)+1, nvf))) then ! case 5
        Vf = truncate_face_N (1, cell%rho, Y, trunc_vol_face) + &
             truncate_face_N (2, cell%rho, Y, trunc_vol_face)
      else                                                                       ! case 2
        Vf = truncate_face_2 (cell%rho, Y, trunc_vol_face)
      end if
    else if ( cell%Rho > trunc_vol_face%MUa(3)) then                                         ! case 3 & 4
      if (cell%Rho < trunc_vol_face%MUa(4)) &          ! case 3
           Vf = - truncate_face_N (4, cell%rho, Y, trunc_vol_face)
      
      Vf = Vf + truncate_face_4 (cell, trunc_vol_face) ! case 4
    end if
  end function truncate_face

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><>

  function truncate_face_2 (cell_rho, Y, trunc_vol_face) result(Q)
    !=======================================================================
    ! PURPOSE - 
    !   Compute an expression (for case "2") that is part of
    !   the volume truncated along a hex face by the plane
    !   X*Normal - Ro = 0
    !=======================================================================

    ! Arguments
    real(r8), intent(in) :: cell_rho
    real(r8), dimension(nvf), intent(IN)  :: Y
    type(truncvol_data), intent(in) :: trunc_vol_face
    real(r8) :: Q 

    ! Local Variables
    integer :: i, k
    real(r8) :: J1, J2, J3, W1, W2, W3, W4, Z1, Z2, Z3

    real(r8), parameter :: one_third     = 1.0_r8 / 3.0_r8
    real(r8), parameter :: one_fourth    = 0.25_r8
    real(r8), parameter :: one_fifth     = 0.2_r8
    real(r8), parameter :: one_sixth     = 1.0_r8 / 6.0_r8
    real(r8), parameter :: one_seventh   = 1.0_r8 / 7.0_r8
    real(r8), parameter :: one_eighth    = 0.125_r8

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Q  = 0.0_r8
    J1 = 0.0_r8
    J2 = 0.0_r8
    J3 = 0.0_r8
    W1 = 0.0_r8
    W2 = 0.0_r8
    W3 = 0.0_r8
    W4 = 0.0_r8
    Z1 = 0.0_r8
    Z2 = 0.0_r8
    Z3 = 0.0_r8

    do k = 1, nvf
      i = 1 + mod(k + 1,nvf)

      ! Scan the four face indicies until we get to the 
      ! 2nd predecessor of MU(b).
      if (trunc_vol_face%MUp(2) == k) then

        ! Z3 is the Y function of the successor of the A vertex
        Z3 = Y(k)

        ! W1 is the difference between: 
        !   the Mu of the 2nd successor of vertex B and the Mu of vertex A.
        W1 = (trunc_vol_face%MUi(i) - trunc_vol_face%MUa(1))
      end if
    end do

    if (ABS(W1) > alittle) then
      W1 = 1.0_r8 / W1
    else
      W1 = 0.0_r8
    end if

    do k = 1, nvf
      i = 1 + mod(k + 1,nvf)

      ! Scan the four face indicies until we get to the 
      ! 2nd predecessor of MU(a).
      if (trunc_vol_face%MUp(1) == k) then

        ! Z1 is the Y function of the A vertex
        Z1 = Y(k)

        ! Z2 is the Y function of the 2nd successor of the A vertex
        Z2 = Y(i) 
        W3 = eps(k)*trunc_vol_face%Nu*W1
        Q  = eps(k)*trunc_vol_face%V1234*0.5_r8*W1*W1 

        ! W2 is the difference between: 
        !   the Mu of the 2nd successor of vertex A and 
        !   the Mu of vertex B.
        W2 = trunc_vol_face%MUi(i) - trunc_vol_face%MUa(2)
      end if
    end do

    if (abs(W2) > alittle) then
      W2 = 1.0_r8 / W2
    else
      W2 = 0.0_r8
    end if

    if ((abs(W3) > 1.0e-2) .and. (W3 /= -1.0_r8)) then
      W4 = 1.0_r8/W3
      J1 = (1.0_r8 - log(abs(1.0_r8 + W3))*W4)*W4
      J2 = (0.5_r8 - J1)*W4
      J3 = (one_third - J2)*W4
    else
      J3 = one_fourth - one_fifth*W3 + one_sixth*W3*W3 - &
           one_seventh*W3*W3*W3 + one_eighth*W3**4
      J2 = one_third - W3*J3
      J1 = 0.5_r8 - W3*J2
    end if

    Q = Q*(J1*(Cell_Rho - trunc_vol_face%MUa(1))**2 -  &
         2.0_r8*(trunc_vol_face%MUa(2) - trunc_vol_face%MUa(1))*&
         (Cell_Rho - trunc_vol_face%MUa(1))*J2 &
         + J3*(trunc_vol_face%MUa(2) - trunc_vol_face%MUa(1))**2)

    Q = Q + W1*Z1*(2.0_r8*Cell_Rho - trunc_vol_face%MUa(1) - &
         trunc_vol_face%MUa(2))/6.0_r8

    Q = Q + (W1*W2*(Z2 - Z3)*(Cell_Rho - trunc_vol_face%MUa(2))**2)/6.0_r8

  end function TRUNCATE_FACE_2

  function TRUNCATE_FACE_4 (cell, trunc_vol_face) result(Q)
    !=======================================================================
    ! PURPOSE - 
    !   Compute an expression (for case "4") that is part of
    !   the volume truncated along a hex face by the plane
    !   X*Normal - Ro = 0
    !=======================================================================
    use hex_types, only: reconstruction_hex

    ! Arguments
    class(reconstruction_hex), intent(in) :: cell
    type(truncvol_data), intent(in) :: trunc_vol_face
    real(r8) :: Q

    ! Local Variables
    real(r8), dimension(3) :: Tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Tmp = 0.25_r8*(trunc_vol_face%X(1,:) + trunc_vol_face%X(2,:)  + &
         &         trunc_vol_face%X(3,:) + trunc_vol_face%X(4,:)) - &
         Cell%Rho*Cell%Normal(:)

    Q = sum( tmp * trunc_vol_face%K ) / 6.0_r8

  end function TRUNCATE_FACE_4

  function TRUNCATE_FACE_N (n, plane_const, Y, trunc_vol_face) result(Q)
    !=======================================================================
    ! PURPOSE - 
    !   Compute an expression (for case "n") that is part of
    !   the volume truncated along a hex face by the plane
    !   X*Normal - Ro = 0
    !=======================================================================

    ! Arguments
    integer,  intent(in) :: n
    real(r8), intent(in) :: plane_const
    real(r8), dimension(nvf), intent(in) :: Y
    type(truncvol_data), intent(in) :: trunc_vol_face
    real(r8)  :: Q

    ! Local Variables
    integer  :: k
    real(r8) :: J1, S, T, W1, W3
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


    T = plane_const - trunc_vol_face%MUa(n)
    W3 = trunc_vol_face%Lambda + trunc_vol_face%Nu*trunc_vol_face%MUa(n)

    if (abs(W3) > alittle) then
      W3 = 1.0_r8/W3
    else
      W3 = 0.0_r8
    end if

    W1 = T*T*W3
    T = trunc_vol_face%Nu * T * W3

    S = 0.0_r8
    if (abs(T) > alittle) S = 1.0_r8/T

    J1 = 0.0_r8
    if (T /= -1.0_r8) J1 = (1.0_r8-log(abs(1.0_r8+T))*S)*S

    if (abs(T) > 1.0e-2) then
      W3 = J1 + (-2.0_r8/3.0_r8 + 2.0_r8*J1 + (J1 - 0.5_r8) * S) * S
    else
      W3 = 1.0_r8/12.0_r8 + T*(-1.0_r8/30.0_r8 + T*(1.0_r8/60.0_r8 + T*(-1.0_r8/105.0_r8 + T*1.0_r8/168.0_r8)))
    end if

    Q = 0.0_r8
    do k = 1, nvf
      if (trunc_vol_face%MUp(n) == k) &
           Q = eps(k)*(Y(k)*W1/6.0_r8 + 0.5_r8*trunc_vol_face%V1234*W3*W1*W1)
    end do

  end function TRUNCATE_FACE_N

  function Y_FUNCTION (cell, trunc_vol_face)
    !=======================================================================
    ! PURPOSE -
    !   THIS ROUTINE COMPUTES THE VARIABLE :
    !   Yi = (Xi'''-Normal.Ro) . (Xi-Normal.Ro) X (Xi'-Normal.Ro)
    !
    ! Note:
    !   Xi     - is the i-th vertex of face "FACE".
    !   Xi'''  - is the predecessor vertex (third successor vertex) to
    !            vertex Xi
    !   Xi'    - is the successor vertex to vertex Xi
    !   Normal - is the normal vector to the material interface
    !   Ro     - is the distance of vertex Xi from the origin along
    !            the interface normal "NORMAL"
    !   X      - the coordinates of the four vertices that define
    !            call face "FACE".
    !   Y      - the value of the Y-FUNCTION for this value of RO.
    !
    ! This routine is called from truncate_face.
    !=======================================================================
    use hex_types, only: reconstruction_hex

    ! Arguments
    class(reconstruction_hex), intent(in) :: cell
    type(truncvol_data), intent(in) :: trunc_vol_face
    real(r8), dimension(nvf) :: Y_FUNCTION

    ! Local Variables
    integer, parameter :: ICMP = 1, JCMP = 2, KCMP = 3
    real(r8), dimension(3) :: R1, R2, R3, R4, S1, S3

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    R4 = Cell%Rho * Cell%Normal

    R1 = trunc_vol_face%X(1,:) - R4
    R2 = trunc_vol_face%X(2,:) - R4
    R3 = trunc_vol_face%X(3,:) - R4
    R4 = trunc_vol_face%X(4,:) - R4

    !  S1 = (X1 - Normal.Ro) X (X2 - Normal.Ro)
    S1(ICMP) = R1(JCMP)*R2(KCMP) - R1(KCMP)*R2(JCMP)
    S1(JCMP) = R1(KCMP)*R2(ICMP) - R1(ICMP)*R2(KCMP)
    S1(KCMP) = R1(ICMP)*R2(JCMP) - R1(JCMP)*R2(ICMP)

    !  S3 = (X3 - Normal.Ro) X (X4 - Normal.Ro)
    S3(ICMP) = R3(JCMP)*R4(KCMP) - R3(KCMP)*R4(JCMP)
    S3(JCMP) = R3(KCMP)*R4(ICMP) - R3(ICMP)*R4(KCMP)
    S3(KCMP) = R3(ICMP)*R4(JCMP) - R3(JCMP)*R4(ICMP)

    !  Y1 = (X4 - Normal.Ro) . (X1 - Normal.Ro) X (X2 - Normal.Ro)
    Y_FUNCTION(1) = R4(ICMP)*S1(ICMP) + R4(JCMP)*S1(JCMP) + R4(KCMP)*S1(KCMP)

    !  Y2 = (X1 - Normal.Ro) . (X2 - Normal.Ro) X (X3 - Normal.Ro)
    !     However, we take advantage of the fact the vector identity:
    !     A . B X C = C . A X B  which means that we don-t need to
    !     compute the S2 vector as above.  We can use the S1 vector:
    !     Y2 = (X3 - Normal.Ro) . (X1 - Normal.Ro) X (X2 - Normal.Ro)
    Y_FUNCTION(2) = R3(ICMP)*S1(ICMP) + R3(JCMP)*S1(JCMP) + R3(KCMP)*S1(KCMP)

    !  Y3 = (X2 - Normal.Ro) . (X3 - Normal.Ro) X (X4 - Normal.Ro)
    Y_FUNCTION(3) = R2(ICMP)*S3(ICMP) + R2(JCMP)*S3(JCMP) + R2(KCMP)*S3(KCMP)

    !  Y4 = (X3 - Normal.Ro) . (X4 - Normal.Ro) X (X1 - Normal.Ro)
    !     However, we take advantage of the fact the vector identity:
    !     A . B X C = C . A X B  which means that we don-t need to
    !     compute the S4 vector as above.  We can use the S3 vector:
    !     Y4 = (X1 - Normal.Ro) . (X3 - Normal.Ro) X (X4 - Normal.Ro)
    Y_FUNCTION(4) = R1(ICMP)*S3(ICMP) + R1(JCMP)*S3(JCMP) + R1(KCMP)*S3(KCMP)

  end function Y_FUNCTION

end module truncate_volume_module
