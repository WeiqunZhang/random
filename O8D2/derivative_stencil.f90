module derivative_stencil_module

  use bl_types

  implicit none

  public

  integer, parameter :: narrow=1, wide=2, S3D=2, chain=3
  integer, save :: stencil, stencil_ng

  integer, save :: D2_order

  ! for 8th-order first derivatives
  real(dp_t),dimension(4),parameter :: D8 = (/ 0.8_dp_t, -0.2_dp_t, 4._dp_t/105._dp_t, -1._dp_t/280._dp_t /)

  ! 6th-order first derivative
  real(dp_t),dimension(3),parameter :: D6 = (/ 0.75_dp_t, -0.15_dp_t, 1._dp_t/60._dp_t /)

  ! coefficients for 8th-order stencil of second derivatives
  ! d(a*du/dx)/dx = H_{i+1/2} - H_{i-1/2},
  ! where H = a.M.u
  !
  ! optimized for more zeros
!  real(dp_t), private, save :: M8_47 = 683._dp_t/10080._dp_t, M8_48 = -1._dp_t/224._dp_t
  ! optimized for leading order truncation error assuming equal weight for the error terms
  real(dp_t), private, save :: M8_47 = 3557._dp_t/44100._dp_t, M8_48 = -2083._dp_t/117600._dp_t
  !
  real(dp_t), save, dimension(8,8) :: M8

  ! coefficients for 6th-order stencil of second derivatives
  !
  ! optimized for more zeros
!  real(dp_t), private, save :: M6_36 = 1._dp_t/90._dp_t
  ! optimized for leading order truncation error assuming equal weight for the error terms
  real(dp_t), private, save :: M6_36 = 281._dp_t/3600._dp_t
  !
  real(dp_t), save, dimension(6,6) :: M6

contains
  
  subroutine stencil_init(order, O8_m47, O8_m48, O6_m36, stencil_type)
    use bl_error_module
    implicit none
    integer , intent(in) :: order
    real(dp_t), intent(in) :: O8_m47, O8_m48, O6_m36
    character (len=*), intent(in) :: stencil_type

    D2_order = order
    stencil_ng = order/2

    if (trim(stencil_type) == "narrow") then
       stencil = narrow
    else if (trim(stencil_type) == "S3D" .or. trim(stencil_type) == "wide") then
       stencil = wide
    else if (trim(stencil_type) == "chain") then
       stencil = chain
    else
       call bl_error("unknow stencil_type")
    end if

    M8_47 = O8_m47
    M8_48 = O8_m48
    M6_36 = O6_m36

    ! 8th-order
    M8(1,1) = 5._dp_t/336._dp_t + M8_48
    M8(2,1) = -11._dp_t/560._dp_t - 2._dp_t*M8_48
    M8(3,1) = -1._dp_t/280._dp_t
    M8(4,1) = 17._dp_t/1680._dp_t + 2._dp_t*M8_48
    M8(5,1) = -M8_48
    M8(6,1) = 0._dp_t
    M8(7,1) = 0._dp_t
    M8(8,1) = 0._dp_t

    M8(1,2) = -83._dp_t/3600._dp_t - M8_47/5._dp_t - 14._dp_t*M8_48/5._dp_t
    M8(2,2) = -31._dp_t/360._dp_t + M8_47 + 3._dp_t*M8_48
    M8(3,2) = 1097._dp_t/5040._dp_t - 2._dp_t*M8_47 + 6._dp_t*M8_48
    M8(4,2) = -319._dp_t/2520._dp_t + 2._dp_t*M8_47 - 8._dp_t*M8_48
    M8(5,2) = -M8_47
    M8(6,2) = -139._dp_t/25200._dp_t + M8_47/5._dp_t + 9._dp_t*M8_48/5._dp_t
    M8(7,2) = 0._dp_t
    M8(8,2) = 0._dp_t

    M8(1,3) = 299._dp_t/50400._dp_t + 2._dp_t*M8_47/5._dp_t + 13._dp_t*M8_48/5._dp_t
    M8(2,3) = 41._dp_t/200._dp_t - 9._dp_t*M8_47/5._dp_t + 4._dp_t*M8_48/5._dp_t
    M8(3,3) = -1349._dp_t/10080._dp_t + 3._dp_t*M8_47 - 12._dp_t*M8_48
    M8(4,3) = -919._dp_t/5040._dp_t - 2._dp_t*M8_47 + 6._dp_t*M8_48
    M8(5,3) = 65._dp_t/224._dp_t + 7._dp_t*M8_48
    M8(6,3) = -467._dp_t/25200._dp_t + 3._dp_t*M8_47/5._dp_t - 18._dp_t*M8_48/5._dp_t
    M8(7,3) = 503._dp_t/50400._dp_t - M8_47/5._dp_t - 4._dp_t*M8_48/5._dp_t
    M8(8,3) = 0._dp_t

    M8(1,4) = 17._dp_t/12600._dp_t - M8_47/5._dp_t - 4._dp_t*M8_48/5._dp_t
    M8(2,4) = -5927._dp_t/50400._dp_t + 4._dp_t*M8_47/5._dp_t - 9._dp_t*M8_48/5._dp_t
    M8(3,4) = -887._dp_t/5040._dp_t - M8_47 + 6._dp_t*M8_48
    M8(4,4) = -445._dp_t/2016._dp_t
    M8(5,4) = -583._dp_t/720._dp_t + M8_47 - 6._dp_t*M8_48
    M8(6,4) = -3613._dp_t/50400._dp_t - 4._dp_t*M8_47/5._dp_t + 9._dp_t*M8_48/5._dp_t
    M8(7,4) = -17._dp_t/600._dp_t + M8_47/5._dp_t + 4._dp_t*M8_48/5._dp_t
    M8(8,4) = -1._dp_t/1120._dp_t

    M8(1,5) = -M8(8,4)
    M8(2,5) = -M8(7,4)
    M8(3,5) = -M8(6,4)
    M8(4,5) = -M8(5,4)
    M8(5,5) = -M8(4,4)
    M8(6,5) = -M8(3,4)
    M8(7,5) = -M8(2,4)
    M8(8,5) = -M8(1,4)

    M8(1,6) = 0._dp_t
    M8(2,6) = -M8(7,3)
    M8(3,6) = -M8(6,3)
    M8(4,6) = -M8(5,3)
    M8(5,6) = -M8(4,3)
    M8(6,6) = -M8(3,3)
    M8(7,6) = -M8(2,3)
    M8(8,6) = -M8(1,3)

    M8(1,7) = 0._dp_t
    M8(2,7) = 0._dp_t
    M8(3,7) = -M8(6,2)
    M8(4,7) = -M8(5,2)
    M8(5,7) = -M8(4,2)
    M8(6,7) = -M8(3,2)
    M8(7,7) = -M8(2,2)
    M8(8,7) = -M8(1,2)

    M8(1,8) = 0._dp_t
    M8(2,8) = 0._dp_t
    M8(3,8) = 0._dp_t
    M8(4,8) = -M8(5,1)
    M8(5,8) = -M8(4,1)
    M8(6,8) = -M8(3,1)
    M8(7,8) = -M8(2,1)
    M8(8,8) = -M8(1,1)

    ! 6th-order
    M6(1,1) = -11._dp_t/180._dp_t + M6_36
    M6(2,1) = 7._dp_t/60._dp_t - 3._dp_t * M6_36
    M6(3,1) = -1._dp_t/15._dp_t + 3._dp_t * M6_36
    M6(4,1) = -M6_36
    M6(5,1) = 0._dp_t
    M6(6,1) = 0._dp_t
    
    M6(1,2) = 1._dp_t/9._dp_t - 2._dp_t * M6_36
    M6(2,2) = -1._dp_t/120._dp_t + 5._dp_t * M6_36
    M6(3,2) = -11._dp_t/60._dp_t - 3._dp_t * M6_36
    M6(4,2) = 83._dp_t/360._dp_t - M6_36
    M6(5,2) = -1._dp_t/90._dp_t + M6_36
    M6(6,2) = 0._dp_t
    
    M6(1,3) = -1._dp_t/18._dp_t + M6_36
    M6(2,3) = -17._dp_t/90._dp_t - 2._dp_t * M6_36
    M6(3,3) = -101._dp_t/360._dp_t 
    M6(4,3) = -137._dp_t/180._dp_t + 2._dp_t * M6_36
    M6(5,3) = -5._dp_t/72._dp_t - M6_36
    M6(6,3) = -1._dp_t/180._dp_t
    
    M6(1,4) = -M6(6,3)
    M6(2,4) = -M6(5,3)
    M6(3,4) = -M6(4,3)
    M6(4,4) = -M6(3,3)
    M6(5,4) = -M6(2,3)
    M6(6,4) = -M6(1,3)
    
    M6(1,5) = 0._dp_t
    M6(2,5) = -M6(5,2)
    M6(3,5) = -M6(4,2)
    M6(4,5) = -M6(3,2)
    M6(5,5) = -M6(2,2)
    M6(6,5) = -M6(1,2)
    
    M6(1,6) = 0._dp_t
    M6(2,6) = 0._dp_t
    M6(3,6) = M6_36
    M6(4,6) = -M6(3,1)
    M6(5,6) = -M6(2,1)
    M6(6,6) = -M6(1,1)

  end subroutine stencil_init

  function first_deriv_8(u) result(du)
    real(dp_t) :: du
    real(dp_t), intent(in) :: u(-4:4)
    du =   D8(1)*(u(1)-u(-1)) &
         + D8(2)*(u(2)-u(-2)) &
         + D8(3)*(u(3)-u(-3)) &
         + D8(4)*(u(4)-u(-4))
  end function first_deriv_8

  function first_deriv_6(u) result(du)
    real(dp_t) :: du
    real(dp_t), intent(in) :: u(-3:3)
    du =   D6(1)*(u(1)-u(-1)) &
         + D6(2)*(u(2)-u(-2)) &
         + D6(3)*(u(3)-u(-3))
  end function first_deriv_6

end module derivative_stencil_module
