module prob_module

  use multifab_module
  use derivative_stencil_module, only : D2_order

  implicit none

  private

  integer, save :: prob_type
  real(dp_t), save :: diffcoef_max

  public :: init_prob, compute_diffcoef, compute_source, compute_norms, diffcoef_max

contains
  
  subroutine init_prob(u,dx,prob_lo,prob_type_in)

    type(multifab) , intent(inout) :: u
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: prob_lo(u%dim)
    integer        , intent(in   ) :: prob_type_in

    ! local variables
    integer :: lo(u%dim), hi(u%dim)
    integer :: dm, ng, i

    real(kind=dp_t), pointer :: du(:,:,:,:)

    prob_type = prob_type_in

    dm = u%dim
    ng = u%ng

    do i=1,nfabs(u)
       du => dataptr(u,i)
       lo = lwb(get_box(u,i))
       hi = upb(get_box(u,i))
       select case(dm)
       case (2)
          call init_prob_2d(du(:,:,1,1), ng, lo, hi, prob_lo, dx)
       case default
          call bl_error('dim must be 2')
       end select
    end do

  end subroutine init_prob

  subroutine init_prob_2d(u, ng, lo, hi, prob_lo, dx)

    integer          :: lo(2), hi(2), ng
    real(dp_t) :: u(lo(1)-ng:,lo(2)-ng:)
    real(dp_t) :: prob_lo(2)
    real(dp_t) :: dx
 
    ! local varables
    integer          :: i,j
    real(dp_t) :: x,y

    !$omp parallel do private(i,j,x,y)
    do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+0.5_dp_t) * dx
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5_dp_t) * dx

            u(i,j) = sin(x) * sin(y)

         end do
      end do
      !$omp end parallel do

    end subroutine init_prob_2d

  subroutine compute_diffcoef(a,b,u,dx,prob_lo)

    type(multifab) , intent(in   ) :: u
    type(multifab) , intent(inout) :: a, b
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: prob_lo(u%dim)

    ! local variables
    integer :: lo(u%dim), hi(u%dim)
    integer :: dm, ng, i

    real(kind=dp_t), pointer :: du(:,:,:,:), da(:,:,:,:), db(:,:,:,:)
    logical, save :: first_call = .true.

    if (first_call) then
       select case(prob_type)
       case (1)
          diffcoef_max = 1.1_dp_t
       case (2)
          diffcoef_max = 1.9_dp_t
       case (3)
          diffcoef_max = 1.9_dp_t
       case default
          call bl_error('Unknown prob_type')
       end select

       first_call = .false.
    else
       if (prob_type .ne. 3) then ! a & b do not change
          return
       end if
    end if

    dm = u%dim
    ng = u%ng

    do i=1,nfabs(u)
       du => dataptr(u,i)
       da => dataptr(a,i)
       db => dataptr(b,i)
       lo = lwb(get_box(u,i))
       hi = upb(get_box(u,i))
       select case(dm)
       case (2)
          call comp_dc_2d(du(:,:,1,1), da(:,:,1,1), db(:,:,1,1), &
               ng, lo, hi, prob_lo, dx)
       case default
          call bl_error('dim must be 2')
       end select
    end do

  end subroutine compute_diffcoef

  subroutine comp_dc_2d(u, a, b, ng, lo, hi, prob_lo, dx)

    integer          :: lo(2), hi(2), ng
    real(dp_t) :: u(lo(1)-ng:,lo(2)-ng:), a(lo(1)-ng:,lo(2)-ng:), b(lo(1)-ng:,lo(2)-ng:)
    real(dp_t) :: prob_lo(2)
    real(dp_t) :: dx
 
    ! local varables
    integer          :: i,j
    real(dp_t) :: x,y,dudx,dudy
    real(dp_t) :: pi = 4._dp_t*atan(1._dp_t)

    !$omp parallel do private(i,j,x,y,dudx,dudy)
    do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+0.5_dp_t) * dx
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5_dp_t) * dx

            if (prob_type .eq. 1) then
               a(i,j) = 1._dp_t + 0.1_dp_t*cos(x)*cos(y)
               b(i,j) = a(i,j)
            else if (prob_type .eq. 2) then
               a(i,j) = 1._dp_t + 0.9_dp_t * cos(2._dp_t*x) * sin(2._dp_t*y + pi/3._dp_t)
               b(i,j) = 1._dp_t + 0.9_dp_t * cos(2._dp_t*x + pi/3._dp_t) * sin(2._dp_t*y)
            else if (prob_type .eq. 3) then
               if (D2_order .eq. 8) then
                  dudx = (       0.8_dp_t*(u(i+1,j)-u(i-1,j)) &
                       &       - 0.2_dp_t*(u(i+2,j)-u(i-2,j)) &
                       + (4._dp_t/105._dp_t)*(u(i+3,j)-u(i-3,j)) &
                       - (1._dp_t/280._dp_t)*(u(i+4,j)-u(i-4,j)) ) / dx
                  dudy = (       0.8_dp_t*(u(i,j+1)-u(i,j-1)) &
                       &       - 0.2_dp_t*(u(i,j+2)-u(i,j-2)) &
                       + (4._dp_t/105._dp_t)*(u(i,j+3)-u(i,j-3)) &
                       - (1._dp_t/280._dp_t)*(u(i,j+4)-u(i,j-4)) ) / dx
               else
                  call bl_error('Only 8th-order is supported')
               end if
               a(i,j) = 1._dp_t + 0.9_dp_t*cos(10._dp_t*dudx)*sin(10._dp_t*y+pi/3._dp_t)
               b(i,j) = 1._dp_t + 0.9_dp_t*cos(10._dp_t*dudy)*sin(10._dp_t*x+pi/3._dp_t)
            else
               call bl_error('Unknow prob_type')
            end if

         end do
      end do
      !$omp end parallel do

    end subroutine comp_dc_2d

    subroutine compute_source(g,time,dx,prob_lo)
      type(multifab) , intent(inout) :: g
      real(kind=dp_t), intent(in   ) :: time, dx
      real(kind=dp_t), intent(in   ) :: prob_lo(g%dim)

      ! local variables
      integer :: lo(g%dim), hi(g%dim)
      integer :: dm, ng, i
      real(kind=dp_t), pointer :: dg(:,:,:,:)

      if (prob_type .ne. 1) return

      dm = g%dim
      ng = g%ng
      
      do i=1,nfabs(g)
         dg => dataptr(g,i)
         lo = lwb(get_box(g,i))
         hi = upb(get_box(g,i))
         select case(dm)
         case (2)
            call comp_src_2d(dg(:,:,1,1), time, ng, lo, hi, prob_lo, dx)
         case default
            call bl_error('dim must be 2')
         end select
      end do

    end subroutine compute_source

    subroutine comp_src_2d(g, time, ng, lo, hi, prob_lo, dx)

      integer          :: lo(2), hi(2), ng
      real(dp_t) :: g(lo(1)-ng:,lo(2)-ng:)
      real(dp_t) :: prob_lo(2)
      real(dp_t) :: dx, time
 
      ! local varables
      integer          :: i,j
      real(dp_t) :: x,y,uex

      !$omp parallel do private(i,j,x,y,uex)
      do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+0.5_dp_t) * dx
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+0.5_dp_t) * dx

            uex = exp(-time) * sin(x) * sin(y)
            g(i,j) = (1._dp_t+4._dp_t*0.1_dp_t*cos(x)*cos(y))*uex

         end do
      end do
      !$omp end parallel do

    end subroutine comp_src_2d

    subroutine compute_norms(u, dx, time, prob_lo, l0, l2)
      type(multifab), intent(in) :: u
      real(kind=dp_t), intent(in) :: dx, time
      real(kind=dp_t), intent(in) :: prob_lo(u%dim)
      real(kind=dp_t), intent(out) :: l0, l2

      type(multifab) :: error

      ! local variables
      integer :: lo(u%dim), hi(u%dim)
      integer :: dm, ng, i
      real(kind=dp_t) :: dvol
      real(kind=dp_t), pointer :: derr(:,:,:,:)

      dm = u%dim
      ng = u%ng

      call multifab_build(error,u%la,1,0)

      do i=1,nfabs(error)
         derr => dataptr(error,i)
         lo = lwb(get_box(error,i))
         hi = upb(get_box(error,i))
         select case(dm)
         case (2)
            call comp_ex_2d(derr(:,:,1,1), time, lo, hi, prob_lo, dx)
         case default
            call bl_error('dim must be 2')
         end select
      end do

      call multifab_sub_sub(error, u)

      l0 = multifab_norm_inf(error)

      l2 = multifab_norm_l2(error)
      dvol = dx**dm
      l2 = l2 * sqrt(dvol)

      call multifab_destroy(error)

    end subroutine compute_norms

    subroutine comp_ex_2d(ex, time, lo, hi, prob_lo, dx)
      integer          :: lo(2), hi(2)
      real(dp_t) :: ex(lo(1):,lo(2):)
      real(dp_t) :: prob_lo(2)
      real(dp_t) :: dx, time
 
      ! local varables
      integer          :: i,j
      real(dp_t) :: x,y

      if (prob_type .eq. 1) then
         !$omp parallel do private(i,j,x,y)
         do j=lo(2),hi(2)
            y = prob_lo(2) + (dble(j)+0.5_dp_t) * dx
            do i=lo(1),hi(1)
               x = prob_lo(1) + (dble(i)+0.5_dp_t) * dx
               ex(i,j) = exp(-time)*sin(x)*sin(y)
            end do
         end do
         !$omp end parallel do
      else
         ex = 0._dp_t
      end if
    end subroutine comp_ex_2d

end module prob_module
