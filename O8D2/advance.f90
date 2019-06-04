module advance_module

  use multifab_module
  use prob_module
  use derivative_stencil_module
  use sdcquad_module

  implicit none

  private

  integer, save :: verbose
  integer, save :: advance_method, sdc_nnodes, sdc_iters
  real(dp_t), save :: sdc_tol_residual
  type(sdcquad), save :: sdc

  public :: advance, compute_dt, advance_init, advance_close

contains
  
  subroutine advance_init(advance_method_in, verbose_in, &
       sdc_nnodes_in, sdc_iters_in, sdc_tol_residual_in)
    integer, intent(in)  :: advance_method_in, verbose_in, sdc_nnodes_in, sdc_iters_in
    real(dp_t), intent(in) :: sdc_tol_residual_in

    advance_method = advance_method_in
    sdc_nnodes = sdc_nnodes_in
    sdc_iters = sdc_iters_in
    sdc_tol_residual = sdc_tol_residual_in

    verbose = verbose_in

    if (advance_method == 2) then
       call create(sdc, 1, sdc_nnodes)
       sdc%iters        = sdc_iters
       sdc%tol_residual = sdc_tol_residual
       call build(sdc)
    end if
  end subroutine advance_init

  subroutine advance_close()
    if (advance_method == 2) then
       call destroy(sdc)
    end if
  end subroutine advance_close

  subroutine advance(u,a,b,g,dx,time,dt,prob_lo,prob_hi)

    type(multifab) , intent(inout) :: u,g,a,b
    real(kind=dp_t), intent(in   ) :: dx, time, dt
    real(kind=dp_t), intent(in   ) :: prob_lo(u%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(u%dim)

    select case(advance_method)
    case(1)
       call advance_rk4(u,a,b,g,dx,time,dt,prob_lo,prob_hi)
    case(2)
       call advance_sdc(u,a,b,g,dx,time,dt,prob_lo,prob_hi)
    case default
       call bl_error("Invalid advance_method.")
    end select

    if (dp_t .eq. 8) then
       if (contains_nan(U)) then
          call bl_error("U contains nan")
       end if
    end if
       
  end subroutine advance


  subroutine advance_sdc(u,a,b,g,dx,time,dt,prob_lo,prob_hi)

    type(multifab) , intent(inout) :: u,g,a,b
    real(kind=dp_t), intent(in   ) :: dx, time, dt
    real(kind=dp_t), intent(in   ) :: prob_lo(u%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(u%dim)

    integer :: k, m, n, ng
    real(dp_t) :: res_proc, res
    type(layout) :: la
    type(multifab) :: uSDC(sdc%nnodes), fSDC(sdc%nnodes), S(sdc%nnodes-1)

    real(dp_t) :: dtsdc(sdc%nnodes-1), tnodes(sdc%nnodes)

    ng = nghost(U)
    la = get_layout(U)

    ! build u and u' multifabs for each node
    do m = 1, sdc%nnodes
       call build(uSDC(m), la, 1, ng)
       call build(fSDC(m), la, 1, 0)
    end do

    ! build S multifab (node to node integrals)
    do m = 1, sdc%nnodes-1
       call build(S(m), la, 1, 0)
    end do

    ! set provisional solution

    call copy(uSDC(1), U)
    call compute_dudt(uSDC(1), a, b, g, fSDC(1), time, dx, prob_lo)

    do m = 2, sdc%nnodes
       call copy(uSDC(m), uSDC(1))
       call copy(fSDC(m), fSDC(1))
    end do

     ! perform sdc iterations
    res = 0.0_dp_t
    dtsdc = dt * (sdc%nodes(2:sdc%nnodes) - sdc%nodes(1:sdc%nnodes-1))

    tnodes = time + dt * sdc%nodes

    do k = 1, sdc%iters

       ! compute integrals (compact forward euler)
       do m = 1, sdc%nnodes-1
          call setval(S(m), 0.0_dp_t)
          do n = 1, sdc%nnodes
             call saxpy(S(m), sdc%smats(m,n,1), fSDC(n))
          end do
       end do

       ! perform sub-step correction
       do m = 1, sdc%nnodes-1

          ! U(m+1) = U(m) + dt dUdt(m) + dt S(m)

          call copy(uSDC(m+1), uSDC(m))
          call saxpy(uSDC(m+1), dtsdc(m), fSDC(m))
          call saxpy(uSDC(m+1), dt, S(m))

          call compute_dudt(uSDC(m+1), a, b, g, fSDC(m+1), tnodes(m+1), dx, prob_lo)

       end do

       ! check residual
       if (sdc%tol_residual > 0._dp_t) then
          res_proc = sdc_residual(uSDC, fSDC, S(1), dt, sdc)
          call parallel_reduce(res, res_proc, MPI_MAX)

          if (parallel_IOProcessor() .and. verbose > 0) then
             print *, "SDC: iter:", k, "residual:", res
          end if

          if (res < sdc%tol_residual) exit
       end if
    end do

    call copy(U, uSDC(sdc%nnodes))

    ! destroy
    do m = 1, sdc%nnodes
       call destroy(uSDC(m))
       call destroy(fSDC(m))
    end do

    do m = 1, sdc%nnodes-1
       call destroy(S(m))
    end do

  end subroutine advance_sdc


  function sdc_residual (uSDC,fSDC,R,dt,sdc) result(res)
    real(dp_t)                      :: res
    type(sdcquad),    intent(in   ) :: sdc
    type(multifab),   intent(inout) :: uSDC(sdc%nnodes), fSDC(sdc%nnodes), R
    real(dp_t),       intent(in   ) :: dt

    integer :: m, n

    ! compute integral
    call copy(R, uSDC(1))

    do m = 1, sdc%nnodes-1
       do n = 1, sdc%nnodes
          call saxpy(R, dt*sdc%smat(m,n), fSDC(n))
       end do
    end do

    call saxpy(R, -1.0_dp_t, uSDC(sdc%nnodes))

    res = norm_inf(R)

  end function sdc_residual


  subroutine advance_rk4(u,a,b,g,dx,time,dt,prob_lo,prob_hi)

    type(multifab) , intent(inout) :: u,g,a,b
    real(kind=dp_t), intent(in   ) :: dx, time, dt
    real(kind=dp_t), intent(in   ) :: prob_lo(u%dim)
    real(kind=dp_t), intent(in   ) :: prob_hi(u%dim)

    ! local variables
    type(multifab) :: dudt, u1, u2, u3
    real(kind=dp_t) :: t

    call multifab_build(dudt,u%la,1,0)
    call multifab_build(u1,u%la,1,u%ng)
    call multifab_build(u2,u%la,1,u%ng)
    call multifab_build(u3,u%la,1,u%ng)

    ! step 1
    t = time

    call compute_dudt(u, a, b, g, dudt, t, dx, prob_lo)

    call multifab_mult_mult_s(dudt, 0.5_dp_t*dt)
    call multifab_copy(u1, u)
    call multifab_plus_plus(u1, dudt)

    ! step 2
    t = time + 0.5_dp_t*dt

    call compute_dudt(u1, a, b, g, dudt, t, dx, prob_lo)

    call multifab_mult_mult_s(dudt, 0.5_dp_t*dt)
    call multifab_copy(u2, u)
    call multifab_plus_plus(u2, dudt)
    
    ! step 3
    t = time + 0.5_dp_t*dt

    call compute_dudt(u2, a, b, g, dudt, t, dx, prob_lo)
    
    call multifab_mult_mult_s(dudt, dt)
    call multifab_copy(u3, u)
    call multifab_plus_plus(u3, dudt)
    
    ! step 4
    t = time + dt

    call compute_dudt(u3, a, b, g, dudt, t, dx, prob_lo)
    
    call multifab_mult_mult_s(dudt, dt/6._dp_t)

    call multifab_mult_mult_s(u,-1._dp_t/3._dp_t)
    call multifab_plus_plus(u, dudt)
    
    call multifab_mult_mult_s(u1,1._dp_t/3._dp_t)
    call multifab_mult_mult_s(u2,2._dp_t/3._dp_t)
    call multifab_mult_mult_s(u3,1._dp_t/3._dp_t)

    call multifab_plus_plus(u, u1)
    call multifab_plus_plus(u, u2)
    call multifab_plus_plus(u, u3)

    call multifab_destroy(dudt)
    call multifab_destroy(u1)
    call multifab_destroy(u2)
    call multifab_destroy(u3)

  end subroutine advance_rk4


  function compute_dt(a, b, dx) result(dt)
    real(dp_t), intent(in) :: dx
    type(multifab) , intent(in) :: a, b
    real(dp_t) :: dt
    real(dp_t) :: abmax, abmax_proc
!    abmax_proc = max(max_val(a), max_val(b))
!    call parallel_reduce(abmax, abmax_proc, MPI_MAX)
    abmax = diffcoef_max
    dt = 0.4_dp_t * dx**2 / (2._dp_t*abmax)
  end function compute_dt


  subroutine compute_dudt(u, a, b, g, dudt, t, dx, prob_lo)
    type(multifab), intent(inout) :: u, a, b, g, dudt
    real(dp_t), intent(in) :: t, dx, prob_lo(u%dim)
    if (stencil .eq. narrow) then
       call compute_dudt_narrow(u, a, b, g, dudt, t, dx, prob_lo)
    else if (stencil .eq. wide) then
       call compute_dudt_wide(u, a, b, g, dudt, t, dx, prob_lo)
    else if (stencil .eq. chain) then
       call compute_dudt_chain(u, a, b, g, dudt, t, dx, prob_lo)
    else
       call bl_error('Unknown stencil type')
    end if
  end subroutine compute_dudt

  subroutine compute_dudt_narrow(u, a, b, g, dudt, t, dx, prob_lo)
    type(multifab), intent(inout) :: u, a, b, g, dudt
    real(dp_t), intent(in) :: t, dx, prob_lo(u%dim)

    integer :: i, dm, ng
    type(layout)   :: la
    type(multifab) :: flux(U%dim)

    dm = u%dim
    ng = nghost(u)
    la = get_layout(u)

    do i=1,dm
       ! flux(i) has one component, zero ghost cells, and is nodal in direction i
       call multifab_build_edge(flux(i),u%la,1,0,i)
    end do

    call multifab_fill_boundary(u)
    call compute_diffcoef(a,b,u,dx,prob_lo)
    call multifab_fill_boundary(a)
    call multifab_fill_boundary(b)
    call compute_source(g,t,dx,prob_lo)

    call compute_flux(flux(1), a, u, 1)
    call compute_flux(flux(2), b, u, 2)

    call add_dudt(dudt, flux, g, dx)

    do i=1,dm
       call multifab_destroy(flux(i))
    end do

  end subroutine compute_dudt_narrow

  subroutine compute_flux(flux, a, u, idir)

    integer, intent(in) :: idir
    type(multifab) , intent(in   ) :: a, u
    type(multifab) , intent(inout) :: flux
    
    real(kind=dp_t), pointer :: du(:,:,:,:), da(:,:,:,:), df(:,:,:,:)
    integer :: lo(u%dim), hi(u%dim)
    integer :: dm, ng, i

    dm = u%dim
    ng = u%ng

    do i=1,nfabs(u)
       df => dataptr(flux,i)
       da => dataptr(a,i)
       du => dataptr(u,i)
       lo = lwb(get_box(u,i))
       hi = upb(get_box(u,i))
       select case(dm)
       case (2)
          call comp_flux_2d(du(:,:,1,1), da(:,:,1,1), df(:,:,1,1), &
               ng, lo, hi, idir)
       case default
          call bl_error('dim must be 2')
       end select       
    end do

  end subroutine compute_flux


  subroutine comp_flux_2d(u, a, f, ng, lo, hi, idir)

    integer          :: lo(2), hi(2), ng, idir
    real(dp_t) :: u(lo(1)-ng:,lo(2)-ng:), a(lo(1)-ng:,lo(2)-ng:) 
    real(dp_t) :: f(lo(1):,lo(2):)
 
    ! local varables
    integer          :: i,j

    if (D2_order .eq. 8) then
       if (idir .eq. 1) then
          !$omp parallel do private(i,j)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                f(i,j) = dot_product(matmul(a(i-4:i+3,j),M8), &
                     &                      u(i-4:i+3,j))
             end do
          end do
          !$omp end parallel do
       else if (idir .eq. 2) then
          !$omp parallel do private(i,j)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                f(i,j) = dot_product(matmul(a(i,j-4:j+3),M8), &
                     &                      u(i,j-4:j+3))
             end do
          end do
          !$omp end parallel do
       else 
          call bl_error("idir = 3???")
       end if
    else
       if (idir .eq. 1) then
          !$omp parallel do private(i,j)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                f(i,j) = dot_product(matmul(a(i-3:i+2,j),M6), &
                     &                      u(i-3:i+2,j))
             end do
          end do
          !$omp end parallel do
       else if (idir .eq. 2) then
          !$omp parallel do private(i,j)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                f(i,j) = dot_product(matmul(a(i,j-3:j+2),M6), &
                     &                      u(i,j-3:j+2))
             end do
          end do
          !$omp end parallel do
       else 
          call bl_error("idir = 3???")
       end if
    end if

  end subroutine comp_flux_2d

  subroutine add_dudt(dudt, flux, g, dx)
    type(multifab), intent(inout) :: dudt
    type(multifab), intent(in   ) :: flux(2), g
    real(kind=dp_t), intent(in) :: dx

    real(kind=dp_t), pointer :: ddudt(:,:,:,:), dfx(:,:,:,:), dfy(:,:,:,:), dg(:,:,:,:)
    integer :: lo(dudt%dim), hi(dudt%dim)
    integer :: dm, ng,i

    dm = dudt%dim
    ng = g%ng

    do i=1,nfabs(dudt)
       dfx => dataptr(flux(1),i)
       dfy => dataptr(flux(2),i)
       ddudt => dataptr(dudt,i)
       dg    => dataptr(g   ,i)
       lo = lwb(get_box(dudt,i))
       hi = upb(get_box(dudt,i))
       select case(dm)
       case (2)
          call add_dudt_2d(ddudt(:,:,1,1), dfx(:,:,1,1), dfy(:,:,1,1), dg(:,:,1,1), &
               ng, lo, hi, dx)
       case default
          call bl_error('dim must be 2')
       end select       
    end do

  end subroutine add_dudt

  subroutine add_dudt_2d(dudt, fx, fy, g, ng, lo, hi, dx)
    integer :: lo(2), hi(2), ng
    real(dp_t) :: dx
    real(dp_t) :: dudt(lo(1):,lo(2):), fx(lo(1):,lo(2):), fy(lo(1):,lo(2):)
    real(dp_t) :: g(lo(1)-ng:,lo(2)-ng:)

    integer :: i,j
    real(dp_t) dx2

    dx2 = 1._dp_t / dx**2

    !$omp parallel do private(i,j)
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          dudt(i,j) = dx2 * (fx(i+1,j) - fx(i,j) + fy(i,j+1) - fy(i,j)) + g(i,j)
       end do
    end do
    !$omp end parallel do
  end subroutine add_dudt_2d


  subroutine compute_dudt_wide(u, a, b, g, dudt, t, dx, prob_lo)
    type(multifab), intent(inout) :: u, a, b, g, dudt
    real(dp_t), intent(in) :: t, dx, prob_lo(u%dim)
    
    integer :: ib, i, j, dm, ng
    type(layout)   :: la
    type(multifab) :: ux, uy, fx, fy
    integer :: lo(u%dim), hi(u%dim)
    real(dp_t) :: dx2
    real(kind=dp_t), pointer :: fxp(:,:,:,:), fyp(:,:,:,:), gp(:,:,:,:), upp(:,:,:,:)

    call multifab_fill_boundary(u)

    dm = u%dim
    ng = nghost(u)
    la = get_layout(u)
    
    call multifab_build(ux, la, 1, ng)
    call multifab_build(uy, la, 1, ng)
    call multifab_build(fx, la, 1, 0)
    call multifab_build(fy, la, 1, 0)

    call compute_diffcoef(a,b,u,dx,prob_lo)
    call compute_source(g,t,dx,prob_lo)

    call compute_xderiv(u, ux)
    call compute_yderiv(u, uy)

    call multifab_mult_mult(ux, a)
    call multifab_mult_mult(uy, b)
    
    call multifab_fill_boundary(ux)
    call multifab_fill_boundary(uy)

    call compute_xderiv(ux, fx)
    call compute_yderiv(uy, fy)

    dx2 = 1._dp_t/dx**2
    do ib=1,nfabs(u)
       fxp => dataptr(fx, ib)
       fyp => dataptr(fy, ib)
       gp  => dataptr(g , ib)
       upp => dataptr(dudt, ib)
       lo = lwb(get_box(u,ib))
       hi = upb(get_box(u,ib))

       !$omp parallel do private(i,j)
       do j=lo(2), hi(2)
          do i=lo(1), hi(1)
             upp(i,j,1,1) = dx2*(fxp(i,j,1,1)+fyp(i,j,1,1)) + gp(i,j,1,1)
          end do
       end do
       !$omp end parallel do
    end do

    call multifab_destroy(ux)
    call multifab_destroy(uy)
    call multifab_destroy(fx)
    call multifab_destroy(fy)

  end subroutine compute_dudt_wide
  
  subroutine compute_xderiv(u, ux)
    type(multifab), intent(in) :: u
    type(multifab), intent(inout) :: ux

    integer :: ib, i, j
    integer :: lo(u%dim), hi(u%dim)
    real(kind=dp_t), pointer :: up(:,:,:,:), uxp(:,:,:,:)

    do ib=1,nfabs(u)
       up  => dataptr(u , ib)
       uxp => dataptr(ux, ib)
       lo = lwb(get_box(u,ib))
       hi = upb(get_box(u,ib))

       !$omp parallel do private(i,j)
       do j=lo(2), hi(2)
          do i=lo(1), hi(1)
             uxp(i,j,1,1) = first_deriv_8(up(i-4:i+4,j,1,1))
          end do
       end do
       !$omp end parallel do
    end do
  end subroutine compute_xderiv

  subroutine compute_yderiv(u, uy)
    type(multifab), intent(in) :: u
    type(multifab), intent(inout) :: uy

    integer :: ib, i, j
    integer :: lo(u%dim), hi(u%dim)
    real(kind=dp_t), pointer :: up(:,:,:,:), uyp(:,:,:,:)

    do ib=1,nfabs(u)
       up  => dataptr(u , ib)
       uyp => dataptr(uy, ib)
       lo = lwb(get_box(u,ib))
       hi = upb(get_box(u,ib))

       !$omp parallel do private(i,j)
       do j=lo(2), hi(2)
          do i=lo(1), hi(1)
             uyp(i,j,1,1) = first_deriv_8(up(i,j-4:j+4,1,1))
          end do
       end do
       !$omp end parallel do
    end do
  end subroutine compute_yderiv


  subroutine compute_dudt_chain(u, a, b, g, dudt, t, dx, prob_lo)
    type(multifab), intent(inout) :: u, a, b, g, dudt
    real(dp_t), intent(in) :: t, dx, prob_lo(u%dim)
    
    integer :: ib, i, j, dm, ng
    type(layout)   :: la
    type(multifab) :: ax, by, ux, uy, fx, fy, one
    integer :: lo(u%dim), hi(u%dim)
    real(dp_t) :: dx2
    real(kind=dp_t), dimension(:,:,:,:), pointer :: fxp, fyp, gp, upp, &
         axp, byp, uxp, uyp, ap, bp

    call multifab_fill_boundary(u)

    dm = u%dim
    ng = nghost(u)
    la = get_layout(u)
    
    call multifab_build(ux, la, 1, 0)
    call multifab_build(uy, la, 1, 0)

    call multifab_build(ax, la, 1, 0)
    call multifab_build(by, la, 1, 0)

    call multifab_build_edge(fx, la, 1, 0, 1)
    call multifab_build_edge(fy, la, 1, 0, 2)

    call multifab_build(one, la, 1, ng)
    call setval(one, 1._dp_t, all=.true.)

    call compute_diffcoef(a,b,u,dx,prob_lo)
    call compute_source(g,t,dx,prob_lo)

    call compute_xderiv(u, ux)
    call compute_yderiv(u, uy)

    call multifab_fill_boundary(a)
    call multifab_fill_boundary(b)

    call compute_xderiv(a, ax)
    call compute_yderiv(b, by)

    call compute_flux(fx, one, u, 1)
    call compute_flux(fy, one, u, 2)

    dx2 = 1._dp_t/dx**2
    do ib=1,nfabs(u)
       fxp => dataptr(fx, ib)
       fyp => dataptr(fy, ib)
       gp  => dataptr(g , ib)
       upp => dataptr(dudt, ib)
       axp => dataptr(ax, ib)
       byp => dataptr(by, ib)
       uxp => dataptr(ux, ib)
       uyp => dataptr(uy, ib)
       ap  => dataptr(a , ib)
       bp  => dataptr(b , ib)

       lo = lwb(get_box(u,ib))
       hi = upb(get_box(u,ib))

       !$omp parallel do private(i,j)
       do j=lo(2), hi(2)
          do i=lo(1), hi(1)
             upp(i,j,1,1) = dx2 * (ap(i,j,1,1)*(fxp(i+1,j,1,1) - fxp(i,j,1,1)) &
                  &              + bp(i,j,1,1)*(fyp(i,j+1,1,1) - fyp(i,j,1,1)) &
                  &           + axp(i,j,1,1)*uxp(i,j,1,1) + byp(i,j,1,1)*uyp(i,j,1,1)) &
                  + gp(i,j,1,1)
          end do
       end do
       !$omp end parallel do
    end do

    call multifab_destroy(ux)
    call multifab_destroy(uy)
    call multifab_destroy(ax)
    call multifab_destroy(by)
    call multifab_destroy(fx)
    call multifab_destroy(fy)
    call multifab_destroy(one)

  end subroutine compute_dudt_chain

end module advance_module

