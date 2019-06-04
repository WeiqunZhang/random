program main

  use boxlib
  use multifab_module
  use bl_IO_module
  use layout_module
  use prob_module
  use write_plotfile_module
  use advance_module
  use derivative_stencil_module, only : stencil_init

  implicit none

  ! stuff you can set with the inputs file (otherwise use default values below)
  integer :: prob_type, D2_order, dim, nsteps, plot_int, n_cell, verbose, max_grid_size
  integer :: advance_method, sdc_nnodes, sdc_iters
  real(dp_t) :: sdc_tol_residual
  real(dp_t) :: O8_m47, O8_m48, O6_m36
  real(dp_t) :: stop_time
  character(len=16) :: stencil_type

  ! dummy indices using for reading in inputs file
  integer :: un, farg, narg
  logical :: need_inputs_file, found_inputs_file
  character(len=128) :: inputs_file_name
  character(len=128) :: fname

  integer, allocatable :: lo(:), hi(:)
  integer :: istep

  real(dp_t), allocatable :: prob_lo(:), prob_hi(:)
  real(dp_t) :: dx, dt, time, l0, l2
  
  logical, allocatable :: is_periodic(:)

  type(box)      :: bx
  type(boxarray) :: ba
  type(layout)   :: la
  type(multifab) :: u, a, b, g
  integer :: ng
  logical :: laststep

  namelist /probin/ prob_type, D2_order, dim, nsteps, plot_int, n_cell, &
       advance_method, sdc_nnodes, sdc_iters, sdc_tol_residual, &
       O8_m47, O8_m48, o6_m36, stop_time, verbose, stencil_type, &
       max_grid_size

  ! if running in parallel, this will print out the number of MPI 
  ! processes and OpenMP threads
  call boxlib_initialize()

  ! default values - will get overwritten by the inputs file
  verbose       = 0
  prob_type     = 1
  dim           = 2
  nsteps        = 1000
  plot_int      = 100
  stop_time     = 1._dp_t
  n_cell        = 20
  max_grid_size = 256

  advance_method = 2
  sdc_nnodes = 3
  sdc_iters = 7
  sdc_tol_residual = -1.e-5_dp_t

  stencil_type = "narrow"
  D2_order = 8

! optimized for more zeros
!  O8_m47 = 683._dp_t/10080._dp_t
!  O8_m48 = -1._dp_t/224._dp_t
! optimized for leading order truncation error assuming equal weight for the error terms
  O8_m47 = 3557._dp_t/44100._dp_t
  O8_m48 = -2083._dp_t/117600._dp_t

  O6_m36 = 1._dp_t/90._dp_t

  ! read inputs file and overwrite any default values
  narg = command_argument_count()
  need_inputs_file = .true.
  farg = 1
  if ( need_inputs_file .AND. narg >= 1 ) then
     call get_command_argument(farg, value = inputs_file_name)
     inquire(file = inputs_file_name, exist = found_inputs_file )
     if ( found_inputs_file ) then
        farg = farg + 1
        un = unit_new()
        open(unit=un, file = inputs_file_name, status = 'old', action = 'read')
        read(unit=un, nml = probin)
        close(unit=un)
        need_inputs_file = .false.
     end if
  end if

  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)
     select case (fname)

     case ('--prob_type')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) prob_type

     case ('--verbose')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) verbose

     case ('--stop_time')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) stop_time

     case ('--n_cell')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) n_cell

     case ('--max_grid_size')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) max_grid_size

     case ('--advance_method')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) advance_method

     case ('--sdc_nnodes')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) sdc_nnodes

     case ('--sdc_iters')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) sdc_iters

     case ('--sdc_tol_residual')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) sdc_tol_residual

     case ('--D2_order')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) D2_order

     case ('--stencil_type')
        farg = farg + 1
        call get_command_argument(farg, value = stencil_type)

     case ('--')
        farg = farg + 1
        exit
        
     case default
        if ( .not. parallel_q() ) then
           write(*,*) 'UNKNOWN option = ', fname
           call bl_error("MAIN")
        end if
     end select
     
     farg = farg + 1
  end do

  call stencil_init(D2_order, O8_m47, O8_m48, O6_m36, stencil_type)
  call advance_init(advance_method, verbose, sdc_nnodes, sdc_iters, sdc_tol_residual)

  ! now that we have dim, we can allocate these
  allocate(lo(dim),hi(dim))
  allocate(is_periodic(dim))
  allocate(prob_lo(dim),prob_hi(dim))

  prob_lo(:) = 0._dp_t
  prob_hi(:) = atan(1._dp_t) * 8._dp_t
  is_periodic(:) = .true.

  ! create a box from (0,0) to (n_cell-1,n_cell-1)
  lo(:) = 0
  hi(:) = n_cell-1
  bx = make_box(lo,hi)

  ! the grid spacing is the same in each direction
  dx = (prob_hi(1)-prob_lo(1)) / n_cell
  prob_lo = prob_lo - 0.5_dp_t*dx
  prob_hi = prob_hi - 0.5_dp_t*dx

  ! initialize the boxarray to be one single box
  call boxarray_build_bx(ba,bx)
  call boxarray_maxsize(ba,max_grid_size)
  call layout_build_ba(la,ba,bx,pmask=is_periodic)

  call destroy(ba)

  ! build multifab with 1 component and ghost cells
  ng = D2_order / 2
  call multifab_build(u,la,1,ng)
  call multifab_build(a,la,1,ng)
  call multifab_build(b,la,1,ng)
  call multifab_build(g,la,1,0)
  
  call multifab_setval(u, 0._dp_t, all=.true.)
  call multifab_setval(a, 0._dp_t, all=.true.)
  call multifab_setval(b, 0._dp_t, all=.true.)
  call multifab_setval(g, 0._dp_t, all=.true.)

  ! initialze phi
  call init_prob(u,dx,prob_lo,prob_type)

  istep = 0
  time = 0._dp_t

  call multifab_fill_boundary(u)
  call compute_diffcoef(a,b,u,dx,prob_lo)
  dt = compute_dt(a, b, dx)

  if ( parallel_IOProcessor() ) then
     print *, 'dt = ', dt
     call flush(6)
  end if

  ! write out plotfile 0
  call write_plotfile(la,u,istep,dx,time,prob_lo,prob_hi)

  laststep = .false.
  do istep=1,nsteps

     if (time + dt .gt. stop_time) then
        dt = stop_time - time
        laststep = .true.
     end if

     call advance(u,a,b,g,dx,time,dt,prob_lo,prob_hi)

     time = time + dt

     if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps .or. laststep) then
        ! write out plotfile
        call write_plotfile(la,u,istep,dx,time,prob_lo,prob_hi)
     end if

     if (laststep) then
        exit
     end if

  end do

  call compute_norms(u, dx, time, prob_lo, l0, l2)

  if (parallel_IOProcessor()) then
     if (prob_type .eq. 1) then
        print *, 'n_cell = ', n_cell
        print *, 'L0 error = ', L0
        print *, 'L2 error = ', L2
        call flush(6)
     end if
  end if

  ! make sure to destroy the multifab or you'll leak memory
  call destroy(u)
  call destroy(a)
  call destroy(b)
  call destroy(g)
  call destroy(la)

  call advance_close()

  deallocate(lo,hi,is_periodic,prob_lo,prob_hi)

  ! deallocate temporary boxarrays and communication mappings
  call layout_flush_copyassoc_cache()

  call boxlib_finalize()

end program main
