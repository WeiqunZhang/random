
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BaseFab.H>

using namespace amrex;

void main_main ()
{
    BL_PROFILE("main");
//    constexpr int N = 7;
    constexpr int N = 27;
    constexpr int BLOCKDIM = 128;
    int ncell = 256;
    ParmParse pp;
    pp.query("ncell", ncell);

    Box box(IntVect(0),IntVect(ncell-1));
    BaseFab<Real> soa_fab(box,N);
    BaseFab<GpuArray<Real,N> > aos_fab(box,1);

    auto const& soa = soa_fab.array();
    auto const& aos = aos_fab.array();

    amrex::ParallelFor(box,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        for (int n = 0; n < N; ++n) {
            soa(i,j,k,n) = 0.0;
            aos(i,j,k)[n] = 0.0;
        }
    });
    Gpu::synchronize();

    {
        BL_PROFILE("soa");
        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            soa(i,j,k,0) = Real(N-1);
            for (int n = 1; n < N; ++n) {
                soa(i,j,k,n) = -1.;
            }
        });
        Gpu::synchronize();
    }

    {
        BL_PROFILE("aos");
        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            aos(i,j,k)[0] = Real(N-1) * 3.2;
            for (int n = 1; n < N; ++n) {
                aos(i,j,k)[n] = -1.*3.3;
            }
        });
        Gpu::synchronize();
    }

    {
        BL_PROFILE("aos-shared");
        int ntotcells = box.numPts(); 
        int nblocks = (ntotcells+BLOCKDIM-1)/BLOCKDIM;
        std::size_t shared_mem_bytes = BLOCKDIM*N*sizeof(Real);
        Real* p = (Real*)aos_fab.dataPtr();
        amrex::launch(nblocks, BLOCKDIM, shared_mem_bytes, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept
        {
            int icell = blockDim.x*blockIdx.x+threadIdx.x;
            Gpu::SharedMemory<Real> gsm;
            Real* const shared = gsm.dataPtr();
            Real* q = shared + threadIdx.x*N;
            if (icell < ntotcells) {
                q[0] = Real(N-1);
                for (int n = 1; n < N; ++n) {
                    q[n] = -1.0;
                }
            }
            __syncthreads();
            q = p + blockDim.x*blockIdx.x*N;
            for (unsigned int m = threadIdx.x,
                     mend = amrex::min(blockDim.x, ntotcells-blockDim.x*blockIdx.x) * N;
                 m < mend; m += blockDim.x) {
                q[m] = shared[m];
            }
        });
        Gpu::synchronize();
    }

    {
        BL_PROFILE("aos-transpose");
        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            soa(i,j,k,0) = Real(N-1);
            for (int n = 1; n < N; ++n) {
                soa(i,j,k,n) = -1.;
            }
        });
        int ntotcells = box.numPts(); 
        int nblocks = (ntotcells+BLOCKDIM-1)/BLOCKDIM;
        std::size_t shared_mem_bytes = BLOCKDIM*N*sizeof(Real);
        Real* pdst = (Real*)aos_fab.dataPtr();
        Real* psrc = soa_fab.dataPtr();
        amrex::launch(nblocks, BLOCKDIM, shared_mem_bytes, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept
        {
            int icell = blockDim.x*blockIdx.x+threadIdx.x;
            Gpu::SharedMemory<Real> gsm;
            Real* const shared = gsm.dataPtr();
            Real* q = shared + threadIdx.x*N;
            if (icell < ntotcells) {
                for (int n = 0; n < N; ++n) {
                    q[n] = psrc[n*ntotcells+icell];
                }
            }
            __syncthreads();
            q = pdst + blockDim.x*blockIdx.x*N;
            for (unsigned int m = threadIdx.x,
                     mend = amrex::min(blockDim.x, ntotcells-blockDim.x*blockIdx.x) * N;
                 m < mend; m += blockDim.x) {
                q[m] = shared[m];
            }
        });
        Gpu::synchronize();
    }

    {
        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (aos(i,j,k)[0] != Real(N-1)) {
                printf("aos(%d,%d,%d)[0] = %g\n",i,j,k,aos(i,j,k)[0]);
            }
            for (int n = 1; n < N; ++n) {
                if (aos(i,j,k)[n] != -1.) {
                    printf("aos(%d,%d,%d)[%d] = %g\n",i,j,k,n,aos(i,j,k)[n]);
                }
            }
        });
        Gpu::synchronize();
    }
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}

