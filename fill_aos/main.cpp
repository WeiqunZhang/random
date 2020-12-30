
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BaseFabUtility.H>

using namespace amrex;

void main_main ()
{
    BL_PROFILE("main");
//    constexpr int N = 7;
    constexpr int N = 27;
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
        amrex::fill<N>(aos_fab,
        [=] AMREX_GPU_HOST_DEVICE (GpuArray<Real,N>& a, int /*i*/, int /*j*/, int /*k*/)
        {
            a[0] = Real(N-1);
            for (int n = 1; n < N; ++n) {
                a[n] = -1.;
            }
        });
        Gpu::synchronize();
    }

    {
        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (aos(i,j,k)[0] != Real(N-1)) {
#if AMREX_DEVICE_COMPILE
                AMREX_DEVICE_PRINTF("aos(%d,%d,%d)[0] = %g\n",i,j,k,aos(i,j,k)[0]);
#endif
            }
            for (int n = 1; n < N; ++n) {
                if (aos(i,j,k)[n] != -1.) {
#if AMREX_DEVICE_COMPILE
                    AMREX_DEVICE_PRINTF("aos(%d,%d,%d)[%d] = %g\n",i,j,k,n,aos(i,j,k)[n]);
#endif
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

