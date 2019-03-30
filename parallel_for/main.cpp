
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <Kokkos_Core.hpp>

using namespace amrex;

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        Kokkos::InitArguments args;
        args.device_id = amrex::Gpu::Device::deviceId();
        Kokkos::initialize(args);

        main_main();

        Kokkos::finalize();
    }

    amrex::Finalize();
}


void main_main ()
{
    int ncell = 128;
    int max_grid_size = 64;
    {
        ParmParse pp;
        pp.query("ncell", ncell);
        pp.query("max_grid_size", max_grid_size);
    }

    BoxArray ba;
    {
        Box domain_box(IntVect(0), IntVect(ncell-1));
        ba.define(domain_box);
        ba.maxSize(max_grid_size);
    }

    MultiFab mf(ba,DistributionMapping{ba},1,0);
    mf.setVal(0.0);

    {
        BL_PROFILE("amrex_for_3d");
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                fab(i,j,k) += 1.;
            });
        }
    }

    {
        BL_PROFILE("kokkos_for_3d");
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            const auto lo = lbound(bx);
            const auto hi = ubound(bx);
            Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3> >({lo.x,lo.y,lo.z},{hi.x+1,hi.y+1,hi.z+1}, {32,4,4}),
            KOKKOS_LAMBDA (int i, int j, int k)
            {
                fab(i,j,k) += 1.;
            });
        }
    }

    {
        BL_PROFILE("amrex_for_1d");
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = mf[mfi];
            Real* AMREX_RESTRICT p = fab.dataPtr();
            const long nitems = fab.box().numPts();
            amrex::ParallelFor(nitems,
            [=] AMREX_GPU_DEVICE (long idx)
            {
                p[idx] += 1.;
            });
        }
    }

    {
        BL_PROFILE("kokkos_for_1d");
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = mf[mfi];
            Real* AMREX_RESTRICT p = fab.dataPtr();
            const long nitems = fab.box().numPts();
            Kokkos::parallel_for(nitems, 
            KOKKOS_LAMBDA (const long& idx)
            {
                p[idx] += 1.;
            });
        }
    }
}
