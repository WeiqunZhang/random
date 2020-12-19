
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}

void main_main ()
{
    int N = 10000000;
    int ncell = 256;
    int max_grid_size = 128;
    int nghost = 0;
    {
        ParmParse pp;
        pp.query("n", N);
        pp.query("ncell", ncell);
        pp.query("max_grid_size", max_grid_size);
        pp.query("nghost", nghost);
    }

    Gpu::DeviceVector<Real> vec(N, 1.0);
    Real* pvec = vec.data();

    MultiFab mf;
    {
        BoxArray ba;
        Box domain_box(IntVect(0), IntVect(ncell-1));
        ba.define(domain_box);
        ba.maxSize(max_grid_size);
        mf.define(ba,DistributionMapping{ba},1,nghost);
        mf.setVal(3.0);
    }

    { // warm up
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(N, reduce_data,
        [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
        {
            return pvec[i];
        });

        Real hv = amrex::get<0>(reduce_data.value());
        amrex::Print().SetPrecision(17) << "vec sum: " << hv << "\n";
    }

    Gpu::synchronize();
    Real t_vec_orig = amrex::second();
    {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(N, reduce_data,
        [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
        {
            return pvec[i];
        });

        Real hv = amrex::get<0>(reduce_data.value());
        amrex::Print().SetPrecision(17) << "vec sum: " << hv << "\n";
    }
    Gpu::synchronize();
    t_vec_orig = amrex::second() - t_vec_orig;

    { // warm up
        DeviceBuffer<Real> da({0.0});
        Real* dp = da.data();
        amrex::ParallelFor(Gpu::KernelInfo().setReduction(true), N,
        [=] AMREX_GPU_DEVICE (int i, Gpu::Handler const& handler) noexcept
        {
            Gpu::deviceReduceSum(dp, pvec[i], handler);
        });
        Real* hp = da.copyToHost();
        amrex::Print().SetPrecision(17) << "vec sum: "  << hp[0] << "\n";
    }

    Gpu::synchronize();
    Real t_vec_new = amrex::second();
    {
        DeviceBuffer<Real> da({0.0});
        Real* dp = da.data();
        amrex::ParallelFor(Gpu::KernelInfo().setReduction(true), N,
        [=] AMREX_GPU_DEVICE (int i, Gpu::Handler const& handler) noexcept
        {
            Gpu::deviceReduceSum(dp, pvec[i], handler);
        });
        Real* hp = da.copyToHost();
        amrex::Print().SetPrecision(17) << "vec sum: "  << hp[0] << "\n";
    }
    Gpu::synchronize();
    t_vec_new = amrex::second() - t_vec_new;

    { // warm up
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.fabbox();
            auto const& fab = mf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return fab(i,j,k);
            });
        }

        Real hv = amrex::get<0>(reduce_data.value());
        amrex::Print().SetPrecision(17) << "mf sum: " << hv << "\n";
    }

    Gpu::synchronize();
    Real t_mf_orig = amrex::second();
    {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.fabbox();
            auto const& fab = mf.array(mfi);
            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return fab(i,j,k);
            });
        }

        Real hv = amrex::get<0>(reduce_data.value());
        amrex::Print().SetPrecision(17) << "mf sum: " << hv << "\n";
    }
    Gpu::synchronize();
    t_mf_orig = amrex::second() - t_mf_orig;

    { // warm up
        DeviceBuffer<Real> da({0.0});
        Real* dp = da.data();
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.fabbox();
            Array4<Real const> const& fab = mf.const_array(mfi);
            amrex::ParallelFor(Gpu::KernelInfo().setReduction(true), bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Gpu::Handler const& handler) noexcept
            {
                Gpu::deviceReduceSum(dp, fab(i,j,k), handler);
            });
        }
        Real* hp = da.copyToHost();
        amrex::Print().SetPrecision(17) << "mf sum: " << hp[0] << "\n";
    }

    Gpu::synchronize();
    Real t_mf_new = amrex::second();
    {
        DeviceBuffer<Real> da({0.0});
        Real* dp = da.data();
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.fabbox();
            Array4<Real const> const& fab = mf.const_array(mfi);
            amrex::ParallelFor(Gpu::KernelInfo().setReduction(true), bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Gpu::Handler const& handler) noexcept
            {
                Gpu::deviceReduceSum(dp, fab(i,j,k), handler);
            });
        }
        Real* hp = da.copyToHost();
        amrex::Print().SetPrecision(17) << "mf sum: " << hp[0] << "\n";
    }
    Gpu::synchronize();
    t_mf_new = amrex::second() - t_mf_new;

    amrex::Print() << "Sum of Vector of " << N << " Reals"
                   << ": ReduceOps time is " << t_vec_orig
                   << ", ParallelFor time is " << t_vec_new << "\n"
                   << "Sum of MultiFab of " << mf.boxArray().numPts() << " cells"
                   << ": ReduceOps time is " << t_mf_orig
                   << ", ParallelFor time is " << t_mf_new << std::endl;
}
