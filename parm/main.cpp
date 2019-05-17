
#include <AMReX.H>
#include <AMReX_Gpu.H>
#include "eos.H"

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        eos::init();

        // launch 1 thread
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int i)
        {
            eos::test();
        });

        eos::finalize();
    }
    amrex::Finalize();
}

