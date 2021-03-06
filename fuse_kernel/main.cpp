
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_FArrayBox.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        Box bx1(IntVect(0),IntVect(15));
        FArrayBox fab1(bx1,1);
        Array4<Real> arr1 = fab1.array();
        auto f1 = [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr1(i,j,k) = 2.0;
        };

        static_assert(std::is_trivially_copyable<decltype(f1)>::value,"xxxxx");

        Box bx2(IntVect(0),IntVect(31));
        FArrayBox fab2(bx2,5);
        Array4<Real> arr2 = fab2.array();
        auto f2 = [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            arr2(i,j,k,n) = 4.0;
        };

        Gpu::DeviceVector<int> vec(1000);
        int* pvec = vec.data();
        auto f3 = [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            pvec[i] = 8;
        };

        Gpu::Fuser fuser;
        fuser.Register(bx1, f1);
        fuser.Register(bx2, fab2.nComp(), f2);
        fuser.Register(vec.size(), f3);
        fuser.Launch();
    }
    amrex::Finalize();
}

