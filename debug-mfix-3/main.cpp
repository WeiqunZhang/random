
#include <AMReX.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

//    for (int i = 0; i < 100; ++i)
    {
        BL_PROFILE("main");
        MyTest mytest;
        for (int i = 0; i < 1; ++i) {
            mytest.solve();
            mytest.writePlotfile();
        }
//        amrex::Print() << "Test # " << i << " done." << std::endl;
    }

    amrex::Finalize();
}
