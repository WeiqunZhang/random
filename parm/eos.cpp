#include "eos.H"
#include <AMReX_Arena.H>
#include <cstdio>

using namespace amrex;

namespace eos {

AMREX_GPU_DEVICE_MANAGED Real eos_scalar;
AMREX_GPU_DEVICE_MANAGED Real* eos_array;
AMREX_GPU_DEVICE_MANAGED int array_size;

void init ()
{
    eos_scalar = 3.14;
    array_size = 5;
    eos_array = static_cast<Real*>(The_Managed_Arena()->alloc(sizeof(Real)*array_size));
    for (int i = 0; i < array_size; ++i) {
        eos_array[i] = i*1.1;
    }
}

void finalize ()
{
    The_Managed_Arena()->free(eos_array);
}

AMREX_GPU_DEVICE
void test ()
{
    printf("eos_scalar = %f\n", eos_scalar);

    for (int i = 0; i < array_size; ++i) {
        printf("eos_array[%d] = %f\n", i, eos_array[i]);
    }
}

}

