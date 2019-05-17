#ifndef EOS_H_
#define EOS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace eos {

void init ();
void finalize ();

AMREX_GPU_DEVICE void test ();

extern AMREX_GPU_DEVICE_MANAGED amrex::Real eos_scalar;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_array;
extern AMREX_GPU_DEVICE_MANAGED int array_size;

}

#endif
