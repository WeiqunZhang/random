
#include "MyTest.H"
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void
MyTest::writePlotfile () const
{
    const int ncomp = AMREX_SPACEDIM + 3;
    Vector<std::string> varname =
        {AMREX_D_DECL("vx", "vy", "vz"), "exact", "error", "vfrc"};

    MultiFab plotmf(grids, dmap, ncomp, 0);
    MultiFab::Copy(plotmf, solution, 0, 0               , AMREX_SPACEDIM, 0);
    MultiFab::Copy(plotmf, exact   , 0, AMREX_SPACEDIM  , 1             , 0);
    MultiFab::Copy(plotmf, solution, flow_dir, AMREX_SPACEDIM+1, 1      , 0);
    MultiFab::Subtract(plotmf, exact, 0, AMREX_SPACEDIM+1, 1, 0);
    const MultiFab& vfrc = factory->getVolFrac();
    MultiFab::Copy(plotmf, vfrc, 0, AMREX_SPACEDIM+2, 1, 0);

    WriteMultiLevelPlotfile("plot", 1, {&plotmf},
                            varname, {geom}, 0.0, {0}, {IntVect(2)});
}
