#include "OutFlowPressure.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        Box box(IntVect(0), IntVect(31));
        BoxArray ba(box);
        ba.maxSize(16);
        DistributionMapping dm{ba};
        MultiFab rho(ba,dm,1,1);
        rho.setVal(1.0);
        MultiFab p(amrex::convert(ba,IntVect::TheNodeVector()),dm,1,0);
        p.setVal(-100.);

        RealBox rb(0.,0.,0.,1.,1.,1.);
        Array<int,AMREX_SPACEDIM> is_periodic = {0,0,0};
        Geometry geom(box, rb, 0, is_periodic);

        Array<bool,AMREX_SPACEDIM*2> is_outflow;
        is_outflow.fill(false);
        is_outflow[Orientation(Direction::y,Orientation::low)] = true;
        is_outflow[Orientation(Direction::y,Orientation::high)] = true;

        OutFlowPressure ofp(geom, 1.0, is_outflow);
        ofp.computePressure(p, rho);
        VisMF::Write(p,"p");
    }
    amrex::Finalize();
}

