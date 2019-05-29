#include "MyTest.H"

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initializeEB();

    initData();
}

//
// Given vel, rhs & sig, this solves Div (sig * Grad phi) = Div vel + rhs.
// On return, vel becomes vel  - sig * Grad phi.
//
void
MyTest::solve ()
{
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            mlmg_lobc[idim] = LinOpBCType::Periodic;
            mlmg_hibc[idim] = LinOpBCType::Periodic;
        } else {
            mlmg_lobc[idim] = LinOpBCType::Neumann;
            mlmg_hibc[idim] = LinOpBCType::Neumann;
        }
    }

    mlmg_lobc[0] = LinOpBCType::Neumann;
    mlmg_hibc[0] = LinOpBCType::Dirichlet;

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);
//    info.setAgglomeration(false);
//    info.setConsolidation(false);

    MLNodeLaplacian mlndlap({geom}, {grids}, {dmap}, info,
                            Vector<EBFArrayBoxFactory const*>{factory.get()});

//    mlndlap.setGaussSeidel(false);

    mlndlap.setDomainBC(mlmg_lobc, mlmg_hibc);

    {
        MultiFab tmp;
        VisMF::Read(tmp, "b");
        sig.ParallelCopy(tmp);

        EB_set_covered(sig, 0.0);
        mlndlap.setSigma(0, sig);
    }

    mlndlap.setNormalizationThreshold(1.e-10);

//    const auto& vfrc = factory[0]->getVolFrac();
//    VisMF::Write(vfrc, "vfrc");

    MLMG mlmg(mlndlap);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);

//    mlmg.setBottomMaxIter(100);

    mlmg.setBottomTolerance(1.e-4);

//    mlmg.setBottomSolver(MLMG::BottomSolver::cg);
//    mlmg.setBottomSolver(MLMG::BottomSolver::cgbicg);
    mlmg.setBottomSolver(MLMG::BottomSolver::bicgcg);
//    mlmg.setBottomSolver(MLMG::BottomSolver::smoother);

//    mlmg.setMaxIter(1);
//    mlmg.setFinalSmooth(200);

    phi.setVal(0.0);
    {
        MultiFab tmp;
        VisMF::Read(tmp, "phi_nd");
        phi.ParallelCopy(tmp);
    }

    rhs.setVal(0.0);
    {
        MultiFab tmp;
        VisMF::Read(tmp, "rhs");
        rhs.ParallelCopy(tmp);
    }

    Real mlmg_err = mlmg.solve({&phi}, {&rhs}, 1.e-11, 0.0);
}

void
MyTest::readParameters ()
{   
    ParmParse pp;
    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("max_coarsening_level", max_coarsening_level);

    ParmParse pp_eb("eb2");
    std::string geom_type;
    pp_eb.get("cylinder_direction", cylinder_direction);
}

void
MyTest::initGrids ()
{
    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(0.096,0.024,0.024)});
    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    Box domain(IntVect{0,0,0}, IntVect{511,127,127});
    geom.define(domain, &rb, 0, is_periodic.data());

    grids.define(domain);
    grids.maxSize(64);
}

void
MyTest::initData ()
{
    dmap.define(grids);
    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
    const EB2::Level& eb_level = eb_is.getLevel(geom);
    factory.reset(new EBFArrayBoxFactory(eb_level, geom, grids, dmap,
                                         {2,2,2}, EBSupport::full));

    phi.define(amrex::convert(grids,IntVect::TheNodeVector()),
               dmap, 1, 1, MFInfo(), *factory);
    rhs.define(amrex::convert(grids,IntVect::TheNodeVector()),
               dmap, 1, 0, MFInfo(), *factory);
    sig.define(grids, dmap, 1, 1, MFInfo(), *factory);
}
