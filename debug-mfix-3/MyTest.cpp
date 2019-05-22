#include "MyTest.H"

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>

#if(AMREX_SPACEDIM == 3)
#include <AMReX_algoim_integrals.H>
#endif

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
//            mlmg_lobc[idim] = LinOpBCType::Dirichlet;
//            mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            mlmg_hibc[idim] = LinOpBCType::Neumann;
        }
    }

    mlmg_lobc[0] = LinOpBCType::Neumann;
//    mlmg_lobc[0] = LinOpBCType::Dirichlet;
    mlmg_hibc[0] = LinOpBCType::Dirichlet;
//    mlmg_hibc[0] = LinOpBCType::Neumann;

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setAgglomeration(false);
    info.setConsolidation(false);

    MLNodeLaplacian mlndlap(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(factory));

    if (sigma) {
        mlndlap.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::Sigma);
    }

//    mlndlap.setGaussSeidel(false);

    mlndlap.setDomainBC(mlmg_lobc, mlmg_hibc);

    AMREX_ALWAYS_ASSERT("max_level == 0");

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        MultiFab tmp;
        VisMF::Read(tmp, "solver_mf22May19/b");
        sig[ilev].ParallelCopy(tmp);
#if 0
        for (MFIter mfi(sig[ilev]); mfi.isValid(); ++mfi) {
            auto const& fab = sig[ilev].array(mfi);
            Box const& bx = mfi.validbox();
            amrex::For(bx, [=] (int i, int j, int k)
            {
                fab(i,j,k) = 0.5 + 10.*amrex::Random();
            });
        }
#endif
        EB_set_covered(sig[ilev], 0.0);
        mlndlap.setSigma(ilev, sig[ilev]);
    }

//    mlndlap.compRHS(amrex::GetVecOfPtrs(rhs), amrex::GetVecOfPtrs(vel), {}, {});

    const auto& vfrc = factory[0]->getVolFrac();
    VisMF::Write(vfrc, "vfrc");

#if 0
#if (AMREX_SPACEDIM == 2)
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(rhs[ilev], "rhs2d");
    }
#else
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(rhs[ilev], "rhs3d");
    }
#endif
#endif

    MLMG mlmg(mlndlap);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);

//    mlmg.setBottomMaxIter(100);

//    mlmg.setBottomTolerance(1.e-5);

//    mlmg.setBottomSolver(MLMG::BottomSolver::cg);
//    mlmg.setBottomSolver(MLMG::BottomSolver::cgbicg);
    mlmg.setBottomSolver(MLMG::BottomSolver::bicgcg);
//    mlmg.setBottomSolver(MLMG::BottomSolver::smoother);

//    mlmg.setMaxIter(1);
//    mlmg.setFinalSmooth(200);

    phi[0].setVal(0.0);
    {
        MultiFab tmp;
        VisMF::Read(tmp, "solver_mf22May19/phi_nd");
        phi[0].ParallelCopy(tmp);
    }

    rhs[0].setVal(0.0);
    {
        MultiFab tmp;
        VisMF::Read(tmp, "solver_mf22May19/rhs");
        rhs[0].ParallelCopy(tmp);
    }

    Real mlmg_err = mlmg.solve(amrex::GetVecOfPtrs(phi), amrex::GetVecOfConstPtrs(rhs),
                               1.e-11, 0.0);

#if 1
#if (AMREX_SPACEDIM == 2)
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(phi[ilev], "phi2d");
    }
#else
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(phi[ilev], "phi3d");
    }
#endif
#endif

//    mlndlap.compRHS(amrex::GetVecOfPtrs(rhs), amrex::GetVecOfPtrs(vel), {}, {});

#if 0
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(rhs[ilev], "rhs"+std::to_string(ilev));
        amrex::Print() << "rhs.norm0() = " << rhs[ilev].norm0() << "\n";
        amrex::Print() << "rhs.norm1()/npoints = " << rhs[ilev].norm1() / grids[0].d_numPts() << "\n";
    }
#endif
}

void
MyTest::writePlotfile ()
{
//    amrex::WriteSingleLevelPlotfile("plot", vel[0], {AMREX_D_DECL("xvel","yvel","zvel")}, geom[0], 0.0, 0);
}

void
MyTest::readParameters ()
{   
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("n_cell_x", n_cell_x);
    pp.query("n_cell_y", n_cell_y);
    pp.query("n_cell_z", n_cell_z);
    pp.query("max_grid_size_x", max_grid_size_x);
    pp.query("max_grid_size_y", max_grid_size_y);
    pp.query("max_grid_size_z", max_grid_size_z);

    pp.query("plot_file", plot_file_name);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("max_coarsening_level", max_coarsening_level);
#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
#endif

    pp.query("sigma", sigma);

#if (AMREX_SPACEDIM == 3)
    ParmParse pp_eb("eb2");
    std::string geom_type;
    pp_eb.get("cylinder_direction", cylinder_direction);
#endif
}

void
MyTest::initGrids ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(0.024,0.006,0.006)});

    // Make the domain periodic at the ends of the cylinder
    if (cylinder_direction == 0)
    {
       std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
       Geometry::Setup(&rb, 0, is_periodic.data());

    } else if (cylinder_direction == 1)
    {
       std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
       Geometry::Setup(&rb, 0, is_periodic.data());

    } else if (cylinder_direction == 2)
    {
       std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
       Geometry::Setup(&rb, 0, is_periodic.data());
    }

    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain);
        domain.refine(ref_ratio);
    }

    std::cout << "GEOM " << geom[0] << std::endl;

    domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        grids[ilev].define(domain);
        grids[ilev].maxSize(IntVect(max_grid_size_x,max_grid_size_y,max_grid_size_z));

        domain.grow(0,-n_cell_x/4);   // fine level cover the middle of the coarse domain
        domain.grow(1,-n_cell_y/4);   // fine level cover the middle of the coarse domain
#if (AMREX_SPACEDIM == 3)
        domain.grow(2,-n_cell_z/4);   // fine level cover the middle of the coarse domain
#endif

        domain.refine(ref_ratio); 
    }
}

void
MyTest::initData ()
{
    int nlevels = max_level + 1;
    dmap.resize(nlevels);
    factory.resize(nlevels);
    phi.resize(nlevels);
    rhs.resize(nlevels);
    sig.resize(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[ilev]);
        factory[ilev].reset(new EBFArrayBoxFactory(eb_level, geom[ilev], grids[ilev], dmap[ilev],
                                                   {2,2,2}, EBSupport::full));

        phi[ilev].define(amrex::convert(grids[ilev],IntVect::TheNodeVector()),
                         dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        rhs[ilev].define(amrex::convert(grids[ilev],IntVect::TheNodeVector()),
                         dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        sig[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
    }
}
