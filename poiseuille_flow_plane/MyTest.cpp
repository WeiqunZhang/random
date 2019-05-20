#include "MyTest.H"

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MLEBTensorOp.H>
#include <AMReX_MLMG.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_VisMF.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.0,1.0,1.0)});

    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    is_periodic[other_dir] = 1;
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    geom.define(domain);

    {
        bool static first = true;
        if (first) {
            first = false;
            RealArray point{0.,0.,0.};
            const Real dx = geom.CellSize(wall_dir);
            point[wall_dir] = dx_wall;
            RealArray normal{0.,0.,0.};
            normal[wall_dir] = -1.0;
            EB2::PlaneIF pf(point, normal);
            point[wall_dir] = 1.0 - dx_wall;
            normal[wall_dir] = 1.0;
            EB2::PlaneIF pf2(point, normal);
            auto gshop = EB2::makeShop(EB2::makeUnion(pf, pf2));
            EB2::Build(gshop, geom, 0, 100, 2);
        }
    }

    initData();
}

void
MyTest::solve ()
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    const Real tol_rel = 1.e-11;
    const Real tol_abs = 0.0;

    MLEBTensorOp ebtensorop({geom}, {grids}, {dmap}, info, {factory.get()});

    ebtensorop.setMaxOrder(linop_maxorder);

    Array<LinOpBCType,AMREX_SPACEDIM> v_lo_bc{AMREX_D_DECL(LinOpBCType::Neumann,
                                                           LinOpBCType::Neumann,
                                                           LinOpBCType::Neumann)};
    Array<LinOpBCType,AMREX_SPACEDIM> v_hi_bc{AMREX_D_DECL(LinOpBCType::Neumann,
                                                           LinOpBCType::Neumann,
                                                           LinOpBCType::Neumann)};

    v_lo_bc[flow_dir] = LinOpBCType::Dirichlet;
    v_lo_bc[other_dir] = LinOpBCType::Periodic;
    v_hi_bc[other_dir] = LinOpBCType::Periodic;

    ebtensorop.setDomainBC({AMREX_D_DECL(v_lo_bc,v_lo_bc,v_lo_bc)},
                           {AMREX_D_DECL(v_hi_bc,v_hi_bc,v_hi_bc)});

    ebtensorop.setLevelBC(0, &solution);

    const Real a = 0.0; // 1.0e6;
    {
        MultiFab tmp(grids, dmap, 1, 0, MFInfo(), *factory);
        tmp.setVal(a);

        ebtensorop.setACoeffs(0, tmp);

        Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(grids, IntVect::TheDimensionVector(idim));
            face_bcoef[idim].define(ba, dmap, 1, 0);
            face_bcoef[idim].setVal(mu);
        }
        ebtensorop.setShearViscosity(0, amrex::GetArrOfConstPtrs(face_bcoef));
        tmp.setVal(mu);
        ebtensorop.setEBShearViscosity(0, tmp);
    }

    MLMG mlmg(ebtensorop);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    mlmg.setBottomTolerance(1.e-4);

    MultiFab rhs(grids, dmap, 3, 0, MFInfo(), *factory);
    rhs.setVal(0.0);
    rhs.setVal(dpdx, flow_dir, 1);
    MultiFab::Saxpy(rhs, a, exact, 0, 0, 1, 0);

//    solution.setVal(0.0);
    mlmg.solve({&solution}, {&rhs}, tol_rel, tol_abs);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        amrex::Print() << "Component " << idim << " min and max "
                       << solution.min(idim) << ", " << solution.max(idim)
                       << "\n";
    }

    MultiFab error(grids, dmap, 1, 0, MFInfo(), *factory);
    MultiFab::Copy(error, solution, flow_dir, 0, 1, 0);
    MultiFab::Subtract(error, exact, 0, 0, 1, 0);
    const MultiFab& vfrc = factory->getVolFrac();
    MultiFab::Multiply(error, vfrc, 0, 0, 1, 0);
    const auto dx = geom.CellSize();
    error.mult(dx[0]*dx[1]*dx[2]);
    amrex::Print() << "1-norm error = " << error.norm1() << std::endl;
}

void
MyTest::readParameters ()
{
    ParmParse pp;
 
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("max_coarsening_level", max_coarsening_level);

    pp.query("mu", mu);
    pp.query("dpdx", dpdx);

    pp.query("flow_dir", flow_dir);
    pp.query("wall_dir", wall_dir);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(flow_dir != wall_dir, "check flow and wall dir");
    other_dir = 3 - flow_dir - wall_dir;

    pp.query("dx_wall", dx_wall);
}

void
MyTest::initData ()
{
    grids.define(geom.Domain());
    grids.maxSize(max_grid_size);
    dmap.define(grids);

    factory = makeEBFabFactory(geom, grids, dmap, {2,2,2}, EBSupport::full);

    solution.define(grids, dmap, 3, 1, MFInfo(), *factory);
    exact.define(grids, dmap, 1, 1, MFInfo(), *factory);

    const auto& dx = geom.CellSizeArray();
    const Real h = 1.0 - 2.*dx_wall;

    int ijk[3];
    for (MFIter mfi(exact); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.fabbox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        const Array4<Real> fab = exact.array(mfi);
        for         (int k = lo.z; k <= hi.z; ++k) {
            ijk[2] = k;
            for     (int j = lo.y; j <= hi.y; ++j) {
                ijk[1] = j;
                for (int i = lo.x; i <= hi.x; ++i) {
                    ijk[0] = i;
                    Real d = (ijk[wall_dir]+0.5)*dx[wall_dir] - dx_wall;
                    fab(i,j,k) = (0.5/mu)*dpdx*d*(h-d);
                }
            }
        }
    }

    amrex::EB_set_covered(exact, 0.0);
    solution.setVal(0.0);
    MultiFab::Copy(solution, exact, 0, flow_dir, 1, 1);
}
