#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLMG.H>

class MyTest
{
public:

    MyTest ();

    void solve ();

private:

    void initializeEB ();
    void readParameters ();
    void initGrids ();
    void initData ();

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 2;
    int max_iter = 100;
    int max_fmg_iter = 0;
    int max_coarsening_level = 30;

    int cylinder_direction = 0;

    amrex::Geometry geom;
    amrex::BoxArray grids;
    amrex::DistributionMapping dmap;
    std::unique_ptr<amrex::EBFArrayBoxFactory>  factory;

    amrex::MultiFab phi;
    amrex::MultiFab rhs;
    amrex::MultiFab sig;
};

#endif
