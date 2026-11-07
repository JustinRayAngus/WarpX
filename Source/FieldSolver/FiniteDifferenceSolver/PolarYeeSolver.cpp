/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FiniteDifferenceSolver.H"

#ifdef WARPX_DIM_RTH
#include "EmbeddedBoundary/Enabled.H"
#include "FiniteDifferenceAlgorithms/PolarYeeAlgorithm.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_Geometry.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>

using namespace amrex::literals;

namespace
{
    PolarYeeAlgorithm GetPolarAlgorithm (int lev = 0)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0,
            "RTH Yee operators currently support only level zero.");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!EB::enabled(),
            "RTH Yee operators do not support embedded boundaries.");
        auto const& geom = WarpX::GetInstance().Geom(lev);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(geom.ProbLo(0) > 0.0_rt && !geom.isPeriodic(0),
            "RTH Yee operators require an annulus with nonperiodic radial boundaries.");
        return {geom.CellSize(0), 1.0_rt/geom.CellSize(0), 1.0_rt/geom.CellSize(1),
                geom.ProbLo(0), geom.Domain().smallEnd(0)};
    }

    // edge=true: Er(0,1), Etheta(1,0), Ez(1,1); B has the complementary staggering.
    void CheckPolarVector (ablastr::fields::VectorField const& field, bool edge,
                           amrex::IntVect const& ngrow)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(ngrow.allGE(amrex::IntVect(0)),
            "RTH Yee operators require nonnegative guard extents.");
        for (int dir = 0; dir < 3; ++dir) {
            amrex::IntVect const nodes{int(edge != (dir == 0)), int(edge != (dir == 1))};
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                field[dir]->ixType().toIntVect() == nodes && field[dir]->nComp() == 1,
                "RTH Yee operators require single-component fields at the RTH staggering.");
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(field[dir]->nGrowVect().allGE(ngrow),
                "RTH Yee operator has insufficient field guard cells.");
        }
    }
}

// `upward` maps physical edge fields to faces; otherwise physical faces to edges.
template <bool upward>
void FiniteDifferenceSolver::ComputeCurlPolar (
    ablastr::fields::VectorField const& out,
    ablastr::fields::VectorField const& in,
    int lev, amrex::IntVect const& ngrow)
{
    auto const algo = GetPolarAlgorithm(lev);
    CheckPolarVector(out, !upward, ngrow);
    CheckPolarVector(in, upward, ngrow + 1);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(algo.m_rmin - ngrow[0]*algo.m_dr > 0.0_rt,
        "RTH curl requires positive radii throughout the output guard cells.");

    for (int dir = 0; dir < 3; ++dir) { out[dir]->setVal(0.0_rt); }
    auto* cost = WarpX::getCosts(lev);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*out[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers) {
            amrex::Gpu::synchronize();
        }
        auto const start = amrex::second();
        auto const fr = in[0]->const_array(mfi);
        auto const ft = in[1]->const_array(mfi);
        auto const fz = in[2]->const_array(mfi);
        auto const cr = out[0]->array(mfi);
        auto const ct = out[1]->array(mfi);
        auto const cz = out[2]->array(mfi);
        // Only boundary tiles grow, so ghost-cell writes do not overlap.
        auto const br = mfi.tilebox(out[0]->ixType().toIntVect(), ngrow);
        auto const bt = mfi.tilebox(out[1]->ixType().toIntVect(), ngrow);
        auto const bz = mfi.tilebox(out[2]->ixType().toIntVect(), ngrow);
        amrex::ParallelFor(br, bt, bz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                cr(i,j,k) = algo.CurlR<upward>(fz, i, j, k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                ct(i,j,k) = algo.CurlTheta<upward>(fz, i, j, k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                cz(i,j,k) = algo.CurlZ<upward>(fr, ft, i, j, k);
            });
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers) {
            amrex::Gpu::synchronize();
            auto const elapsed = static_cast<amrex::Real>(amrex::second() - start);
            amrex::HostDevice::Atomic::Add(&(*cost)[mfi.index()], elapsed);
        }
    }
}

template void FiniteDifferenceSolver::ComputeCurlPolar<true>(
    ablastr::fields::VectorField const&, ablastr::fields::VectorField const&,
    int, amrex::IntVect const&);
template void FiniteDifferenceSolver::ComputeCurlPolar<false>(
    ablastr::fields::VectorField const&, ablastr::fields::VectorField const&,
    int, amrex::IntVect const&);

void FiniteDifferenceSolver::EvolveBPolar (
    ablastr::fields::VectorField const& Bfield,
    ablastr::fields::VectorField const& Efield,
    amrex::MultiFab const* Gfield, int lev, amrex::Real dt)
{
    auto const algo = GetPolarAlgorithm(lev);
    CheckPolarVector(Bfield, false, amrex::IntVect(0));
    CheckPolarVector(Efield, true, amrex::IntVect(1));
    if (Gfield) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            Gfield->ixType().cellCentered() && Gfield->nGrowVect().allGE(amrex::IntVect(1)),
            "RTH magnetic cleaning requires cell-centered G with one guard cell.");
    }
    auto* cost = WarpX::getCosts(lev);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Bfield[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers) {
            amrex::Gpu::synchronize();
        }
        auto const start = amrex::second();
        auto const er = Efield[0]->const_array(mfi);
        auto const et_cov = Efield[1]->const_array(mfi);
        auto const ez = Efield[2]->const_array(mfi);
        auto const br = Bfield[0]->array(mfi);
        auto const bt = Bfield[1]->array(mfi);
        auto const bz = Bfield[2]->array(mfi);
        auto const tr = mfi.tilebox(Bfield[0]->ixType().toIntVect());
        auto const tt = mfi.tilebox(Bfield[1]->ixType().toIntVect());
        auto const tz = mfi.tilebox(Bfield[2]->ixType().toIntVect());
        amrex::ParallelFor(tr, tt, tz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                br(i,j,k) -= dt*algo.CurlR<true>(ez, i, j, k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                bt(i,j,k) -= dt*algo.CurlTheta<true>(ez, i, j, k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                bz(i,j,k) -= dt*algo.CurlZ<true, true>(er, et_cov, i, j, k);
            });
        if (Gfield) {
            auto const g = Gfield->const_array(mfi);
            amrex::ParallelFor(tr, tt,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    br(i,j,k) += dt*algo.Dr<false>(g, i, j, k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    bt(i,j,k) += dt*algo.Dtheta<false>(g, i, j, k)/algo.Radius(i, false);
                });
        }
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers) {
            amrex::Gpu::synchronize();
            auto const elapsed = static_cast<amrex::Real>(amrex::second() - start);
            amrex::HostDevice::Atomic::Add(&(*cost)[mfi.index()], elapsed);
        }
    }
}

void FiniteDifferenceSolver::EvolveEPolar (
    ablastr::fields::VectorField const& Efield,
    ablastr::fields::VectorField const& Bfield,
    ablastr::fields::VectorField const& Jfield,
    amrex::MultiFab const* Ffield, int lev, amrex::Real dt)
{
    auto const algo = GetPolarAlgorithm(lev);
    CheckPolarVector(Efield, true, amrex::IntVect(0));
    CheckPolarVector(Bfield, false, amrex::IntVect(1));
    CheckPolarVector(Jfield, true, amrex::IntVect(0));
    if (Ffield) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            Ffield->ixType().nodeCentered() && Ffield->nGrowVect().allGE(amrex::IntVect(1)),
            "RTH electric cleaning requires nodal F with one guard cell.");
    }
    auto* cost = WarpX::getCosts(lev);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Efield[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers) {
            amrex::Gpu::synchronize();
        }
        auto const start = amrex::second();
        auto const br = Bfield[0]->const_array(mfi);
        auto const bt = Bfield[1]->const_array(mfi);
        auto const bz = Bfield[2]->const_array(mfi);
        auto const jr = Jfield[0]->const_array(mfi);
        auto const jt = Jfield[1]->const_array(mfi);
        auto const jz = Jfield[2]->const_array(mfi);
        auto const er = Efield[0]->array(mfi);
        auto const et = Efield[1]->array(mfi);
        auto const ez = Efield[2]->array(mfi);
        auto const tr = mfi.tilebox(Efield[0]->ixType().toIntVect());
        auto const tt = mfi.tilebox(Efield[1]->ixType().toIntVect());
        auto const tz = mfi.tilebox(Efield[2]->ixType().toIntVect());
        amrex::ParallelFor(tr, tt, tz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                er(i,j,k) += PhysConst::c2*dt*(algo.CurlR<false>(bz, i, j, k)
                                             - PhysConst::mu0*jr(i,j,k));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                et(i,j,k) += PhysConst::c2*dt*(algo.CurlTheta<false>(bz, i, j, k)
                                             - PhysConst::mu0*jt(i,j,k));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                ez(i,j,k) += PhysConst::c2*dt*(algo.CurlZ<false>(br, bt, i, j, k)
                                             - PhysConst::mu0*jz(i,j,k));
            });
        if (Ffield) {
            auto const f = Ffield->const_array(mfi);
            amrex::ParallelFor(tr, tt,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    er(i,j,k) += PhysConst::c2*dt*algo.Dr<true>(f, i, j, k);
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    et(i,j,k) += PhysConst::c2*dt*algo.Dtheta<true>(f, i, j, k)
                        / algo.Radius(i, true);
                });
        }
        if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers) {
            amrex::Gpu::synchronize();
            auto const elapsed = static_cast<amrex::Real>(amrex::second() - start);
            amrex::HostDevice::Atomic::Add(&(*cost)[mfi.index()], elapsed);
        }
    }
}

void FiniteDifferenceSolver::ComputeDivEPolar (
    ablastr::fields::VectorField const& Efield, amrex::MultiFab& divEfield)
{
    auto const algo = GetPolarAlgorithm();
    CheckPolarVector(Efield, true, amrex::IntVect(1));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(divEfield.ixType().nodeCentered(),
        "RTH div E must be nodal.");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(divEfield, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const er = Efield[0]->const_array(mfi);
        auto const et = Efield[1]->const_array(mfi);
        auto const div = divEfield.array(mfi);
        amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                div(i,j,k) = algo.Div<false>(er, et, i, j, k);
            });
    }
}

void FiniteDifferenceSolver::EvolveFPolar (
    amrex::MultiFab* Ffield, ablastr::fields::VectorField const& Efield,
    amrex::MultiFab const* rhofield, int rho_comp, amrex::Real dt)
{
    auto const algo = GetPolarAlgorithm();
    CheckPolarVector(Efield, true, amrex::IntVect(1));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        Ffield->ixType().nodeCentered() && rhofield->ixType().nodeCentered() &&
        rho_comp >= 0 && rho_comp < rhofield->nComp(),
        "RTH electric cleaning requires nodal F and rho with a valid charge component.");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Ffield, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const er = Efield[0]->const_array(mfi);
        auto const et = Efield[1]->const_array(mfi);
        auto const rho = rhofield->const_array(mfi);
        auto const f = Ffield->array(mfi);
        amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                f(i,j,k) += dt*(algo.Div<false>(er, et, i, j, k)
                                - rho(i,j,k,rho_comp)/PhysConst::epsilon_0);
            });
    }
}

void FiniteDifferenceSolver::EvolveGPolar (
    amrex::MultiFab* Gfield, ablastr::fields::VectorField const& Bfield, amrex::Real dt)
{
    auto const algo = GetPolarAlgorithm();
    CheckPolarVector(Bfield, false, amrex::IntVect(1));
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(Gfield->ixType().cellCentered(),
        "RTH magnetic cleaning requires cell-centered G.");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*Gfield, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const br = Bfield[0]->const_array(mfi);
        auto const bt = Bfield[1]->const_array(mfi);
        auto const g = Gfield->array(mfi);
        amrex::ParallelFor(mfi.tilebox(),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                g(i,j,k) += PhysConst::c2*dt*algo.Div<true>(br, bt, i, j, k);
            });
    }
}
#endif
