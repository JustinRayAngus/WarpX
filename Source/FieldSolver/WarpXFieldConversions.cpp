/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "Fields.H"
#include "Utils/TextMsg.H"

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;
using warpx::fields::FieldType;

#ifdef WARPX_DIM_RTH
namespace
{
    // PICNIC's covariant edge field is (Er, r*Etheta, Ez).
    // Radius is evaluated at each component's own staggering.
    template <bool to_covariant>
    void ConvertPolarE (
        ablastr::fields::VectorField const& dst,
        ablastr::fields::VectorField const& src,
        amrex::Geometry const& geom)
    {
        const amrex::Real dr = geom.CellSize(0);
        const amrex::Real rmin = geom.ProbLo(0);
        const int ilo = geom.Domain().smallEnd(0);
        for (int dir = 0; dir < 3; ++dir) {
            if (dir != 1) {
                amrex::MultiFab::Copy(*dst[dir], *src[dir], 0, 0,
                                     src[dir]->nComp(), src[dir]->nGrowVect());
                continue;
            }
            const amrex::Real offset = src[dir]->ixType().nodeCentered(0) ? 0.0_rt : 0.5_rt;
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                rmin + (offset - src[dir]->nGrowVect()[0])*dr > 0.0_rt,
                "RTH E conversion requires positive radii throughout the radial guard cells.");
            // Process each FAB once: its full allocation includes physical-boundary guards.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(*dst[dir]); mfi.isValid(); ++mfi) {
                const auto in = src[dir]->const_array(mfi);
                const auto out = dst[dir]->array(mfi);
                amrex::ParallelFor(mfi.fabbox(), dst[dir]->nComp(),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                        const amrex::Real r = rmin + (i - ilo + offset)*dr;
                        if constexpr (to_covariant) {
                            out(i,j,k,n) = r*in(i,j,k,n);
                        } else {
                            out(i,j,k,n) = in(i,j,k,n)/r;
                        }
                    });
            }
        }
    }
}

void WarpX::ConvertPhysicalEToCovariant (int lev, PatchType patch_type)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0 && patch_type == PatchType::fine,
        "RTH E conversion currently supports only the level-zero fine patch.");
    ConvertPolarE<true>(
        m_fields.get_alldirs(FieldType::Efield_cov_fp, lev),
        m_fields.get_alldirs(FieldType::Efield_fp, lev), Geom(lev));
}

void WarpX::ConvertCovariantEToPhysical (int lev, PatchType patch_type)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0 && patch_type == PatchType::fine,
        "RTH E conversion currently supports only the level-zero fine patch.");
    ConvertPolarE<false>(
        m_fields.get_alldirs(FieldType::Efield_fp, lev),
        m_fields.get_alldirs(FieldType::Efield_cov_fp, lev), Geom(lev));
}
#endif
