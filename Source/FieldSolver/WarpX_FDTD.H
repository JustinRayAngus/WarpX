/* Copyright 2019-2025 The WarpX Community
 * Axel Huebl, David Grote, Eric Clark
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FDTD_H_
#define WARPX_FDTD_H_

#ifdef WARPX_DIM_RTH
#include "FiniteDifferenceSolver/FiniteDifferenceAlgorithms/PolarYeeAlgorithm.H"
#endif

#include <AMReX.H>
#include <AMReX_FArrayBox.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void warpx_computedivb(
    int i, int j, int k, int dcomp,
    amrex::Array4<amrex::Real> const& divB,
    amrex::Array4<amrex::Real const> const& Bx,
    amrex::Array4<amrex::Real const> const& By,
    amrex::Array4<amrex::Real const> const& Bz,
    amrex::Real dxinv,
    amrex::Real dyinv,
    amrex::Real dzinv,
#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
    amrex::Real const rmin,
#elif defined(WARPX_DIM_RTH)
    amrex::Real const rmin, int const radial_domain_lo,
#endif
    bool collocated = false)
{
    // TODO: Update to use the templated operators in the FiniteDifference object
    using namespace amrex;

#if defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
    const amrex::Real ru = 1._rt + 0.5_rt/(rmin*dxinv + i + 0.5_rt);
    const amrex::Real rd = 1._rt - 0.5_rt/(rmin*dxinv + i + 0.5_rt);
#endif

    if (collocated) {
        //Stencil is across 2 cells.
#if defined WARPX_DIM_3D
        divB(i,j,k,dcomp) = 0.5_rt*(Bx(i+1,j  ,k  ) - Bx(i-1,j,k))*dxinv
                +        0.5_rt*(By(i  ,j+1,k  ) - By(i,j-1,k))*dyinv
                +        0.5_rt*(Bz(i  ,j  ,k+1) - Bz(i,j,k-1))*dzinv;
#elif defined WARPX_DIM_XZ
        divB(i,j,0,dcomp) = 0.5_rt*(Bx(i+1,j  ,0) - Bx(i-1,j,0))*dxinv
                +        0.5_rt*(Bz(i  ,j+1,0) - Bz(i,j-1,0))*dzinv;
        amrex::ignore_unused(k, By, dyinv);
#elif defined WARPX_DIM_1D_Z
        divB(i,0,0,dcomp) = 0.5_rt*(Bz(i+1,0  ,0) - Bz(i-1,0,0))*dzinv;
        amrex::ignore_unused(j, Bx, dxinv);
        amrex::ignore_unused(k, By, dyinv);
#elif defined(WARPX_DIM_RZ) || defined(WARPX_DIM_RCYLINDER) || defined(WARPX_DIM_RSPHERE)
        amrex::Abort("Collocated non-Cartesian grids are not a currently supported mode for computing div B.");
#elif defined(WARPX_DIM_RTH)
        amrex::Abort("RTH div B requires a Yee staggered grid.");
#endif
    } else {
#if defined WARPX_DIM_3D
        divB(i,j,k,dcomp) = (Bx(i+1,j  ,k  ) - Bx(i,j,k))*dxinv
                +        (By(i  ,j+1,k  ) - By(i,j,k))*dyinv
                +        (Bz(i  ,j  ,k+1) - Bz(i,j,k))*dzinv;
#elif defined WARPX_DIM_XZ
        divB(i,j,0,dcomp) = (Bx(i+1,j  ,0) - Bx(i,j,0))*dxinv
                +        (Bz(i  ,j+1,0) - Bz(i,j,0))*dzinv;
        amrex::ignore_unused(k, By, dyinv);
#elif defined WARPX_DIM_1D_Z
        divB(i,0,0,dcomp) = (Bz(i+1,0  ,0) - Bz(i,0,0))*dzinv;
        amrex::ignore_unused(j, Bx, dxinv);
        amrex::ignore_unused(k, By, dyinv);
#elif defined(WARPX_DIM_RCYLINDER)
        divB(i,0,0,dcomp) = (ru*Bx(i+1,0  ,0) - rd*Bx(i,0,0))*dxinv;
        amrex::ignore_unused(j, By, dyinv);
        amrex::ignore_unused(k, Bz, dzinv);
#elif defined(WARPX_DIM_RSPHERE)
        divB(i,0,0,dcomp) = (ru*ru*Bx(i+1,0  ,0) - rd*rd*Bx(i,0,0))*dxinv;
        amrex::ignore_unused(j, By, dyinv);
        amrex::ignore_unused(k, Bz, dzinv);
#elif defined WARPX_DIM_RZ
        divB(i,j,0,dcomp) = (ru*Bx(i+1,j,0) - rd*Bx(i,j,0))*dxinv
                           + (Bz(i,j+1,0) - Bz(i,j,0))*dzinv;
        amrex::ignore_unused(k, By, dyinv);
#elif defined(WARPX_DIM_RTH)
        PolarYeeAlgorithm const algo{1.0_rt/dxinv, dxinv, dyinv, rmin, radial_domain_lo};
        divB(i,j,k,dcomp) = algo.Div<true>(Bx, By, i, j, k);
        amrex::ignore_unused(Bz, dzinv);
#endif
    }
}

#endif
