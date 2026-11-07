/* Copyright 2026 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FINITE_DIFFERENCE_ALGORITHM_POLAR_YEE_H_
#define WARPX_FINITE_DIFFERENCE_ALGORITHM_POLAR_YEE_H_

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

/** Yee differences on an annular (r, theta) mesh, with vector order (r, theta, z).
 * `upward` differences go from nodes to cells; downward differences go from
 * cells to nodes. All results are physical components. Inputs are physical
 * except when CurlZ is explicitly called with a covariant angular component.
 */
struct PolarYeeAlgorithm
{
    amrex::Real m_dr;
    amrex::Real m_inv_dr;
    amrex::Real m_inv_dtheta;
    amrex::Real m_rmin;
    int m_ilo;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real Radius (int i, bool nodal) const
    {
        using namespace amrex::literals;
        return m_rmin + (i - m_ilo + (nodal ? 0.0_rt : 0.5_rt))*m_dr;
    }

    template <bool upward>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real Dr (amrex::Array4<amrex::Real const> const& f,
                    int i, int j, int k) const
    {
        constexpr int hi = upward ? 1 : 0;
        constexpr int lo = upward ? 0 : -1;
        return m_inv_dr*(f(i+hi,j,k) - f(i+lo,j,k));
    }

    template <bool upward>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real Dtheta (amrex::Array4<amrex::Real const> const& f,
                        int i, int j, int k) const
    {
        constexpr int hi = upward ? 1 : 0;
        constexpr int lo = upward ? 0 : -1;
        return m_inv_dtheta*(f(i,j+hi,k) - f(i,j+lo,k));
    }

    // d_r(r*f): evaluate each input radius at its own staggering.
    template <bool upward>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real Drr (amrex::Array4<amrex::Real const> const& f,
                     int i, int j, int k) const
    {
        constexpr int hi = upward ? 1 : 0;
        constexpr int lo = upward ? 0 : -1;
        return m_inv_dr*(Radius(i+hi, upward)*f(i+hi,j,k)
                      - Radius(i+lo, upward)*f(i+lo,j,k));
    }

    // (curl f)_r = d_theta(f_z)/r; radial output is nodal for curl E.
    template <bool upward>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real CurlR (amrex::Array4<amrex::Real const> const& fz,
                       int i, int j, int k) const
    {
        return Dtheta<upward>(fz, i, j, k)/Radius(i, upward);
    }

    // (curl f)_theta = -d_r(f_z).
    template <bool upward>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real CurlTheta (amrex::Array4<amrex::Real const> const& fz,
                           int i, int j, int k) const
    {
        return -Dr<upward>(fz, i, j, k);
    }

    // (curl f)_z = [d_r(r*f_theta) - d_theta(f_r)]/r.
    // EvolveB already receives r*Etheta, while standalone curls take physical fields.
    template <bool upward, bool covariant = false>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real CurlZ (amrex::Array4<amrex::Real const> const& fr,
                       amrex::Array4<amrex::Real const> const& ft,
                       int i, int j, int k) const
    {
        amrex::Real radial;
        if constexpr (covariant) {
            radial = Dr<upward>(ft, i, j, k);
        } else {
            radial = Drr<upward>(ft, i, j, k);
        }
        return (radial - Dtheta<upward>(fr, i, j, k))/Radius(i, !upward);
    }

    // div f = [d_r(r*f_r) + d_theta(f_theta)]/r.
    // Upward gives cell-centered div B; downward gives nodal div E.
    template <bool upward>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real Div (amrex::Array4<amrex::Real const> const& fr,
                     amrex::Array4<amrex::Real const> const& ft,
                     int i, int j, int k) const
    {
        return (Drr<upward>(fr, i, j, k) + Dtheta<upward>(ft, i, j, k))
            / Radius(i, !upward);
    }
};

#endif
