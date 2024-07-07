/* Copyright 2019 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_
#define WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_

#include "Utils/WarpXConst.H"

#include <AMReX_Math.H>
#include <AMReX_Random.H>

#include <cmath>  // isnan() isinf()
#include <limits> // numeric_limits<float>::min()

/* \brief Update particle velocities according to
 *        F. Perez et al., Phys.Plasmas.19.083104 (2012),
 *        which is based on Nanbu's method, PhysRevE.55.4642 (1997).
 *        @param[in] bmax is max(Debye length, minimal interparticle distance).
 *        @param[in] L is the Coulomb log. A fixed L will be used if L > 0,
 *        otherwise L will be calculated based on the algorithm.
 *        @param[in] n12 = max(w1,w2)*min(N1,N2)/dV is the effective density used for s12
 *        @param[in] sigma_max is the maximum cross section based on mfp = atomic spacing
 *        used for the normalized scattering length s12
 *        To see if there are nan or inf updated velocities,
 *        compile with USE_ASSERTION=TRUE.
 *
 * Updates and corrections to the original publication are documented in
 * https://github.com/ECP-WarpX/WarpX/issues/429
 * https://github.com/ECP-WarpX/WarpX/files/3799803/main.pdf
 */

template <typename T_PR, typename T_R>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void UpdateMomentumPerezElastic (
    T_PR& u1x, T_PR& u1y, T_PR& u1z, T_PR& u2x, T_PR& u2y, T_PR& u2z,
    T_PR const n12, T_PR const sigma_max,
    T_PR const q1, T_PR const m1, T_PR const w1,
    T_PR const q2, T_PR const m2, T_PR const w2,
    T_R const dt, T_PR const L, T_PR const bmax,
    amrex::RandomEngine const& engine)
{

    T_PR const diffx = amrex::Math::abs(u1x-u2x);
    T_PR const diffy = amrex::Math::abs(u1y-u2y);
    T_PR const diffz = amrex::Math::abs(u1z-u2z);
    T_PR const diffm = std::sqrt(diffx*diffx+diffy*diffy+diffz*diffz);
    T_PR const summm = std::sqrt(u1x*u1x+u1y*u1y+u1z*u1z) + std::sqrt(u2x*u2x+u2y*u2y+u2z*u2z);
    // If g = u1 - u2 = 0, do not collide.
    // Or if the relative difference is less than 1.0e-10.
    if ( diffm < std::numeric_limits<T_PR>::min() || diffm/summm < 1.0e-10 ) { return; }

    T_PR constexpr inv_c2 = T_PR(1.0) / ( PhysConst::c * PhysConst::c );

    // Compute Lorentz factor gamma
    T_PR const g1 = std::sqrt( T_PR(1.0) + (u1x*u1x+u1y*u1y+u1z*u1z)*inv_c2 );
    T_PR const g2 = std::sqrt( T_PR(1.0) + (u2x*u2x+u2y*u2y+u2z*u2z)*inv_c2 );

    // Compute momenta
    T_PR const p1x = u1x * m1;
    T_PR const p1y = u1y * m1;
    T_PR const p1z = u1z * m1;
    T_PR const ptotx = p1x + u2x * m2;
    T_PR const ptoty = p1y + u2y * m2;
    T_PR const ptotz = p1z + u2z * m2;

    // Compute center-of-mass (COM) velocity and gamma
    T_PR const mass_g = (g1 * m1 + g2 * m2); // Etot/c^2
    T_PR const vcmx   = ptotx / mass_g;
    T_PR const vcmy   = ptoty / mass_g;
    T_PR const vcmz   = ptotz / mass_g;
    T_PR const vcmsq  = vcmx*vcmx + vcmy*vcmy + vcmz*vcmz;
    T_PR const gcm    = T_PR(1.0) / std::sqrt(T_PR(1.0) - vcmsq*inv_c2);

    // Compute gbetacm dot gbeta1 and gbeta2
    T_PR const gbcDgb1 = gcm*(vcmx*u1x + vcmy*u1y + vcmz*u1z)*inv_c2;
    T_PR const gbcDgb2 = gcm*(vcmx*u2x + vcmy*u2y + vcmz*u2z)*inv_c2;

    // Compute p1* (i.e., transform p1 from lab to COM frame)
    T_PR const lorentz_transform_factor = gcm*(gbcDgb1/(T_PR(1.0) + gcm) - g1)*m1;
    T_PR const p1stx = p1x + vcmx*lorentz_transform_factor;
    T_PR const p1sty = p1y + vcmy*lorentz_transform_factor;
    T_PR const p1stz = p1z + vcmz*lorentz_transform_factor;
    T_PR const p1stm = std::sqrt( p1stx*p1stx + p1sty*p1sty + p1stz*p1stz );

    // Compute gamma star for particles 1 and 2
    T_PR const g1st = gcm*g1 - gbcDgb1;
    T_PR const g2st = gcm*g2 - gbcDgb2;
    
    // Compute variant relative velocity in cm frame
    T_PR const muRst = g1st*m1*g2st*m2/(g1st*m1 + g2st*m2);
    T_PR const vrelst = p1stm/muRst;
    
    // Compute invariant relative velocity (frame independent)
    T_PR const denom = T_PR(1.0) + p1stm*p1stm/(m1*g1st*m2*g2st)*inv_c2; // (1.0 - v1st*v2st/c^2)
    T_PR const vrelst_invar = vrelst/denom;

    // Compute s12
    T_PR s12 = 0;
    if (vrelst > std::numeric_limits<T_PR>::min()) {

        // s12 is non-zero (i.e. particles scatter) only if the relative
        // motion between particles is not negligible
            
        // Writing b0 = b90 in a form that is directly analagous to the well-known
        // non-relativistic form. Note that b0 on Eq. 22 of Perez POP 19 (2012) is bmin = b0/2,
        // Note: there is a typo in Eq 22 of Perez, the last square is incorrect!
        // See the SMILEI documentation: https://smileipic.github.io/Smilei/Understand/collisions.html
        // and https://github.com/ECP-WarpX/WarpX/files/3799803/main.pdf from GitHub #429
        T_PR const b0 = amrex::Math::abs(q1*q2) /
                        (T_PR(2.0)*MathConst::pi*PhysConst::ep0*muRst*vrelst*vrelst_invar);

        // Compute the Coulomb log lnLmd first
        T_PR lnLmd;
        if ( L > T_PR(0.0) ) { lnLmd = L; }
        else
        {
            // Compute the minimum impact parameter from quantum (bmin = lambda_DB/(4*pi))
            // Note that there is an error in Perez 2012 and in Lee and More 1984, where
            // they incorrectly use h instead of hbar. See NRL formulary. Also see "An
            // introduction to the physics of the Coulomb logarithm, with emphasis on
            // quantum-mechanical effects", J. Plasma Phys. vol. 85 (2019). by J.A. Krommes
            const T_PR bmin_qm = static_cast<T_PR>(PhysConst::hbar*0.5/p1stm);

            // Set the minimum impact parameter
            const T_PR bmin = amrex::max(bmin_qm, T_PR(0.5)*b0);

            // Compute the Coulomb log lnLmd
            lnLmd = amrex::max( T_PR(2.0),
                    T_PR(0.5)*std::log(T_PR(1.0) + bmax*bmax/(bmin*bmin)) );
        }

        // Compute s12 (with sigma limited by sigma_max where mfp = atomic spacing)
        const auto sigma_eff = amrex::min(MathConst::pi * b0 * b0 * lnLmd,sigma_max);
        s12 = sigma_eff * n12 * dt * vrelst * g1st*g2st/(g1*g2);

    }

    // Only modify momenta if s12 is non-zero
    if (s12 > std::numeric_limits<T_PR>::min()) {

        // Get random numbers
        T_PR r = amrex::Random(engine);

        // Compute scattering angle
        T_PR cosXs;
        T_PR sinXs;
        if ( s12 <= T_PR(0.1466) )
        {
            T_PR const Ainv = static_cast<T_PR>(s12*(1.0 - s12/2.0 + s12*s12/6.0));
            r = amrex::Random(engine);
            cosXs = T_PR(1.0) + Ainv*std::log(T_PR(1.0) - r*(T_PR(1.0)-std::exp(-T_PR(2.0)/Ainv)));
        }
        else if ( s12 > T_PR(0.1466) && s12 <= T_PR(3.0) )
        {
            T_PR const s12sq = s12*s12;
            T_PR const s12cu = s12*s12sq;
            T_PR const Ainv = static_cast<T_PR>(
                0.0056958 + 0.9560202*s12 - 0.508139*s12sq +
                0.47913906*s12cu - 0.12788975*s12sq*s12sq + 0.02389567*s12sq*s12cu);
            cosXs = T_PR(1.0) + Ainv*std::log(T_PR(1.0) - r*(T_PR(1.0)-std::exp(-T_PR(2.0)/Ainv)));
        }
        else if ( s12 > T_PR(3.0) && s12 <= T_PR(6.0) )
        {
            T_PR const Ainv = T_PR(1.0)/(T_PR(3.0) * std::exp(-s12));
            cosXs = T_PR(1.0) + Ainv*std::log(T_PR(1.0) - r*(T_PR(1.0)-std::exp(-T_PR(2.0)/Ainv)));
        }
        else
        {
            cosXs = T_PR(2.0) * r - T_PR(1.0);
        }
        sinXs = std::sqrt(T_PR(1.0) - cosXs*cosXs);

        // Get random azimuthal angle
        T_PR const phis = amrex::Random(engine) * T_PR(2.0) * MathConst::pi;
        T_PR const cosphis = std::cos(phis);
        T_PR const sinphis = std::sin(phis);

        // Compute post-collision momenta pfs in COM
        T_PR p1fstx;
        T_PR p1fsty;
        T_PR p1fstz;
        T_PR const p1stp = std::sqrt( p1stx*p1stx + p1sty*p1sty );
        if ( p1stp > std::numeric_limits<T_PR>::min() )
        {
            p1fstx = ( p1stx*p1stz/p1stp ) * sinXs*cosphis +
                     (-p1sty*p1stm/p1stp ) * sinXs*sinphis +
                     ( p1stx           ) * cosXs;
            p1fsty = ( p1sty*p1stz/p1stp ) * sinXs*cosphis +
                     ( p1stx*p1stm/p1stp ) * sinXs*sinphis +
                     ( p1sty           ) * cosXs;
            p1fstz = (-p1stp           ) * sinXs*cosphis +
                     ( p1stz           ) * cosXs;
        }
        else
        {
            // If p1stp is almost zero
            p1fstx = p1stm * sinXs*cosphis;
            p1fsty = p1stm * sinXs*sinphis;
            p1fstz = p1stm * cosXs;
        }

        // Transform p1fst from COM to lab frame
        T_PR const ucDp1fs = gcm*(vcmx*p1fstx + vcmy*p1fsty + vcmz*p1fstz);
        T_PR const factor1 = gcm*(ucDp1fs*inv_c2/(T_PR(1.0)+gcm) + m1*g1st);
        T_PR const p1fx = p1fstx + vcmx * factor1;
        T_PR const p1fy = p1fsty + vcmy * factor1;
        T_PR const p1fz = p1fstz + vcmz * factor1;

        // Set final momentum of 2 using conservation of momentum
        T_PR const p2fx = ptotx - p1fx;
        T_PR const p2fy = ptoty - p1fy;
        T_PR const p2fz = ptotz - p1fz;

        // Rejection method
        r = amrex::Random(engine);
        const auto wmax = amrex::max(w1, w2);
        if ( w2 > r*wmax)
        {
            u1x  = p1fx / m1;
            u1y  = p1fy / m1;
            u1z  = p1fz / m1;
        }
        r = amrex::Random(engine);
        if ( w1 > r*wmax )
        {
            u2x  = p2fx / m2;
            u2y  = p2fy / m2;
            u2z  = p2fz / m2;
        }
#ifndef AMREX_USE_DPCPP
        AMREX_ASSERT(!std::isnan(u1x+u1y+u1z+u2x+u2y+u2z));
        AMREX_ASSERT(!std::isinf(u1x+u1y+u1z+u2x+u2y+u2z));
#endif

    } // if s12 > std::numeric_limits<T_PR>::min()

}

#endif // WARPX_PARTICLES_COLLISION_UPDATE_MOMENTUM_PEREZ_ELASTIC_H_
