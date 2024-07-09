/* Copyright 2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "CollisionHandler.H"

#include "Particles/Collision/BackgroundMCC/BackgroundMCCCollision.H"
#include "Particles/Collision/BackgroundStopping/BackgroundStopping.H"
#include "Particles/Collision/BinaryCollision/BinaryCollision.H"
#include "Particles/Collision/BinaryCollision/Coulomb/PairWiseCoulombCollisionFunc.H"
#include "Particles/Collision/BinaryCollision/DSMC/DSMCFunc.H"
#include "Particles/Collision/BinaryCollision/DSMC/SplitAndScatterFunc.H"
#include "Particles/Collision/BinaryCollision/NuclearFusion/NuclearFusionFunc.H"
#include "Particles/Collision/BinaryCollision/ParticleCreationFunc.H"
#include "Utils/TextMsg.H"

#include <AMReX_ParmParse.H>

#include <vector>

CollisionHandler::CollisionHandler(MultiParticleContainer const * const mypc)
{

    // Read in collision input
    const amrex::ParmParse pp_collisions("collisions");
    pp_collisions.queryarr("collision_names", collision_names);

    // Create instances based on the collision type
    auto const ncollisions = collision_names.size();
    collision_types.resize(ncollisions);
    allcollisions.resize(ncollisions);
    for (int i = 0; i < static_cast<int>(ncollisions); ++i) {
        const amrex::ParmParse pp_collision_name(collision_names[i]);

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(WarpX::n_rz_azimuthal_modes==1,
        "RZ mode `warpx.n_rz_azimuthal_modes` must be 1 when using the binary collision module.");

        // For legacy, pairwisecoulomb is the default
        std::string type = "pairwisecoulomb";

        pp_collision_name.query("type", type);
        collision_types[i] = type;

        if (type == "pairwisecoulomb") {
            allcollisions[i] =
               std::make_unique<BinaryCollision<PairWiseCoulombCollisionFunc>>(
                    collision_names[i], mypc
                );
            m_contains_coulomb_type = true;
        }
        else if (type == "background_mcc") {
            allcollisions[i] = std::make_unique<BackgroundMCCCollision>(collision_names[i]);
        }
        else if (type == "background_stopping") {
            allcollisions[i] = std::make_unique<BackgroundStopping>(collision_names[i]);
        }
        else if (type == "dsmc") {
            allcollisions[i] =
                std::make_unique<BinaryCollision<DSMCFunc, SplitAndScatterFunc>>(
                    collision_names[i], mypc
                );
        }
        else if (type == "nuclearfusion") {
            allcollisions[i] =
               std::make_unique<BinaryCollision<NuclearFusionFunc, ParticleCreationFunc>>(
                    collision_names[i], mypc
                );
        }
        else{
            WARPX_ABORT_WITH_MESSAGE("Unknown collision type.");
        }

    }

}

void CollisionHandler::initCollisions ()
{

    if (m_contains_coulomb_type) { // Define lambda_debye MultiFab
        const WarpX& warpx = WarpX::GetInstance();
        const int num_levels = warpx.maxLevel() + 1;
        const int ncomp = 1;
        m_lambda_debye.resize(num_levels); // size is number of levels
        for (int lev = 0; lev < num_levels; ++lev) {
            const amrex::BoxArray ba = warpx.boxArray(lev);
            const amrex::DistributionMapping dm = warpx.DistributionMap(lev);
            m_lambda_debye[lev] = std::make_unique<amrex::MultiFab>( ba, dm, ncomp,
                                                                    amrex::IntVect::TheZeroVector() );
            m_lambda_debye[lev]->setVal(0.0);
        }

    }

}

/** Perform all collisions
 *
 * @param cur_time Current time
 * @param dt time step size
 * @param mypc MultiParticleContainer calling this method
 *
 */
void CollisionHandler::doCollisions ( amrex::Real cur_time, amrex::Real dt, MultiParticleContainer* mypc)
{
    // Compute the plasma Debye length for Coulomb collisions
    // LDe^{-2} = sum_{species} LDe_s^{-2}, LDe_s^{-2} = ns*qs^2/(kTs*ep0)
    if (m_contains_coulomb_type) {
        const auto nSpecies = mypc->nSpecies();
        for (int lev = 0; lev < m_lambda_debye.size(); ++lev) {
            m_lambda_debye[lev]->setVal(0.0);
        }
        for (int sp = 0; sp<nSpecies; sp++) {
            auto& species = mypc->GetParticleContainer(sp);
            const amrex::ParticleReal charge = species.getCharge();
            if (charge==0.0) { continue; }
            // Here is where I need to compute species Debye length
        }
    }

    for (auto& collision : allcollisions) {
        int const ndt = collision->get_ndt();
        if ( int(std::floor(cur_time/dt)) % ndt == 0 ) {
            collision->doCollisions(cur_time, dt*ndt, mypc);
        }
    }

}
