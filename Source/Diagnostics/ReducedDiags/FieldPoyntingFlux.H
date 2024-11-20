/* Copyright 2019-2020
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPOYTINGFLUX_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPOYTINGFLUX_H_

#include "ReducedDiags.H"

#include <string>

/**
 * \brief This class mainly contains a function that computes the field Poynting flux.
 */
class FieldPoyntingFlux : public ReducedDiags
{
public:

    /**
     * \brief Constructor
     *
     * \param[in] rd_name reduced diags names
     */
    FieldPoyntingFlux(const std::string& rd_name);

    /**
     * \brief This function computes the electromagnetic Poynting flux,
     * obtained by integrating the electromagnetic Poynting flux density g = eps0 * (E x B)
     * on the surface of the domain.
     *
     * \param[in] step current time step
     */
    void ComputeDiags(int step) final;
};

#endif
