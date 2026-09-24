#!/usr/bin/env python3

import sys

import numpy as np
import yt
from scipy.constants import e, epsilon_0

# Compare charge density with Gauss's law on the final diagnostic step.
# Exclude the cells adjacent to the radial boundaries.
pltdir = sys.argv[1]
ds = yt.load(pltdir)
data = ds.covering_grid(
    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
)

rho = (
    (data[("boxlib", "rho_electrons")] + data[("boxlib", "rho_ions")])
    .to_ndarray()
    .squeeze()
)
divE = data[("boxlib", "divE")].to_ndarray().squeeze()

assert rho.ndim == 1 and divE.ndim == 1, "Expected 1D RSPHERE diagnostic fields"
assert rho.size > 2, "Need interior cells after excluding boundary-adjacent cells"

charge_error = rho[1:-1] - epsilon_0 * divE[1:-1]  # ignore values next to boundary.
# Normalize to the characteristic charge density e*n0 from the input deck.
n0 = 1.0e22
error_rms = np.sqrt(np.mean((charge_error / (e * n0)) ** 2))
tolerance = 1.0e-14

print(f"final-step interior RMS charge-conservation error: {error_rms:.6e}")
print(f"tolerance: {tolerance:.6e}")
assert error_rms < tolerance
