"""Illustrative single-scattering slab response, not calibrated Cassini photometry."""
import numpy as np


def slab(tau, sun_elevation_deg=25.0, observer_elevation_deg=30.0):
    tau = np.asarray(tau, dtype=float)
    if not np.isfinite(tau).all() or np.any(tau < 0):
        raise ValueError("Optical depth must be finite and nonnegative")
    if not 0 < sun_elevation_deg <= 90 or not 0 < observer_elevation_deg <= 90:
        raise ValueError("Elevations must be in (0, 90] degrees")
    mu0, mu = np.sin(np.deg2rad([sun_elevation_deg, observer_elevation_deg]))
    return mu0 / (mu + mu0) * -np.expm1(-tau * (1 / mu + 1 / mu0))


def shadow_length(height_m, sun_elevation_deg):
    if not np.isfinite(height_m) or height_m < 0 or not 0 < sun_elevation_deg < 90:
        raise ValueError("Need a nonnegative height and a solar elevation in (0,90)")
    return height_m / np.tan(np.deg2rad(sun_elevation_deg))
