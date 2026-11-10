# Model and validation

## Restricted encounter

The rotating frame has radial coordinate `x`, azimuthal coordinate `y`, and vertical coordinate `z`. The reference orbit is circular, at 136,505 km from Saturn. Positive `x` is outward. Outer-ring material moves in negative `y`. Inner-ring material has the opposite symmetry for a circular, coplanar moon.

The adopted nominal constants are `G = 6.67430e-11 m^3 kg^-1 s^-2`, Saturn `GM = 3.7931208e16 m^3 s^-2`, and Daphnis mass `7.7e13 kg`. They are model inputs, not a new mass determination. The mean motion is `n = sqrt(GM/a^3)`. The Hill radius is `rH = (Gm/(3n^2))^(1/3)`.

The equations in physical units are:

```text
xddot - 2 n ydot - 3 n^2 x = -Gm (x-xm) / distance^3
yddot + 2 n xdot           = -Gm (y-ym) / distance^3
zddot + n^2 z              = -Gm (z-zm) / distance^3
```

The moon follows a prescribed free epicycle:

```text
xm = -ae cos(nt + phase)
ym = 2 ae sin(nt + phase)
zm = ai sin(nt + phase)
```

The eccentric and inclined control uses `ae = 4.5 km` and `ai = 8.6 km`. The phase sweep ties the radial and vertical phases together. This samples one family of geometries; it does not sample all possible arguments of periapse and node. There is no Saturn J2, orbital precession, moon back-reaction, grain self-gravity, or collision in this calculation. It is a local Hill model, not an ephemeris of either image epoch. Six local orbital periods represent travel downstream, not six encounters with the moon. The next synodic encounter takes much longer.

At dimensionless time `s = -120`, each particle has `x = b`, `y = -1.5 b s`, `vx = vz = z = 0`, and `vy = -1.5 b`. Each `b` has a different upstream distance. DOP853 integrates in Hill units with relative tolerance `2e-10`. A terminal surface event rejects a trajectory that reaches the adopted 3.9 km spherical moon surface. Such phases have `NaN` orbital results and an explicit impact flag in the phase CSV. They are not scattering predictions.

The osculating Hill guiding center is `X = 2(vy/n + 2x)`. This approximates semimajor-axis offset from the reference orbit to first order in local displacement. The radial epicycle amplitude is `A = sqrt((x-X)^2 + (vx/n)^2)`. The vertical oscillator amplitude is `Z = sqrt(z^2 + (vz/n)^2)`. These are not the instantaneous radius, visible edge, or RMS height of an ensemble.

For a stationary moon the Jacobi invariant in Hill units is:

```text
C = 3 x^2 - z^2 + 6 mass_factor / distance - (vx^2 + vy^2 + vz^2)
```

Far from a circular moon, `C = 3 X^2/4 - A^2 - Z^2`. The tests verify this identity, its finite-distance potential correction, mirror symmetry, free shear, and planar symmetry. A tighter solver tolerance and a larger initial separation test different errors. The report retains both.

For a circular moon, the streamline map `(b, s) -> (x, y)` gives a surface-density ratio before crossing:

```text
Sigma / Sigma0 = 1.5 |b| / |det(d(x,y)/d(b,s))|
```

A vanishing or reversed Jacobian identifies a fold where the collisionless solution fails as a dense-ring description. No smoothing, pressure, or damping is added to hide that failure. This density diagnostic is not applied to a moving eccentric/inclined moon.

## Grain box

REBOUND 4.4.10 supplies its symplectic epicycle integrator (SEI), shear-periodic horizontal boundaries, a gravity tree, and hard-sphere collisions. Each particle is one grain. Mass is `4 pi rho R^3/3`. The radius sweep holds surface density and material density fixed, so `tau = 3 Sigma/(4 rho R)` apart from integer particle-count rounding. No particle mass, collision radius, or gravity multiplier is enlarged to save compute.

The base box is 128 by 128 m. Surface density is 300 kg/m² and grain material density is 450 kg/m³. Radius is 2 m. There are 326 particles. The initial horizontal dispersions are 0.02 m/s radially and 0.01 m/s azimuthally after subtraction of Keplerian shear. Vertical position has a 5 m Gaussian scale; vertical speed has dispersion `n * 5 m`. Initial placements reject sphere overlap, including horizontal periodic images.

The normal restitution coefficient is a constant 0.5. Smooth spheres have no spin, tangential friction, fragmentation, adhesion, size distribution, or persistent-contact model. Elastic and collisionless cases are controls. The constant is an experimental assumption, not a calibrated ice law. A separate control uses the Bridges form from REBOUND’s ring example: `e = min(1, 0.32 (100 |v_normal|)^(-0.234))`, with speed in m/s and `e = 1` at zero speed. This laboratory-based law is also not a calibration of the observed ring edge. The optional self-gravity calculation uses no force softening. One layer of horizontal ghost boxes approximates neighboring mass; the extended suite tests two layers and a tighter opening angle. It does not calculate an exact infinite periodic gravity sum. The gravity-accuracy controls cover the first half orbit, through the initial heating peak; the baseline gravity run covers three orbits. Late-time gravity convergence is not established.

The vertical extent is 64 root cells. REBOUND's shear boundary also wraps vertically. We therefore provide a much taller box and reject outputs whose vertical oscillator amplitude approaches one quarter of that extent. Zero vertical ghost boxes alone do not disable vertical wrapping. This is covered by a regression test.

An additional high-energy control uses radial dispersion 0.09 m/s, the Bridges law, and 1 s and 0.5 s steps. This input speed is about one third of the 20 km encounter’s radial epicycle speed. The input remains random and homogeneous; it is not a coupled encounter simulation.

Time steps are 2 s in the base case, 1 s for the 1 m grains and the time-step control, and 0.5 s for the elastic control. At each output, the code rejects `maximum peculiar speed * dt / radius > 0.3`. This is a sampled diagnostic, not a proof that every impact is resolved. A half-step run tests the resulting statistics. Close packing or sustained contacts would require a soft-sphere solver. The current results cannot validate that regime.

RMS height subtracts the mean height. The 99th percentile uses absolute height above the reference plane. Vertical specific oscillator energy is `mean(vz^2 + n^2 z^2)/2`. In a collisionless control it stays constant even when phases mix. This distinguishes energy loss from a change in the appearance of an ensemble. With self-gravity, this quantity excludes mutual gravitational potential energy and is not a closed energy budget.

The imposed horizontal motion is random and homogeneous. It is not the coherent, kilometer-scale streamline crossing produced by Daphnis. The box can show that collisions convert horizontal energy into vertical motion and later remove energy. It cannot measure the downstream damping length of the gap edge.

## Optical illustration

For a uniform slab viewed and illuminated from the same side, the single-scattering response is proportional to:

```text
mu0 / (mu + mu0) * [1 - exp(-tau (1/mu + 1/mu0))]
```

`mu0` and `mu` are sines of Sun and observer elevations. A shared albedo and phase-function factor is omitted. The illustration uses 25° solar and 30° observer elevations. These are chosen values, not recovered image geometry. Opaque clumps and empty space have a different mean response than a uniform slab of the same mass. No ray tracing, mutual shadowing, multiple scattering, or Cassini camera response is calculated.
