"""Restricted Hill encounters; input lengths use SI, integration uses Hill units."""

from dataclasses import asdict, dataclass
import numpy as np
from scipy.integrate import solve_ivp

G = 6.67430e-11
SATURN_GM = 3.7931208e16
A_MOON = 136_505_000.0
MOON_MASS = 7.7e13
OMEGA = np.sqrt(SATURN_GM / A_MOON**3)
PERIOD = 2 * np.pi / OMEGA
HILL = (G * MOON_MASS / (3 * OMEGA**2)) ** (1 / 3)


@dataclass(frozen=True)
class Encounter:
    mass_factor: float = 1.0
    radial_epicycle_m: float = 0.0
    vertical_amplitude_m: float = 0.0
    phase: float = 0.0
    start: float = -120.0
    end: float = 12 * np.pi
    samples: int = 1601
    rtol: float = 2e-10

    def __post_init__(self):
        values = np.array(list(asdict(self).values()), dtype=float)
        if (
            not np.isfinite(values).all()
            or self.mass_factor < 0
            or self.end <= self.start
        ):
            raise ValueError("Invalid encounter parameters")
        if self.samples < 2 or self.rtol <= 0:
            raise ValueError("Need at least two samples and a positive tolerance")


def moon(t, config):
    p = t + config.phase
    e = config.radial_epicycle_m / HILL
    i = config.vertical_amplitude_m / HILL
    return np.array([-e * np.cos(p), 2 * e * np.sin(p), i * np.sin(p)])


class SurfaceImpact(ValueError):
    pass


def rhs(t, flat, config):
    q = flat.reshape(6, -1)
    d = q[:3] - moon(t, config)[:, None]
    r = np.sqrt(np.sum(d * d, axis=0))
    f = -3 * config.mass_factor * d / r**3
    return np.vstack(
        (q[3:], 2 * q[4] + 3 * q[0] + f[0], -2 * q[3] + f[1], -q[2] + f[2])
    ).ravel()


def integrate(impact_m, config=Encounter()):
    b = np.asarray(impact_m, dtype=float) / HILL
    if b.ndim != 1 or len(b) == 0 or not np.isfinite(b).all() or np.any(abs(b) < 2):
        raise ValueError(
            "Impact parameters must be a finite vector outside two Hill radii"
        )
    q = np.zeros((6, len(b)))
    q[0] = b
    q[1] = -1.5 * b * config.start
    q[4] = -1.5 * b
    t = np.linspace(config.start, config.end, config.samples)

    def surface(t, flat, config):
        delta = flat.reshape(6, -1)[:3] - moon(t, config)[:, None]
        return np.min(np.linalg.norm(delta, axis=0)) * HILL - 3900

    surface.terminal = True
    surface.direction = -1
    sol = solve_ivp(
        rhs,
        (t[0], t[-1]),
        q.ravel(),
        args=(config,),
        t_eval=t,
        events=surface,
        method="DOP853",
        rtol=config.rtol,
        atol=config.rtol * 0.01,
    )
    if sol.t_events[0].size:
        raise SurfaceImpact("Trajectory intersects the adopted Daphnis surface")
    if not sol.success:
        raise RuntimeError(sol.message)
    return t, sol.y.reshape(6, len(b), -1).transpose(2, 1, 0)


def elements(q):
    x, y, z, vx, vy, vz = np.moveaxis(q, -1, 0)
    guiding = 2 * (vy + 2 * x)
    radial = np.hypot(x - guiding, vx)
    vertical = np.hypot(z, vz)
    return guiding, radial, vertical


def jacobi(q, mass_factor=1):
    x, y, z, vx, vy, vz = np.moveaxis(q, -1, 0)
    return (
        3 * x * x
        - z * z
        + 6 * mass_factor / np.sqrt(x * x + y * y + z * z)
        - vx * vx
        - vy * vy
        - vz * vz
    )


def crowding(t, impact_m, q):
    # Flux conservation applies before folds; it does not resolve collisional caustics.
    b = np.asarray(impact_m) / HILL
    db = np.gradient(q[..., :2], b, axis=1, edge_order=2)
    dt = np.gradient(q[..., :2], t, axis=0, edge_order=2)
    det = db[..., 0] * dt[..., 1] - db[..., 1] * dt[..., 0]
    return (1.5 * np.abs(b))[None, :] / np.abs(det), det
