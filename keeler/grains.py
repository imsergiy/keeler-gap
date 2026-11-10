"""Physical-grain, homogeneous shearing box; no gap edge is represented."""
from dataclasses import dataclass, asdict
import time
import numpy as np
import rebound
from .encounter import G, OMEGA, PERIOD


@dataclass(frozen=True)
class Box:
    radius_m: float = 2.0
    size_m: float = 128.0
    surface_density: float = 300.0
    material_density: float = 450.0
    horizontal_sigma: float = .02
    height_m: float = 5.0
    orbits: float = 3.0
    dt_s: float = 2.0
    restitution: float | str = .5
    collisions: bool = True
    self_gravity: bool = False
    seed: int = 42
    opening_angle: float = .3
    ghosts: int = 1

    def __post_init__(self):
        positive = [self.radius_m, self.size_m, self.surface_density, self.material_density,
                    self.height_m, self.orbits, self.dt_s, self.opening_angle]
        if not np.isfinite(positive).all() or min(positive) <= 0:
            raise ValueError("Box scales must be finite and positive")
        if not np.isfinite(self.horizontal_sigma) or self.horizontal_sigma < 0:
            raise ValueError("Velocity dispersion must be finite and nonnegative")
        valid_cor = self.restitution == "bridges" if isinstance(self.restitution,str) else np.isfinite(self.restitution) and 0 <= self.restitution <= 1
        if not valid_cor or self.ghosts < 1:
            raise ValueError("Invalid collision or ghost parameters")
        if self.size_m < 8 * self.radius_m:
            raise ValueError("Box must span at least four grain diameters")

    @property
    def grain_mass(self):
        return 4 * np.pi / 3 * self.material_density * self.radius_m**3

    @property
    def count(self):
        return int(round(self.surface_density * self.size_m**2 / self.grain_mass))

    @property
    def optical_depth(self):
        return self.count * np.pi * self.radius_m**2 / self.size_m**2


def initial_particles(config):
    rng = np.random.default_rng(config.seed)
    n = config.count
    if n < 2 or n > 2_000_000:
        raise ValueError(f"Requested {n} grains; use 2 to 2,000,000 and budget before running")
    q = np.zeros((n, 6))
    cells = {}
    L, r = config.size_m, config.radius_m
    nc = int(L / (2 * r))
    cell = L / nc
    accepted = 0
    for _ in range(1000 * n):
        p = np.array([rng.uniform(-L/2, L/2), rng.uniform(-L/2, L/2),
                      rng.normal(0, config.height_m)])
        key = np.floor((p + np.array([L/2,L/2,0])) / cell).astype(int)
        near = []
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for dz in (-1,0,1):
                    near.extend(cells.get(((key[0]+dx)%nc, (key[1]+dy)%nc, key[2]+dz), []))
        if near:
            d = q[near, :3] - p
            d[:, :2] -= L * np.round(d[:, :2] / L)
            if np.any(np.sum(d*d, axis=1) < (2*r)**2):
                continue
        q[accepted, :3] = p
        cells.setdefault(tuple(key), []).append(accepted)
        accepted += 1
        if accepted == n:
            break
    if accepted != n:
        raise RuntimeError("Cannot place nonoverlapping grains at this density and height")
    q[:, 3] = rng.normal(0, config.horizontal_sigma, n)
    q[:, 4] = -1.5 * OMEGA * q[:, 0] + rng.normal(0, config.horizontal_sigma/2, n)
    q[:, 5] = rng.normal(0, OMEGA * config.height_m, n)
    q[:, 3] -= np.mean(q[:, 3])
    q[:, 4] -= np.mean(q[:, 4] + 1.5 * OMEGA * q[:, 0])
    q[:, 5] -= np.mean(q[:, 5])
    return q


def coefficient(law, normal_speed):
    if law == "bridges":
        speed=abs(normal_speed)
        return min(1., .32*(100*speed)**(-.234)) if speed else 1.
    return law


def setup(config):
    sim = rebound.Simulation()
    sim.G = G
    sim.integrator = "sei"
    sim.ri_sei.OMEGA = OMEGA
    sim.ri_sei.OMEGAZ = OMEGA
    sim.boundary = "shear"
    sim.configure_box(config.size_m, 1, 1, 64)
    sim.N_ghost_x = sim.N_ghost_y = config.ghosts
    sim.N_ghost_z = 0
    sim.gravity = "tree" if config.self_gravity else "none"
    sim.collision = "tree" if config.collisions else "none"
    sim.collision_resolve = "hardsphere"
    sim.opening_angle2 = config.opening_angle**2
    sim.softening = 0.0
    sim.rand_seed = config.seed
    sim.dt = config.dt_s
    def restitution(sim_pointer, speed):
        return coefficient(config.restitution,speed)
    sim.coefficient_of_restitution = restitution
    q = initial_particles(config)
    for x,y,z,vx,vy,vz in q:
        sim.add(m=config.grain_mass, r=config.radius_m,
                x=x,y=y,z=z,vx=vx,vy=vy,vz=vz)
    return sim, restitution


def state(sim):
    q = np.empty((sim.N, 6))
    for i,p in enumerate(sim.particles):
        q[i] = p.x,p.y,p.z,p.vx,p.vy,p.vz
    return q


def diagnostics(q):
    vx, vy, vz = q[:,3], q[:,4] + 1.5*OMEGA*q[:,0], q[:,5]
    return [np.std(q[:,2]), np.quantile(np.abs(q[:,2]), .99),
            np.std(vx), np.std(vy), np.std(vz),
            .5*np.mean(vz*vz + OMEGA**2*q[:,2]**2),
            np.max(np.linalg.norm(np.column_stack((vx,vy,vz)), axis=1))]


def run(config, samples=121, progress=False):
    sim, callback = setup(config)
    rows, snapshots = [], []
    started = time.monotonic()
    nsteps = int(np.ceil(config.orbits * PERIOD / config.dt_s))
    output_steps = np.unique(np.linspace(0, nsteps, samples).astype(int))
    previous = 0
    max_crossing = 0.0
    for index,step in enumerate(output_steps):
        sim.steps(int(step-previous))
        previous = step
        q = state(sim)
        if sim.N != config.count or not np.isfinite(q).all():
            raise RuntimeError("Particle loss or nonfinite state")
        if np.max(np.hypot(q[:,2], q[:,5]/OMEGA)) > 16*config.size_m:
            raise RuntimeError("Vertical excursion approaches the periodic boundary; increase vertical extent")
        d = diagnostics(q)
        crossing = d[-1] * config.dt_s / config.radius_m
        max_crossing = max(max_crossing, crossing)
        if crossing > .3:
            raise RuntimeError(f"Step too large: speed*dt/r={crossing:.3g}; reduce dt_s")
        rows.append([sim.t/PERIOD, *d])
        if progress and index % 20 == 0:
            print(f"  orbit={sim.t/PERIOD:.3f} height_rms_m={d[0]:.3f}",flush=True)
        if step in (output_steps[0], output_steps[len(output_steps)//2], output_steps[-1]):
            snapshots.append(q)
    return np.asarray(rows), np.asarray(snapshots), {
        "config": asdict(config), "particles": sim.N, "tau": config.optical_depth,
        "sigma_actual": sim.N*config.grain_mass/config.size_m**2,
        "seconds": time.monotonic()-started, "max_speed_dt_over_radius": max_crossing,
        "rebound_version": rebound.__version__,
        "columns": ["orbits", "height_rms_m", "height_abs_p99_m", "vx_std_m_s",
                    "vy_peculiar_std_m_s", "vz_std_m_s", "vertical_energy_J_kg", "max_peculiar_speed_m_s"],
        "meaning": "Homogeneous excited patch; no Daphnis, radial edge, or photometric fit",
    }
