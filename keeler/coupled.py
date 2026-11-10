"""Direct 3D moon encounter of finite, physical-grain packets."""

import argparse
from dataclasses import asdict, dataclass, replace
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import rebound
from scipy.spatial import cKDTree

from .encounter import G, HILL, MOON_MASS, OMEGA, PERIOD, elements
from .grains import coefficient, state

MAX_SIDE = 64
CORE_SIDE = 8
ROOT_M = 4_000_000.0
MOON_RADIUS_M = 3900.0


@dataclass(frozen=True)
class Packet:
    count: int = 128
    radius_m: float = 2.0
    material_density: float = 450.0
    surface_density: float = 300.0
    center_m: float = 21_250.0
    packet_centres_xy_m: tuple | None = None
    height_m: float = 5.0
    upstream_orbits: float = 2.0
    downstream_orbits: float = 3.0
    dt_s: float = 2.0
    seed: int = 42
    moon_vertical_amplitude_m: float = 0.0
    moon_vertical_phase_rad: float = 0.0
    moon_gravity: bool = True
    self_gravity: bool = True
    collisions: bool = True
    restitution: float | str = "bridges"

    def __post_init__(self):
        if not isinstance(self.count, int) or self.count < 8:
            raise ValueError("Count must be twice an even square, at least eight")
        side = int(np.sqrt(self.count // 2))
        if self.count != 2 * side * side or side % 2 or side > MAX_SIDE:
            raise ValueError(
                "Count must be twice an even square, with at most 8192 grains"
            )
        values = [
            self.radius_m,
            self.material_density,
            self.surface_density,
            self.center_m,
            self.height_m,
            self.upstream_orbits,
            self.downstream_orbits,
            self.dt_s,
        ]
        if not np.isfinite(values).all() or min(values) <= 0:
            raise ValueError("Physical scales must be finite and positive")
        if self.cell_m <= 2 * self.radius_m:
            raise ValueError(
                "These initial cells require center spacing greater than a grain diameter"
            )
        centers = (
            np.array([[self.center_m, 0], [-self.center_m, 0]])
            if self.packet_centres_xy_m is None
            else np.asarray(self.packet_centres_xy_m, dtype=float)
        )
        if centers.shape != (2, 2) or not np.isfinite(centers).all():
            raise ValueError("Exactly two finite packet centres are required")
        if np.min(abs(centers[:, 0])) - self.side_length_m / 2 <= 2 * np.sqrt(3) * HILL:
            raise ValueError(
                "This packet experiment starts outside the close Hill-encounter region"
            )
        if self.restitution != "bridges" and not (
            isinstance(self.restitution, (float, int)) and 0 <= self.restitution <= 1
        ):
            raise ValueError("Restitution must be bridges or a number in [0,1]")
        if (
            not np.isfinite(
                [self.moon_vertical_amplitude_m, self.moon_vertical_phase_rad]
            ).all()
            or self.moon_vertical_amplitude_m < 0
        ):
            raise ValueError("Moon vertical amplitude must be finite and nonnegative")
        if self.moon_vertical_amplitude_m and not self.self_gravity:
            raise ValueError(
                "This inclined-moon experiment requires all grain gravity active"
            )
        if not isinstance(self.seed, int) or self.seed < 0:
            raise ValueError("Seed must be a nonnegative integer")

    @property
    def grain_mass(self):
        return 4 * np.pi / 3 * self.material_density * self.radius_m**3

    @property
    def cell_m(self):
        return np.sqrt(self.grain_mass / self.surface_density)

    @property
    def side_length_m(self):
        return np.sqrt(self.count / 2) * self.cell_m


def initial_packet(config):
    side = int(np.sqrt(config.count // 2))
    rows = []
    ids = []
    sides = []
    core = []
    jitter = 0.4 * (config.cell_m - 2 * config.radius_m)
    cold_speed = OMEGA * config.height_m
    centers = (
        [
            (
                sign * config.center_m,
                sign * 1.5 * OMEGA * config.center_m * config.upstream_orbits * PERIOD,
            )
            for sign in [1, -1]
        ]
        if config.packet_centres_xy_m is None
        else config.packet_centres_xy_m
    )
    for side_index, sign in enumerate([1, -1]):
        for i in range(-side // 2, side // 2):
            for j in range(-side // 2, side // 2):
                rng = np.random.default_rng(
                    np.random.SeedSequence(
                        [config.seed, side_index, i % 2**32, j % 2**32]
                    )
                )
                x = (
                    centers[side_index][0]
                    + (i + 0.5) * config.cell_m
                    + rng.uniform(-jitter, jitter)
                )
                y = (
                    centers[side_index][1]
                    + (j + 0.5) * config.cell_m
                    + rng.uniform(-jitter, jitter)
                )
                z = rng.normal(0, config.height_m)
                rows.append(
                    [
                        x,
                        y,
                        z,
                        rng.normal(0, cold_speed),
                        -1.5 * OMEGA * x + rng.normal(0, cold_speed / 2),
                        rng.normal(0, cold_speed),
                    ]
                )
                ids.append(
                    1
                    + side_index * MAX_SIDE**2
                    + (i + MAX_SIDE // 2) * MAX_SIDE
                    + j
                    + MAX_SIDE // 2
                )
                sides.append(sign)
                core.append(
                    -CORE_SIDE // 2 <= i < CORE_SIDE // 2
                    and -CORE_SIDE // 2 <= j < CORE_SIDE // 2
                )
    return {
        "state": np.asarray(rows),
        "ids": np.asarray(ids),
        "side": np.asarray(sides),
        "core": np.asarray(core),
    }


def install_collision_law(sim, config):
    contacts = {
        "approaching_calls": 0,
        "dissipated_J": 0.0,
        "max_normal_speed_m_s": 0.0,
    }

    def restitution(pointer, speed):
        e = coefficient(config.restitution, speed)
        if speed < -1e-14:
            contacts["approaching_calls"] += 1
            contacts["max_normal_speed_m_s"] = max(
                contacts["max_normal_speed_m_s"], abs(speed)
            )
            contacts["dissipated_J"] += (
                config.grain_mass * 0.25 * (1 - e * e) * speed * speed
            )
        return e

    sim.coefficient_of_restitution = restitution
    return contacts


def setup(config):
    layout = initial_packet(config)
    sim = rebound.Simulation()
    sim.G = G
    sim.integrator = "sei"
    sim.ri_sei.OMEGA = sim.ri_sei.OMEGAZ = OMEGA
    sim.dt = config.dt_s
    sim.boundary = "none"
    sim.configure_box(ROOT_M)
    sim.N_ghost_x = sim.N_ghost_y = sim.N_ghost_z = 0
    sim.softening = 0
    sim.rand_seed = config.seed
    if config.self_gravity:
        sim.gravity = "basic"
        sim.N_active = -1
    elif config.moon_gravity:
        sim.gravity = "basic"
        sim.N_active = 1
        sim.testparticle_type = 0
    else:
        sim.gravity = "none"
    # Contact-tree modes must be enabled before adding bodies to populate the index.
    sim.collision = "tree" if config.collisions else "none"
    if config.moon_gravity:
        sim.add(
            m=MOON_MASS,
            r=0,
            hash=0,
            z=config.moon_vertical_amplitude_m * np.sin(config.moon_vertical_phase_rad),
            vz=OMEGA
            * config.moon_vertical_amplitude_m
            * np.cos(config.moon_vertical_phase_rad),
        )
    for identity, row in zip(layout["ids"], layout["state"]):
        x, y, z, vx, vy, vz = row
        sim.add(
            m=config.grain_mass,
            r=config.radius_m,
            hash=int(identity),
            x=x,
            y=y,
            z=z,
            vx=vx,
            vy=vy,
            vz=vz,
        )
    contacts = install_collision_law(sim, config)
    sim.collision_resolve = "hardsphere"
    sim.minimum_collision_velocity = 0
    return sim, layout, contacts


def read_frame(sim, config, layout):
    identities = [p.hash.value for p in sim.particles]
    expected = np.r_[0, layout["ids"]] if config.moon_gravity else layout["ids"]
    if len(set(identities)) != len(expected) or set(identities) != set(expected):
        raise RuntimeError(
            "Grain identities changed; no removal or merger is permitted"
        )
    indices = {identity: i for i, identity in enumerate(identities)}
    ordered = state(sim)[[indices[int(identity)] for identity in expected]]
    if config.moon_gravity:
        if not config.self_gravity and identities[0] != 0:
            raise RuntimeError("The active moon moved out of the library active slot")
        return ordered[1:], ordered[0].copy()
    return ordered, np.zeros(6)


def energy(sim, config, q, moon):
    tidal = config.grain_mass * np.sum(
        -1.5 * OMEGA**2 * q[:, 0] ** 2 + 0.5 * OMEGA**2 * q[:, 2] ** 2
    )
    if config.self_gravity:
        if config.moon_gravity:
            tidal += MOON_MASS * (
                -1.5 * OMEGA**2 * moon[0] ** 2 + 0.5 * OMEGA**2 * moon[2] ** 2
            )
        return sim.energy() + tidal
    result = 0.5 * config.grain_mass * np.sum(q[:, 3:] ** 2) + tidal
    if config.moon_gravity:
        result -= (
            G
            * MOON_MASS
            * config.grain_mass
            * np.sum(1 / np.linalg.norm(q[:, :3] - moon[:3], axis=1))
        )
    return result


def deep_overlap_fraction(q, radius):
    pairs = cKDTree(q[:, :3]).query_pairs(1.8 * radius, output_type="ndarray")
    return float(len(np.unique(pairs)) / len(q)) if len(pairs) else 0.0


def packet_step_ratio(q, side, dt, radius):
    maximum = 0.0
    for sign in [1, -1]:
        velocities = q[side == sign, 3:]
        residual = velocities - np.mean(velocities, axis=0)
        maximum = max(maximum, float(np.max(np.linalg.norm(residual, axis=1))))
    return maximum * dt / radius


COLUMNS = [
    "orbits",
    "core_packet_a_guiding_shift_m",
    "core_packet_b_guiding_shift_m",
    "core_height_rms_m",
    "core_abs_height_p99_m",
    "core_mean_z_m",
    "core_vertical_oscillator_energy_J_kg",
    "contacts",
    "dissipated_J",
    "rotating_energy_J",
    "corrected_energy_error_J",
    "deep_overlap_fraction",
    "packet_speed_dt_over_radius",
]


def run(config, samples=201):
    sim, layout, contacts = setup(config)
    initial_guiding = (
        elements(
            layout["state"]
            / np.array([HILL, HILL, HILL, OMEGA * HILL, OMEGA * HILL, OMEGA * HILL])
        )[0]
        * HILL
    )
    steps = np.unique(
        np.linspace(
            0,
            int(
                np.ceil(
                    (config.upstream_orbits + config.downstream_orbits)
                    * PERIOD
                    / config.dt_s
                )
            ),
            samples,
        ).astype(int)
    )
    rows = []
    snapshots = []
    moons = []
    before = 0
    initial_energy = None
    started = time.monotonic()
    for index, step in enumerate(steps):
        sim.steps(int(step - before))
        before = step
        q, moon = read_frame(sim, config, layout)
        if (
            not np.isfinite(q).all()
            or not np.isfinite(moon).all()
            or max(np.max(abs(q[:, :3])), np.max(abs(moon[:3]))) > 0.45 * ROOT_M
        ):
            raise RuntimeError(
                "Nonfinite state or approach to the gravity-tree boundary"
            )
        if (
            config.moon_gravity
            and np.min(np.linalg.norm(q[:, :3] - moon[:3], axis=1))
            <= MOON_RADIUS_M + config.radius_m
        ):
            raise RuntimeError("Moon impact leaves this experiment domain")
        overlap = deep_overlap_fraction(q, config.radius_m)
        if config.collisions and overlap > 0.02:
            raise RuntimeError(
                f"Deep overlap in {overlap:.1%} of grains; hard-sphere step/contact model is unresolved"
            )
        crossing = packet_step_ratio(q, layout["side"], config.dt_s, config.radius_m)
        if crossing > 0.3:
            raise RuntimeError(f"Contact step bound exceeded: {crossing:g}; reduce dt")
        core = q[layout["core"]]
        guiding = 2 * (q[:, 4] / OMEGA + 2 * q[:, 0])
        shifts = guiding - initial_guiding
        outer = layout["core"] & (layout["side"] == 1)
        inner = layout["core"] & (layout["side"] == -1)
        hamiltonian = energy(sim, config, q, moon)
        if initial_energy is None:
            initial_energy = hamiltonian
        rows.append(
            [
                sim.t / PERIOD,
                np.mean(shifts[outer]),
                np.mean(shifts[inner]),
                np.std(core[:, 2]),
                np.quantile(abs(core[:, 2]), 0.99),
                np.mean(core[:, 2]),
                0.5 * np.mean(core[:, 5] ** 2 + OMEGA**2 * core[:, 2] ** 2),
                contacts["approaching_calls"],
                contacts["dissipated_J"],
                hamiltonian,
                hamiltonian + contacts["dissipated_J"] - initial_energy,
                overlap,
                crossing,
            ]
        )
        snapshots.append(q.copy())
        moons.append(moon)
        if index % 20 == 0:
            print(
                f"  t={sim.t / PERIOD:.2f} orbits H={rows[-1][3]:.3f}m contacts={contacts['approaching_calls']} balance={rows[-1][10]:.4g}J",
                flush=True,
            )
    rows = np.asarray(rows)
    meta = {
        "config": asdict(config),
        "grain_mass_kg": config.grain_mass,
        "side_length_m": config.side_length_m,
        "actual_surface_density_kg_m2": config.count
        * config.grain_mass
        / (2 * config.side_length_m**2),
        "normal_optical_depth": np.pi * config.radius_m**2 / config.cell_m**2,
        "initial_thermal_speed_m_s": OMEGA * config.height_m,
        "core_grains": int(layout["core"].sum()),
        "columns": COLUMNS,
        "elapsed_seconds": time.monotonic() - started,
        "rebound_version": rebound.__version__,
        "gravity_solver": sim.gravity,
        "contact_max_normal_speed_m_s": contacts["max_normal_speed_m_s"],
        "contact_tree_extent_m": ROOT_M,
        "model": "Finite open packets; Hill planetary tide; massive moon and grains; physical hard spheres; no target waves or heating kicks",
        "max_sampled_step_ratio": float(rows[:, -1].max()),
        "max_sampled_deep_overlap_fraction": float(rows[:, -2].max()),
        "energy_diagnostic_note": "For a vertically moving moon, total energy is dominated by the moon oscillator; its roundoff limits grain-scale budget interpretation. Use trajectory and step checks as well."
        if config.moon_vertical_amplitude_m
        else "Rotating energy includes moon and grain potentials; collision loss is added back.",
        "max_corrected_energy_error_J": float(abs(rows[:, 10]).max()),
    }
    return rows, np.asarray(snapshots), np.asarray(moons), layout, meta


def save_case(out, name, config):
    if (out / f"{name}.json").exists():
        raise FileExistsError("Use a new case/output directory")
    print("RUN", name, config.count, flush=True)
    rows, q, moon, layout, meta = run(config)
    np.savetxt(
        out / f"{name}.csv", rows, delimiter=",", header=",".join(COLUMNS), comments=""
    )
    np.savez_compressed(
        out / f"{name}.npz",
        orbits=rows[:, 0],
        state=q,
        moon=moon,
        ids=layout["ids"],
        side=layout["side"],
        core=layout["core"],
        initial_state=layout["state"],
    )
    (out / f"{name}.json").write_text(
        json.dumps(meta, indent=2, allow_nan=False) + "\n"
    )
    print("DONE", name, json.dumps(meta), flush=True)
    return rows, q, moon, layout, meta


def plan_stream_meeting():
    from scipy.integrate import solve_ivp
    from scipy.optimize import root
    from .encounter import Encounter, rhs

    impact = np.array([21_000.0, 23_000.0])
    event_orbits = 4.0
    end = event_orbits * 2 * np.pi

    def final(y):
        q = np.zeros((6, 2))
        q[0] = impact / HILL
        q[1] = y / HILL
        q[4] = -1.5 * q[0]
        solution = solve_ivp(
            rhs,
            [0, end],
            q.ravel(),
            args=(Encounter(),),
            method="DOP853",
            rtol=2e-10,
            atol=2e-12,
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        units = np.array([HILL, HILL, HILL, OMEGA * HILL, OMEGA * HILL, OMEGA * HILL])
        return solution.y[:, -1].reshape(6, 2) * units[:, None]

    def difference(y):
        q = final(y)
        return q[:2, 0] - q[:2, 1]

    guess = 1.5 * impact * end - 450_000.0
    solution = root(difference, guess, tol=1e-8)
    residual = float(np.linalg.norm(difference(solution.x)))
    if not solution.success or residual > 0.01 or np.any(solution.x < 100_000):
        raise RuntimeError("No valid upstream stream meeting found")
    q = final(solution.x)
    free = np.column_stack((impact, solution.x - 1.5 * impact * end))
    centers = tuple((float(x), float(y)) for x, y in zip(impact, solution.x))
    return centers, {
        "initial_centres_xy_m": centers,
        "predicted_meeting_orbits": event_orbits,
        "predicted_meeting_xy_m": q[:2].mean(axis=1).tolist(),
        "point_track_residual_m": residual,
        "predicted_relative_speed_m_s": float(np.linalg.norm(q[3:, 0] - q[3:, 1])),
        "no_moon_point_track_separation_m": float(np.linalg.norm(free[0] - free[1])),
        "selection": "Two initially circular upstream regions in the same outer strip, selected from gravity-predicted crossing orbits. No post-encounter velocity is supplied.",
    }


def animate_pair(out, count, comparison="moon"):
    if comparison not in {"moon", "inclination", "streams"}:
        raise ValueError("Comparison must be moon, inclination, or streams")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import PillowWriter

    datasets = []
    metadata = []
    names = (
        [f"full_{count}", f"no_moon_{count}"]
        if comparison == "moon"
        else (
            [f"inclined_{count}", f"full_{count}"]
            if comparison == "inclination"
            else [f"crossing_full_{count}", f"crossing_no_moon_{count}"]
        )
    )
    for name in names:
        datasets.append(np.load(out / f"{name}.npz"))
        metadata.append(json.loads((out / f"{name}.json").read_text()))
    a, b = datasets
    if not np.array_equal(a["ids"], b["ids"]) or not np.array_equal(
        a["initial_state"], b["initial_state"]
    ):
        raise RuntimeError("The paired runs must start with identical real grains")
    curves = [
        np.loadtxt(out / f"{name}.csv", delimiter=",", skiprows=1) for name in names
    ]
    core = a["core"]
    side = a["side"]
    outer = core & (side == 1)
    colors = ["#60d2ff", "#ffaf65"]
    labels = (
        ["Moon gravity on", "Moon gravity off"]
        if comparison != "inclination"
        else ["Tilted moon", "Coplanar moon"]
    )
    positions = []
    centers = []
    for data in datasets:
        q = data["state"]
        center = q[:, outer, :3].mean(axis=1)
        positions.append(
            np.column_stack(
                (
                    q[:, outer, 0].ravel() - np.repeat(center[:, 0], outer.sum()),
                    q[:, outer, 2].ravel(),
                )
            ).reshape(len(q), outer.sum(), 2)
        )
        centers.append(
            np.stack(
                [q[:, core & (side == sign), :3].mean(axis=1) for sign in [1, -1]],
                axis=1,
            )
        )
    limit_x = max(10, float(np.max([abs(p[:, :, 0]).max() for p in positions])) * 1.1)
    limit_z = max(15, float(np.max([abs(p[:, :, 1]).max() for p in positions])) * 1.1)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.8))
    fig.subplots_adjust(
        left=0.09, right=0.98, top=0.86, bottom=0.16, hspace=0.48, wspace=0.28
    )
    fig.patch.set_facecolor("#071726")
    for ax in axes.flat:
        ax.set_facecolor("#071726")
        ax.tick_params(colors="#dbe8f1")
        ax.xaxis.label.set_color("#dbe8f1")
        ax.yaxis.label.set_color("#dbe8f1")
        for sp in ax.spines.values():
            sp.set_color("#536b7b")
        ax.grid(color="#254252", linewidth=0.4)
    fig.text(
        0.5,
        0.97,
        f"{count:,} physical grains; grain radius = {metadata[0]['config']['radius_m']:g} m",
        ha="center",
        color="white",
        fontsize=15,
    )
    clock = fig.text(0.5, 0.92, "", ha="center", color="white", fontsize=10)
    overview, local, orbit, height = axes.flat
    overview_span = max(
        800, 100 * np.ceil(max(abs(c[:, :, 1]).max() for c in centers) / 1000 / 100)
    )
    overview.set(
        xlim=(-overview_span, overview_span),
        ylim=(-28, 28),
        xlabel="Along orbit (km)",
        ylabel="Radial offset (km)",
    )
    overview.set_title(
        "Centres of the tracked grain packets", color="white", fontsize=10
    )
    overview.scatter([0], [0], s=40, color="#ffd166")
    overview.annotate(
        "Moon reference",
        (0, 0),
        xytext=(6, 6),
        textcoords="offset points",
        color="#ffd166",
        fontsize=8,
    )
    local.set(
        xlim=(-limit_x, limit_x),
        ylim=(-limit_z, limit_z),
        xlabel="Radial distance from packet centre (m)",
        ylabel="Height above ring plane (m)",
    )
    local.set_title(
        f"Side view: same {outer.sum()} packet-A core grains",
        color="white",
        fontsize=10,
    )
    local.axhline(0, color="#8398a6", linewidth=0.5)
    times = curves[0][:, 0] - metadata[0]["config"]["upstream_orbits"]
    for ax in [orbit, height]:
        ax.set_xlim(times[0], times[-1])
        ax.set_xlabel("Orbits relative to the reference event")
        ax.axvline(0, color="#a7b6c1", linestyle=":", linewidth=0.8)
    signs = [np.sign(a["initial_state"][side == value, 0].mean()) for value in [1, -1]]
    guiding = (
        [0.5 * (signs[0] * c[:, 1] + signs[1] * c[:, 2]) for c in curves]
        if comparison != "inclination"
        else [d["state"][:, outer, 2].mean(axis=1) for d in datasets]
    )
    extrema = np.concatenate(guiding)
    orbit.set_ylim(min(-5, float(extrema.min()) - 5), max(10, float(extrema.max()) + 5))
    orbit.set_ylabel("Mean orbit shift away from moon orbit (m)")
    orbit.set_title(
        "Orbit estimate (transient near encounter)", color="white", fontsize=10
    )
    thickness = (
        [c[:, 3] for c in curves]
        if comparison != "inclination"
        else [np.std(d["state"][:, outer, 2], axis=1) for d in datasets]
    )
    if comparison == "inclination":
        orbit.set_ylabel("Mean outer packet height (m)")
        orbit.set_title("Coherent vertical lift", color="white", fontsize=10)
    height.set_ylim(0, max(float(h.max()) for h in thickness) * 1.2)
    height.set_ylabel("Core RMS height spread (m)")
    height.set_title(
        "Vertical motion in the same core grains", color="white", fontsize=10
    )
    dots = []
    trails = []
    orbit_lines = []
    height_lines = []
    for color, label in zip(colors, labels):
        dots.append(local.scatter([], [], s=8, color=color, alpha=0.8, label=label))
        trails.append(
            [overview.plot([], [], color=color, linewidth=1)[0] for _ in [1, -1]]
        )
        orbit_lines.append(orbit.plot([], [], color=color, label=label)[0])
        height_lines.append(height.plot([], [], color=color, label=label)[0])
    local.legend(fontsize=7, loc="upper right", facecolor="#203240", labelcolor="white")
    height.legend(fontsize=7, facecolor="#203240", labelcolor="white")
    fig.text(
        0.5,
        0.075,
        "3D dynamics; two projected views. Packets have free edges. No imposed waves or added heating impulse.",
        ha="center",
        color="white",
        fontsize=9,
    )
    fig.text(
        0.5,
        0.035,
        "Count adds neighboring grains; the tracked core stays fixed. Axes use different scales. Lighting is not modeled.",
        ha="center",
        color="white",
        fontsize=9,
    )
    movie_name = (
        f"coupled_{count}"
        if comparison == "moon"
        else f"inclination_{count}"
        if comparison == "inclination"
        else f"crossing_{count}"
    )
    event_name = (
        "predicted stream meeting"
        if comparison == "streams"
        else "nominal moon encounter"
    )
    with PillowWriter(fps=10).saving(
        fig, str(out / f"{movie_name}.gif"), dpi=90
    ) as writer:
        for frame in range(0, len(times), 2):
            for k in [0, 1]:
                dots[k].set_offsets(positions[k][frame])
                for group, line in enumerate(trails[k]):
                    line.set_data(
                        centers[k][: frame + 1, group, 1] / 1000,
                        centers[k][: frame + 1, group, 0] / 1000,
                    )
                orbit_lines[k].set_data(times[: frame + 1], guiding[k][: frame + 1])
                height_lines[k].set_data(times[: frame + 1], thickness[k][: frame + 1])
            clock.set_text(
                f"Time: {times[frame]:+.2f} orbits from {event_name} | Contacts: {int(curves[0][frame, 7])} / {int(curves[1][frame, 7])}"
            )
            writer.grab_frame()
            if frame == 120:
                fig.savefig(
                    out / f"{movie_name}.png",
                    dpi=120,
                    facecolor=fig.get_facecolor(),
                )
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--counts", nargs="+", type=int, default=[128, 288, 512])
    parser.add_argument("--controls", action="store_true")
    parser.add_argument(
        "--inclined",
        action="store_true",
        help="Also test a documented 8.6 km moon vertical amplitude at the smallest count",
    )
    args = parser.parse_args()
    if args.counts != sorted(set(args.counts)):
        parser.error("Counts must be increasing and distinct")
    for n in args.counts:
        Packet(count=n)
        if n < 128:
            parser.error(
                "The series tracks the same 128 core grains; use at least 128 total"
            )
    if args.out.exists() and any(args.out.iterdir()):
        parser.error("Use a new output directory")
    args.out.mkdir(parents=True, exist_ok=True)
    identity = {
        "commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    for n in args.counts:
        cfg = Packet(count=n)
        save_case(args.out, f"full_{n}", cfg)
        save_case(args.out, f"no_moon_{n}", replace(cfg, moon_gravity=False))
        animate_pair(args.out, n)
    if args.inclined:
        tilted = Packet(count=args.counts[0], moon_vertical_amplitude_m=8600)
        save_case(args.out, f"inclined_{tilted.count}", tilted)
        save_case(
            args.out,
            f"inclined_half_step_{tilted.count}",
            replace(tilted, dt_s=tilted.dt_s / 2),
        )
        animate_pair(args.out, tilted.count, comparison="inclination")
    if args.controls:
        cfg = Packet(count=args.counts[0])
        for name, c in [
            ("half_step", replace(cfg, dt_s=cfg.dt_s / 2)),
            ("no_grain_gravity", replace(cfg, self_gravity=False)),
            ("force_only", replace(cfg, self_gravity=False, collisions=False)),
        ]:
            save_case(args.out, f"{name}_{cfg.count}", c)
    (args.out / "manifest.json").write_text(json.dumps(identity, indent=2) + "\n")


if __name__ == "__main__":
    main()
