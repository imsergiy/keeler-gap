import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time
import numpy as np
from .encounter import (
    SurfaceImpact,
    Encounter,
    integrate,
    elements,
    jacobi,
    crowding,
    HILL,
    PERIOD,
    OMEGA,
)
from .grains import Box, run


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def encounters(out):
    impact = np.linspace(13_000, 35_000, 89)
    summary = {}
    for name, config in {
        "circular": Encounter(),
        "inclined": Encounter(vertical_amplitude_m=8600),
        "eccentric_inclined": Encounter(
            radial_epicycle_m=4500, vertical_amplitude_m=8600
        ),
        "half_mass": Encounter(mass_factor=0.5),
        "no_moon": Encounter(mass_factor=0),
    }.items():
        t, q = integrate(impact, config)
        guiding, radial, vertical = elements(q)
        stationary = config.radial_epicycle_m == config.vertical_amplitude_m == 0
        density_fields = {}
        if stationary:
            density, det = crowding(t, impact, q)
            density_fields = {"crowding": density, "determinant": det}
        np.savez_compressed(
            out / f"{name}.npz",
            t=t,
            impact_m=impact,
            state_hill=q,
            guiding_m=guiding * HILL,
            radial_epicycle_m=radial * HILL,
            vertical_amplitude_m=vertical * HILL,
            **density_fields,
        )
        j = 28
        summary[name] = {
            "config": asdict(config),
            "impact_m": float(impact[j]),
            "guiding_shift_m": float((guiding[-1, j] - guiding[0, j]) * HILL),
            "radial_epicycle_m": float(radial[-1, j] * HILL),
            "vertical_amplitude_m": float(vertical[-1, j] * HILL),
            "fold_present": bool(np.any(det > 0)) if stationary else None,
        }
        if config.radial_epicycle_m == config.vertical_amplitude_m == 0:
            c = jacobi(q, config.mass_factor)
            summary[name]["max_jacobi_absolute_drift"] = float(np.max(abs(c - c[0])))
        print(name, summary[name], flush=True)
    phases = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    phase_b = np.array([13e3, 15e3, 20e3, 21e3])
    phase_rows = []
    for phase in phases:
        cfg = Encounter(
            radial_epicycle_m=4500,
            vertical_amplitude_m=8600,
            phase=float(phase),
            samples=501,
        )
        for impact_value in phase_b:
            try:
                _, q = integrate([impact_value], cfg)
            except SurfaceImpact:
                phase_rows.append([phase, impact_value, np.nan, np.nan, np.nan, 1])
                continue
            g, a, z = elements(q)
            phase_rows.append(
                [
                    phase,
                    impact_value,
                    (g[-1, 0] - g[0, 0]) * HILL,
                    a[-1, 0] * HILL,
                    z[-1, 0] * HILL,
                    0,
                ]
            )
    np.savetxt(
        out / "phase_sweep.csv",
        phase_rows,
        delimiter=",",
        header="phase_rad,impact_m,guiding_shift_m,radial_epicycle_m,vertical_amplitude_m,surface_impact_1_else_0",
        comments="",
    )
    checks = {}
    b = np.array([15e3, 20e3, 25e3])
    for tol in [2e-8, 2e-10, 2e-12]:
        _, q = integrate(b, Encounter(rtol=tol, samples=501))
        checks[str(tol)] = (elements(q)[0][-1] * HILL - b).tolist()
    for start in [-30.0, -60.0, -120.0, -240.0, -480.0]:
        _, q = integrate(b, Encounter(start=start, samples=501))
        checks[f"start_{start}"] = (elements(q)[0][-1] * HILL - b).tolist()
    for end in [12 * np.pi, 24 * np.pi, 48 * np.pi]:
        _, q = integrate(b, Encounter(start=-480, end=end, samples=501))
        checks[f"end_{end}"] = (elements(q)[0][-1] * HILL - b).tolist()
    write_json(out / "encounters.json", summary)
    write_json(out / "encounter_convergence.json", checks)


def boxes(out, extended, selected):
    baseline = Box()
    configs = {
        "no_collisions": replace(baseline, collisions=False),
        "inelastic": baseline,
        "elastic": replace(baseline, restitution=1, dt_s=0.5),
        "self_gravity": replace(baseline, self_gravity=True),
        "radius_1m": replace(baseline, radius_m=1, dt_s=1),
        "radius_4m": replace(baseline, radius_m=4),
        "half_step": replace(baseline, dt_s=1),
        "larger_box": replace(baseline, size_m=256),
        "seed_43": replace(baseline, seed=43),
        "seed_44": replace(baseline, seed=44),
        "bridges": replace(baseline, restitution="bridges"),
        "encounter_energy": replace(
            baseline, restitution="bridges", horizontal_sigma=0.09, dt_s=1
        ),
    }
    if extended:
        configs.update(
            {
                "gravity_ghosts2": replace(
                    baseline, self_gravity=True, ghosts=2, orbits=0.5
                ),
                "gravity_tighter": replace(
                    baseline, self_gravity=True, opening_angle=0.15, orbits=0.5
                ),
                "gravity_half_step": replace(
                    baseline, self_gravity=True, dt_s=1, orbits=0.5
                ),
                "radius_4m_large": replace(baseline, radius_m=4, size_m=256),
                "radius_1m_seed43": replace(baseline, radius_m=1, dt_s=1, seed=43),
                "radius_4m_seed43": replace(baseline, radius_m=4, seed=43),
                "encounter_energy_half_step": replace(
                    baseline, restitution="bridges", horizontal_sigma=0.09, dt_s=0.5
                ),
            }
        )
    if selected:
        unknown = set(selected) - set(configs)
        if unknown:
            raise ValueError(f"Unknown box cases: {sorted(unknown)}")
        configs = {k: v for k, v in configs.items() if k in selected}
    for name, c in configs.items():
        target = out / f"box_{name}.json"
        if target.exists():
            raise FileExistsError(f"{target}: use a new output directory")
        print("Running", name, "N=", c.count, flush=True)
        rows, snapshots, meta = run(c, progress=True)
        np.savetxt(
            out / f"box_{name}.csv",
            rows,
            delimiter=",",
            header=",".join(meta["columns"]),
            comments="",
        )
        np.savez_compressed(out / f"box_{name}.npz", snapshots=snapshots)
        write_json(target, meta)
        print("Done", name, round(meta["seconds"], 2), "s", flush=True)


def budget(args):
    area = args.width_km * args.length_km * 1e6
    n = args.tau * area / (np.pi * args.radius_m**2)
    print(
        json.dumps(
            {
                "area_m2": area,
                "particles": int(np.ceil(n)),
                "particle_array_lower_bound_GiB": n * 128 / 2**30,
                "note": "Tree, collision, ghost and output storage are additional. No runtime prediction.",
            },
            indent=2,
        )
    )


def main():
    p = argparse.ArgumentParser(
        description="Daphnis encounters and local physical-grain controls"
    )
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("encounters", "boxes", "all"):
        s = sub.add_parser(name)
        s.add_argument("--out", type=Path, required=True)
        if name != "encounters":
            s.add_argument("--extended", action="store_true")
            s.add_argument(
                "--case",
                action="append",
                help="Run only a named box case; repeat to select several",
            )
    s = sub.add_parser("budget")
    s.add_argument("--radius-m", type=float, default=1)
    s.add_argument("--tau", type=float, default=0.5)
    s.add_argument("--width-km", type=float, default=10)
    s.add_argument("--length-km", type=float, default=1000)
    args = p.parse_args()
    if args.command == "budget":
        if not all(
            np.isfinite(x) and x > 0
            for x in [args.radius_m, args.tau, args.width_km, args.length_km]
        ):
            p.error("Budget inputs must be finite and positive")
        budget(args)
        return
    sha = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain"], capture_output=True, text=True, check=True
    ).stdout
    source_hashes = {
        str(f): hashlib.sha256(f.read_bytes()).hexdigest()
        for f in sorted(Path("keeler").glob("*.py"))
    }
    args.out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    if args.command in ("encounters", "all"):
        encounters(args.out)
    if args.command in ("boxes", "all"):
        boxes(args.out, args.extended, args.case)
    write_json(
        args.out / "manifest.json",
        {
            "commit": sha,
            "working_tree_dirty": bool(dirty),
            "python": platform.python_version(),
            "source_sha256": source_hashes,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "seconds": time.monotonic() - start,
            "hill_radius_m": HILL,
            "omega_per_s": OMEGA,
            "period_s": PERIOD,
            "files": {
                f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                for f in sorted(args.out.iterdir())
                if f.is_file() and f.name != "manifest.json"
            },
        },
    )


if __name__ == "__main__":
    main()
