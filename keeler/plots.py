"""Rebuild figures and small review tables from saved numerical outputs."""

import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .encounter import HILL
from .optics import slab


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    d = args.directory
    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    c = np.load(d / "circular.npz")
    t, b, q = c["t"], c["impact_m"], c["state_hill"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for name, label in [
        ("circular", "Full moon mass"),
        ("half_mass", "Half moon mass"),
    ]:
        x = np.load(d / f"{name}.npz")
        axes[0, 0].plot(b / 1000, (x["guiding_m"][-1] - b) / 1000, label=label)
    axes[0, 0].set(
        xlabel="Initial offset from moon (km)",
        ylabel="Guiding-center shift (km)",
        title="Collisionless scattering after six orbits",
    )
    axes[0, 0].legend()
    for impact in (13e3, 15e3, 20e3, 21e3):
        j = np.argmin(abs(b - impact))
        mask = (t > 0) & (t < 12 * np.pi)
        axes[0, 1].plot(
            -q[mask, j, 1] * HILL / 1000,
            q[mask, j, 0] * HILL / 1000,
            label=f"{impact / 1000:g} km",
        )
        axes[1, 0].plot(
            t / (2 * np.pi),
            (c["guiding_m"][:, j] - b[j]) / 1000,
            label=f"{impact / 1000:g} km",
        )
    axes[0, 1].set(
        xlabel="Downstream distance (km)",
        ylabel="Radial offset (km)",
        title="Test-particle tracks; crossings are not resolved",
    )
    axes[0, 1].legend()
    axes[1, 0].set(
        xlim=(-1, 6),
        xlabel="Time / orbital period (nominal encounter at zero)",
        ylabel="Guiding-center shift (km)",
        title="Guiding center differs from radial oscillation",
    )
    phase = np.loadtxt(d / "phase_sweep.csv", delimiter=",", skiprows=1)
    for impact in np.unique(phase[:, 1]):
        rows = phase[phase[:, 1] == impact]
        axes[1, 1].plot(
            rows[:, 0] * 180 / np.pi,
            rows[:, 4] / 1000,
            ".-",
            label=f"{impact / 1000:g} km",
        )
    axes[1, 1].axhspan(
        0.5, 1.5, alpha=0.13, color="gray", label="Reported height range (not RMS)"
    )
    axes[1, 1].set(
        xlabel="Assumed moon epicycle phase (degrees)",
        ylabel="Final test-particle vertical amplitude (km)",
        title="Eccentric + inclined moon; no collisions",
    )
    axes[1, 1].legend(fontsize=8)
    fig.savefig(d / "encounters.png", dpi=160)
    plt.close(fig)
    names = [
        ("no_collisions", "No collisions"),
        ("inelastic", "Inelastic e=0.5"),
        ("elastic", "Elastic e=1"),
        ("self_gravity", "Inelastic + self-gravity"),
        ("bridges", "Velocity-dependent ice law"),
    ]
    if all((d / f"box_{name}.csv").exists() for name, _ in names):
        fig, ax = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
        for name, label in names:
            a = np.loadtxt(d / f"box_{name}.csv", delimiter=",", skiprows=1)
            ax[0].plot(a[:, 0], a[:, 1], label=label)
            ax[1].plot(a[:, 0], a[:, 6], label=label)
        ax[0].set(
            xlabel="Orbits",
            ylabel="RMS height (m)",
            title="Physical grains in an excited, homogeneous patch",
        )
        ax[0].legend(fontsize=8)
        ax[1].set(
            xlabel="Orbits",
            ylabel="Vertical oscillator energy (J/kg)",
            yscale="log",
            title="Energy transfer and loss; no imposed damping",
        )
        fig.savefig(d / "collisions.png", dpi=160)
        plt.close(fig)
    names = [
        ("radius_1m", "1 m radius"),
        ("inelastic", "2 m radius"),
        ("radius_4m", "4 m radius"),
        ("half_step", "2 m; half time step"),
        ("larger_box", "2 m; double box width"),
    ]
    if all((d / f"box_{name}.csv").exists() for name, _ in names):
        fig, ax = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
        for name, label in names:
            a = np.loadtxt(d / f"box_{name}.csv", delimiter=",", skiprows=1)
            ax[0].plot(a[:, 0], a[:, 1], label=label)
        for name in ["inelastic", "seed_43", "seed_44"]:
            if (d / f"box_{name}.csv").exists():
                a = np.loadtxt(d / f"box_{name}.csv", delimiter=",", skiprows=1)
                ax[1].plot(a[:, 0], a[:, 1], label=name)
        ax[0].set(
            xlabel="Orbits",
            ylabel="RMS height (m)",
            title="Fixed surface density; optical depth changes with radius",
        )
        ax[1].set(
            xlabel="Orbits",
            ylabel="RMS height (m)",
            title="Finite sample variation in the 2 m run",
        )
        ax[0].legend(fontsize=8)
        ax[1].legend(fontsize=8)
        fig.savefig(d / "sensitivity.png", dpi=160)
        plt.close(fig)
    energy_cases = [
        ("encounter_energy", "1 s step"),
        ("encounter_energy_half_step", "0.5 s step"),
    ]
    if all((d / f"box_{name}.csv").exists() for name, _ in energy_cases):
        fig, ax = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
        for name, label in energy_cases:
            a = np.loadtxt(d / f"box_{name}.csv", delimiter=",", skiprows=1)
            ax[0].plot(a[:, 0], a[:, 1] / 1000, label=label)
            ax[1].plot(a[:, 0], a[:, 2] / 1000, label=label)
        ax[0].set(
            xlabel="Orbits",
            ylabel="RMS height (km)",
            title="Imposed radial dispersion 0.09 m/s; Bridges law",
        )
        ax[1].axhspan(0.5, 1.5, alpha=0.13, color="gray", label="Reported wall heights")
        ax[1].set(
            xlabel="Orbits",
            ylabel="99th percentile of absolute altitude (km)",
            title="Particle percentile; no optical visibility model",
        )
        ax[0].legend(fontsize=8)
        ax[1].legend(fontsize=8)
        fig.savefig(d / "height_scale.png", dpi=160)
        plt.close(fig)
        grid = np.linspace(2, 3, 401)
        means = {}
        for name, _ in energy_cases:
            a = np.loadtxt(d / f"box_{name}.csv", delimiter=",", skiprows=1)
            means[name] = {
                key: float(np.trapezoid(np.interp(grid, a[:, 0], a[:, column]), grid))
                for key, column in [("rms_height_m", 1), ("vertical_energy_J_kg", 6)]
            }
        comparison = {
            "interval_orbits": [2, 3],
            "time_means": means,
            "percent_change_when_halving_step": {
                key: 100
                * (
                    means["encounter_energy_half_step"][key]
                    / means["encounter_energy"][key]
                    - 1
                )
                for key in means["encounter_energy"]
            },
        }
        (d / "energy_step_comparison.json").write_text(
            json.dumps(comparison, indent=2, allow_nan=False) + "\n"
        )

    tau = np.linspace(0, 1.5, 300)
    f = 0.2
    mean_tau = 0.5
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    ax[0].plot(tau, slab(tau))
    ax[0].set(
        xlabel="Normal optical depth",
        ylabel="Relative slab radiance",
        title="Single scattering; fixed geometry and grain optics",
    )
    ax[1].bar(
        ["Clumps + empty space", "Uniform same mass"],
        [f * slab(mean_tau / f), slab(mean_tau)],
        color=["#777777", "#297eb5"],
    )
    ax[1].set(
        ylabel="Relative slab radiance",
        title="Illustration: spreading opaque clumps can brighten material",
    )
    fig.savefig(d / "optics.png", dpi=160)
    plt.close(fig)
    rows = []
    for file in sorted(d.glob("box_*.json")):
        a = np.loadtxt(file.with_suffix(".csv"), delimiter=",", skiprows=1)
        meta = json.loads(file.read_text())
        peak = int(np.argmax(a[:, 1]))
        rows.append(
            [
                file.stem.removeprefix("box_"),
                meta["particles"],
                meta["tau"],
                a[peak, 0],
                a[peak, 1],
                a[-1, 1],
                a[-1, 2],
                a[-1, 6],
                meta["seconds"],
            ]
        )
    header = "case,particles,tau,peak_orbits,peak_rms_height_m,final_rms_height_m,final_abs_p99_m,final_vertical_energy_J_kg,seconds"
    (d / "box_summary.csv").write_text(
        header + "\n" + "\n".join(",".join(map(str, r)) for r in rows) + "\n"
    )
    j = [np.argmin(abs(b - x)) for x in [13e3, 14e3, 15e3, 20e3, 21e3, 25e3]]
    np.savetxt(
        d / "encounter_summary.csv",
        np.column_stack(
            (b[j], c["guiding_m"][-1, j] - b[j], c["radial_epicycle_m"][-1, j])
        ),
        delimiter=",",
        header="impact_m,guiding_shift_m,radial_epicycle_m",
        comments="",
    )


if __name__ == "__main__":
    main()
