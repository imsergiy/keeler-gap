"""Refine the sampling of a finite ring packet; all markers are independent."""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
import numpy as np
from scipy.stats import qmc

from .encounter import Encounter, HILL, PERIOD, elements, integrate, jacobi

BMIN = 21_000.0
BMAX = 31_000.0
MOVIE_ORBITS = np.linspace(0, 6, 97)
PROFILE_FRAME = 64
PROFILE_BINS_M = np.arange(100_000, 900_001, 40_000)


def packet_sample(count):
    """One side; reflection supplies the circular moon's identical other side."""
    u = qmc.Halton(2, scramble=False).random(count)
    impact = np.sqrt(BMIN**2 + u[:, 0] * (BMAX**2 - BMIN**2))
    # The incoming mass flux is proportional to b, so uniform time gives uniform area density.
    initial_age = (u[:, 1] - 1) * 8 * np.pi
    return impact, initial_age


def interpolate_positions(t, q, ages):
    if ages.shape[1] != q.shape[1] or ages.min() < t[0] or ages.max() > t[-1]:
        raise ValueError("Requested movie ages must lie within the integrated tracks")
    positions = np.empty((*ages.shape, 2))
    for j in range(q.shape[1]):
        for k in range(2):
            positions[:, j, k] = np.interp(ages[:, j], t, q[:, j, k]) * HILL
    return positions


def new_tracks(impact, initial_age, config, batch_size=128):
    positions = np.empty((len(MOVIE_ORBITS), len(impact), 2))
    shifts = np.empty(len(impact))
    max_jacobi_drift = 0.0
    for start in range(0, len(impact), batch_size):
        stop = min(start + batch_size, len(impact))
        t, q = integrate(impact[start:stop], config)
        ages = 2 * np.pi * MOVIE_ORBITS[:, None] + initial_age[None, start:stop]
        positions[:, start:stop] = interpolate_positions(t, q, ages)
        guiding = elements(q[[0, -1]])[0]
        shifts[start:stop] = (guiding[-1] - guiding[0]) * HILL
        invariant = jacobi(q)
        max_jacobi_drift = max(
            max_jacobi_drift, float(np.max(abs(invariant - invariant[0])))
        )
        print(f"  integrated {stop}/{len(impact)} new outer-side tracks", flush=True)
    return positions, shifts, max_jacobi_drift


def profile(positive_positions):
    downstream = -positive_positions[:, 1]
    radial = positive_positions[:, 0]
    counts, edge = [], []
    for lo, hi in zip(PROFILE_BINS_M[:-1], PROFILE_BINS_M[1:]):
        selected = radial[(downstream >= lo) & (downstream < hi)]
        counts.append(len(selected))
        edge.append(float(np.quantile(selected, 0.1)) if len(selected) else np.nan)
    return np.asarray(counts), np.asarray(edge)


def animate(path, positive, impact, initial_age, count):
    plt.rcParams.update({"font.size": 10})
    fig, axes = plt.subplots(2, 1, figsize=(10, 6.8), sharex=True, sharey=True)
    fig.subplots_adjust(top=0.88, bottom=0.14, left=0.09, right=0.98, hspace=0.23)
    fig.patch.set_facecolor("#071726")
    fig.text(
        0.5,
        0.96,
        f"{count:,} independent orbit markers",
        ha="center",
        color="white",
        fontsize=15,
    )
    clock = fig.text(0.5, 0.915, "", ha="center", color="white")
    dots = []
    for ax, title in zip(
        axes, ["Saturn gravity only — comparison", "Saturn + Daphnis gravity"]
    ):
        ax.set_facecolor("#071726")
        ax.set_title(title, color="white", fontsize=11)
        ax.set(xlim=(-1800, 1800), ylim=(-36, 36), ylabel="Radial offset (km)")
        ax.tick_params(colors="#d5e5ef")
        ax.yaxis.label.set_color("#d5e5ef")
        for spine in ax.spines.values():
            spine.set_color("#456173")
        ax.grid(color="#244050", linewidth=0.4)
        dots.append(
            ax.scatter([], [], s=3, color="#79d5ff", alpha=0.8, edgecolors="none")
        )
        ax.axhline(0, color="#456173", linewidth=0.5)
        ax.annotate(
            "",
            xy=(-700, 32),
            xytext=(-300, 32),
            arrowprops={"arrowstyle": "->", "color": "#d5e5ef"},
        )
        ax.annotate(
            "",
            xy=(700, -32),
            xytext=(300, -32),
            arrowprops={"arrowstyle": "->", "color": "#d5e5ef"},
        )
    axes[0].scatter([0], [0], marker="+", s=30, color="#8b9aa5")
    axes[1].scatter([0], [0], s=55, color="#ffd166", edgecolor="none", zorder=5)
    axes[1].annotate(
        "Daphnis", (0, 0), xytext=(8, 8), textcoords="offset points", color="#ffd166"
    )
    axes[1].set_xlabel("Along orbit relative to the moon (km)", color="#d5e5ef")
    fig.text(
        0.5,
        0.065,
        "View follows the moon. All markers orbit Saturn. No collisions or mutual particle gravity.",
        ha="center",
        color="white",
        fontsize=9,
    )
    fig.text(
        0.5,
        0.03,
        "Radial direction stretched. Dots and moon symbol are not physical sizes or brightness.",
        ha="center",
        color="white",
        fontsize=9,
    )
    writer = PillowWriter(fps=10)
    with writer.saving(fig, str(path), dpi=90):
        for i, orbit in enumerate(MOVIE_ORBITS):
            age = 2 * np.pi * orbit + initial_age
            free = np.column_stack((impact, -1.5 * impact * age))
            for artist, p in zip(dots, [free, positive[i]]):
                paired = np.concatenate((p, -p))
                artist.set_offsets(paired[:, [1, 0]] / 1000)
            clock.set_text(
                f"Elapsed time: {orbit:.2f} Saturn orbits ({orbit * PERIOD / 3600:.1f} hours)"
            )
            writer.grab_frame()
            if i == PROFILE_FRAME:
                fig.savefig(
                    path.with_suffix(".png"), dpi=120, facecolor=fig.get_facecolor()
                )
    plt.close(fig)


def interpolation_check(out):
    impact = np.array([BMIN, 25_000.0, BMAX])
    initial_age = np.array([-7.0, -15.0, -22.0])
    ages = 2 * np.pi * MOVIE_ORBITS[:, None] + initial_age
    values = []
    for samples in (2401, 4801):
        t, q = integrate(impact, Encounter(samples=samples))
        values.append(interpolate_positions(t, q, ages))
    result = {
        "method": "Same three tracks, 2401 versus 4801 interpolation samples",
        "max_position_difference_m": float(
            np.max(np.linalg.norm(values[1] - values[0], axis=-1))
        ),
    }
    (out / "interpolation_check.json").write_text(json.dumps(result, indent=2) + "\n")
    if result["max_position_difference_m"] > 10:
        raise RuntimeError("Movie interpolation error exceeds 10 m")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--counts", type=int, nargs="+", default=[250, 1000, 4000, 16000]
    )
    args = parser.parse_args()
    if args.counts != sorted(set(args.counts)) or any(
        n < 10 or n % 2 or n > 64000 for n in args.counts
    ):
        parser.error(
            "Counts must be distinct, increasing, even integers from 10 to 64000"
        )
    if args.out.exists() and any(args.out.iterdir()):
        parser.error("Use a new output directory")
    args.out.mkdir(parents=True, exist_ok=True)
    config = Encounter(samples=2401)
    sample_b, sample_age = packet_sample(args.counts[-1] // 2)
    all_positions = np.empty((len(MOVIE_ORBITS), 0, 2))
    all_shifts = np.empty(0)
    summaries = []
    previous_profile = None
    previous_mean = None
    source_identity = {
        "commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "working_tree_dirty": bool(
            subprocess.check_output(["git", "status", "--porcelain"], text=True)
        ),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    started = time.monotonic()
    for count in args.counts:
        stage = time.monotonic()
        old = all_positions.shape[1]
        new = count // 2
        print(
            f"RUN {count} markers; adding {new - old} outer tracks and their reflected inner tracks",
            flush=True,
        )
        positions, shifts, drift = new_tracks(
            sample_b[old:new], sample_age[old:new], config
        )
        all_positions = np.concatenate((all_positions, positions), axis=1)
        all_shifts = np.concatenate((all_shifts, shifts))
        cells, edge = profile(all_positions[PROFILE_FRAME])
        valid = np.isfinite(edge)
        change = None
        if previous_profile is not None:
            common = valid & np.isfinite(previous_profile)
            change = (
                float(np.sqrt(np.mean((edge[common] - previous_profile[common]) ** 2)))
                if common.any()
                else None
            )
        mean = float(np.mean(all_shifts))
        summary = {
            "markers": count,
            "independent_outer_tracks": new,
            "mean_outward_guiding_shift_m": mean,
            "shift_epoch": "Six orbits after each track nominal conjunction, not one common movie time",
            "profile_movie_orbits": 4,
            "mean_shift_change_percent": None
            if previous_mean is None
            else float(100 * abs(mean / previous_mean - 1)),
            "edge_profile_rms_change_m": change,
            "occupied_profile_bins": int(valid.sum()),
            "minimum_markers_per_profile_bin": int(cells.min()),
            "max_new_track_jacobi_absolute_drift": drift,
            "gif": f"wave_{count}.gif",
            "model": "Circular coplanar moon; independent massless orbit markers",
            "config": asdict(config),
        }
        np.savetxt(
            args.out / f"profile_{count}.csv",
            np.column_stack(
                ((PROFILE_BINS_M[:-1] + PROFILE_BINS_M[1:]) / 2, cells, edge)
            ),
            delimiter=",",
            header="downstream_bin_center_m,outer_markers,radial_position_percentile_10_m",
            comments="",
        )
        np.savez_compressed(
            args.out / f"markers_{count}.npz",
            orbits=MOVIE_ORBITS,
            outer_xy_m=all_positions,
            impact_m=sample_b[:new],
            initial_age=sample_age[:new],
            outward_shift_m=all_shifts,
        )
        animate(
            args.out / f"wave_{count}.gif",
            all_positions,
            sample_b[:new],
            sample_age[:new],
            count,
        )
        summary["elapsed_seconds"] = time.monotonic() - stage
        summaries.append(summary)
        (args.out / f"wave_{count}.json").write_text(
            json.dumps(summary, indent=2, allow_nan=False) + "\n"
        )
        (args.out / "resolution.json").write_text(
            json.dumps(summaries, indent=2, allow_nan=False) + "\n"
        )
        previous_profile = edge
        previous_mean = mean
        print("DONE", count, json.dumps(summary), flush=True)
    interpolation = interpolation_check(args.out)
    fig, ax = plt.subplots(figsize=(9, 4), layout="constrained")
    for n in args.counts:
        a = np.loadtxt(args.out / f"profile_{n}.csv", delimiter=",", skiprows=1)
        ax.plot(a[:, 0] / 1000, a[:, 2] / 1000, ".-", label=f"{n:,} markers")
    ax.set(
        xlabel="Outer-ring downstream distance (km)",
        ylabel="Radial position, lower 10th percentile (km)",
        title="Same packet after four orbits; fixed 40 km bins",
    )
    ax.legend()
    fig.savefig(args.out / "profiles.png", dpi=160)
    plt.close(fig)
    manifest = {
        **source_identity,
        "elapsed_seconds": time.monotonic() - started,
        "interpolation_check": interpolation,
        "initial_gap_width_m": 2 * BMIN,
        "radial_strip_width_m": BMAX - BMIN,
        "profile_definition": "Lower tenth percentile of outer radial position in fixed downstream bins, at four orbits. Not a physical or optical sharp edge.",
        "sampling": "Nested Halton points, b density proportional to b; identical reflected inner strip. Equal weight 1/N per marker. Finite packet, no reinjection.",
        "files": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(args.out.iterdir())
            if p.is_file()
        },
    }
    (args.out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n"
    )
    print("COMPLETE", json.dumps(summaries), flush=True)


if __name__ == "__main__":
    main()
