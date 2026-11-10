"""Offline check for deeply overlapping grains, including sheared images."""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from .encounter import OMEGA, PERIOD


def deep_overlap_fraction(q, size, radius, time_s):
    xyz = q[:, :3]
    tree = cKDTree(xyz)
    affected = set()
    for ix in (-1, 0, 1):
        dy = (-1.5 * OMEGA * size * ix * time_s + size / 2) % size - size / 2
        for iy in (-1, 0, 1):
            neighbors = tree.query_ball_point(
                xyz + np.array([ix * size, dy + iy * size, 0]), 1.8 * radius
            )
            for i, near in enumerate(neighbors):
                for j in near:
                    if i != j:
                        affected.update((i, j))
    return len(affected) / len(q)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("directory", type=Path)
    args = p.parse_args()
    out = {}
    for path in sorted(args.directory.glob("box_*.npz")):
        meta = json.loads(path.with_suffix(".json").read_text())
        rows = np.loadtxt(path.with_suffix(".csv"), delimiter=",", skiprows=1)
        q = np.load(path)["snapshots"]
        times = rows[[0, len(rows) // 2, -1], 0] * PERIOD
        config = meta["config"]
        out[path.stem] = [
            deep_overlap_fraction(x, config["size_m"], config["radius_m"], t)
            for x, t in zip(q, times)
        ]
    (args.directory / "contact_checks.json").write_text(
        json.dumps(
            {
                "meaning": "Fraction of grains with a center separation below 0.9 times the sum of radii, at start/middle/end. Three snapshots do not bound every contact.",
                "cases": out,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
