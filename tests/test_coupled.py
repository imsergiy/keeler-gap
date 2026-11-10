from dataclasses import replace
import numpy as np
import pytest
import rebound
from scipy.integrate import solve_ivp
from scipy.spatial.distance import pdist

from keeler.coupled import (
    Packet,
    initial_packet,
    setup,
    read_frame,
    deep_overlap_fraction,
    install_collision_law,
)
from keeler.encounter import Encounter, HILL, OMEGA, rhs


def test_larger_packets_retain_identical_central_grains_and_density():
    a, b = Packet(count=128), Packet(count=512)
    x, y = initial_packet(a), initial_packet(b)
    lookup = {identity: i for i, identity in enumerate(y["ids"])}
    indices = [lookup[identity] for identity in x["ids"]]
    np.testing.assert_array_equal(x["state"], y["state"][indices])
    np.testing.assert_array_equal(x["ids"], y["ids"][y["core"]])
    assert a.count * a.grain_mass / (2 * a.side_length_m**2) == pytest.approx(
        a.surface_density
    )
    assert b.count * b.grain_mass / (2 * b.side_length_m**2) == pytest.approx(
        a.surface_density
    )
    assert pdist(x["state"][:, :3]).min() > 2 * a.radius_m
    q = x["state"]
    assert np.max(abs(q[:, 3])) < 0.01
    assert np.max(abs(q[:, 4] + 1.5 * OMEGA * q[:, 0])) < 0.01


def test_gravity_tree_is_populated_and_agrees_with_direct_forces():
    cfg = Packet()
    tree, layout, _ = setup(cfg)
    tree.gravity = "tree"
    tree.opening_angle2 = 0.1**2
    direct, _, _ = setup(cfg)
    direct.gravity = "basic"
    tree.steps(1)
    direct.steps(1)

    def acceleration(sim):
        return {p.hash.value: np.array([p.ax, p.ay, p.az]) for p in sim.particles}

    a, b = acceleration(tree), acceleration(direct)
    actual = np.array([a[int(i)] for i in layout["ids"]])
    ref = np.array([b[int(i)] for i in layout["ids"]])
    assert np.linalg.norm(ref) > 1e-8
    assert np.linalg.norm(actual - ref) / np.linalg.norm(ref) < 0.005


def test_tree_contacts_match_direct_search_in_an_advecting_packet():
    results = []
    for method in ["tree", "direct"]:
        sim, layout, contacts = setup(Packet())
        sim.collision = method
        sim.steps(1500)
        q, _ = read_frame(sim, Packet(), layout)
        assert contacts["approaching_calls"] > 0
        assert deep_overlap_fraction(q, 2) == 0
        results.append((q, contacts))
    np.testing.assert_allclose(results[0][0], results[1][0], rtol=1e-9, atol=1e-6)
    assert results[0][1]["dissipated_J"] == pytest.approx(
        results[1][1]["dissipated_J"], rel=1e-9
    )


def test_point_moon_force_matches_existing_hill_equations():
    cfg = Packet(count=8, self_gravity=False, collisions=False)
    sim, layout, _ = setup(cfg)
    units = np.array([HILL, HILL, HILL, OMEGA * HILL, OMEGA * HILL, OMEGA * HILL])
    initial = layout["state"][0] / units
    duration = 1000.0
    sol = solve_ivp(
        rhs,
        [0, duration * OMEGA],
        initial,
        args=(Encounter(),),
        method="DOP853",
        rtol=1e-12,
        atol=1e-14,
    )
    sim.steps(int(duration / cfg.dt_s))
    q, _ = read_frame(sim, cfg, layout)
    np.testing.assert_allclose(q[0, :3], sol.y[:3, -1] * HILL, rtol=0, atol=0.001)


def test_contact_counter_records_actual_energy_loss():
    cfg = Packet(count=8, moon_gravity=False, self_gravity=False, restitution=0.5)
    sim = rebound.Simulation()
    sim.gravity = "none"
    sim.integrator = "leapfrog"
    sim.collision = "direct"
    sim.collision_resolve = "hardsphere"
    sim.dt = 0.001
    contacts = install_collision_law(sim, cfg)
    sim.add(m=cfg.grain_mass, r=2, x=-2.1, vx=1)
    sim.add(m=cfg.grain_mass, r=2, x=2.1, vx=-1)
    sim.steps(150)
    assert contacts["approaching_calls"] == 1
    assert contacts["dissipated_J"] == pytest.approx(0.75 * cfg.grain_mass)
    assert sum(p.m * p.vx for p in sim.particles) == pytest.approx(0, abs=1e-10)


def test_inclined_moon_drives_the_same_three_dimensional_equations():
    cfg = Packet(count=8, collisions=False, moon_vertical_amplitude_m=8600)
    sim, layout, _ = setup(cfg)
    while sim.N > 2:
        sim.remove(index=sim.N - 1)
    sim.particles[1].y = 0
    p = sim.particles[1]
    initial = np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz])
    units = np.array([HILL, HILL, HILL, OMEGA * HILL, OMEGA * HILL, OMEGA * HILL])
    duration = 1000.0
    reference = solve_ivp(
        rhs,
        [0, duration * OMEGA],
        initial / units,
        args=(Encounter(vertical_amplitude_m=8600),),
        method="DOP853",
        rtol=1e-12,
        atol=1e-14,
    )
    sim.steps(int(duration / cfg.dt_s))
    p = sim.particles[1]
    np.testing.assert_allclose(
        [p.x, p.y, p.z], reference.y[:3, -1] * HILL, rtol=0, atol=0.001
    )
    moon = sim.particles[0]
    assert moon.z == pytest.approx(8600 * np.sin(duration * OMEGA), abs=1e-5)


def test_direct_packet_energy_is_conservative_without_contacts():
    from keeler.coupled import energy

    cfg = Packet(count=8, collisions=False)
    sim, layout, _ = setup(cfg)
    q, moon = read_frame(sim, cfg, layout)
    before = energy(sim, cfg, q, moon)
    sim.steps(100)
    q, moon = read_frame(sim, cfg, layout)
    after = energy(sim, cfg, q, moon)
    assert abs(after - before) < 1e-4
