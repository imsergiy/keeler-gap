import numpy as np

from keeler.encounter import Encounter, integrate
from keeler.waves import (
    BMAX,
    BMIN,
    MOVIE_ORBITS,
    PROFILE_BINS_M,
    interpolate_positions,
    packet_sample,
    profile,
)


def test_packet_samples_are_nested_and_cover_the_same_physical_domain():
    b, age = packet_sample(125)
    fine_b, fine_age = packet_sample(500)
    np.testing.assert_array_equal(b, fine_b[:125])
    np.testing.assert_array_equal(age, fine_age[:125])
    assert np.all((b >= BMIN) & (b < BMAX))
    assert np.all((age >= -8 * np.pi) & (age < 0))


def test_incoming_sampling_has_uniform_area_density():
    b, _ = packet_sample(2000)
    cumulative = np.sort((b * b - BMIN**2) / (BMAX**2 - BMIN**2))
    expected = np.arange(len(b)) / len(b)
    assert np.max(abs(cumulative - expected)) < 0.002


def test_unforced_interpolation_preserves_exact_shear():
    b = np.array([21_000.0, 25_000.0, 31_000.0])
    phase = np.array([-20.0, -10.0, -3.0])
    t, q = integrate(b, Encounter(mass_factor=0, samples=201))
    ages = MOVIE_ORBITS[:, None] * 2 * np.pi + phase
    positions = interpolate_positions(t, q, ages)
    np.testing.assert_allclose(
        positions[:, :, 0], np.broadcast_to(b, ages.shape), atol=1e-7
    )
    np.testing.assert_allclose(positions[:, :, 1], -1.5 * b * ages, atol=1e-6)


def test_profile_reports_empty_bins_without_inventing_an_edge():
    center = (PROFILE_BINS_M[0] + PROFILE_BINS_M[1]) / 2
    p = np.array([[21000.0, -center], [22000.0, -center], [23000.0, -center]])
    counts, edge = profile(p)
    assert counts[0] == 3
    assert edge[0] == 21200
    assert np.all(counts[1:] == 0)
    assert np.isnan(edge[1:]).all()


def test_profile_window_stays_behind_the_slowest_packet_front():
    b = np.array([BMIN])
    t, q = integrate(b, Encounter(samples=2401))
    front = interpolate_positions(t, q, np.array([[8 * np.pi]]))[0, 0]
    assert PROFILE_BINS_M[-1] < -front[1] - 40_000
