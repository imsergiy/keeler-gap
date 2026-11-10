import numpy as np
import pytest

from keeler.encounter import Encounter, integrate
from keeler.waves import (
    BMAX,
    BMIN,
    MOVIE_ORBITS,
    PROFILE_BINS_M,
    interpolate_positions,
    packet_sample,
    profile,
    profile_difference,
    new_tracks,
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


def test_profile_comparison_requires_every_measurement_bin():
    a = np.array([2.0, 0.0])
    b = np.array([0.0, 0.0])
    assert profile_difference(a, b) == np.sqrt(2)
    assert profile_difference(np.array([0.0, np.nan]), b) is None


@pytest.mark.parametrize(
    "config", [Encounter(radial_epicycle_m=4500), Encounter(vertical_amplitude_m=8600)]
)
def test_packet_reuse_rejects_a_time_dependent_moon(config):
    with pytest.raises(ValueError, match="stationary"):
        new_tracks(np.array([BMIN]), np.array([-10.0]), config)


@pytest.mark.parametrize("mass_factor", [0.0, 0.5])
def test_track_invariant_uses_the_configured_moon_mass(mass_factor):
    _, _, drift = new_tracks(
        np.array([BMIN]),
        np.array([-10.0]),
        Encounter(mass_factor=mass_factor, samples=201),
    )
    assert drift < 1e-7
