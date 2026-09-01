#!/usr/bin/env python3
"""Unit checks for the Boozer symmetry-breaking diagnostic."""

import importlib.util
from pathlib import Path
import sys

import numpy as np


SCRIPT = Path(__file__).parents[1] / "scripts" / \
    "plot_boozer_symmetry_breaking.py"
SPEC = importlib.util.spec_from_file_location("boozer_symmetry", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def main():
    surfaces = 3
    ntheta = 24
    nzeta = 20
    theta = MODULE.PERIOD * np.arange(ntheta) / ntheta
    zeta = MODULE.PERIOD * np.arange(nzeta) / nzeta
    s = np.array([0.2, 0.6, 1.0])
    first = np.array([0.02, 0.04, 0.08])
    second = np.array([0.03, 0.01, 0.06])
    b00 = np.array([1.0, 1.2, 1.4])
    phase_first = theta[None, :] - zeta[:, None]
    phase_second = 2.0 * theta[None, :] + zeta[:, None]
    b = (
        b00[:, None, None]
        + first[:, None, None] * np.cos(phase_first)[None, :, :]
        + second[:, None, None] * np.cos(phase_second)[None, :, :]
        + 0.1 * np.cos(theta)[None, None, :]
    )
    mode_m = np.array([0, 1, 2, 0, 1, 2], dtype=np.int32)
    mode_n = np.array([0, 0, 0, 1, 1, 1], dtype=np.int32)
    zeros = np.zeros((surfaces, mode_m.size))
    data = MODULE.BoozerField(
        source_ns=surfaces + 1,
        first_surface=1,
        ntheta=ntheta,
        nzeta=nzeta,
        s=s,
        mode_m=mode_m,
        mode_n=mode_n,
        numnsc=zeros,
        numncs=zeros,
        b=b,
    )
    profile = MODULE.symmetry_breaking_profile(data, helicity=0)
    expected_maximum = np.maximum(first, second) / b00
    expected_rms = np.sqrt(first * first + second * second) / b00
    np.testing.assert_allclose(profile.maximum, expected_maximum, atol=1e-14)
    np.testing.assert_allclose(profile.rms, expected_rms, atol=1e-14)
    print("Boozer symmetry-breaking checks passed")


if __name__ == "__main__":
    main()
