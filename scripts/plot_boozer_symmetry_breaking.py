#!/usr/bin/env python3
"""Compare symmetry-breaking |B_mn| profiles from Boozer-v2 results."""

import argparse
from dataclasses import dataclass
from pathlib import Path
import struct

import numpy as np


MAGIC = b"MCBOOZ02"
VERSION = 2
PERIOD = 2.0 * np.pi
FAMILY_NAMES = (
    "rmncc", "rmnss", "zmnsc", "zmncs", "numnsc", "numncs",
)


@dataclass(frozen=True)
class BoozerField:
    """The Boozer-v2 arrays required for a |B| spectrum."""

    source_ns: int
    first_surface: int
    ntheta: int
    nzeta: int
    s: np.ndarray
    mode_m: np.ndarray
    mode_n: np.ndarray
    numnsc: np.ndarray
    numncs: np.ndarray
    b: np.ndarray


@dataclass(frozen=True)
class SymmetryBreakingProfile:
    """Normalized maximum and RMS nonsymmetric Fourier amplitudes."""

    s: np.ndarray
    maximum: np.ndarray
    rms: np.ndarray
    mmax: int
    nmax: int


class _BinaryReader:
    def __init__(self, path):
        self.path = Path(path)
        self.stream = self.path.open("rb")

    def close(self):
        self.stream.close()

    def read(self, count):
        value = self.stream.read(count)
        if len(value) != count:
            raise ValueError(f"truncated Boozer file: {self.path}")
        return value

    def i32(self):
        return struct.unpack("<i", self.read(4))[0]

    def f64(self):
        return struct.unpack("<d", self.read(8))[0]

    def string(self):
        count = self.i32()
        if not 0 <= count < 1 << 24:
            raise ValueError(f"invalid string length in {self.path}")
        return self.read(count).decode("utf-8")

    def array(self, dtype, count):
        size = np.dtype(dtype).itemsize * count
        return np.frombuffer(self.read(size), dtype=dtype).copy()


def load_boozer_binary(path):
    """Load the subset of a cuMES Boozer-v2 binary needed here."""
    reader = _BinaryReader(path)
    try:
        if reader.read(8) != MAGIC or reader.i32() != VERSION:
            raise ValueError(f"unsupported Boozer binary format: {path}")
        reader.string()  # coordinate convention
        reader.string()  # Fourier convention
        reader.string()  # source path
        header = [reader.i32() for _ in range(13)]
        (_, source_ns, _, _, _, _, _, first_surface, ntheta, nzeta,
         mmax, nmax, _) = header
        reader.f64()  # resonance tolerance
        surfaces = source_ns - first_surface
        modes = (mmax + 1) * (nmax + 1)
        if surfaces < 1 or ntheta < 4 or nzeta < 2 or mmax < 0 or nmax < 0:
            raise ValueError(f"invalid Boozer dimensions in {path}")
        s = reader.array("<f8", surfaces)
        reader.array("<f8", surfaces)  # iota
        mode_m = reader.array("<i4", modes)
        mode_n = reader.array("<i4", modes)
        families = {
            name: reader.array("<f8", surfaces * modes).reshape(
                surfaces, modes)
            for name in FAMILY_NAMES
        }
        b = reader.array("<f8", surfaces * nzeta * ntheta).reshape(
            surfaces, nzeta, ntheta)
    finally:
        reader.close()

    if np.any(np.diff(s) <= 0.0) or np.any(b <= 0.0):
        raise ValueError(f"nonphysical Boozer field in {path}")
    if np.any(mode_m < 0) or np.any(mode_n < 0):
        raise ValueError(f"negative stored mode number in {path}")
    return BoozerField(
        source_ns, first_surface, ntheta, nzeta, s, mode_m, mode_n,
        families["numnsc"], families["numncs"], b,
    )


def symmetry_breaking_profile(data, helicity=0):
    """Compute max and RMS |B_mn|/B_00 outside ``n = helicity*m``.

    ``n`` is the field-period toroidal integer stored by Boozer-v2. The file
    uses a grid uniform in ``theta_b`` and source ``zeta``. Since
    ``zeta_b = zeta + nu``, integration includes
    ``d zeta_b / d zeta = 1 + partial_zeta nu``.
    """
    theta = PERIOD * np.arange(data.ntheta) / data.ntheta
    zeta = PERIOD * np.arange(data.nzeta) / data.nzeta
    mt = data.mode_m[:, None] * theta[None, :]
    nz = data.mode_n[:, None] * zeta[None, :]
    cos_mt, sin_mt = np.cos(mt), np.sin(mt)
    cos_nz, sin_nz = np.cos(nz), np.sin(nz)
    nu = (
        np.einsum("sm,mz,mt->szt", data.numnsc, cos_nz, sin_mt,
                  optimize=True)
        + np.einsum("sm,mz,mt->szt", data.numncs, sin_nz, cos_mt,
                    optimize=True)
    )
    dnu_dzeta = (
        np.einsum(
            "sm,mz,mt->szt", -data.mode_n[None, :] * data.numnsc,
            sin_nz, sin_mt, optimize=True)
        + np.einsum(
            "sm,mz,mt->szt", data.mode_n[None, :] * data.numncs,
            cos_nz, cos_mt, optimize=True)
    )
    jacobian = 1.0 + dnu_dzeta
    if np.any(jacobian <= 0.0) or not np.all(np.isfinite(jacobian)):
        raise ValueError("source-zeta to Boozer-zeta map is not monotone")
    zeta_b = zeta[None, :, None] + nu
    b00 = np.mean(data.b * jacobian, axis=(1, 2))
    if np.any(b00 <= 0.0) or not np.all(np.isfinite(b00)):
        raise ValueError("invalid B_00 normalization")

    mmax = int(np.max(data.mode_m))
    nmax = int(np.max(data.mode_n))
    amplitudes = []
    for m in range(mmax + 1):
        for n in range(-nmax, nmax + 1):
            # cos(0*theta-n*zeta) duplicates the -n mode, so retain only
            # positive n when m=0, matching the standard Boozer mode set.
            if m == 0 and n < 0:
                continue
            if n == helicity * m:
                continue
            phase = m * theta[None, None, :] - n * zeta_b
            cosine = 2.0 * np.mean(
                data.b * np.cos(phase) * jacobian, axis=(1, 2))
            sine = 2.0 * np.mean(
                data.b * np.sin(phase) * jacobian, axis=(1, 2))
            amplitudes.append(np.hypot(cosine, sine) / b00)
    if not amplitudes:
        raise ValueError("the selected truncation has no symmetry-breaking modes")
    amplitudes = np.asarray(amplitudes).T
    return SymmetryBreakingProfile(
        data.s.copy(), np.max(amplitudes, axis=1),
        np.sqrt(np.sum(amplitudes * amplitudes, axis=1)), mmax, nmax,
    )


def plot_comparison(profiles, labels, output, helicity=0):
    """Write a two-panel initial/final symmetry-breaking comparison."""
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharex=True)
    for profile, label in zip(profiles, labels):
        axes[0].semilogy(profile.s, profile.maximum, linewidth=2, label=label)
        axes[1].semilogy(profile.s, profile.rms, linewidth=2, label=label)
    axes[0].set_title(r"Maximum symmetry-breaking $B_{m,n}$")
    axes[1].set_title(r"RMS symmetry-breaking $B_{m,n}$")
    axes[0].set_ylabel(r"$max |B_{m,n}|/B_{0,0}$")
    axes[1].set_ylabel(r"$\sqrt{\sum |B_{m,n}|^2}/B_{0,0}$")
    for axis in axes:
        axis.set_xlabel(r"Normalized toroidal flux $s$")
        axis.set_xlim(0.0, 1.0)
        axis.grid(True, which="both", alpha=0.25)
        axis.legend()
    mmax = min(profile.mmax for profile in profiles)
    nmax = min(profile.nmax for profile in profiles)
    symmetric = "n=0" if helicity == 0 else f"n={helicity}m"
    figure.suptitle(
        "Boozer |B| symmetry breaking "
        f"(symmetric family: {symmetric}; 0≤m≤{mmax}, |n|≤{nmax})")
    figure.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=220)
    plt.close(figure)


def _parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("initial", type=Path, help="initial Boozer-v2 .bin")
    parser.add_argument("final", type=Path, help="final Boozer-v2 .bin")
    parser.add_argument("--output", "-o", type=Path, required=True)
    parser.add_argument(
        "--helicity", type=int, default=0,
        help="symmetric field-period family n=helicity*m (QA default: 0)")
    parser.add_argument(
        "--labels", nargs=2, default=("Initial", "Final"),
        metavar=("INITIAL", "FINAL"))
    return parser.parse_args()


def main():
    arguments = _parse_arguments()
    fields = [
        load_boozer_binary(arguments.initial),
        load_boozer_binary(arguments.final),
    ]
    profiles = [
        symmetry_breaking_profile(field, arguments.helicity)
        for field in fields
    ]
    plot_comparison(
        profiles, arguments.labels, arguments.output, arguments.helicity)
    for label, profile in zip(arguments.labels, profiles):
        print(
            f"{label}: edge max={profile.maximum[-1]:.12e}, "
            f"edge RMS={profile.rms[-1]:.12e}")
    print(f"saved {arguments.output}")


if __name__ == "__main__":
    main()
