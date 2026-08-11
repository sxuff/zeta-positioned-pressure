#!/usr/bin/env python3
"""Build interval tables for the positioned-pressure certificate.

Trust base: Python, mpmath.iv interval arithmetic, IEEE-754 binary64, and the
short formulas below.  The repository-native Arb verifier should be used as an
independent reproduction before treating the result as publication-ready.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import struct
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import numpy as np
from mpmath import iv

WINDOW_DEN = 1_000_000_000
WINDOW_NUM = (
    1_000_000_000,
    3_322_500,
    -7_609_135,
    1_190_194,
    -731_476,
    -1_680_572,
    1_141_360,
)
TARGET = Fraction(51_063, 10_000_000)
PRESSURE_DEN = 2_300_000_000
PRESSURE_NUM = (831_522, 1_096_590, 1_071_888, 1_071_888, 1_096_590, 831_522)
SERIES_TERMS = 24
SERIES_RADIUS = 0.75
SERIES_TAIL = 1e-60


def _configure(precision: int) -> None:
    iv.dps = precision


def _sinc_derivatives(z):
    """Interval enclosures for sinc(z), sinc'(z), sinc''(z)."""
    upper = max(abs(float(z.a)), abs(float(z.b)))
    if upper < SERIES_RADIUS:
        z2 = z * z
        value, first, second = iv.mpf(1), iv.mpf(0), iv.mpf(0)
        power, sign = iv.mpf(1), -1
        for n in range(1, SERIES_TERMS):
            fact = math.factorial(2 * n + 1)
            term0 = power * z2 / fact
            term1 = (2 * n) * power * z / fact
            term2 = (2 * n) * (2 * n - 1) * power / fact
            if sign > 0:
                value += term0
                first += term1
                second += term2
            else:
                value -= term0
                first -= term1
                second -= term2
            sign = -sign
            power *= z2
        # At |z| <= 3/4 the largest first omitted second-derivative term is
        # below 7e-66; 1e-60 is therefore a conservative common tail radius.
        tail = iv.mpf([-SERIES_TAIL, SERIES_TAIL])
        return value + tail, first + tail, second + tail

    sine, cosine = iv.sin(z), iv.cos(z)
    z2 = z * z
    return (
        sine / z,
        (z * cosine - sine) / z2,
        ((2 - z2) * sine - 2 * z * cosine) / (z2 * z),
    )


def _constants():
    coeff = [iv.mpf(n) / WINDOW_DEN for n in WINDOW_NUM]
    omega = [iv.sqrt(2)] + [2 * j * iv.pi for j in range(1, 7)]
    return coeff, omega


def _kernel_derivatives(x, coeff, omega):
    value, first, second = iv.mpf(0), iv.mpf(0), iv.mpf(0)
    for c, w in zip(coeff, omega):
        left = w / 2 - iv.pi * x
        right = w / 2 + iv.pi * x
        lv, l1, l2 = _sinc_derivatives(left)
        rv, r1, r2 = _sinc_derivatives(right)
        value += c * (lv + rv) / 2
        first += c * iv.pi * (r1 - l1) / 2
        second += c * iv.pi**2 * (l2 + r2) / 2
    return value, first, second


def _coarse_worker(task):
    start, stop, grid, precision = task
    _configure(precision)
    coeff, omega = _constants()
    k0 = _kernel_derivatives(iv.mpf(0), coeff, omega)[0]
    k0sq = k0 * k0
    lower = np.empty(stop - start, dtype=np.float64)
    second_lower = np.empty(stop - start, dtype=np.float64)
    for out, index in enumerate(range(start, stop)):
        x = iv.mpf([index / grid, (index + 1) / grid])
        k, k1, k2 = _kernel_derivatives(x, coeff, omega)
        w = k * k / k0sq
        w2 = 2 * (k1 * k1 + k * k2) / k0sq
        lo = max(0.0, float(w.a))
        lower[out] = math.nextafter(lo, -math.inf) if lo > 0 else 0.0
        second_lower[out] = math.nextafter(float(w2.a), -math.inf)
    return start, lower, second_lower


def _midpoint_worker(task):
    start, stop, grid, precision = task
    _configure(precision)
    coeff, omega = _constants()
    k0 = _kernel_derivatives(iv.mpf(0), coeff, omega)[0]
    k0sq = k0 * k0
    result = np.empty((4, stop - start), dtype=np.float64)
    for out, index in enumerate(range(start, stop)):
        x = iv.mpf(index) / (2 * grid)
        k, k1, _ = _kernel_derivatives(x, coeff, omega)
        w = k * k / k0sq
        derivative = 2 * k * k1 / k0sq
        result[0, out] = math.nextafter(float(w.a), -math.inf)
        result[1, out] = math.nextafter(float(w.b), math.inf)
        result[2, out] = math.nextafter(float(derivative.a), -math.inf)
        result[3, out] = math.nextafter(float(derivative.b), math.inf)
    return start, result


def _chunks(length: int, workers: int, grid: int, precision: int):
    chunk = math.ceil(length / workers)
    return [(start, min(length, start + chunk), grid, precision)
            for start in range(0, length, chunk)]


def _sha256(values: np.ndarray) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(struct.pack(">d", float(value)))
    return digest.hexdigest()


def _write(path: Path, values: np.ndarray) -> None:
    values.astype("<f8").tofile(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("tables"))
    parser.add_argument("--grid", type=int, default=4_000)
    parser.add_argument("--precision", type=int, default=45)
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--extra-cells", type=int, default=33)
    parser.add_argument("--smoke-cells", type=int, default=0,
                        help="build only this many coarse cells for a code-path smoke test")
    args = parser.parse_args()

    min_pressure = min(Fraction(n, PRESSURE_DEN) for n in PRESSURE_NUM)
    exact_cutoff = TARGET * args.grid / min_pressure
    cell_count = (exact_cutoff.numerator + exact_cutoff.denominator - 1) // exact_cutoff.denominator
    cell_count += args.extra_cells
    if args.smoke_cells:
        cell_count = args.smoke_cells
    midpoint_count = 2 * cell_count + 1
    workers = max(1, min(args.workers, cell_count))
    args.output.mkdir(parents=True, exist_ok=True)

    context = mp.get_context("spawn")
    with context.Pool(workers) as pool:
        coarse_parts = pool.map(
            _coarse_worker,
            _chunks(cell_count, workers, args.grid, args.precision),
        )
    lower = np.empty(cell_count)
    second = np.empty(cell_count)
    for start, values, values2 in coarse_parts:
        lower[start:start + len(values)] = values
        second[start:start + len(values2)] = values2

    midpoint_workers = max(1, min(args.workers, midpoint_count))
    with context.Pool(midpoint_workers) as pool:
        midpoint_parts = pool.map(
            _midpoint_worker,
            _chunks(midpoint_count, midpoint_workers, args.grid, args.precision),
        )
    mid = np.empty((4, midpoint_count))
    for start, values in midpoint_parts:
        mid[:, start:start + values.shape[1]] = values

    files = {
        "w_lower.bin": lower,
        "w_second_lower.bin": second,
        "w_mid_lower.bin": mid[0],
        "w_mid_upper.bin": mid[1],
        "w_prime_mid_lower.bin": mid[2],
        "w_prime_mid_upper.bin": mid[3],
    }
    manifest = {
        "grid": args.grid,
        "precision_decimal_digits": args.precision,
        "target": f"{TARGET.numerator}/{TARGET.denominator}",
        "coarse_cells": cell_count,
        "midpoint_values": midpoint_count,
        "files": {},
    }
    for name, values in files.items():
        _write(args.output / name, values)
        manifest["files"][name] = {
            "length": len(values),
            "sha256_big_endian_float_stream": _sha256(values),
        }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
