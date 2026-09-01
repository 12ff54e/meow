#!/usr/bin/env python3
"""Convert the supported fixed-boundary VMEC namelist subset to cuMES JSON."""

import argparse
import json
from pathlib import Path
import re


HARMONIC = re.compile(
    r"\b(RBC|ZBS)\s*\(\s*([+-]?\d+)\s*,\s*(\d+)\s*\)\s*=\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?)",
    re.IGNORECASE,
)

INTEGER_KEYS = {"mpol", "ntor", "nfp", "ncurr", "nvacskip"}
FLOAT_KEYS = {
    "delt", "phiedge", "pres_scale", "gamma", "spres_ped", "bloat",
    "curtor", "tcon0",
}
BOOLEAN_KEYS = {"lasym", "lfreeb"}
ARRAY_KEYS = {
    "ns_array", "niter_array", "ftol_array", "am", "ac", "ai", "aphi",
}
STRING_KEYS = {"pmass_type", "pcurr_type", "piota_type"}
AXIS_KEYS = {"raxis_cc": "raxis_c", "zaxis_cs": "zaxis_s"}
IGNORED_KEYS = {"niter", "nstep", "precon_type", "prec2d_threshold"}
SUPPORTED_KEYS = (
    INTEGER_KEYS | FLOAT_KEYS | BOOLEAN_KEYS | ARRAY_KEYS | STRING_KEYS
    | set(AXIS_KEYS) | IGNORED_KEYS
)


def _tokens(text):
    return [token for token in re.split(r"[\s,]+", text.strip()) if token]


def _number(token):
    return float(token.replace("D", "e").replace("d", "e"))


def _logical(token):
    normalized = token.strip().strip(".").lower()
    if normalized in ("t", "true"):
        return True
    if normalized in ("f", "false"):
        return False
    raise ValueError(f"invalid Fortran logical value: {token}")


def _string(token):
    token = token.strip()
    if len(token) >= 2 and token[0] == token[-1] and token[0] in "'\"":
        return token[1:-1].lower()
    return token.lower()


def _read_assignments(text):
    harmonics = {"rbc": [], "zbs": []}
    assignments = {}
    current_key = None
    current_value = []

    def finish_assignment():
        nonlocal current_key, current_value
        if current_key is not None:
            assignments[current_key] = " ".join(current_value).strip(" ,")
        current_key = None
        current_value = []

    for raw_line in text.splitlines():
        line = raw_line.split("!", 1)[0]
        for match in HARMONIC.finditer(line):
            family, n, m, value = match.groups()
            harmonics[family.lower()].append({
                "n": int(n), "m": int(m), "value": _number(value),
            })
        line = HARMONIC.sub("", line).strip(" ,\t")
        if not line or line.startswith("&") or line == "/":
            continue
        match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)$", line)
        if match:
            finish_assignment()
            current_key = match.group(1).lower()
            current_value = [match.group(2)]
        elif current_key is not None:
            current_value.append(line)
        else:
            raise ValueError(f"cannot parse namelist line: {raw_line}")
    finish_assignment()
    return assignments, harmonics


def convert_text(text, minimum_ftol=None, minimum_niter=None):
    assignments, harmonics = _read_assignments(text)
    unknown = sorted(set(assignments) - SUPPORTED_KEYS)
    if unknown:
        raise ValueError("unsupported active namelist keys: " + ", ".join(unknown))

    output = {}
    for key in INTEGER_KEYS:
        if key in assignments:
            values = _tokens(assignments[key])
            if len(values) != 1:
                raise ValueError(f"{key} must be a scalar")
            output[key] = int(_number(values[0]))
    for key in FLOAT_KEYS:
        if key in assignments:
            values = _tokens(assignments[key])
            if len(values) != 1:
                raise ValueError(f"{key} must be a scalar")
            output[key] = _number(values[0])
    for key in BOOLEAN_KEYS:
        if key in assignments:
            output[key] = _logical(assignments[key])
    for key in ARRAY_KEYS:
        if key in assignments:
            values = [_number(token) for token in _tokens(assignments[key])]
            if key in ("ns_array", "niter_array"):
                values = [int(value) for value in values]
            output[key] = values
    if minimum_ftol is not None and "ftol_array" in output:
        if minimum_ftol <= 0.0:
            raise ValueError("minimum_ftol must be positive")
        output["ftol_array"] = [
            max(value, minimum_ftol) for value in output["ftol_array"]
        ]
    if minimum_niter is not None and "niter_array" in output:
        if minimum_niter <= 0:
            raise ValueError("minimum_niter must be positive")
        output["niter_array"] = [
            max(value, minimum_niter) for value in output["niter_array"]
        ]
    for key in STRING_KEYS:
        if key in assignments:
            output[key] = _string(assignments[key])

    if "ntor" not in output:
        raise ValueError("ntor is required before axis arrays can be expanded")
    axis_size = output["ntor"] + 1
    for source, destination in AXIS_KEYS.items():
        if source in assignments:
            values = [_number(token) for token in _tokens(assignments[source])]
            if len(values) > axis_size:
                raise ValueError(f"{source} has more than ntor+1 entries")
            output[destination] = values + [0.0] * (axis_size - len(values))

    output.update(harmonics)
    if not output["rbc"] or not output["zbs"]:
        raise ValueError("both RBC and ZBS boundary families are required")
    if output.get("lasym", False):
        raise ValueError("cuMES does not support LASYM=T")
    if output.get("lfreeb", False):
        raise ValueError("this converter is limited to fixed-boundary inputs")
    return output


def convert_file(input_path, output_path, minimum_ftol=None,
                 minimum_niter=None):
    converted = convert_text(
        Path(input_path).read_text(encoding="utf-8"), minimum_ftol,
        minimum_niter)
    Path(output_path).write_text(
        json.dumps(converted, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")


def read_wout_axis(wout_path, ntor):
    try:
        from scipy.io import netcdf_file
    except ImportError as error:
        raise RuntimeError("--wout-axis requires scipy") from error
    with netcdf_file(str(wout_path), "r", mmap=False) as wout:
        required = ("raxis_cc", "zaxis_cs")
        missing = [name for name in required if name not in wout.variables]
        if missing:
            raise ValueError("wout is missing axis variables: "
                             + ", ".join(missing))
        raxis = [float(value) for value in wout.variables["raxis_cc"].data]
        zaxis = [float(value) for value in wout.variables["zaxis_cs"].data]
    expected = ntor + 1
    if not raxis or not zaxis:
        raise ValueError("wout axis coefficient arrays must not be empty")
    # An iteration-zero wout may have been written after a continuation driver
    # lowered ntor below the value in the source namelist. Modes absent from
    # that lower-resolution axis predictor are identically zero.
    raxis = raxis[:expected] + [0.0] * max(0, expected - len(raxis))
    zaxis = zaxis[:expected] + [0.0] * max(0, expected - len(zaxis))
    return raxis, zaxis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--minimum-ftol", type=float,
        help="explicitly clamp VMEC stage tolerances to this cuMES floor")
    parser.add_argument(
        "--minimum-niter", type=int,
        help="explicitly raise VMEC stage iteration caps for cuMES")
    parser.add_argument(
        "--wout-axis", type=Path,
        help="replace VMEC's zero auto-axis input with an axis from this wout")
    arguments = parser.parse_args()
    converted = convert_text(
        arguments.input.read_text(encoding="utf-8"), arguments.minimum_ftol,
        arguments.minimum_niter)
    if arguments.wout_axis is not None:
        converted["raxis_c"], converted["zaxis_s"] = read_wout_axis(
            arguments.wout_axis, converted["ntor"])
    arguments.output.write_text(
        json.dumps(converted, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
