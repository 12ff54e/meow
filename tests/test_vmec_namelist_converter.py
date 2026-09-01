#!/usr/bin/env python3

from pathlib import Path
import json
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from vmec_namelist_to_cumes_json import convert_text, read_wout_axis


class NamelistConverterTest(unittest.TestCase):
    def test_fixed_boundary_subset_and_axis_padding(self):
        converted = convert_text("""
&INDATA
  NS_ARRAY = 5 11
  NITER_ARRAY = 100 200
  FTOL_ARRAY = 1d-10 1d-12
  NFP = 2
  MPOL = 4
  NTOR = 2
  NCURR = 1
  CURTOR = 0
  LASYM = F
  LFREEB = F
  AM = 0
  RAXIS_CC = 1
  ZAXIS_CS = 0
  RBC(0,0) = 1.0, ZBS(0,0) = 0.0
  RBC(-1,1) = 0.2, ZBS(-1,1) = -0.3
/
""")
        self.assertEqual(converted["ns_array"], [5, 11])
        self.assertEqual(converted["ftol_array"], [1e-10, 1e-12])
        self.assertEqual(converted["raxis_c"], [1.0, 0.0, 0.0])
        self.assertEqual(converted["rbc"][1],
                         {"n": -1, "m": 1, "value": 0.2})
        self.assertEqual(converted["zbs"][1]["value"], -0.3)

        clamped = convert_text("""
&INDATA
  MPOL=2
  NTOR=0
  NITER_ARRAY=100 200
  FTOL_ARRAY=1e-17 2e-16
  RBC(0,0)=1
  ZBS(0,0)=0
/
""", minimum_ftol=1e-16, minimum_niter=150)
        self.assertEqual(clamped["ftol_array"], [1e-16, 2e-16])
        self.assertEqual(clamped["niter_array"], [150, 200])

    def test_unknown_active_key_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "unsupported active"):
            convert_text("""
&INDATA
  MPOL=2
  NTOR=0
  MYSTERY=1
  RBC(0,0)=1
  ZBS(0,0)=0
/
""")

    def test_lower_resolution_wout_axis_is_zero_padded(self):
        try:
            from scipy.io import netcdf_file
            import numpy as np
        except ImportError:
            self.skipTest("scipy and numpy are required for wout-axis test")

        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "wout.nc"
            with netcdf_file(path, "w") as wout:
                wout.createDimension("axis_mode", 2)
                raxis = wout.createVariable(
                    "raxis_cc", "d", ("axis_mode",))
                zaxis = wout.createVariable(
                    "zaxis_cs", "d", ("axis_mode",))
                raxis[:] = np.array([1.0, 0.2])
                zaxis[:] = np.array([0.0, -0.2])

            raxis, zaxis = read_wout_axis(path, ntor=3)
            self.assertEqual(raxis, [1.0, 0.2, 0.0, 0.0])
            self.assertEqual(zaxis, [0.0, -0.2, 0.0, 0.0])

    def test_checked_in_landreman_analytic_boundaries_are_sparse(self):
        example_directory = (Path(__file__).resolve().parents[1]
                             / "examples" / "landreman")
        qa = json.loads((example_directory / "qa_analytic.json").read_text())
        qh = json.loads((example_directory / "qh_analytic.json").read_text())

        self.assertEqual(qa["nfp"], 2)
        self.assertEqual(qa["rbc"], [
            {"m": 0, "n": 0, "value": 1.0},
            {"m": 1, "n": 0, "value": 0.2},
        ])
        self.assertEqual(qa["zbs"], [
            {"m": 0, "n": 0, "value": 0.0},
            {"m": 1, "n": 0, "value": 0.2},
        ])

        self.assertEqual(qh["nfp"], 4)
        self.assertEqual(qh["rbc"], [
            {"m": 0, "n": 0, "value": 1.0},
            {"m": 1, "n": 0, "value": 0.2},
            {"m": 0, "n": 1, "value": 0.2},
        ])
        self.assertEqual(qh["zbs"], [
            {"m": 0, "n": 0, "value": 0.0},
            {"m": 1, "n": 0, "value": 0.2},
            {"m": 0, "n": 1, "value": 0.2},
        ])


if __name__ == "__main__":
    unittest.main()
