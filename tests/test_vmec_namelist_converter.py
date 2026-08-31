#!/usr/bin/env python3

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from vmec_namelist_to_cumes_json import convert_text


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


if __name__ == "__main__":
    unittest.main()
