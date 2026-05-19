import math
import unittest

from snom_web.calculator import calculate_spectrum, parse_numeric_table


class ParseNumericTableTests(unittest.TestCase):
    def test_accepts_mixed_delimiters(self) -> None:
        table = parse_numeric_table(
            "1650,2.2,0.1\n1651 2.3 0.2\n1652\t2.4\t0.3\n",
            minimum_columns=3,
            label="Permittivity file",
        )
        self.assertEqual(table.shape, (3, 3))
        self.assertAlmostEqual(table[1, 1], 2.3)


class CalculatorSmokeTests(unittest.TestCase):
    def test_fast_bulk_calculation_returns_expected_columns(self) -> None:
        payload = {
            "modelType": "bulk",
            "algorithm": "fast",
            "tip": {
                "radius": 25,
                "length": 600,
                "amplitude": 40,
                "sampleGap": 2,
                "referenceGap": 2,
                "epsilonReal": -1000,
                "epsilonImag": 800,
            },
            "farField": {
                "enabled": False,
                "coefficient": 1,
                "angleDegrees": 45,
            },
            "reference": {
                "epsilonReal": -5000,
                "epsilonImag": 1000,
            },
            "harmonicNumber": 3,
            "frequencyRange": {
                "start": 1650,
                "end": 1660,
                "points": 4,
            },
            "sample": {
                "inputMethod": "builtIn",
                "builtInMaterial": "pmma",
                "uploadedPermittivity": None,
            },
            "experimentalSpectrum": {
                "enabled": False,
                "content": None,
            },
            "discretization": {
                "N": 50,
                "M": 20,
            },
        }

        result = calculate_spectrum(payload)
        self.assertEqual(len(result["frequency"]), 4)
        self.assertEqual(len(result["sigmaReal"]), 4)
        self.assertTrue(all(math.isfinite(value) for value in result["amplitude"]))
        self.assertTrue(all(math.isfinite(value) for value in result["phase"]))


if __name__ == "__main__":
    unittest.main()
