import math
import unittest
from unittest import mock

import numpy as np

from snom_web.calculator import ValidationError, calculate_spectrum, parse_numeric_table
from snom_web.materials import generate_epsilon_drude_lorentz
from snom_web.server import _run_calculation_with_timeout


def build_payload() -> dict[str, object]:
    return {
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
            "drudeLorentz": None,
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


def build_layered_payload(structure: str = "single") -> dict[str, object]:
    payload = build_payload()
    payload["modelType"] = "layered"
    payload["sample"] = None
    payload["frequencyRange"] = {
        "start": 900,
        "end": 920,
        "points": 3,
    }
    payload["layered"] = {
        "structure": structure,
        "imageCharges": 2,
        "layer1": {
            "thickness": 20,
            "material": {
                "inputMethod": "builtIn",
                "builtInMaterial": "pmma",
                "uploadedPermittivity": None,
                "drudeLorentz": None,
            },
        },
        "layer2": None,
        "substrate": {
            "inputMethod": "builtIn",
            "builtInMaterial": "si",
            "uploadedPermittivity": None,
            "drudeLorentz": None,
        },
    }
    if structure == "double":
        payload["layered"]["layer2"] = {
            "thickness": 50,
            "material": {
                "inputMethod": "builtIn",
                "builtInMaterial": "sic_3c",
                "uploadedPermittivity": None,
                "drudeLorentz": None,
            },
        }
    return payload


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
        payload = build_payload()

        result = calculate_spectrum(payload)
        self.assertEqual(len(result["frequency"]), 4)
        self.assertEqual(len(result["sigmaReal"]), 4)
        self.assertTrue(all(math.isfinite(value) for value in result["amplitude"]))
        self.assertTrue(all(math.isfinite(value) for value in result["phase"]))
        self.assertEqual(result["metadata"]["sampleDiscretizationN"], 30)
        self.assertEqual(result["metadata"]["sampleDiscretizationM"], 30)
        self.assertEqual(result["metadata"]["referenceDiscretizationN"], 30)

    def test_fractional_harmonic_is_rejected(self) -> None:
        payload = build_payload()
        payload["harmonicNumber"] = 3.5

        with self.assertRaisesRegex(ValidationError, "positive integer"):
            calculate_spectrum(payload)

    def test_drude_lorentz_generator_matches_formula(self) -> None:
        frequency = np.array([1000.0, 1200.0])
        parameters = {
            "epsilonInfinity": 2.0,
            "useDrude": True,
            "plasmaFrequency": 500.0,
            "drudeDamping": 25.0,
            "oscillators": [
                {"strength": 1.2, "resonanceFrequency": 1500.0, "damping": 30.0},
            ],
        }

        epsilon = generate_epsilon_drude_lorentz(frequency, parameters)
        expected = (
            2.0
            - 500.0**2 / (frequency**2 + 1j * 25.0 * frequency)
            + 1.2 * 1500.0**2 / (1500.0**2 - frequency**2 - 1j * 30.0 * frequency)
        )
        np.testing.assert_allclose(epsilon, expected)

    def test_drude_lorentz_sample_calculates(self) -> None:
        payload = build_payload()
        payload["sample"] = {
            "inputMethod": "drudeLorentz",
            "builtInMaterial": None,
            "uploadedPermittivity": None,
            "drudeLorentz": {
                "epsilonInfinity": 1.0,
                "useDrude": False,
                "plasmaFrequency": 0.0,
                "drudeDamping": 0.0,
                "oscillators": [
                    {"strength": 1.0, "resonanceFrequency": 1700.0, "damping": 20.0},
                ],
            },
        }

        result = calculate_spectrum(payload)
        self.assertEqual(result["metadata"]["materialLabel"], "Drude-Lorentz model")
        self.assertEqual(len(result["epsilonReal"]), 4)

    def test_server_calculation_returns_more_than_pipe_buffer_threshold(self) -> None:
        payload = build_payload()
        payload["frequencyRange"] = {
            "start": 1650,
            "end": 1800,
            "points": 130,
        }

        with mock.patch("snom_web.server.TIMEOUT_SECONDS", 20):
            result = _run_calculation_with_timeout(payload, "large-payload-test")

        self.assertEqual(len(result["frequency"]), 130)

    def test_single_layer_calculation_returns_layer_permittivities(self) -> None:
        result = calculate_spectrum(build_layered_payload("single"))

        self.assertEqual(result["metadata"]["modelType"], "layered")
        self.assertEqual(result["metadata"]["layeredStructure"], "single")
        self.assertEqual(result["metadata"]["imageChargeCount"], 2)
        self.assertEqual(len(result["frequency"]), 3)
        self.assertEqual(len(result["permittivitySeries"]), 2)
        self.assertTrue(all(math.isfinite(value) for value in result["amplitude"]))

    def test_double_layer_calculation_returns_all_permittivities(self) -> None:
        result = calculate_spectrum(build_layered_payload("double"))

        self.assertEqual(result["metadata"]["layeredStructure"], "double")
        self.assertEqual(len(result["permittivitySeries"]), 3)
        self.assertEqual(result["permittivitySeries"][0]["key"], "layer1")
        self.assertEqual(result["permittivitySeries"][1]["key"], "layer2")
        self.assertEqual(result["permittivitySeries"][2]["key"], "substrate")

    def test_layered_image_charge_limit_is_rejected(self) -> None:
        payload = build_layered_payload("single")
        payload["layered"]["imageCharges"] = 51

        with self.assertRaisesRegex(ValidationError, "K must be between 1 and 50"):
            calculate_spectrum(payload)


if __name__ == "__main__":
    unittest.main()
