from __future__ import annotations

import math
import re
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss

from .materials import (
    generate_epsilon_drude_lorentz,
    get_builtin_permittivity,
    get_material_label,
    interpolate_complex_data,
)


MAX_DISCRETIZATION_N = 500
MAX_DISCRETIZATION_M = 100
MAX_LAYER_COUNT = 10
MAX_FREQUENCY_POINTS = 1000
DEFAULT_FAST_DISCRETIZATION = 10
GAUSS_LEGENDRE_ORDER = 64


class ValidationError(ValueError):
    """User-facing request validation error."""


@dataclass(frozen=True)
class GeometrySolver:
    a: float
    c: float
    harmonic_weights: np.ndarray
    eigenvalues: np.ndarray
    modal_numerators: np.ndarray
    tip_offset: complex

    def solve(self, beta: complex, field_multiplier: complex) -> complex:
        denominator = 1.0 - beta * self.eigenvalues
        dipole_first_component = beta * field_multiplier * np.sum(
            self.modal_numerators / denominator, axis=1
        )
        dipole_profile = 2.0 * self.a * self.c * dipole_first_component + self.tip_offset * field_multiplier
        return np.sum(dipole_profile * self.harmonic_weights)


def calculate_spectrum(payload: dict[str, object]) -> dict[str, object]:
    config = _validate_payload(payload)

    frequency_grid = np.linspace(
        config["frequency_start"],
        config["frequency_end"],
        config["frequency_points"],
        dtype=float,
    )

    sample_permittivity = _resolve_sample_permittivity(config, frequency_grid)
    experimental_spectrum = _resolve_experimental_spectrum(config)

    if config["algorithm"] == "fast":
        sample_discretization_n = DEFAULT_FAST_DISCRETIZATION
        sample_discretization_m = DEFAULT_FAST_DISCRETIZATION
    else:
        sample_discretization_n = config["discretization_n"]
        sample_discretization_m = config["discretization_m"]

    reference_solver = _build_geometry_solver(
        radius=config["tip_radius"],
        length=config["tip_length"],
        amplitude=config["tip_amplitude"],
        minimum_gap=config["tip_reference_gap"],
        tip_permittivity=complex(config["tip_epsilon_real"], config["tip_epsilon_imag"]),
        harmonic_number=config["harmonic_number"],
        discretization_n=DEFAULT_FAST_DISCRETIZATION,
        matrix_size=DEFAULT_FAST_DISCRETIZATION,
    )
    sample_solver = _build_geometry_solver(
        radius=config["tip_radius"],
        length=config["tip_length"],
        amplitude=config["tip_amplitude"],
        minimum_gap=config["tip_sample_gap"],
        tip_permittivity=complex(config["tip_epsilon_real"], config["tip_epsilon_imag"]),
        harmonic_number=config["harmonic_number"],
        discretization_n=sample_discretization_n,
        matrix_size=sample_discretization_m,
    )

    reference_permittivity = complex(
        config["reference_epsilon_real"], config["reference_epsilon_imag"]
    )
    reference_beta = (reference_permittivity - 1.0) / (reference_permittivity + 1.0)
    reference_field = _far_field_multiplier(
        epsilon=reference_permittivity,
        enabled=config["far_field_enabled"],
        coefficient=config["far_field_coefficient"],
        angle_degrees=config["far_field_angle"],
    )
    reference_harmonic = reference_solver.solve(reference_beta, reference_field)

    if abs(reference_harmonic) < 1e-18:
        raise ValidationError("The selected reference material produced a zero normalization signal.")

    sigma_values = np.empty(frequency_grid.shape[0], dtype=np.complex128)
    for index, epsilon in enumerate(sample_permittivity):
        beta = (epsilon - 1.0) / (epsilon + 1.0)
        field_multiplier = _far_field_multiplier(
            epsilon=epsilon,
            enabled=config["far_field_enabled"],
            coefficient=config["far_field_coefficient"],
            angle_degrees=config["far_field_angle"],
        )
        sample_harmonic = sample_solver.solve(beta, field_multiplier)
        sigma_values[index] = sample_harmonic / reference_harmonic

    phase_values = _unwrap_phase(sigma_values)
    amplitude_values = np.abs(sigma_values)

    return {
        "frequency": frequency_grid.tolist(),
        "epsilonReal": sample_permittivity.real.tolist(),
        "epsilonImag": sample_permittivity.imag.tolist(),
        "sigmaReal": sigma_values.real.tolist(),
        "sigmaImag": sigma_values.imag.tolist(),
        "amplitude": amplitude_values.tolist(),
        "phase": phase_values.tolist(),
        "experimental": experimental_spectrum,
        "metadata": {
            "modelType": config["model_type"],
            "algorithm": config["algorithm"],
            "harmonicNumber": config["harmonic_number"],
            "farFieldEnabled": config["far_field_enabled"],
            "sampleInputMethod": config["sample_input_method"],
            "materialLabel": config["sample_material_label"],
            "frequencyStart": config["frequency_start"],
            "frequencyEnd": config["frequency_end"],
            "frequencyPoints": config["frequency_points"],
            "sampleDiscretizationN": sample_discretization_n,
            "sampleDiscretizationM": sample_discretization_m,
            "referenceDiscretizationN": DEFAULT_FAST_DISCRETIZATION,
            "referenceDiscretizationM": DEFAULT_FAST_DISCRETIZATION,
        },
        "suggestedFilename": _build_suggested_filename(config),
    }


def parse_numeric_table(
    text: str,
    *,
    minimum_columns: int,
    label: str,
) -> np.ndarray:
    rows: list[list[float]] = []

    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        columns = [token for token in re.split(r"[\s,;]+", line) if token]
        if len(columns) < minimum_columns:
            raise ValidationError(
                f"{label}: line {line_number} must contain at least {minimum_columns} numeric columns."
            )

        try:
            rows.append([float(columns[index]) for index in range(minimum_columns)])
        except ValueError as error:
            raise ValidationError(
                f"{label}: line {line_number} contains a non-numeric value."
            ) from error

    if len(rows) < 2:
        raise ValidationError(f"{label}: at least two data rows are required.")

    return np.asarray(rows, dtype=float)


def _validate_payload(payload: dict[str, object]) -> dict[str, object]:
    if not isinstance(payload, dict):
        raise ValidationError("Request payload must be a JSON object.")

    model_type = payload.get("modelType")
    if model_type != "bulk":
        raise ValidationError("Layered sample mode will be added in the next stage.")

    algorithm = payload.get("algorithm")
    if algorithm not in {"fast", "accurate"}:
        raise ValidationError("Calculation algorithm must be either Fast or Accurate.")

    tip = _require_mapping(payload, "tip")
    far_field = _require_mapping(payload, "farField")
    reference = _require_mapping(payload, "reference")
    frequency_range = _require_mapping(payload, "frequencyRange")
    sample = _require_mapping(payload, "sample")
    experimental = _optional_mapping(payload.get("experimentalSpectrum"))
    discretization = _require_mapping(payload, "discretization")

    tip_radius = _require_positive_number(tip, "radius", "R must be greater than 0.")
    tip_length = _require_positive_number(tip, "length", "L must be greater than 0.")
    tip_amplitude = _require_positive_number(tip, "amplitude", "A must be greater than 0.")
    tip_sample_gap = _require_positive_number(
        tip, "sampleGap", "H0,s must be greater than 0."
    )
    tip_reference_gap = _require_positive_number(
        tip, "referenceGap", "H0,r must be greater than 0."
    )

    if tip_length <= 2.0 * tip_radius:
        raise ValidationError("L must be greater than 2R so the spheroid geometry remains valid.")

    tip_epsilon_real = _require_number(tip, "epsilonReal", "Tip permittivity real part is required.")
    tip_epsilon_imag = _require_number(tip, "epsilonImag", "Tip permittivity imaginary part is required.")

    far_field_enabled = bool(far_field.get("enabled", False))
    far_field_coefficient = _require_number(
        far_field,
        "coefficient",
        "Far-field coefficient is required.",
    )
    far_field_angle = _require_number(
        far_field,
        "angleDegrees",
        "Illumination angle is required.",
    )
    if far_field_angle < 0 or far_field_angle > 90:
        raise ValidationError("Illumination angle must be between 0 and 90 degrees.")

    reference_epsilon_real = _require_number(
        reference,
        "epsilonReal",
        "Reference permittivity real part is required.",
    )
    reference_epsilon_imag = _require_number(
        reference,
        "epsilonImag",
        "Reference permittivity imaginary part is required.",
    )

    harmonic_number = _require_positive_integer(payload, "harmonicNumber")

    frequency_start = _require_number(
        frequency_range,
        "start",
        "Start frequency is required.",
    )
    frequency_end = _require_number(
        frequency_range,
        "end",
        "End frequency is required.",
    )
    if frequency_start >= frequency_end:
        raise ValidationError("Start frequency must be smaller than end frequency.")

    frequency_points = _require_positive_integer(frequency_range, "points")
    if frequency_points > MAX_FREQUENCY_POINTS:
        raise ValidationError(
            f"Number of frequency points must be between 1 and {MAX_FREQUENCY_POINTS}."
        )

    sample_input_method = sample.get("inputMethod")
    if sample_input_method not in {"builtIn", "upload", "drudeLorentz"}:
        raise ValidationError("Sample material input method is invalid.")
    if sample_input_method == "drudeLorentz":
        raise ValidationError("Drude-Lorentz material input will be added in the next stage.")

    sample_material_id = sample.get("builtInMaterial") if sample_input_method == "builtIn" else None
    sample_material_label = "Uploaded file"
    if sample_input_method == "builtIn":
        if not isinstance(sample_material_id, str) or not sample_material_id:
            raise ValidationError("Please choose a built-in material.")
        try:
            sample_material_label = get_material_label(sample_material_id)
        except ValueError as error:
            raise ValidationError("Unknown built-in material selected.") from error

    uploaded_permittivity = sample.get("uploadedPermittivity")
    if sample_input_method == "upload" and not isinstance(uploaded_permittivity, str):
        raise ValidationError("Please upload a permittivity file.")

    experimental_enabled = bool(experimental.get("enabled", False))
    experimental_content = experimental.get("content")
    if experimental_enabled and not isinstance(experimental_content, str):
        raise ValidationError("Please upload an experimental spectrum file.")

    discretization_n = _require_positive_integer(discretization, "N")
    discretization_m = _require_positive_integer(discretization, "M")
    if discretization_n < 1 or discretization_n > MAX_DISCRETIZATION_N:
        raise ValidationError(f"N must be between 1 and {MAX_DISCRETIZATION_N}.")
    if discretization_m < 1 or discretization_m > MAX_DISCRETIZATION_M:
        raise ValidationError(f"M must be between 1 and {MAX_DISCRETIZATION_M}.")

    return {
        "model_type": model_type,
        "algorithm": algorithm,
        "tip_radius": tip_radius,
        "tip_length": tip_length,
        "tip_amplitude": tip_amplitude,
        "tip_sample_gap": tip_sample_gap,
        "tip_reference_gap": tip_reference_gap,
        "tip_epsilon_real": tip_epsilon_real,
        "tip_epsilon_imag": tip_epsilon_imag,
        "far_field_enabled": far_field_enabled,
        "far_field_coefficient": far_field_coefficient,
        "far_field_angle": far_field_angle,
        "reference_epsilon_real": reference_epsilon_real,
        "reference_epsilon_imag": reference_epsilon_imag,
        "harmonic_number": harmonic_number,
        "frequency_start": frequency_start,
        "frequency_end": frequency_end,
        "frequency_points": frequency_points,
        "sample_input_method": sample_input_method,
        "sample_material_id": sample_material_id,
        "sample_material_label": sample_material_label,
        "uploaded_permittivity": uploaded_permittivity,
        "experimental_enabled": experimental_enabled,
        "experimental_content": experimental_content,
        "discretization_n": discretization_n,
        "discretization_m": discretization_m,
    }


def _resolve_sample_permittivity(
    config: dict[str, object], frequency_grid: np.ndarray
) -> np.ndarray:
    input_method = config["sample_input_method"]

    if input_method == "builtIn":
        return get_builtin_permittivity(str(config["sample_material_id"]), frequency_grid)

    if input_method == "upload":
        table = parse_numeric_table(
            str(config["uploaded_permittivity"]),
            minimum_columns=3,
            label="Permittivity file",
        )
        epsilon_values = table[:, 1] + 1j * table[:, 2]
        return interpolate_complex_data(
            table[:, 0],
            epsilon_values,
            frequency_grid,
            range_error_message="The selected frequency range is outside the uploaded permittivity data range.",
        )

    generate_epsilon_drude_lorentz(frequency_grid, {})
    raise ValidationError("Unsupported sample material input method.")


def _resolve_experimental_spectrum(config: dict[str, object]) -> dict[str, list[float]] | None:
    if not config["experimental_enabled"]:
        return None

    table = parse_numeric_table(
        str(config["experimental_content"]),
        minimum_columns=3,
        label="Experimental spectrum file",
    )
    order = np.argsort(table[:, 0])
    table = table[order]
    return {
        "frequency": table[:, 0].tolist(),
        "amplitude": table[:, 1].tolist(),
        "phase": table[:, 2].tolist(),
    }


def _build_suggested_filename(config: dict[str, object]) -> str:
    material_token = str(config["sample_material_label"]).lower().replace(" ", "-")
    return (
        f"s-snom-{material_token}-{config['algorithm']}-"
        f"n{config['harmonic_number']}-{int(config['frequency_start'])}-{int(config['frequency_end'])}.csv"
    )


def _unwrap_phase(sigma_values: np.ndarray) -> np.ndarray:
    phase = np.angle(sigma_values)
    phase = np.where(phase < 0, phase + 2.0 * np.pi, phase)
    for index in range(phase.size - 1):
        delta = phase[index] - phase[index + 1]
        if abs(delta) > 3.0:
            phase[index + 1 :] += math.copysign(2.0 * np.pi, delta)
    return phase


def _far_field_multiplier(
    *, epsilon: complex, enabled: bool, coefficient: float, angle_degrees: float
) -> complex:
    if not enabled:
        return 1.0 + 0.0j

    phi = math.radians(angle_degrees)
    cosine_phi = math.cos(phi)
    if abs(cosine_phi) < 1e-12:
        cosine_phi = 1e-12

    kz1 = cosine_phi
    kz2 = np.sqrt(epsilon - math.sin(phi) ** 2)
    rp = (epsilon / kz2 - 1.0 / kz1) / (epsilon / kz2 + 1.0 / kz1)
    return (1.0 + coefficient * rp) ** 2


@lru_cache(maxsize=32)
def _build_geometry_solver(
    *,
    radius: float,
    length: float,
    amplitude: float,
    minimum_gap: float,
    tip_permittivity: complex,
    harmonic_number: int,
    discretization_n: int,
    matrix_size: int,
) -> GeometrySolver:
    a = length / 2.0
    c = math.sqrt(a * (a - radius))
    axis_ratio = a / c

    nodes, weights = leggauss(GAUSS_LEGENDRE_ORDER)
    p_nodes = _legendre_p_table(matrix_size, nodes)
    weighted_p_nodes = p_nodes[1:] * weights[np.newaxis, :]

    p_axis = _legendre_p_table(matrix_size, np.array([axis_ratio], dtype=float))[:, 0]
    q_axis = _legendre_q_table(matrix_size, np.array([axis_ratio], dtype=float))[:, 0]

    denominator_terms = np.empty(matrix_size + 1, dtype=np.complex128)
    denominator_terms[0] = np.nan + 0j
    for order in range(1, matrix_size + 1):
        denominator_terms[order] = tip_permittivity * q_axis[order] / p_axis[order] - (
            axis_ratio * q_axis[order] - q_axis[order - 1]
        ) / (axis_ratio * p_axis[order] - p_axis[order - 1])

    reference_constant = (
        tip_permittivity * q_axis[1] - axis_ratio * q_axis[0] + a * a / (a * a - c * c)
    )
    tip_offset = -(tip_permittivity - 1.0) * a * c * c / (3.0 * reference_constant)

    prefactors = (tip_permittivity - 1.0) * (2.0 * np.arange(1, matrix_size + 1) + 1.0) / 2.0
    geometry_denominator = p_axis[1:] * denominator_terms[1:]

    eigenvalues = np.empty((discretization_n + 1, matrix_size), dtype=np.complex128)
    modal_numerators = np.empty_like(eigenvalues)
    harmonic_weights = _harmonic_weights(harmonic_number, discretization_n)

    h0 = minimum_gap + amplitude
    for k_index in range(discretization_n + 1):
        psi = math.pi * k_index / discretization_n
        distance = h0 - amplitude * math.cos(psi)
        z_shift = 2.0 * a / c + 2.0 * distance / c
        numerator = z_shift - axis_ratio * nodes
        base = 1.0 + numerator * numerator + (a * a / (c * c) - 1.0) * (1.0 - nodes * nodes)
        radical = np.maximum(base * base - 4.0 * numerator * numerator, 0.0)
        radial_argument = np.sqrt((base + np.sqrt(radical)) / 2.0)
        transformed_argument = numerator / radial_argument

        p_argument = _legendre_p_table(matrix_size, transformed_argument)
        q_radial = _legendre_q_table(matrix_size, radial_argument)
        right_matrix = p_argument[1:] * q_radial[1:]
        j_matrix = weighted_p_nodes @ right_matrix.T

        source_vector = (
            -(tip_permittivity - 1.0) ** 2
            * c
            / 4.0
            * j_matrix[0, :]
            / (p_axis[1:] * reference_constant * denominator_terms[1:])
        )
        system_matrix = (j_matrix.T / geometry_denominator[:, np.newaxis]) * prefactors[np.newaxis, :]

        spectral_values, spectral_vectors = np.linalg.eig(system_matrix)
        modal_projection = np.linalg.solve(spectral_vectors, source_vector)

        eigenvalues[k_index, :] = spectral_values
        modal_numerators[k_index, :] = spectral_vectors[0, :] * modal_projection

    return GeometrySolver(
        a=a,
        c=c,
        harmonic_weights=harmonic_weights,
        eigenvalues=eigenvalues,
        modal_numerators=modal_numerators,
        tip_offset=tip_offset,
    )


def _harmonic_weights(harmonic_number: int, discretization_n: int) -> np.ndarray:
    weights = np.empty(discretization_n + 1, dtype=np.complex128)
    for k_index in range(discretization_n + 1):
        if k_index in {0, discretization_n}:
            weights[k_index] = np.exp(-1j * np.pi * k_index * harmonic_number / discretization_n) / (
                2.0 * discretization_n
            )
        else:
            weights[k_index] = math.cos(np.pi * k_index * harmonic_number / discretization_n) / discretization_n
    return weights


def _legendre_p_table(max_order: int, x_values: np.ndarray) -> np.ndarray:
    x_array = np.asarray(x_values, dtype=np.complex128)
    table = np.empty((max_order + 1, x_array.size), dtype=np.complex128)
    table[0, :] = 1.0
    if max_order >= 1:
        table[1, :] = x_array
    for order in range(1, max_order):
        table[order + 1, :] = (
            ((2.0 * order + 1.0) * x_array * table[order, :]) - order * table[order - 1, :]
        ) / (order + 1.0)
    return table


def _legendre_q_table(max_order: int, x_values: np.ndarray) -> np.ndarray:
    x_array = np.asarray(x_values, dtype=np.complex128)
    p_table = _legendre_p_table(max_order, x_array)
    table = np.empty((max_order + 1, x_array.size), dtype=np.complex128)
    table[0, :] = 0.5 * np.log(np.abs((x_array + 1.0) / (x_array - 1.0)))
    if max_order >= 1:
        table[1, :] = table[0, :] * x_array - 1.0
    for order in range(2, max_order + 1):
        table[order, :] = _legendre_q_single(order, x_array, p_table[order, :], p_table[order - 1, :])
    return table


def _legendre_q_single(
    order: int,
    x_values: np.ndarray,
    p_order: np.ndarray,
    p_previous: np.ndarray,
) -> np.ndarray:
    epsilon = 1e-15
    tiny = 1e-30
    h_n = np.full(x_values.shape, tiny, dtype=np.complex128)
    c_previous = h_n.copy()
    d_previous = np.zeros(x_values.shape, dtype=np.complex128)

    a_j = 1.0
    b_j = (2.0 + 1.0 / order) * x_values
    d_current = b_j + a_j * d_previous
    d_current = np.where(d_current == 0, tiny, d_current)
    c_current = b_j + a_j / c_previous
    c_current = np.where(c_current == 0, tiny, c_current)
    d_current = 1.0 / d_current
    delta = c_current * d_current
    h_n *= delta
    d_previous = d_current
    c_previous = c_current

    iteration = 0
    while np.max(np.abs(delta - 1.0)) >= epsilon:
        a_j = -(1.0 + 1.0 / (order + iteration))
        b_j = (2.0 + 1.0 / (order + iteration + 1.0)) * x_values
        d_current = b_j + a_j * d_previous
        d_current = np.where(d_current == 0, tiny, d_current)
        c_current = b_j + a_j / c_previous
        c_current = np.where(c_current == 0, tiny, c_current)
        d_current = 1.0 / d_current
        delta = c_current * d_current
        h_n *= delta
        d_previous = d_current
        c_previous = c_current
        iteration += 1
        if iteration > 512:
            break

    return h_n / order / (p_order - h_n * p_previous)


def _require_mapping(payload: dict[str, object], key: str) -> dict[str, object]:
    value = payload.get(key)
    if not isinstance(value, dict):
        raise ValidationError(f"'{key}' must be an object.")
    return value


def _optional_mapping(value: object) -> dict[str, object]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValidationError("'experimentalSpectrum' must be an object.")
    return value


def _require_number(payload: dict[str, object], key: str, message: str) -> float:
    value = payload.get(key)
    if value is None:
        raise ValidationError(message)
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValidationError(message) from error
    if not math.isfinite(number):
        raise ValidationError(message)
    return number


def _require_positive_number(payload: dict[str, object], key: str, message: str) -> float:
    value = _require_number(payload, key, message)
    if value <= 0:
        raise ValidationError(message)
    return value


def _require_positive_integer(payload: dict[str, object], key: str) -> int:
    value = payload.get(key)
    try:
        integer = int(value)
    except (TypeError, ValueError) as error:
        raise ValidationError(f"'{key}' must be a positive integer.") from error

    if float(value) != integer or integer <= 0:
        raise ValidationError(f"'{key}' must be a positive integer.")
    return integer
