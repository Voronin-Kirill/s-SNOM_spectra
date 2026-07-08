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
MAX_DRUDE_LORENTZ_OSCILLATORS = 10
MAX_IMAGE_CHARGES = 50
DEFAULT_FAST_DISCRETIZATION = 30
GAUSS_LEGENDRE_ORDER = 64
LAYERED_REFERENCE_GAUSS_ORDER = 48
LAYERED_SAMPLE_GAUSS_ORDER = 64


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


@dataclass(frozen=True)
class DirectGeometry:
    a: float
    c: float
    harmonic_weights: np.ndarray
    j_tables: np.ndarray
    image_tables: np.ndarray
    row_scale: np.ndarray
    col_scale: np.ndarray
    const_c: complex
    const_pz: complex


def calculate_spectrum(payload: dict[str, object]) -> dict[str, object]:
    config = _validate_payload(payload)

    frequency_grid = np.linspace(
        config["frequency_start"],
        config["frequency_end"],
        config["frequency_points"],
        dtype=float,
    )

    if config["model_type"] == "layered":
        return _calculate_layered_spectrum(config, frequency_grid)
    return _calculate_bulk_spectrum(config, frequency_grid)


def _calculate_bulk_spectrum(
    config: dict[str, object], frequency_grid: np.ndarray
) -> dict[str, object]:
    try:
        sample_permittivity = _resolve_sample_permittivity(config, frequency_grid)
    except ValueError as error:
        raise ValidationError(str(error)) from error
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

    return _build_result_payload(
        config=config,
        frequency_grid=frequency_grid,
        sigma_values=sigma_values,
        experimental_spectrum=experimental_spectrum,
        permittivity_series=[
            {
                "key": "sample",
                "label": str(config["sample_material_label"]),
                "epsilon": sample_permittivity,
            }
        ],
        metadata={
            "sampleInputMethod": config["sample_input_method"],
            "materialLabel": config["sample_material_label"],
            "sampleDiscretizationN": sample_discretization_n,
            "sampleDiscretizationM": sample_discretization_m,
            "referenceDiscretizationN": DEFAULT_FAST_DISCRETIZATION,
            "referenceDiscretizationM": DEFAULT_FAST_DISCRETIZATION,
        },
    )


def _calculate_layered_spectrum(
    config: dict[str, object], frequency_grid: np.ndarray
) -> dict[str, object]:
    layered = dict(config["layered"])
    layer1 = dict(layered["layer1"])
    layer1_material = dict(layer1["material"])
    layer2 = dict(layered["layer2"]) if layered["layer2"] is not None else None
    layer2_material = dict(layer2["material"]) if layer2 is not None else None
    substrate_material = dict(layered["substrate"])

    try:
        layer1_permittivity = _resolve_material_permittivity(
            layer1_material, frequency_grid, "Layer 1"
        )
        layer2_permittivity = (
            _resolve_material_permittivity(layer2_material, frequency_grid, "Layer 2")
            if layer2_material is not None
            else None
        )
        substrate_permittivity = _resolve_material_permittivity(
            substrate_material, frequency_grid, "Substrate"
        )
    except ValueError as error:
        raise ValidationError(str(error)) from error

    experimental_spectrum = _resolve_experimental_spectrum(config)

    if config["algorithm"] == "fast":
        sample_discretization_n = DEFAULT_FAST_DISCRETIZATION
        sample_discretization_m = DEFAULT_FAST_DISCRETIZATION
    else:
        sample_discretization_n = config["discretization_n"]
        sample_discretization_m = config["discretization_m"]

    reference_geometry = _build_direct_geometry(
        radius=config["tip_radius"],
        length=config["tip_length"],
        amplitude=config["tip_amplitude"],
        minimum_gap=config["tip_reference_gap"],
        tip_permittivity=complex(config["tip_epsilon_real"], config["tip_epsilon_imag"]),
        harmonic_number=config["harmonic_number"],
        discretization_n=DEFAULT_FAST_DISCRETIZATION,
        matrix_size=DEFAULT_FAST_DISCRETIZATION,
        image_shifts=(),
        gauss_order=LAYERED_REFERENCE_GAUSS_ORDER,
    )
    reference_permittivity = complex(
        config["reference_epsilon_real"], config["reference_epsilon_imag"]
    )
    reference_harmonic = _solve_direct_harmonic(
        reference_geometry,
        beta=(reference_permittivity - 1.0) / (reference_permittivity + 1.0),
        field_multiplier=_far_field_multiplier(
            epsilon=reference_permittivity,
            enabled=config["far_field_enabled"],
            coefficient=config["far_field_coefficient"],
            angle_degrees=config["far_field_angle"],
        ),
    )
    if abs(reference_harmonic) < 1e-18:
        raise ValidationError("The selected reference material produced a zero normalization signal.")

    image_charge_count = int(layered["image_charge_count"])
    if layered["structure"] == "single":
        d1 = float(layer1["thickness"])
        image_shifts = tuple(2.0 * d1 * charge_index for charge_index in range(1, image_charge_count + 1))
    else:
        d1 = float(layer1["thickness"])
        d2 = float(layer2["thickness"])
        image_shifts = (2.0 * d1,) + tuple(
            2.0 * d1 + 2.0 * min(d1, d2) * charge_index
            for charge_index in range(1, image_charge_count)
        )

    sample_geometry = _build_direct_geometry(
        radius=config["tip_radius"],
        length=config["tip_length"],
        amplitude=config["tip_amplitude"],
        minimum_gap=config["tip_sample_gap"],
        tip_permittivity=complex(config["tip_epsilon_real"], config["tip_epsilon_imag"]),
        harmonic_number=config["harmonic_number"],
        discretization_n=sample_discretization_n,
        matrix_size=sample_discretization_m,
        image_shifts=image_shifts,
        gauss_order=LAYERED_SAMPLE_GAUSS_ORDER,
    )

    sigma_values = np.empty(frequency_grid.shape[0], dtype=np.complex128)
    for index, frequency in enumerate(frequency_grid):
        epsilon_layer1 = layer1_permittivity[index]
        epsilon_substrate = substrate_permittivity[index]
        beta1 = (epsilon_layer1 - 1.0) / (epsilon_layer1 + 1.0)

        if layered["structure"] == "single":
            image_charges = _single_layer_image_charges(
                1.0 + 0.0j,
                epsilon_layer1,
                epsilon_substrate,
                d1,
                image_charge_count,
            )
            field_multiplier = _single_layer_far_field_multiplier(
                epsilon_layer=epsilon_layer1,
                epsilon_substrate=epsilon_substrate,
                thickness=d1,
                frequency=frequency,
                enabled=config["far_field_enabled"],
                coefficient=config["far_field_coefficient"],
                angle_degrees=config["far_field_angle"],
            )
        else:
            epsilon_layer2 = layer2_permittivity[index]
            image_charges = _double_layer_image_charges(
                1.0 + 0.0j,
                epsilon_layer1,
                epsilon_layer2,
                epsilon_substrate,
                d1,
                d2,
                image_charge_count,
            )
            field_multiplier = _double_layer_far_field_multiplier(
                epsilon_layer1=epsilon_layer1,
                epsilon_layer2=epsilon_layer2,
                epsilon_substrate=epsilon_substrate,
                thickness1=d1,
                thickness2=d2,
                frequency=frequency,
                enabled=config["far_field_enabled"],
                coefficient=config["far_field_coefficient"],
                angle_degrees=config["far_field_angle"],
            )

        sample_harmonic = _solve_layered_harmonic(
            sample_geometry,
            beta1=beta1,
            image_charges=image_charges,
            field_multiplier=field_multiplier,
        )
        sigma_values[index] = sample_harmonic / reference_harmonic

    permittivity_series = [
        {
            "key": "layer1",
            "label": f"Layer 1 ({layer1_material['material_label']})",
            "epsilon": layer1_permittivity,
        }
    ]
    if layer2_permittivity is not None:
        permittivity_series.append(
            {
                "key": "layer2",
                "label": f"Layer 2 ({layer2_material['material_label']})",
                "epsilon": layer2_permittivity,
            }
        )
    permittivity_series.append(
        {
            "key": "substrate",
            "label": f"Substrate ({substrate_material['material_label']})",
            "epsilon": substrate_permittivity,
        }
    )

    structure_label = "1 layer on substrate" if layered["structure"] == "single" else "2 layers on substrate"
    return _build_result_payload(
        config=config,
        frequency_grid=frequency_grid,
        sigma_values=sigma_values,
        experimental_spectrum=experimental_spectrum,
        permittivity_series=permittivity_series,
        metadata={
            "layeredStructure": layered["structure"],
            "materialLabel": structure_label,
            "imageChargeCount": image_charge_count,
            "sampleDiscretizationN": sample_discretization_n,
            "sampleDiscretizationM": sample_discretization_m,
            "referenceDiscretizationN": DEFAULT_FAST_DISCRETIZATION,
            "referenceDiscretizationM": DEFAULT_FAST_DISCRETIZATION,
        },
    )


def _build_result_payload(
    *,
    config: dict[str, object],
    frequency_grid: np.ndarray,
    sigma_values: np.ndarray,
    experimental_spectrum: dict[str, list[float]] | None,
    permittivity_series: list[dict[str, object]],
    metadata: dict[str, object],
) -> dict[str, object]:
    phase_values = _unwrap_phase(sigma_values)
    amplitude_values = np.abs(sigma_values)
    primary_permittivity = np.asarray(permittivity_series[0]["epsilon"], dtype=np.complex128)
    series_payload = []
    for series in permittivity_series:
        epsilon = np.asarray(series["epsilon"], dtype=np.complex128)
        series_payload.append(
            {
                "key": series["key"],
                "label": series["label"],
                "epsilonReal": epsilon.real.tolist(),
                "epsilonImag": epsilon.imag.tolist(),
            }
        )

    base_metadata = {
        "modelType": config["model_type"],
        "algorithm": config["algorithm"],
        "harmonicNumber": config["harmonic_number"],
        "farFieldEnabled": config["far_field_enabled"],
        "frequencyStart": config["frequency_start"],
        "frequencyEnd": config["frequency_end"],
        "frequencyPoints": config["frequency_points"],
    }
    base_metadata.update(metadata)

    return {
        "frequency": frequency_grid.tolist(),
        "epsilonReal": primary_permittivity.real.tolist(),
        "epsilonImag": primary_permittivity.imag.tolist(),
        "permittivitySeries": series_payload,
        "sigmaReal": sigma_values.real.tolist(),
        "sigmaImag": sigma_values.imag.tolist(),
        "amplitude": amplitude_values.tolist(),
        "phase": phase_values.tolist(),
        "experimental": experimental_spectrum,
        "metadata": base_metadata,
        "suggestedFilename": _build_suggested_filename(config),
    }


def _solve_direct_harmonic(
    geometry: DirectGeometry, *, beta: complex, field_multiplier: complex
) -> complex:
    return _solve_harmonic_from_wtabs(
        geometry,
        field_multiplier=field_multiplier,
        wtab_builder=lambda height_index: -beta * geometry.j_tables[height_index],
    )


def _solve_layered_harmonic(
    geometry: DirectGeometry,
    *,
    beta1: complex,
    image_charges: np.ndarray,
    field_multiplier: complex,
) -> complex:
    def build_wtab(height_index: int) -> np.ndarray:
        wtab = -beta1 * geometry.j_tables[height_index]
        for charge_index, charge in enumerate(image_charges):
            wtab = wtab + charge * geometry.image_tables[height_index, charge_index]
        return wtab

    return _solve_harmonic_from_wtabs(
        geometry,
        field_multiplier=field_multiplier,
        wtab_builder=build_wtab,
    )


def _solve_harmonic_from_wtabs(
    geometry: DirectGeometry,
    *,
    field_multiplier: complex,
    wtab_builder,
) -> complex:
    matrix_size = geometry.row_scale.size
    identity = np.eye(matrix_size, dtype=np.complex128)
    dipole_profile = np.empty(geometry.harmonic_weights.size, dtype=np.complex128)

    for height_index in range(geometry.harmonic_weights.size):
        wtab = np.asarray(wtab_builder(height_index), dtype=np.complex128)
        right_side = geometry.const_c * field_multiplier * (wtab[0, :] * geometry.row_scale)
        system_matrix = identity + (wtab.T * geometry.col_scale[np.newaxis, :]) * geometry.row_scale[:, np.newaxis]
        coefficients = np.linalg.solve(system_matrix, right_side)
        dipole_profile[height_index] = 2.0 * geometry.a * geometry.c * coefficients[0] - geometry.const_pz * field_multiplier

    return np.sum(dipole_profile * geometry.harmonic_weights)


@lru_cache(maxsize=16)
def _build_direct_geometry(
    *,
    radius: float,
    length: float,
    amplitude: float,
    minimum_gap: float,
    tip_permittivity: complex,
    harmonic_number: int,
    discretization_n: int,
    matrix_size: int,
    image_shifts: tuple[float, ...],
    gauss_order: int,
) -> DirectGeometry:
    a = length / 2.0
    c = math.sqrt(a * (a - radius))
    axis_ratio = a / c

    nodes, weights = leggauss(gauss_order)
    p_nodes = _legendre_p_table(matrix_size, nodes)[1:]
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
    row_scale = 1.0 / (p_axis[1:] * denominator_terms[1:])
    col_scale = (tip_permittivity - 1.0) * (2.0 * np.arange(1, matrix_size + 1) + 1.0) / 2.0
    const_c = (tip_permittivity - 1.0) ** 2 * c / 4.0 / reference_constant
    const_pz = (tip_permittivity - 1.0) * a * c * c / 3.0 / reference_constant
    harmonic_weights = _harmonic_weights(harmonic_number, discretization_n)

    j_tables = np.empty((discretization_n + 1, matrix_size, matrix_size), dtype=float)
    image_tables = np.empty(
        (discretization_n + 1, len(image_shifts), matrix_size, matrix_size),
        dtype=float,
    )
    h0 = minimum_gap + amplitude
    for height_index in range(discretization_n + 1):
        psi = math.pi * height_index / discretization_n
        distance = h0 - amplitude * math.cos(psi)
        j_tables[height_index] = _integral_table_for_shift(
            p_nodes,
            nodes,
            weights,
            matrix_size,
            a,
            c,
            distance,
            0.0,
        ).real
        for image_index, image_shift in enumerate(image_shifts):
            image_tables[height_index, image_index] = _integral_table_for_shift(
                p_nodes,
                nodes,
                weights,
                matrix_size,
                a,
                c,
                distance,
                image_shift,
            ).real

    return DirectGeometry(
        a=a,
        c=c,
        harmonic_weights=harmonic_weights,
        j_tables=j_tables,
        image_tables=image_tables,
        row_scale=row_scale,
        col_scale=col_scale,
        const_c=const_c,
        const_pz=const_pz,
    )


def _integral_table_for_shift(
    p_nodes: np.ndarray,
    nodes: np.ndarray,
    weights: np.ndarray,
    matrix_size: int,
    a: float,
    c: float,
    distance: float,
    image_shift: float,
) -> np.ndarray:
    shifted_position = 2.0 * a / c + 2.0 * distance / c + image_shift / c - a / c * nodes
    base = 1.0 + shifted_position * shifted_position + (a * a / (c * c) - 1.0) * (1.0 - nodes * nodes)
    radical = np.maximum(base * base - 4.0 * shifted_position * shifted_position, 0.0)
    eta = np.sqrt((base + np.sqrt(radical)) / 2.0)
    transformed_argument = shifted_position / eta

    p_argument = _legendre_p_table(matrix_size, transformed_argument)[1:]
    q_radial = _legendre_q_table(matrix_size, eta)[1:]
    return (p_nodes * weights[np.newaxis, :]) @ (p_argument * q_radial).T


def _single_layer_image_charges(
    epsilon_air: complex,
    epsilon_layer: complex,
    epsilon_substrate: complex,
    thickness: float,
    image_charge_count: int,
) -> np.ndarray:
    beta1 = (epsilon_layer - epsilon_air) / (epsilon_layer + epsilon_air)
    beta2 = (epsilon_substrate - epsilon_layer) / (epsilon_substrate + epsilon_layer)
    distances = 2.0 * thickness * np.arange(1, image_charge_count + 1, dtype=float)
    wave_numbers = np.linspace(0.0, 10.0 / thickness, 1000)
    delta_k = wave_numbers[1] - wave_numbers[0]
    reflection = beta1 - (beta1 + beta2 * np.exp(-2.0 * wave_numbers * thickness)) / (
        1.0 + beta1 * beta2 * np.exp(-2.0 * wave_numbers * thickness)
    )
    return _fit_image_charges(distances, delta_k, reflection)


def _double_layer_image_charges(
    epsilon_air: complex,
    epsilon_layer1: complex,
    epsilon_layer2: complex,
    epsilon_substrate: complex,
    thickness1: float,
    thickness2: float,
    image_charge_count: int,
) -> np.ndarray:
    beta1 = (epsilon_layer1 - epsilon_air) / (epsilon_layer1 + epsilon_air)
    beta2 = (epsilon_layer2 - epsilon_layer1) / (epsilon_layer2 + epsilon_layer1)
    beta3 = (epsilon_substrate - epsilon_layer2) / (epsilon_substrate + epsilon_layer2)
    min_thickness = min(thickness1, thickness2)
    distances = np.array(
        [2.0 * thickness1]
        + [
            2.0 * thickness1 + 2.0 * min_thickness * charge_index
            for charge_index in range(1, image_charge_count)
        ],
        dtype=float,
    )
    wave_numbers = np.linspace(0.0, 10.0 / min_thickness, 1000)
    delta_k = wave_numbers[1] - wave_numbers[0]
    r1 = -(
        beta2 + beta3 * np.exp(-2.0 * wave_numbers * thickness2)
    ) / (1.0 + beta2 * beta3 * np.exp(-2.0 * wave_numbers * thickness2))
    reflection = beta1 - (beta1 - r1 * np.exp(-2.0 * wave_numbers * thickness1)) / (
        1.0 - beta1 * r1 * np.exp(-2.0 * wave_numbers * thickness1)
    )
    return _fit_image_charges(distances, delta_k, reflection)


def _fit_image_charges(distances: np.ndarray, delta_k: float, reflection: np.ndarray) -> np.ndarray:
    powers = np.arange(reflection.size, dtype=float)[:, np.newaxis]
    basis = np.exp(-distances * delta_k)[np.newaxis, :] ** powers
    charges, *_ = np.linalg.lstsq(basis, reflection, rcond=None)
    return charges


def _single_layer_far_field_multiplier(
    *,
    epsilon_layer: complex,
    epsilon_substrate: complex,
    thickness: float,
    frequency: float,
    enabled: bool,
    coefficient: float,
    angle_degrees: float,
) -> complex:
    if not enabled:
        return 1.0 + 0.0j

    phi = math.radians(angle_degrees)
    kz1 = _safe_cosine(phi)
    kz2 = np.sqrt(epsilon_layer - math.sin(phi) ** 2)
    kz3 = np.sqrt(epsilon_substrate - math.sin(phi) ** 2)
    wave_number = 2.0 * math.pi * frequency / 1e7
    rp1 = (epsilon_layer / kz2 - 1.0 / kz1) / (epsilon_layer / kz2 + 1.0 / kz1)
    rp2 = (epsilon_substrate / kz3 - epsilon_layer / kz2) / (
        epsilon_substrate / kz3 + epsilon_layer / kz2
    )
    phase = np.exp(2.0j * kz2 * thickness * wave_number)
    reflection = (rp1 + rp2 * phase) / (1.0 + rp1 * rp2 * phase)
    return (1.0 + coefficient * reflection) ** 2


def _double_layer_far_field_multiplier(
    *,
    epsilon_layer1: complex,
    epsilon_layer2: complex,
    epsilon_substrate: complex,
    thickness1: float,
    thickness2: float,
    frequency: float,
    enabled: bool,
    coefficient: float,
    angle_degrees: float,
) -> complex:
    if not enabled:
        return 1.0 + 0.0j

    phi = math.radians(angle_degrees)
    kz1 = _safe_cosine(phi)
    kz2 = np.sqrt(epsilon_layer1 - math.sin(phi) ** 2)
    kz3 = np.sqrt(epsilon_layer2 - math.sin(phi) ** 2)
    kz4 = np.sqrt(epsilon_substrate - math.sin(phi) ** 2)
    wave_number = 2.0 * math.pi * frequency / 1e7
    rp1 = (epsilon_layer1 / kz2 - 1.0 / kz1) / (epsilon_layer1 / kz2 + 1.0 / kz1)
    rp2 = (epsilon_layer2 / kz3 - epsilon_layer1 / kz2) / (
        epsilon_layer2 / kz3 + epsilon_layer1 / kz2
    )
    rp3 = (epsilon_substrate / kz4 - epsilon_layer2 / kz3) / (
        epsilon_substrate / kz4 + epsilon_layer2 / kz3
    )
    phase2 = np.exp(2.0j * kz3 * thickness2 * wave_number)
    nested_reflection = (rp2 + rp3 * phase2) / (1.0 + rp2 * rp3 * phase2)
    phase1 = np.exp(2.0j * kz2 * thickness1 * wave_number)
    reflection = (rp1 + nested_reflection * phase1) / (
        1.0 + rp1 * nested_reflection * phase1
    )
    return (1.0 + coefficient * reflection) ** 2


def _safe_cosine(angle_radians: float) -> float:
    cosine = math.cos(angle_radians)
    if abs(cosine) < 1e-12:
        return 1e-12
    return cosine


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
    if model_type not in {"bulk", "layered"}:
        raise ValidationError("Model type must be either Bulk sample or Layered sample.")

    algorithm = payload.get("algorithm")
    if algorithm not in {"fast", "accurate"}:
        raise ValidationError("Calculation algorithm must be either Fast or Accurate.")

    tip = _require_mapping(payload, "tip")
    far_field = _require_mapping(payload, "farField")
    reference = _require_mapping(payload, "reference")
    frequency_range = _require_mapping(payload, "frequencyRange")
    sample = _require_mapping(payload, "sample") if model_type == "bulk" else {}
    layered = _require_mapping(payload, "layered") if model_type == "layered" else None
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

    sample_input_method = None
    sample_material_id = None
    sample_material_label = None
    uploaded_permittivity = None
    drude_lorentz_parameters = None
    if model_type == "bulk":
        sample_input_method = sample.get("inputMethod")
        if sample_input_method not in {"builtIn", "upload", "drudeLorentz"}:
            raise ValidationError("Sample material input method is invalid.")

        sample_material_id = sample.get("builtInMaterial") if sample_input_method == "builtIn" else None
        sample_material_label = "Uploaded file"
        if sample_input_method == "builtIn":
            if not isinstance(sample_material_id, str) or not sample_material_id:
                raise ValidationError("Please choose a built-in material.")
            try:
                sample_material_label = get_material_label(sample_material_id)
            except ValueError as error:
                raise ValidationError("Unknown built-in material selected.") from error
        if sample_input_method == "drudeLorentz":
            sample_material_label = "Drude-Lorentz model"

        uploaded_permittivity = sample.get("uploadedPermittivity")
        if sample_input_method == "upload" and not isinstance(uploaded_permittivity, str):
            raise ValidationError("Please upload a permittivity file.")

        if sample_input_method == "drudeLorentz":
            drude_lorentz_parameters = _validate_drude_lorentz_parameters(
                _require_mapping(sample, "drudeLorentz")
            )

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

    layered_config = _validate_layered_config(layered) if layered is not None else None

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
        "drude_lorentz_parameters": drude_lorentz_parameters,
        "layered": layered_config,
        "experimental_enabled": experimental_enabled,
        "experimental_content": experimental_content,
        "discretization_n": discretization_n,
        "discretization_m": discretization_m,
    }


def _validate_layered_config(layered: dict[str, object]) -> dict[str, object]:
    structure = layered.get("structure")
    if structure not in {"single", "double"}:
        raise ValidationError("Layered structure must be either 1 layer or 2 layers on substrate.")

    image_charge_count = _require_positive_integer(layered, "imageCharges")
    if image_charge_count > MAX_IMAGE_CHARGES:
        raise ValidationError(f"K must be between 1 and {MAX_IMAGE_CHARGES}.")

    layer1 = _validate_layer_payload(_require_mapping(layered, "layer1"), "Layer 1", "layer1")
    layer2 = None
    if structure == "double":
        layer2 = _validate_layer_payload(_require_mapping(layered, "layer2"), "Layer 2", "layer2")
    substrate = _validate_material_input(
        _require_mapping(layered, "substrate"),
        label="Substrate",
        role="substrate",
    )

    return {
        "structure": structure,
        "image_charge_count": image_charge_count,
        "layer1": layer1,
        "layer2": layer2,
        "substrate": substrate,
    }


def _validate_layer_payload(
    payload: dict[str, object], label: str, role: str
) -> dict[str, object]:
    thickness = _require_positive_number(payload, "thickness", f"{label} thickness must be positive.")
    material = _validate_material_input(
        _require_mapping(payload, "material"),
        label=label,
        role=role,
    )
    return {
        "thickness": thickness,
        "material": material,
    }


def _validate_material_input(
    payload: dict[str, object], *, label: str, role: str
) -> dict[str, object]:
    input_method = payload.get("inputMethod")
    if input_method not in {"builtIn", "upload", "drudeLorentz"}:
        raise ValidationError(f"{label} material input method is invalid.")

    material_id = None
    material_label = f"{label} uploaded file"
    if input_method == "builtIn":
        material_id = payload.get("builtInMaterial")
        if not isinstance(material_id, str) or not material_id:
            raise ValidationError(f"Please choose a built-in material for {label}.")
        try:
            material_label = get_material_label(material_id, role)
        except ValueError as error:
            raise ValidationError(f"Unknown built-in material selected for {label}.") from error

    uploaded_permittivity = payload.get("uploadedPermittivity")
    if input_method == "upload" and not isinstance(uploaded_permittivity, str):
        raise ValidationError(f"Please upload a permittivity file for {label}.")

    drude_lorentz_parameters = None
    if input_method == "drudeLorentz":
        material_label = f"{label} Drude-Lorentz model"
        drude_lorentz_parameters = _validate_drude_lorentz_parameters(
            _require_mapping(payload, "drudeLorentz")
        )

    return {
        "input_method": input_method,
        "material_id": material_id,
        "material_label": material_label,
        "uploaded_permittivity": uploaded_permittivity,
        "drude_lorentz_parameters": drude_lorentz_parameters,
        "role": role,
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

    if input_method == "drudeLorentz":
        return generate_epsilon_drude_lorentz(
            frequency_grid, dict(config["drude_lorentz_parameters"])
        )

    raise ValidationError("Unsupported sample material input method.")


def _resolve_material_permittivity(
    material: dict[str, object], frequency_grid: np.ndarray, label: str
) -> np.ndarray:
    input_method = material["input_method"]

    if input_method == "builtIn":
        return get_builtin_permittivity(
            str(material["material_id"]),
            frequency_grid,
            role=str(material["role"]),
        )

    if input_method == "upload":
        table = parse_numeric_table(
            str(material["uploaded_permittivity"]),
            minimum_columns=3,
            label=f"{label} permittivity file",
        )
        epsilon_values = table[:, 1] + 1j * table[:, 2]
        return interpolate_complex_data(
            table[:, 0],
            epsilon_values,
            frequency_grid,
            range_error_message=f"The selected frequency range is outside the uploaded {label} permittivity data range.",
        )

    if input_method == "drudeLorentz":
        return generate_epsilon_drude_lorentz(
            frequency_grid, dict(material["drude_lorentz_parameters"])
        )

    raise ValidationError(f"Unsupported material input method for {label}.")


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
    if config["model_type"] == "layered":
        layered = dict(config["layered"])
        material_token = "one-layer" if layered["structure"] == "single" else "two-layer"
    else:
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


def _validate_drude_lorentz_parameters(parameters: dict[str, object]) -> dict[str, object]:
    epsilon_infinity = _require_positive_number(
        parameters,
        "epsilonInfinity",
        "High-frequency permittivity must be greater than 0.",
    )

    use_drude = bool(parameters.get("useDrude", False))
    plasma_frequency = 0.0
    drude_damping = 0.0
    if use_drude:
        plasma_frequency = _require_positive_number(
            parameters,
            "plasmaFrequency",
            "Drude plasma frequency must be greater than 0.",
        )
        drude_damping = _require_positive_number(
            parameters,
            "drudeDamping",
            "Drude damping must be greater than 0.",
        )

    oscillators_value = parameters.get("oscillators")
    if not isinstance(oscillators_value, list):
        raise ValidationError("Lorentz oscillators must be provided as a list.")
    if len(oscillators_value) > MAX_DRUDE_LORENTZ_OSCILLATORS:
        raise ValidationError(
            f"At most {MAX_DRUDE_LORENTZ_OSCILLATORS} Lorentz oscillators are allowed."
        )
    oscillators: list[dict[str, float]] = []
    for index, oscillator_value in enumerate(oscillators_value, start=1):
        if not isinstance(oscillator_value, dict):
            raise ValidationError(f"Oscillator {index} must be an object.")
        oscillators.append(
            {
                "strength": _require_positive_number(
                    oscillator_value,
                    "strength",
                    f"Oscillator {index} strength must be greater than 0.",
                ),
                "resonanceFrequency": _require_positive_number(
                    oscillator_value,
                    "resonanceFrequency",
                    f"Oscillator {index} resonance frequency must be greater than 0.",
                ),
                "damping": _require_positive_number(
                    oscillator_value,
                    "damping",
                    f"Oscillator {index} damping must be greater than 0.",
                ),
            }
        )

    return {
        "epsilonInfinity": epsilon_infinity,
        "useDrude": use_drude,
        "plasmaFrequency": plasma_frequency,
        "drudeDamping": drude_damping,
        "oscillators": oscillators,
    }


def _require_number(payload: dict[str, object], key: str, message: str) -> float:
    value = payload.get(key)
    if value is None or isinstance(value, bool):
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
    if isinstance(value, bool):
        raise ValidationError(f"'{key}' must be a positive integer.")
    try:
        integer = int(value)
    except (TypeError, ValueError) as error:
        raise ValidationError(f"'{key}' must be a positive integer.") from error

    if float(value) != integer or integer <= 0:
        raise ValidationError(f"'{key}' must be a positive integer.")
    return integer
