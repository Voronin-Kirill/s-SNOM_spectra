from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MATLAB_APP_ROOT = PROJECT_ROOT / "MATLAB App project"
LAYERED_SIO2_PATH = PROJECT_ROOT / "Eps_SiO2.txt"


@dataclass(frozen=True)
class BuiltInMaterial:
    identifier: str
    label: str
    description: str


BUILT_IN_MATERIALS = (
    BuiltInMaterial(
        identifier="pmma",
        label="PMMA",
        description="Tabulated dielectric permittivity from the original MATLAB app.",
    ),
    BuiltInMaterial(
        identifier="sio2",
        label="SiO2",
        description="Anisotropic Lorentz model from the original MATLAB app.",
    ),
    BuiltInMaterial(
        identifier="sic_4h",
        label="4H-SiC",
        description="Anisotropic Lorentz-Drude model from the original MATLAB app.",
    ),
    BuiltInMaterial(
        identifier="sic_6h",
        label="6H-SiC",
        description="Anisotropic Lorentz-Drude model from the original MATLAB app.",
    ),
)

LAYER_BUILT_IN_MATERIALS = (
    BuiltInMaterial(
        identifier="pmma",
        label="PMMA",
        description="Tabulated dielectric permittivity from the original MATLAB app.",
    ),
    BuiltInMaterial(
        identifier="sio2",
        label="SiO2",
        description="Tabulated isotropic permittivity for layered-sample calculations.",
    ),
    BuiltInMaterial(
        identifier="sic_3c",
        label="3C-SiC",
        description="Lorentz model used by the layered MATLAB functions.",
    ),
)

SUBSTRATE_BUILT_IN_MATERIALS = (
    BuiltInMaterial(
        identifier="si",
        label="Si",
        description="Constant substrate permittivity, epsilon = 12.",
    ),
    BuiltInMaterial(
        identifier="sio2",
        label="SiO2",
        description="Tabulated isotropic permittivity for layered-sample calculations.",
    ),
    BuiltInMaterial(
        identifier="sic_3c",
        label="3C-SiC",
        description="Lorentz model used by the layered MATLAB functions.",
    ),
)


def list_materials() -> list[dict[str, str]]:
    return _serialize_materials(BUILT_IN_MATERIALS)


def list_layer_materials() -> list[dict[str, str]]:
    return _serialize_materials(LAYER_BUILT_IN_MATERIALS)


def list_substrate_materials() -> list[dict[str, str]]:
    return _serialize_materials(SUBSTRATE_BUILT_IN_MATERIALS)


def _serialize_materials(materials: tuple[BuiltInMaterial, ...]) -> list[dict[str, str]]:
    return [
        {
            "id": material.identifier,
            "label": material.label,
            "description": material.description,
        }
        for material in materials
    ]


def get_material_label(material_id: str, role: str = "bulk") -> str:
    for material in _materials_for_role(role):
        if material.identifier == material_id:
            return material.label
    raise ValueError("Unknown built-in material.")


def _materials_for_role(role: str) -> tuple[BuiltInMaterial, ...]:
    if role == "bulk":
        return BUILT_IN_MATERIALS
    if role in {"layer", "layer1", "layer2"}:
        return LAYER_BUILT_IN_MATERIALS
    if role == "substrate":
        return SUBSTRATE_BUILT_IN_MATERIALS
    raise ValueError("Unknown material role.")


def generate_epsilon_drude_lorentz(
    frequency_grid: np.ndarray, parameters: dict[str, object]
) -> np.ndarray:
    omega = np.asarray(frequency_grid, dtype=float)
    epsilon = np.full(
        omega.shape,
        float(parameters["epsilonInfinity"]),
        dtype=np.complex128,
    )

    if parameters.get("useDrude", False):
        plasma_frequency = float(parameters["plasmaFrequency"])
        damping = float(parameters["drudeDamping"])
        epsilon -= plasma_frequency**2 / (omega**2 + 1j * damping * omega)

    for oscillator in parameters.get("oscillators", []):
        oscillator_parameters = dict(oscillator)
        strength = float(oscillator_parameters["strength"])
        resonance = float(oscillator_parameters["resonanceFrequency"])
        damping = float(oscillator_parameters["damping"])
        epsilon += strength * resonance**2 / (
            resonance**2 - omega**2 - 1j * damping * omega
        )

    return epsilon


def get_builtin_permittivity(
    material_id: str, frequency_grid: np.ndarray, role: str = "bulk"
) -> np.ndarray:
    if role in {"layer", "layer1", "layer2", "substrate"}:
        return _get_layered_builtin_permittivity(material_id, frequency_grid, role)

    if material_id == "pmma":
        frequencies, epsilon_values = _load_pmma_dataset()
        return interpolate_complex_data(frequencies, epsilon_values, frequency_grid)

    if material_id == "sio2":
        return _sio2_permittivity(frequency_grid)

    if material_id == "sic_4h":
        return _sic_permittivity(
            frequency_grid,
            eb_par=6.78,
            eb_per=6.56,
            w_to_par=782.0,
            w_to_per=797.0,
            w_lo_par=967.0,
            w_lo_per=971.0,
            gamma_par=3.3,
            gamma_per=3.3,
            wp_par=220.0,
            wp_per=275.0,
            gamma_plasma_par=450.0,
            gamma_plasma_per=450.0,
        )

    if material_id == "sic_6h":
        return _sic_permittivity(
            frequency_grid,
            eb_par=6.72,
            eb_per=6.56,
            w_to_par=788.0,
            w_to_per=797.0,
            w_lo_par=964.0,
            w_lo_per=970.0,
            gamma_par=5.5,
            gamma_per=5.9,
            wp_par=120.0,
            wp_per=230.0,
            gamma_plasma_par=250.0,
            gamma_plasma_per=500.0,
        )

    raise ValueError("Unknown built-in material.")


def _get_layered_builtin_permittivity(
    material_id: str, frequency_grid: np.ndarray, role: str
) -> np.ndarray:
    if material_id == "pmma" and role != "substrate":
        frequencies, epsilon_values = _load_pmma_dataset()
        return interpolate_complex_data(frequencies, epsilon_values, frequency_grid)

    if material_id == "si" and role == "substrate":
        return np.full(np.asarray(frequency_grid, dtype=float).shape, 12.0 + 0.0j, dtype=np.complex128)

    if material_id == "sio2":
        frequencies, epsilon_values = _load_layered_sio2_dataset()
        return interpolate_complex_data(frequencies, epsilon_values, frequency_grid)

    if material_id == "sic_3c":
        return _sic_3c_permittivity(frequency_grid)

    raise ValueError("Unknown built-in material.")


def interpolate_complex_data(
    frequencies: np.ndarray,
    epsilon_values: np.ndarray,
    frequency_grid: np.ndarray,
    range_error_message: str | None = None,
) -> np.ndarray:
    source_frequencies = np.asarray(frequencies, dtype=float)
    source_values = np.asarray(epsilon_values, dtype=np.complex128)
    grid = np.asarray(frequency_grid, dtype=float)

    order = np.argsort(source_frequencies)
    source_frequencies = source_frequencies[order]
    source_values = source_values[order]

    if np.any(np.diff(source_frequencies) <= 0):
        raise ValueError("Frequency values must be strictly increasing.")

    if grid[0] < source_frequencies[0] or grid[-1] > source_frequencies[-1]:
        raise ValueError(
            range_error_message
            or "The selected frequency range is outside the available permittivity data range."
        )

    real_part = np.interp(grid, source_frequencies, source_values.real)
    imaginary_part = np.interp(grid, source_frequencies, source_values.imag)
    return real_part + 1j * imaginary_part


@lru_cache(maxsize=1)
def _load_pmma_dataset() -> tuple[np.ndarray, np.ndarray]:
    real_data = np.loadtxt(MATLAB_APP_ROOT / "real_eps_PMMA.txt")
    imaginary_data = np.loadtxt(MATLAB_APP_ROOT / "imag_eps_PMMA.txt")

    frequencies = real_data[:, 0]
    epsilon_values = real_data[:, 1] + 1j * imaginary_data[:, 1]
    return frequencies, epsilon_values


@lru_cache(maxsize=1)
def _load_layered_sio2_dataset() -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(LAYERED_SIO2_PATH)
    return data[:, 0], data[:, 1] + 1j * data[:, 2]


def _positive_imaginary_sqrt(values: np.ndarray) -> np.ndarray:
    roots = np.sqrt(np.asarray(values, dtype=np.complex128))
    negative_imaginary_branch = roots.imag < 0
    roots[negative_imaginary_branch] *= -1
    return roots


def _sio2_permittivity(frequency_grid: np.ndarray) -> np.ndarray:
    omega = np.asarray(frequency_grid, dtype=float)

    w_lo_par = np.array([1223, 1222, 794, 540, 514, 559, 413], dtype=float)
    w_to_par = np.array([1220, 1080, 778, 539, 509, 495, 364], dtype=float)
    gamma_par = np.array([183, 7.5, 78, 22, 7.1, 4.5, 5.1], dtype=float)

    w_lo_per = np.array([1229, 1165, 1215, 815, 700, 522, 421], dtype=float)
    w_to_per = np.array([1227, 1163, 1072, 797, 697, 450, 394], dtype=float)
    gamma_per = np.array([135, 7.0, 7.6, 7.2, 8.4, 4.0, 2.8], dtype=float)

    epsilon_par = np.ones_like(omega, dtype=np.complex128)
    epsilon_per = np.ones_like(omega, dtype=np.complex128)

    for index in range(w_lo_par.size):
        epsilon_par += (w_lo_par[index] ** 2 - w_to_par[index] ** 2) / (
            w_to_par[index] ** 2 - omega**2 - 1j * gamma_par[index] * omega
        )
        epsilon_per += (w_lo_per[index] ** 2 - w_to_per[index] ** 2) / (
            w_to_per[index] ** 2 - omega**2 - 1j * gamma_per[index] * omega
        )

    epsilon_par *= 2.383
    epsilon_per *= 2.356
    return _positive_imaginary_sqrt(epsilon_par * epsilon_per)


def _sic_3c_permittivity(frequency_grid: np.ndarray) -> np.ndarray:
    omega = np.asarray(frequency_grid, dtype=float)
    epsilon_infinity = 6.6
    w_to = 797.0
    w_lo = 973.0
    damping = 6.0
    return epsilon_infinity * (omega**2 - w_lo**2 + 1j * damping * omega) / (
        omega**2 - w_to**2 + 1j * damping * omega
    )


def _sic_permittivity(
    frequency_grid: np.ndarray,
    *,
    eb_par: float,
    eb_per: float,
    w_to_par: float,
    w_to_per: float,
    w_lo_par: float,
    w_lo_per: float,
    gamma_par: float,
    gamma_per: float,
    wp_par: float,
    wp_per: float,
    gamma_plasma_par: float,
    gamma_plasma_per: float,
) -> np.ndarray:
    omega = np.asarray(frequency_grid, dtype=float)
    epsilon_par = eb_par * (
        1
        + (w_lo_par**2 - w_to_par**2) / (w_to_par**2 - omega**2 - 1j * gamma_par * omega)
        - wp_par**2 / (omega**2 + 1j * omega * gamma_plasma_par)
    )
    epsilon_per = eb_per * (
        1
        + (w_lo_per**2 - w_to_per**2) / (w_to_per**2 - omega**2 - 1j * gamma_per * omega)
        - wp_per**2 / (omega**2 + 1j * omega * gamma_plasma_per)
    )
    return _positive_imaginary_sqrt(epsilon_par * epsilon_per)
