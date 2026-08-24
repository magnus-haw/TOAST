import numpy as np

from toast.materials import Material
from toast.postprocess import (
    analytic_layered_temperature,
    layered_wall_heat_flux,
    layered_wall_resistance,
)
from toast.regions import Layer


def test_layered_analytic_endpoints_and_resistance():
    materials = {
        "a": Material(k=0.10, rho=1.0, cp=1.0),
        "b": Material(k=0.20, rho=1.0, cp=1.0),
    }
    layers = (
        Layer("layer_a", "a", 1, 0.02, 2),
        Layer("layer_b", "b", 2, 0.04, 2),
    )
    h_ext = 10.0
    h_int = 5.0
    T_ext = 273.15
    T_int = 293.15

    expected_R = 1 / h_ext + 0.02 / 0.10 + 0.04 / 0.20 + 1 / h_int
    assert np.isclose(
        layered_wall_resistance(layers, materials, h_ext, h_int), expected_R
    )

    q = layered_wall_heat_flux(
        layers, materials, h_ext, h_int, T_ext, T_int
    )
    assert np.isclose(q, (T_int - T_ext) / expected_R)

    x = np.array([0.0, sum(layer.thickness for layer in layers)])
    T = analytic_layered_temperature(
        x, layers, materials, h_ext, h_int, T_ext, T_int
    )
    assert np.isclose(T[0], T_ext + q / h_ext)
    assert np.isclose(T[-1], T_int - q / h_int)
