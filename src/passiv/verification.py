"""
Verification utilities for analytical/manufactured heat-transfer cases.
"""

import numpy as np
import ufl


def manufactured_temperature_numpy(
    x,
    *,
    lx,
    ly,
    t0,
    amplitude,
):
    return (
        t0
        + amplitude
        * np.sin(np.pi * x[0] / lx)
        * np.sin(np.pi * x[1] / ly)
    )


def manufactured_temperature_ufl(
    domain,
    *,
    lx,
    ly,
    t0,
    amplitude,
):
    x = ufl.SpatialCoordinate(domain)

    return (
        t0
        + amplitude
        * ufl.sin(np.pi * x[0] / lx)
        * ufl.sin(np.pi * x[1] / ly)
    )


def manufactured_source_ufl(
    domain,
    *,
    lx,
    ly,
    conductivity,
    amplitude,
):
    x = ufl.SpatialCoordinate(domain)

    shape = (
        ufl.sin(np.pi * x[0] / lx)
        * ufl.sin(np.pi * x[1] / ly)
    )

    return (
        conductivity
        * amplitude
        * np.pi**2
        * (
            1.0 / lx**2
            + 1.0 / ly**2
        )
        * shape
    )




    