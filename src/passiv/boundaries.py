"""Thermal boundary-condition descriptions."""

from dataclasses import dataclass


@dataclass(frozen=True)
class ConvectionBoundary:
    """Robin boundary: -k grad(T).n = h (T - T_inf)."""

    tag: int
    h: float
    temperature: float
    name: str = ""
