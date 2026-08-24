"""Semantic geometry-region descriptions."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Region:
    """A tagged material region in a mesh."""

    name: str
    material: str
    tag: int


@dataclass(frozen=True)
class Layer(Region):
    """A 1-D-through-thickness region used by a layered-wall geometry."""

    thickness: float
    n_through: int = 2
