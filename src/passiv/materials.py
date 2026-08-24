"""Thermophysical material definitions."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Material:
    """Constant thermal properties for an isotropic material."""

    k: float       # W/(m K)
    rho: float     # kg/m^3
    cp: float      # J/(kg K)


# Values preserved from the original working materials.py.
MATERIALS = {    
    "fiber_cement": Material(k=1.0, rho=1800.0, cp=840.0),
    "air": Material(k=0.144, rho=1.2, cp=1005.0),
    "mineral_wool_ext": Material(k=0.034, rho=80.0, cp=1400.0),
    "mineral_wool_int": Material(k=0.035, rho=100.0, cp=1400.0),
    "plywood": Material(k=0.13, rho=520.0, cp=1600.0),
    "wood": Material(k=0.12, rho=450.0, cp=1600.0),
    "drywall": Material(k=0.25, rho=800.0, cp=1090.0),
    "concrete": Material(k=2.1,rho=2500.0,cp=950.0,),
    "soil": Material(k=2.0,rho=1500.0,cp=900.0,),
}
