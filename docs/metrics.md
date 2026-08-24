# Metrics, post-processing, and visualization

## Sign convention

`metrics.boundary_heat_flow()` integrates the conductive heat flow **outward** from the computational domain:

```text
q_out = -k grad(T) . n
```

For a steady source-free calculation, signed outward heat flows should sum to zero.

For convection boundaries, `convection_boundary_heat_flow()` evaluates the equivalent Robin flux directly:

```text
q_out = h (T - T_inf)
```

This is the preferred global energy-conservation diagnostic for the current P1 temperature formulation because the raw elementwise gradient is not an H(div)-conforming reconstructed flux.

## Global heat-flow utilities

### `boundary_heat_flow(...)`

Integrates the gradient-based conductive flux over a tagged boundary.

### `convection_boundary_heat_flow(...)`

Integrates the applied Robin heat flux over a tagged convection boundary.

### `boundary_measure(...)`

Returns tagged boundary length for a 2-D model.

### `energy_balance_error(*heat_flows)`

Returns

```text
abs(sum(Q_i)) / sum(abs(Q_i))
```

for signed outward heat flows.

## U-value and resistance utilities

`metrics.py` includes:

- `series_thermal_resistance()`
- `series_u_value()`
- `multilayer_thermal_resistance()`
- `multilayer_u_value()`
- `parallel_path_u_value()`
- `effective_u_value()`

These are useful for comparing FEM results with one-dimensional engineering reference models.

For a simple planar assembly, `compute_assembly_metrics()` combines exterior/interior Robin heat flows, energy balance, boundary length, temperature difference, and effective U-value into an `AssemblyMetrics` object. It expects exterior and interior reference boundary lengths to be equal and will reject geometries where that assumption is false.

For junctions such as a wall/floor/foundation detail, use explicit heat-flow and coupling/psi metrics rather than forcing the result into a planar U-value.

## 2-D thermal-bridge quantities

### Thermal coupling coefficient

```python
L2D = thermal_coupling_coefficient(Q, delta_T)
```

For a 2-D cross-section, this has units W/(m K).

### Linear thermal transmittance

```python
psi = linear_thermal_transmittance(
    heat_flow,
    delta_temperature,
    reference_u_lengths,
)
```

The implementation uses:

```text
psi = |Q|/|delta_T| - sum(U_i L_i)
```

where the reference `U_i L_i` terms must be chosen consistently with the intended building-physics convention.

## Analytic layered-wall tools

`postprocess.py` provides:

- `layered_wall_resistance()`
- `layered_wall_heat_flux()`
- `analytic_layered_temperature()`

These are used by the steady layered-wall verification test and are useful for regression checks on nominally 1-D assemblies.

## Point and profile sampling

`evaluate_function_at_points()` evaluates a scalar DOLFINx function at supplied physical points. It is intended primarily for serial post-processing.

`sample_horizontal_profile()` samples temperature through a domain along a horizontal line, which is useful for comparing a layered FEM solution with a piecewise-linear analytic profile.

## Visualization

`visualization.py` includes:

- `plot_temperature_field()`
- `plot_material_regions()`
- `plot_temperature_profile()`
- `plot_mesh()`
- `plot_multiple_temperature_profiles()`

PyVista is used for 2-D mesh/field views; Matplotlib is used for profiles.

On macOS, visualization has previously been more fragile than the FEM solve itself. A useful debugging pattern is to complete all solver and metric calculations before opening interactive PyVista windows. The wall/floor/foundation example follows this pattern.
