"""Tests for the finite-volume ideal MHD solver."""

from __future__ import annotations

import numpy as np
import pytest

from mhd_sandbox import Grid, IdealMHD2DSolver, MultiFluidIdealMHD2DSolver


def _make_grid(num_gridpoints: tuple[int, int] = (24, 24)) -> Grid:
    """Build a periodic unit-square grid for MHD tests."""
    return Grid(
        lower_bound=(0.0, 0.0),
        upper_bound=(1.0, 1.0),
        num_gridpoints=num_gridpoints,
    )


def test_primitive_conservative_roundtrip() -> None:
    """Primitive variables survive conversion into conservative form."""
    grid = _make_grid((8, 8))
    solver = IdealMHD2DSolver(grid)

    rho = np.full(grid.num_gridpoints, 1.2)
    pressure = np.full(grid.num_gridpoints, 0.7)
    velocity = np.zeros(grid.num_gridpoints + (3,))
    velocity[..., 0] = 0.3
    velocity[..., 1] = -0.2
    magnetic_field = np.zeros(grid.num_gridpoints + (3,))
    magnetic_field[..., 0] = 0.1
    magnetic_field[..., 1] = 0.4

    solver.set_primitive(
        density=rho,
        velocity=velocity,
        magnetic_field=magnetic_field,
        pressure=pressure,
    )
    primitives = solver.primitive_variables()

    assert np.allclose(primitives["density"], rho)
    assert np.allclose(primitives["pressure"], pressure)
    assert np.allclose(primitives["velocity"], velocity)
    assert np.allclose(primitives["magnetic_field"], magnetic_field)


def test_constant_state_is_preserved() -> None:
    """A uniform periodic MHD state has zero finite-volume update."""
    grid = _make_grid((16, 12))
    solver = IdealMHD2DSolver(grid)

    rho = np.full(grid.num_gridpoints, 1.0)
    pressure = np.full(grid.num_gridpoints, 1.0)
    velocity = np.zeros(grid.num_gridpoints + (3,))
    velocity[..., 0] = 0.25
    velocity[..., 1] = -0.1
    magnetic_field = np.zeros(grid.num_gridpoints + (3,))
    magnetic_field[..., 0] = 0.5
    magnetic_field[..., 1] = 0.2

    solver.set_primitive(
        density=rho,
        velocity=velocity,
        magnetic_field=magnetic_field,
        pressure=pressure,
    )
    state0 = solver.state.copy()
    solver.step(dt=0.001)

    assert np.allclose(solver.state, state0)


def test_orszag_tang_initial_state_is_divergence_free() -> None:
    """The analytic Orszag-Tang magnetic field starts with div(B)=0."""
    grid = _make_grid((32, 32))
    solver = IdealMHD2DSolver(grid)
    solver.set_state(IdealMHD2DSolver.orszag_tang_initial_state(grid))

    assert solver.max_abs_divergence_b() < 1e-12


def test_short_orszag_tang_run_stays_physical() -> None:
    """A short Orszag-Tang run keeps density and pressure above floors."""
    grid = _make_grid((24, 24))
    solver = IdealMHD2DSolver(grid, cfl=0.25)
    solver.set_state(IdealMHD2DSolver.orszag_tang_initial_state(grid))
    result = solver.solve(t_end=0.01)
    primitives = solver.primitive_variables()

    assert result.steps > 0
    assert result.time == 0.01
    assert result.max_divergence_b < 1e-10
    assert float(np.min(primitives["density"])) > solver.density_floor
    assert float(np.min(primitives["pressure"])) > solver.pressure_floor


def test_multifluid_solver_uses_most_restrictive_timestep() -> None:
    """Multi-fluid timestep follows the smallest per-fluid CFL timestep."""
    grid = _make_grid((12, 12))
    solver = MultiFluidIdealMHD2DSolver(
        grid=grid,
        fluid_names=("electrons", "ions"),
        gammas={"electrons": 1.4, "ions": 5.0 / 3.0},
        cfl=0.3,
        project_divergence=False,
    )

    density = np.full(grid.num_gridpoints, 1.0)
    velocity = np.zeros(grid.num_gridpoints + (3,))
    magnetic_field = np.zeros(grid.num_gridpoints + (3,))
    magnetic_field[..., 0] = 0.2
    magnetic_field[..., 1] = 0.1

    solver.set_primitive(
        "electrons",
        density=density,
        velocity=velocity,
        magnetic_field=magnetic_field,
        pressure=np.full(grid.num_gridpoints, 1.5),
    )
    solver.set_primitive(
        "ions",
        density=density,
        velocity=velocity,
        magnetic_field=magnetic_field,
        pressure=np.full(grid.num_gridpoints, 0.2),
    )

    electrons_dt = solver.fluid_solver("electrons").compute_dt()
    ions_dt = solver.fluid_solver("ions").compute_dt()
    expected_dt = min(electrons_dt, ions_dt)

    assert solver.compute_dt() == pytest.approx(expected_dt)
    used_dt = solver.step()
    assert used_dt == pytest.approx(expected_dt)
    assert solver.time == pytest.approx(expected_dt)
    assert solver.steps == 1
    assert solver.fluid_solver("electrons").time == pytest.approx(expected_dt)
    assert solver.fluid_solver("ions").time == pytest.approx(expected_dt)


def test_multifluid_solver_advances_all_fluids_to_target_time() -> None:
    """Multi-fluid solve() keeps all fluid solvers synchronized."""
    grid = _make_grid((8, 8))
    solver = MultiFluidIdealMHD2DSolver(
        grid=grid,
        fluid_names=("e", "i"),
        cfl=0.25,
        project_divergence=False,
    )

    density = np.full(grid.num_gridpoints, 1.0)
    velocity = np.zeros(grid.num_gridpoints + (3,))
    magnetic_field = np.zeros(grid.num_gridpoints + (3,))
    pressure = np.full(grid.num_gridpoints, 0.5)

    solver.set_primitive("e", density, velocity, magnetic_field, pressure)
    solver.set_primitive("i", density, velocity, magnetic_field, pressure)
    result = solver.solve(t_end=0.01, dt=0.002)

    assert result.time == pytest.approx(0.01)
    assert result.steps == 5
    assert set(result.states.keys()) == {"e", "i"}
    assert set(result.max_divergence_b.keys()) == {"e", "i"}
    assert all(
        state.shape == grid.num_gridpoints + (8,) for state in result.states.values()
    )
    assert solver.fluid_solver("e").time == pytest.approx(0.01)
    assert solver.fluid_solver("i").time == pytest.approx(0.01)
