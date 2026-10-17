"""Finite-volume solvers for the two-dimensional ideal MHD equations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mhd_sandbox.grid import Grid


@dataclass
class MHDResult:
    """Result returned by an MHD time integration."""

    state: np.ndarray
    time: float
    steps: int
    dt_history: list[float]
    max_divergence_b: float


@dataclass
class MultiFluidMHDResult:
    """Result returned by a multi-fluid MHD integration."""

    states: dict[str, np.ndarray]
    time: float
    steps: int
    dt_history: list[float]
    max_divergence_b: dict[str, float]


class IdealMHD2DSolver:
    """Explicit finite-volume solver for the two-dimensional ideal MHD equations.

    The conservative state is stored as
    ``[rho, rho*vx, rho*vy, rho*vz, E, Bx, By, Bz]`` on a periodic,
    cell-centered grid.  Fluxes are computed with the local Lax-Friedrichs
    (Rusanov) numerical flux and the time step is forward Euler, both first
    order, with an optional periodic magnetic-field projection for the
    div(B)=0 constraint: the simplest scheme that works, as a starting point
    for higher-order reconstructions, better Riemann solvers and other
    div(B) treatments.
    """

    NUM_COMPONENTS = 8
    RHO = 0
    MOM_X = 1
    MOM_Y = 2
    MOM_Z = 3
    ENERGY = 4
    B_X = 5
    B_Y = 6
    B_Z = 7

    def __init__(
        self,
        grid: Grid,
        gamma: float = 5.0 / 3.0,
        cfl: float = 0.35,
        density_floor: float = 1e-8,
        pressure_floor: float = 1e-8,
        project_divergence: bool = True,
    ) -> None:
        """Initialize the MHD solver."""
        if grid.ndim != 2:
            raise ValueError("IdealMHD2DSolver requires a two-dimensional grid")
        if gamma <= 1.0:
            raise ValueError("gamma must be greater than one")
        if cfl <= 0.0:
            raise ValueError("cfl must be positive")

        self._grid = grid
        self._gamma = gamma
        self._cfl = cfl
        self._density_floor = density_floor
        self._pressure_floor = pressure_floor
        self._project_divergence = project_divergence

        nx, ny = grid.num_gridpoints
        lx = grid.upper_bound[0] - grid.lower_bound[0]
        ly = grid.upper_bound[1] - grid.lower_bound[1]
        self._dx = lx / nx
        self._dy = ly / ny

        self._state = np.zeros((nx, ny, self.NUM_COMPONENTS), dtype=float)
        self._time = 0.0
        self._steps = 0
        self._dt_history: list[float] = []

    @staticmethod
    def _as_three_components(vector: np.ndarray, name: str) -> np.ndarray:
        """Return a vector field with three components."""
        arr = np.asarray(vector, dtype=float)
        if arr.shape[-1] == 3:
            return arr
        if arr.shape[-1] == 2:
            zeros = np.zeros(arr.shape[:-1] + (1,), dtype=arr.dtype)
            return np.concatenate((arr, zeros), axis=-1)
        raise ValueError(f"{name} must have two or three components")

    @classmethod
    def primitive_to_conservative(
        cls,
        density: np.ndarray,
        velocity: np.ndarray,
        magnetic_field: np.ndarray,
        pressure: np.ndarray,
        gamma: float = 5.0 / 3.0,
    ) -> np.ndarray:
        """Convert primitive variables to conservative MHD state variables."""
        rho = np.asarray(density, dtype=float)
        pressure = np.asarray(pressure, dtype=float)
        velocity = cls._as_three_components(velocity, "velocity")
        magnetic_field = cls._as_three_components(magnetic_field, "magnetic_field")

        if velocity.shape[:-1] != rho.shape:
            raise ValueError("velocity shape must match density shape plus components")
        if magnetic_field.shape[:-1] != rho.shape:
            raise ValueError(
                "magnetic_field shape must match density shape plus components",
            )
        if pressure.shape != rho.shape:
            raise ValueError("pressure shape must match density shape")

        state = np.zeros(rho.shape + (cls.NUM_COMPONENTS,), dtype=float)
        state[..., cls.RHO] = rho
        state[..., cls.MOM_X : cls.MOM_Z + 1] = rho[..., None] * velocity
        state[..., cls.B_X : cls.B_Z + 1] = magnetic_field

        kinetic_energy = 0.5 * rho * np.sum(velocity * velocity, axis=-1)
        magnetic_energy = 0.5 * np.sum(magnetic_field * magnetic_field, axis=-1)
        state[..., cls.ENERGY] = pressure / (gamma - 1.0)
        state[..., cls.ENERGY] += kinetic_energy + magnetic_energy
        return state

    @classmethod
    def orszag_tang_initial_state(
        cls,
        grid: Grid,
        gamma: float = 5.0 / 3.0,
    ) -> np.ndarray:
        """Build the standard Orszag-Tang vortex initial condition.

        The domain is interpreted as periodic and cell-centered.  The returned
        state uses normalized MHD units with magnetic permeability set to one.
        """
        if grid.ndim != 2:
            raise ValueError("Orszag-Tang initial state requires a 2-D grid")

        nx, ny = grid.num_gridpoints
        x0, y0 = grid.lower_bound
        x1, y1 = grid.upper_bound
        dx = (x1 - x0) / nx
        dy = (y1 - y0) / ny
        x = x0 + (np.arange(nx) + 0.5) * dx
        y = y0 + (np.arange(ny) + 0.5) * dy
        xx, yy = np.meshgrid(x, y, indexing="ij")

        density = np.full((nx, ny), 25.0 / (36.0 * np.pi), dtype=float)
        pressure = np.full((nx, ny), 5.0 / (12.0 * np.pi), dtype=float)

        velocity = np.zeros((nx, ny, 3), dtype=float)
        velocity[..., 0] = -np.sin(2.0 * np.pi * yy)
        velocity[..., 1] = np.sin(2.0 * np.pi * xx)

        magnetic_field = np.zeros((nx, ny, 3), dtype=float)
        magnetic_scale = np.sqrt(4.0 * np.pi)
        magnetic_field[..., 0] = -np.sin(2.0 * np.pi * yy) / magnetic_scale
        magnetic_field[..., 1] = np.sin(4.0 * np.pi * xx) / magnetic_scale

        return cls.primitive_to_conservative(
            density=density,
            velocity=velocity,
            magnetic_field=magnetic_field,
            pressure=pressure,
            gamma=gamma,
        )

    def _raw_pressure(self, state: np.ndarray) -> np.ndarray:
        """Return pressure without applying the pressure floor."""
        rho = np.maximum(state[..., self.RHO], self.density_floor)
        momentum = state[..., self.MOM_X : self.MOM_Z + 1]
        velocity = momentum / rho[..., None]
        magnetic_field = state[..., self.B_X : self.B_Z + 1]
        kinetic_energy = 0.5 * rho * np.sum(velocity * velocity, axis=-1)
        magnetic_energy = 0.5 * np.sum(magnetic_field * magnetic_field, axis=-1)
        return (self.gamma - 1.0) * (
            state[..., self.ENERGY] - kinetic_energy - magnetic_energy
        )

    def _enforce_floors(self) -> None:
        """Apply density and pressure floors in-place."""
        self._state[..., self.RHO] = np.maximum(
            self._state[..., self.RHO],
            self.density_floor,
        )

        primitives = self.primitive_variables(self._state)
        rho = primitives["density"]
        velocity = primitives["velocity"]
        magnetic_field = primitives["magnetic_field"]
        raw_pressure = self._raw_pressure(self._state)

        kinetic_energy = 0.5 * rho * np.sum(velocity * velocity, axis=-1)
        magnetic_energy = 0.5 * np.sum(magnetic_field * magnetic_field, axis=-1)
        floored_energy = (
            self.pressure_floor / (self.gamma - 1.0) + kinetic_energy + magnetic_energy
        )
        self._state[..., self.ENERGY] = np.where(
            raw_pressure < self.pressure_floor,
            floored_energy,
            self._state[..., self.ENERGY],
        )

    def set_state(self, state: np.ndarray) -> None:
        """Set the conservative MHD state."""
        expected_shape = self.grid.num_gridpoints + (self.NUM_COMPONENTS,)
        if state.shape != expected_shape:
            raise ValueError(
                f"Expected state shape {expected_shape}, got {state.shape}",
            )
        self._state = np.asarray(state, dtype=float).copy()
        self._enforce_floors()
        if self.project_divergence:
            self.project_magnetic_field()

    def set_primitive(
        self,
        density: np.ndarray,
        velocity: np.ndarray,
        magnetic_field: np.ndarray,
        pressure: np.ndarray,
    ) -> None:
        """Set the state from primitive MHD variables."""
        self.set_state(
            self.primitive_to_conservative(
                density=density,
                velocity=velocity,
                magnetic_field=magnetic_field,
                pressure=pressure,
                gamma=self.gamma,
            ),
        )

    def primitive_variables(
        self, state: np.ndarray | None = None
    ) -> dict[str, np.ndarray]:
        """Return primitive variables for a conservative state."""
        if state is None:
            state = self.state
        return self.primitive_variables_from_state(
            state=state,
            gamma=self.gamma,
            density_floor=self.density_floor,
            pressure_floor=self.pressure_floor,
        )

    @classmethod
    def primitive_variables_from_state(
        cls,
        state: np.ndarray,
        gamma: float = 5.0 / 3.0,
        density_floor: float = 1e-8,
        pressure_floor: float = 1e-8,
    ) -> dict[str, np.ndarray]:
        """Return primitive variables from a conservative MHD state."""
        rho = np.maximum(state[..., cls.RHO], density_floor)
        momentum = state[..., cls.MOM_X : cls.MOM_Z + 1]
        velocity = momentum / rho[..., None]
        magnetic_field = state[..., cls.B_X : cls.B_Z + 1]

        kinetic_energy = 0.5 * rho * np.sum(velocity * velocity, axis=-1)
        magnetic_energy = 0.5 * np.sum(magnetic_field * magnetic_field, axis=-1)
        pressure = (gamma - 1.0) * (
            state[..., cls.ENERGY] - kinetic_energy - magnetic_energy
        )
        pressure = np.maximum(pressure, pressure_floor)

        return {
            "density": rho,
            "velocity": velocity,
            "magnetic_field": magnetic_field,
            "pressure": pressure,
        }

    def _flux(self, state: np.ndarray, axis: int) -> np.ndarray:
        """Return physical MHD fluxes along x (axis=0) or y (axis=1)."""
        primitives = self.primitive_variables(state)
        rho = primitives["density"]
        velocity = primitives["velocity"]
        magnetic_field = primitives["magnetic_field"]
        pressure = primitives["pressure"]

        vx = velocity[..., 0]
        vy = velocity[..., 1]
        vz = velocity[..., 2]
        bx = magnetic_field[..., 0]
        by = magnetic_field[..., 1]
        bz = magnetic_field[..., 2]

        magnetic_energy = 0.5 * (bx * bx + by * by + bz * bz)
        total_pressure = pressure + magnetic_energy
        v_dot_b = vx * bx + vy * by + vz * bz
        energy = state[..., self.ENERGY]

        flux = np.zeros_like(state)
        if axis == 0:
            flux[..., self.RHO] = rho * vx
            flux[..., self.MOM_X] = rho * vx * vx + total_pressure - bx * bx
            flux[..., self.MOM_Y] = rho * vx * vy - bx * by
            flux[..., self.MOM_Z] = rho * vx * vz - bx * bz
            flux[..., self.ENERGY] = (energy + total_pressure) * vx - bx * v_dot_b
            flux[..., self.B_X] = 0.0
            flux[..., self.B_Y] = vx * by - vy * bx
            flux[..., self.B_Z] = vx * bz - vz * bx
        elif axis == 1:
            flux[..., self.RHO] = rho * vy
            flux[..., self.MOM_X] = rho * vy * vx - by * bx
            flux[..., self.MOM_Y] = rho * vy * vy + total_pressure - by * by
            flux[..., self.MOM_Z] = rho * vy * vz - by * bz
            flux[..., self.ENERGY] = (energy + total_pressure) * vy - by * v_dot_b
            flux[..., self.B_X] = vy * bx - vx * by
            flux[..., self.B_Y] = 0.0
            flux[..., self.B_Z] = vy * bz - vz * by
        else:
            raise ValueError("axis must be 0 or 1")

        return flux

    def _fast_magnetosonic_speed(self, state: np.ndarray, axis: int) -> np.ndarray:
        """Return the fast magnetosonic speed along one coordinate axis."""
        primitives = self.primitive_variables(state)
        rho = primitives["density"]
        pressure = primitives["pressure"]
        magnetic_field = primitives["magnetic_field"]

        bx = magnetic_field[..., 0]
        by = magnetic_field[..., 1]
        bz = magnetic_field[..., 2]
        b_normal = bx if axis == 0 else by

        sound_speed_sq = self.gamma * pressure / rho
        alfven_speed_sq = (bx * bx + by * by + bz * bz) / rho
        normal_alfven_speed_sq = b_normal * b_normal / rho
        sum_speeds = sound_speed_sq + alfven_speed_sq
        discriminant = sum_speeds * sum_speeds
        discriminant -= 4.0 * sound_speed_sq * normal_alfven_speed_sq
        discriminant = np.maximum(discriminant, 0.0)
        fast_sq = 0.5 * (sum_speeds + np.sqrt(discriminant))
        return np.sqrt(np.maximum(fast_sq, 0.0))

    def _characteristic_speed(self, state: np.ndarray, axis: int) -> np.ndarray:
        """Return the local maximum wave speed along one axis."""
        primitives = self.primitive_variables(state)
        velocity = primitives["velocity"]
        normal_velocity = velocity[..., axis]
        return np.abs(normal_velocity) + self._fast_magnetosonic_speed(state, axis)

    def _rusanov_flux(
        self, left_state: np.ndarray, right_state: np.ndarray, axis: int
    ) -> np.ndarray:
        """Return the Rusanov numerical flux at cell interfaces."""
        left_flux = self._flux(left_state, axis=axis)
        right_flux = self._flux(right_state, axis=axis)
        speed = np.maximum(
            self._characteristic_speed(left_state, axis=axis),
            self._characteristic_speed(right_state, axis=axis),
        )
        return 0.5 * (left_flux + right_flux) - 0.5 * speed[..., None] * (
            right_state - left_state
        )

    def compute_dt(self) -> float:
        """Return a CFL-limited timestep for the current state."""
        speed_x = self._characteristic_speed(self.state, axis=0)
        speed_y = self._characteristic_speed(self.state, axis=1)
        max_speed = max(float(np.max(speed_x)), float(np.max(speed_y)))
        if max_speed <= 0.0:
            return self.cfl * min(self.dx, self.dy)
        return self.cfl * min(self.dx, self.dy) / max_speed

    def step(self, dt: float | None = None) -> float:
        """Advance the MHD state by one timestep and return the timestep used."""
        if dt is None:
            dt = self.compute_dt()
        if dt <= 0.0:
            raise ValueError("dt must be positive")

        state = self.state
        right_x = np.roll(state, shift=-1, axis=0)
        right_y = np.roll(state, shift=-1, axis=1)

        flux_x = self._rusanov_flux(state, right_x, axis=0)
        flux_y = self._rusanov_flux(state, right_y, axis=1)

        divergence_flux_x = flux_x - np.roll(flux_x, shift=1, axis=0)
        divergence_flux_y = flux_y - np.roll(flux_y, shift=1, axis=1)

        self._state = state - (dt / self.dx) * divergence_flux_x
        self._state -= (dt / self.dy) * divergence_flux_y
        self._enforce_floors()
        if self.project_divergence:
            self.project_magnetic_field()

        self._time += dt
        self._steps += 1
        self._dt_history.append(float(dt))
        return float(dt)

    def solve(
        self,
        t_end: float,
        max_steps: int = 100_000,
        dt: float | None = None,
        verbose: bool = False,
    ) -> MHDResult:
        """Advance the state until ``t_end`` and return the final state.

        If ``dt`` is omitted, each step uses a CFL-limited timestep.  When
        ``dt`` is provided, it is clipped only on the last step so the solver
        lands exactly on ``t_end``.
        """
        if t_end < self.time:
            raise ValueError("t_end must be greater than or equal to the current time")
        if max_steps <= 0:
            raise ValueError("max_steps must be positive")

        stop_step = self.steps + max_steps
        while self.time < t_end:
            if self.steps >= stop_step:
                raise RuntimeError(
                    f"IdealMHD2DSolver reached max_steps={max_steps} before t_end",
                )

            next_dt = self.compute_dt() if dt is None else dt
            next_dt = min(next_dt, t_end - self.time)
            self.step(next_dt)

            if verbose and self.steps % 25 == 0:
                print(
                    f"MHD step {self.steps:5d}, "
                    f"time={self.time:.6e}, dt={next_dt:.3e}, "
                    f"max|div B|={self.max_abs_divergence_b():.3e}",
                )

        return MHDResult(
            state=self.state.copy(),
            time=self.time,
            steps=self.steps,
            dt_history=self.dt_history.copy(),
            max_divergence_b=self.max_abs_divergence_b(),
        )

    def project_magnetic_field(self) -> None:
        """Project the in-plane magnetic field onto a periodic div(B)=0 subspace.

        The projection uses the Fourier symbol of the centered derivative used
        by ``divergence_b``.  Total energy is then rebuilt so density, velocity,
        and gas pressure remain unchanged by the projection.
        """
        primitives = self.primitive_variables()
        rho = primitives["density"]
        velocity = primitives["velocity"]
        pressure = primitives["pressure"]

        nx, ny = self.grid.num_gridpoints
        qx = np.sin(2.0 * np.pi * np.fft.fftfreq(nx)[:, None]) / self.dx
        qy = np.sin(2.0 * np.pi * np.fft.fftfreq(ny)[None, :]) / self.dy
        q2 = qx * qx + qy * qy

        bx_hat = np.fft.fftn(self._state[..., self.B_X])
        by_hat = np.fft.fftn(self._state[..., self.B_Y])
        q_dot_b = qx * bx_hat + qy * by_hat

        correction = np.zeros_like(q_dot_b)
        mask = q2 > 0.0
        correction[mask] = q_dot_b[mask] / q2[mask]

        bx_hat -= qx * correction
        by_hat -= qy * correction
        self._state[..., self.B_X] = np.real(np.fft.ifftn(bx_hat))
        self._state[..., self.B_Y] = np.real(np.fft.ifftn(by_hat))

        magnetic_field = self._state[..., self.B_X : self.B_Z + 1]
        kinetic_energy = 0.5 * rho * np.sum(velocity * velocity, axis=-1)
        magnetic_energy = 0.5 * np.sum(magnetic_field * magnetic_field, axis=-1)
        self._state[..., self.ENERGY] = (
            pressure / (self.gamma - 1.0) + kinetic_energy + magnetic_energy
        )
        self._enforce_floors()

    def divergence_b(self) -> np.ndarray:
        """Return a centered finite-difference estimate of div(B)."""
        bx = self.state[..., self.B_X]
        by = self.state[..., self.B_Y]
        div_x = np.roll(bx, shift=-1, axis=0) - np.roll(bx, shift=1, axis=0)
        div_x /= 2.0 * self.dx
        div_y = np.roll(by, shift=-1, axis=1) - np.roll(by, shift=1, axis=1)
        div_y /= 2.0 * self.dy
        return div_x + div_y

    def max_abs_divergence_b(self) -> float:
        """Return the maximum absolute centered div(B) diagnostic."""
        return float(np.max(np.abs(self.divergence_b())))

    @property
    def grid(self) -> Grid:
        """Return the grid."""
        return self._grid

    @property
    def gamma(self) -> float:
        """Return the adiabatic index."""
        return self._gamma

    @property
    def cfl(self) -> float:
        """Return the CFL number."""
        return self._cfl

    @property
    def density_floor(self) -> float:
        """Return the density floor."""
        return self._density_floor

    @property
    def pressure_floor(self) -> float:
        """Return the pressure floor."""
        return self._pressure_floor

    @property
    def project_divergence(self) -> bool:
        """Return whether periodic magnetic divergence projection is enabled."""
        return self._project_divergence

    @property
    def dx(self) -> float:
        """Return the cell width along x."""
        return self._dx

    @property
    def dy(self) -> float:
        """Return the cell width along y."""
        return self._dy

    @property
    def time(self) -> float:
        """Return the current simulation time."""
        return self._time

    @property
    def steps(self) -> int:
        """Return the number of completed timesteps."""
        return self._steps

    @property
    def dt_history(self) -> list[float]:
        """Return the timestep history."""
        return self._dt_history

    @property
    def state(self) -> np.ndarray:
        """Return the conservative state array."""
        return self._state


class MultiFluidIdealMHD2DSolver:
    """Coordinate multiple ideal-MHD fluid solvers on one grid."""

    def __init__(
        self,
        grid: Grid,
        fluid_names: list[str] | tuple[str, ...],
        gammas: dict[str, float] | float = 5.0 / 3.0,
        cfl: float = 0.35,
        density_floor: float = 1e-8,
        pressure_floor: float = 1e-8,
        project_divergence: bool = True,
    ) -> None:
        """Initialize the multi-fluid solver."""
        if grid.ndim != 2:
            raise ValueError("MultiFluidIdealMHD2DSolver requires a 2-D grid")
        for name in fluid_names:
            if not isinstance(name, str) or not name:
                raise ValueError("fluid_names entries must be non-empty strings")
        fluid_names = list(fluid_names)
        if len(fluid_names) == 0:
            raise ValueError("fluid_names must contain at least one fluid")
        if len(set(fluid_names)) != len(fluid_names):
            raise ValueError("fluid_names must be unique")

        if isinstance(gammas, dict):
            missing = [name for name in fluid_names if name not in gammas]
            if missing:
                raise ValueError(
                    "Missing gamma values for fluids: " + ", ".join(sorted(missing)),
                )
            gamma_by_fluid = {name: float(gammas[name]) for name in fluid_names}
        else:
            gamma_by_fluid = {name: float(gammas) for name in fluid_names}

        self._grid = grid
        self._fluid_names = fluid_names
        self._solvers = {
            name: IdealMHD2DSolver(
                grid=grid,
                gamma=gamma_by_fluid[name],
                cfl=cfl,
                density_floor=density_floor,
                pressure_floor=pressure_floor,
                project_divergence=project_divergence,
            )
            for name in fluid_names
        }
        self._time = 0.0
        self._steps = 0
        self._dt_history: list[float] = []

    def _get_solver(self, fluid_name: str) -> IdealMHD2DSolver:
        """Return one registered fluid solver."""
        try:
            return self._solvers[fluid_name]
        except KeyError as exc:
            available = ", ".join(self._fluid_names)
            raise KeyError(
                f"Unknown fluid {fluid_name!r}. Available fluids: {available}",
            ) from exc

    def set_state(self, fluid_name: str, state: np.ndarray) -> None:
        """Set conservative state for one fluid."""
        self._get_solver(fluid_name).set_state(state)

    def set_primitive(
        self,
        fluid_name: str,
        density: np.ndarray,
        velocity: np.ndarray,
        magnetic_field: np.ndarray,
        pressure: np.ndarray,
    ) -> None:
        """Set primitive state for one fluid."""
        self._get_solver(fluid_name).set_primitive(
            density=density,
            velocity=velocity,
            magnetic_field=magnetic_field,
            pressure=pressure,
        )

    def primitive_variables(
        self,
        fluid_name: str,
        state: np.ndarray | None = None,
    ) -> dict[str, np.ndarray]:
        """Return primitive variables for one fluid."""
        return self._get_solver(fluid_name).primitive_variables(state=state)

    def fluid_solver(self, fluid_name: str) -> IdealMHD2DSolver:
        """Return the underlying one-fluid solver by name."""
        return self._get_solver(fluid_name)

    def compute_dt(self) -> float:
        """Return the most restrictive CFL timestep across all fluids."""
        return min(solver.compute_dt() for solver in self._solvers.values())

    def step(self, dt: float | None = None) -> float:
        """Advance all fluids by one synchronized timestep."""
        if dt is None:
            dt = self.compute_dt()
        if dt <= 0.0:
            raise ValueError("dt must be greater than zero")

        for solver in self._solvers.values():
            solver.step(dt)

        reference_solver = self._solvers[self._fluid_names[0]]
        for name, solver in self._solvers.items():
            if abs(float(solver.time) - float(reference_solver.time)) > 1e-15 or int(
                solver.steps,
            ) != int(reference_solver.steps):
                raise RuntimeError(
                    f"Fluid solver {name!r} became desynchronized after step()",
                )
        self._time = float(reference_solver.time)
        self._steps = int(reference_solver.steps)
        self._dt_history.append(float(dt))
        return float(dt)

    def solve(
        self,
        t_end: float,
        max_steps: int = 100_000,
        dt: float | None = None,
        verbose: bool = False,
    ) -> MultiFluidMHDResult:
        """Advance all fluid states until ``t_end``."""
        if t_end < self.time:
            raise ValueError("t_end must be greater than or equal to the current time")
        if max_steps <= 0:
            raise ValueError("max_steps must be positive")

        stop_step = self.steps + max_steps
        while self.time < t_end:
            if self.steps >= stop_step:
                raise RuntimeError(
                    "MultiFluidIdealMHD2DSolver reached "
                    f"max_steps={max_steps} before t_end",
                )

            next_dt = self.compute_dt() if dt is None else dt
            next_dt = min(next_dt, t_end - self.time)
            self.step(next_dt)

            if verbose and self.steps % 25 == 0:
                max_div = max(self.max_abs_divergence_b().values())
                print(
                    f"Multi-fluid MHD step {self.steps:5d}, "
                    f"time={self.time:.6e}, dt={next_dt:.3e}, "
                    f"max|div B|={max_div:.3e}",
                )

        return MultiFluidMHDResult(
            states={
                name: solver.state.copy() for name, solver in self._solvers.items()
            },
            time=self.time,
            steps=self.steps,
            dt_history=self.dt_history.copy(),
            max_divergence_b=self.max_abs_divergence_b(),
        )

    def max_abs_divergence_b(self) -> dict[str, float]:
        """Return max |div(B)| by fluid."""
        return {
            name: solver.max_abs_divergence_b()
            for name, solver in self._solvers.items()
        }

    @property
    def grid(self) -> Grid:
        """Return the shared simulation grid."""
        return self._grid

    @property
    def fluid_names(self) -> list[str]:
        """Return all registered fluid names."""
        return self._fluid_names

    @property
    def time(self) -> float:
        """Return the current simulation time."""
        return self._time

    @property
    def steps(self) -> int:
        """Return the number of completed timesteps."""
        return self._steps

    @property
    def dt_history(self) -> list[float]:
        """Return timestep history."""
        return self._dt_history

    @property
    def state(self) -> dict[str, np.ndarray]:
        """Return conservative state by fluid name."""
        return {name: solver.state for name, solver in self._solvers.items()}
