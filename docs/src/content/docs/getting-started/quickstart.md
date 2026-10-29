---
title: Quickstart
description: Run the Orszag-Tang vortex with mhd-sandbox.
---

First, ensure that `mhd-sandbox` is [installed](/mhd-sandbox/getting-started/installation/).

## Run the Orszag-Tang vortex

```python
import matplotlib.pyplot as plt

from mhd_sandbox import Grid, IdealMHD2DSolver

grid = Grid(lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0), num_gridpoints=(128, 128))
solver = IdealMHD2DSolver(grid, cfl=0.3)
solver.set_state(IdealMHD2DSolver.orszag_tang_initial_state(grid))
result = solver.solve(t_end=0.5)

density = solver.primitive_variables()["density"]
plt.imshow(density.T, origin="lower", extent=(0, 1, 0, 1), cmap="magma")
plt.colorbar(label="density")
plt.show()
```

The state is conservative, `[rho, rho vx, rho vy, rho vz, E, Bx, By, Bz]` per cell, in units
with the magnetic permeability set to one. `primitive_variables()` returns density, velocity,
magnetic field and pressure.

The [examples](/mhd-sandbox/examples/) run the classic test problems, and the
[API reference](/mhd-sandbox/api/mhd_sandbox/) documents every class and function.
