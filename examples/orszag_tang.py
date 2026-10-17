"""Run the Orszag-Tang vortex.

The standard 2D ideal MHD test: smooth, periodic vortices in the velocity and magnetic field
steepen into shocks that interact. The figure shows the density at t = 0.5 on a 128 x 128 grid,
with the magnetic field lines.
"""

import matplotlib.pyplot as plt
import numpy as np

from mhd_sandbox import Grid, IdealMHD2DSolver

grid = Grid(lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0), num_gridpoints=(128, 128))
solver = IdealMHD2DSolver(grid, cfl=0.3)
solver.set_state(IdealMHD2DSolver.orszag_tang_initial_state(grid))
result = solver.solve(t_end=0.5)

primitives = solver.primitive_variables()
nx, ny = grid.num_gridpoints
x = (np.arange(nx) + 0.5) * solver.dx
y = (np.arange(ny) + 0.5) * solver.dy

fig, ax = plt.subplots(figsize=(6, 5), layout="constrained")
image = ax.imshow(
    primitives["density"].T, origin="lower", extent=(0, 1, 0, 1), cmap="magma"
)
ax.streamplot(
    x,
    y,
    primitives["magnetic_field"][..., 0].T,
    primitives["magnetic_field"][..., 1].T,
    color="white",
    density=1.0,
    linewidth=0.4,
    arrowsize=0.5,
)
fig.colorbar(image, ax=ax, label="density")
ax.set(
    xlabel="x",
    ylabel="y",
    title=f"Orszag-Tang vortex, t = {result.time:.2f} ({result.steps} steps)",
)
