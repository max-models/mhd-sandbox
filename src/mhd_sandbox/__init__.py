"""A sandbox for learning ideal MHD: finite-volume solvers on a periodic 2D grid."""

from mhd_sandbox.grid import Grid
from mhd_sandbox.solver import (
    IdealMHD2DSolver,
    MHDResult,
    MultiFluidIdealMHD2DSolver,
    MultiFluidMHDResult,
)

__version__ = "0.1.0"

__all__ = [
    "Grid",
    "IdealMHD2DSolver",
    "MHDResult",
    "MultiFluidIdealMHD2DSolver",
    "MultiFluidMHDResult",
]
