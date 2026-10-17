"""A periodic, cell-centred rectangular grid."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Grid:
    """A box ``lower_bound`` to ``upper_bound`` split into ``num_gridpoints`` cells.

    The values live at the cell centres; every axis is periodic.
    """

    lower_bound: tuple[float, ...]
    upper_bound: tuple[float, ...]
    num_gridpoints: tuple[int, ...]

    def __post_init__(self) -> None:
        """Normalise the bounds and cell counts to tuples and check them."""
        lower = tuple(float(v) for v in self.lower_bound)
        upper = tuple(float(v) for v in self.upper_bound)
        cells = tuple(int(n) for n in self.num_gridpoints)
        if not len(lower) == len(upper) == len(cells):
            raise ValueError("bounds and num_gridpoints must have the same length")
        if any(u <= lo for lo, u in zip(lower, upper, strict=True)):
            raise ValueError("upper_bound must exceed lower_bound on every axis")
        if any(n < 1 for n in cells):
            raise ValueError("num_gridpoints must be positive")
        object.__setattr__(self, "lower_bound", lower)
        object.__setattr__(self, "upper_bound", upper)
        object.__setattr__(self, "num_gridpoints", cells)

    @property
    def ndim(self) -> int:
        """Return the number of axes."""
        return len(self.num_gridpoints)

    @property
    def spacing(self) -> tuple[float, ...]:
        """Return the cell widths."""
        return tuple(
            (u - lo) / n
            for lo, u, n in zip(
                self.lower_bound, self.upper_bound, self.num_gridpoints, strict=True
            )
        )
