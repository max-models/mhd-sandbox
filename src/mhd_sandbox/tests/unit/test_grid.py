import pytest

from mhd_sandbox import Grid


def test_grid_spacing_and_validation() -> None:
    grid = Grid(lower_bound=(0, -1), upper_bound=(2, 1), num_gridpoints=(4, 8))
    assert grid.ndim == 2
    assert grid.num_gridpoints == (4, 8)
    assert grid.spacing == pytest.approx((0.5, 0.25))
    with pytest.raises(ValueError, match="same length"):
        Grid(lower_bound=(0, 0), upper_bound=(1, 1), num_gridpoints=(4,))
    with pytest.raises(ValueError, match="exceed"):
        Grid(lower_bound=(0, 1), upper_bound=(1, 1), num_gridpoints=(4, 4))
