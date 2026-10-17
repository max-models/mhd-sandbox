# mhd-sandbox


<!-- README.md is generated from README.qmd: edit the .qmd and run `make readme`. -->

[![Tests](https://github.com/max-models/mhd-sandbox/actions/workflows/test_pytest.yml/badge.svg)](https://github.com/max-models/mhd-sandbox/actions/workflows/test_pytest.yml)
[![Static
analysis](https://github.com/max-models/mhd-sandbox/actions/workflows/static_analysis.yml/badge.svg)](https://github.com/max-models/mhd-sandbox/actions/workflows/static_analysis.yml)
[![Docs](https://github.com/max-models/mhd-sandbox/actions/workflows/docs.yml/badge.svg)](https://max-models.github.io/mhd-sandbox/)
[![codecov](https://codecov.io/gh/max-models/mhd-sandbox/branch/main/graph/badge.svg)](https://codecov.io/gh/max-models/mhd-sandbox)
[![PyPI](https://img.shields.io/pypi/v/mhd-sandbox.png)](https://pypi.org/project/mhd-sandbox/)
[![Python](https://img.shields.io/pypi/pyversions/mhd-sandbox.png)](https://pypi.org/project/mhd-sandbox/)

A sandbox for learning ideal magnetohydrodynamics by writing and
improving a solver. It holds a deliberately simple finite-volume code
for the 2D ideal MHD equations on a periodic grid (first-order Rusanov
fluxes, forward Euler in time, a Fourier projection that keeps div B =
0) and notebooks that run the classic test problems with it: the
Orszag-Tang vortex, the MHD rotor and a magnetized blast wave. It needs
only NumPy.

The solver started in [picnic](https://github.com/max-models/picnic),
where it was frozen; it lives here so that it can change freely without
touching the max-models simulation codes.

Documentation: <https://max-models.github.io/mhd-sandbox/>

``` python
from mhd_sandbox import Grid, IdealMHD2DSolver

grid = Grid(lower_bound=(0.0, 0.0), upper_bound=(1.0, 1.0), num_gridpoints=(128, 128))
solver = IdealMHD2DSolver(grid, cfl=0.3)
solver.set_state(IdealMHD2DSolver.orszag_tang_initial_state(grid))
result = solver.solve(t_end=0.5)
print(result.steps, result.max_divergence_b)
```

## Where to go from here

Each step is a standard lesson in MHD numerics with a check you can run:

1.  **Brio-Wu shock tube** (1D): the MHD Riemann problem and its
    compound waves, against published reference solutions.
2.  **Higher order**: MUSCL reconstruction with a limiter and SSP-RK2 in
    time; measure the order of convergence on a circularly polarized
    Alfvén wave, an exact nonlinear solution.
3.  **Better fluxes**: Rusanov, then HLL, then HLLD, and how sharply
    each resolves contact and Alfvén discontinuities.
4.  **div B**: projection against constrained transport and GLM cleaning
    on Orszag-Tang.
5.  **Physics**: Kelvin-Helmholtz with a guide field, and reconnection
    once resistivity is in.

[Athena++](https://github.com/PrincetonUniversity/athena) is the usual
reference for the results.

## Install

With [uv](https://docs.astral.sh/uv/):

``` bash
make install    # uv sync --extra dev, plus the pre-commit hooks
```

Or with pip, in a Python 3.10 or newer environment:

``` bash
pip install -e ".[dev]"
```

The `test`, `docs` and `dev` extras install the test runner, the
documentation tooling (Matplotlib and Jupyter, for the tutorials) and
the linters.

## Development

Formatting and linting use [ruff](https://docs.astral.sh/ruff/), type
checking [pyright](https://microsoft.github.io/pyright/), run by
[pre-commit](https://pre-commit.com/) and in CI:

``` bash
make lint     # ruff check, ruff format --check, pyright
make test     # pytest with coverage
```

Commit messages follow [Conventional
Commits](https://www.conventionalcommits.org/); see
[CONTRIBUTING.md](CONTRIBUTING.md).

## Build docs

The documentation in `docs/` is an Astro + Starlight site: hand-written
pages, the notebooks in `tutorials/` executed and published as pages,
complete Python scripts in `examples/` published with their plots, and
the API reference generated from the docstrings with
[starlight-pydocs](https://ewels.github.io/starlight-pydocs/). It needs
Node 22 or newer.

``` bash
make docs-install     # npm packages and the Python docs extra
make docs-notebooks   # execute tutorials/*.ipynb and convert them to pages
make docs-examples    # run examples/*.py and publish their plots
make docs-dev         # live preview at http://localhost:4321/mhd-sandbox/
make docs-build       # the static site in docs/dist
```

## Build the README

`README.md` is rendered from `README.qmd` with
[Quarto](https://quarto.org/):

``` bash
make readme
```

## Releases

Before merging a release to `main`, update the version in
`pyproject.toml`, `src/mhd_sandbox/__init__.py` and `CITATION.cff`
(including its release date), and add the release notes to
`CHANGELOG.md`. The push to `main` creates a GitHub release with a
`vX.Y.Z` tag and publishes the package to PyPI with trusted publishing
(OIDC). The one-time PyPI and GitHub configuration is described in the
[publishing
tutorial](https://max-models.github.io/mhd-sandbox/development/publishing/).
