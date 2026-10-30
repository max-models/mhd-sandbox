---
title: Installation
description: Create a Python environment and install mhd-sandbox.
---

Create and activate a Python environment:

```bash
python -m venv env
source env/bin/activate
pip install --upgrade pip
```

Install the package and its requirements with pip:

```bash
pip install -e .
```

Check it with:

```bash
python -c "import mhd_sandbox; print(mhd_sandbox.__version__)"
```

The examples also need Matplotlib and Jupyter, which the `docs` extra installs.

## Optional extras

| Extra  | Installs                                                |
| ------ | ------------------------------------------------------- |
| `test` | `pytest` and `coverage`                                 |
| `docs` | Matplotlib, the notebook runner and griffe for the API reference |
| `dev`  | formatters and linters, plus the `test` and `docs` extras |

```bash
pip install -e ".[dev]"
```
