# Mineye-Terranigma — Agent Guide

## Project Overview

Bayesian segmentation of satellite imagery + probabilistic 3D geological modeling and geophysical joint inversion with **GemPy**. Developed by **Terranigma Solutions GmbH**.

Two core pipelines:
- **Bayesian Segmentation** — EnMAP satellite imagery → lithological classification with BaySeg
- **Probabilistic Geological Modeling** — 3D implicit modeling + geophysical inversion (gravity, magnetics, EnMAP labels) via GemPy + Pyro

---

## Essential Commands

| Command | Purpose |
|---|---|
| `pip install -e .` | Install package in editable mode (uses `setup.py` with `setuptools_scm` versioning) |
| `python -m pytest tests/ -v --tb=short` | Run all tests |
| `python -m pytest tests/test_data_processing.py -v` | Run specific test file |
| `python -m pytest tests/tests_inversions/` | Run inversion tests (slow, need data) |
| `cd docs && make html` | Build Sphinx documentation |
| `cd docs && python -m sphinx -D plot_gallery=0 -b html source build/html` | Build docs without executing gallery plots (faster) |
| `bash run_docs.sh` | Build docs from project root using venv `~/.venv/2025` |
| `python -m pytest --co` | List collected tests (dry run) |
| `python -m pytest -k "test_name"` | Run by test name pattern |

**Test configuration** (`pytest.ini` / `pyproject.toml`): `-v --tb=short`, `console_output_style = progress`. Test paths include both `tests/` and `examples/`.

---

## Code Organization

```
mineye/                          # Main package
├── GeoModel/                    # Geological modeling & inversion
│   ├── model_one/               # Tharsis Plutonites model
│   │   ├── model_setup.py       # read_gravity(), setup_geomodel(), baseline()
│   │   ├── probabilistic_model.py       # Priors, orientation modifiers, normalization
│   │   ├── probabilistic_model_likelihoods.py  # Likelihood functions (gravity, enmap, magnetic)
│   │   ├── joint_probabilistic_model.py        # Joint inversion model
│   │   ├── deterministics.py    # @dataclass deterministic functions
│   │   ├── inference_diagnostics.py  # check_mcmc_quality(), check_likelihood_balance()
│   │   └── visualization.py     # Gravity uncertainty plots, 3D viz
│   ├── geophysics.py            # compute/align_forward_to_observed() — quantile alignment
│   ├── helper_methods.py        # Data processing (orientations, simplification, boundary removal)
│   ├── helper_plotter.py        # Basic plotting helpers
│   ├── plotting/
│   │   └── probabilistic_analysis.py  # Advanced comparison/uncertainty plots
│   └── Old_Scripts/             # DEPRECATED — legacy code, don't use
├── BayesianSegmentation/        # Satellite imagery processing
│   ├── EnMap_feature_extraction.py  # EnMAP band parsing/wavelengths
│   ├── full_workflow.py         # BaySeg integration
│   └── prepare_data.py          # GeoTIFF cropping/sampling
├── GisHelpers/                  # GIS utilities
│   ├── crs_utils.py
│   ├── extractPointsFromMap.py
│   ├── raster2obj.py
│   └── shapefiles2numpy.py
└── config/
    ├── paths.py                 # Path resolution (data, geomodel, geophysical dirs)
    └── example_parameters.py    # TharsisModelConfig, SoricomModelConfig
examples/
├── 01_basic_examples/           # Simple deterministic models
│   ├── 01_simple_tharsis_model.py
│   ├── 02_complex_tharsis_model.py
│   ├── 03_gravity_forward_model.py
│   └── 04_soricom_fault_model.py
├── 02_probabilistic_modeling/   # Bayesian inversion examples
│   ├── 04_gravity_inversion.py      # Gold-standard template (1368 lines)
│   ├── 05_magnetics_inversion.py
│   ├── 06_enmap_inversion.py
│   └── 07_joint_inversion.py
├── 03_segmentation/             # EnMAP lithology segmentation
└── Segmentation/                # Standalone (deprecated) segmentation scripts
tests/
├── conftest.py                  # Empty (avoids pytest import issues)
├── test_data_processing.py      # Recreates temp files from raw data
├── memory_consumption_test.py   # BaySeg memory profiling
├── test_benchmarking/
│   └── test_benchmark_I.py      # Basic model compute benchmark
└── tests_inversions/            # Integration tests for inversion
    ├── conftest.py              # Key fixtures: simple_geo_model, topography_dir, etc.
    └── Model 1/
        ├── test_structural_model.py
        ├── test_geophysics.py
        ├── test_gravity_inversion.py      # NUTS inference
        ├── test_gravity_inversion_vi.py   # Variational inference (normalizing flows)
        ├── test_magnetics_inversion.py
        ├── test_enmap_inversion.py
        ├── test_enmap_preprocess.py
        ├── test_enmap_residuals.py
        ├── test_enmap_likelihood.py
        ├── test_joint_inversion.py
        └── test_error_propagation*.py
docs/                            # Sphinx documentation (adapted from GemPy)
```

---

## Application Architecture

### Inversion Workflow (gravity → template for all others)

```
1. Read gravity data (GeoJSON) → subset via K-means clustering
2. Setup GeoModel  → PyTorch backend → centered grid → compute forward model
3. Baseline forward model → print statistics
4. Compute alignment params (quantile mapping from forward→observed)
5. Define priors (dips: Normal, density: Normal, .to_event(1))
6. Define deterministics (raw gravity, aligned gravity, mean, max)
7. Define likelihood function (diagonal / hierarchical per-station / hierarchical spatial)
8. Build probabilistic model via gempy_probability.make_gempy_pyro_model[_extended]()
9. Run prior predictive checks
10. Run MCMC inference (NUTS) or Variational Inference (SVI + AutoIAFNormal)
11. Posterior analysis → trace plots, uncertainty maps, cross-sections
12. Save results to NetCDF via arviz
```

### Control Flow

`make_gempy_pyro_model` orchestrates: sample from priors → `set_interp_input_fn(samples, geo_model)` → GemPy forward compute → `likelihood_fn(solutions)` → Pyro distribution → HMC/NUTS gradients.

### Data Flow

- **Paths** are resolved via `mineye.config.paths` (get_data_dir, get_tmp_dir, get_geophysical_dir, etc.)
- Model input data: `examples/Data/Model_Input_Data/Simple-Models/` (points_mod.csv, orientiations_mod.csv)
- Geophysical data: `examples/Data/General_Input_Data/Geophysical_Cleaned_Data/`
- Pre-computed inversion results: `.nc` files in test directories (loaded via `az.from_netcdf`)

---

## Key Patterns & Conventions

### Code Style
- Google-style NumPy docstrings with Parameters/Returns/Notes sections
- `# noinspection PyUnusedImport` comments above unused imports
- F-string formatting throughout
- Dataclass-based deterministic functions (`@dataclass(frozen=True)` with `__call__`)
- Config classes use inner classes for grouping (TharsisModelConfig.TharsisGravityConfig)
- Sphinx Gallery formatting: `# %%` cell separators, `# Commentary text as documentation`

### Probabilistic Programming Patterns
- **Seed**: `pyro.set_rng_seed(4003)`, `torch.manual_seed(4003)` — always set both
- **Priors**: `dist.Normal(loc=..., scale=..., validate_args=True).to_event(1)` for vector-valued
- **Gravity alignment**: `-solutions.gravity` (negate sign), then `align_forward_to_observed()`
- **Normalization methods**: `"align_to_reference"` (mean/std, simpler) or `"quantile_align"` (robust)
- **Likelihood functions** live in `probabilistic_model_likelihoods.py`, are generated via factory functions
- **Deterministics** are lambdas passed as dict: `{"name": lambda samples, gm, sol: ...}`
- **NUTSConfig**: `step_size=0.0001, target_accept_prob=0.65, max_tree_depth=5, init_strategy='median'`
- **PyTorch backend**: Always call `BackendTensor.change_backend_gempy(engine_backend=PYTORCH)` before forward models
- **Gradient preservation**: use `torch.where()` instead of boolean indexing, `torch.clamp()` instead of clipping

### Geophysical Data
- Gravity field name: `VALU_BOU267` (in **mGal**), immediately converted to μGal (`* 1000`)
- Reserved density contrast convention: `density_plutonites - 2.67` (reduction to Bouguer plate)
- Magnetics: susceptibility in SI units, field name `VALUE`
- Spatially distributed subset via `KMeans(n_clusters=20)` + nearest-neighbor to cluster centers

### Testing Patterns
- Tests are **class-based** (`class TestProbabilisticInversion:`)
- Shared fixture in `tests/tests_inversions/conftest.py`: `simple_geo_model`, `topography_dir`, `geophysical_dir`, etc.
- Some tests **skip** if data files missing: `if not os.path.exists(...): pytest.skip("...")`
- Inversion tests save costly results: `data.to_netcdf(filename)` — reuse with `az.from_netcdf(filename)`
- Integration tests are long-running — mark as slow if adding more

---

## Important Gotchas

1. **Dual conftest files**: `tests/conftest.py` is empty (prevents import issues). Real fixtures live in `tests/tests_inversions/conftest.py`.
2. **gempy_probability implicit dependency**: needs `dotenv` loaded (`dotenv.load_dotenv()`) — done automatically in the inversion conftest.
3. **Legacy path resolution**: Some path functions accept optional `base_dir` for backward compatibility — prefer the parameterless versions.
4. **bayseg package**: Referenced in `full_workflow.py` and examples but **not in this repository** — it's an external dependency.
5. **Setuptools SCM**: Version is auto-generated by `setuptools_scm` from git tags. Falls back to `0.1.0`.
6. **Custom vs Centered grids**: Gravity uses `CENTRED` grid, EnMAP uses `CUSTOM` grid with extracted coordinates. Grid type switching is explicit.
7. **Sigmoid slope**: `geo_model.interpolation_options.sigmoid_slope = 100` is set in inversion tests for sharper boundaries.
8. **Octree vs Full resolution**: Models use octree refinement (level 4-5). For inversion, explicitly set `mesh_extraction = False`.
9. **Soricom model variants**: Three config classes — simple (no fault), with fault, erosive. Each has different extent/stratigraphy.
10. **`pyproject.toml` vs `setup.py`**: Both exist. `setup.py` has the authoritative install config including `use_scm_version`. `pyproject.toml` has test config.
11. **Gravity sign convention**: Solutions.gravity is negated: `-solutions.gravity` before alignment.
12. **Formation colors**: Hardcoded dict in `TharsisDataProcessingConfig.FORMATION_COLORS`.
13. **`gempy` API**: Uses `gp.data.GeoModel`, `gp.data.AvailableBackends`, solutions from `gp.compute_model(gempy_model)`.
14. **Python version**: `setup.py` requires `>=3.10`, `pyproject.toml` says `>=3.9`. Actual minimum is 3.10 due to PyTorch.