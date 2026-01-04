# pyDots3DMP

Python codes for dots3DMP experiments modelling and analysis.

### Basic Usage

#### Environment setup 🔧

- Option A — **recommended**: install `uv` and use `uv sync` to install dependencies directly from `pyproject.toml`:
  1. Install `uv` if needed:
     ```bash
     curl -LsSf https://astral.sh/uv/install.sh | sh  # macOS / Linux
     ```
     Visit [Getting Started](https://docs.astral.sh/uv/getting-started/installation/) for more information.
  2. From the project root, sync and install dependencies - uv will automatically create a venv in the project workspace.
     ```bash
     uv sync
     ```

- Option B — Standard venv creation and pip installation:
    ```bash
    python -m venv <env_name>
    source <env_name>/bin/activate  # macOS / Linux
    pip install -r requirements.txt
    ```

- Option C — conda-based environment (alternative):
  ```bash
  conda create --name <env_name> --file environment.yml
  # if that fails:
  conda env create -f environment.yml
  conda activate <env_name>
  ```

---

#### Project structure 📁

- `behavior/` — behavioral models and helpers (e.g., `selfmotionddm.py`, `Accumulator.py`, `preprocessing.py`)
- `neural/` — neural analyses and helpers
- `SJscripts/` — example wrapper scripts and utilities
- `scraps/` — experimental/testing scripts
- `archive/` — deprecated code

---

#### Using `SelfMotionDDM` (behavior/selfmotionddm.py) ✅

`SelfMotionDDM` implements the 3DMP accumulator model with optional confidence/wager readouts. Example workflow:

1. Instantiate the model:
```python
import numpy as np
from behavior.selfmotionddm import SelfMotionDDM

# set diffusion grid resolution
grid_vec = np.arange(-3, 0, 0.05)
time_vec = np.arange(0, 2, 0.01)

# define initial parameters
init_params = {
    'kmult': [0.6, 0.6],
    'bound': [1.0, 1.0, 1.0],
    'non_dec_time': [0.3],
    'wager_thr': [1, 1, 1],
    'wager_alpha': [0.05],
}
ddm = SelfMotionDDM(grid_vec=grid_vec, tvec=time_vec, **init_params,
                      stim_scaling=False, return_wager=False)
```

2. Prepare data (use provided helpers in `behavior.preprocessing`):
- `data = data_cleanup(filepath)`
- `data = format_onetargconf(data, remove_one_targ=True)`

Data shapes expected:
- X DataFrame: columns `['modality', 'coherence', 'delta', 'heading']`
- y DataFrame: columns `['choice', 'PDW', 'RT']` (RT is used for RT likelihood computations)

3. Fit the model:
```python
# optionally specify parameters to hold fixed during fitting
accum.fit(X, y, fixed_params=['kmult'])
```
- `fixed_params` is an optional list of parameter names to keep constant during optimization - these do not get passed to the optimization call but are used by the `predict` method.
- By default, optimization uses `pybads.BADS` (requires `pybads`); the code also contains a hook to use `scipy.optimize.minimize`.

4. Generate predictions:
```python
y_pred, y_pred_samp = ddm.predict(X, n_samples=1, cache_accumulators=True, seed=1)
```
- `predict` returns a DataFrame with predicted likelihoods for `choice`, `PDW`, and `RT` and an optional sampled predictions DataFrame when `n_samples > 0`.

Key parameters (see `behavior/selfmotionddm.py` for defaults):
- `kmult` — k multipliers per modality (e.g., `[k_ves, k_vis]`)
- `bound` — bounds per modality (e.g., `[ves, vis, comb]`)
- `non_dec_time` — non-decision times per modality
- `wager_thr` — wager threshold(s)
- `wager_alpha` — wager mapping alpha(s)
- `return_wager` — whether to compute wager/confidence predictions
- `stim_scaling` — whether to compute stimulus-driven urgency signals (or pass tuple of urgency arrays)

---

#### Notes & roadmap

- The `Accumulator` class and helpers (`behavior/Accumulator.py`, `behavior/moi.py`) provides the low-level method-of-images computations (CDF/PDF, RT distributions, and log posterior odds) that `SelfMotionDDM` uses.
- The older `ddm_2d` codebase is deprecated; you may still find legacy code under `archive/` but active analyses should use `SelfMotionDDM` and `Accumulator`.

---

## TO DO

1. Some experimental features (cue-combination strategies, different confidence mappings) are skeletons and not fully implemented or tested.
2. Split larger functions into smaller responsibilities (e.g., accumulator setup vs prediction)
3. add unit tests.
4. Improve and expand documentation.
5. Improve and expand diagnostic visualizations of accumulators.
6. Improve logging/saving of json results and iteration history for checks and resuming optimization runs.
