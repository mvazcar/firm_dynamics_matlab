# Hopenhayn (1992) Stationary Equilibrium

Implementation of the stationary equilibrium from Hopenhayn (1992), with calibration following Hopenhayn, Neira, and Singhania (2022). Firms use labor as the only input, face idiosyncratic productivity shocks, and make binary exit decisions.

## Workflow

### Entry point

- `main.m` -- Sets parameters, builds the productivity grid and startup distribution, calls `stationary()`, and displays results.

### Functions called by `main.m`

1. `tauchen.m` -- Discretizes the AR(1) productivity process into a finite Markov chain using the Tauchen (1986) method.
2. `stationary.m` -- Solves for the stationary equilibrium:
   - Calls `vfi()` to get the value function and exit threshold.
   - Computes the free-entry cost `ce`.
   - Solves the stationary firm distribution via matrix inversion.
   - Computes aggregates (employment, output, job turnover, exit rate).
3. `vfi.m` -- Value function iteration over the 1D productivity grid. Firms choose to continue or exit each period.
4. `n.m` -- Static optimal labor demand from the firm FOC: `n(s) = ((p * exp(s) * theta) / w) ^ (1/(1-theta))`.
5. `profit.m` -- Static flow profit net of the fixed operating cost `cf`.

### Input parameters (set in `main.m`)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `beta` | 1/1.04 | Discount factor |
| `theta` | 0.64 | Returns to scale |
| `ns` | 100 | Productivity grid points |
| `rho` | 0.9842 | AR(1) persistence |
| `mu` | -1.4361 | AR(1) long-run mean |
| `sigma` | 0.2456 | AR(1) innovation std dev |
| `cf` | 24.308 | Fixed operating cost (labor units) |
| `mu0` | -4.3444 | Startup distribution mean |
| `sigma0` | 1.3311 | Startup distribution std dev |
| `pstar` | 1 | Output price (normalized) |
| `wstar` | 1 | Wage (normalized) |
| `mstar` | 100 | Mass of entrants |

### Output

The `stationary()` function returns a struct with:

| Field | Description |
|-------|-------------|
| `ce` | Free-entry cost (labor units) |
| `exit_rate` | Fraction of incumbents exiting per period |
| `startup_rate` | Entry rate relative to incumbent mass |
| `avg_fsize_n` | Average firm size (total workers) |
| `avg_stsize_n` | Average startup size (total workers) |
| `avg_productivity` | Output-weighted average productivity |
| `job_turnover` | (Job creation + destruction) / (2 * production employment) |
| `N` | Total employment (production + fixed cost + entry) |
| `Y` | Total output |
| `v` | Value function (ns x 1) |
| `mustar` | Stationary distribution (ns x 1) |
| `sstar_ind` | Exit threshold index |
