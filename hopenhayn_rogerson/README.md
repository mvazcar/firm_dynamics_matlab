# Hopenhayn-Rogerson (1993) Comparative Statics

Extension of the Hopenhayn (1992) model with firing costs, following Hopenhayn and Rogerson (1993). The state space is 2D: productivity `s` and current employment `n`. Firms face a per-worker firing cost `tau` when reducing employment, making labor a dynamic choice.

## Workflow

### Entry point

- `main.m` -- Sets parameters, builds the productivity grid, startup distribution, and employment grid, then solves the stationary equilibrium for `tau = 0, 0.1, 0.5` and displays comparative statics.

### Functions called by `main.m`

1. `tauchen.m` -- Discretizes the AR(1) productivity process into a finite Markov chain (Tauchen 1986).
2. `stationary.m` -- Solves for the stationary equilibrium:
   - Calls `vfi()` to get the 2D value function and employment policy.
   - Computes the free-entry cost `ce` (entrants start with `n = 0`).
   - Iterates the 2D distribution `mu(s, n)` until convergence.
   - Computes aggregates (employment, output, job turnover, exit rate).
3. `vfi.m` -- Value function iteration over the 2D grid `(s, n)`. Firms choose next-period employment `nprime` from the grid, facing firing costs `phi(n, nprime) = tau * w * max(0, n - nprime)`.
4. `n.m` -- Static optimal labor demand (used to set the employment grid bounds).
5. `profit.m` -- Static flow profit net of the fixed operating cost `cf`.
6. `phi.m` -- Firing cost function: `phi(nprime, n, tau, w) = tau * w * max(0, n - nprime)`.

### Input parameters (set in `main.m`)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `beta` | 1/1.04 | Discount factor |
| `theta` | 0.64 | Returns to scale |
| `ns` | 100 | Productivity grid points |
| `nn` | 300 | Employment grid points |
| `rho` | 0.9842 | AR(1) persistence |
| `mu` | -1.4361 | AR(1) long-run mean |
| `sigma` | 0.2456 | AR(1) innovation std dev |
| `cf` | 24.308 | Fixed operating cost (labor units) |
| `mu0` | -4.3444 | Startup distribution mean |
| `sigma0` | 1.3311 | Startup distribution std dev |
| `pstar` | 1 | Output price (normalized) |
| `wstar` | 1 | Wage (normalized) |
| `mstar` | 100 | Mass of entrants |
| `tau` | 0, 0.1, 0.5 | Firing cost per worker |

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
| `v` | Value function (ns x nn) |
| `npol_ind` | Employment policy indices (ns x nn) |
| `mustar` | Stationary distribution (ns x nn) |
| `tau` | Firing cost parameter used |
