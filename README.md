# Firm Dynamics

MATLAB implementations of firm dynamics models with heterogeneous firms, entry/exit, and labor adjustment costs.

Started with Hopenhayn-Neira-Singhania (2022) replication codes. Adapted to my taste on notation and procedure.

Note: Still unclear on whether to use the FOC for the policies or, as we do now, just divide profits and factor payments.

## Models

- [`hopenhayn/`](hopenhayn/) -- Hopenhayn (1992) stationary equilibrium with calibration from Hopenhayn, Neira, and Singhania (2022). 1D state space (productivity only), static labor demand, binary exit decision.
- [`hopenhayn_rogerson/`](hopenhayn_rogerson/) -- Hopenhayn-Rogerson (1993) with firing costs. 2D state space (productivity, employment), dynamic labor choice on a discrete grid, comparative statics across firing cost levels.

## Results

### Hopenhayn (1992) Stationary Equilibrium

| | |
|---|---|
| Entry cost (ce) | 0.0119 |
| Exit rate | 0.0976 |
| Startup rate | 0.0976 |
| Avg firm size | 256.5217 |
| Avg startup size | 36.5263 |
| Avg productivity | 1.5625 |
| Job turnover | 0.2711 |
| Total employment | 15.7443 |
| Output | 20.5838 |

### Hopenhayn-Rogerson (1993) Comparative Statics

| | tau=0 | tau=0.1 | tau=0.5 |
|---|---|---|---|
| Entry cost (ce) | 0.0120 | 0.0110 | 0.0085 |
| Exit rate | 0.0979 | 0.0981 | 0.0788 |
| Startup rate | 0.0981 | 0.0984 | 0.0792 |
| Avg firm size | 250.7925 | 232.5851 | 211.7189 |
| Avg startup size | 36.5306 | 34.8554 | 34.2341 |
| Avg productivity | 1.5624 | 1.5997 | 1.6860 |
| Job turnover | 0.2613 | 0.2216 | 0.1421 |
| Total employment | 15.3490 | 14.1850 | 10.4187 |
| Output | 19.9713 | 18.7454 | 14.2799 |
