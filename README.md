# Factor-Revealing LPs for the Unsplittable CVRP

## Overview

The results in the paper use the $40\times20$ parameter grid (`N=40`).

- `single_lp.py`: Solves the single-setting LP (B) for each parameter pair on a grid and returns the smallest value and its parameter pair (Theorem 1).
- `mixed_lp.py`: Solves the mixed-parameter LP (C), combining all parameter pairs on the same grid (Theorem 2).
- `results/Grid_40x20/`: Contains the four result files for the grid used in the paper.

Results for three additional grids, $10\times5$, $20\times10$, and $30\times15$, are also included for comparison.

## Parameters

The two LP scripts take two command-line arguments: `<generalflag> <N>`.

- `generalflag`: `1` = general vehicle capacity; `0` = fixed vehicle capacity.
- `N`: Grid resolution, which must be an even integer at least 2. The parameter pairs are $p=i/N$ and $\delta=j/N$, with $i=1,\ldots,N$ and $j=1,\ldots,N/2$.

Thus, `N=40` gives the $40\times20$ grid used in the paper.

For `single_lp.py`, $\theta=\delta$ for each parameter pair. For `mixed_lp.py`, $\theta=1/N$. In both scripts, $\sigma=\theta$ for general capacity and $\sigma=0$ for fixed capacity; the rounding parameter is $\gamma=-\ln p$.


## Python Files

### `single_lp.py`

The function `singlelp(p, delta, generalflag, N)` builds and solves one LP for the specified parameter pair and returns its optimal value. The grid loop is in `main`, which calls this function for each pair and retains the smallest value.

The final output contains `generalflag`, `N`, `value`, `p`, `delta`, and `runtime`. The runtime covers the entire grid search, including model construction and solution.

### `mixed_lp.py`

The function `mixlp(generalflag, N)` builds and solves one LP containing all parameter pairs. Demand-case variables are shared across these pairs.

The final output contains `generalflag`, `N`, `value`, and `runtime`. The runtime includes model construction and solution.

In both LP scripts, `value` is the additive term beyond $\alpha$, rather than the full approximation ratio.

## Requirements

- Python 3 (Python 3.11 was used in the experiments).
- Gurobi Optimizer and a compatible `gurobipy` installation, with a valid Gurobi license, for the two LP scripts.

## Computational Environment

The LPs were implemented in Python 3.11 and solved using Gurobi Optimizer 12.0.3 (build `v12.0.3rc0`).

All experiments were conducted on a server with an Intel Xeon Gold 6226R CPU @ 2.90 GHz and 512 GB RAM, running Ubuntu Server 20.04.4 LTS. Gurobi used 4 threads.

## Gurobi Parameters

Both LP scripts use the following settings:

```python
m.setParam('Threads', 4)
m.setParam('Method', 2)
m.setParam('BarOrder', 1)
m.setParam('BarHomogeneous', 1)
m.setParam('NumericFocus', 2)
m.setParam('Crossover', 0)
```

## How to Run

Run the following commands from this directory to reproduce the computations on the $40\times20$ grid used in the paper:

```bash
# Single-setting LPs: general and fixed capacity
python single_lp.py 1 40
python single_lp.py 0 40

# Mixed-parameter LP: general and fixed capacity
python mixed_lp.py 1 40
python mixed_lp.py 0 40
```

To run the additional experiments on the smaller grids, replace `40` by `10`, `20`, or `30`.


### Example Output

The following examples use general capacity and `N=40`; intermediate solver output is omitted.

For `python single_lp.py 1 40`:

```text
generalflag: 1, N: 40, value: 1.600754219895, p: 0.85, delta: 0.15, runtime: 7.801785 seconds
```

For `python mixed_lp.py 1 40`:

```text
generalflag: 1, N: 40, value: 1.5909756650354863, runtime: 121.26805901527405 seconds
```

Runtime and the final numerical digits may vary across runs and environments.

## Results

The four files in `results/Grid_40x20/` correspond to the results in the paper:

| Capacity | Single-setting LP | Mixed-parameter LP |
| --- | ---: | ---: |
| General | 1.600754219895 | 1.5909756650354863 |
| Fixed | 1.485628627078 | 1.48346238729408 |

The corresponding additive terms stated in the paper are 1.601 and 1.486 for the single-setting results, and 1.591 and 1.484 for the mixed-parameter results, respectively.

For comparison, results for the three smaller grids are provided in `results/Grid_10x5/`, `results/Grid_20x10/`, and `results/Grid_30x15/`. These additional experiments are not reported in the paper.


### Summary of All Four Grids

The table below summarizes the 16 LP runs. The parameter pair $(p,\delta)$ is the best pair found by `single_lp.py`; the mixed LP uses all pairs on the corresponding grid. LP values are the additive terms beyond $\alpha$ and are rounded to nine decimal places. Runtimes are in seconds and include model construction and solution.

| Grid | Capacity | Best $(p,\delta)$ | Single LP value | Single runtime (s) | Mixed LP value | Mixed runtime (s) |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| $10\times5$ | General | $(0.7,0.1)$ | 1.612230493 | 0.276 | 1.594827543 | 0.145 |
| $10\times5$ | Fixed | $(0.6,0.2)$ | 1.485825618 | 0.287 | 1.484887085 | 0.154 |
| $20\times10$ | General | $(0.85,0.15)$ | 1.600754220 | 1.483 | 1.591612266 | 2.254 |
| $20\times10$ | Fixed | $(0.6,0.2)$ | 1.485825618 | 1.467 | 1.484299582 | 2.700 |
| $30\times15$ | General | $(0.8,\frac{2}{15})$ | 1.600066622 | 3.950 | 1.591419322 | 22.397 |
| $30\times15$ | Fixed | $(0.6,0.2)$ | 1.485825618 | 3.803 | 1.483916652 | 24.153 |
| $40\times20$ | General | $(0.85,0.15)$ | 1.600754220 | 7.802 | 1.590975665 | 121.268 |
| $40\times20$ | Fixed | $(0.625,0.2)$ | 1.485628627 | 7.582 | 1.483462387 | 124.349 |
