# GILA Integrator

Friedmann equation implementation and RK4 integration for GILA (exponential) and GR in Fortran 2008.

## Featured in...

- TOFILL

## How to use

### Requirements

- A Unix system (for the moment `src/gila.f90` uses the unix version of `mkdir`)
- 'GCC' (tested on `gfortran 14.2.1 20250207`)
- [Fortran Package Manager ('fpm')](https://fpm.fortran-lang.org/)
- `conda` (`conda-forge`) for environment and dependency management (make sure
  `channel_priority: strict` in `.condarc`)
    - 'python'
        - `matplotlib`
        - `numpy`
        - `scipy`
        - [`ford`](https://github.com/Fortran-FOSS-Programmers/ford) (optional - doc creation)

### Compilation & execution

Run

    conda env create
    conda activate gila-integrator

to activate the environment and download all dependencies.

Run

    conda deactivate

or close your terminal when finished to close it.

#### Fortran source

To compile using `gcc` and run,

    fpm run <app name> --profile release

where `app name` can be

- `find_solutions_initial`: simulates example cases
- `gr_integration`: simulates ΛCDM
- `sr_variance_normalized`: generates data for ϵ_(1,2,3) heatmaps
- `find_solutions`: simulates selected cases

#### Python scripts

    python scripts/<script name>
    
where `script name` can be any of

- `universe.py` creates plots for the introduction
- `initial_plots.py` plots data from `find_solutions_initial`
- `heatmap.py` generates the heatmap from `sr_variance_normalized` data
- `selected_plots.py` plots data from `find_solutions`

#### Docs

    ford doc.md
    
## To be fixed

- Switch directory creation to OS-agnostic implementation (e.g. Python script)
- Make parallelization possible
- Fix output formatting for other compilers (output format not usable when
  compiling with 'ifort')
