# GILA Integrator

Friedmann equation implementation and RK4 integration for GILA (exponential) and GR in Fortran 2008.

## Featured in...

- TOFILL

## How to use

### Requirements

- A Unix system (for the moment `src/gila.f90` uses the unix version of `mkdir`)
- Fortran 2008 compiler (tested on `gfortran 14.2.1 20250207` and `ifx 2025.1.0 20250317` (conda))
- `conda` (`conda-forge`) for environment and dependency management (make sure
  `channel_priority: strict` in `.condarc`)
    - 'python'
        - `matplotlib`
        - `numpy`
        - `scipy`
    - [Fortran Package Manager ('fpm')](https://fpm.fortran-lang.org/)
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

To compile using `gcc`,

    fpm run sr_invariance
    fpm run gr_integration --flag '-std=f2008' --profile release
    fpm run curve_sandwich --flag '-std=f2008' --profile release

For other compatible compilers, replace `-std=f2008` witht the apropiate flag for 2008 standard
compliance (e.g. `-std08` for `ifx`).

#### Scripts

If using `bash` or similar, you may use

    for SCRIPT in $(ls scripts) ; do echo --$SCRIPT-- && python scripts/$SCRIPT ; done

#### Docs

    ford doc.md
    
## To be fixed

-  Make parallelization possible
