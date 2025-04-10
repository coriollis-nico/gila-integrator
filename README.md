# GILA Integrator

Friedmann equation implementation and RK4 integration for GILA (exponential) and GR in Fortran 2008.

## Featured in...

- TOFILL

## How to use

### Requirements

- Fortran 2008 compiler (e.g. `gfortran`, tested on `14.2.1 20250207`)
- `conda` for environment and dependency management
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

    fpm run gr_integration --flag -std=f2008 --profile release
    fpm run curve_sandwich --flag -std=f2008 --profile release

For other compatible compilers, replace `-std=f2008` witht the apropiate flag for 2008 standard
compliance.

#### Scripts

If using `bash` or similar, you may use

    for SCRIPT in $(ls scripts) ; do echo --$SCRIPT-- && python scripts/$SCRIPT ; done

#### Docs

    ford doc.md

