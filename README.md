# Paraboloid Curvature

This is a fork of [Pececillo](https://gitlab.com/truchas/pececillo) which reproduces the results for the paper "A Paraboloid Fitting Technique for Calculating Curvature from Piecewise-Linear Interface Reconstructions on 3D Unstructured Meshes" by Z. Jibben, N. N. Carlson, and M. M. Francois, 2017.

## Building

This depends on the libraries necessary for building Pececillo. A superbuild may be found [here](https://gitlab.com/truchas/pececillo-tpl/tags/v1). Tests successfully run with the Intel compilers v17.0.4.

```bash
$ mkdir build
$ cd build
$ cmake -DCMAKE_PREFIX_PATH=/your/pececillo-tpl/install/path ..
$ make
```

## Running

```bash
$ cd run
$ ../build/src/parabfit/pececillo-parabfit
```

By default, this will calculate the curvature error norms for the ellipsoid test case on the coarsest tetrahedral mesh, which should take 5-20 minutes, depending on your hardware.

To reproduce all results in the paper, you must first download all tetrahedral and distorted hexahedral meshes from [here](https://zenodo.org/record/1135295/files/meshes.tar.gz) and unpack them in the run directory. Then, execute `$ ../build/src/parabfit/pececillo-parabfit full`. This will take days to weeks to calculate on a single workstation, due to the high cost of initializing volume fraction fields, and to a lesser extent calculating normal vectors via LVIRA. You can speed this up by instead downloading volume fraction fields and normal vectors from [here](https://zenodo.org/record/1135295/files/fields.tar.gz), unpacking them in the run directory. By reading in these fields, all paper results will be calculated in a few hours.

Both the fields and meshes tarballs are identified by [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1135295.svg)](https://doi.org/10.5281/zenodo.1135295).

## Paraboloid Fitting Algorithm in the Source

To find the source code for the paraboloid fitting algorithm, look at the following in the `src/lib` directory:

- `surface_type.F90`, function `local_patch` returns a collection of interface reconstruction polygons surrounding a target cell.
- `interface_patch_type.F90`, function `curvature_from_patch` takes a collection of polygons and drives other routines to fit a paraboloid and calculate curvature.
- `paraboloid_type.F90`, subroutine `volumetricFit` fits a paraboloid to a given set of polygons.
- `paraboloid_type.F90`, function `curvature` calculates curvature at a point on the paraboloid.
