# one (Regridding data files for CCAM)

One is used to convert data from lat/lon to cubic grid so that it can be
used by the Conformal Cubic Atmospheric Model (CCAM).  Typically one is used
to prepare sea surface temperatures (SSTs) or related datasets for reading
with CCAM.

## Website

For documentation, see our website at

[https://research.csiro.au/ccam/]

## Dependencies

One requires the NetCDF C library.

## Building one

To build one with intel, gnu and cray fortran compiler use

```
make
make GFORTRAN=yes
make CRAY=yes
```

Debugging is also enabled with

```
make TEST=yes
```
