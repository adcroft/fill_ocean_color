# Fill SeaWIFS monthly climatological dataset

- The SeaWIFS data has many gaps where the satellite imagery was obscured.
- To use ocean color in a model we need a complete dataset and so use OI
  to fill in missing data.
- We fill in missing data on the original SeaWIFS grid by generating an
  ocean/land mask from GEBCO bathymetry.

## Usage:
```
make
```
does everything (at GFDL).

```
./fill_and_join_chlor_a.py -h
```
gives usage information.
