# Data

No data files are included in this repository.

## Monte Carlo simulation

`src/heatwave_tvar_sim.py` uses **no external data**. All inputs are illustrative, synthetic parameters (w = 1, M = 0.3, CRRA 3, λ0 = 0.1, q = 0.9, κ = 0.1) and a fixed random seed. Outputs go to `results/` and `reports/figures/`.

## ERA5 reanalysis (preliminary empirical statistics)

The preliminary heat-wave statistics quoted in the main README (e.g. the 0.84 mean inter-regional correlation of annual heat-wave days) come from:

- **Source:** ECMWF ERA5 reanalysis, available from the Copernicus Climate Data Store (https://cds.climate.copernicus.eu/). Free to download after registration, subject to the Copernicus licence.
- **Coverage used:** a representative point in each of Korea's 17 provinces, 1991–2025, daily maximum temperature.
- **Heat-wave day definition:** daily maximum at or above that region's 1991–2020 June–September 95th percentile.

**Why it is not included:** the raw ERA5 downloads are large (NetCDF/GRIB) and are excluded via `.gitignore`. The statistics are preliminary and have not yet been validated against KMA ASOS station data.

**Code:** the ERA5 download and processing scripts are not in this repository yet. They will be added once the station-data validation is done.

## Planned data (not yet available)

- KMA ASOS station records (for bias correction of ERA5).
- Economically Active Population Survey microdata and Construction Workers Mutual Aid Association workday records (requested), for the loss coefficient.
