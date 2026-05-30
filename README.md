# MC Propagator

Monte Carlo orbit propagation tool for small Solar System bodies. The project queries JPL Horizons/SBDB for nominal orbital elements and uncertainties, generates Monte Carlo clones, propagates them through GMAT, merges the resulting ephemerides, and exports averaged trajectories with standard deviations.

## Features

* Query osculating orbital elements from JPL Horizons
* Query orbital-element uncertainties from JPL SBDB
* Generate Monte Carlo Keplerian-element clones
* Run parallel GMAT simulations using Python multiprocessing
* Export individual clone ephemerides as CSV files
* Merge uneven-length GMAT reports using cubic-spline interpolation
* Save a full ephemeris CSV with mean and sigma columns
* Plot heliocentric and geocentric orbital evolution

## Project structure

```text
MC Propagator/
├── main.py              # Main Monte Carlo execution script
├── MC_handler.py        # Monte Carlo clone generation, GMAT run wrapper, data merging
├── gmat_handler.py      # JPL Horizons query, GMAT script parsing/loading/execution
├── load_gmat.py         # GMAT Python API bootstrap configuration
├── miscfuncs.py         # Report loading, orbital conversions, uncertainty math helpers
├── config.py            # Shared paths, date formats, constants, report configuration
├── make_plots.py        # Plot generation from saved full ephemeris CSV
├── analysis.py          # Placeholder for analysis code
├── test.py              # Experimental plotting/test script
├── debug.ipynb          # Debug notebook
├── ssd.jpl.nasa.gov.crt # Certificate used for JPL requests
└── GmatAPILog.txt       # GMAT API log
```

## Requirements

This project requires Python 3.10 and a working GMAT installation with the GMAT Python API enabled.

Python packages used by the project include:

```text
numpy
pandas
scipy
astropy
astroquery
uncertainties
matplotlib
```

GMAT is not installable through `pip`; it must be installed separately.

## Setup

Clone or unzip the project, then create and activate a virtual environment:

```bash
python -m venv .venv
```

On Windows:

```bash
.venv\Scripts\activate
```

On macOS/Linux:

```bash
source .venv/bin/activate
```

Install dependencies:

```bash
pip install numpy pandas scipy astropy astroquery uncertainties matplotlib
```

## GMAT configuration

Edit `load_gmat.py` so `GmatInstall` points to your local GMAT installation:

```python
GmatInstall = "E:\\Projects\\MEng Project\\GMAT"
```

The script expects the GMAT startup file here:

```text
<GmatInstall>/bin/api_startup_file.txt
```

If GMAT is configured correctly, importing `load_gmat.py` should add the GMAT binary folder to `sys.path` and initialize `gmatpy`.

## Required folders

Before running simulations, create the expected output folders:

```bash
mkdir scripts
mkdir data
mkdir -p MCWorkspace/Scripts
mkdir -p MCWorkspace/Data
```

On Windows PowerShell:

```powershell
New-Item -ItemType Directory -Force scripts
New-Item -ItemType Directory -Force data
New-Item -ItemType Directory -Force MCWorkspace/Scripts
New-Item -ItemType Directory -Force MCWorkspace/Data
```

## GMAT template script

`main.py` expects a GMAT template script at:

```text
scripts/template.script
```

The template should contain GMAT objects that `gmat_handler.py` can parse, including:

* `Spacecraft`
* `ForceModel`
* `Propagator`
* `CoordinateSystem`
* `ReportFile`
* `BeginMissionSequence`

The parser expects GMAT script lines in formats such as:

```text
Create Spacecraft <name>;
GMAT <name>.<property> = <value>;
BeginMissionSequence;
Propagate <PropagatorName>(<SpacecraftName>) {...};
```

## Usage

Configure the target and simulation settings in `main.py`:

```python
obj_name = '2022 NX1'
start_epoch = '2021-07-01'
step = 1 * u.day
duration = 5 * u.year
n = 500
```

Run the Monte Carlo propagation:

```bash
python main.py
```

The script will:

1. Query JPL Horizons for the target’s initial orbital elements.
2. Query JPL SBDB for orbital-element uncertainties.
3. Generate `n` Monte Carlo clones.
4. Run each clone through GMAT in parallel.
5. Save individual clone CSV files in `MCWorkspace/Data`.
6. Merge the results onto a common time axis.
7. Save the final averaged ephemeris to:

```text
data/Full ephemeris <object name>.csv
```

Example:

```text
data/Full ephemeris 2022 NX1.csv
```

## Plotting results

After a full ephemeris CSV has been generated, update the object name in `make_plots.py` if needed:

```python
obj_name = '2024 PT5'
```

Then run:

```bash
python make_plots.py
```

This generates:

```text
data/<object name> heliocentric plots.png
data/<object name> geocentric plots.png
```

## Output columns

The default reported quantities are defined in `config.py`:

```python
default_params = [
    "UTCGregorian",
    "Earth.C3Energy",
    "Earth.Altitude",
    "Earth.SMA",
    "Earth.ECC",
    "EarthMJ2000Eq.INC",
    "EarthMJ2000Eq.RAAN",
    "EarthMJ2000Eq.AOP",
    "Earth.TA",
    "Sun.SMA",
    "Sun.ECC",
    "Heliocentric.INC",
    "Heliocentric.RAAN",
    "Heliocentric.AOP",
    "Sun.TA"
]
```

The final CSV contains mean values for each propagated quantity plus corresponding `_sigma` standard-deviation columns.

## Notes

* `main.py` uses all available CPU cores through `multiprocessing`.
* GMAT report files are deleted after they are converted to CSV to save space.
* Clones that enter Earth or forbidden regions are skipped.
* The JPL certificate path is set in `gmat_handler.py` as `ssd.jpl.nasa.gov.crt`.
* The project currently hard-codes several paths and object names inside scripts; update them before each run.
* The uploaded project archive does not include `scripts/template.script`, so a compatible GMAT template must be supplied before running `main.py`.

## Troubleshooting

### `Cannot find api_startup_file.txt`

Update `GmatInstall` in `load_gmat.py` and confirm that this file exists:

```text
<GmatInstall>/bin/api_startup_file.txt
```

### `FileNotFoundError: scripts/template.script`

Create or copy a compatible GMAT template script to:

```text
scripts/template.script
```

### JPL query errors

Confirm that the target object name is valid in JPL Horizons/SBDB and that the certificate file exists:

```text
ssd.jpl.nasa.gov.crt
```

### Missing output folders

Create the required folders listed in the setup section before running `main.py`.

## License

No license file is included in this project.
