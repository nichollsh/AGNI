# Line-by-line radiative transfer

AGNI includes an interface to the
[Reference Forward Model](https://eodg.atm.ox.ac.uk/RFM/) (RFM) for line-by-line (LbL)
radiative transfer calculations. This provides a way to validate and benchmark the
correlated-k SOCRATES results.

!!! note
    The RFM source code is proprietary; AGNI packages it as a pre-compiled binary.

## Setup

Obtain a HITRAN-formatted `.par` line-list file. The file can contain absorption from
multiple species and can be downloaded from [hitran.org](https://hitran.org/lbl/).
Alternatively, fetch the parfiles stored on Zenodo:
```bash
./src/get_data.sh parfiles
```

## Configuration

Set the following parameters in the `[files]` and `[execution]` tables of your
configuration file:

```toml
[files]
rfm_parfile = "path/to/file.par"

[execution]
rfm_wn_min = 0      # minimum wavenumber [cm-1]
rfm_wn_max = 50000  # maximum wavenumber [cm-1]
```

The wavenumber resolution is fixed at 1 cm⁻¹.

## Output

LbL results are saved alongside all other output in the NetCDF file, as the arrays
`rfm_wn` (wavenumber axis) and `rfm_fl` (spectral flux).
