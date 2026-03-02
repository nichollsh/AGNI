# Using FastChem

[FastChem](https://newstrangeworlds.github.io/FastChem/) is a fast numerical code for
computing gas-phase thermochemical equilibrium. Given a set of elemental abundances and a
temperature-pressure profile, it solves for the mixing ratios of hundreds of gas-phase
species simultaneously, based on the principle of minimising the total Gibbs free energy
of the system.

When FastChem is enabled in AGNI, the following pipeline runs at each solver step:
1. **Surface saturation** — condensable partial pressures at the surface are adjusted
   to the saturation vapour pressure, updating the pressure grid accordingly.
2. **FastChem (gas-phase chemistry)** — the metallicity derived from the current surface
   composition is passed to FastChem alongside the T-P profile. FastChem returns
   gas-phase equilibrium mixing ratios at every model level.
3. **Rainout aloft** — AGNI's own condensation scheme then operates on the post-FastChem
   composition to handle super-saturation, cold-trapping, and latent heat transport.

Steps 2 and 3 are therefore **compatible and complementary**: FastChem handles gas-phase
speciation while AGNI handles condensation independently. The condensation scheme
(`composition.condensates`) can be used alongside FastChem.

## Installation

Run the installation script from the AGNI root directory:
```bash
./src/get_fastchem.sh
```

Then set the `FC_DIR` environment variable to the FastChem installation folder.
Adding this to your shell rc file is recommended:
```bash
export FC_DIR=/path/to/fastchem/
```

## Enabling FastChem in a configuration file

Set `physics.chemistry = true` in the `[physics]` table of your configuration file.
Provide the atmospheric composition as elemental mass abundance ratios relative to
hydrogen via `composition.metallicities`.

For example:
```toml
[physics]
chemistry = true

[composition]
p_surf = 100.0   # bar
metallicities = { H = 1.0, O = 0.5, C = 0.1 }
```

See the [Configuration reference](@ref) for all available composition parameters.

!!! warning
    The `FC_DIR` environment variable must be set before running AGNI with chemistry
    enabled, or the model will error.
