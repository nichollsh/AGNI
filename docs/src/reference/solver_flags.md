# Solver and output flags

During execution, AGNI emits flags in the log output to indicate what operations were
performed and whether any issues arose. These flags appear on the solver status line,
separated by hyphens (`-`).

## Gas species flags

These flags appear in the composition table printed at startup, alongside each gas species:

| Flag        | Meaning |
|:------------|:--------|
| `EOS_[XX]`  | Using the `[XX]` equation of state (e.g. `EOS_IDEAL`, `EOS_AQUA`). |
| `NO_OPACITY`| No opacity data available; gas can still contribute to thermodynamics. |
| `NO_THERMO` | No thermodynamic data available; gas is treated as a diatomic ideal gas. |
| `COND`      | This gas is allowed to condense. |

## Solver step flags

These flags appear at each solver iteration:

| Flag          | Meaning |
|:--------------|:--------|
| `Cs` / `Cf`   | Chemistry and condensation either (s)ucceeded or (f)ailed. |
| `Gg`          | Radiative transfer performed with double-grey scheme. |
| `M` / `Mr`    | Convective fluxes are being modulated for stability. |
| `C2` / `C4`   | Central finite-difference scheme used (at 2nd or 4th order). |
| `F2` / `F4`   | Forward finite-difference scheme used (at 2nd or 4th order). |
| `Ls`          | A linesearch method was applied to scale the update step. |
| `P`           | Step was forcibly extrapolated because the solver is not making progress. |
| `U`           | The atmosphere has become gravitationally unbound. |
