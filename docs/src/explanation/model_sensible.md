# Sensible heat flux

## Turbulent kinetic energy

Surface-atmosphere turbulent heat exchange is parameterised using a bulk aerodynamic formula with Monin-Obukhov similarity theory (see [nicholson_exchange_2006](@cite), Equation 9; and [pierrehumbert_book_2010](@cite)):
```math
F_{\text{sens}} = \rho c_p C_d u (T_{\text{surf}} - T_{\text{atm}})
```
where $\rho$ is air density, $c_p$ is specific heat capacity, $u$ is wind speed, $T_{\text{surf}}$ is surface temperature, $T_{\text{atm}}$ is the temperature of the lowest atmospheric layer, and $C_d$ is the drag coefficient:
```math
C_d = \left(\frac{\kappa}{\ln(h/z_0)}\right)^2
```
where $\kappa$ is the von Kármán constant (0.4), $h$ is the height of the lowest atmospheric layer, and $z_0$ is the surface roughness length [hogstrom_karman_1988](@citep).

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
