# Height structure

## Hydrostatic and gravity equations

The atmosphere is assumed to be hydrostatically supported.

The radius at each pressure level is obtained by integrating from the surface upwards using a fourth order Runge-Kutta method [press_numerical_2007](@citep). This solves the coupled system:
```math
\frac{dr}{dp} = -\frac{1}{\rho a}
```
```math
\frac{dg}{dr} = -\frac{G M(r)}{r^2}
```
```math
a = g - r ( 2 \pi  \cos(\theta) / d )^2
```
where  $r$ is radial distance from planet centre, $p$ is pressure, $\rho$ is density, $g$ is local gravitational acceleration, $G$ is the gravitational constant, $M(r)$ is the mass enclosed within radius $r$, $\theta$ is latitude, and $d$ is the planet's axial rotation period (day length). This includes self-gravitational attraction, and makes AGNI applicable as an atmospheric structure model.

Densities are evaluated using the scheme described in [Gas densities](@ref).

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
