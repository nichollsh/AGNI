# Advective heating

## Deep heating

In addition to the energy transport terms, some heat production/loss may occur within a given layer of the atmosphere; e.g. from advective heat transport [cronin_advective_2016](@citep) or ohmic dissipation [matt_angular_2015](@citep). AGNI includes a parameterisation of this 'deep' heating through a Gaussian energy deposition profile in log-pressure space, centred at a user-specified centre and width:
```math
\frac{dF}{dP} = \frac{F_{\text{total}}}{\sqrt{2\pi} \sigma P} \exp\left(-\frac{(\ln P - \ln P_0)^2}{2\sigma^2}\right)
```
where $F_{\text{total}}$ is the total integrated heating rate, $P_0$ is the centre pressure of the Gaussian profile, $\sigma$ is the width in log-pressure space, and $P$ is pressure. This is a representation of energy sources/sinks from otherwise unmodelled physics. The 'deep heating' functionality was first introduced to AGNI by [Cheng An Hsieh](https://didymos65803.github.io/).

The deposited flux $F_{\text{total}}$ is set by the `physics.deep_heating.power_mode` configuration option, and can be specified in one of two ways:
* `abs`: a fixed absolute flux, $F_{\text{total}} = F_\text{abs}$ [$\mathrm{W\ m^{-2}}$].
* `rel`: a fraction of the instellation, $F_{\text{total}} = \epsilon \, F^\text{ins}$, where $\epsilon$ is the `flux_rel` parameter and $F^\text{ins}$ is defined as in [Stellar irradiation](@ref).

Setting `power_mode = "off"` disables deep heating entirely.

Two normalisation conventions are available for distributing $F_{\text{total}}$ across levels, selected via `norm_method`:
* `pressure`: the deposition rate $dF/dP$ follows the Gaussian profile directly, as written above, and is integrated over pressure from the TOA downwards.
* `mass`: the profile is instead normalised so that the flux deposited in each layer is weighted by that layer's column mass $dm = dp/g$, ensuring that $\sum_l \varepsilon_l \, dm_l = F_{\text{total}}$ exactly regardless of the local gravity profile. This is the more physically appropriate choice for a heating rate expressed per unit mass (e.g. tidal or radiogenic heating), since it is insensitive to the level spacing chosen in pressure.

If the requested deposition pressure $P_0$ lies outside the model domain (i.e. below the surface or above the TOA), the `domain` option controls how this is handled: `clamp` restricts $P_0$ to lie within $[p_\text{toa}, p_\text{boa}]$, while `boundary_flux` instead applies $F_{\text{total}}$ as a uniform flux entering through the base of the domain, bypassing the Gaussian profile altogether. This is the mechanism used to inject a prescribed heat-redistribution flux into each column of the [multicolumn mode](@ref "Multicolumn and global mode"), where the redistribution flux is applied as a lower-boundary source rather than distributed with height.

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
