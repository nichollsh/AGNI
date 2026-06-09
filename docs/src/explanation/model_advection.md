# Advective heating

## Deep heating

In addition to the energy transport terms, some heat production/loss may occur within a given layer of the atmosphere; e.g. from advective heat transport [cronin_advective_2016](@citep) or ohmic dissipation [matt_angular_2015](@citep). AGNI includes a parameterisation of this 'deep' heating through a Gaussian energy deposition profile in log-pressure space, centred at a user-specified centre and width:
```math
\frac{dF}{dP} = \frac{F_{\text{total}}}{\sqrt{2\pi} \sigma P} \exp\left(-\frac{(\ln P - \ln P_0)^2}{2\sigma^2}\right)
```
where $F_{\text{total}}$ is the total integrated heating rate, $P_0$ is the centre pressure of the Gaussian profile, $\sigma$ is the width in log-pressure space, and $P$ is pressure. This is a representation of energy sources/sinks from otherwise unmodelled physics []. The 'deep heating' functionality was first introduced to AGNI by [Cheng An Hsieh](https://didymos65803.github.io/).

## Bibliography for this page

```@bibliography
Pages = [@__FILE__]
Canonical = false
```
