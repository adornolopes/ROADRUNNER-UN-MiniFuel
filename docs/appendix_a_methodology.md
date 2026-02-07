# Appendix A: Uncertainty Propagation Methodology

## Overview

The sensitivity of FGR predictions to input parameter variations was quantified using first-order uncertainty propagation based on Taylor series expansion.

## Formulation

For a correlation of the form *f(T, BU, ρ)*, the propagated variance is:

$$\sigma^2(\text{FGR}) = \left(\frac{\partial f}{\partial T}\right)^2 \Delta T^2 + \left(\frac{\partial f}{\partial BU}\right)^2 \Delta BU^2 + \left(\frac{\partial f}{\partial \rho}\right)^2 \Delta \rho^2$$

where the partial derivatives represent the local sensitivity of the FGR prediction to each input variable, and ΔT, ΔBU, and Δρ are the assumed input uncertainties.

## Numerical Implementation

Partial derivatives are evaluated via central finite differences:

$$\frac{\partial f}{\partial x} \approx \frac{f(x + h) - f(x - h)}{2h}$$

with step sizes:
- *h* = 1 K for temperature
- *h* = 0.01 %FIMA for burnup
- *h* = 0.1 %TD for density

## Input Uncertainties

| Parameter | Symbol | Value | Basis |
|-----------|--------|-------|-------|
| Temperature | ΔT | ±30 K | Thermal design margin + SiC thermometry resolution |
| Burnup | ΔBU | ±0.5 %FIMA | Neutronic calculation uncertainty |
| Density | Δρ | ±2 %TD | As-fabricated density variability within each target group |

## Variance Decomposition

The variance contribution from each input variable is calculated as:

$$\text{Variance contribution (\%)} = \frac{\sigma_i^2}{\sigma_{\text{total}}^2} \times 100$$

This decomposition identifies which input parameters dominate the uncertainty in the FGR prediction for each sample and temperature regime.

## FGR Correlations

### Storms Correlation
Function of temperature *T* (K), burnup *BU* (%FIMA), and density *ρ* (%TD):

$$\text{FGR} = \frac{100}{\exp\left[0.0025 \times \left(90 \times \frac{\rho^{0.77}}{BU^{0.09}} - T\right)\right] + 1}$$

### Rogozkin Correlation
Function of temperature *T* (°C) and burnup *BU* (%FIMA):

$$\text{FGR} = 3.05 \times BU^{1.92} \times \exp\left(\frac{-2086}{T}\right)$$

## Key Findings

- At **low temperatures (~873 K)**: propagated uncertainty is small (~±0.2% for Storms, ~±0.7% for Rogozkin)
- At **intermediate temperatures (~1,173 K)**: uncertainty increases to ~±0.5% and ~±1.7%, respectively
- At **high temperatures (~1,473 K)**: Storms uncertainty remains bounded at ~±0.9%, while Rogozkin grows to ~±2.9%
- **Inter-model discrepancies substantially exceed propagated uncertainty bounds**, demonstrating structural rather than input-driven divergence

### Variance Decomposition by Regime

**Storms Correlation**: Density contributes ~56% of total variance, temperature ~29%, burnup ~15% (stable across regimes due to explicit density dependence).

**Rogozkin Correlation**: Temperature and burnup contribute approximately equally at low temperature. Burnup dominates at intermediate (~82%) and high temperatures (~93%).
