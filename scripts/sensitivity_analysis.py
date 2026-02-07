#!/usr/bin/env python3
"""
Sensitivity Analysis: First-order uncertainty propagation and variance
decomposition for the Storms and Rogozkin FGR correlations applied to
the ROADRUNNER UN MiniFuel irradiation matrix.

Methodology (Appendix A):
    σ²(FGR) = (∂f/∂T)² · ΔT² + (∂f/∂BU)² · ΔBU² + (∂f/∂ρ)² · Δρ²

    Partial derivatives evaluated via central finite differences:
        ∂f/∂x ≈ [f(x + h) − f(x − h)] / 2h

Correlations:
    Storms:   FGR = 100 / {exp[0.0025 × (90 × ρ^0.77 / BU^0.09 − T)] + 1}
    Rogozkin: FGR = 3.05 × BU^1.92 × exp(−2086 / T_°C)

Input uncertainties:
    ΔT  = ±30 K    (thermal design margin + SiC thermometry resolution)
    ΔBU = ±0.5 %FIMA (neutronic calculation uncertainty)
    Δρ  = ±2 %TD   (as-fabricated density variability)

Reference:
    Adorno Lopes et al., "ROADRUNNER Uranium Nitride MiniFuel: Experimental
    Design, Fabrication and Pre-irradiation Baseline Characterization for
    Accelerated Burnup Testing", Nuclear Engineering and Design (2025/2026).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ===================================================================
# Data: ROADRUNNER specimen conditions
# ===================================================================
data = {
    "Sample_ID": [
        "RRN01-6", "RRN01-5", "RRN01-4", "RRN01-3", "RRN01-2", "RRN01-1",
        "RRN02-6", "RRN02-5", "RRN02-4", "RRN02-3", "RRN02-2", "RRN02-1",
        "RRN03-6", "RRN03-5", "RRN03-4", "RRN03-3", "RRN03-2", "RRN03-1",
        "RRN04-6", "RRN04-5", "RRN04-4", "RRN04-3", "RRN04-2", "RRN04-1",
        "RRN05-6", "RRN05-5", "RRN05-4", "RRN05-3", "RRN05-2", "RRN05-1",
        "RRN06-6", "RRN06-5", "RRN06-4", "RRN06-3", "RRN06-2", "RRN06-1",
    ],
    "Target": [1]*6 + [2]*6 + [3]*6 + [4]*6 + [5]*6 + [6]*6,
    "T_K": [
        1190, 1190, 1181, 1183, 1163, 1182,
        877, 881, 1183, 1479, 1487, 1181,
        1184, 1187, 1181, 1172, 1187, 1181,
        875, 881, 1172, 1483, 1490, 1175,
        1191, 1183, 1183, 1179, 1182, 1175,
        877, 881, 1183, 1479, 1487, 1181,
    ],
    "BU_FIMA": [
        7.058, 7.470, 7.657, 7.839, 7.730, 7.576,
        7.338, 7.726, 7.908, 8.081, 8.020, 7.716,
        3.498, 3.681, 3.741, 3.852, 3.842, 3.740,
        5.936, 6.229, 6.399, 6.517, 6.451, 6.237,
        5.485, 5.751, 5.985, 5.990, 5.983, 5.863,
        3.757, 4.027, 4.096, 4.217, 4.163, 3.958,
    ],
    "TD_pct": [
        92.89, 94.21, 93.43, 94.73, 95.21, 95.76,
        94.79, 94.36, 93.19, 92.00, 95.01, 95.49,
        89.09, 86.27, 89.02, 95.21, 95.91, 94.50,
        94.33, 95.45, 95.59, 93.48, 96.05, 95.60,
        87.30, 88.57, 88.19, 93.27, 94.07, 96.32,
        94.10, 95.43, 95.63, 95.14, 96.14, 95.01,
    ],
}

df = pd.DataFrame(data)

# ===================================================================
# FGR Correlations
# ===================================================================

def storms_fgr(T_K, BU, TD):
    """Storms FGR correlation. T in K, BU in %FIMA, TD in %TD."""
    return 100.0 / (np.exp(0.0025 * (90.0 * TD**0.77 / BU**0.09 - T_K)) + 1.0)


def rogozkin_fgr(T_K, BU):
    """Rogozkin FGR correlation. T in K (converted to °C internally)."""
    T_C = T_K - 273.15
    return 3.05 * BU**1.92 * np.exp(-2086.0 / T_C)


# ===================================================================
# Uncertainty Propagation Parameters
# ===================================================================
DELTA_T = 30.0    # K
DELTA_BU = 0.5    # %FIMA
DELTA_TD = 2.0    # %TD

# Central finite difference step sizes
H_T = 1.0         # K
H_BU = 0.01       # %FIMA
H_TD = 0.1        # %TD


# ===================================================================
# Compute FGR, uncertainties, and variance decomposition
# ===================================================================
results = []

for _, row in df.iterrows():
    T = row["T_K"]
    BU = row["BU_FIMA"]
    TD = row["TD_pct"]

    # --- Storms ---
    fgr_s = storms_fgr(T, BU, TD)

    dfdT_s = (storms_fgr(T + H_T, BU, TD) - storms_fgr(T - H_T, BU, TD)) / (2 * H_T)
    dfdBU_s = (storms_fgr(T, BU + H_BU, TD) - storms_fgr(T, BU - H_BU, TD)) / (2 * H_BU)
    dfdTD_s = (storms_fgr(T, BU, TD + H_TD) - storms_fgr(T, BU, TD - H_TD)) / (2 * H_TD)

    var_T_s = (dfdT_s * DELTA_T) ** 2
    var_BU_s = (dfdBU_s * DELTA_BU) ** 2
    var_TD_s = (dfdTD_s * DELTA_TD) ** 2
    var_total_s = var_T_s + var_BU_s + var_TD_s
    sigma_s = np.sqrt(var_total_s)

    # --- Rogozkin ---
    fgr_r = rogozkin_fgr(T, BU)

    dfdT_r = (rogozkin_fgr(T + H_T, BU) - rogozkin_fgr(T - H_T, BU)) / (2 * H_T)
    dfdBU_r = (rogozkin_fgr(T, BU + H_BU) - rogozkin_fgr(T, BU - H_BU)) / (2 * H_BU)

    var_T_r = (dfdT_r * DELTA_T) ** 2
    var_BU_r = (dfdBU_r * DELTA_BU) ** 2
    var_total_r = var_T_r + var_BU_r
    sigma_r = np.sqrt(var_total_r)

    results.append({
        "Sample_ID": row["Sample_ID"],
        "Target": row["Target"],
        "T_K": T,
        "BU_FIMA": BU,
        "TD_pct": TD,
        # Storms
        "FGR_Storms": fgr_s,
        "sigma_Storms": sigma_s,
        "Storms_var_T_pct": 100 * var_T_s / var_total_s if var_total_s > 0 else 0,
        "Storms_var_BU_pct": 100 * var_BU_s / var_total_s if var_total_s > 0 else 0,
        "Storms_var_TD_pct": 100 * var_TD_s / var_total_s if var_total_s > 0 else 0,
        # Rogozkin
        "FGR_Rogozkin": fgr_r,
        "sigma_Rogozkin": sigma_r,
        "Rogozkin_var_T_pct": 100 * var_T_r / var_total_r if var_total_r > 0 else 0,
        "Rogozkin_var_BU_pct": 100 * var_BU_r / var_total_r if var_total_r > 0 else 0,
        # Discrepancy
        "Abs_Discrepancy": abs(fgr_s - fgr_r),
    })

res_df = pd.DataFrame(results)

# ===================================================================
# Print summary table (Table A1 equivalent)
# ===================================================================
print("=" * 110)
print("Table A1: FGR Predictions, Propagated Uncertainties, and Model Discrepancy")
print(f"Input uncertainties: ΔT = ±{DELTA_T} K, ΔBU = ±{DELTA_BU} %FIMA, Δρ = ±{DELTA_TD} %TD")
print("=" * 110)
print(f"{'Sample':<12} {'T(K)':>7} {'BU':>7} {'ρ(%TD)':>8} "
      f"{'FGR_S(%)':>9} {'σ_S(%)':>7} {'FGR_R(%)':>9} {'σ_R(%)':>7} {'|Δ|(%)':>7}")
print("-" * 110)

for _, r in res_df.iterrows():
    print(f"{r['Sample_ID']:<12} {r['T_K']:7.1f} {r['BU_FIMA']:7.2f} {r['TD_pct']:8.2f} "
          f"{r['FGR_Storms']:9.2f} {r['sigma_Storms']:7.2f} "
          f"{r['FGR_Rogozkin']:9.2f} {r['sigma_Rogozkin']:7.2f} "
          f"{r['Abs_Discrepancy']:7.2f}")

# ===================================================================
# Variance Decomposition Summary
# ===================================================================
print("\n" + "=" * 80)
print("Variance Decomposition Summary (% of total variance)")
print("=" * 80)

# Group by approximate temperature regime
def temp_regime(T):
    if T < 1000:
        return "Low (~873 K)"
    elif T < 1350:
        return "Intermediate (~1,173 K)"
    else:
        return "High (~1,473 K)"

res_df["Regime"] = res_df["T_K"].apply(temp_regime)

for regime in ["Low (~873 K)", "Intermediate (~1,173 K)", "High (~1,473 K)"]:
    sub = res_df[res_df["Regime"] == regime]
    if len(sub) == 0:
        continue
    print(f"\n{regime} (n={len(sub)} specimens):")
    print(f"  Storms:   T = {sub['Storms_var_T_pct'].mean():.1f}%, "
          f"BU = {sub['Storms_var_BU_pct'].mean():.1f}%, "
          f"ρ = {sub['Storms_var_TD_pct'].mean():.1f}%")
    print(f"  Rogozkin: T = {sub['Rogozkin_var_T_pct'].mean():.1f}%, "
          f"BU = {sub['Rogozkin_var_BU_pct'].mean():.1f}%")
    print(f"  Avg |Δ| = {sub['Abs_Discrepancy'].mean():.2f}%")

# ===================================================================
# Plot: Variance Decomposition Bar Chart (Figure A1)
# ===================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# --- Storms variance decomposition ---
ax = axes[0]
x = np.arange(len(res_df))
width = 0.8
ax.bar(x, res_df["Storms_var_T_pct"], width, label="Temperature", color="#2196F3")
ax.bar(x, res_df["Storms_var_BU_pct"], width,
       bottom=res_df["Storms_var_T_pct"], label="Burnup", color="#FF9800")
ax.bar(x, res_df["Storms_var_TD_pct"], width,
       bottom=res_df["Storms_var_T_pct"] + res_df["Storms_var_BU_pct"],
       label="Density", color="#4CAF50")
ax.set_xlabel("Sample Index")
ax.set_ylabel("Variance Contribution (%)")
ax.set_title("(a) Storms Correlation")
ax.legend(fontsize=9)
ax.set_ylim(0, 105)

# --- Rogozkin variance decomposition ---
ax = axes[1]
ax.bar(x, res_df["Rogozkin_var_T_pct"], width, label="Temperature", color="#2196F3")
ax.bar(x, res_df["Rogozkin_var_BU_pct"], width,
       bottom=res_df["Rogozkin_var_T_pct"], label="Burnup", color="#FF9800")
ax.set_xlabel("Sample Index")
ax.set_ylabel("Variance Contribution (%)")
ax.set_title("(b) Rogozkin Correlation")
ax.legend(fontsize=9)
ax.set_ylim(0, 105)

plt.tight_layout()

# Save
output_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
os.makedirs(output_dir, exist_ok=True)
fig.savefig(os.path.join(output_dir, "figure_A1_variance_decomposition.png"), dpi=300)
plt.show()

# Also save results to CSV
csv_path = os.path.join(os.path.dirname(__file__), "..", "data", "sensitivity_results_computed.csv")
res_df.to_csv(csv_path, index=False, float_format="%.4f")

print(f"\nFigure A1 saved to docs/figure_A1_variance_decomposition.png")
print(f"Results saved to data/sensitivity_results_computed.csv")
