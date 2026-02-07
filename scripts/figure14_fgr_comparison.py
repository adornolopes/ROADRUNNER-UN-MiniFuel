#!/usr/bin/env python3
"""
Figure 14: Predicted fission gas release (FGR) for individual ROADRUNNER UN
MiniFuel specimens as a function of (a) burnup and (b) temperature, calculated
using the empirical correlations of Storms (solid symbols) and Rogozkin
(open symbols). Error bars represent propagated uncertainties arising from
representative input variability in temperature, burnup, and as-fabricated
density (Appendix A).

Reference:
    Adorno Lopes et al., "ROADRUNNER Uranium Nitride MiniFuel: Experimental
    Design, Fabrication and Pre-irradiation Baseline Characterization for
    Accelerated Burnup Testing", Nuclear Engineering and Design (2025/2026).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ---------------------------------------------------------------------------
# Data: ROADRUNNER specimen conditions (from Tables 2 & 3)
# ---------------------------------------------------------------------------
data = {
    "Target": [1]*6 + [2]*6 + [3]*6 + [4]*6 + [5]*6 + [6]*6,
    "Sample_ID": [
        "RRN01-6", "RRN01-5", "RRN01-4", "RRN01-3", "RRN01-2", "RRN01-1",
        "RRN02-6", "RRN02-5", "RRN02-4", "RRN02-3", "RRN02-2", "RRN02-1",
        "RRN03-6", "RRN03-5", "RRN03-4", "RRN03-3", "RRN03-2", "RRN03-1",
        "RRN04-6", "RRN04-5", "RRN04-4", "RRN04-3", "RRN04-2", "RRN04-1",
        "RRN05-6", "RRN05-5", "RRN05-4", "RRN05-3", "RRN05-2", "RRN05-1",
        "RRN06-6", "RRN06-5", "RRN06-4", "RRN06-3", "RRN06-2", "RRN06-1",
    ],
    "BU": [
        7.058, 7.470, 7.657, 7.839, 7.730, 7.576,
        7.338, 7.726, 7.908, 8.081, 8.020, 7.716,
        3.498, 3.681, 3.741, 3.852, 3.842, 3.740,
        5.936, 6.229, 6.399, 6.517, 6.451, 6.237,
        5.485, 5.751, 5.985, 5.990, 5.983, 5.863,
        3.757, 4.027, 4.096, 4.217, 4.163, 3.958,
    ],
    "T": [
        1190, 1190, 1181, 1183, 1163, 1182,
        877, 881, 1183, 1479, 1487, 1181,
        1184, 1187, 1181, 1172, 1187, 1181,
        875, 881, 1172, 1483, 1490, 1175,
        1191, 1183, 1183, 1179, 1182, 1175,
        877, 881, 1183, 1479, 1487, 1181,
    ],
    "TD": [
        92.89, 94.21, 93.43, 94.73, 95.21, 95.76,
        94.79, 94.36, 93.19, 92.00, 95.01, 95.49,
        89.09, 86.27, 89.02, 95.21, 95.91, 94.50,
        94.33, 95.45, 95.59, 93.48, 96.05, 95.60,
        87.30, 88.57, 88.19, 93.27, 94.07, 96.32,
        94.10, 95.43, 95.63, 95.14, 96.14, 95.01,
    ],
}

df = pd.DataFrame(data)

# ---------------------------------------------------------------------------
# FGR Correlations
# ---------------------------------------------------------------------------

def storms_fgr(T_K, BU, TD):
    """Storms FGR correlation (T in K, BU in %FIMA, TD in %TD)."""
    return 100.0 / (np.exp(0.0025 * (90.0 * TD**0.77 / BU**0.09 - T_K)) + 1.0)


def rogozkin_fgr(T_K, BU):
    """Rogozkin FGR correlation (T in °C internally, BU in %FIMA)."""
    T_C = T_K - 273.15
    return 3.05 * BU**1.92 * np.exp(-2086.0 / T_C)


# ---------------------------------------------------------------------------
# Uncertainty Propagation (first-order Taylor expansion)
# ---------------------------------------------------------------------------
# Input uncertainties
dT = 30.0      # K
dBU = 0.5      # %FIMA
dTD = 2.0      # %TD

# Step sizes for central finite differences
h_T = 1.0      # K
h_BU = 0.01    # %FIMA
h_TD = 0.1     # %TD


def propagate_storms(T, BU, TD):
    """Return FGR and propagated sigma for Storms correlation."""
    fgr = storms_fgr(T, BU, TD)
    df_dT = (storms_fgr(T + h_T, BU, TD) - storms_fgr(T - h_T, BU, TD)) / (2 * h_T)
    df_dBU = (storms_fgr(T, BU + h_BU, TD) - storms_fgr(T, BU - h_BU, TD)) / (2 * h_BU)
    df_dTD = (storms_fgr(T, BU, TD + h_TD) - storms_fgr(T, BU, TD - h_TD)) / (2 * h_TD)
    sigma = np.sqrt((df_dT * dT)**2 + (df_dBU * dBU)**2 + (df_dTD * dTD)**2)
    return fgr, sigma


def propagate_rogozkin(T, BU):
    """Return FGR and propagated sigma for Rogozkin correlation."""
    fgr = rogozkin_fgr(T, BU)
    df_dT = (rogozkin_fgr(T + h_T, BU) - rogozkin_fgr(T - h_T, BU)) / (2 * h_T)
    df_dBU = (rogozkin_fgr(T, BU + h_BU) - rogozkin_fgr(T, BU - h_BU)) / (2 * h_BU)
    sigma = np.sqrt((df_dT * dT)**2 + (df_dBU * dBU)**2)
    return fgr, sigma


# Calculate FGR predictions with uncertainties
storms_results = df.apply(lambda r: propagate_storms(r["T"], r["BU"], r["TD"]), axis=1)
rogozkin_results = df.apply(lambda r: propagate_rogozkin(r["T"], r["BU"]), axis=1)

df["FGR_Storms"] = [r[0] for r in storms_results]
df["sigma_Storms"] = [r[1] for r in storms_results]
df["FGR_Rogozkin"] = [r[0] for r in rogozkin_results]
df["sigma_Rogozkin"] = [r[1] for r in rogozkin_results]

# ---------------------------------------------------------------------------
# Plotting: Figure 14(a) FGR vs Burnup, Figure 14(b) FGR vs Temperature
# ---------------------------------------------------------------------------
target_labels = {1: "RRN01", 2: "RRN02", 3: "RRN03",
                 4: "RRN04", 5: "RRN05", 6: "RRN06"}
markers = {1: "o", 2: "s", 3: "^", 4: "D", 5: "v", 6: "H"}
colors = {1: "#1f77b4", 2: "#ff7f0e", 3: "#2ca02c",
          4: "#d62728", 5: "#9467bd", 6: "#8c564b"}

fig, axes = plt.subplots(2, 1, figsize=(9, 11), sharex=False)

for ax_idx, (x_col, x_label) in enumerate([("BU", "Burnup (%FIMA)"),
                                             ("T", "Temperature (K)")]):
    ax = axes[ax_idx]
    panel = "(a)" if ax_idx == 0 else "(b)"

    for target_id in sorted(df["Target"].unique()):
        subset = df[df["Target"] == target_id]
        c = colors[target_id]
        m = markers[target_id]
        lbl = target_labels[target_id]

        # Storms (solid)
        ax.errorbar(
            subset[x_col], subset["FGR_Storms"], yerr=subset["sigma_Storms"],
            fmt=m, color=c, markersize=8, capsize=3, capthick=1,
            label=f"Storms – {lbl}", markeredgecolor="k", markeredgewidth=0.5,
        )
        # Rogozkin (open)
        ax.errorbar(
            subset[x_col], subset["FGR_Rogozkin"], yerr=subset["sigma_Rogozkin"],
            fmt=m, color=c, markersize=8, capsize=3, capthick=1,
            markerfacecolor="none", markeredgecolor=c, markeredgewidth=1.5,
            label=f"Rogozkin – {lbl}", linestyle="none",
        )

    ax.set_ylabel("Fission Gas Release (%)", fontsize=12)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_title(f"{panel} FGR vs {x_label.split('(')[0].strip()}", fontsize=12)
    ax.grid(True, alpha=0.3)

# Shared legend below plots
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=7,
           bbox_to_anchor=(0.5, -0.02), framealpha=0.9)

plt.tight_layout(rect=[0, 0.06, 1, 1])

# Save
output_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
os.makedirs(output_dir, exist_ok=True)
fig.savefig(os.path.join(output_dir, "figure14_fgr_comparison.png"), dpi=300,
            bbox_inches="tight")
plt.show()
print("Figure 14 saved to docs/figure14_fgr_comparison.png")
