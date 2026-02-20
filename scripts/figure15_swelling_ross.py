#!/usr/bin/env python3
"""
Figure 15: Predicted Swelling for ROADRUNNER UN MiniFuel Samples
================================================================

Predicted volumetric swelling (ΔV/V %) for individual ROADRUNNER samples
as a function of burnup (%FIMA) and temperature (K), based on the
empirical correlation developed by Ross et al. [59]:

    ΔV/V (%) = 4.7e-11 * T_av^3.12 * Bu^0.83 * ρ^0.5

where:
    T_av = volume-averaged fuel temperature (K)
    Bu   = fuel burnup (at% ≈ %FIMA for UN)
    ρ    = as-fabricated fuel density (%TD)

Reference:
    Ross, S.B. and El-Genk, M.S., "Uranium nitride fuel swelling
    correlation," Journal of Nuclear Materials, 170 (1990) 169–177.

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
    "BU": [  # Burnup in %FIMA
        7.058, 7.470, 7.657, 7.839, 7.730, 7.576,
        7.338, 7.726, 7.908, 8.081, 8.020, 7.716,
        3.498, 3.681, 3.741, 3.852, 3.842, 3.740,
        5.936, 6.229, 6.399, 6.517, 6.451, 6.237,
        5.485, 5.751, 5.985, 5.990, 5.983, 5.863,
        3.757, 4.027, 4.096, 4.217, 4.163, 3.958,
    ],
    "T": [  # Temperature in K
        1190, 1190, 1181, 1183, 1163, 1182,
        877, 881, 1183, 1479, 1487, 1181,
        1184, 1187, 1181, 1172, 1187, 1181,
        875, 881, 1172, 1483, 1490, 1175,
        1191, 1183, 1183, 1179, 1182, 1175,
        877, 881, 1183, 1479, 1487, 1181,
    ],
    "TD": [  # As-fabricated density in %TD
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
# Ross Swelling Correlation
# ---------------------------------------------------------------------------

def ross_swelling(T_K, BU, TD):
    """
    Ross et al. (1990) volumetric swelling correlation for UN fuel.

    Parameters
    ----------
    T_K : float
        Volume-averaged fuel temperature (K)
    BU : float
        Fuel burnup (at% ≈ %FIMA for UN)
    TD : float
        As-fabricated fuel density (%TD)

    Returns
    -------
    float
        Volumetric swelling ΔV/V (%)
    """
    return 4.7e-11 * T_K**3.12 * BU**0.83 * TD**0.5


# ---------------------------------------------------------------------------
# Calculate swelling predictions
# ---------------------------------------------------------------------------
df["Swelling_Ross"] = df.apply(
    lambda r: ross_swelling(r["T"], r["BU"], r["TD"]), axis=1
)

# ---------------------------------------------------------------------------
# Plotting: Figure 15(a) Swelling vs Burnup, Figure 15(b) Swelling vs Temperature
# ---------------------------------------------------------------------------
target_labels = {1: "RNN01", 2: "RNN02", 3: "RNN03",
                 4: "RNN04", 5: "RNN05", 6: "RNN06"}
markers = {1: "o", 2: "s", 3: "^", 4: "D", 5: "v", 6: "H"}
colors = {1: "#1f77b4", 2: "#ff7f0e", 3: "#2ca02c",
          4: "#d62728", 5: "#9467bd", 6: "#8c564b"}

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))

# --- Panel (a): Swelling vs Burnup ---
for tgt in sorted(df["Target"].unique()):
    sub = df[df["Target"] == tgt]
    ax1.scatter(
        sub["BU"], sub["Swelling_Ross"],
        marker=markers[tgt], color=colors[tgt],
        s=120, edgecolors="black", linewidths=0.5,
        label=target_labels[tgt], zorder=3
    )

ax1.set_xlabel("Burnup (%FIMA)", fontsize=12)
ax1.set_ylabel(r"$\Delta V/V$ (%)", fontsize=12)
ax1.legend(title="Target", fontsize=9, title_fontsize=10,
           loc="upper left", frameon=True, edgecolor="grey")
ax1.grid(axis="both", linestyle="--", alpha=0.4)
ax1.set_xlim(left=3.0)
ax1.set_ylim(bottom=0)
ax1.text(0.02, 0.95, "(a)", transform=ax1.transAxes,
         fontsize=14, fontweight="bold", va="top")

# --- Panel (b): Swelling vs Temperature ---
for tgt in sorted(df["Target"].unique()):
    sub = df[df["Target"] == tgt]
    ax2.scatter(
        sub["T"], sub["Swelling_Ross"],
        marker=markers[tgt], color=colors[tgt],
        s=120, edgecolors="black", linewidths=0.5,
        label=target_labels[tgt], zorder=3
    )

ax2.set_xlabel("Temperature (K)", fontsize=12)
ax2.set_ylabel(r"$\Delta V/V$ (%)", fontsize=12)
ax2.legend(title="Target", fontsize=9, title_fontsize=10,
           loc="upper left", frameon=True, edgecolor="grey")
ax2.grid(axis="both", linestyle="--", alpha=0.4)
ax2.set_xlim(left=800)
ax2.set_ylim(bottom=0)
ax2.text(0.02, 0.95, "(b)", transform=ax2.transAxes,
         fontsize=14, fontweight="bold", va="top")

plt.tight_layout()

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
output_dir = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(output_dir, "figure15_swelling_ross.png"),
            dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(output_dir, "figure15_swelling_ross.pdf"),
            bbox_inches="tight")
plt.show()

# ---------------------------------------------------------------------------
# Print summary table
# ---------------------------------------------------------------------------
print("\n" + "="*70)
print("Figure 15: Ross Swelling Predictions for ROADRUNNER Samples")
print("="*70)
print(f"{'Sample':<12} {'T (K)':>8} {'BU (%FIMA)':>12} {'TD (%TD)':>10} {'ΔV/V (%)':>10}")
print("-"*70)
for _, row in df.iterrows():
    print(f"{row['Sample_ID']:<12} {row['T']:>8.0f} {row['BU']:>12.3f} "
          f"{row['TD']:>10.2f} {row['Swelling_Ross']:>10.2f}")
print("="*70)
