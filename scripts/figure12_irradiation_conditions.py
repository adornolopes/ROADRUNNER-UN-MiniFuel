#!/usr/bin/env python3
"""
Figure 12: Summary of individual ROADRUNNER sample irradiation conditions
and material properties, including burnup, average irradiation temperature
(TAVA), and fuel density (%TD).

Reference:
    Adorno Lopes et al., "ROADRUNNER Uranium Nitride MiniFuel: Experimental
    Design, Fabrication and Pre-irradiation Baseline Characterization for
    Accelerated Burnup Testing", Nuclear Engineering and Design (2025/2026).
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# ---------------------------------------------------------------------------
# Data: ROADRUNNER specimen-specific irradiation conditions (from Table 2)
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
    "Burnup_FIMA": [
        7.058, 7.470, 7.657, 7.839, 7.730, 7.576,
        7.338, 7.726, 7.908, 8.081, 8.020, 7.716,
        3.498, 3.681, 3.741, 3.852, 3.842, 3.740,
        5.936, 6.229, 6.399, 6.517, 6.451, 6.237,
        5.485, 5.751, 5.985, 5.990, 5.983, 5.863,
        3.757, 4.027, 4.096, 4.217, 4.163, 3.958,
    ],
    "Temperature_K": [
        1190, 1190, 1181, 1183, 1163, 1182,
        877, 881, 1183, 1479, 1487, 1181,
        1184, 1187, 1181, 1172, 1187, 1181,
        875, 881, 1172, 1483, 1490, 1175,
        1191, 1183, 1183, 1179, 1182, 1175,
        877, 881, 1183, 1479, 1487, 1181,
    ],
    "Density_TD": [
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
# Plot configuration
# ---------------------------------------------------------------------------
target_labels = {1: "RRN01", 2: "RRN02", 3: "RRN03",
                 4: "RRN04", 5: "RRN05", 6: "RRN06"}
markers = {1: "o", 2: "s", 3: "^", 4: "D", 5: "v", 6: "H"}

fig, ax = plt.subplots(figsize=(10, 6))

for target_id in sorted(df["Target"].unique()):
    subset = df[df["Target"] == target_id]
    sc = ax.scatter(
        subset["Burnup_FIMA"],
        subset["Temperature_K"],
        c=subset["Density_TD"],
        cmap="viridis",
        marker=markers[target_id],
        s=100,
        label=target_labels[target_id],
        edgecolor="k",
        linewidths=0.5,
        vmin=86,
        vmax=97,
    )

cbar = fig.colorbar(sc, ax=ax, pad=0.02)
cbar.set_label("As-fabricated Density (%TD)", fontsize=11)

ax.set_xlabel("Burnup (%FIMA)", fontsize=12)
ax.set_ylabel("TAVA Irradiation Temperature (K)", fontsize=12)
ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
ax.grid(True, alpha=0.3)
ax.set_xlim(2.5, 9.0)
ax.set_ylim(750, 1600)

plt.tight_layout()

# Save
output_dir = os.path.join(os.path.dirname(__file__), "..", "docs")
os.makedirs(output_dir, exist_ok=True)
fig.savefig(os.path.join(output_dir, "figure12_irradiation_conditions.png"), dpi=300)
plt.show()
print("Figure 12 saved to docs/figure12_irradiation_conditions.png")
