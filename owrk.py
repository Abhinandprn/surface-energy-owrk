"""Surface free energy of a solid from contact angles (OWRK method).

Owens, Wendt, Rabel and Kaelble split the surface energy into a dispersive (d)
and a polar (p) part. For each test liquid L on the solid S:

    gamma_L (1 + cos theta) / (2 sqrt(gamma_L^d)) = sqrt(gamma_S^p) * sqrt(gamma_L^p / gamma_L^d) + sqrt(gamma_S^d)

Plotting y = gamma_L (1 + cos theta) / (2 sqrt(gamma_L^d)) against
x = sqrt(gamma_L^p / gamma_L^d) gives a straight line:
    slope     = sqrt(gamma_S^p)
    intercept = sqrt(gamma_S^d)

Usage:
    python owrk.py examples/example_angles.csv
    python owrk.py examples/example_angles.csv --plot owrk_plot.png
"""

import argparse
import sys

import numpy as np
import pandas as pd

# Surface tension of common test liquids in mN/m at 20 °C: (total, dispersive, polar)
# Values as commonly used with the OWRK method (Stroem et al. 1987 / KRUESS tables).
LIQUIDS = {
    "water": (72.8, 21.8, 51.0),
    "diiodomethane": (50.8, 50.8, 0.0),
    "ethylene glycol": (47.7, 26.4, 21.3),
    "formamide": (58.0, 39.0, 19.0),
    "glycerol": (63.4, 37.0, 26.4),
}


def owrk(liquids, angles_deg):
    """Return (total, dispersive, polar) surface energy of the solid in mN/m,
    plus the x and y values used for the fit and the R^2 of the line."""
    if len(liquids) < 2:
        raise ValueError("OWRK needs at least two test liquids.")

    x, y = [], []
    for name, theta in zip(liquids, angles_deg):
        key = name.strip().lower()
        if key not in LIQUIDS:
            raise ValueError(f"Unknown liquid '{name}'. Known: {', '.join(LIQUIDS)}")
        g_total, g_d, g_p = LIQUIDS[key]
        x.append(np.sqrt(g_p / g_d))
        y.append(g_total * (1 + np.cos(np.radians(theta))) / (2 * np.sqrt(g_d)))

    x, y = np.array(x), np.array(y)
    slope, intercept = np.polyfit(x, y, 1)

    gamma_p = max(slope, 0.0) ** 2
    gamma_d = max(intercept, 0.0) ** 2

    y_fit = slope * x + intercept
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return gamma_d + gamma_p, gamma_d, gamma_p, x, y, slope, intercept, r2


def main():
    parser = argparse.ArgumentParser(description="OWRK surface free energy from contact angles")
    parser.add_argument("csv", help="CSV with columns: sample, liquid, angle_deg")
    parser.add_argument("--plot", help="save an OWRK plot to this PNG file")
    args = parser.parse_args()

    data = pd.read_csv(args.csv)
    required = {"sample", "liquid", "angle_deg"}
    if not required.issubset(data.columns):
        sys.exit(f"CSV needs the columns {sorted(required)}")

    # average repeated drops of the same liquid on the same sample
    mean = data.groupby(["sample", "liquid"], sort=False)["angle_deg"].agg(["mean", "std", "count"]).reset_index()

    results = []
    fits = {}
    for sample, group in mean.groupby("sample", sort=False):
        total, d, p, x, y, slope, intercept, r2 = owrk(group["liquid"], group["mean"])
        results.append({"sample": sample, "total_mN_m": total, "dispersive_mN_m": d,
                        "polar_mN_m": p, "polar_fraction_%": 100 * p / total, "R2": r2})
        fits[sample] = (x, y, slope, intercept, list(group["liquid"]))

    print("\nMean contact angles (deg):")
    print(mean.round(1).to_string(index=False))
    print("\nSurface free energy (OWRK):")
    print(pd.DataFrame(results).round(2).to_string(index=False))

    if args.plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4.5))
        for sample, (x, y, slope, intercept, names) in fits.items():
            pts = ax.plot(x, y, "o", label=f"{sample}")[0]
            xs = np.linspace(0, max(x) * 1.1, 50)
            ax.plot(xs, slope * xs + intercept, "-", color=pts.get_color(), alpha=0.7)
        ax.set_xlabel(r"$\sqrt{\gamma_L^p / \gamma_L^d}$")
        ax.set_ylabel(r"$\gamma_L (1+\cos\theta) / (2\sqrt{\gamma_L^d})$  [$\sqrt{mN/m}$]")
        ax.set_title("OWRK plot: slope = $\\sqrt{\\gamma_S^p}$, intercept = $\\sqrt{\\gamma_S^d}$")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.plot, dpi=150)
        print(f"\nPlot saved to {args.plot}")


if __name__ == "__main__":
    main()
