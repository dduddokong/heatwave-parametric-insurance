"""Build reports/figures/hero.png from results/results.csv (no re-simulation).

Run from the repo root:  python reports/make_figures.py
Numbers come only from results/results.csv, written by src/heatwave_tvar_sim.py.
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "results" / "results.csv"
OUT = ROOT / "reports" / "figures" / "hero.png"

# Inline copy of the shared portfolio style (pstyle.py)
COLORS = {
    "primary": "#1F4E79",
    "accent": "#E07A1F",
    "muted": "#9AA5B1",
    "light": "#D6E2EE",
    "text": "#1F2933",
    "note": "#6B7785",
}
plt.rcParams.update({
    "figure.figsize": (10, 4.2),
    "figure.dpi": 160,
    "savefig.dpi": 160,
    "savefig.bbox": "tight",
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.labelcolor": COLORS["text"],
    "axes.edgecolor": "#C9D1D9",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#EDF0F3",
    "grid.linewidth": 0.8,
    "xtick.color": COLORS["text"],
    "ytick.color": COLORS["text"],
    "legend.frameon": False,
    "axes.axisbelow": True,
})


def load():
    with open(CSV, encoding="utf-8") as fh:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(fh)]
    return rows


def main():
    rows = load()
    rho = [r["rho"] for r in rows]
    mult = [r["tvar_mult"] for r in rows]
    kstar = [r["kstar"] for r in rows]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.2))
    fig.suptitle("Same mean payout, but at ρ = 0.84 the tail is 3.45x the mean "
                 "and κ* falls from 1.03 to 0.13",
                 x=0.01, ha="left", fontsize=13, fontweight="bold",
                 color=COLORS["text"])

    # Left: TVaR90 / E[A] per rho
    labels = [f"{x:.2f}" for x in rho]
    cols = [COLORS["accent"] if abs(x - 0.84) < 1e-9 else COLORS["primary"]
            for x in rho]
    bars = a1.bar(labels, mult, color=cols, width=0.65)
    for b, m in zip(bars, mult):
        a1.text(b.get_x() + b.get_width() / 2, m + 0.06, f"{m:.2f}x",
                ha="center", va="bottom", fontsize=9, color=COLORS["text"])
    a1.axhline(1.0, color=COLORS["muted"], lw=1, ls="--")
    a1.set_ylim(0, 4.3)
    a1.set_xlabel("latent correlation ρ")
    a1.set_ylabel("TVaR₉₀ / E[payout]")
    a1.set_title("1-in-10-year tail of per-capita payout", fontsize=11)
    a1.grid(axis="x", visible=False)

    # Right: kappa*(rho)
    a2.plot(rho, kstar, "o-", color=COLORS["primary"], lw=2)
    i84 = rho.index(0.84)
    a2.plot([0.84], [kstar[i84]], "o", color=COLORS["accent"], ms=9, zorder=3)
    a2.annotate(f"κ*(0) = {kstar[0]:.3f}", (0, kstar[0]), xytext=(0.12, 0.98),
                fontsize=9.5, color=COLORS["text"],
                arrowprops=dict(arrowstyle="-", color=COLORS["muted"]))
    a2.annotate(f"κ*(0.84) = {kstar[i84]:.3f}\n(0.84 = prelim. ERA5 heat-day-count\n"
                "correlation, not latent ρ)",
                (0.84, kstar[i84]), xytext=(0.42, 0.42), fontsize=9.5,
                color=COLORS["accent"],
                arrowprops=dict(arrowstyle="-", color=COLORS["accent"]))
    a2.set_ylim(0, 1.15)
    a2.set_xlim(-0.05, 1.0)
    a2.set_xlabel("latent correlation ρ")
    a2.set_ylabel("κ*  (tail-capital cost at which\ninsurance demand vanishes)")
    a2.set_title("Market-extinction threshold κ*(ρ)", fontsize=11)

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.text(0.01, -0.02,
             "Source: results/results.csv (src/heatwave_tvar_sim.py, seed 20261001). "
             "Illustrative synthetic parameters (w=1, M=0.3, CRRA 3, λ0=0.1, q=0.9, κ=0.1); "
             "loss map h assumed. Not a Korean calibration.",
             fontsize=8, color=COLORS["note"], ha="left", va="top")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    print(f"saved {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
