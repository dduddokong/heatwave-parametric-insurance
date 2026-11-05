"""
Correlated heat risk -> per-capita payout -> TVaR capital loading.

ILLUSTRATIVE, SYNTHETIC PARAMETERS. NOT a Korean calibration.

This script reproduces the *structure* of "Model A" in the working paper
(static optimal contract under correlated heat risk):

    Z_i = sqrt(rho) * S + sqrt(1 - rho) * eps_i,   S, eps_i ~ iid N(0, 1)
    X_i = h(Z_i) in [0, M]                         (seasonal income loss)
    A(S) = E[I(X) | S]                             (per-capita payout, N -> inf)
    Cost(I) = (1 + lambda0) E[A] + kappa * TVaR_q(A)
    m~(x) = 1 + lambda0 + kappa * pi(x) / (1 - q),  pi(x) = P(S > z_q | X = x)

Because Z is standard normal for every rho, the marginal loss distribution is
held fixed; only the dependence across workers changes with rho.

What is computed
  1. Monte Carlo over the common heat shock S (seed fixed), with A(S) evaluated
     by Gauss-Hermite quadrature over the idiosyncratic shock -> E[A], sd(A),
     TVaR_q(A) and the effective full-insurance loading lambda_full(rho).
  2. A finite-pool check (N workers simulated explicitly) showing that pooling
     removes idiosyncratic risk but not the common-shock tail.
  3. The extinction threshold kappa*(rho) from the closed form (eq. 7 of the
     draft), via numerical integration, plus lambda_full evaluated at kappa*.

Assumption that is NOT in the draft: the loss map h. The draft's theory
section does not publish h, so we use a capped linear ramp
    h(z) = M * clip((z - z0) / c, 0, 1)
("lost work-days x daily wage, capped at the season maximum") and calibrate
(z0, c) to two numbers the draft does report for its example:
E[X] = 0.075 and kappa*(0) = 1.032. Other values therefore differ from the
draft's tables; only the qualitative pattern is meant to carry over.

Dependencies: numpy (matplotlib optional, for the figure).
"""


import math
from pathlib import Path

import numpy as np

# ----------------------------------------------------------------------------
# Example parameters (same as the draft's illustrative example)
# ----------------------------------------------------------------------------
SEED = 20261001
W = 1.0          # seasonal income
M = 0.3          # maximum seasonal loss (bounded, M < W)
GAMMA = 3.0      # CRRA risk aversion
LAMBDA0 = 0.1    # base (correlation-free) loading
Q = 0.9          # tail level: TVaR over the worst 10% of years
KAPPA = 0.1      # cost of tail capital used in the main table
TARGET_EX = 0.075
TARGET_KSTAR0 = 1.032

RHOS = [0.0, 0.2, 0.4, 0.6, 0.8, 0.84, 0.95]
N_YEARS = 200_000          # Monte Carlo draws of the common shock S
N_GH = 201                 # Gauss-Hermite nodes for E[. | S]
FINITE_N = 1_000           # workers in the finite-pool check
FINITE_YEARS = 20_000

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"              # small summary outputs
FIG_DIR = ROOT / "reports" / "figures"

_erf = np.frompyfunc(math.erf, 1, 1)


def Phi(x):
    x = np.asarray(x, dtype=float)
    return 0.5 * (1.0 + _erf(x / math.sqrt(2.0)).astype(float))


Z_Q = 1.2815515655446004   # Phi^{-1}(0.9)

# Fine grid over Z ~ N(0,1) for numerical integration (marginal quantities)
ZG = np.linspace(-8.0, 8.0, 16001)
FG = np.exp(-ZG**2 / 2.0)
FG /= FG.sum()


def h(z, z0, c):
    return M * np.clip((z - z0) / c, 0.0, 1.0)


def u_prime(c):
    return c ** (-GAMMA)


# ----------------------------------------------------------------------------
# Calibration of h (the only non-draft ingredient)
# ----------------------------------------------------------------------------
def _solve_c(z0):
    """c such that E[X] = TARGET_EX for a given z0 (E[X] decreases in c)."""
    lo, hi = 1e-4, 100.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if (h(ZG, z0, mid) * FG).sum() > TARGET_EX:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def calibrate():
    # kappa*(0) = u'(w-M) / E[u'(w-X)] - (1 + lambda0)   (Prop. 3(c))
    target_eu = u_prime(W - M) / (1.0 + LAMBDA0 + TARGET_KSTAR0)

    def gap(z0):
        x = h(ZG, z0, _solve_c(z0))
        return (u_prime(W - x) * FG).sum() - target_eu

    lo, hi = -1.0, 0.6     # E[u'] is increasing in z0 on this bracket
    g_lo = gap(lo)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        g_mid = gap(mid)
        if g_mid * g_lo > 0:
            lo, g_lo = mid, g_mid
        else:
            hi = mid
    z0 = 0.5 * (lo + hi)
    return z0, _solve_c(z0)


# ----------------------------------------------------------------------------
# Model quantities
# ----------------------------------------------------------------------------
def pi_of_z(z, rho):
    """P(S > z_q | Z = z). S | Z=z ~ N(sqrt(rho) z, 1 - rho)."""
    if rho == 0.0:
        return np.full_like(z, 1.0 - Q)
    return 1.0 - Phi((Z_Q - math.sqrt(rho) * z) / math.sqrt(1.0 - rho))


def lambda_full_closed(rho, kappa, z0, c):
    """Eq. (4): lambda_full = lambda0 + kappa * E[X pi(X)] / ((1-q) E[X])."""
    x = h(ZG, z0, c)
    pi = pi_of_z(ZG, rho)
    return LAMBDA0 + kappa * (x * pi * FG).sum() / ((1.0 - Q) * (x * FG).sum())


def kappa_star(rho, z0, c, n_t=600):
    """Eq. (7): extinction threshold of tail-capital cost."""
    x = h(ZG, z0, c)
    up = u_prime(W - x)
    eu = (up * FG).sum()
    pi = pi_of_z(ZG, rho)
    best = 0.0
    for t in np.linspace(0.0, M * 0.9999, n_t):
        ind = x > t
        p_tail = (ind * FG).sum()
        if p_tail < 1e-14:
            continue
        num = (up * ind * FG).sum() - (1.0 + LAMBDA0) * eu * p_tail
        den = eu * (ind * pi * FG).sum() / (1.0 - Q)
        best = max(best, num / den)
    return best


def payout_given_S(s, rho, z0, c, indemnity):
    """A(S) = E[I(h(Z)) | S] via Gauss-Hermite over the idiosyncratic shock."""
    nodes, weights = np.polynomial.hermite_e.hermegauss(N_GH)
    weights = weights / weights.sum()
    z = math.sqrt(rho) * s[:, None] + math.sqrt(1.0 - rho) * nodes[None, :]
    return (indemnity(h(z, z0, c)) * weights[None, :]).sum(axis=1)


def tvar(a, q=Q):
    var = np.quantile(a, q)
    return a[a >= var].mean(), var


def full_insurance(x):
    return x


def layer_contract(x):
    # fixed layer used in the draft's multi-period section: min(max(x-0.1,0),0.1)
    return np.minimum(np.maximum(x - 0.1, 0.0), 0.1)


# ----------------------------------------------------------------------------
def main():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)

    z0, c = calibrate()
    ex = (h(ZG, z0, c) * FG).sum()
    p_loss = (h(ZG, z0, c) > 0).astype(float) @ FG
    p_cap = (h(ZG, z0, c) >= M).astype(float) @ FG
    print("ILLUSTRATIVE, SYNTHETIC PARAMETERS -- not a Korean calibration")
    print(f"seed={SEED}  w={W} M={M} gamma={GAMMA} lambda0={LAMBDA0} q={Q} kappa={KAPPA}")
    print(f"calibrated h: z0={z0:.4f}, c={c:.4f}  ->  E[X]={ex:.4f}, "
          f"P(X>0)={p_loss:.3f}, P(X=M)={p_cap:.3f}")
    print()

    s = rng.standard_normal(N_YEARS)     # common heat shock, one per year
    rows = []
    payouts = {}
    for rho in RHOS:
        a = payout_given_S(s, rho, z0, c, full_insurance)
        t, v = tvar(a)
        cost = (1.0 + LAMBDA0) * a.mean() + KAPPA * t
        lam_mc = cost / a.mean() - 1.0
        lam_cf = lambda_full_closed(rho, KAPPA, z0, c)
        ks = kappa_star(rho, z0, c)
        lam_at_ks = lambda_full_closed(rho, ks, z0, c)
        a_layer = payout_given_S(s, rho, z0, c, layer_contract)
        t_layer, _ = tvar(a_layer)
        rows.append(dict(rho=rho, EA=a.mean(), sdA=a.std(), VaR=v, TVaR=t,
                         tvar_mult=t / a.mean(), cap_load=KAPPA * t,
                         lam_mc=lam_mc, lam_cf=lam_cf, kstar=ks,
                         lam_at_ks=lam_at_ks,
                         layer_tvar_mult=t_layer / a_layer.mean()))
        payouts[rho] = a

    hdr = ("rho", "E[A]", "sd(A)", "VaR90", "TVaR90", "TVaR/E[A]",
           "k*TVaR", "lam_full MC", "lam_full eq4", "kappa*", "lam_full@k*",
           "layer TVaR/E")
    fmt = "{:>5} {:>7} {:>7} {:>7} {:>7} {:>9} {:>7} {:>11} {:>12} {:>7} {:>11} {:>12}"
    lines = [fmt.format(*hdr)]
    for r in rows:
        lines.append(fmt.format(
            f"{r['rho']:.2f}", f"{r['EA']:.4f}", f"{r['sdA']:.4f}",
            f"{r['VaR']:.4f}", f"{r['TVaR']:.4f}", f"{r['tvar_mult']:.2f}",
            f"{r['cap_load']:.4f}", f"{r['lam_mc']:.3f}", f"{r['lam_cf']:.3f}",
            f"{r['kstar']:.3f}", f"{r['lam_at_ks']:.2f}",
            f"{r['layer_tvar_mult']:.2f}"))
    table = "\n".join(lines)
    print("Full insurance I(x)=x, per-capita annual payout A(S), kappa=0.1:")
    print(table)
    print()

    # Lemma 2(i): mean price coefficient does not depend on rho
    for rho in (0.0, 0.6, 0.95):
        mbar = 1 + LAMBDA0 + KAPPA * (pi_of_z(ZG, rho) * FG).sum() / (1 - Q)
        print(f"Lemma 2(i) check rho={rho:.2f}: E[m~(X)] = {mbar:.4f} "
              f"(theory {1 + LAMBDA0 + KAPPA:.4f})")
    print()

    # Finite pool: pooling removes idiosyncratic risk, not the common shock
    print(f"Finite-pool check: N={FINITE_N} workers, {FINITE_YEARS} simulated years")
    fin_lines = []
    for rho in (0.0, 0.6, 0.84):
        s_f = rng.standard_normal(FINITE_YEARS)
        avg = np.empty(FINITE_YEARS)
        for k in range(0, FINITE_YEARS, 2000):
            e = rng.standard_normal((min(2000, FINITE_YEARS - k), FINITE_N))
            zz = math.sqrt(rho) * s_f[k:k + 2000, None] + math.sqrt(1 - rho) * e
            avg[k:k + 2000] = h(zz, z0, c).mean(axis=1)
        t_f, _ = tvar(avg)
        ln = (f"  rho={rho:.2f}: mean={avg.mean():.4f}  sd={avg.std():.5f}  "
              f"TVaR90={t_f:.4f}  TVaR/mean={t_f / avg.mean():.3f}")
        print(ln)
        fin_lines.append(ln)
    print()

    # Mapping note: corr of losses X_i, X_j vs latent rho
    print("Loss correlation corr(X_i, X_j) implied by latent rho (two workers, same S):")
    corr_lines = []
    s2 = rng.standard_normal(400_000)
    e2 = rng.standard_normal((400_000, 2))
    for rho in (0.6, 0.84, 0.9, 0.95):
        zz = math.sqrt(rho) * s2[:, None] + math.sqrt(1 - rho) * e2
        xx = h(zz, z0, c)
        r_x = np.corrcoef(xx[:, 0], xx[:, 1])[0, 1]
        ln = f"  latent rho={rho:.2f} -> corr(X_i,X_j)={r_x:.3f}"
        print(ln)
        corr_lines.append(ln)

    # Save text output
    with open(OUT / "results.txt", "w", encoding="utf-8") as fh:
        fh.write("ILLUSTRATIVE, SYNTHETIC PARAMETERS -- not a Korean calibration\n")
        fh.write(f"seed={SEED} w={W} M={M} gamma={GAMMA} lambda0={LAMBDA0} "
                 f"q={Q} kappa={KAPPA}\n")
        fh.write(f"h(z)=M*clip((z-z0)/c,0,1), z0={z0:.4f}, c={c:.4f}, E[X]={ex:.4f}\n\n")
        fh.write(table + "\n\n")
        fh.write("Finite pool:\n" + "\n".join(fin_lines) + "\n\n")
        fh.write("Loss correlation vs latent rho:\n" + "\n".join(corr_lines) + "\n")
    with open(OUT / "results.csv", "w", encoding="utf-8") as fh:
        keys = list(rows[0].keys())
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(f"{r[k]:.6g}" for k in keys) + "\n")

    # Figure (optional)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping figure")
        return

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    colors = {0.0: "#9aa0a6", 0.6: "#4c78a8", 0.84: "#d1495b"}
    bins = np.linspace(0, 0.3, 61)
    for rho, col in colors.items():
        ax[0].hist(payouts[rho], bins=bins, density=True, histtype="step",
                   lw=1.8, color=col, label=f"rho = {rho:.2f}")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("per-capita annual payout A(S)")
    ax[0].set_ylabel("density (log)")
    ax[0].set_title("Same marginal loss, fatter aggregate tail")
    ax[0].legend(frameon=False)

    rr = [r["rho"] for r in rows]
    ax[1].plot(rr, [r["lam_mc"] for r in rows], "o-", color="#4c78a8",
               label="Monte Carlo (TVaR of simulated A)")
    ax[1].plot(rr, [r["lam_cf"] for r in rows], "x--", color="#222222",
               label="closed form, eq. (4)")
    ax[1].set_xlabel("latent correlation rho")
    ax[1].set_ylabel("effective loading lambda_full")
    ax[1].set_title("Full-insurance loading, kappa = 0.1")
    ax[1].legend(frameon=False)

    ax[2].plot(rr, [r["kstar"] for r in rows], "o-", color="#d1495b")
    ax[2].set_xlabel("latent correlation rho")
    ax[2].set_ylabel("kappa*(rho)")
    ax[2].set_title("Extinction threshold of tail-capital cost")
    for a_ in ax[1:]:
        a_.axvline(0.84, color="#999999", lw=0.8, ls=":")
    ax[2].text(0.83, ax[2].get_ylim()[1] * 0.9, "0.84 (ERA5, prelim.) ", ha="right",
               fontsize=8, color="#666666")
    fig.suptitle("Illustrative, synthetic parameters (w=1, M=0.3, CRRA 3, "
                 "lambda0=0.1, q=0.9) - not a Korean calibration", fontsize=10)
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "tvar_loading_vs_correlation.png", dpi=150)
    print("\nsaved reports/figures/tvar_loading_vs_correlation.png")


if __name__ == "__main__":
    main()
