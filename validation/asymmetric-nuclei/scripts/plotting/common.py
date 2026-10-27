"""Shared loader / helpers for the asymmetric-nuclei 2p2h validation plots (see ../SPEC.md).

Environment overrides (used for development on mock data):
  PLOTS_DATA : directory holding flat_*.npz and tot*.txt          (default S/data)
  PLOTS_FIGS : output directory for PNGs, summary.json, INDEX.md  (default S/plots/figs)
"""
import os
import re
import json
import glob
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- paths
S = Path(os.environ.get("PLOTS_S", Path(__file__).resolve().parents[2]))  # validation/asymmetric-nuclei
DATA = Path(os.environ.get("PLOTS_DATA", S / "data"))
AUX = S / "data" / "aux"
HARNESS = S / "harness"
RUNS = S / "runs"
FIGS = Path(os.environ.get("PLOTS_FIGS", S / "plots" / "figs"))
CACHE = S / "plots" / "cache"
FIGS.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------- physics tables
NUCS = ["C12", "O16", "O18", "Ar40", "Ca40", "Ca44", "Pb208"]
A_OF = {"C12": 12, "O16": 16, "O18": 18, "Ar40": 40, "Ca40": 40, "Ca44": 44, "Pb208": 208}
Z_OF = {"C12": 6, "O16": 8, "O18": 8, "Ar40": 18, "Ca40": 20, "Ca44": 20, "Pb208": 82}
N_OF = {"C12": 6, "O16": 8, "O18": 10, "Ar40": 22, "Ca40": 20, "Ca44": 24, "Pb208": 126}
BEAMS = ["nu", "nubar"]
E_MAX = 2000.0
BIN_E = 50.0

# slot -> final state (N1',N2') and initial pair
SLOT_FINAL = {"nu": ["pp", "np", "pn"], "nubar": ["nn", "pn", "np"]}
SLOT_INIT = {"nu": ["np", "nn", "nn"], "nubar": ["np", "pp", "pp"]}
BEAM_LABEL = {"nu": r"$\nu$", "nubar": r"$\bar\nu$"}
BEAM_LS = {"nu": "-", "nubar": "--"}
BEAM_MARK = {"nu": "o", "nubar": "s"}

# --------------------------------------------------------------------------- colours (colour-blind safe)
# nuclei: Paul Tol "bright"-like; slots: Okabe-Ito; pieces: Okabe-Ito remainder
NUC_COL = {"C12": "#4477AA", "O16": "#66CCEE", "O18": "#228833", "Ar40": "#CCBB44",
           "Ca40": "#EE6677", "Ca44": "#AA3377", "Pb208": "#222222"}
SLOT_COL = ["#D55E00", "#0072B2", "#009E73"]
TOTAL_COL = "#000000"
PIECE_COL = {"rho": "#E69F00", "pi": "#56B4E9", "int": "#CC79A7"}
PIECE_LABEL = {"rho": r"$\rho$-transverse", "pi": r"$\pi$-longitudinal", "int": r"$\rho\pi$ interference"}
RD_COL = {0: "#000000", 1: "#0072B2", 2: "#D55E00"}
RD_LABEL = {0: r"Re$\Sigma_\Delta$=0", 1: r"Re$\Sigma_\Delta$=$-$50$\,\rho/\rho_0$ MeV",
            2: r"Re$\Sigma_\Delta$=+30$\,\rho/\rho_0$ MeV"}


def slot_label(beam, k, with_beam=False):
    s = SLOT_FINAL[beam][k]
    return f"{'nu' if beam == 'nu' else 'nubar'} {s}" if with_beam else s


def mirror_label(k):
    """slot k label showing both beams, e.g. 'pp | nn' (nu | nubar)."""
    return f"{SLOT_FINAL['nu'][k]} ($\\nu$) / {SLOT_FINAL['nubar'][k]} ($\\bar\\nu$)"


# --------------------------------------------------------------------------- matplotlib defaults
plt.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 150, "font.size": 10, "axes.titlesize": 11,
    "axes.labelsize": 10.5, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "legend.fontsize": 9, "legend.frameon": False, "axes.grid": True, "grid.alpha": 0.25,
    "lines.linewidth": 1.6, "figure.constrained_layout.use": False, "mathtext.default": "regular",
    "axes.spines.top": False, "axes.spines.right": False,
})

# --------------------------------------------------------------------------- bookkeeping
MISSING = []


def note_missing(what):
    if what not in MISSING:
        MISSING.append(what)
        print(f"  [missing input] {what}", flush=True)


def update_summary(d):
    """Deep-merge dict d into FIGS/summary.json."""
    p = FIGS / "summary.json"
    try:
        cur = json.loads(p.read_text())
    except Exception:
        cur = {}

    def merge(a, b):
        for k, v in b.items():
            if isinstance(v, dict) and isinstance(a.get(k), dict):
                merge(a[k], v)
            else:
                a[k] = v
        return a
    merge(cur, _jsonable(d))
    p.write_text(json.dumps(cur, indent=1, sort_keys=True))


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if np.isfinite(v) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def save_fig(fig, name, caption, extra=None):
    """Save FIGS/<name>.png at dpi 150 and register the caption (+ optional numbers) in summary.json."""
    out = FIGS / f"{name}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    d = {"captions": {name: caption}}
    if extra:
        d[name] = extra
    update_summary(d)
    print(f"  wrote {out}", flush=True)


def no_data(ax, text="data missing"):
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", fontsize=10,
            color="#666666", style="italic")
    ax.set_xticks([])
    ax.set_yticks([])


def stamp_missing(fig, items):
    """Footnote listing samples that were not (yet) available."""
    if items:
        fig.text(0.005, 0.003, "not yet available: " + ", ".join(items), fontsize=9, color="#a00000",
                 ha="left", va="bottom")


# --------------------------------------------------------------------------- flat samples
class Sample:
    """One flat-flux event sample (npz).  Arrays are loaded lazily and cached as float64."""

    def __init__(self, nuc, beam, rd=0):
        self.nuc, self.beam, self.rd = nuc, beam, rd
        self.tag = f"flat_{nuc}_{beam}_RD{rd}"
        self.path = DATA / f"{self.tag}.npz"
        self._z = np.load(self.path)
        c = np.asarray(self._z["constants"], dtype=float)
        self.nr_events, self.last_E, self.mlep, self.C3, self.Q, self.n_throws = c[:6]
        self.E_lo = self.mlep + self.Q
        self.A, self.Z, self.N = A_OF[nuc], Z_OF[nuc], N_OF[nuc]
        # sigma = sum(w) * fac
        self.fac = self.C3 / 1000.0 * 2.0 * np.pi / self.n_throws
        self._c = {}

    # raw arrays --------------------------------------------------------------
    def a(self, name):
        if name not in self._c:
            self._c[name] = np.asarray(self._z[name], dtype=np.float64)
        return self._c[name]

    def free(self, *names):
        for n in names:
            self._c.pop(n, None)

    @property
    def n(self):
        return len(self.a("w0"))

    def w(self, k=None):
        """weight of slot k (0..2) or the slot-summed weight (k=None)."""
        if k is None:
            if "wtot" not in self._c:
                self._c["wtot"] = self.a("w0") + self.a("w1") + self.a("w2")
            return self._c["wtot"]
        return self.a(f"w{k}")

    def wp(self, piece, k=None):
        if k is None:
            return sum(self.a(f"w{i}_{piece}") for i in range(3))
        return self.a(f"w{k}_{piece}")

    # derived variables -------------------------------------------------------
    @property
    def E(self):
        return self.a("energy_neutrino")

    @property
    def q0(self):
        return self.a("q0")

    @property
    def qz(self):
        return self.a("qz")

    @property
    def omega(self):
        return self.a("q0") - self.Q

    @property
    def Q2(self):
        """Q^2 = |q|^2 - q0^2 in GeV^2"""
        return (self.a("qz") ** 2 - self.a("q0") ** 2) * 1e-6

    @property
    def Emu(self):
        return self.E - self.q0

    @property
    def cth(self):
        return np.cos(self.a("theta_lepton"))

    # cross-section helpers -----------------------------------------------------
    def sigma_tot(self, k=None):
        return float(self.w(k).sum() * self.fac)

    def sigma_tot_err(self, k=None):
        return float(np.sqrt((self.w(k) ** 2).sum()) * self.fac)

    def ebins(self):
        return ebins(self.E_lo)

    def sigma_E(self, wts, edges=None, mask=None):
        """sigma(E_nu) per bin (1e-38 cm^2 per nucleus) and MC error."""
        edges = self.ebins() if edges is None else edges
        b = Binner(self.E, edges, mask)
        s, e = b.sums(wts)
        f = self.fac * (E_MAX - self.E_lo) / np.diff(edges)
        return s * f, e * f

    def dsig(self, x, wts, edges, mask=None):
        """flux-averaged d sigma/dx (1e-38 cm^2 per nucleus per unit x) and MC error."""
        b = Binner(x, edges, mask)
        s, e = b.sums(wts)
        f = self.fac / np.diff(edges)
        return s * f, e * f


class Binner:
    def __init__(self, x, edges, mask=None):
        edges = np.asarray(edges, dtype=float)
        self.nb = len(edges) - 1
        idx = np.searchsorted(edges, x, side="right") - 1
        # right-most edge inclusive
        idx[x == edges[-1]] = self.nb - 1
        ok = (idx >= 0) & (idx < self.nb)
        if mask is not None:
            ok &= mask
        self.ok = ok
        self.idx = idx[ok]

    def sums(self, w):
        w = w[self.ok]
        s = np.bincount(self.idx, weights=w, minlength=self.nb)
        e = np.sqrt(np.bincount(self.idx, weights=w * w, minlength=self.nb))
        return s, e

    def count(self, extra=None):
        if extra is None:
            return np.bincount(self.idx, minlength=self.nb)
        return np.bincount(self.idx, weights=extra[self.ok].astype(float), minlength=self.nb)


def ebins(E_lo):
    """50 MeV bins on the multiples-of-50 grid, first edge above E_lo, up to 2000 MeV."""
    first = BIN_E * np.ceil(E_lo / BIN_E)
    return np.arange(first, E_MAX + 1e-9, BIN_E)


def centres(edges):
    return 0.5 * (edges[1:] + edges[:-1])


def load_sample(nuc, beam, rd=0, quiet=False):
    p = DATA / f"flat_{nuc}_{beam}_RD{rd}.npz"
    if not p.exists():
        note_missing(str(p.name))
        return None
    try:
        return Sample(nuc, beam, rd)
    except Exception as ex:  # partially written file etc.
        note_missing(f"{p.name} (unreadable: {ex})")
        return None


def available(nucs=NUCS, beams=BEAMS, rd=0):
    out = []
    for nuc in nucs:
        for beam in beams:
            if (DATA / f"flat_{nuc}_{beam}_RD{rd}.npz").exists():
                out.append((nuc, beam))
    return out


def missing_list(nucs=NUCS, beams=BEAMS, rd=0):
    return [f"{n}_{b}_RD{rd}" for n in nucs for b in beams
            if not (DATA / f"flat_{n}_{b}_RD{rd}.npz").exists()]


# --------------------------------------------------------------------------- total-mode files
def load_totals(kind, beam, rd, nuc="C12"):
    """kind in {'tot3p3h','tot2p2h'}. returns dict E -> (mean, sem or nan, n) from all files present."""
    out = {}
    pat = str(DATA / f"{kind}_{nuc}_{beam}_RD{rd}_E*.txt")
    for fn in sorted(glob.glob(pat)):
        m = re.search(r"_E(\d+)\.txt$", fn)
        if not m:
            continue
        rows = []
        for line in open(fn):
            sp = line.split()
            if len(sp) >= 2:
                try:
                    rows.append((float(sp[0]), float(sp[1])))
                except ValueError:
                    pass
        if not rows:
            continue
        v = np.array([r[1] for r in rows])
        E = float(m.group(1))
        sem = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan
        out[E] = (float(v.mean()), float(sem), len(v))
    return out


# --------------------------------------------------------------------------- auxiliary densities / luminosities
def load_dens():
    import pandas as pd
    p = AUX / "densities.csv"
    if not p.exists():
        note_missing("aux/densities.csv")
        return None
    return pd.read_csv(p)


def pair_luminosities(dens=None):
    """L_ab = 4 pi int r^2 rho_a rho_b dr (fm^-3 ... dimensionless counts-like), per nucleus."""
    dens = load_dens() if dens is None else dens
    if dens is None:
        return {}
    out = {}
    for nuc, d in dens.groupby("nucleus"):
        r = d["r_fm"].to_numpy()
        rp, rn = d["rho_p_fm3"].to_numpy(), d["rho_n_fm3"].to_numpy()
        f = 4 * np.pi * r ** 2
        tr = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
        out[nuc] = {"np": float(tr(f * rp * rn, r)), "nn": float(tr(f * rn * rn, r)),
                    "pp": float(tr(f * rp * rp, r)), "tot": float(tr(f * (rp + rn) ** 2, r)),
                    "Np": float(tr(f * rp, r)), "Nn": float(tr(f * rn, r))}
    return out


def slot_pair(beam, k):
    """initial-pair type (np, nn or pp) feeding slot k."""
    return SLOT_INIT[beam][k]


# --------------------------------------------------------------------------- misc numeric helpers
def pos_only(y):
    y = np.array(y, dtype=float)
    y[~(y > 0)] = np.nan
    return y


def ratio_err(a, ea, b, eb):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = a / b
        er = np.abs(r) * np.sqrt((ea / a) ** 2 + (eb / b) ** 2)
    return r, er


def step_band(ax, edges, y, e, color, ls="-", label=None, alpha=0.25, lw=1.6, band=True, **kw):
    """step line over bin edges with optional +-e band, using bin centres (plain line) for readability."""
    xc = centres(edges)
    ln, = ax.plot(xc, y, color=color, ls=ls, lw=lw, label=label, **kw)
    if band and e is not None:
        ax.fill_between(xc, y - e, y + e, color=color, alpha=alpha, lw=0)
    return ln


def finish(name):
    if MISSING:
        print(f"[{name}] missing inputs: {len(MISSING)}")
    else:
        print(f"[{name}] all inputs present")


def common_edges(edge_list):
    """50 MeV bin edges common to all samples (they all sit on the multiples-of-50 grid): start at the largest first edge."""
    edge_list = [e for e in edge_list if e is not None]
    if not edge_list:
        return None
    first = max(e[0] for e in edge_list)
    return np.arange(first, E_MAX + 1e-9, BIN_E)


def sigma_window(s, k=None, lo=150.0, hi=2000.0):
    """flux-averaged sigma over the common window [lo,hi] (flat flux): (value, MC error) in 1e-38 cm^2.
    All nuclei share the same window, so Q-dependent E_lo differences do not bias cross-nucleus comparisons."""
    m = (s.E >= lo) & (s.E < hi)
    w = s.w(k)[m]
    f = s.fac * (E_MAX - s.E_lo) / (hi - lo)
    return float(w.sum() * f), float(np.sqrt((w * w).sum()) * f)
