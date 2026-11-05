#!/usr/bin/env python3
"""
plot_ogt.py - genome-wide copy number plot from a Nexus OGT file.

Draws a log2 ratio panel and a B-allele frequency panel along the genome,
with chromosome boundaries and centromeres taken from a UCSC cytoband file.
Optionally overlays segments (CNVkit .cns / .call.cns, a SEG file, or any
tab-delimited file with chromosome / start / end / log2 columns) and colours
each probe by its segment's call against a gain cutoff and a loss cutoff,
which can be set separately (e.g. +0.2 for gains, -0.25 for losses).

Input (tab-delimited, header row, BAF may be blank):
    Chromosome  Position  Position  Log R Ratio  B-Allele Frequency
    1           2329601   2329721   0.237995
    1           3638584   3638771   -0.16126     0.24165

Examples
    # genome-wide, no segments
    python plot_ogt.py sample.txt --cytoband hg19_cytoBand.txt.gz

    # with segments and a symmetric +/-0.25 call cutoff, saved as PDF
    python plot_ogt.py sample.txt --cytoband hg19_cytoBand.txt.gz \\
        --segments sample.cns --cutoff 0.25 --out sample_cn.pdf

    # separate cutoffs: gain at >= 0.2, loss at <= -0.25
    python plot_ogt.py sample.txt --cytoband hg19_cytoBand.txt.gz \\
        --segments sample.cns --gain-cutoff 0.2 --loss-cutoff -0.25

    # zoom to chromosomes 9 and 17
    python plot_ogt.py sample.txt --cytoband hg19_cytoBand.txt.gz --chrom 9 17

Requires: numpy, pandas, matplotlib.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# Colours (validated for colour-vision deficiency against a light background)
C_NEUTRAL = "#8a8984"   # probes in neutral segments / no segments
C_LOSS = "#2a78d6"      # probes in loss segments
C_GAIN = "#e34948"      # probes in gain segments
C_SEGMENT = "#0b0b0b"   # segment lines
C_CENTRO = "#4a3aa7"    # centromere marker
C_CUTOFF = "#1baf7a"    # +/- cutoff lines
C_GRID = "#6e6d68"      # chromosome dividers
C_BAF = "#52514e"       # BAF points

CHROM_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]


def norm_chrom(s):
    s = pd.Series(s).astype(str).str.strip()
    s = s.str.replace(r"^chr", "", regex=True, case=False)
    return s.replace({"23": "X", "24": "Y"})


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------
def read_ogt(path):
    """Read a Nexus OGT export. Columns are taken by position, because the
    header repeats 'Position' and the BAF column is often blank."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str,
                     skip_blank_lines=True, engine="python")
    if df.shape[1] < 4:
        sys.exit(f"{path}: expected >= 4 tab-separated columns, found {df.shape[1]}")
    # Drop the header row(s): anything whose start position isn't numeric
    df = df[pd.to_numeric(df[1], errors="coerce").notna()].copy()
    out = pd.DataFrame({
        "chrom": norm_chrom(df[0]).values,
        "start": pd.to_numeric(df[1]).values,
        "end": pd.to_numeric(df[2], errors="coerce").values,
        "lrr": pd.to_numeric(df[3], errors="coerce").values,
        "baf": (pd.to_numeric(df[4], errors="coerce").values
                if df.shape[1] > 4 else np.nan),
    })
    out["end"] = out["end"].fillna(out["start"])
    out["pos"] = (out["start"] + out["end"]) / 2
    return out[out["chrom"].isin(CHROM_ORDER)].reset_index(drop=True)


def read_cytoband(path):
    cb = pd.read_csv(path, sep="\t", header=None, comment="#",
                     names=["chrom", "start", "end", "band", "stain"],
                     usecols=range(5))
    cb["chrom"] = norm_chrom(cb["chrom"]).values
    cb = cb[cb["chrom"].isin(CHROM_ORDER)]
    lengths = cb.groupby("chrom")["end"].max()
    acen = cb[cb["stain"] == "acen"].groupby("chrom").agg(
        cen_start=("start", "min"), cen_end=("end", "max"))
    if acen.empty:
        sys.exit(f"{path}: no 'acen' bands found - is this a UCSC cytoBand file?")
    return lengths, acen


def read_segments(path):
    """Accept CNVkit .cns/.call.cns, SEG, or a generic TSV; find the columns
    by name."""
    seg = pd.read_csv(path, sep="\t", comment="#")
    cols = {c.lower(): c for c in seg.columns}

    def pick(options, what):
        for o in options:
            if o in cols:
                return cols[o]
        sys.exit(f"{path}: couldn't find a {what} column (tried {options}); "
                 f"columns are {list(seg.columns)}")

    c = pick(["chromosome", "chrom", "chr"], "chromosome")
    s = pick(["start", "loc.start", "startpos", "start_position"], "start")
    e = pick(["end", "loc.end", "endpos", "end_position"], "end")
    v = pick(["log2", "seg.mean", "segment_mean", "mean", "log2ratio", "value"],
             "log2 value")
    out = pd.DataFrame({"chrom": norm_chrom(seg[c]).values,
                        "start": seg[s].values, "end": seg[e].values,
                        "log2": pd.to_numeric(seg[v], errors="coerce").values})
    return out[out["chrom"].isin(CHROM_ORDER)].dropna().reset_index(drop=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def assign_segments(probes, segs):
    """Segment log2 for each probe (NaN if outside every segment)."""
    val = np.full(len(probes), np.nan)
    for ch, idx in probes.groupby("chrom").indices.items():
        s = segs[segs["chrom"] == ch].sort_values("start")
        if s.empty:
            continue
        pos = probes["pos"].values[idx]
        k = np.searchsorted(s["start"].values, pos, side="right") - 1
        ok = (k >= 0)
        ok[ok] = pos[ok] <= s["end"].values[k[ok]]
        v = np.full(len(pos), np.nan)
        v[ok] = s["log2"].values[k[ok]]
        val[idx] = v
    return val


def qc_metrics(probes):
    """Overall spread and derivative log ratio spread (probe-to-probe noise,
    robust to real copy number changes)."""
    lrr = probes["lrr"].dropna()
    d = probes.dropna(subset=["lrr"]).groupby("chrom", sort=False)["lrr"].diff().dropna()
    return {
        "probes": len(lrr),
        "SD": lrr.std(),
        "MAD": 1.4826 * np.median(np.abs(lrr - lrr.median())),
        "DLRS": d.std() / np.sqrt(2) if len(d) > 1 else np.nan,
    }


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ogt", help="Nexus OGT file")
    ap.add_argument("--cytoband", required=True,
                    help="UCSC cytoBand.txt(.gz) for the same genome build")
    ap.add_argument("--segments", help="optional segment file (.cns, .seg, TSV)")
    ap.add_argument("--cutoff", type=float, default=None,
                    help="symmetric call cutoff, sets both gain and loss to +/-cutoff "
                         "(default 0.25 when --segments is given)")
    ap.add_argument("--gain-cutoff", type=float, default=None,
                    help="segments with log2 >= this are gains (overrides --cutoff)")
    ap.add_argument("--loss-cutoff", type=float, default=None,
                    help="segments with log2 <= this are losses, e.g. -0.25 "
                         "(a positive value is taken as negative; overrides --cutoff)")
    ap.add_argument("--chrom", nargs="+", help="plot only these chromosomes")
    ap.add_argument("--sample", help="title [file name]")
    ap.add_argument("--out", help="output .png/.pdf/.svg [<sample>.png]")
    ap.add_argument("--ylim", type=float, default=2.0,
                    help="log2 axis limit; points beyond it are drawn as triangles [2]")
    ap.add_argument("--no-baf", action="store_true", help="omit the BAF panel")
    ap.add_argument("--no-qc", action="store_true", help="omit QC metrics from the title")
    ap.add_argument("--width", type=float, default=None, help="figure width in inches")
    ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()

    sample = a.sample or os.path.basename(a.ogt).split(".")[0]
    out = a.out or f"{sample}.png"

    probes = read_ogt(a.ogt)
    lengths, acen = read_cytoband(a.cytoband)

    chroms = [c for c in CHROM_ORDER if c in set(probes["chrom"]) and c in lengths.index]
    if a.chrom:
        want = list(norm_chrom(a.chrom))
        missing = [c for c in want if c not in chroms]
        if missing:
            sys.exit(f"Not in data or cytoband: {', '.join(missing)}")
        chroms = [c for c in CHROM_ORDER if c in want]
    probes = probes[probes["chrom"].isin(chroms)].reset_index(drop=True)
    if probes.empty:
        sys.exit("No probes to plot.")

    # Genome coordinates: chromosomes laid end to end at true length
    offsets, x = {}, 0
    for c in chroms:
        offsets[c] = x
        x += lengths[c]
    genome_len = x
    probes["gx"] = probes["pos"] + probes["chrom"].map(offsets)

    # Cutoffs: --gain-cutoff / --loss-cutoff override the symmetric --cutoff
    base = a.cutoff
    if base is None and a.segments:
        base = 0.25
    gain_cut = a.gain_cutoff if a.gain_cutoff is not None else (
        abs(base) if base is not None else None)
    loss_cut = a.loss_cutoff if a.loss_cutoff is not None else (
        -abs(base) if base is not None else None)
    if loss_cut is not None:
        loss_cut = -abs(loss_cut)
    if gain_cut is not None and gain_cut < 0:
        sys.exit("--gain-cutoff should be positive")
    cutoff = gain_cut is not None or loss_cut is not None

    # Segment calls
    segs = None
    if a.segments:
        segs = read_segments(a.segments)
        segs = segs[segs["chrom"].isin(chroms)]
        probes["seg"] = assign_segments(probes, segs)
        is_gain = probes["seg"] >= gain_cut if gain_cut is not None else False
        is_loss = probes["seg"] <= loss_cut if loss_cut is not None else False
        state = np.where(is_gain, "gain", np.where(is_loss, "loss", "neutral"))
    else:
        state = np.full(len(probes), "neutral")
    probes["state"] = state

    # Figure
    single = len(chroms) == 1
    width = a.width or (8 if single else max(12, min(20, 4 + 0.6 * len(chroms))))
    show_baf = not a.no_baf and probes["baf"].notna().any()
    heights = [4, 1.2] if show_baf else [4]
    fig, axes = plt.subplots(len(heights), 1, figsize=(width, sum(heights) + 1.2),
                             sharex=True, gridspec_kw={"height_ratios": heights,
                                                       "hspace": 0.08})
    axes = np.atleast_1d(axes)
    ax = axes[0]
    lim = a.ylim
    ms = 3 if len(probes) > 20000 else 6

    # Probes, clamped to +/-lim
    colours = {"neutral": C_NEUTRAL, "loss": C_LOSS, "gain": C_GAIN}
    lrr = probes["lrr"]
    inside = lrr.abs() <= lim
    for st in ("neutral", "loss", "gain"):            # neutral underneath
        m = (probes["state"] == st) & inside
        ax.scatter(probes.loc[m, "gx"], lrr[m], s=ms, c=colours[st],
                   alpha=0.45 if st == "neutral" else 0.7,
                   linewidths=0, rasterized=True, zorder=2)
        for sign, mk in ((1, "^"), (-1, "v")):
            mm = (probes["state"] == st) & (lrr * sign > lim)
            if mm.any():
                ax.scatter(probes.loc[mm, "gx"], np.full(mm.sum(), sign * lim * 0.97),
                           marker=mk, s=18, c=colours[st], linewidths=0,
                           zorder=3, clip_on=False)

    # Segments: dark line with a thin white halo so it reads over the points
    if segs is not None and not segs.empty:
        for _, s in segs.iterrows():
            y = float(np.clip(s["log2"], -lim, lim))
            x0, x1 = s["start"] + offsets[s["chrom"]], s["end"] + offsets[s["chrom"]]
            ax.plot([x0, x1], [y, y], color="white", lw=3.2, solid_capstyle="butt", zorder=4)
            ax.plot([x0, x1], [y, y], color=C_SEGMENT, lw=1.8, solid_capstyle="butt", zorder=5)

    # Cutoff and zero lines
    ax.axhline(0, color="#52514e", lw=0.8, zorder=1)
    for y in (gain_cut, loss_cut):
        if y is not None:
            ax.axhline(y, color=C_CUTOFF, lw=1, ls="--", zorder=1)

    # Chromosome dividers, centromeres, labels
    for axx in axes:
        for c in chroms[1:]:
            axx.axvline(offsets[c], color=C_GRID, lw=1.0, zorder=0)
        for c in chroms:
            if c in acen.index:
                cx = (acen.loc[c, "cen_start"] + acen.loc[c, "cen_end"]) / 2 + offsets[c]
                axx.axvline(cx, color=C_CENTRO, lw=0.8, ls="--", alpha=0.3, zorder=1)

    ax.set_ylim(-lim, lim)
    ax.set_ylabel("Log2 copy ratio")
    ax.set_xlim(0, genome_len)

    # BAF panel
    if show_baf:
        bx = axes[1]
        b = probes.dropna(subset=["baf"])
        bx.scatter(b["gx"], b["baf"], s=ms, c=C_BAF, alpha=0.6, linewidths=0,
                   rasterized=True, zorder=2)
        bx.axhline(0.5, color=C_CUTOFF, lw=1, ls="--", zorder=1)
        bx.set_ylim(-0.02, 1.02)
        bx.set_yticks([0, 0.5, 1])
        bx.set_ylabel("BAF")

    # X axis: chromosome names genome-wide, Mb when zoomed to one chromosome
    bottom = axes[-1]
    if single:
        c = chroms[0]
        ticks = np.arange(0, lengths[c] + 1, 10e6 if lengths[c] > 60e6 else 5e6)
        bottom.set_xticks(ticks)
        bottom.set_xticklabels([f"{t / 1e6:.0f}" for t in ticks])
        bottom.set_xlabel(f"Chromosome {c} position (Mb)")
    else:
        mids = [offsets[c] + lengths[c] / 2 for c in chroms]
        bottom.set_xticks(mids)
        bottom.set_xticklabels(chroms, fontsize=8)
        bottom.tick_params(axis="x", length=0)
        bottom.set_xlabel("Chromosome")

    for axx in axes:
        for side in ("top", "right"):
            axx.spines[side].set_visible(False)
        axx.tick_params(labelsize=8)

    # Title with QC metrics
    title = sample
    if not a.no_qc:
        q = qc_metrics(probes)
        title += (f"\nprobes: {q['probes']:,}   SD: {q['SD']:.3f}   "
                  f"MAD: {q['MAD']:.3f}   DLRS: {q['DLRS']:.3f}")
        if gain_cut is not None and loss_cut is not None and gain_cut == -loss_cut:
            title += f"   cutoff: ±{gain_cut:g}"
        else:
            if gain_cut is not None:
                title += f"   gain: ≥{gain_cut:g}"
            if loss_cut is not None:
                title += f"   loss: ≤{loss_cut:g}"
    ax.set_title(title, fontsize=10, loc="left")

    # Legend
    handles = [Line2D([], [], marker="o", ls="", color=C_NEUTRAL, label="Neutral")]
    if segs is not None:
        handles += [Line2D([], [], marker="o", ls="", color=C_LOSS, label="Loss"),
                    Line2D([], [], marker="o", ls="", color=C_GAIN, label="Gain"),
                    Line2D([], [], color=C_SEGMENT, lw=1.8, label="Segment")]
    if cutoff:
        handles.append(Line2D([], [], color=C_CUTOFF, lw=1, ls="--", label="Cutoff"))
    handles.append(Line2D([], [], color=C_CENTRO, lw=0.8, ls="--", alpha=0.3,
                          label="Centromere"))
    ax.legend(handles=handles, loc="upper right", fontsize=7, frameon=False,
              ncol=len(handles), bbox_to_anchor=(1, 1.13))

    fig.savefig(out, dpi=a.dpi, bbox_inches="tight")
    print(f"Wrote {out}  ({len(probes):,} probes, {len(chroms)} chromosomes)")


if __name__ == "__main__":
    main()
