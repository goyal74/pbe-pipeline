#!/usr/bin/env python3
"""
PBE Pipeline - Stage 1 (in-silico) for the lncRNA-RBP scaffold-disruption platform
==================================================================================
Map the protein-binding element (PBE) on a disease lncRNA and nominate candidate
steric-block / decoy disruptors -- entirely from public data and standard tools.

This is the "Target and interface nomination" step described in the Phase I proposal.
It is intentionally built from EXISTING resources (no new algorithm): the only custom
piece is the lightweight glue that overlays an RBP footprint onto the predicted/probed
structure and tiles candidate disruptors across the PBE.

WHAT YOU SUPPLY (all public):
  --fasta   FASTA of the lncRNA transcript (one sequence). Anchor example: MANTIS
            (lncRNA n342419). Get the transcript sequence from Ensembl/RefSeq/UCSC.
  --peaks   BED of RBP-binding intervals IN TRANSCRIPT COORDINATES (0-based start, end).
            Derive from ENCODE eCLIP for the partner RBP (e.g., SMARCA4/BRG1) or POSTAR3:
              * Download ENCODE eCLIP narrowPeak (genomic, e.g. SMARCA4 ... eCLIP).
              * Lift genomic peaks to transcript coordinates with the transcript GTF
                (e.g., bedtools intersect + a genomic->transcript map, or `gffread`/
                custom mapping). A helper note is in README.md.
            Columns: chrom(=transcript_id)  start  end  [name]  [score]
  --shape   (optional) per-nucleotide SHAPE reactivities (one value per line, or
            '<pos>\\t<reactivity>'); used as folding constraints if provided.

DEPENDENCIES (all free / standard):
  pip install ViennaRNA biopython numpy
  (RNAfold/ViennaRNA does the folding; SHAPE constraints use Deigan et al. method.)
  Optional, for the uniqueness check (run separately): NCBI BLAST+ (blastn).

OUTPUTS (written to --outdir):
  pbe_summary.tsv          PBE window(s): coordinates, mean RBP coverage, structuredness
  disruptor_candidates.tsv ranked antisense disruptors (seq, Tm, GC%, target accessibility,
                           self-structure dG)
  structure.json           sequence, MFE structure, per-base unpaired prob, RBP coverage,
                           PBE windows, top disruptors  (drives the figure)

NOTE: This nominates candidates in silico. Confirm uniqueness with blastn against the
transcriptome, and confirm the interaction experimentally (RIP/CLIP, SHAPE) in Phase II.
"""

import argparse, json, sys, os
import numpy as np

# ---------- I/O ----------
def read_fasta(path):
    seq, name = [], None
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            seq.append(line.upper().replace("T", "U"))
    return name or "lncRNA", "".join(seq)

def read_peaks(path):
    """BED in transcript coords -> list of (start, end, score)."""
    peaks = []
    for line in open(path):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split()
        s, e = int(f[1]), int(f[2])
        sc = float(f[4]) if len(f) > 4 and _isnum(f[4]) else 1.0
        peaks.append((min(s, e), max(s, e), sc))
    return peaks

def read_shape(path, n):
    react = [-999.0] * n  # -999 = no data
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        f = line.split()
        if len(f) == 1:
            # one value per line, 1-based implicit
            pass
        elif len(f) >= 2:
            i = int(f[0]) - 1
            v = float(f[1])
            if 0 <= i < n:
                react[i] = v
    # if one-per-line form was used, re-read positionally
    vals = [l.strip() for l in open(path) if l.strip()]
    if vals and len(vals[0].split()) == 1:
        for i, v in enumerate(vals):
            if i < n:
                react[i] = float(v)
    return react

def _isnum(x):
    try:
        float(x); return True
    except ValueError:
        return False

# ---------- motif scanning (use a known RBP recognition motif instead of CLIP) ----------
IUPAC = {"A": "A", "C": "C", "G": "G", "U": "U", "T": "U",
         "R": "[AG]", "Y": "[CU]", "S": "[GC]", "W": "[AU]", "K": "[GU]",
         "M": "[AC]", "B": "[CGU]", "D": "[AGU]", "H": "[ACU]", "V": "[ACG]", "N": "[ACGU]"}
def motif_to_regex(motif):
    import re
    return re.compile("".join(IUPAC.get(b.upper(), b) for b in motif.replace("T", "U")))

def motif_scan(seq, motif):
    rx = motif_to_regex(motif)
    sites = []
    for m in rx.finditer(seq):
        sites.append((m.start(), m.end()))
    return sites

# ---------- local folding accessibility (RNAplfold-style; scales to long lncRNAs) ----------
def local_accessibility(seq, W=200, L=150, cutoff=1e-4):
    """Per-base unpaired probability from local partition folding (RNA.pfl_fold)."""
    import RNA
    n = len(seq)
    W = min(W, n); L = min(L, W)
    plist = RNA.pfl_fold(seq, W, L, cutoff)
    paired = np.zeros(n)
    for e in plist:
        try:
            i, j, p = e.i, e.j, e.p
        except AttributeError:
            i, j, p = e[0], e[1], e[2]
        paired[i - 1] += p
        paired[j - 1] += p
    return np.clip(1.0 - paired, 0.0, 1.0)

def nominate_pbe_motif(sites, unpaired, flank=8, top=5, min_gap=10):
    """Each motif hit -> a PBE window; rank by accessibility (binding-competent sites)."""
    n = len(unpaired)
    cands = []
    for (s, e) in sites:
        ws, we = max(0, s - flank), min(n, e + flank)
        acc = float(unpaired[ws:we].mean())
        cands.append((ws, we, acc, s))
    cands.sort(key=lambda c: -c[2])
    picked = []
    for ws, we, acc, s in cands:
        if all(abs(s - ps) >= min_gap for _, _, _, ps in picked):
            picked.append((ws, we, acc, s))
        if len(picked) >= top:
            break
    return [(ws, we, acc) for ws, we, acc, _ in picked]

# ---------- structure ----------
def fold(seq, shape=None):
    """Return (mfe_structure, unpaired_probability[per base]) using ViennaRNA."""
    import RNA
    fc = RNA.fold_compound(seq)
    if shape is not None:
        # Deigan SHAPE pseudo-energies (standard m=1.8, b=-0.6)
        fc.sc_add_SHAPE_deigan(shape, 1.8, -0.6)
    mfe_struct, _ = fc.mfe()
    fc.pf()
    bpp = fc.bpp()  # base-pair probabilities (1-indexed matrix)
    n = len(seq)
    paired = np.zeros(n)
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            p = bpp[i][j]
            paired[i - 1] += p
            paired[j - 1] += p
    unpaired = np.clip(1.0 - paired, 0.0, 1.0)
    return mfe_struct, unpaired

# ---------- RBP coverage ----------
def coverage(peaks, n):
    cov = np.zeros(n)
    for s, e, sc in peaks:
        s = max(0, s); e = min(n, e)
        cov[s:e] += sc
    if cov.max() > 0:
        cov = cov / cov.max()
    return cov

# ---------- PBE nomination ----------
def nominate_pbe(cov, unpaired, win=30, top=3, min_gap=20):
    """A PBE is a window of high RBP coverage that is structured (low unpaired prob)."""
    n = len(cov)
    score = np.zeros(n)
    structured = 1.0 - unpaired
    for i in range(0, n - win + 1):
        c = cov[i:i + win].mean()
        s = structured[i:i + win].mean()
        # require real RBP signal; reward structured, RBP-covered windows
        score[i] = c * (0.5 + 0.5 * s) if c > 0 else 0.0
    picks = []
    order = np.argsort(score)[::-1]
    for i in order:
        if score[i] <= 0:
            break
        if all(abs(i - p) >= min_gap for p, _ in picks):
            picks.append((int(i), float(score[i])))
        if len(picks) >= top:
            break
    return [(i, i + win, sc) for i, sc in picks]

# ---------- disruptor design ----------
COMP = {"A": "U", "U": "A", "G": "C", "C": "G"}
def revcomp(s):
    return "".join(COMP.get(b, "N") for b in reversed(s))

def tm_wallace_nn(seq):
    """Approximate duplex Tm (nearest-neighbor SantaLucia, DNA-like). Modified
       chemistries (2'-MOE/LNA) raise Tm further; treat as a relative ranking."""
    s = seq.replace("U", "T")
    nn = {  # dH(kcal/mol), dS(cal/mol/K)
        "AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3), "CA": (-8.5, -22.7),
        "GT": (-8.4, -22.4), "CT": (-7.8, -21.0), "GA": (-8.2, -22.2), "CG": (-10.6, -27.2),
        "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "AC": (-8.4, -22.4), "TC": (-8.2, -22.2),
        "TG": (-8.5, -22.7), "AG": (-7.8, -21.0), "TT": (-7.9, -22.2), "CC": (-8.0, -19.9),
    }
    dH, dS = 0.2, -5.7  # initiation
    for i in range(len(s) - 1):
        d = nn.get(s[i:i + 2])
        if d:
            dH += d[0]; dS += d[1]
    C = 0.25e-6  # 0.25 uM strand
    R = 1.987
    try:
        tm = (dH * 1000) / (dS + R * np.log(C / 4)) - 273.15
    except Exception:
        tm = float("nan")
    return round(tm, 1)

def gc(seq):
    return round(100.0 * sum(b in "GC" for b in seq) / len(seq), 1)

def self_struct_dg(seq):
    try:
        import RNA
        _, dg = RNA.fold(seq)
        return round(dg, 1)
    except Exception:
        return None

def design_disruptors(seq, unpaired, pbe_windows, lens=(16, 18, 20, 22), step=2):
    cands = []
    for (ps, pe, _) in pbe_windows:
        for L in lens:
            for i in range(ps, pe - L + 1, step):
                target = seq[i:i + L]
                acc = float(unpaired[i:i + L].mean())  # higher = more accessible target
                dis = revcomp(target)                  # antisense steric-block oligo
                cands.append({
                    "pbe": f"{ps}-{pe}",
                    "target_start": i, "target_end": i + L, "len": L,
                    "target_seq_5to3": target,
                    "disruptor_antisense_5to3": dis,
                    "Tm_C_approx": tm_wallace_nn(dis),
                    "GC_pct": gc(dis),
                    "target_accessibility": round(acc, 3),
                    "self_struct_dG": self_struct_dg(dis),
                })
    # rank: accessible target, Tm in 50-70 window, low self-structure
    def keyfn(c):
        tm = c["Tm_C_approx"] or 0
        tm_pen = abs(tm - 60)
        ss = c["self_struct_dG"] if c["self_struct_dG"] is not None else 0
        return (-c["target_accessibility"], tm_pen, -(ss))  # ss less negative = better
    cands.sort(key=keyfn)
    return cands

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="lncRNA PBE mapping + disruptor nomination")
    ap.add_argument("--fasta", required=True, help="lncRNA transcript FASTA (RNA or DNA)")
    ap.add_argument("--peaks", help="RBP-binding BED in TRANSCRIPT coords")
    ap.add_argument("--motif", help="RBP recognition motif, IUPAC (e.g. Pumilio PRE = UGUANAUA); used instead of --peaks")
    ap.add_argument("--local", action="store_true", help="local folding for accessibility (auto for >1500 nt)")
    ap.add_argument("--flank", type=int, default=8, help="nt flank around a motif hit to define the PBE window")
    ap.add_argument("--shape", help="optional SHAPE reactivities")
    ap.add_argument("--win", type=int, default=30, help="PBE window size (nt; --peaks mode)")
    ap.add_argument("--top", type=int, default=3, help="number of PBE windows")
    ap.add_argument("--ncand", type=int, default=15, help="disruptor candidates to report")
    ap.add_argument("--outdir", default="pbe_out")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    name, seq = read_fasta(args.fasta)
    n = len(seq)
    if n == 0:
        sys.exit("ERROR: empty sequence")
    sites = None
    if args.motif:
        sites = motif_scan(seq, args.motif)
        if not sites:
            sys.exit(f"ERROR: motif '{args.motif}' not found in transcript")
        peaks = [(s, e, 1.0) for (s, e) in sites]
        src = f"{len(sites)} '{args.motif}' motif site(s)"
    elif args.peaks:
        peaks = read_peaks(args.peaks)
        if not peaks:
            sys.exit("ERROR: no RBP peaks parsed (need transcript-coordinate BED)")
        src = f"{len(peaks)} RBP peak(s)"
    else:
        sys.exit("ERROR: provide --peaks (CLIP) or --motif (RBP recognition motif)")
    shape = read_shape(args.shape, n) if args.shape else None

    use_local = args.local or n > 1500
    print(f"[{name}] length={n} nt, {src}"
          + (", SHAPE ON" if shape else "") + (", local folding" if use_local else ""))
    if use_local:
        mfe_struct, unpaired = "", local_accessibility(seq)
    else:
        mfe_struct, unpaired = fold(seq, shape)
    cov = coverage(peaks, n)
    if sites is not None:
        pbe = nominate_pbe_motif(sites, unpaired, flank=args.flank, top=args.top)
    else:
        pbe = nominate_pbe(cov, unpaired, win=args.win, top=args.top)
    if not pbe:
        sys.exit("ERROR: no PBE nominated")
    disruptors = design_disruptors(seq, unpaired, pbe)[: args.ncand]

    # write outputs
    with open(os.path.join(args.outdir, "pbe_summary.tsv"), "w") as fh:
        fh.write("pbe_start\tpbe_end\tmean_RBP_coverage\tmean_structuredness\tscore\n")
        for s, e, sc in pbe:
            fh.write(f"{s}\t{e}\t{cov[s:e].mean():.3f}\t{(1-unpaired[s:e].mean()):.3f}\t{sc:.3f}\n")

    cols = ["pbe", "target_start", "target_end", "len", "target_seq_5to3",
            "disruptor_antisense_5to3", "Tm_C_approx", "GC_pct",
            "target_accessibility", "self_struct_dG"]
    with open(os.path.join(args.outdir, "disruptor_candidates.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for c in disruptors:
            fh.write("\t".join(str(c[k]) for k in cols) + "\n")

    json.dump({
        "name": name, "length": n, "sequence": seq, "mfe_structure": mfe_struct,
        "unpaired_prob": [round(float(x), 3) for x in unpaired],
        "rbp_coverage": [round(float(x), 3) for x in cov],
        "pbe_windows": [{"start": s, "end": e, "score": round(sc, 3)} for s, e, sc in pbe],
        "top_disruptors": disruptors,
    }, open(os.path.join(args.outdir, "structure.json"), "w"), indent=1)

    print(f"PBE window(s): " + "; ".join(f"{s}-{e}" for s, e, _ in pbe))
    print(f"Top disruptor: {disruptors[0]['disruptor_antisense_5to3']} "
          f"(Tm~{disruptors[0]['Tm_C_approx']}C, target {disruptors[0]['target_start']}-{disruptors[0]['target_end']})")
    print(f"Wrote results to {args.outdir}/  (then run blastn for uniqueness; see README)")

if __name__ == "__main__":
    main()
