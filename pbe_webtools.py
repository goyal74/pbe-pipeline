#!/usr/bin/env python3
"""
pbe_webtools.py - helper utilities for the PBE Pipeline web app (Epigenuity LLC)
================================================================================
Two capabilities used by app.py:

1. fetch_by_name(symbol, species)  -> (label, sequence)
   Resolve a gene / lncRNA *name* (e.g. "MALAT1", "NEAT1", "NORAD") to its
   transcript sequence automatically, so the user never has to paste a sequence.
   Tries Ensembl REST first, then NCBI E-utilities as a fallback. Runs on the
   server that hosts the app (normal outbound HTTPS to public bioinformatics APIs).

2. structure_figure(subseq, hl_start, hl_end, title) -> (fig, structure, mfe)
   Fold a window with ViennaRNA and draw an RNAfold-style 2-D secondary-structure
   diagram, with the candidate disruptor target site highlighted.
"""
import requests

NAVY, BLUE, RED, GREY = "#1F3864", "#2E6DB4", "#C0504D", "#BBBBBB"

SPECIES_COMMON = {
    "homo_sapiens": "Homo sapiens",
    "mus_musculus": "Mus musculus",
    "rattus_norvegicus": "Rattus norvegicus",
}


def _parse_fasta(text):
    name, seq = "transcript", []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            seq.append(line.upper().replace("T", "U"))
    return name, "".join(seq)


def fetch_by_name(symbol, species="homo_sapiens", timeout=25):
    """Resolve a gene/lncRNA symbol to a transcript sequence. Returns (label, seq)."""
    symbol = (symbol or "").strip()
    if not symbol:
        raise ValueError("Enter a gene / lncRNA name (e.g. MALAT1).")

    # --- 1) Ensembl REST: symbol -> gene -> canonical transcript cDNA ---
    try:
        base = "https://rest.ensembl.org"
        r = requests.get(f"{base}/lookup/symbol/{species}/{symbol}",
                         params={"expand": "1"},
                         headers={"Content-Type": "application/json"}, timeout=timeout)
        if r.ok:
            g = r.json()
            tid = (g.get("canonical_transcript") or "").split(".")[0]
            if not tid:
                trs = g.get("Transcript") or []
                if trs:
                    tid = trs[0].get("id", "")
            if tid:
                s = requests.get(f"{base}/sequence/id/{tid}",
                                 params={"type": "cdna"},
                                 headers={"Content-Type": "text/x-fasta"}, timeout=timeout)
                if s.ok and s.text.lstrip().startswith(">"):
                    _, seq = _parse_fasta(s.text)
                    if seq:
                        return f"{symbol} [Ensembl {tid}]", seq
    except Exception:
        pass

    # --- 2) NCBI E-utilities fallback ---
    try:
        eu = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        org = SPECIES_COMMON.get(species, "Homo sapiens")
        for term in (
            f'{symbol}[gene] AND "{org}"[orgn] AND srcdb_refseq[prop] AND biomol_ncrna[prop]',
            f'{symbol}[gene] AND "{org}"[orgn] AND srcdb_refseq[prop]',
        ):
            es = requests.get(f"{eu}/esearch.fcgi",
                              params={"db": "nuccore", "term": term, "retmax": "1",
                                      "retmode": "json"}, timeout=timeout)
            ids = es.json().get("esearchresult", {}).get("idlist", []) if es.ok else []
            if ids:
                ef = requests.get(f"{eu}/efetch.fcgi",
                                  params={"db": "nuccore", "id": ids[0],
                                          "rettype": "fasta", "retmode": "text"}, timeout=timeout)
                if ef.ok and ef.text.lstrip().startswith(">"):
                    _, seq = _parse_fasta(ef.text)
                    if seq:
                        return f"{symbol} [RefSeq {ids[0]}]", seq
    except Exception:
        pass

    raise RuntimeError(
        f"Could not auto-fetch '{symbol}'. Check the spelling and species, "
        f"or paste the sequence directly below.")


def structure_figure(subseq, hl_start, hl_end, title):
    """Fold subseq and draw an RNAfold-style 2-D diagram with the target site highlighted.
    Returns (matplotlib_figure, dot_bracket_structure, mfe)."""
    import RNA
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fc = RNA.fold_compound(subseq)
    structure, mfe = fc.mfe()
    xy = RNA.simple_xy_coordinates(structure)
    xs = [c.X for c in xy]
    ys = [c.Y for c in xy]
    pt = RNA.ptable(structure)

    fig, ax = plt.subplots(figsize=(6.6, 6.6))
    ax.plot(xs, ys, "-", color=GREY, lw=1.1, zorder=1)               # backbone
    for i in range(1, len(pt)):                                       # base pairs
        j = pt[i]
        if j > i:
            ax.plot([xs[i - 1], xs[j - 1]], [ys[i - 1], ys[j - 1]],
                    "-", color=BLUE, lw=0.8, alpha=0.55, zorder=1)
    for k, (x, y) in enumerate(zip(xs, ys)):                          # nucleotides
        hl = hl_start <= k < hl_end
        ax.scatter([x], [y], s=54 if hl else 22,
                   c=(RED if hl else "white"),
                   edgecolors=(RED if hl else "#999999"),
                   linewidths=1.0, zorder=3 if hl else 2)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=12, color=NAVY, pad=10)
    ax.legend([Line2D([0], [0], marker="o", color="w", markerfacecolor=RED, markersize=9),
               Line2D([0], [0], color=BLUE, lw=2)],
              ["disruptor target site", "base pair"],
              loc="lower right", fontsize=9, frameon=False)
    fig.text(0.02, 0.02, f"MFE {mfe:.1f} kcal/mol  |  {len(subseq)} nt window",
             fontsize=8, color="#666666")
    return fig, structure, mfe
