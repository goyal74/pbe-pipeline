#!/usr/bin/env python3
"""
PBE Pipeline — Streamlit web wrapper (optional Phase II/III dissemination tool)
===============================================================================
Lets users run the lncRNA PBE-mapping + disruptor-nomination pipeline in a browser,
with no coding: paste/upload a transcript, give a partner-RBP motif (or a CLIP BED),
and get the PBE map, accessibility, and ranked disruptor candidates.

Run locally:        streamlit run app.py
Deploy (free):      push the PBE_Pipeline folder to GitHub, then deploy on
                    Streamlit Community Cloud (share.streamlit.io) pointing at app.py.

Depends on pbe_pipeline.py (same folder) for the analysis functions.
"""
import io
import numpy as np
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pbe_pipeline as pp

st.set_page_config(page_title="lncRNA PBE Pipeline", layout="wide")
NAVY, BLUE, AMBER, TEAL, RED = "#1F3864", "#2E6DB4", "#E0922E", "#2E9B8F", "#C0504D"

st.title("lncRNA–RBP interface mapping & disruptor nomination")
st.caption("Map the protein-binding element (PBE) on a disease lncRNA and nominate "
           "steric-block / decoy disruptors — from public data and standard tools. "
           "In-silico nomination only; confirm experimentally.")

with st.sidebar:
    st.header("Input")
    fa = st.file_uploader("lncRNA transcript FASTA", type=["fa", "fasta", "txt"])
    seq_text = st.text_area("…or paste sequence (RNA/DNA)", height=120,
                            placeholder=">my_lncRNA\nACGU...")
    st.divider()
    mode = st.radio("Define the partner-RBP footprint by:",
                    ["Recognition motif", "CLIP BED (transcript coords)"])
    motif = st.text_input("RBP motif (IUPAC)", value="UGUANAUA",
                          help="e.g. Pumilio PRE = UGUANAUA") if mode == "Recognition motif" else None
    peaks_file = st.file_uploader("CLIP BED (id  start  end  [name]  [score])",
                                  type=["bed", "tsv", "txt"]) if mode != "Recognition motif" else None
    st.divider()
    top = st.slider("PBE windows to nominate", 1, 12, 6)
    flank = st.slider("Flank around motif (nt)", 0, 30, 8)
    ncand = st.slider("Disruptor candidates", 4, 30, 15)
    local = st.checkbox("Local folding (auto for >1500 nt)", value=True)
    go = st.button("Run analysis", type="primary")

def parse_fasta_text(t):
    name, s = "lncRNA", []
    for line in t.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            s.append(line.upper().replace("T", "U"))
    return name, "".join(s)

def get_sequence():
    if fa is not None:
        return parse_fasta_text(fa.getvalue().decode("utf-8", "ignore"))
    if seq_text.strip():
        return parse_fasta_text(seq_text)
    return None, None

if go:
    name, seq = get_sequence()
    if not seq:
        st.error("Provide a sequence (upload a FASTA or paste one)."); st.stop()
    n = len(seq)
    use_local = local or n > 1500
    with st.spinner(f"Folding {n} nt and mapping the interface…"):
        # footprint
        sites = None
        if mode == "Recognition motif":
            sites = pp.motif_scan(seq, motif)
            if not sites:
                st.error(f"Motif '{motif}' not found."); st.stop()
            peaks = [(s, e, 1.0) for (s, e) in sites]
        else:
            if peaks_file is None:
                st.error("Upload a CLIP BED or switch to motif mode."); st.stop()
            import tempfile, os
            tmp = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
            tmp.write(peaks_file.getvalue().decode("utf-8", "ignore")); tmp.close()
            peaks = pp.read_peaks(tmp.name); os.unlink(tmp.name)
            if not peaks:
                st.error("No peaks parsed."); st.stop()
        # accessibility
        unpaired = pp.local_accessibility(seq) if use_local else pp.fold(seq)[1]
        cov = pp.coverage(peaks, n)
        pbe = (pp.nominate_pbe_motif(sites, unpaired, flank=flank, top=top)
               if sites is not None else pp.nominate_pbe(cov, unpaired, top=top))
        disr = pp.design_disruptors(seq, unpaired, pbe)[:ncand]

    c1, c2, c3 = st.columns(3)
    c1.metric("Transcript length", f"{n} nt")
    c2.metric("Footprint sites", f"{len(sites) if sites is not None else len(peaks)}")
    c3.metric("PBE windows", f"{len(pbe)}")

    # figure
    x = np.arange(n)
    minw = max(1, int(n * 0.006))
    fig, axes = plt.subplots(2, 1, figsize=(11, 3.6), sharex=True,
                             gridspec_kw={"hspace": 0.15})
    def intervals(a):
        out, i = [], 0
        while i < len(a):
            if a[i] > 0:
                j = i
                while j < len(a) and a[j] > 0: j += 1
                out.append((i, j)); i = j
            else: i += 1
        return out
    for w in pbe:
        for ax in axes:
            ax.add_patch(Rectangle((w["start"] if isinstance(w, dict) else w[0], 0),
                                   max((w["end"] if isinstance(w, dict) else w[1]) -
                                       (w["start"] if isinstance(w, dict) else w[0]), minw),
                                   1, transform=ax.get_xaxis_transform(),
                                   color=AMBER, alpha=0.3, zorder=0))
    for (s, e) in intervals(cov):
        axes[0].add_patch(Rectangle((s, 0), max(e - s, minw), 1.0,
                          transform=axes[0].get_xaxis_transform(), color=RED, alpha=0.9, lw=0))
    axes[0].set_ylabel("RBP\nfootprint", fontsize=9, color=NAVY); axes[0].set_yticks([0, 1])
    axes[1].fill_between(x, unpaired, color=BLUE, alpha=0.75, lw=0)
    axes[1].set_ylabel("Accessibility", fontsize=9, color=NAVY); axes[1].set_yticks([0, 1])
    axes[1].set_xlabel("Transcript position (nt)", fontsize=9, color=NAVY)
    for c in disr:
        axes[1].add_patch(Rectangle((c["target_start"], -0.18), max(c["len"], minw), 0.11,
                          transform=axes[1].get_xaxis_transform(), clip_on=False, color=TEAL))
    fig.subplots_adjust(bottom=0.18, left=0.08, right=0.98, top=0.95)
    st.pyplot(fig)

    st.subheader("Top disruptor candidates")
    import pandas as pd
    cols = ["pbe", "target_start", "target_end", "len", "target_seq_5to3",
            "disruptor_antisense_5to3", "Tm_C_approx", "GC_pct", "target_accessibility"]
    df = pd.DataFrame(disr)[cols]
    st.dataframe(df, use_container_width=True)
    st.download_button("Download candidates (TSV)", df.to_csv(sep="\t", index=False),
                       file_name=f"{name}_disruptors.tsv")
    st.info("In-silico nomination only. Confirm uniqueness by BLAST and confirm the "
            "interaction/disruption experimentally (RIP/CLIP, SHAPE, EMSA/MST, cell assays).")
else:
    st.write("Configure inputs in the sidebar and click **Run analysis**. "
             "Example: paste the NORAD transcript and use motif **UGUANAUA** (Pumilio PRE).")
