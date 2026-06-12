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
import pbe_webtools as wt

st.set_page_config(page_title="PBE Pipeline · Epigenuity LLC", layout="wide")
NAVY, BLUE, AMBER, TEAL, RED = "#1F3864", "#2E6DB4", "#E0922E", "#2E9B8F", "#C0504D"

st.title("RNA–protein interface mapping & disruptor nomination")
st.markdown("**Developed by Ravi Goyal, MD, PhD · Epigenuity LLC** &nbsp;·&nbsp; "
            "[github.com/goyal74/pbe-pipeline](https://github.com/goyal74/pbe-pipeline) "
            "&nbsp;·&nbsp; DOI [10.5281/zenodo.20615550](https://doi.org/10.5281/zenodo.20615550)")
st.caption("Map the protein-binding element (PBE) on a disease RNA (lncRNA or toxic "
           "structured/repeat RNA) and nominate steric-block / decoy disruptors — from "
           "public data and standard tools. In-silico nomination only; confirm experimentally.")

with st.sidebar:
    st.header("Input")
    name_in = st.text_input("lncRNA / gene name (auto-fetch)",
                            placeholder="e.g. MALAT1, NEAT1, NORAD",
                            help="Type a name and the sequence is fetched automatically "
                                 "from Ensembl / NCBI — no need to paste it.")
    species = st.selectbox("Species", ["homo_sapiens", "mus_musculus", "rattus_norvegicus"],
                           index=0, format_func=lambda s: s.replace("_", " ").title())
    st.caption("…or provide a sequence directly:")
    fa = st.file_uploader("Transcript FASTA", type=["fa", "fasta", "txt"])
    seq_text = st.text_area("…or paste sequence (RNA/DNA)", height=100,
                            placeholder=">my_RNA\nACGU...")
    st.divider()
    mode = st.radio("Find the protein-binding site by:",
                    ["Auto — scan known RBP motifs (no input needed)",
                     "One recognition motif",
                     "CLIP BED (transcript coords)"])
    if mode.startswith("Auto"):
        st.caption(f"Scans the transcript against {len(pp.RBP_MOTIFS)} known RBP "
                   "recognition motifs (PUM, MBNL, HuR, PTBP1, TDP-43, hnRNPA1, QKI, "
                   "Nova, CELF1, …) and reports which RBP each candidate site matches.")
    motif = st.text_input("RBP motif (IUPAC)", value="UGUANAUA",
                          help="e.g. Pumilio PRE = UGUANAUA") if mode == "One recognition motif" else None
    peaks_file = st.file_uploader("CLIP BED (id  start  end  [name]  [score])",
                                  type=["bed", "tsv", "txt"]) if mode.startswith("CLIP") else None
    st.divider()
    top = st.slider("PBE windows to nominate", 1, 12, 6)
    flank = st.slider("Flank around motif (nt)", 0, 30, 8)
    ncand = st.slider("Disruptor candidates", 4, 30, 15)
    local = st.checkbox("Local folding (auto for >1500 nt)", value=True)
    go = st.button("Run analysis", type="primary")
    st.divider()
    st.caption("© 2026 Epigenuity LLC · MIT License")

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
    if name_in.strip():
        return wt.fetch_by_name(name_in, species)   # may raise; handled by caller
    if fa is not None:
        return parse_fasta_text(fa.getvalue().decode("utf-8", "ignore"))
    if seq_text.strip():
        return parse_fasta_text(seq_text)
    return None, None

if go:
    try:
        with st.spinner("Fetching sequence…" if name_in.strip() else "Reading sequence…"):
            name, seq = get_sequence()
    except Exception as e:
        st.error(str(e)); st.stop()
    if not seq:
        st.error("Provide an input: type a gene/lncRNA name, upload a FASTA, or paste a sequence.")
        st.stop()
    if name_in.strip():
        st.success(f"Fetched **{name}** — {len(seq)} nt.")
    n = len(seq)
    use_local = local or n > 1500
    with st.spinner(f"Folding {n} nt and mapping the interface…"):
        # footprint
        sites = None
        hits_by_rbp = {}
        if mode.startswith("Auto"):
            all_sites, hits_by_rbp = pp.scan_rbp_library(seq)
            if not all_sites:
                st.error("No known RBP recognition-motif sites found in this transcript.")
                st.stop()
            peaks = [(s, e, 1.0) for (s, e) in all_sites]
        elif mode == "One recognition motif":
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
        # annotate each PBE window with the RBP(s) whose motif lands there (auto mode)
        pbe_rbp = {}
        if hits_by_rbp:
            for (ps, pe, _) in pbe:
                pbe_rbp[f"{ps}-{pe}"] = ", ".join(pp.rbp_labels_for_window(ps, pe, hits_by_rbp)) or "—"

    c1, c2, c3 = st.columns(3)
    c1.metric("Transcript length", f"{n} nt")
    c2.metric("RBPs detected" if hits_by_rbp else "Footprint sites",
              f"{len(hits_by_rbp)}" if hits_by_rbp else f"{len(sites) if sites is not None else len(peaks)}")
    c3.metric("PBE windows", f"{len(pbe)}")

    if hits_by_rbp:
        summary = sorted(((rbp, len(s)) for rbp, s in hits_by_rbp.items()),
                         key=lambda t: -t[1])
        st.markdown("**Candidate RBPs found (motif hits across the transcript):** "
                    + " · ".join(f"{rbp} ({c})" for rbp, c in summary))

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
    df = pd.DataFrame(disr)
    if hits_by_rbp:
        df.insert(1, "candidate_RBP", df["pbe"].map(lambda p: pbe_rbp.get(p, "—")))
    cols = (["pbe", "candidate_RBP"] if hits_by_rbp else ["pbe"]) + [
        "target_start", "target_end", "len", "target_seq_5to3",
        "disruptor_antisense_5to3", "Tm_C_approx", "GC_pct", "target_accessibility"]
    df = df[cols]
    st.dataframe(df, use_container_width=True)
    st.download_button("Download candidates (TSV)", df.to_csv(sep="\t", index=False),
                       file_name=f"{name}_disruptors.tsv")

    # --- RNAfold-style structure diagram around the top disruption site ---
    if disr:
        top = disr[0]
        center = (top["target_start"] + top["target_end"]) // 2
        half = 70
        w0, w1 = max(0, center - half), min(n, center + half)
        sub = seq[w0:w1]
        hl0, hl1 = max(0, top["target_start"] - w0), min(len(sub), top["target_end"] - w0)
        st.subheader("Secondary structure at the disruption site (RNAfold)")
        with st.spinner("Folding the disruption-site window…"):
            sfig, sstruct, smfe = wt.structure_figure(
                sub, hl0, hl1,
                f"Predicted structure around the top disruptor site (nt {w0 + 1}–{w1})")
        c1, c2 = st.columns([3, 2])
        c1.pyplot(sfig)
        top_rbp = pbe_rbp.get(top["pbe"], "") if hits_by_rbp else ""
        c2.markdown(
            f"**Top disruptor**\n\n"
            + (f"- Candidate RBP: **{top_rbp}**\n" if top_rbp and top_rbp != "—" else "")
            + f"- Target site: **nt {top['target_start'] + 1}–{top['target_end']}**\n"
            f"- Target (5′→3′): `{top['target_seq_5to3']}`\n"
            f"- Antisense disruptor: `{top['disruptor_antisense_5to3']}`\n"
            f"- Predicted Tm: **{top['Tm_C_approx']:.1f} °C** · GC {top['GC_pct']:.0f}%\n"
            f"- Window MFE: **{smfe:.1f} kcal/mol**\n\n"
            f"Red nucleotides mark where the steric-block / decoy disruptor binds, "
            f"occluding the protein-binding element without degrading the RNA.")

    st.info("In-silico nomination only. Confirm uniqueness by BLAST and confirm the "
            "interaction/disruption experimentally (RIP/CLIP, SHAPE, EMSA/MST, cell assays).")
else:
    st.write("**Just type an lncRNA / gene name** in the sidebar (e.g. **MALAT1**, **NEAT1**, "
             "**NORAD**) and click **Run analysis** — the app fetches the sequence, scans it "
             "against a library of known RBP motifs to find the protein-binding sites itself, "
             "ranks the best-disruptable site, designs the antisense disruptor, and draws the "
             "RNAfold structure with the disruption site marked. No motif or sequence needed.")

st.divider()
st.caption("PBE Pipeline · Developed by **Ravi Goyal, MD, PhD, Epigenuity LLC** (Tucson, AZ) · "
           "© 2026 Epigenuity LLC · MIT License · Cite: doi.org/10.5281/zenodo.20615550")
