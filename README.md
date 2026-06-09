# PBE Pipeline — in-silico target & interface nomination (Stage 1)

Maps the protein-binding element (PBE) on a disease lncRNA and nominates candidate
steric-block / decoy disruptors, **using only public data and standard tools**. This
generates the kind of preliminary in-silico result that supports the Phase I feasibility
case, and is the open-source artifact for the dissemination plan.

Nothing here is a new algorithm — it orchestrates existing public resources
(POSTAR3 / ENCODE eCLIP + ViennaRNA/RNAfold) with a light overlay-and-tile script.

---

## 1. Install (free)

```bash
pip install ViennaRNA biopython numpy
# Optional, for the uniqueness check: NCBI BLAST+ (blastn) and bedtools
```

## 2. Get the public inputs (anchor example: MANTIS / SMARCA4–BRG1)

**(a) lncRNA transcript sequence** → `mantis.fa`
- MANTIS = lncRNA *n342419* (antisense to *GLCCI1*). Download the transcript FASTA from
  Ensembl, RefSeq, or the UCSC Table Browser. Save as a single-record FASTA.

**(b) Partner-RBP binding intervals, in TRANSCRIPT coordinates** → `smarca4_peaks.bed`
- Download SMARCA4 (BRG1) binding sites from **POSTAR3** (https://postar3.ncrnalab.org)
  or **ENCODE eCLIP** (https://www.encodeproject.org, search "SMARCA4 eCLIP";
  download the `bed narrowPeak`). These are in **genomic** coordinates.
- Convert genomic peaks → transcript coordinates (one-time preprocessing):
  ```bash
  # keep peaks overlapping the lncRNA locus, then map genome->transcript:
  bedtools intersect -a smarca4_genomic.narrowPeak -b mantis.bed12 -wa > hits.bed
  # map to transcript offsets using the transcript exon structure (BED12 / GTF).
  # Output columns required by the pipeline:  transcript_id  start  end  name  score
  ```
  (If you prefer, any tool that gives you RBP-bound *transcript* intervals works —
  the pipeline just needs `transcript_id  start  end  [name]  [score]`.)

**(c) (optional) SHAPE reactivities** → `mantis.shape`
- If you have in-vitro/in-cell SHAPE for the region, pass it to fold with experimental
  constraints (Deigan method). One value per nucleotide, or `position<TAB>reactivity`.
  Phase II generates this experimentally; it is optional for the in-silico nomination.

## 3. Run

```bash
python3 pbe_pipeline.py \
    --fasta mantis.fa \
    --peaks smarca4_peaks.bed \
    --shape mantis.shape \      # optional
    --win 30 --top 3 --ncand 15 \
    --outdir mantis_out
```

## 4. Outputs (`mantis_out/`)

| File | Contents |
|------|----------|
| `pbe_summary.tsv` | nominated PBE window(s): coords, mean RBP coverage, structuredness, score |
| `disruptor_candidates.tsv` | ranked antisense disruptors: sequence, approx Tm, GC%, target accessibility, self-structure ΔG |
| `structure.json` | sequence, MFE structure, per-base unpaired probability, RBP coverage track, PBE windows, top disruptors — drives the figure |

## 5. Finish the nomination (standard, separate steps)

- **Uniqueness:** BLAST each candidate against the transcriptome and reject those above
  the identity ceiling:
  ```bash
  blastn -query candidates.fa -db human_transcriptome -task blastn-short \
         -word_size 11 -outfmt 6 > blast.tsv
  ```
- **Chemistry note:** reported Tm is a DNA-like nearest-neighbor *ranking* value; 2'-MOE/LNA
  modifications raise the true Tm. Use it to rank, not as an absolute.

## 6. Scope / honesty

This step nominates candidates **in silico**. It does **not** prove the interaction or
the disruption — that is the Phase II wet-lab program (RIP/CLIP, SHAPE, EMSA/MST, cell
assays). Present pipeline output as *preliminary in-silico mapping*, not as validated data.

## 6b. Worked example included: NORAD – PUMILIO (motif mode, real data)

When the partner RBP has a **defined recognition motif**, you don't need a CLIP file at
all — scan for the motif directly. This repo ships a genuine worked example:

- `norad.fa` — the real NORAD transcript (NCBI **NR_027451.1**, ~5.3 kb).
- Pumilio binds the **PRE** motif `UGUANAUA`. Run:
  ```bash
  python3 pbe_pipeline.py --fasta norad.fa --motif UGUANAUA --flank 8 \
        --top 6 --ncand 15 --outdir norad_out
  python3 make_pbe_figure.py norad_out/structure.json norad_fig.png \
        "In-silico PRE mapping & disruptor nomination - NORAD (Pumilio sites)"
  ```
- Result: the pipeline recovers **15 PRE sites** (consistent with the ~17 PREs reported
  for NORAD), local-folds the 5.3 kb transcript for accessibility, ranks the most
  accessible PREs as candidate PBEs, and designs disruptors against them. Recovering the
  **independently known** binding elements is a clean validation of the Stage-1 workflow.

NOTE on framing: NORAD is a *protective* PUMILIO sponge, so this is a **method-validation /
benchmark** demonstration (does the pipeline find the right elements on real data?), not a
claim that NORAD should be disrupted therapeutically. PREs are AU-rich, so the reported
DNA-like Tm is low; 2'-MOE/LNA chemistry raises the real Tm (rank, don't read absolutely).

Modes summary:
- `--peaks file.bed`  use a CLIP footprint (transcript coords)
- `--motif UGUANAUA`  use a known RBP recognition motif (no CLIP file needed)
- `--local`           local folding for accessibility (auto for transcripts > 1500 nt)

## 7. Dissemination

A GitHub repository containing this script + the design rules satisfies the Challenge's
public-availability requirement for Phase I. A simple web wrapper (e.g., Streamlit/Galaxy)
that lets users run the pipeline without coding is an optional Phase II/III dissemination
deliverable — not required for the submission.
