# nf-core/cutandrun workshop demo

## Data

K562 CUT&RUN from GEO **GSE145187**. These are the files nf-core uses for its own full-size cutandrun test.

| Sample | Mark | Class | Control |
|---|---|---|---|
| K562_H3K4me3_rep1 / rep2 | H3K4me3 | histone_narrow | K562_IgG_rep1 / rep2 |
| K562_H3K27me3_rep1 / rep2 | H3K27me3 | histone_broad | K562_IgG_rep1 / rep2 |
| K562_IgG_rep1 / rep2 | IgG | control | — |

- **Reads:** 25 bp paired-end, about 1.6–2.1M read pairs per sample.
- **Genome:** human GRCh38, GENCODE 50.
- **Location:** shared, read-only, at `/grid/bsr/data/data/bsr_readable_data/CodeSpringFlow/examples/cutrun_nfcore_GSE145187`. To point the app at another copy, set `CSL_CUTRUN_EXAMPLE_ROOT`.

## Before the workshop

1. Each attendee needs the app running from a CodeSpringApp checkout in their own home folder. They also need read access to the shared folder above.
2. Do one test run from an account outside the `bioinformatics` group.
3. Make `openxlsx` and `DT` available in the R library that attendees use.

## Live run (attendees)

1. Go to **New project** and set the analysis type to **CUT&RUN**.
2. Click **Use Example Dataset**. This fills in the FASTQ folder, the design matrix, human genome and the nf-core backend.
3. Choose a project name and click **Create project**.
4. Optionally, open **Samples & Design** and show the IgG matching: rep1 targets use IgG rep1, rep2 targets use IgG rep2.
5. Open **Run pipeline**, tick MACS2 if you want both peak callers, and click **Run**.
6. Watch **Progress**: stages go from Preprocessing through Alignment, Dedup & QC, Signal, Peak calling, Peak QC and Reports.

The reference run took **about 71 minutes** on `cpuq`: 264 tasks, from 15:05 to 16:16. Attendees' runs will take longer if they all submit at once.

## Presenting the finished run

The presenter's account has the project **`workshop_cutrun_demo`**, which is already completed. Its results are in `/grid/bsr/data/data/gajri/csl_results/workshop_cutrun_demo/data/nfcore_cutandrun/results`.

Things to point out:

- **Alignment (Results Explorer → QC):** targets align at 95–98%. IgG rep1 aligns at 85%.
- **Duplication:** targets are at 1–7% duplicates. IgG rep1 is at **86%** and rep2 at 36%. IgG libraries are low-complexity by design. This is a good moment to explain why a control can't be judged with the same rules as a target.
- **Peaks:** SEACR (stringent, against IgG) finds about 4–6k H3K4me3 peaks and 6–11k H3K27me3 peaks per replicate. MACS2 narrowPeak finds about 8–11k H3K4me3 peaks. Compare narrow H3K4me3 promoter peaks with broad H3K27me3 domains in the **Genome Browser**.
- **MultiQC:** this has the fragment-length distribution (look for the nucleosome ladder), FRiP and the peak-reproducibility plots.
- **One task "failed" but the run completed.** Preseq's library-complexity extrapolation failed for H3K27me3 rep2, and nf-core tolerates this error. The Progress tab shows the stage as Complete with 1 in the failed column. Expect the same thing in attendee runs, because Preseq often fails on small inputs.

# Differential peaks example

nf-core/cutandrun stops at peak calling. Differential binding is a separate DiffBind/DESeq2 step that you start by hand in the app. The K562 example has no second condition, so a second example exists to show this step.

## Data

NUT carcinoma cell line 10-15 from GEO **GSE228533** (BioProject PRJNA950247). It is H3K27me3 CUT&RUN after treatment with DMSO or tazemetostat, an EZH2 inhibitor.

| Sample | Condition | Class | Control | SRA run |
|---|---|---|---|---|
| NC1015_H3K27me3_DMSO_rep1 / rep2 | DMSO | histone_broad | NC1015_IgG_DMSO_rep1 / rep2 | SRR24013104 / SRR24013113 |
| NC1015_H3K27me3_Taz_rep1 / rep2 | Taz | histone_broad | NC1015_IgG_Taz_rep1 / rep2 | SRR24013119 / SRR24013110 |
| NC1015_IgG_DMSO_rep1 / rep2 | DMSO | control | — | SRR24013121 / SRR24013112 |
| NC1015_IgG_Taz_rep1 / rep2 | Taz | control | — | SRR24013118 / SRR24013109 |

- **Reads:** the first 1,000,000 read pairs of each run, 51 bp paired-end. The folder is 410 MB in total.
- **Genome:** human GRCh38.
- **Location:** shared, read-only, at `/grid/bsr/data/data/bsr_readable_data/CodeSpringFlow/examples/cutrun_nfcore_GSE228533`. To point the app at another copy, set `CSL_CUTRUN_DIFF_EXAMPLE_ROOT`. The folder's `README.txt` records where the data came from.

## Running it

1. Go to **New project**, set the analysis type to **CUT&RUN** and choose **Differential peaks** under the example choices.
2. Click **Use Example Dataset**, name the project and click **Create project**.
3. Run nf-core as for the K562 example. The reference run took **about 82 minutes** on `cpuq`.
4. Open the **Differential Peaks** panel and set:
   - **Reference condition:** `DMSO`
   - **Peak source:** nf-core SEACR stringent
   - **Minimum called peaks per sample:** leave it at 100.

   The panel shows one eligible comparison, `NC1015 — H3K27me3 — Taz vs DMSO`, with 2 + 2 replicates.
5. Click **Submit all eligible comparisons**. The Slurm job took 2 min 20 s.
6. Optionally, refresh the **Differential Peak Summary** panel.

## What to point out

The presenter's account has the completed project **`workshop_cutrun_diff_demo`**, with nf-core and DiffBind already run. Its DiffBind results are in `/grid/bsr/data/data/gajri/csl_results/workshop_cutrun_diff_demo/data/cutrun_diffbind/nfcore_seacr_stringent__NC1015__H3K27me3__Taz_vs_DMSO`.

- **Peak calling already shows the effect.** SEACR finds 18.5k and 19.4k H3K27me3 peaks in the DMSO replicates, but only 2.0k and 1.4k in the tazemetostat replicates.
- **DiffBind tests 27,373 consensus peaks.** It finds **551 significant peaks** (FDR < 0.05): **538 lose H3K27me3** with tazemetostat and 13 gain it. Loss of the mark is what you expect when EZH2, the enzyme that writes it, is inhibited.
- **Normalisation.** The tazemetostat libraries have about half as many reads in peaks as DMSO. `normalization_factors.tsv` shows that DiffBind normalised on genomic background bins, not on reads in peaks. Normalising on reads in peaks would hide a global loss like this one, so this is a good place to discuss the choice. Spike-in normalisation is the more robust option when a whole mark is expected to drop.
- **Plots:** the volcano, MA, PCA and correlation heatmap PNGs are in the results folder. HOMER annotation of the significant peaks is in `significant_differential_peaks_annotated.tsv`.
- **Small input.** Each library is only 1M read pairs, so DiffBind has limited power. The full-depth data would find many more differential peaks.
