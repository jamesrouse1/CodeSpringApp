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
