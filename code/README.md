# TCR Overlap Analysis — Code for Publication

Analysis code accompanying the Treg manuscript by Cerosaletti et al. This code reproduces the
TCR repertoire overlap analyses comparing regulatory T cells (Treg) and conventional
T cells (Tconv) against published antigen-specific TCR reference datasets.

## Contents

| File | Description |
|------|-------------|
| `tcr_overlap_analysis.py` | Main analysis script (Tasks 1–5, see below) |
| `README.md` | This file |

## Analyses performed

1. **Removal of invariant TCRs.** iNKT and MAIT cells, which carry semi-invariant TCRs that
   would confound antigen-specific overlap analyses, are identified and removed from the S10,
   Michels, and Maki Nakayama datasets.
   - **iNKT**: TRAV10 + TRAJ18 **and** CDR3a exactly matches the canonical invariant sequence
     `CVVSDRGSTLGRLYF`. This rule is applied identically to all three datasets. Any deviation
     from the canonical CDR3 indicates non-germline junctional diversity at the V-J joint —
     i.e. a coincidental TRAV10+TRAJ18 rearrangement from an ordinary (non-invariant) T cell,
     not the true semi-invariant iNKT clone.
   - **MAIT**: TRAV1-2 + TRAJ33 / TRAJ20 / TRAJ12 (gene identity only, no CDR3 requirement).
     MAIT junctional diversity is broader than iNKT's, so no single canonical CDR3 is
     required.
   - Gene names are normalized before comparison (IMGT/committee prefix `TRAV`/`TRAJ`/`TRBV`/
     `TRBJ` and allele suffix `*01` stripped), because the Michels reference stores bare
     numeric/committee codes (e.g. `10`, `1-2`) while S10 and Maki use IMGT-prefixed names
     (e.g. `TRAV10`, `TRAV1-2`).

2. **Restrict Michels to T1D-associated antigens.** The Michels reference is a compiled
   dataset that includes both islet-autoantigen (T1D) and viral-antigen TCR entries. Downstream
   comparisons (Task 4) use only the T1D-associated subset (`T1D or VIR` == `T1D`).

3. **Treg vs Tconv TCR overlap within S10.** TCR chains shared between Treg and Tconv cells
   (within the iNKT/MAIT-filtered S10 dataset) are identified at two levels: (a) identical full
   nucleotide sequence, and (b) identical V-CDR3-J amino acid identity. TRA and TRB chains are
   analyzed separately.

4. **CDR3 sharing, S10 vs Michels (T1D-only).** Exact amino acid matching of CDR3 (junction)
   sequences of TRA and TRB chains between (filtered) S10, split by Treg / Tconv, and the
   (filtered, T1D-only) Michels islet-antigen TCR dataset. Results are reported as the number
   and percentage of unique S10 CDR3 sequences shared with the reference.

5. **CDR3 sharing, S10 vs Maki Nakayama.** Same as Task 4, but against the (filtered) Maki
   Nakayama pancreatic TCR dataset, using only in-frame primary-chain junctions from Maki.

## Conventions

- **Tconv** ("Tconventional") throughout this analysis denotes Tconventional and Tscm cells
  merged.
- **Chain-sharing matching (Task 3)** is exact string identity on the full nucleotide
  sequence, or on the V-CDR3-J amino acid string.
- **CDR3-sharing matching (Tasks 4–5)** is exact string identity on the CDR3 amino acid
  sequence alone — V/J gene usage is deliberately ignored, since a shared CDR3 arising from
  different V/J genes (convergent recombination) may still reflect antigen-relevant sequence
  convergence, and the reference datasets do not share a common V/J gene naming convention
  with S10.
- The **"unique CDR3" denominators** used in Tasks 4–5 count bare CDR3 amino acid sequences
  (not V-CDR3-J bioidentity), consistent with the CDR3-only matching key. This matters:
  counting by CDR3 alone collapses convergent recombination events (same CDR3, different V/J
  genes) into one entry, which is *not* the right denominator for a V-CDR3-J-level comparison
  — see the corresponding note in the script if adapting this code for a bioidentity-level
  comparison instead.
- `""` / `NA` / common placeholder strings (`Not reported`, `ND`, `N/A`, `None`) are treated as
  missing, case-insensitively.
- Column names are validated on load (`require_columns`) and the script raises an explicit
  error if an expected column is missing, rather than silently producing zero matches. This
  guards against the two classes of bugs found during development of this pipeline: a
  reference-file column-naming mismatch, and a gene-format mismatch between datasets that use
  different naming conventions for the same gene (both previously failed silently).

## Data sources

### S10 (this study)
The S10 T cell receptor sequences are provided as **supplemental data table 10** with the
manuscript. Data begin on spreadsheet row 3 (row 2 is the header).

Relevant columns (referenced by name, not position):
`libid`, `donorID`, `labID`, `cellType` (Treg / Tconventional / Tscm), `chain` (TRA/TRB/TRG),
`v_gene`, `j_gene`, `junction` (CDR3 amino acid), `V_CDR3_J`, `full_nt_sequence`.

### Michels islet-antigen TCR dataset (published; not included in this manuscript)
Islet- and viral-antigen-specific human TCR sequences compiled from the published literature.
This dataset is from **Mitchell et al. 2023, *Science Advances*** (see Citations) and is
**previously published**; it is **not** redistributed with this manuscript. Obtain it from the
original publication/supplementary materials or the corresponding authors.

Relevant columns (referenced by name): `TRAV`, `CDR3a`, `TRAJ`, `TRBV`, `CDR3b`, `TRBJ`,
`T1D or VIR`, `Antigen`, `Epitope`.

> **Note:** `TRAV`/`TRAJ`/`TRBV`/`TRBJ` in this dataset are bare numeric/committee codes
> (e.g. `10`, `1-2`, `18`), not IMGT-prefixed names — this is handled by gene-name
> normalization in the script, not by column position.

### Maki Nakayama pancreatic TCR dataset (published; not included in this manuscript)
Pancreas-infiltrating T cell (PIT) TCR sequences from healthy, at-risk, and T1D organ donors.
This dataset is from **Linsley et al. 2024, *Nature Communications*** (see Citations) and is
**previously published**; it is **not** redistributed with this manuscript. Obtain it from the
original publication/supplementary materials.

Relevant columns (referenced by name):

| Column name | Field |
|-------------|-------|
| `Group` | donor group (AAb, non-T1D, T1D) |
| `Vgene...6`, `Jgene...7`, `Junction...8` | primary alpha V gene / J gene / CDR3 junction |
| `frame...9` | primary alpha reading frame |
| `Vgene...11`, `Jgene...12`, `Junction...13` | secondary alpha V gene / J gene / CDR3 junction |
| `frame...14` | secondary alpha reading frame |
| `Vgene...16`, `Jgene...17`, `Junction...18` | primary beta V gene / J gene / CDR3 junction |
| `frame...19` | primary beta reading frame |

Both alpha slots (primary and secondary) are screened for iNKT/MAIT removal (Task 1). Only the
in-frame primary chains (`frame...9` == `in-frame`, `frame...19` == `in-frame`) are used for
CDR3 sharing (Task 5).

## Requirements

- Python >= 3.8
- pandas
- openpyxl

```bash
pip install pandas openpyxl
```

## Usage

1. Place the three input files in a working directory (S10 supplemental table; Michels and
   Maki Nakayama datasets obtained from their original sources).
2. Edit the `CONFIG` section at the top of `tcr_overlap_analysis.py` to point to the input
   files and desired output directory.
3. Run:

```bash
python tcr_overlap_analysis.py
```

### Outputs

| File | Content |
|------|---------|
| `S10_iNKT-MAIT_filtered.xlsx` | S10 with iNKT/MAIT cells removed |
| `S10_iNKT-MAIT_removed_rows.xlsx` | the removed S10 rows (all chains, for QC) |
| `Michels_iNKT-MAIT_filtered.xlsx` | Michels with iNKT/MAIT entries removed |
| `Michels_iNKT-MAIT_removed_rows.xlsx` | the removed Michels rows |
| `Maki_iNKT-MAIT_filtered.xlsx` | Maki Nakayama with iNKT/MAIT entries removed |
| `Maki_iNKT-MAIT_removed_rows.xlsx` | the removed Maki rows |
| `Michels_iNKT-MAIT_filtered_T1Donly.xlsx` | filtered Michels restricted to T1D antigens |
| `S10_Treg_vs_Tconv_overlap.xlsx` | Treg–Tconv shared chains (identical NT; shared V-CDR3-J), TRA and TRB |
| `CDR3_sharing_S10_vs_Michels_T1Donly.xlsx` | S10 CDR3 sequences shared with Michels (T1D-only), with antigen/epitope |
| `CDR3_sharing_S10_vs_Maki.xlsx` | S10 CDR3 sequences shared with Maki Nakayama |

Match summaries (unique matched / denominator / % of S10, by Treg/Tconv and TRA/TRB) are
printed to standard output for Tasks 3–5.

## Notes on reproducibility

- No random sampling is performed; results are deterministic.
- This revised pipeline was validated by reproducing three independently-generated prior-session
  reference files exactly (S10: 24→21 after adopting the exact-CDR3-for-all-three rule; Michels:
  9; Maki: 143 removed), and by cross-checking derived unique-chain counts (bioidentity vs bare
  CDR3) against a previously-computed antigen-stratified summary table.

## Citations

The reference TCR datasets used in this analysis are previously published and must be
obtained from their original sources:

- **Michels islet-antigen TCR dataset:**
  Mitchell AM, Baschal EE, McDaniel KA, Fleury T, Choi H, Pyle L, Yu L, Rewers MJ,
  Nakayama M, Michels AW. Tracking DNA-based antigen-specific T cell receptors during
  progression to type 1 diabetes. *Science Advances*. 2023 Dec 8;9(49):eadj6975.
  doi:10.1126/sciadv.adj6975. PMID: 38064552.

- **Maki pancreatic TCR dataset:**
  Linsley PS, Nakayama M, Balmas E, Chen J, Barahmand-pour-Whitman F, Bansal S,
  Bottorff T, Serti E, Speake C, Pugliese A, Cerosaletti K. Germline-like TCR-α chains
  shared between autoreactive T cells in blood and pancreas. *Nature Communications*.
  2024 Jun 11;15:4971. doi:10.1038/s41467-024-48833-w. PMID: 38871688.
