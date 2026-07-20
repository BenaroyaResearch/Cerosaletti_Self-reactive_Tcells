#!/usr/bin/env python3
"""
TCR overlap analysis for Treg vs Tconventional cells (Treg manuscript, JCI Insight).
Revised analysis.

This script performs five tasks:

    1. Remove invariant iNKT and MAIT TCRs from the S10, Michels, and Maki datasets.
    2. Restrict the Michels reference dataset to T1D-associated antigens only
       (excludes the viral-antigen entries also present in that compiled dataset).
    3. Treg vs Tconv TCR overlap within (filtered) S10: identical full nucleotide
       sequence, and identical V-CDR3-J amino acid identity. TRA and TRB separately.
    4. CDR3 amino acid sharing between (filtered) S10 (Treg / Tconv) and (filtered,
       T1D-only) Michels, TRA and TRB separately.
    5. CDR3 amino acid sharing between (filtered) S10 (Treg / Tconv) and (filtered)
       Maki Nakayama pancreatic TCRs, TRA and TRB separately.

Conventions
    - Tconv ("Tconventional") = Tconventional + Tscm merged.
    - iNKT: TRAV10 + TRAJ18 AND CDR3a exactly matches the canonical invariant sequence
      CVVSDRGSTLGRLYF. This rule is applied identically to all three datasets. Any
      deviation from the canonical CDR3 indicates non-germline junctional diversity at
      the V-J joint, i.e. a coincidental TRAV10+TRAJ18 rearrangement from an ordinary
      (non-invariant) T cell rather than the true semi-invariant iNKT clone.
    - MAIT: TRAV1-2 + TRAJ33 / TRAJ20 / TRAJ12 (gene identity only). Unlike iNKT, MAIT
      CDR3 is not held to a single exact sequence, since MAIT junctional diversity is
      documented to be broader than iNKT's.
    - Gene names are normalized before comparison (IMGT/committee prefix TRAV/TRAJ/
      TRBV/TRBJ and allele suffix *01 stripped), because the Michels reference stores
      bare numeric/committee codes (e.g. "10", "1-2") while S10 and Maki use
      IMGT-prefixed names (e.g. "TRAV10", "TRAV1-2").
    - Chain-sharing (Task 3) matching is exact string identity on full nucleotide
      sequence, or on the V-CDR3-J amino acid string. CDR3-sharing (Tasks 4-5) matching
      is exact string identity on the CDR3 amino acid sequence alone -- V/J gene usage
      is deliberately ignored there, since a shared CDR3 arising from different V/J
      genes (convergent recombination) may still reflect antigen-relevant sequence
      convergence, and because the reference datasets do not share a common V/J gene
      naming convention with S10.
    - "" / "NA" / common placeholder strings ("Not reported", "ND", "N/A", "None") are
      treated as missing, case-insensitively.
    - The "unique CDR3" denominators used in Tasks 4-5 count bare CDR3 amino acid
      sequences (not V-CDR3-J bioidentity), consistent with the CDR3-only matching key.

Requirements
    python >= 3.8, pandas, openpyxl

Usage
    Edit the CONFIG section below, then:  python tcr_overlap_analysis.py
"""

import re
import sys
import pandas as pd

# ============================================================================
# CONFIG -- edit paths for your environment
# ============================================================================
DATA_DIR = "."                      # directory containing the input files
OUT_DIR  = "."                      # directory for output files

S10_FILE     = f"{DATA_DIR}/Table S10_revised 062326.xlsx"
MICHELS_FILE = f"{DATA_DIR}/Combined islet antigen TCR list 08_13_2025_Aaron_Michels.xlsx"
MAKI_FILE    = f"{DATA_DIR}/maki's_pancreatic_TCRs_filtered_060926.csv"

CANONICAL_INKT_CDR3 = "CVVSDRGSTLGRLYF"

# ============================================================================
# Shared helpers
# ============================================================================
_PLACEHOLDERS = {"", "NA", "NOT REPORTED", "ND", "N/A", "NONE"}


def _clean(x):
    """Normalize a cell value to a stripped string; '' for missing/placeholder values."""
    if pd.isna(x):
        return ""
    s = str(x).strip()
    return "" if s.upper() in _PLACEHOLDERS else s


def normalize_gene(gene):
    """Strip an optional TRAV/TRAJ/TRBV/TRBJ prefix and any allele suffix (*01).

    Handles both IMGT-prefixed strings ('TRAV10', 'TRAJ18*01') and bare
    numeric/committee codes ('10', '1-2', '18'), so the same invariant-gene
    check works across datasets that use either naming convention.
    """
    s = _clean(gene)
    s = re.sub(r"\*.*$", "", s)
    s = re.sub(r"^(TRAV|TRAJ|TRBV|TRBJ)", "", s, flags=re.IGNORECASE)
    return s


def is_mait(v_gene, j_gene):
    """MAIT invariant alpha chain: TRAV1-2 + TRAJ33/20/12 (gene identity only)."""
    return normalize_gene(v_gene) == "1-2" and normalize_gene(j_gene) in ("33", "20", "12")


def invariant_type(v_gene, j_gene, cdr3):
    """Classify an alpha-chain rearrangement as 'iNKT', 'MAIT', or '' (not invariant).

    iNKT requires both the canonical genes (TRAV10+TRAJ18) AND the canonical CDR3.
    MAIT requires only the canonical genes.
    """
    v, j = normalize_gene(v_gene), normalize_gene(j_gene)
    if v == "10" and j == "18" and _clean(cdr3) == CANONICAL_INKT_CDR3:
        return "iNKT"
    if is_mait(v_gene, j_gene):
        return "MAIT"
    return ""


def cell_group(cell_type):
    """Map cellType to Treg or Tconv (Tconv = Tconventional + Tscm)."""
    return "Treg" if _clean(cell_type) == "Treg" else "Tconv"


def require_columns(df, columns, dataset_name):
    """Fail loudly if expected columns are missing, instead of silently no-op'ing.

    Two earlier bugs in this pipeline (a Maki column-name mismatch, and a Michels
    gene-format mismatch) each silently produced zero matches with no error. This
    check converts that failure mode into an immediate, explicit error.
    """
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(
            f"{dataset_name}: expected column(s) not found: {missing}. "
            f"Available columns: {list(df.columns)}"
        )


def export_xlsx(path, df):
    df.to_excel(path, index=False)
    print(f"  Wrote {path}  ({len(df)} rows)")


def export_sheets(path, sheets):
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in sheets.items():
            (df if isinstance(df, pd.DataFrame) else pd.DataFrame(df)) \
                .to_excel(xw, sheet_name=name[:31], index=False)
    print(f"  Wrote {path}")


# ============================================================================
# TASK 1 -- Remove iNKT and MAIT TCRs from S10, Michels, and Maki
# ============================================================================
def remove_invariant_s10(path):
    """Return (kept, removed) DataFrames for S10 with iNKT/MAIT cells removed.

    A cell (libid) is dropped if any of its TRA rows is an invariant iNKT or MAIT
    rearrangement. All rows belonging to a dropped libid are removed (both chains).
    Header row is spreadsheet row 2 (0-based header=1); data begins row 3.
    """
    df = pd.read_excel(path, header=1)
    df.columns = [str(c).strip() for c in df.columns]
    require_columns(df, ["libid", "chain", "v_gene", "j_gene", "junction"], "S10")

    invariant_libids = set()
    for _, row in df.iterrows():
        if _clean(row["chain"]) == "TRA":
            if invariant_type(row["v_gene"], row["j_gene"], row["junction"]):
                invariant_libids.add(_clean(row["libid"]))

    is_removed = df["libid"].map(_clean).isin(invariant_libids)
    kept, removed = df[~is_removed].copy(), df[is_removed].copy()
    print(f"S10: removed {len(invariant_libids)} iNKT/MAIT cells "
          f"({len(removed)} rows across all chains); {len(kept)} rows remain (of {len(df)}).")
    return kept, removed


def remove_invariant_michels(path):
    """Return (kept, removed) DataFrames for Michels with iNKT/MAIT entries removed."""
    df = pd.read_excel(path, header=0)
    df.columns = [str(c).strip() for c in df.columns]
    require_columns(df, ["TRAV", "TRAJ", "CDR3a"], "Michels")

    types = df.apply(lambda r: invariant_type(r["TRAV"], r["TRAJ"], r["CDR3a"]), axis=1)
    is_removed = types != ""
    kept, removed = df[~is_removed].copy(), df[is_removed].copy()
    print(f"Michels: removed {int(is_removed.sum())} iNKT/MAIT entries; "
          f"{len(kept)} rows remain (of {len(df)}).")
    return kept, removed


def remove_invariant_maki(path):
    """Return (kept, removed) DataFrames for Maki with iNKT/MAIT entries removed.

    Both alpha slots are screened (primary: Vgene...6/Jgene...7/Junction...8;
    secondary: Vgene...11/Jgene...12/Junction...13).
    """
    df = pd.read_csv(path)
    df.columns = [str(c).strip() for c in df.columns]
    require_columns(df, ["Vgene...6", "Jgene...7", "Junction...8",
                          "Vgene...11", "Jgene...12", "Junction...13"], "Maki")

    slots = [("Vgene...6", "Jgene...7", "Junction...8"),
             ("Vgene...11", "Jgene...12", "Junction...13")]

    def row_is_invariant(r):
        for vcol, jcol, jncol in slots:
            if invariant_type(r[vcol], r[jcol], r[jncol]):
                return True
        return False

    is_removed = df.apply(row_is_invariant, axis=1)
    kept, removed = df[~is_removed].copy(), df[is_removed].copy()
    print(f"Maki: removed {int(is_removed.sum())} iNKT/MAIT entries; "
          f"{len(kept)} rows remain (of {len(df)}).")
    return kept, removed


# ============================================================================
# TASK 2 -- Restrict Michels to T1D-associated antigens only
# ============================================================================
def restrict_michels_t1d(michels):
    """Return Michels rows where 'T1D or VIR' == 'T1D' (excludes viral entries)."""
    require_columns(michels, ["T1D or VIR"], "Michels")
    t1d = michels[michels["T1D or VIR"].map(_clean) == "T1D"].copy()
    print(f"Michels: restricted to T1D antigens; {len(t1d)} rows remain (of {len(michels)}).")
    return t1d


# ============================================================================
# TASK 3 -- Treg vs Tconv TCR overlap within S10 (NT identical; V-CDR3-J identical)
# ============================================================================
def s10_treg_tconv_overlap(s10):
    """Identify TCR chains shared between Treg and Tconv cells within S10.

    Two definitions of "shared":
      (A) identical full nucleotide sequence, and
      (B) identical V-CDR3-J amino acid identity.
    TRA and TRB are analyzed separately.
    """
    require_columns(s10, ["chain", "full_nt_sequence", "V_CDR3_J", "cellType",
                          "libid", "labID"], "S10")

    results = {}
    for chain in ("TRA", "TRB"):
        sub = s10[s10["chain"].map(_clean) == chain].copy()
        sub["grp"] = sub["cellType"].map(cell_group)
        sub["nt"] = sub["full_nt_sequence"].map(_clean)
        sub["vcdj"] = sub["V_CDR3_J"].map(_clean)

        nt_treg = set(sub.loc[(sub.grp == "Treg") & (sub.nt != ""), "nt"])
        nt_tconv = set(sub.loc[(sub.grp == "Tconv") & (sub.nt != ""), "nt"])
        shared_nt = nt_treg & nt_tconv
        results[f"NT_{chain}"] = (
            sub[sub.nt.isin(shared_nt)][["libid", "labID", "grp", "nt", "vcdj"]]
            .sort_values(["nt", "grp"])
        )

        aa_treg = set(sub.loc[(sub.grp == "Treg") & (sub.vcdj != ""), "vcdj"])
        aa_tconv = set(sub.loc[(sub.grp == "Tconv") & (sub.vcdj != ""), "vcdj"])
        shared_aa = aa_treg & aa_tconv
        results[f"VCDJ_{chain}"] = (
            sub[sub.vcdj.isin(shared_aa)][["libid", "labID", "grp", "vcdj"]]
            .sort_values(["vcdj", "grp"])
        )

        print(f"S10 Treg-Tconv sharing [{chain}]: "
              f"{len(shared_nt)} identical-NT, {len(shared_aa)} identical V-CDR3-J.")
    return results


# ============================================================================
# TASKS 4-5 -- CDR3 sharing: S10 vs Michels (T1D-only) and S10 vs Maki
# ============================================================================
def _s10_cdr3_by_group(s10):
    """Return {(grp, chain): set(unique CDR3 aa)} and per-cell rows for S10."""
    require_columns(s10, ["chain", "junction", "cellType", "libid", "labID", "donorID"], "S10")
    denom, rows = {}, []
    for _, r in s10.iterrows():
        chain, cdr3 = _clean(r["chain"]), _clean(r["junction"])
        if chain not in ("TRA", "TRB") or cdr3 == "":
            continue
        grp = cell_group(r["cellType"])
        denom.setdefault((grp, chain), set()).add(cdr3)
        rows.append({"chain": chain, "cell_group": grp, "CDR3": cdr3,
                      "libid": _clean(r["libid"]), "labID": _clean(r["labID"]),
                      "donor": _clean(r["donorID"])})
    return denom, pd.DataFrame(rows)


def _michels_cdr3_sets(michels):
    """CDR3 amino acid sets for Michels: {'TRA': set, 'TRB': set}, plus antigen/epitope lookup."""
    require_columns(michels, ["CDR3a", "CDR3b", "Antigen", "Epitope"], "Michels")
    sets = {"TRA": set(), "TRB": set()}
    annot = {"TRA": {}, "TRB": {}}
    for _, r in michels.iterrows():
        for chain, ccol in (("TRA", "CDR3a"), ("TRB", "CDR3b")):
            c = _clean(r[ccol])
            if c:
                sets[chain].add(c)
                annot[chain].setdefault(c, set()).add((_clean(r["Antigen"]), _clean(r["Epitope"])))
    return sets, annot


def _maki_cdr3_sets(maki):
    """CDR3 amino acid sets for Maki, in-frame primary chains only: {'TRA': set, 'TRB': set}."""
    require_columns(maki, ["Junction...8", "frame...9", "Junction...18", "frame...19"], "Maki")
    sets = {"TRA": set(), "TRB": set()}
    spec = [("TRA", "Junction...8", "frame...9"), ("TRB", "Junction...18", "frame...19")]
    for _, r in maki.iterrows():
        for chain, jcol, fcol in spec:
            if _clean(r[fcol]) == "in-frame":
                j = _clean(r[jcol])
                if j:
                    sets[chain].add(j)
    return sets


def _print_match_summary(match_df, denom):
    for grp in ("Treg", "Tconv"):
        for chain in ("TRA", "TRB"):
            den = len(denom.get((grp, chain), set()))
            uniq = 0 if match_df.empty else match_df[
                (match_df.cell_group == grp) & (match_df.chain == chain)]["CDR3"].nunique()
            pct = (100.0 * uniq / den) if den else 0.0
            print(f"  {grp:5s} {chain}: {uniq:4d} / {den:4d}  ({pct:5.2f}%)")


def cdr3_sharing_s10_vs_michels(s10, michels_t1d):
    """CDR3 sharing of TRA/TRB between S10 (Treg/Tconv) and Michels (T1D-only)."""
    denom, s10_rows = _s10_cdr3_by_group(s10)
    ref_sets, annot = _michels_cdr3_sets(michels_t1d)

    matches = []
    for _, r in s10_rows.iterrows():
        chain, cdr3 = r["chain"], r["CDR3"]
        if cdr3 in ref_sets[chain]:
            ags = sorted(annot[chain][cdr3])
            matches.append({
                **r,
                "Michels_antigens": "; ".join(sorted({a for a, _ in ags if a})),
                "Michels_epitopes": "; ".join(sorted({e for _, e in ags if e})),
            })
    match_df = pd.DataFrame(matches)
    print("\nS10 vs Michels (T1D-only) -- unique CDR3 shared (% of S10):")
    _print_match_summary(match_df, denom)
    return match_df


def cdr3_sharing_s10_vs_maki(s10, maki):
    """CDR3 sharing of TRA/TRB between S10 (Treg/Tconv) and Maki (in-frame primary chains)."""
    denom, s10_rows = _s10_cdr3_by_group(s10)
    ref_sets = _maki_cdr3_sets(maki)

    matches = []
    for _, r in s10_rows.iterrows():
        chain, cdr3 = r["chain"], r["CDR3"]
        if cdr3 in ref_sets[chain]:
            matches.append(dict(r))
    match_df = pd.DataFrame(matches)
    print("\nS10 vs Maki -- unique CDR3 shared (% of S10):")
    _print_match_summary(match_df, denom)
    return match_df


# ============================================================================
# Main
# ============================================================================
def main():
    # --- Task 1: remove iNKT/MAIT from all three datasets --------------------
    print("=== Task 1: remove invariant iNKT/MAIT TCRs ===")
    s10, s10_removed = remove_invariant_s10(S10_FILE)
    michels, michels_removed = remove_invariant_michels(MICHELS_FILE)
    maki, maki_removed = remove_invariant_maki(MAKI_FILE)

    export_xlsx(f"{OUT_DIR}/S10_iNKT-MAIT_filtered.xlsx", s10)
    export_xlsx(f"{OUT_DIR}/S10_iNKT-MAIT_removed_rows.xlsx", s10_removed)
    export_xlsx(f"{OUT_DIR}/Michels_iNKT-MAIT_filtered.xlsx", michels)
    export_xlsx(f"{OUT_DIR}/Michels_iNKT-MAIT_removed_rows.xlsx", michels_removed)
    export_xlsx(f"{OUT_DIR}/Maki_iNKT-MAIT_filtered.xlsx", maki)
    export_xlsx(f"{OUT_DIR}/Maki_iNKT-MAIT_removed_rows.xlsx", maki_removed)

    # --- Task 2: restrict Michels to T1D-associated antigens -----------------
    print("\n=== Task 2: restrict Michels to T1D antigens ===")
    michels_t1d = restrict_michels_t1d(michels)
    export_xlsx(f"{OUT_DIR}/Michels_iNKT-MAIT_filtered_T1Donly.xlsx", michels_t1d)

    # --- Task 3: Treg vs Tconv overlap within S10 -----------------------------
    print("\n=== Task 3: Treg vs Tconv overlap within S10 (NT; V-CDR3-J) ===")
    overlap = s10_treg_tconv_overlap(s10)
    export_sheets(f"{OUT_DIR}/S10_Treg_vs_Tconv_overlap.xlsx", {
        "Identical NT - TRA": overlap["NT_TRA"],
        "Identical NT - TRB": overlap["NT_TRB"],
        "Shared V_CDR3_J - TRA": overlap["VCDJ_TRA"],
        "Shared V_CDR3_J - TRB": overlap["VCDJ_TRB"],
    })

    # --- Task 4: CDR3 sharing S10 vs Michels (T1D-only) -----------------------
    print("\n=== Task 4: CDR3 sharing, S10 vs Michels (T1D-only) ===")
    m_matches = cdr3_sharing_s10_vs_michels(s10, michels_t1d)
    export_sheets(f"{OUT_DIR}/CDR3_sharing_S10_vs_Michels_T1Donly.xlsx",
                  {"CDR3 matches": m_matches})

    # --- Task 5: CDR3 sharing S10 vs Maki -------------------------------------
    print("\n=== Task 5: CDR3 sharing, S10 vs Maki ===")
    k_matches = cdr3_sharing_s10_vs_maki(s10, maki)
    export_sheets(f"{OUT_DIR}/CDR3_sharing_S10_vs_Maki.xlsx",
                  {"CDR3 matches": k_matches})


if __name__ == "__main__":
    try:
        main()
    except ValueError as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)
