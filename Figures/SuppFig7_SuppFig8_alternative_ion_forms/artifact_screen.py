"""
Pairwise offset + retention-time co-elution screen for spectral-library artifacts.

Idea: in a library built by collecting every spectrum that satisfies a diagnostic-ion
query (e.g. acylcarnitine head-group ions), some entries are not novel modifications but
analytical artifacts of OTHER entries in the same library: isotope peaks, alkali/alkaline-
earth adducts, and in-source fragments. Each such artifact is related to a "parent" entry
in the library by a known mass offset AND co-elutes with it (same file, same RT), because
it is physically the same molecule entering the source together.

This module does NOT need a designated parent. It self-joins the library: for every pair of
entries in the SAME file whose precursor masses differ by a known artifact offset (within a
ppm tolerance) AND whose retention times match (within an RT tolerance), it flags the
artifact member of the pair.

Conventions
-----------
- For each candidate offset we treat the LOWER-mass entry as `i` and look for a HIGHER-mass
  partner `j` at  precmz_i + offset.
- For isotopes and adducts, the HIGHER-mass member is the artifact (parent is the lighter ion).
- For in-source-fragment neutral losses, the LOWER-mass member is the artifact (parent is the
  heavier, intact ion). This is encoded per offset via `artifact_member`.
"""

import numpy as np
import pandas as pd


# (name, exact mass offset in Da, which member of the pair is the artifact)
# Higher member = artifact  -> isotopes & adducts.
ISOTOPE_ADDUCT_OFFSETS = [
    ("13C isotope (M+1)",   1.003355, "higher"),
    ("13C isotope (M+1, 2+)", 0.501678, "higher"),
    ("Na adduct (Na-H)",   21.981944, "higher"),
    ("Ca adduct (Ca-2H)",  37.946941, "higher"),
    ("K adduct (K-H)",     37.955882, "higher"),
    ("NH4 adduct (NH4-H)", 17.026549, "higher"),
]

# Optional: common in-source neutral losses. Lower member = artifact (the fragment).
# These are heuristic; rigorous ISF detection is better done via MS/MS fragment matching.
# Append to the offset list if you want them included.
ISF_NEUTRAL_LOSS_OFFSETS = [
    ("ISF: -H2O",                 18.010565, "lower"),
    ("ISF: -2H2O",                36.021130, "lower"),
    ("ISF: -3H2O",                54.031695, "lower"),
    ("ISF: -4H2O",                72.042260, "lower"),
    ("ISF: -CO",                  27.994915, "lower"),
    ("ISF: -CO2",                 43.989830, "lower"),
]


def screen_artifacts(df,
                     offsets=ISOTOPE_ADDUCT_OFFSETS,
                     ppm_tol=20.0,
                     rt_tol_sec=10.0,
                     mz_col="precmz",
                     rt_col="rt",
                     file_col="original_path",
                     scan_col="original_scan",
                     rt_in_minutes=True):
    """
    Screen a library table for offset + co-elution artifact pairs.

    Returns
    -------
    pairs : DataFrame, one row per flagged co-eluting pair, columns:
        original_path, offset_name, offset_da, observed_delta, ppm_error, rt_diff_sec,
        parent_idx, parent_scan, parent_precmz, parent_rt,
        artifact_idx, artifact_scan, artifact_precmz, artifact_rt
    annotated : copy of `df` with added columns:
        is_artifact (bool), artifact_types (str, ';'-joined), n_pairs (int),
        n_files_flagged (int),
        parent_scans (str, ';'-joined original_scan(s) of the parent entry/entries
                      that this row was classified against)
    `parent_*`/`artifact_*` indices refer to df's original index.
    """
    rt_scale = 60.0 if rt_in_minutes else 1.0
    records = []

    for file_path, g in df.groupby(file_col, sort=False):
        idx = g.index.to_numpy()                       # original df indices
        mz = g[mz_col].to_numpy(dtype=float)
        rt = g[rt_col].to_numpy(dtype=float) * rt_scale  # -> seconds
        scan = g[scan_col].to_numpy()

        order = np.argsort(mz)
        mz_sorted = mz[order]

        for name, off, artifact_member in offsets:
            targets = mz + off                          # expected higher partner for each row
            tol = targets * ppm_tol * 1e-6              # ppm tol referenced to the heavier ion
            lo = np.searchsorted(mz_sorted, targets - tol, side="left")
            hi = np.searchsorted(mz_sorted, targets + tol, side="right")

            for i in range(len(mz)):
                for k in range(lo[i], hi[i]):
                    j = order[k]                        # j is the higher-mass partner
                    if j == i:
                        continue
                    if abs(rt[i] - rt[j]) > rt_tol_sec:
                        continue

                    observed = mz[j] - mz[i]
                    ppm_err = (observed - off) / targets[i] * 1e6

                    # i = lower mass, j = higher mass
                    if artifact_member == "higher":
                        a, p = j, i                     # artifact = higher, parent = lower
                    else:
                        a, p = i, j                     # artifact = lower,  parent = higher

                    records.append({
                        "original_path":   file_path,
                        "offset_name":     name,
                        "offset_da":       off,
                        "observed_delta":  observed,
                        "ppm_error":       ppm_err,
                        "rt_diff_sec":     abs(rt[i] - rt[j]),
                        "parent_idx":      idx[p],
                        "parent_scan":     scan[p],
                        "parent_precmz":   mz[p],
                        "parent_rt":       rt[p] / rt_scale,
                        "artifact_idx":    idx[a],
                        "artifact_scan":   scan[a],
                        "artifact_precmz": mz[a],
                        "artifact_rt":     rt[a] / rt_scale,
                    })

    pairs = pd.DataFrame.from_records(records)

    annotated = df.copy()
    annotated["is_artifact"] = False
    annotated["artifact_types"] = ""
    annotated["n_pairs"] = 0
    annotated["n_files_flagged"] = 0
    annotated["parent_scans"] = ""

    if not pairs.empty:
        grp = pairs.groupby("artifact_idx")
        types = grp["offset_name"].apply(lambda s: ";".join(sorted(set(s))))
        npairs = grp.size()
        nfiles = grp["original_path"].nunique()
        parent_scans = grp["parent_scan"].apply(
            lambda s: ";".join(str(v) for v in dict.fromkeys(s))  # unique, order-preserving
        )

        annotated.loc[types.index, "artifact_types"] = types
        annotated.loc[npairs.index, "n_pairs"] = npairs
        annotated.loc[nfiles.index, "n_files_flagged"] = nfiles
        annotated.loc[parent_scans.index, "parent_scans"] = parent_scans
        annotated.loc[types.index, "is_artifact"] = True

    return pairs, annotated


# ----------------------------------------------------------------------------- #
# Self-test on synthetic data
# ----------------------------------------------------------------------------- #
if __name__ == "__main__":
    rng = np.random.default_rng(0)

    rows = []
    # File A: a parent acylcarnitine at m/z 400.300, rt 5.00 min, with a Na adduct and an M+1
    rows.append(dict(precmz=400.3000, rt=5.000, original_path="fileA.mzML", original_scan=101))
    rows.append(dict(precmz=400.3000 + 21.981944 + 400.3 * 5e-6,  # +0.5 ppm off, co-elutes
                     rt=5.001, original_path="fileA.mzML", original_scan=102))  # Na adduct
    rows.append(dict(precmz=400.3000 + 1.003355,
                     rt=4.998, original_path="fileA.mzML", original_scan=103))  # M+1 isotope

    # A genuine analog: +14.0157 (CH2), different RT -> must NOT be flagged
    rows.append(dict(precmz=400.3000 + 14.0157,
                     rt=6.200, original_path="fileA.mzML", original_scan=104))

    # A Na-offset partner but at a DIFFERENT RT -> must NOT be flagged (not co-eluting)
    rows.append(dict(precmz=500.4000, rt=7.00, original_path="fileA.mzML", original_scan=105))
    rows.append(dict(precmz=500.4000 + 21.981944,
                     rt=7.50, original_path="fileA.mzML", original_scan=106))  # +0.5 min = 30 s

    # File B: same parent appears again with a K adduct co-eluting
    rows.append(dict(precmz=400.3000, rt=5.10, original_path="fileB.mzML", original_scan=201))
    rows.append(dict(precmz=400.3000 + 37.955882,
                     rt=5.11, original_path="fileB.mzML", original_scan=202))  # K adduct

    df = pd.DataFrame(rows)

    pairs, annotated = screen_artifacts(df, ppm_tol=20, rt_tol_sec=10)

    pd.set_option("display.width", 200, "display.max_columns", 30)
    print("=== FLAGGED PAIRS ===")
    print(pairs[["original_path", "offset_name", "ppm_error", "rt_diff_sec",
                 "parent_scan", "artifact_scan"]].to_string(index=False))
    print("\n=== PER-ROW ANNOTATION ===")
    print(annotated[["original_path", "original_scan", "precmz", "rt",
                     "is_artifact", "artifact_types", "parent_scans",
                     "n_files_flagged"]].to_string(index=False))

    # Assertions encoding the expected behavior
    flagged_scans = set(annotated.loc[annotated.is_artifact, "original_scan"])
    assert 102 in flagged_scans, "Na adduct should be flagged"
    assert 103 in flagged_scans, "M+1 isotope should be flagged"
    assert 202 in flagged_scans, "K adduct (file B) should be flagged"
    assert 104 not in flagged_scans, "Genuine +CH2 analog must NOT be flagged"
    assert 106 not in flagged_scans, "Na-offset at 30 s apart must NOT be flagged (no co-elution)"
    print("\nAll assertions passed.")