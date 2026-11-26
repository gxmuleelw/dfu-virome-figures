"""Calculate Figure 1b–f summaries and write input JSON for figure1_reflow.py.

Inputs are the retained catalogue's CheckV annotations, vConTACT2 overview,
geNomad MSL39 rank assignments, unmasked mapping matrices and sample groups.
Matrices are TSVs with vOTUs in rows and sample IDs as column names. Breadth
is a fraction. All quantification cells must be observed. An observed sample
with no confirmed detection stays in the detection denominator; its relative
composition is undefined. Panel a is a study workflow drawn separately.
"""

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np


def read_table(path, key, fields=(), delimiter="\t"):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not {key, *fields} <= set(header):
            raise ValueError(f"Missing required columns or duplicate headers in {path.name}")
        records = {}
        for row in reader:
            if None in row or None in row.values() or not row[key].strip() or row[key] in records:
                raise ValueError(f"Incomplete, oversized or duplicate record in {path.name}")
            records[row[key]] = row
        if not records:
            raise ValueError(f"Empty input table: {path.name}")
    return header, records


def read_matrix(path):
    header, rows = read_table(path, "vOTU")
    if header[0] != "vOTU" or len(header) < 2 or any(not s.strip() for s in header[1:]):
        raise ValueError(f"Expected vOTU followed by sample IDs in {path.name}")
    samples = header[1:]
    matrix = np.array([[row[s] for s in samples] for row in rows.values()], dtype=float)
    if not np.isfinite(matrix).all() or (matrix < 0).any():
        raise ValueError(f"Missing, non-finite or negative measurement in {path.name}")
    return list(rows), samples, matrix


def category_counts(rows, field, levels):
    counts = Counter(row[field] for row in rows.values())
    if not set(counts) <= set(levels):
        raise ValueError(f"Unknown or missing {field}")
    return [counts[level] for level in levels]


def summarize(args):
    _, info = read_table(args.quality, "vOTU", ["checkv_quality"])
    _, vc = read_table(args.vc, "Genome", ["VC Status", "VC"], delimiter=",")
    _, taxonomy = read_table(args.taxonomy, "vOTU", ["assigned_rank_MSL39"])
    ids, samples, raw = read_matrix(args.counts)
    rpkm_ids, rpkm_samples, raw_rpkm = read_matrix(args.rpkm)
    breadth_ids, breadth_samples, breadth = read_matrix(args.breadth)
    _, meta = read_table(args.metadata, "sample", ["severity"])
    if not (ids == rpkm_ids == breadth_ids and samples == rpkm_samples == breadth_samples):
        raise ValueError("Count, RPKM and breadth matrix labels must match in order")
    if set(samples) != set(meta) or set(ids) != set(info) or set(ids) != set(taxonomy) or not set(ids) <= set(vc):
        raise ValueError("Sample or catalogue identifiers differ between inputs")
    if (raw != np.floor(raw)).any() or (breadth > 1).any():
        raise ValueError("Expected integer counts and breadth fractions from zero to one")

    # b–d: CheckV's quality category, study members of shared VCs, deepest rank.
    quality = category_counts(info, "checkv_quality", ["Complete", "High-quality"])
    ranks = category_counts(taxonomy, "assigned_rank_MSL39", ["class", "order", "family", "unclassified"])
    clustered = []
    for vid in ids:
        status = vc[vid]["VC Status"]
        if status.strip().lower() in ("", "na", "nan"):
            raise ValueError("Missing VC Status for a catalogue member")
        if status == "Clustered":
            cluster = vc[vid]["VC"].strip()
            if cluster.lower() in ("", "na", "nan"):
                raise ValueError("A clustered catalogue member requires a VC identifier")
            clustered.append(cluster)

    # e: count positive, breadth-confirmed observations over all participants.
    confirmed = breadth >= 0.50
    detected = (raw > 0) & confirmed
    positives = detected.sum(axis=1)
    n = len(samples)
    histogram = np.bincount(positives, minlength=n + 1)
    majority = int((positives > n / 2).sum())

    # f: the ANCOM-BC2 prevalence threshold is based on the composition cohort.
    composition = detected.any(axis=0)
    if not np.array_equal(composition, ((raw_rpkm > 0) & confirmed).any(axis=0)):
        raise ValueError("Confirmed count and RPKM totals disagree on sample eligibility")
    if not composition.any():
        raise ValueError("Composition is undefined for all samples")
    minimum = math.ceil(0.20 * int(composition.sum()))
    retained = detected[:, composition].sum(axis=1) >= minimum
    ancom = composition & detected[retained].any(axis=0)
    groups = ("Wagner_1_3", "Wagner_4")
    severity = np.array([meta[s]["severity"] for s in samples])
    if not set(severity) <= set(groups):
        raise ValueError("Expected Wagner_1_3 or Wagner_4 in severity")
    cohorts = [[int((severity[mask] == g).sum()) for g in groups]
               for mask in (np.ones(n, dtype=bool), composition, ancom)]

    return {
        "quality": quality, "clusters": [len(clustered), len(ids) - len(clustered)], "ranks": ranks,
        "cluster_label": f"{len(set(clustered))} VCs", "cohorts": cohorts,
        "frequency": histogram[1:].tolist(), "undetected_votus": int(histogram[0]),
        "frequency_note": f"{histogram[1]} detected once\n{majority} detected in >50%",
        "detection_xticks": sorted({1, *range(10, n - 5, 10), n}),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, description in [
        ("quality", "TSV with vOTU, checkv_quality"),
        ("vc", "vConTACT2 CSV with Genome, VC Status, VC; may contain reference genomes"),
        ("taxonomy", "TSV with vOTU, assigned_rank_MSL39"),
        ("counts", "Unmasked mapped counts TSV"), ("rpkm", "Unmasked RPKM TSV"),
        ("breadth", "Observed breadth fraction TSV"), ("metadata", "TSV with sample, severity"),
        ("output", "Output JSON for the quantitative panels"),
    ]:
        parser.add_argument("--" + name, required=True, type=Path, help=description)
    args = parser.parse_args()
    result = summarize(args)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
