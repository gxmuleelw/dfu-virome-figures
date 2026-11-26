"""Describe Figure 5 from existing annotations, reviewed split tables and a finished tree.

Run python -m main.figure5_analysis --help from the repository root. Each command
reads externally supplied inputs and writes a JSON summary to standard output.
"""
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path


def tsv_rows(path, columns):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(columns) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns {sorted(missing)}")
        yield from reader


def marker_span(annotations, target, first, last):
    if first < 1 or first > last:
        raise ValueError("Expected an inclusive range of positive gene indices")
    positions = {}
    columns = ["vOTU", "genomad_gene_index", "start_1based", "end_1based"]
    for row in tsv_rows(annotations, columns):
        if row["vOTU"] != target:
            continue
        index = int(row["genomad_gene_index"])
        if first <= index <= last:
            start, end = int(row["start_1based"]), int(row["end_1based"])
            if index in positions or not 1 <= start <= end:
                raise ValueError(f"Duplicate index or invalid coordinate at index {index}")
            positions[index] = (start, end)
    if len(positions) != last - first + 1:
        raise ValueError("Incomplete marker range for the supplied target")
    return {"target": target, "marker_count": len(positions),
            "start_1based": min(x[0] for x in positions.values()),
            "end_1based": max(x[1] for x in positions.values())}


def tip_categories(labels):
    categories = {}
    for row in tsv_rows(labels, ["sequence_id", "category"]):
        tip, category = row["sequence_id"], row["category"]
        if not tip or not category or tip in categories:
            raise ValueError("Each tip needs one nonempty identifier and category")
        categories[tip] = category
    return categories


def marker_consistency(splits, membership, labels, target, jobs):
    if not jobs or len(set(jobs)) != len(jobs):
        raise ValueError("Supply distinct marker jobs")
    tips = {job: set() for job in jobs}
    for row in tsv_rows(membership, ["target", "fasta_file", "sequence_id"]):
        job, tip = Path(row["fasta_file"]).stem, row["sequence_id"]
        if row["target"] == target and job in tips:
            if not tip or tip in tips[job]:
                raise ValueError(f"Empty or repeated tip in {job}")
            tips[job].add(tip)
    common = set.intersection(*tips.values())
    if target not in common:
        raise ValueError("The target must occur in each selected marker's membership")
    sides, supports = {}, {}
    columns = ["job", "target_side_rank", "target_side_members", "target_side_size",
               "SH_aLRT", "UFBoot"]
    for row in tsv_rows(splits, columns):
        job = row["job"]
        if job not in tips or int(row["target_side_rank"]) != 1:
            continue
        members = row["target_side_members"].split(";")
        side = set(members)
        if (job in sides or target not in side or not side <= common
                or len(side) != len(members) or len(side) != int(row["target_side_size"])):
            raise ValueError(f"Inconsistent first-ranked target side for {job}")
        sh, uf = float(row["SH_aLRT"]), float(row["UFBoot"])
        if not (80 <= sh <= 100 and 95 <= uf <= 100):
            raise ValueError(f"First-ranked split lacks the required support for {job}")
        sides[job] = side
        supports[job] = {"target_side_size": len(side), "SH_aLRT": sh, "UFBoot": uf}
    if sides.keys() != tips.keys():
        raise ValueError("Missing first-ranked target side for a requested marker")
    same = all(side == sides[jobs[0]] for side in sides.values())
    counts = None
    if same:
        category = tip_categories(labels)
        side = sides[jobs[0]]
        if not side <= category.keys() or category.get(target) != "Target":
            raise ValueError("The common target side requires matching tip annotations")
        counts = dict(sorted(Counter(category[tip] for tip in side).items()))
    return {"common_tip_count": len(common), "reviewed_marker_supports": supports,
            "identical_first_ranked_sides": same, "shared_side_category_counts": counts}


def tree_neighbourhood(newick, labels, target, nearest, anchors):
    from Bio import Phylo
    from table_analysis.tree_summary import neighbourhood

    tree = Phylo.read(newick, "newick")
    names = [tip.name for tip in tree.get_terminals()]
    category = tip_categories(labels)
    if (len(names) != len(set(names)) or set(names) != category.keys()
            or category.get(target) != "Target" or min(nearest, anchors) < 0):
        raise ValueError("Check tree tips, label identifiers, target and display counts")
    for clade in tree.find_clades():
        if clade is tree.root:
            continue
        length = clade.branch_length
        if length is None or not math.isfinite(length) or length < 0:
            raise ValueError("A distance summary requires finite, nonnegative branch lengths")
    _, rows = neighbourhood(tree, target, {k: {"category": v} for k, v in category.items()},
                            nearest, anchors)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="panel", required=True)
    a = commands.add_parser("span", help="5a inclusive span of a contiguous annotated marker range")
    a.add_argument("annotations", help="coordinate mapping TSV")
    a.add_argument("target")
    a.add_argument("first", type=int)
    a.add_argument("last", type=int)
    b = commands.add_parser("markers", help="5b shared membership and supplied split supports")
    b.add_argument("splits", help="reviewed target-side TSV with ranked supported splits")
    b.add_argument("membership", help="marker sequence membership TSV")
    b.add_argument("labels", help="tip labels TSV from one selected marker")
    b.add_argument("target")
    b.add_argument("jobs", nargs="+", help="marker job names matching the membership filenames' stems")
    c = commands.add_parser("neighbourhood", help="5c patristic distances and display tip selection")
    c.add_argument("newick", help="finished tree with branch lengths")
    c.add_argument("labels", help="matching tip labels TSV")
    c.add_argument("target")
    c.add_argument("--nearest", type=int, default=20)
    c.add_argument("--anchors", type=int, default=3)
    args = vars(parser.parse_args())
    command = args.pop("panel")
    result = {"span": marker_span, "markers": marker_consistency,
              "neighbourhood": tree_neighbourhood}[command](**args)
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == "__main__":
    if not __package__:
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    main()
