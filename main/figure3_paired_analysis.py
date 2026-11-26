"""Figure 3 co-abundance analysis from patient-aligned viral and MAG tables.

Call analyze() with pandas inputs. Host_genus_correlations, Member_correlations
and Cross_pairing_exclusion supply the figure renderer. The matrix tables
retain every tested pair for inspection.
"""

import numpy as np
import pandas as pd
from scipy.stats import rankdata
def unit_ranks(x, axis=-1):
    r = rankdata(x, axis=axis, method='average')
    r -= np.mean(r, axis=axis, keepdims=True)
    norm = np.sqrt(np.sum(r*r, axis=axis, keepdims=True))
    return np.divide(r, norm, out=np.full_like(r, np.nan), where=norm > 0)

def rho(x, y):
    return float(np.sum(unit_ranks(x)*unit_ranks(y)))

def bh(p, size):
    p = np.array(p, float)
    out = np.full(len(p), np.nan)
    order = np.flatnonzero(np.isfinite(p))
    order = order[np.argsort(p[order], kind='stable')]
    if len(order) > size:
        raise ValueError("More finite P values than the family size")
    if len(order):
        out[order] = np.minimum(1, np.minimum.accumulate(
            (p[order]*size/np.arange(1, len(order)+1))[::-1])[::-1])
    return out

def permutation(x, y, B, seed):
    """Rows are matched pairs, columns patients; shared bacterial permutations."""
    xr, yr = unit_ranks(x), unit_ranks(y)
    obs = np.sum(xr*yr, axis=1)
    if not np.isfinite(obs).all():
        raise ValueError("A selected member or genus has an undefined correlation")
    rng = np.random.default_rng(seed)
    count = np.zeros(len(x), dtype=np.int64)
    for start in range(0, B, 500):
        indices = np.array([rng.permutation(x.shape[1]) for _ in range(min(500, B-start))])
        simulated = np.einsum('in,ibn->ib', xr, yr[:, indices], optimize=False)
        count += (np.abs(simulated) >= np.abs(obs[:, None])-1e-12).sum(axis=1)
    return obs, (count+1)/(B+1), count

def boot_ci(x, y, B, seed):
    rng = np.random.default_rng(seed)
    values = []
    for start in range(0, B, 500):
        indices = rng.integers(0, x.shape[1], (min(500, B-start), x.shape[1]))
        values.append(np.sum(unit_ranks(x[:, indices])*unit_ranks(y[:, indices]), axis=2))
    values = np.concatenate(values, axis=1)
    valid = np.isfinite(values).sum(axis=1)
    lo, hi = np.nanquantile(values, [.025, .975], axis=1)
    lo[valid < .95*B] = np.nan
    hi[valid < .95*B] = np.nan
    return lo, hi, valid


def _check_matrix(frame, name, counts=False):
    if frame.empty or not frame.index.is_unique or not frame.columns.is_unique:
        raise ValueError(f"{name} needs unique feature and patient labels and at least one entry")
    if frame.index.hasnans or frame.columns.hasnans:
        raise ValueError(f"{name} contains missing labels")
    values = frame.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"{name} must contain finite nonnegative values")
    if counts and (values != np.floor(values)).any():
        raise ValueError(f"{name} requires integer read counts")


def _selection(rpkm, counts, hosts, bacteria, samples, minimum):
    totals = rpkm.loc[:, samples].sum(axis=0)
    if (totals <= 0).any():
        raise ValueError("Retained viral catalogue has zero total for a paired patient")
    rel = rpkm.loc[:, samples].div(totals)
    genera, members = [], []
    viral_groups = {}
    assigned = hosts[~hosts.isin(["Unassigned", "Ambiguous"])]
    for genus, group in assigned.groupby(assigned, sort=True):
        ids = group.index
        positive = int(counts.loc[ids].gt(0).any(axis=0).sum())
        if positive < minimum:
            continue
        x = rel.loc[ids].sum(axis=0)
        viral_groups[genus] = x
        y = bacteria.loc[genus, samples] if genus in bacteria.index else None
        bacterial_positive = int(y.gt(0).sum()) if y is not None else np.nan
        if y is None:
            status = "no_exact_local_genus"
        elif bacterial_positive < minimum:
            status = "low_bacterial_detection"
        elif x.nunique() <= 1 or y.nunique() <= 1:
            status = "constant_vector"
        else:
            status = "estimated"
            for member in sorted(ids):
                support = int(rel.loc[member].gt(0).sum())
                if support >= minimum:
                    members.append((member, genus, support))
        genera.append((genus, len(ids), positive, bacterial_positive, status,
                       rho(x, y) if y is not None else np.nan))
    genus_table = pd.DataFrame(genera, columns=[
        "host_genus", "n_viral_members", "viral_positive", "bacterial_positive",
        "status", "rho"])
    member_table = pd.DataFrame(members, columns=["vOTU", "host_genus", "n_positive"])
    return rel, pd.DataFrame(viral_groups).T, genus_table, member_table


def _cross_matrix(rel, bacteria, genera, members, B, seed):
    """Cross each eligible member with each estimable bacterial genus."""
    columns = ["bacterium", "vOTU", "predicted_host", "is_predicted_host",
               "n", "rho", "p_permutation", "q_BH_matrix", "exceedances"]
    if members.empty:
        return pd.DataFrame(columns=columns)
    xr = unit_ranks(rel.loc[members.vOTU].to_numpy())
    yr = unit_ranks(bacteria.loc[genera, rel.columns].to_numpy())
    observed = np.einsum("vn,gn->gv", xr, yr, optimize=False)
    if not np.isfinite(observed).all():
        raise ValueError("A selected matrix vector has an undefined correlation")
    rng = np.random.default_rng(seed)
    exceed = np.zeros(observed.shape, dtype=np.int64)
    for start in range(0, B, 500):
        indices = np.array([rng.permutation(rel.shape[1]) for _ in range(min(500, B-start))])
        simulated = np.einsum("vn,gbn->gvb", xr, yr[:, indices], optimize=False)
        exceed += (np.abs(simulated) >= np.abs(observed[:, :, None])-1e-12).sum(axis=2)
    p = (exceed.ravel()+1)/(B+1)
    rows = np.repeat(genera, len(members))
    predicted = np.tile(members.host_genus, len(genera))
    return pd.DataFrame(dict(bacterium=rows, vOTU=np.tile(members.vOTU, len(genera)),
                             predicted_host=predicted, is_predicted_host=rows == predicted,
                             n=rel.shape[1], rho=observed.ravel(), p_permutation=p,
                             q_BH_matrix=bh(p, len(p)), exceedances=exceed.ravel()))


def _cross_counts(matrix, scenario):
    result = {"scenario": scenario, "matrix_BH_family_n": len(matrix)}
    for label, flag in [("matched", True), ("nonmatched", False)]:
        subset = matrix[matrix.is_predicted_host.eq(flag)]
        hits = int(subset.q_BH_matrix.lt(.05).sum())
        result[f"matrix_{label}_n"] = len(subset)
        result[f"matrix_{label}_hits"] = hits
        result[f"matrix_{label}_percent"] = 100*hits/len(subset) if len(subset) else np.nan
    return result


def analyze(rpkm, counts, host_groups, mag_counts, mag_genus, *, seed, bootstrap_seed,
            excluded=None, min_positive=9, permutations=99999, bootstrap_replicates=9999):
    """Compute matched genus and member tests, and matched/nonmatched matrix counts.

    rpkm and counts are breadth-filtered vOTU-by-patient DataFrames for the
    full catalogue. host_groups is the indexed Series of unique highest-score
    genus assignments at score >=90, with Unassigned and Ambiguous labels.
    mag_counts contains all catalogue MAGs by patient. mag_genus is indexed by
    MAG, with the exact genus label or an empty/missing label for unnamed MAGs.

    min_positive counts patients. All patients contribute to viral detection.
    Only patients with positive total viral RPKM enter the correlations.
    Supply the permutation and bootstrap seeds. When supplied, excluded
    lists the vOTU IDs meeting the CheckV host_genes > viral_genes criterion.
    The second cross matrix uses the original paired patients and bacterial
    denominator, recomputes viral fractions and reapplies support filters.
    """
    for frame, name, require_counts in [(rpkm, "rpkm", False), (counts, "counts", True),
                                        (mag_counts, "mag_counts", True)]:
        _check_matrix(frame, name, require_counts)
    if (set(counts.index) != set(rpkm.index) or set(counts.columns) != set(rpkm.columns)
            or set(mag_counts.columns) != set(rpkm.columns)):
        raise ValueError("Feature or patient labels disagree between input matrices")
    for series, index in [(host_groups, rpkm.index), (mag_genus, mag_counts.index)]:
        if not series.index.is_unique or set(series.index) != set(index):
            raise ValueError("Annotation labels must match their feature matrix")
    if host_groups.isna().any() or host_groups.eq("").any():
        raise ValueError("Use explicit Unassigned or Ambiguous host labels")
    for n in [min_positive, permutations, bootstrap_replicates]:
        if not isinstance(n, (int, np.integer)) or n < 1:
            raise ValueError("Support and resampling counts must be positive integers")
    counts = counts.loc[rpkm.index, rpkm.columns]
    if not np.array_equal(counts.to_numpy() > 0, rpkm.to_numpy() > 0):
        raise ValueError("Filtered counts and RPKM disagree on detection")
    hosts = host_groups.loc[rpkm.index]
    mag = mag_counts.loc[:, rpkm.columns]
    tax = mag_genus.loc[mag.index].fillna("")
    mag_total = mag.sum(axis=0)
    if (mag_total <= 0).any():
        raise ValueError("Every patient requires a positive MAG count total")
    named = tax.ne("")
    # The matched and matrix analyses use these respective summation orders.
    bacteria = mag.div(mag_total).loc[named].groupby(tax[named]).sum()
    cross_bacteria = mag.loc[named].groupby(tax[named]).sum().div(mag_total)
    samples = rpkm.columns[rpkm.sum(axis=0).gt(0)]
    if len(samples) < 2:
        raise ValueError("At least two patients need a defined viral composition")
    rel, viral_groups, genus, member = _selection(rpkm, counts, hosts, bacteria, samples, min_positive)
    eligible = genus.status.eq("estimated")
    genera = genus.loc[eligible, "host_genus"].tolist()
    for col in ["p_permutation", "q_BH_genus", "permutation_exceedances",
                "ci_lower", "ci_upper", "bootstrap_valid"]:
        genus[col] = np.nan
    for col in ["rho", "p_permutation", "permutation_exceedances", "q_BH_members"]:
        member[col] = np.nan
    genus["family_n"] = len(genus)
    member["member_family_n"] = len(member)
    if genera:
        x = np.vstack([viral_groups.loc[genera].to_numpy(), rel.loc[member.vOTU].to_numpy()])
        y = bacteria.loc[genera + member.host_genus.tolist(), samples].to_numpy()
        obs, p, count = permutation(x, y, permutations, seed)
        n = len(genera)
        lo, hi, valid = boot_ci(x[:n], y[:n], bootstrap_replicates, bootstrap_seed)
        genus.loc[eligible, ["p_permutation", "permutation_exceedances", "ci_lower", "ci_upper",
                             "bootstrap_valid"]] = np.column_stack((p[:n], count[:n], lo, hi, valid))
        genus["q_BH_genus"] = bh(genus.p_permutation, len(genus))
        member["rho"], member["p_permutation"] = obs[n:], p[n:]
        member["permutation_exceedances"], member["q_BH_members"] = count[n:], bh(p[n:], len(member))
    matrix = _cross_matrix(rel, cross_bacteria, genera, member, permutations, seed)
    summary = [_cross_counts(matrix, "original")]
    results = {"Host_genus_correlations": genus, "Member_correlations": member,
               "Original_matrix": matrix}
    if excluded is not None:
        excluded = list(excluded)
        if len(excluded) != len(set(excluded)) or set(excluded).difference(rpkm.index):
            raise ValueError("Excluded vOTUs must be unique catalogue IDs")
        retained = rpkm.drop(index=excluded)
        reduced, _, selected, tested = _selection(retained, counts.loc[retained.index],
                                                   hosts.loc[retained.index], cross_bacteria,
                                                   samples, min_positive)
        genera = selected.loc[selected.status.eq("estimated"), "host_genus"].tolist()
        sensitivity = _cross_matrix(reduced, cross_bacteria, genera, tested, permutations, seed)
        results["Exclusion_matrix"] = sensitivity
        summary.append(_cross_counts(sensitivity, "exclusion"))
    results["Cross_pairing_exclusion"] = pd.DataFrame(summary)
    return results
