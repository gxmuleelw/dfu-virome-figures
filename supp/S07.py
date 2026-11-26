from supp_style import *

import matplotlib
def draw(rows, tick_indices):
    hosts = {}
    for row in rows:
        name, host = row["vOTU"], row["predicted_host"]
        assert hosts.setdefault(name, host) == host
    bacteria = sorted({r["bacterium"] for r in rows})
    row_index = {name: i for i, name in enumerate(bacteria)}
    viruses = sorted(hosts, key=lambda v: (row_index[hosts[v]], int(v.removeprefix("vOTU"))))
    column_index = {v: j for j, v in enumerate(viruses)}
    rho = np.full((len(bacteria), len(viruses)), np.nan)
    q = np.full_like(rho, np.nan)
    for row in rows:
        i, j = row_index[row["bacterium"]], column_index[row["vOTU"]]
        assert np.isnan(rho[i, j])
        assert (row["is_predicted_host"] == "True") == (row["bacterium"] == row["predicted_host"])
        rho[i, j], q[i, j] = float(row["rho"]), float(row["q_BH_matrix"])
    assert np.isfinite(rho).all() and np.isfinite(q).all()
    assert (np.abs(rho) <= 1).all() and ((0 <= q) & (q <= 1)).all()
    groups = []
    for name in bacteria:
        positions = [j for j, v in enumerate(viruses) if hosts[v] == name]
        if positions:
            assert positions == list(range(min(positions), max(positions)+1))
            groups.append((name, min(positions), max(positions)+1))

    plt.rcParams.update({"font.family": ["Arial", "DejaVu Sans"], "font.size": 7.3,
                         "svg.fonttype": "none", "pdf.fonttype": 42,
                         "ps.fonttype": 42, "axes.unicode_minus": False})
    fig = plt.figure(figsize=(176/25.4, 180/25.4))
    handles = [Rectangle((0, 0), 1, 1, facecolor="none", edgecolor="black", lw=0.7,
                         label="Predicted-host match"),
               Line2D([], [], marker="o", color="black", markerfacecolor="white", markersize=3,
                      markeredgewidth=0.4, linestyle="none", label="Matrix q < 0.05")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.018, 0.922),
               frameon=False, ncol=2, handlelength=1, handletextpad=0.5,
               columnspacing=1.5, fontsize=7.3)
    ax = fig.add_axes([0.31, 0.19, 0.65, 0.686])
    mesh = ax.pcolormesh(np.arange(len(viruses)+1), np.arange(len(bacteria)+1), rho, cmap=plt.get_cmap("RdBu_r", 257),
                         vmin=-1, vmax=1, shading="flat", rasterized=False)
    ax.set(xlim=(0, len(viruses)), ylim=(len(bacteria), 0), yticks=np.arange(len(bacteria))+0.5,
           yticklabels=[f"{i+1:02d}  {name}" for i, name in enumerate(bacteria)],
           xticks=[v-0.5 for v in tick_indices],
           xticklabels=[str(v) for v in tick_indices],
           xlabel="Viral member column index")
    ax.tick_params(axis="y", length=0, pad=4, labelsize=7.3)
    ax.tick_params(axis="x", length=2.5, width=0.5, labelsize=7.3)
    for name, label in zip(bacteria, ax.get_yticklabels()):
        if name in ("Escherichia", "Collinsella"):
            label.set_color("#0072B2" if name == "Escherichia" else "#E69F00")
            label.set_fontweight("bold")
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    for name, start, end in groups:
        i = row_index[name]
        # Dual strokes keep the original predicted-host cells legible on either side of the color scale.
        for edge, width in (("white", 0.8), ("black", 0.35)):
            ax.add_patch(Rectangle((start, i), end-start, 1, fill=False, edgecolor=edge, lw=width, zorder=4))
        if start:
            ax.axvline(start, color="white", lw=0.25, alpha=0.65, zorder=2)
    significant_i, significant_j = np.where(q < 0.05)
    ax.scatter(significant_j+0.5, significant_i+0.5, s=0.9, facecolors="white",
               edgecolors="#202020", linewidths=0.2, zorder=3)
    cax = fig.add_axes([0.47, 0.065, 0.35, 0.014])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal", ticks=[-1, -0.5, 0, 0.5, 1])
    cb.solids.set_rasterized(False)
    cb.set_label("Spearman rho", fontsize=7.5, labelpad=3)
    cb.ax.tick_params(length=2.5, width=0.5, labelsize=7.3)
    cb.outline.set_linewidth(0.5)
    return fig
