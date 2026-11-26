"""Bar chart and count labels for the catalogue summaries in Figure 1b–d."""


def bars(ax, labs, vals, cols, ylim, yticks, ylabel, label_offset):
    positions = range(len(vals))
    ax.bar(positions, vals, width=0.57, color=cols, linewidth=0)
    ax.set_xticks(positions, labs)
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_ylabel(ylabel, labelpad=1.5)
    ax.grid(axis="y", color="#E5ECEF", lw=0.45)
    ax.set_axisbelow(True)
    for index, value in enumerate(vals):
        ax.text(index, value + label_offset, str(value), ha="center", fontsize=7.2)
