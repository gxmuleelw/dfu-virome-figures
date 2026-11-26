from collections import Counter
from Bio.Phylo.BaseTree import Clade
from matplotlib.patches import Polygon
import math
from supp_style import *

STYLES = {'Target': ('#D55E00', '*'), 'ICTV reference': ('#0072B2', 'D'), 'Public, not VMR': ('#666666', 's'), 'MGV': ('#009E73', '^'), 'GPD': ('#A74F87', 'o')}
def descendants(tree):
    result = {}
    for c in tree.find_clades(order='postorder'):
        result[c] = frozenset([c.name]) if c.is_terminal() else frozenset().union(*(result[x] for x in c.clades))
    return result
def numeric_support(raw, order):
    if not raw:
        return {}
    parts = raw.split('/')
    if len(parts) != 2:
        return {}
    return dict(zip(order, map(float, parts)))
def draw(tree, target, labels, order, title, view, chosen=None, caution=''):
    desc = descendants(tree)
    depth = tree.depths()
    all_tips = tree.get_terminals()
    maximal = max(depth.values())
    extent = max(maximal, 1e-08)
    visible = {}
    collapsed = set()
    sibling_groups = {}
    leaves = []
    flattened_zero_edges = set()

    def display_children(c):
        children = []
        for child in c.clades:
            if chosen is not None and (not child.is_terminal()) and ((child.branch_length or 0) == 0) and (not child.name):
                flattened_zero_edges.add(child)
                children.extend(display_children(child))
            else:
                children.append(child)
        return children

    def visit(c):
        if chosen is not None and (not desc[c] & chosen):
            collapsed.add(c)
            leaves.append(c)
            return
        if c.is_terminal():
            leaves.append(c)
            return
        children = display_children(c)
        omitted = [child for child in children if chosen is not None and (not desc[child] & chosen)]
        if len(omitted) > 1:
            group = Clade()
            desc[group] = frozenset().union(*(desc[x] for x in omitted))
            depth[group] = depth[c]
            sibling_groups[group] = len(omitted)
            children = [x for x in children if x not in omitted] + [group]
        visible[c] = children
        for child in children:
            if child in sibling_groups:
                collapsed.add(child)
                leaves.append(child)
            else:
                visit(child)
    visit(tree.root)
    ypos = {c: float(i) for i, c in enumerate(leaves)}

    def calc_y(c):
        if c in ypos:
            return ypos[c]
        child_y = [calc_y(child) for child in visible[c]]
        ypos[c] = (child_y[0] + child_y[-1]) / 2
        return ypos[c]
    calc_y(tree.root)
    h = max(7, len(leaves) * (0.15 if view == 'full' else 0.23) + 2.0)
    fig, ax = plt.subplots(figsize=(17, h))
    fig.subplots_adjust(left=0.045, right=0.985, bottom=1.0 / h, top=1 - (1.22 if caution else 0.95) / h)
    label_x = maximal + extent * 0.025
    terminal_labels = []
    for c, children in visible.items():
        ys = [ypos[x] for x in children]
        ax.plot([depth[c], depth[c]], [min(ys), max(ys)], color='#555555', lw=0.65, zorder=1)
        for child in children:
            color = '#D55E00' if child.is_terminal() and child.name == target else '#555555'
            ax.plot([depth[c], depth[child]], [ypos[child]] * 2, color=color, lw=1.4 if color == '#D55E00' else 0.65, zorder=1)
    for c in leaves:
        y = ypos[c]
        if c in collapsed:
            counts = Counter((labels[n]['category'] for n in desc[c]))
            if c in sibling_groups:
                x2 = depth[c]
                ax.scatter([x2], [y], marker='s', s=24, facecolors='white', edgecolors='#9BA4AE', linewidths=0.7, zorder=3)
                text = f'Omitted at node: {len(desc[c])} references in {sibling_groups[c]} sister subtrees'
            else:
                x2 = max((depth[t] for t in all_tips if t.name in desc[c]))
                ax.add_patch(Polygon([(depth[c], y), (x2, y - 0.32), (x2, y + 0.32)], closed=True, facecolor='#D8DDE2', edgecolor='#9BA4AE', lw=0.6))
                text = f'Collapsed: {len(desc[c])} references'
            text += ' (' + ', '.join((f'{k}={v}' for k, v in sorted(counts.items()))) + ')'
            xend = x2
            textcolor = '#555555'
            font = 7
        else:
            meta = labels[c.name]
            color, marker = STYLES[meta['category']]
            ax.scatter([depth[c]], [y], c=color, marker=marker, s=65 if c.name == target else 18, zorder=3, linewidths=0.5)
            text = meta['label']
            xend = depth[c]
            textcolor = '#D55E00' if c.name == target else '#20252B'
            font = 8 if view != 'full' else 7
            if meta['identical_AA_group_size'] > 1:
                text += f"  [identical AA: n={meta['identical_AA_group_size']}]"
        ax.plot([xend, label_x], [y, y], color='#D2D7DC', lw=0.4, ls=':', zorder=0)
        terminal_labels.append(ax.text(label_x, y, text, va='center', ha='left', fontsize=font, color=textcolor, weight='bold' if c.is_terminal() and c.name == target else 'normal'))
    for c in visible:
        if c is tree.root:
            continue
        raw = c.name or ''
        values = numeric_support(raw, order)
        strong = values.get('SH-aLRT', -1) >= 80 and values.get('UFBoot', -1) >= 95
        if strong:
            ax.scatter([depth[c]], [ypos[c]], s=8, color='#111111', zorder=4)
        if view != 'full' and raw:
            ax.annotate(raw, (depth[c], ypos[c]), xytext=(-2, 4), textcoords='offset points', fontsize=6.5, ha='right', va='bottom', color='#353A40', bbox={'fc': 'white', 'ec': 'none', 'alpha': 0.82, 'pad': 0.1})
    ax.set_xlim(-extent * 0.018, maximal + extent * 1.25)
    ax.set_ylim(len(leaves) - 0.4, -0.6)
    ax.set_yticks([])
    ax.set_xlabel('Branch length (amino-acid substitutions per site)', fontsize=8)
    for side in ['left', 'right', 'top']:
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_bounds(0, maximal)
    ticks = [x for x in ax.get_xticks() if 0 <= x <= maximal]
    ax.set_xticks(ticks)
    if title: fig.text(0.045, 1 - 0.24 / h, title, fontsize=13, weight='bold', va='top')
    legend = [Line2D([0], [0], marker=m, color='none', markerfacecolor=c, markeredgecolor=c, label=k, markersize=7) for k, (c, m) in STYLES.items() if any((r['category'] == k for r in labels.values()))]
    legend.append(Line2D([0], [0], marker='o', color='none', markerfacecolor='black', label='SH-aLRT >=80 and UFBoot >=95', markersize=3))
    fig.legend(handles=legend, loc='upper left', bbox_to_anchor=(0.04, 1 - 0.48 / h), ncol=4, frameon=False, fontsize=8)
    if caution:
        fig.text(0.045, 1 - 0.94 / h, caution, fontsize=8, color='#A64020', weight='bold', va='top')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [t.get_window_extent(renderer) for t in terminal_labels]
    assert all((b.width > 0 and b.height > 0 for b in boxes))
    assert all((not boxes[i].overlaps(boxes[i + 1]) for i in range(len(boxes) - 1))), 'Overlapping tip labels.'
    return fig
