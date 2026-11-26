from supp_style import *

def draw(rows):
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig = plt.figure(figsize=(176 / 25.4, 115 / 25.4))
    ax = fig.add_axes([0.04, 0.25, 0.92, 0.4])
    ax.axis('off')
    headers = ['Marker', 'Tree tips', 'Trimmed AA', 'Model', 'SH-aLRT (%)', 'UFBoot (%)']
    values = [[r['job'].split('_')[1], r['n_sequences'], r['trimmed_sites'], r['selected_model'], f"{float(r['SH_aLRT']):.1f}", f"{float(r['UFBoot']):.0f}"] for r in rows]
    table = ax.table(cellText=values, colLabels=headers, colWidths=[0.12, 0.12, 0.15, 0.2, 0.22, 0.19], cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor('#CDD3D8')
        cell.set_linewidth(0.5)
        cell.set_facecolor('#EDF1F4' if r == 0 else 'white')
        if r == 0:
            cell.set_text_props(weight='bold')
    return fig
