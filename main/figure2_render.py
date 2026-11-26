from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
BLUE, ORANGE, GREY, MAGENTA = '#0072B2', '#E69F00', '#697D88', '#CC79A7'

plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
 'axes.titlesize':9,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.65,
 'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':2.5,
 'ytick.major.size':2.5,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
 'axes.unicode_minus':False,'savefig.facecolor':'white'})

def canvas(height_mm):
    return plt.figure(figsize=(175/25.4,height_mm/25.4))

def table(_,name):
    return pd.read_csv(data/(name+'.tsv'),sep='\t')
GROUPS = ['Wagner_1_3', 'Wagner_4']
COLORS = [BLUE, ORANGE]
MARKERS = ['o', 's']
CATS = ['lytic', 'temperate', 'unclassified']
CAT_LABELS = ['Predicted lytic', 'Temperate-associated', 'Unclassified']
def label_panel(ax, letter, title):
    # Panel letters sit in the margin, never on top of the data.
    ax.text(-.14, 1.12, letter, transform=ax.transAxes, fontsize=11,
            fontweight='bold', va='bottom')
    ax.set_title(title, loc='left', pad=9, fontweight='bold', fontsize=9)

def group_legend(fig, y=.96, x=.52):
    h=[Line2D([],[],marker=m,color='none',markerfacecolor=c,markeredgecolor='white',
              markeredgewidth=.35,markersize=4.8,label=l)
       for m,c,l in zip(MARKERS,COLORS,cfg['group_labels'])]
    fig.legend(handles=h, loc='center', bbox_to_anchor=(x,y), frameon=False,
               ncol=2, columnspacing=2, handletextpad=.4, fontsize=8)

def patient_box(ax, df, value, ylim=None):
    rng=np.random.default_rng(cfg['jitter_seed'])
    arrays=[df.loc[df.severity==g,value].astype(float).to_numpy() for g in GROUPS]
    b=ax.boxplot(arrays, positions=[0,1], widths=.46, patch_artist=True,
                 showfliers=False, medianprops={'color':'#222222','linewidth':1.05},
                 boxprops={'linewidth':.8}, whiskerprops={'linewidth':.7},
                 capprops={'linewidth':.7})
    for patch,c in zip(b['boxes'], COLORS):
        patch.set(facecolor=c,alpha=.13,edgecolor=c)
    for i,(a,c,m) in enumerate(zip(arrays,COLORS,MARKERS)):
        ax.scatter(i+rng.uniform(-.15,.15,len(a)),a,s=13,c=c,marker=m,
                   edgecolors='white',linewidths=.25,zorder=3,alpha=.86)
    ax.set_xticks([0,1],['Wagner 1–3','Wagner 4'])
    ax.set_xlim(-.55,1.55)
    if ylim: ax.set_ylim(*ylim)
    ax.grid(axis='y', color='#e9ecef',lw=.4,zorder=0)
    ax.set_axisbelow(True)

def forest_marks(ax, rows, ys, color=BLUE):
    ax.axvline(0,color='#8a8a8a',lw=.7,ls='--',zorder=0)
    for r,y in zip(rows,ys):
        e,lo,hi=map(float,[r['estimate'],r['ci95_lower'],r['ci95_upper']])
        ax.plot([lo,hi],[y,y],color=color,lw=1.2)
        ax.plot([lo,lo],[y-.045,y+.045],color=color,lw=.8)
        ax.plot([hi,hi],[y-.045,y+.045],color=color,lw=.8)
        ax.scatter([e],[y],s=24,marker='D',color=color,zorder=3)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])

def main_figure():
    fig=canvas(183)
    group_legend(fig,y=.936)
    samples=[]
    for i in range(3):
        d=table(None,'samples_'+CATS[i])
        samples.append(d)
        ax=fig.add_axes([.09+i*.306,.611,.235,.25])
        patient_box(ax,d,'relative_abundance_percent',(0,105))
        label_panel(ax,chr(97+i),CAT_LABELS[i].replace('Temperate-associated','Temperate-\nassociated'))
        if i==0:ax.set_ylabel('Relative abundance (%)')
        else:ax.set_yticklabels([])
        q=float(d.model_bh_q.iloc[0])
        ax.text(.5,.965,f'q={q:.4f}',transform=ax.transAxes,ha='center',va='top',fontsize=8)
    h=table(None,'homology_counts')
    ax=fig.add_axes([.09,.211,.235,.267])
    label_panel(ax,'d','Same-genus MAG match')
    ax.bar(np.arange(3),h.percent,width=.58,color=[BLUE,GREY,'#aebbc3'])
    for j,r in h.iterrows():
        ax.text(j,r.percent+1.4,f'{r["count"]}/{r.denominator}\n{r.percent:.1f}%',
                ha='center',va='bottom',fontsize=7.5)
    ax.set_ylim(0,49);ax.set_ylabel('vOTUs with >50% coverage (%)')
    ax.set_xticks([0,1,2],['Lytic','Temperate','Unclassified'],rotation=30,ha='right')
    ax.tick_params(axis='x',labelsize=7.5)
    ax.grid(axis='y',color='#e9ecef',lw=.4);ax.set_axisbelow(True)
    saved=table(None,'sensitivity').fillna({'BHq':'NA'}); frows=saved.to_dict('records'); labels=saved.display_scenario.tolist()
    ax=fig.add_axes([.65,.205,.30,.285]);ys=np.arange(4,-1,-1)
    forest_marks(ax,frows,ys);ax.set_ylim(-.55,4.55);ax.set_xlim(-40,8)
    ax.set_xticks([-40,-20,0]);ax.set_xlabel('Wagner 4 minus 1–3\n(percentage points)',labelpad=6)
    fig.text(.38,.522,'e',fontsize=11,fontweight='bold')
    fig.text(.415,.522,'Predicted lytic: sensitivity',fontsize=9,fontweight='bold')
    for i,(r,y,label) in enumerate(zip(frows,ys,labels)):
        fy=.205+(y+.55)/5.1*.285
        q=r['BHq'];stat=f'q={float(q):.4f}' if str(q)!='NA' else f'P={float(r["raw_p"]):.4f}; q not computed'
        fig.text(.39,fy+.006,label,ha='left',va='bottom',fontsize=7.5)
        fig.text(.39,fy-.011,stat,ha='left',va='top',fontsize=7.3,color='#52616b')
    return fig

if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('data',type=Path); p.add_argument('config',type=Path); p.add_argument('output',type=Path); a=p.parse_args(); data=a.data; cfg=json.loads(a.config.read_text(encoding='utf-8-sig')); a.output.mkdir(parents=True,exist_ok=True)
    fig=main_figure()
    for ext in ['pdf','svg','png']: fig.savefig(a.output/('Figure_2.'+ext),dpi=600)
