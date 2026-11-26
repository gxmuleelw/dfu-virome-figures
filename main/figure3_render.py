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
def numeric(data, columns):
    data = data.copy()
    for column in columns:
        data[column] = pd.to_numeric(data[column], errors='coerce')
    return data


def label(fig, ax, letter, title):
    pos = ax.get_position()
    fig.text(pos.x0, pos.y1 + .025, letter, fontsize=11, fontweight='bold')
    fig.text(pos.x0 + .035, pos.y1 + .025, title, fontsize=9, fontweight='bold')


def qtext(q):
    return f'{q:.2g}' if q < .001 else f'{q:.4f}'


def host_overview():
    genus = numeric(table(7, 'Host_genus_correlations'), ['rho','ci_lower','ci_upper','q_BH_genus'])
    genus = genus[genus.status.eq('estimated')].sort_values('rho', ascending=False).reset_index(drop=True)
    member = numeric(table(7, 'Member_correlations'), ['rho','q_BH_members']).sort_values('rho').reset_index(drop=True)
    cross = table(7, 'Cross_pairing_exclusion')
    qcounts = [int((member.q_BH_members < .05).sum()), int(((member.q_BH_members >= .05)&(member.q_BH_members < .25)).sum()), int((member.q_BH_members >= .25).sum())]

    fig = canvas(193)
    left = fig.add_axes([.24,.115,.275,.80])
    upper = fig.add_axes([.66,.63,.305,.285])
    lower = fig.add_axes([.66,.22,.305,.26])
    ypos = np.arange(len(genus))
    for y, row in genus.iterrows():
        focal = row.host_genus in ['Escherichia','Collinsella','Parabacteroides']
        color = BLUE if focal else GREY
        left.errorbar(row.rho,y,xerr=[[row.rho-row.ci_lower],[row.ci_upper-row.rho]],fmt='o',color=color,
                      markersize=3.4 if focal else 2.7,elinewidth=.8,capsize=1.5)
    left.set_yticks(ypos,genus.host_genus,fontsize=7.5,fontstyle='italic')
    left.invert_yaxis(); left.set_ylim(len(genus)-.2,-.8)
    for tick, name in zip(left.get_yticklabels(),genus.host_genus):
        if name in ['Escherichia','Collinsella','Parabacteroides']:
            tick.set_color(BLUE);tick.set_fontweight('bold')
    left.set_xlim(-.5,1.03); left.set_xticks([-.5,0,.5,1]);left.axvline(0,color='#B9C3C8',ls=':',lw=.7,zorder=0)
    left.set_xlabel('Spearman ρ (95% CI)')
    left.tick_params(axis='y',length=0)
    fig.text(.035,.94,'a',fontsize=11,fontweight='bold')
    fig.text(.073,.94,'Host-genus collections',fontsize=9,fontweight='bold')

    colors = np.where(member.q_BH_members < .05, BLUE, np.where(member.q_BH_members < .25,GREY,'#BFC9CE'))
    upper.scatter(np.arange(1,len(member)+1),member.rho,c=colors,s=12,linewidths=0)
    upper.axhline(0,color='#B9C3C8',ls=':',lw=.7,zorder=0)
    upper.set_xlim(-2,len(member)+3);upper.set_ylim(-.35,1.03);upper.set_xticks(cfg['member_ticks'])
    upper.set_xlabel('vOTUs ordered by ρ');upper.set_ylabel('Spearman ρ')
    label(fig,upper,'b','Individual vOTUs')
    handles=[Line2D([],[],marker='o',linestyle='',color=c,markersize=4,label=t) for c,t in
             [(BLUE,f'q < 0.05: {qcounts[0]}'),(GREY,f'0.05 ≤ q < 0.25: {qcounts[1]}'),('#BFC9CE',f'q ≥ 0.25: {qcounts[2]}')]]
    upper.legend(handles=handles,loc='lower right',frameon=False,handletextpad=.4,borderpad=.2,fontsize=7.5)

    xpos=np.array([0,1,3,4]); vals=[]; counts=[]
    for _, row in cross.iterrows():
        vals.extend([row.matrix_matched_percent,row.matrix_nonmatched_percent])
        counts.extend([(int(row.matrix_matched_hits),int(row.matrix_matched_n)),(int(row.matrix_nonmatched_hits),int(row.matrix_nonmatched_n))])
    lower.bar(xpos,vals,width=.7,color=[BLUE,'#AEBEC7',BLUE,'#AEBEC7'])
    for x,y,(a,b) in zip(xpos,vals,counts):
        lower.text(x,y+3,f'{a}/{b}\n{y:.1f}%',ha='center',va='bottom',fontsize=7.5)
    lower.set_ylim(0,104);lower.set_xlim(-.65,4.65)
    lower.set_xticks(xpos,['Match','Other','Match','Other'],fontsize=7.5)
    lower.set_yticks([0,25,50,75,100]);lower.set_ylabel('Pairs with q < 0.05 (%)')
    label(fig,lower,'c','Host specificity')
    lower.text(.215,-.19,'Original',transform=lower.transAxes,ha='center',fontsize=8)
    lower.text(.78,-.19,cfg['scenario_label'],transform=lower.transAxes,ha='center',fontsize=8)
    fig.text(.66,.115,'Match: predicted host genus\nOther: nonmatched bacterial genus',fontsize=7.5,linespacing=1.4)
    return fig

if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('data',type=Path); p.add_argument('config',type=Path); p.add_argument('output',type=Path); a=p.parse_args(); data=a.data; cfg=json.loads(a.config.read_text(encoding='utf-8-sig')); a.output.mkdir(parents=True,exist_ok=True)
    fig=host_overview()
    for ext in ['pdf','svg','png']: fig.savefig(a.output/('Figure_3.'+ext),dpi=600)
