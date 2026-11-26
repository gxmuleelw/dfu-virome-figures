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


def focal_hosts():
    genus=numeric(table(7,'Host_genus_correlations'),['rho','ci_lower','ci_upper','q_BH_genus'])
    clinical=numeric(table(8,'Target_host_comparisons'),['estimate','lower','upper','p','q'])
    ranks=table(None,'ranks')
    lines=cfg['descriptive_lines']
    fig=canvas(185)
    axes=[fig.add_axes([.12,.61,.35,.30]),fig.add_axes([.61,.61,.35,.30])]
    for ax, host, letter in zip(axes,['Escherichia','Collinsella'],['a','b']):
        part=ranks[ranks.host_genus.eq(host)]
        row=genus[genus.host_genus.eq(host)].iloc[0]
        for severity, color, marker in [('Wagner_1_3',BLUE,'o'),('Wagner_4',ORANGE,'s')]:
            sub=part[part.severity.eq(severity)]
            ax.scatter(sub.bacterial_rank,sub.viral_rank,s=15,color=color,marker=marker,linewidths=.35,edgecolor='white',zorder=3)
        line=next(x for x in lines if x['host_genus']==host)
        xs=np.array(line['x_span_observed_ranks'])
        ax.plot(xs,line['slope_for_drawing_only']*xs+line['intercept_for_drawing_only'],lw=.9,color='#8596A0',zorder=1)
        ax.set_xlim(0,44);ax.set_ylim(0,44);ax.set_xticks([0,10,20,30,40]);ax.set_yticks([0,10,20,30,40])
        ax.set_xlabel('Bacterial abundance rank');ax.set_ylabel('Viral-collection abundance rank')
        label(fig,ax,letter,host)
        stats=f'ρ = {row.rho:.3f}\n95% CI {row.ci_lower:.3f}–{row.ci_upper:.3f}\nq = {qtext(row.q_BH_genus)}'
        if host == 'Escherichia':
            ax.text(.02,.98,stats,transform=ax.transAxes,va='top',fontsize=7.5)
        else:
            ax.text(.98,.06,stats,transform=ax.transAxes,va='bottom',ha='right',fontsize=7.5)
        ax.text(.01,-.29,f'{int(row.n_viral_members)} vOTUs / {int(row.n_local_MAGs)} bacterial MAG'+('s' if row.n_local_MAGs!=1 else ''),transform=ax.transAxes,fontsize=7.5)
    handles=[Line2D([],[],marker=m,color=c,linestyle='',markersize=4,label=t) for m,c,t in [('o',BLUE,'Wagner 1–3'),('s',ORANGE,'Wagner 4')]]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.54,.53),ncol=2,frameon=False,fontsize=8)

    detect=fig.add_axes([.12,.14,.35,.28])
    for i,host in enumerate(['Escherichia','Collinsella']):
        r=clinical[(clinical.host.eq(host))&(clinical.test.eq('detection_severity'))].iloc[0]
        for x,pos,n,c in [(i-.19,int(r.positive_ref),int(r.n_ref),BLUE),(i+.19,int(r.positive_case),int(r.n_case),ORANGE)]:
            pct=pos/n*100
            detect.bar(x,pct,width=.32,color=c)
            detect.text(x,pct+2.8,f'{pos}/{n}',ha='center',fontsize=8)
        detect.text(i,107,f'q = {r.q:.3f}',ha='center',fontsize=7.5)
    detect.set_ylim(0,120);detect.set_yticks([0,25,50,75,100]);detect.set_ylabel('Patients with detected collection (%)')
    detect.set_xticks([0,1],['Escherichia','Collinsella'],fontstyle='italic')
    label(fig,detect,'c',cfg['detection_title'])

    effect=fig.add_axes([.66,.30,.30,.12])
    es=clinical[clinical.host.eq('Escherichia')]
    for y,test,c,m,title in [(1,'positive_abundance_severity',GREY,'o','Age + sex'),(0,'host_adjusted_severity',BLUE,'s','+ E. coli')]:
        r=es[es.test.eq(test)].iloc[0]
        effect.errorbar(r.estimate,y,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt=m,color=c,markersize=4,elinewidth=1,capsize=2)
        effect.text(.98,y/2+.45,f'β = {r.estimate:.2f}; q = {r.q:.4f}',transform=effect.transAxes,ha='right',fontsize=7.5)
    effect.set_yticks([1,0],['Age + sex','+ E. coli'],fontsize=7.5)
    effect.set_ylim(-.6,1.85);effect.set_xlim(-1.2,5.2);effect.set_xticks([0,2,4]);effect.axvline(0,color='#B9C3C8',lw=.7,ls=':')
    effect.set_xlabel('Severity difference (log2 fraction)')
    label(fig,effect,'d','Escherichia severity effect')

    score=fig.add_axes([.66,.105,.30,.08])
    r=es[es.test.eq('positive_abundance_score')].iloc[0]
    score.errorbar(r.estimate,0,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt='o',color=GREY,markersize=4,elinewidth=1,capsize=2)
    score.set_xlim(0,3);score.set_ylim(-.7,.7);score.set_xticks([0,1,2,3]);score.set_yticks([])
    score.axvline(0,color='#B9C3C8',lw=.7,ls=':')
    score.set_xlabel('log2 fraction / original score unit')
    score.text(.98,.92,f'β = {r.estimate:.3f}; q = {qtext(r.q)}',transform=score.transAxes,ha='right',fontsize=7.5)
    label(fig,score,'e','Bacterial composite score')
    return fig

if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('data',type=Path); p.add_argument('config',type=Path); p.add_argument('output',type=Path); a=p.parse_args(); data=a.data; cfg=json.loads(a.config.read_text(encoding='utf-8-sig')); a.output.mkdir(parents=True,exist_ok=True)
    fig=focal_hosts()
    for ext in ['pdf','svg','png']: fig.savefig(a.output/('Figure_4.'+ext),dpi=600)
