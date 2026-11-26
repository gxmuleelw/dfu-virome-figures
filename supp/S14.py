from supp_style import *

from matplotlib.colors import to_rgba
def draw(model,patients,group_labels):
    spec=[('ARG carriers','CARD_overall','b',40),('AMG carriers','AMG total','a',70),('VF carriers','VFDB_overall','c',40)]
    fig=canvas(180)
    gs=fig.add_gridspec(2,3,left=.10,right=.98,bottom=.12,top=.84,
                        height_ratios=[1,1.25],hspace=.78,wspace=.55)
    rng=np.random.default_rng(20261001)
    for i,(title,feature,oldpanel,limit) in enumerate(spec):
        ax=fig.add_subplot(gs[0,i])
        values=patients.loc[patients.panel==oldpanel].sort_values('sample')
        row=model.loc[(model.family_id=='function_totals')&(model.feature==feature)].iloc[0]
        for pos,(group,color,marker) in enumerate([('Wagner_1_3',BLUE,'o'),('Wagner_4',ORANGE,'s')],1):
            y=values.loc[values.severity==group,'relative_abundance_percent'].astype(float).to_numpy()
            bp=ax.boxplot([y],positions=[pos],widths=.48,whis=1.5,showfliers=False,patch_artist=True,
                          medianprops={'color':'black','linewidth':.9},
                          boxprops={'edgecolor':color,'linewidth':.8},
                          whiskerprops={'color':color,'linewidth':.7},capprops={'color':color,'linewidth':.7})
            bp['boxes'][0].set_facecolor(to_rgba(color,.13))
            ax.scatter(pos+rng.uniform(-.15,.15,len(y)),y,s=12,c=color,marker=marker,
                       linewidth=.25,edgecolors='white',alpha=.9,zorder=3,clip_on=False)
        ax.set_xlim(.55,2.45);ax.set_ylim(0,limit)
        ax.set_yticks([0,20,40,60] if limit==70 else [0,10,20,30,40])
        ax.set_xticks([1,2],group_labels)
        ax.set_ylabel('Relative abundance (%)')
        ax.grid(axis='y',color='#E8ECEF',linewidth=.5);ax.set_axisbelow(True)
        ax.text(-.24,1.23,'abc'[i],transform=ax.transAxes,fontsize=11,fontweight='bold')
        ax.set_title(title,loc='left',fontsize=9,fontweight='bold',pad=35)
        ax.text(0,1.035,f'Model P = {float(row.raw_p):.4f}\nBH q = {float(row.bh_q_recalculated):.4f}',
                transform=ax.transAxes,fontsize=7.5,va='bottom')

    cats=model.loc[model.family_id=='amg_categories'].reset_index(drop=True)
    sub=gs[1,:].subgridspec(1,3,width_ratios=[1.9,2.2,1.0],wspace=.06)
    labels=fig.add_subplot(sub[0,0]);forest=fig.add_subplot(sub[0,1]);numbers=fig.add_subplot(sub[0,2])
    for ax in (labels,forest,numbers):ax.set_ylim(len(cats)-.4,-.6)
    labels.set_axis_off();numbers.set_axis_off()
    labels.text(-.15,1.08,'d',transform=labels.transAxes,fontsize=11,fontweight='bold')
    labels.text(0,1.08,'AMG functional categories',transform=labels.transAxes,fontsize=9,fontweight='bold')
    numbers.text(.5,1.02,'P / BH q',transform=numbers.transAxes,ha='center',fontsize=7.5,fontweight='bold')
    forest.axvline(0,color=GREY,ls='--',lw=.8,zorder=1)
    for i,row in cats.iterrows():
        label=row.feature.replace('Nucleotide / information systems','Nucleotide / information\nsystems')
        labels.text(0,i,label,va='center',fontsize=7.5)
        est,lo,hi=map(float,[row.estimate,row.ci95_lower,row.ci95_upper])
        forest.errorbar(est,i,xerr=[[est-lo],[hi-est]],fmt='o',color=GREY,markersize=4,
                        elinewidth=1,capsize=2.5,zorder=3)
        numbers.text(.5,i,f'{float(row.raw_p):.4f} / {float(row.bh_q_recalculated):.4f}',
                     ha='center',va='center',fontsize=7.5)
    forest.set_xlim(-14,25);forest.set_xticks([-10,0,10,20]);forest.set_yticks([])
    forest.set_xlabel('Adjusted difference (percentage points)\nWagner 4 minus Wagner 1–3',fontsize=7.5)
    forest.spines['left'].set_visible(False)
    return fig
