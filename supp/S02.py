from supp_style import *

def label_panel(ax, letter, title):
    # Panel letters sit in the margin, never on top of the data.
    ax.text(-.14, 1.12, letter, transform=ax.transAxes, fontsize=11,
            fontweight='bold', va='bottom')
    ax.set_title(title, loc='left', pad=9, fontweight='bold', fontsize=9)

def draw(loo, contribution, original_estimate):
    loo=loo.sort_values('beta').reset_index(drop=True)
    contribution=contribution.sort_values('mean_share_pct',ascending=False).reset_index(drop=True)
    fig=canvas(186)
    ax=fig.add_axes([.13,.601,.80,.253]);x=np.arange(1,len(loo)+1)
    ax.vlines(x,loo.ci_low,loo.ci_high,color='#b7ceda',lw=.4)
    ax.scatter(x,loo.beta,s=5,color=BLUE,zorder=3)
    ax.axhline(0,color='#8a8a8a',lw=.7,ls='--')
    ax.axhline(float(original_estimate),color=ORANGE,lw=1,ls='--')
    ax.set_xlim(0,len(loo)+1);ax.set_ylim(min(loo.ci_low)-1,2)
    ax.set_ylabel('Wagner 4 minus 1–3 (pp)')
    ax.set_xlabel('Excluded vOTU, ranked by the saved point estimate')
    label_panel(ax,'a',f"{len(loo)} leave-one-out results")
    ax=fig.add_axes([.22,.183,.70,.253])
    top=contribution.head(12).iloc[::-1]
    colors=[ORANGE if x=='Escherichia' else BLUE for x in top.predicted_host_genus]
    ax.barh(np.arange(len(top)),top.mean_share_pct,color=colors,height=.65)
    ax.set_yticks(np.arange(len(top)),top.vOTU,fontsize=7.5)
    ax.set_xlim(0,float(top.mean_share_pct.max())*1.22)
    for i,(_,r) in enumerate(top.iterrows()):
        ax.text(float(r.mean_share_pct)+.12,i,f'{r.mean_share_pct:.2f}%',fontsize=7.3,va='center')
    ax.set_xlabel('Mean share of predicted-lytic abundance (%)')
    label_panel(ax,'b','Largest 12 member contributions')
    return fig
