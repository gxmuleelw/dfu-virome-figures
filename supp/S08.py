from supp_style import *

def label(fig, ax, letter, title):
    pos = ax.get_position()
    fig.text(pos.x0, pos.y1 + .025, letter, fontsize=11, fontweight='bold')
    fig.text(pos.x0 + .035, pos.y1 + .025, title, fontsize=9, fontweight='bold')
def draw(models, corr, loo, scenarios, names, highlighted_ids):
    loo=loo.sort_values('delta_rho').reset_index(drop=True)
    fig=canvas(188)
    severity=fig.add_axes([.20,.64,.32,.28]);rho=fig.add_axes([.69,.64,.27,.28])
    for i,scenario in enumerate(scenarios):
        for delta,test,c,m in [(-.13,'positive_abundance_severity',GREY,'o'),(.13,'host_adjusted_severity',BLUE,'s')]:
            r=models[(models.variant.eq(scenario))&(models.test.eq(test))].iloc[0]
            severity.errorbar(r.estimate,i+delta,xerr=[[r.estimate-r.lower],[r.upper-r.estimate]],fmt=m,color=c,markersize=4,elinewidth=1,capsize=2)
        r=corr[corr.scenario.eq(scenario)].iloc[0]
        rho.errorbar(r.rho,i,xerr=[[r.rho-r.lower],[r.upper-r.rho]],fmt='o',color=GREY,markersize=4,elinewidth=1,capsize=2)
        rho.text(.75,i-.22,f'{r.rho:.3f}',fontsize=7.5)
    severity.set_yticks(range(4),names);severity.set_ylim(3.5,-.5);severity.set_xlim(-1.2,6);severity.set_xticks([0,2,4,6])
    severity.axvline(0,color='#B9C3C8',ls=':',lw=.7);severity.set_xlabel('Severity difference (log2 fraction)')
    label(fig,severity,'a','Escherichia severity effect')
    rho.set_ylim(3.5,-.5);rho.set_yticks([]);rho.set_xlim(.73,1.01);rho.set_xticks([.75,.85,.95]);rho.set_xlabel('Spearman ρ (95% CI)')
    label(fig,rho,'b','Host co-abundance')
    handles=[Line2D([],[],marker=m,color=c,linestyle='',markersize=4,label=t) for m,c,t in [('o',GREY,'Age + sex'),('s',BLUE,'+ E. coli')]]
    severity.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,-.29),ncol=2,frameon=False,fontsize=7.5)

    ax=fig.add_axes([.20,.10,.76,.38])
    colors=[BLUE if x in highlighted_ids else GREY for x in loo.removed_vOTU]
    ax.hlines(np.arange(len(loo)),0,loo.delta_rho,color='#CFD8DD',lw=.7)
    ax.scatter(loo.delta_rho,np.arange(len(loo)),c=colors,s=15,zorder=3)
    ax.set_yticks(np.arange(len(loo)),loo.removed_vOTU,fontsize=7.5);ax.invert_yaxis();ax.set_ylim(len(loo)-.2,-.8)
    ax.axvline(0,color='#B9C3C8',ls=':',lw=.8);ax.tick_params(axis='y',length=0)
    ax.set_xlim(-.0075,.0045);ax.set_xticks([-.006,-.003,0,.003]);ax.set_xlabel('Change in ρ after removing one member')
    for tick,name in zip(ax.get_yticklabels(),loo.removed_vOTU):
        if name in highlighted_ids:tick.set_color(BLUE);tick.set_fontweight('bold')
    label(fig,ax,'c','Leave-one-member-out sensitivity')
    return fig
