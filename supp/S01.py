from supp_style import *

def draw(qc, catalogue, special_label, special_xy, special_label_xy):
    qc=qc.sort_values('sample').reset_index(drop=True)
    length=catalogue.length_bp.astype(float).to_numpy()
    vg=catalogue.viral_genes.astype(float).to_numpy()
    hg=catalogue.host_genes.astype(float).to_numpy()
    host_dominant=hg>vg
    fig=canvas(178)
    gs=fig.add_gridspec(2,2,left=.11,right=.97,bottom=.12,top=.855,wspace=.40,hspace=.64)
    ax=[fig.add_subplot(gs[i//2,i%2]) for i in range(4)]
    sample_x=np.arange(1,len(qc)+1)
    groups=[('Wagner_1_3',BLUE,'o'),('Wagner_4',ORANGE,'s')]
    for group,color,marker in groups:
        mask=qc.severity.eq(group).to_numpy()
        raw=qc.catalogue_mapping_fraction.astype(float).to_numpy()*100
        confirmed=qc.breadth_confirmed_mapping_fraction.astype(float).to_numpy()*100
        ax[0].vlines(sample_x[mask],confirmed[mask],raw[mask],color=color,alpha=.38,lw=.65)
        ax[0].scatter(sample_x[mask],raw[mask],facecolors='white',edgecolors=color,marker=marker,s=15,lw=.7,zorder=2)
        ax[0].scatter(sample_x[mask],confirmed[mask],c=color,marker=marker,s=15,lw=.3,zorder=3)
        ax[1].scatter(sample_x[mask],qc.loc[mask,'nonzero_vOTUs'].astype(float),c=color,marker=marker,s=19,lw=.3,zorder=3)
    for a in ax[:2]:
        a.set_xlim(0,len(qc)+1);a.set_xticks([1,*range(10,len(qc),10),len(qc)]);a.set_xlabel('Participant (ordered by ID)')
        a.grid(axis='y',color='#E8ECEF',lw=.5);a.set_axisbelow(True)
    panel(ax[0],'a','Catalogue mapping',x=-.21)
    ax[0].set_ylabel('Mapped reads (%)');ax[0].set_ylim(bottom=-1)
    ax[0].legend(handles=[Line2D([],[],marker='o',markerfacecolor='white',markeredgecolor=GREY,color=GREY,ls='',label='All mapped'),
                          Line2D([],[],marker='o',color=GREY,ls='',label='Breadth-confirmed')],
                 loc='upper right',frameon=False,fontsize=7.3,handletextpad=.4)
    panel(ax[1],'b','Per-participant detection',x=-.21)
    ax[1].set_ylabel('Detected vOTUs');ax[1].set_ylim(bottom=-2)
    ax[1].annotate(special_label,xy=special_xy,xytext=special_label_xy,arrowprops={'arrowstyle':'-','lw':.6,'color':GREY},fontsize=7.5)
    panel(ax[2],'c','Sequence length',x=-.21)
    ax[2].hist(length/1000,bins=np.arange(0,751,25),color=GREY,edgecolor='white',linewidth=.35)
    ax[2].set_yscale('log');ax[2].set_ylim(.8,600);ax[2].set_xlim(0,750)
    ax[2].set_yticks([1,10,100],['1','10','100'])
    ax[2].set_xlabel('Sequence length (kb)');ax[2].set_ylabel('vOTUs (log scale)')
    ax[2].set_xticks([0,200,400,600]);ax[2].axvline(200,ls='--',lw=.8,color=MAGENTA)
    ax[2].text(.97,.95,f"{int((length>200000).sum())} vOTUs >200 kb",transform=ax[2].transAxes,ha='right',va='top',fontsize=8)
    panel(ax[3],'d','CheckV gene evidence',x=-.21)
    for mask,color,label in [(~host_dominant,'#A7B3BA',f"host ≤ viral ({int((~host_dominant).sum())})"),(host_dominant,MAGENTA,f"host > viral ({int(host_dominant.sum())})")]:
        ax[3].scatter(vg[mask],hg[mask],s=13,c=color,edgecolors='none',alpha=.7,label=label,zorder=3)
    upper=max(hg.max(),vg.max())*1.25
    ax[3].plot([0,upper],[0,upper],color=GREY,ls='--',lw=.7,zorder=1)
    ax[3].set_xscale('symlog',linthresh=1);ax[3].set_yscale('symlog',linthresh=1)
    ax[3].set_xlim(-.1,upper);ax[3].set_ylim(-.1,upper)
    ticks=[0,1,10,100];ax[3].set_xticks(ticks,[str(t) for t in ticks]);ax[3].set_yticks(ticks,[str(t) for t in ticks])
    ax[3].set_xlabel('Viral genes (zero-preserving log scale)');ax[3].set_ylabel('Host genes (same scale)')
    handles, labels=ax[3].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.77,.017),frameon=False,
               fontsize=7.3,handletextpad=.3,ncol=2,columnspacing=.8)
    fig.legend(handles=[Line2D([],[],marker='o',color=BLUE,ls='',label='Wagner 1–3'),
                        Line2D([],[],marker='s',color=ORANGE,ls='',label='Wagner 4')],
               loc='upper center',bbox_to_anchor=(.56,.985),ncol=2,frameon=False,title='Participant panels a–b',title_fontsize=8)
    return fig
