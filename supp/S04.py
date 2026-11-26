from supp_style import *

GROUPS=['Wagner_1_3','Wagner_4']
COLORS=[BLUE,ORANGE]
MARKERS=['o','s']
def label_panel(ax, letter, title):
    # Panel letters sit in the margin, never on top of the data.
    ax.text(-.14, 1.12, letter, transform=ax.transAxes, fontsize=11,
            fontweight='bold', va='bottom')
    ax.set_title(title, loc='left', pad=9, fontweight='bold', fontsize=9)
def group_legend(fig, group_labels, y=.96, x=.52):
    h=[Line2D([],[],marker=m,color='none',markerfacecolor=c,markeredgecolor='white',
              markeredgewidth=.35,markersize=4.8,label=l)
       for m,c,l in zip(MARKERS,COLORS,group_labels)]
    fig.legend(handles=h, loc='center', bbox_to_anchor=(x,y), frameon=False,
               ncol=2, columnspacing=2, handletextpad=.4, fontsize=8)
def patient_box(ax, df, value, ylim=None):
    rng=np.random.default_rng(20261112)
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
def draw(alpha, ar, ordinations, group_labels):
    fig=canvas(187)
    group_legend(fig,group_labels,y=.937)


    ax=fig.add_axes([.115,.602,.34,.251]);patient_box(ax,alpha,'Shannon',(.7,5.3))
    ax.set_ylabel('Shannon diversity');label_panel(ax,'a','Within-sample diversity')
    ax.text(.02,.97,f'Adjusted difference {ar.estimate:.3f}\n95% CI {ar.lower:.3f} to {ar.upper:.3f}\nP={ar.p:.4f}',
            transform=ax.transAxes,fontsize=7.5,va='top',bbox={'facecolor':'white','edgecolor':'none','alpha':.8,'pad':1})
    records=[]
    for idd,(name,title,letter,pos) in enumerate([
        ('Figure3b','Bray–Curtis','b',[.605,.602,.34,.251]),
        ('Figure3c','Robust Aitchison','c',[.115,.187,.34,.251])]):
        d, perm, disp, geom = ordinations[idd]
        pr=perm[perm.term=='severity'].iloc[0]
        dp=disp[disp.term=='Groups'].iloc[0]
        ax=fig.add_axes(pos)
        for g,c,m in zip(GROUPS,COLORS,MARKERS):
            sub=d[d.severity==g]
            ax.scatter(sub.PCoA1,sub.PCoA2,s=16,c=c,marker=m,edgecolors='white',linewidths=.3,zorder=3)
            el=next(v for v in geom if v['group']==g)
            ax.add_patch(Ellipse(el['center'],2*el['semiaxes'][0],2*el['semiaxes'][1],
                                 angle=el['angle_degrees'],facecolor='none',edgecolor=c,lw=.85))
        ax.axhline(0,color='#d5d9dc',lw=.5,zorder=0);ax.axvline(0,color='#d5d9dc',lw=.5,zorder=0)
        ax.autoscale_view();ax.margins(.08)
        ax.set_xlabel(f'PCoA 1 ({d.axis1_percent.iloc[0]:.1f}%)')
        ax.set_ylabel(f'PCoA 2 ({d.axis2_percent.iloc[0]:.1f}%)')
        label_panel(ax,letter,title)
        ax.text(.02,.97,f'R²={pr.R2:.4f}; P={pr.p:.4f}\nPERMDISP P={dp.p:.4f}',transform=ax.transAxes,
                fontsize=7.5,va='top',bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':1})
        records.append(dict(distance=title,R2=float(pr.R2),P=float(pr.p),PERMDISP_P=float(dp.p),axis1=float(d.axis1_percent.iloc[0]),axis2=float(d.axis2_percent.iloc[0])))
    ax=fig.add_axes([.605,.187,.34,.251]);label_panel(ax,'d','Severity-associated variation')
    for i,r in enumerate(records):
        y=1-i
        ax.barh(y,100*r['R2'],color=[GREY,BLUE][i],height=.36)
        ax.text(100*r['R2']+.2,y,f'{100*r["R2"]:.2f}%\nP={r["P"]:.4f}',fontsize=7.5,va='center')
    ax.set_yticks([1,0],['Bray–Curtis','Robust\nAitchison'])
    ax.set_ylim(-.65,1.65);ax.set_xlim(0,9);ax.set_xlabel('Severity R² (%)')
    ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    return fig
