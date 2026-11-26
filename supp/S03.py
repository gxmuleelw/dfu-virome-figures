from supp_style import *

def label_panel(ax, letter, title):
    # Panel letters sit in the margin, never on top of the data.
    ax.text(-.14, 1.12, letter, transform=ax.transAxes, fontsize=11,
            fontweight='bold', va='bottom')
    ax.set_title(title, loc='left', pad=9, fontweight='bold', fontsize=9)
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
def draw(S6, SCENARIOS, SC_LABELS):
    CATS=['lytic','temperate','unclassified']
    CAT_LABELS=['Predicted lytic','Temperate-associated','Unclassified']
    def row(scenario, category='lytic', model='age_sex'):
        d = S6[(S6.scenario==scenario)&(S6.category==category)&(S6.model==model)]
        assert len(d)==1
        return d.iloc[0]
    fig=canvas(192)
    ys=np.arange(3,-1,-1)
    for j,cat in enumerate(CATS):
        ax=fig.add_axes([.29+j*.225,.583,.19,.273]);rr=[row(s,cat) for s in SCENARIOS]
        forest_marks(ax,rr,ys,color=[BLUE,GREY,'#8b7693'][j])
        ax.set_ylim(-.6,3.7);ax.set_xlim(-38,40);ax.set_xticks([-30,0,30])
        ax.set_xlabel('Difference (pp)')
        label_panel(ax,chr(97+j),CAT_LABELS[j].replace('Temperate-associated','Temperate-\nassociated'))
        for r,y in zip(rr,ys):
            ax.text(.5,(y+.6)/4.3-.09,f'q={float(r.BHq):.4f}',transform=ax.transAxes,fontsize=7.3,ha='center')
    for label,y in zip(SC_LABELS,ys):
        fig.text(.07,.583+(y+.6)/4.3*.273,label,fontsize=7.5,va='center')
    fig.text(.07,.496,'d',fontsize=11,fontweight='bold')
    fig.text(.105,.496,'Predicted lytic after additional E. coli adjustment',fontsize=9,fontweight='bold')
    ax=fig.add_axes([.33,.212,.47,.239]);rr=[row(s,model='age_sex_ecoli_z') for s in SCENARIOS]
    forest_marks(ax,rr,ys);ax.set_xlim(-38,12);ax.set_ylim(-.5,3.5)
    ax.set_xlabel('Wagner 4 minus 1–3 (percentage points)');ax.set_xticks([-30,-20,-10,0,10])
    for label,r,y in zip(SC_LABELS,rr,ys):
        fy=.212+(y+.5)/4*.239
        fig.text(.07,fy,label,fontsize=7.5,va='center')
        fig.text(.83,fy,f'P={float(r.raw_p):.4f}',fontsize=7.5,va='center')
    return fig
