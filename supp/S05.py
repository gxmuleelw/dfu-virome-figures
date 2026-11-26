from supp_style import *

def draw(allv,c,sample_n):
    c=c.copy()
    pa=c['iPHoP host genus (score90)'].eq('Parabacteroides')
    fig=canvas(205)
    a=fig.add_axes([.12,.72,.82,.22])
    hit=allv.q<.25
    a.scatter(allv.loc[~hit,'estimate'],-np.log10(allv.loc[~hit,'q'].astype(float)),s=17,color='#BBC4CA',edgecolor='white',linewidth=.3,zorder=2)
    other=c.loc[~pa];pc=c.loc[pa]
    a.scatter(other['Effect estimate'],-np.log10(other.q.astype(float)),s=24,color=BLUE,zorder=3)
    a.scatter(pc['Effect estimate'],-np.log10(pc.q.astype(float)),s=36,marker='D',facecolors='none',edgecolors=MAGENTA,linewidth=1,zorder=4)
    a.axhline(-np.log10(.25),color=GREY,ls='--',lw=.8)
    a.axvline(0,color='#C7CCD1',lw=.7)
    a.set(xlabel='ANCOM-BC2 effect (natural-log abundance)',ylabel='−log10(BH q)',ylim=(-.035,1.5))
    a.text(.98,.93,f"{len(allv)} vOTUs; {sample_n} patients\n{len(c)} candidates at q < 0.25",transform=a.transAxes,ha='right',va='top',fontsize=8)
    a.text(.99,-np.log10(.25)+.02,'q = 0.25',transform=a.get_yaxis_transform(),ha='right',fontsize=7.5)
    a.legend(handles=[Line2D([],[],marker='D',mfc='none',mec=MAGENTA,ls='',label=f"Parabacteroides ({int(pa.sum())})"),Line2D([],[],marker='o',color=BLUE,ls='',label=f"Other candidates ({int((~pa).sum())})")],loc='upper left',frameon=False,handletextpad=.4)
    panel(a,'a','Single-vOTU severity associations',x=-.11)

    c=c.sort_values('Effect estimate',ascending=False).reset_index(drop=True)
    b=fig.add_axes([.43,.14,.33,.47])
    b.axvline(0,color=GREY,lw=.8,ls='--')
    b.set_ylim(len(c)-.2,-1.2)
    b.set_yticks([])
    b.set_xlabel('Effect and 95% CI\n(natural-log abundance)',labelpad=6.1875)
    for i,r in c.iterrows():
        color=MAGENTA if r['iPHoP host genus (score90)']=='Parabacteroides' else BLUE
        b.plot([r['95% CI lower'],r['95% CI upper']],[i,i],color=color,lw=1)
        b.scatter(r['Effect estimate'],i,s=22,marker='D' if color==MAGENTA else 'o',facecolor='white' if color==MAGENTA else color,edgecolor=color,lw=.8,zorder=3)
        y=b.transData.transform((0,i))[1]/fig.bbox.height
        fig.text(.035,y,r.vOTU,fontsize=7.5,va='center')
        fig.text(.15,y,str(r['iPHoP host genus (score90)']),fontsize=7.5,va='center',style='italic')
        fig.text(.835,y,str(int(r.host_genes)),fontsize=7.5,ha='center',va='center',color=ORANGE if r.host_genes>r.viral_genes else 'black')
        fig.text(.935,y,str(int(r.viral_genes)),fontsize=7.5,ha='center',va='center')
    fig.text(.035,.657,'b',fontsize=11,fontweight='bold')
    fig.text(.07,.657,'Candidate effects and CheckV gene counts',fontsize=9,fontweight='bold')
    for x,t,ha in ((.035,'vOTU','left'),(.15,'Predicted host genus','left'),(.835,'Host\ngenes','center'),(.935,'Viral\ngenes','center')):
        fig.text(x,.625,t,fontsize=7.5,fontweight='bold',ha=ha,va='bottom')
    return fig
