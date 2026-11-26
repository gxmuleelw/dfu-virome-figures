draw <- function(genes, backbone, labels, lysis, lysis_label, notes, titles, palette) {
    library(ggplot2)
p <- ggplot() +
  geom_segment(data=backbone,aes(x=0,xend=length_kb,y=0,yend=0),colour='#737B85',linewidth=.3) +
  gggenes::geom_gene_arrow(data=genes,show.legend=TRUE,
    aes(xmin=left/1000,xmax=right/1000,y=track,forward=forward,fill=plot_group),
    arrowhead_width=grid::unit(1.2,'mm'),arrowhead_height=grid::unit(3.4,'mm'),
    arrow_body_height=grid::unit(2.5,'mm'),colour='white',linewidth=.15) +
  geom_segment(data=labels,aes(x=midpoint,xend=midpoint,y=track,yend=label_y-.12),
    colour='#8A929C',linewidth=.25) +
  geom_segment(data=labels,aes(x=midpoint,xend=label_x,y=label_y-.12,yend=label_y-.04),
    colour='#8A929C',linewidth=.25) +
  geom_text(data=labels,aes(x=label_x,y=label_y,label=label),size=2.65,lineheight=.95,
    colour='#17202A',family='Arial') +
  geom_segment(data=lysis,aes(x=left,xend=right,y=.55,yend=.55),colour='#D55E00',linewidth=.45) +
  geom_segment(data=lysis,aes(x=left,xend=left,y=.48,yend=.55),colour='#D55E00',linewidth=.45) +
  geom_segment(data=lysis,aes(x=right,xend=right,y=.48,yend=.55),colour='#D55E00',linewidth=.45) +
  geom_text(data=lysis,aes(x=(left+right)/2,y=.80),label=lysis_label,
    size=2.65,lineheight=.95,colour='#17202A',family='Arial') +
  geom_text(data=notes,aes(x=0,y=1.47,label=note),hjust=0,size=2.7,colour='#525B65',family='Arial') +
  facet_wrap(~contig,ncol=1,labeller=as_labeller(titles),axes='all_x',drop=FALSE) +
  scale_fill_manual(values=palette,drop=FALSE,name=NULL) +
  scale_x_continuous(breaks=seq(0,35,5),expand=expansion(mult=0)) +
  coord_cartesian(xlim=c(-.2,37.0),ylim=c(-.42,1.57),clip='on') +
  labs(x='Position on viral representative (kb)',y=NULL) +
  guides(fill=guide_legend(nrow=3,byrow=TRUE)) +
  theme_classic(base_size=9,base_family='Arial') +
  theme(axis.line.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
    axis.line.x=element_line(linewidth=.25,colour='#737B85'),axis.ticks.x=element_line(linewidth=.25),
    axis.text.x=element_text(size=8,colour='#35404A'),axis.title.x=element_text(size=9,margin=margin(t=5)),
    strip.background=element_blank(),strip.text=element_text(size=10.5,face='bold',hjust=0,margin=margin(b=4)),
    panel.spacing.y=grid::unit(6,'mm'),legend.position='bottom',legend.text=element_text(size=8),
    legend.key.size=grid::unit(3.5,'mm'),legend.margin=margin(t=0,b=0),legend.box.margin=margin(t=0,b=0),
    plot.caption=element_text(size=7.1,hjust=0,lineheight=1.15,colour='#525B65',margin=margin(t=6)),
    plot.margin=margin(8,10,8,10))
return(p)
}
