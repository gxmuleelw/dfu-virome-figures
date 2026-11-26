import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Rectangle
BLUE, ORANGE, GREY, MAGENTA = '#0072B2', '#E69F00', '#697D88', '#CC79A7'
SOURCE_RECORDS = {}
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,
 'axes.titlesize':9,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':0.65,
 'xtick.major.width':0.6,'ytick.major.width':0.6,'xtick.major.size':2.5,
 'ytick.major.size':2.5,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
 'axes.unicode_minus':False,'savefig.facecolor':'white'})
def canvas(height_mm):
    return plt.figure(figsize=(175/25.4,height_mm/25.4))
def panel(ax, letter, title, x=-.12, y=1.08):
    ax.text(x,y,letter,transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom',ha='left')
    ax.set_title(title,loc='left',pad=10,fontweight='bold')
