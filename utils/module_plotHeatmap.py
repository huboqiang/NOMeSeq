from __future__ import division
import re,sys,os
import cPickle as pickle
import numpy as np
import pandas as pd
import matplotlib
import subprocess,time
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import brewer2mpl
import scipy.stats
from collections import Counter

import seaborn as sns
sns.set(style="ticks")

cdictRG = {'red':((0.0, 1.0, 1.0),
                  (0.3, 0.8, 0.8),
                  (0.7, 0.5, 0.5),
                  (1.0, 0.0, 0.0)),
         'green':((0.0, 1.0, 1.0),
                  (0.3, 0.9, 0.9),
                  (0.7, 0.9, 0.9),
                  (1.0, 0.8, 0.8)),
         'blue':( (0.0, 1.0, 1.0),
                  (0.3, 0.8, 0.8),
                  (0.7, 0.5, 0.5),
                  (1.0, 0.0, 0.0))} 

cdictRR = {'red':((0.0, 1.0, 1.0),
                  (0.3, 0.9, 0.9),
                  (0.7, 1.0, 1.0),
                  (1.0, 1.0, 1.0)),
         'green':((0.0, 1.0, 1.0),
                  (0.3, 0.8, 0.8),
                  (0.7, 0.5, 0.5),
                  (1.0, 0.0, 0.0)),
         'blue':( (0.0, 1.0, 1.0),
                  (0.3, 0.8, 0.8),
                  (0.7, 0.5, 0.5),
                  (1.0, 0.0, 0.0))}
cdictTest = {'red':( (0.00, 0.07, 0.07),
                  (0.25, 0.09, 0.09),
                  (0.50, 1.00, 1.00),
                  (0.75, 0.98, 0.98),
                  (1.00, 0.91, 0.91)),
         'green':((0.00, 0.05, 0.05),
                  (0.25, 0.36, 0.36),
                  (0.50, 1.00, 1.00),
                  (0.75, 0.71, 0.71),
                  (1.00, 0.37, 0.37)),
         'blue':( (0.00, 0.36, 0.36),
                  (0.25, 0.76, 0.76),
                  (0.50, 1.00, 1.00),
                  (0.75, 0.12, 0.12),
                  (1.00, 0.02, 0.02))}

cdictRKG= {'red':((0.0, 1.0, 0.0),
                  (0.5, 0.0, 0.0),
                  (1.0, 1.0, 0.0)),
         'green':((0.0, 0.0, 1.0),
                  (0.5, 0.0, 0.0),
                  (1.0, 0.0, 1.0)),  
         'blue': ((0.0, 0.0, 0.5),
                  (0.5, 0.0, 0.0),
                  (1.0, 0.0, 0.0))}

cdictRKB= {'red':((0.0, 1.00, 0.00),

                  (1.0, 1.00, 0.00)),
         'green':((0.0, 0.00, 0.50),

                  (1.0, 0.00, 0.00)),
         'blue': ((0.0, 0.00, 1.00),

                  (1.0, 0.00, 1.00))}

my_cmap1 = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictRKG, 256)
my_cmap2 = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictRKB, 256)
my_cmap3 = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdictTest, 256)


def zscore(data):
    mean = np.mean(data,axis=1)
    dev  = np.var(data,axis=1)
    print np.shape(data),np.shape(mean),np.shape(dev)
    z = (data-mean[:,None])/(dev[:,None])**0.5
    mask = np.isnan((data - mean[:,None])/dev[:,None]**0.5)+np.isinf((data - mean[:,None])/dev[:,None]**0.5)
    z[ mask ] = 0.0
    return z, mean, dev


def get_highest_stage(infile, mat_specific="specific.rev.v2.mat"):
    M_posInfo = {}
    pd_sam = pd.read_csv(mat_specific, sep="\t")
    l_samp = pd_sam['tag']
    
    outFile = "%s.highest_corr.pair" % ( ".".join(infile.split(".")[:-1]) )
    outFileSorted = "%s.highest_corr.sort.pair" % ( ".".join(infile.split(".")[:-1]) )
    f_outFile = open(outFile, "w")
    with open(infile, "r") as f_infile:
        for line in f_infile:
            line = line.strip()
            f = line.split()
            pos = f[-1]
            corr = float(f[1])
            if pos not in M_posInfo:
                M_posInfo[pos] = {"corr" : corr, "line" : line}
            else:
                if (corr > 0 and corr > M_posInfo[pos]) or (corr < 0 and corr < M_posInfo[pos]):
                    M_posInfo[pos]['corr'] = corr
                    M_posInfo[pos]['line'] = line
    for pos in M_posInfo:
        print >>f_outFile, M_posInfo[pos]['line']
    
    f_outFile.close()
    pd_infile = pd.read_csv(outFile, header=None, sep="\t", index_col=[4])
    pd_out = pd_infile.sort([0, 1, 3], ascending=[0, 0, 0])
    
    l_idx_new = []
    for samp in l_samp:
        l_idx_new += list(pd_out[ pd_out[0]== samp].index)

    pd_out = pd_out.loc[l_idx_new]
    
    pd_out[4] = pd_out.index
    pd_out.to_csv(outFileSorted, sep="\t", header=None, index=False)


def read_pdInput(infile, l_order, mat_specific="specific.rev.v2.mat", idx_name="chrpos__gene"):
    pd_sam = pd.read_csv(mat_specific, sep="\t", index_col=[0])
    pd_Exp = pd.read_csv(infile, sep="\t")
    l_sam = pd_sam.columns
    pd_Exp.index = [ inf.split("__")[0] for inf in pd_Exp[idx_name] ]
    pd_Exp = pd_Exp.drop(idx_name, 1)
    pd_Exp = pd_Exp.loc[l_order]
    pd_Exp = pd_Exp[l_sam]
    return pd_Exp


def plot_matrix( fig, M_matrix_info ):
    main_mat    = M_matrix_info['main']['matrix']
    main_mat_pos= M_matrix_info['main']['pos']
    main_y_grid = M_matrix_info['main']['grid']
    main_cmap   = M_matrix_info['main']['cmap']
    main_yanno  = M_matrix_info['main']['anno']
    bottom_mat     = M_matrix_info['bottom']['matrix']
    bottom_mat_pos = M_matrix_info['bottom']['pos']
    bottom_x_grid  = M_matrix_info['bottom']['grid']
    bottom_cmap    = M_matrix_info['bottom']['cmap']
    bottom_samp    = M_matrix_info['bottom']['samp']
    color_pos      = M_matrix_info['cmap']['pos']
    color_tick     = M_matrix_info['cmap']['tick']
    color_label    = M_matrix_info['cmap']['label']
    
    if bottom_mat_pos is None:
        bottom_tag_pos_y = main_mat_pos[1] - 0.025
        bottom_mat_pos   = [ main_mat_pos[0],     bottom_tag_pos_y,main_mat_pos[2],0.005 ]
    
    if color_pos is None:
        color_tag_pos_x  = main_mat_pos[0] + main_mat_pos[2]/3 
        upper_space      = 1 - main_mat_pos[3] - main_mat_pos[1]
        color_tag_pos_y  = main_mat_pos[1] + main_mat_pos[3] + upper_space/2
        color_pos        = [ color_tag_pos_x,color_tag_pos_y, 0.08,      0.02  ]

             
    main_matrix = fig.add_axes( main_mat_pos )
    main_matrix.grid(True,color='black',linestyle='-',linewidth=3)
#    main_matrix.grid(True,color='white',linestyle='-',linewidth=5)
    cax01 = main_matrix.imshow(main_mat, aspect='auto', cmap=main_cmap,interpolation='nearest')
    
    main_matrix.get_xaxis().set_ticks([])
    main_matrix.get_yaxis().set_ticks(main_y_grid)
    main_matrix.get_yaxis().set_ticklabels([])
    main_matrix.get_xaxis().set_ticklabels([])
    for sp in main_matrix.spines.values():
        sp.set_visible(False)

    if main_yanno is not None:
        for i in range( 0,len(main_yanno) ):
           main_matrix.text(main_mat.shape[1], i+0.5, '%s' % (main_yanno[i]), ha = "left", fontsize=12)
    
    bottom_tag_matrix = fig.add_axes( bottom_mat_pos )
    bottom_tag_matrix.grid(True,color='white',linestyle='-',linewidth=3)
    cax2 = bottom_tag_matrix.imshow( bottom_mat, aspect='auto', cmap=bottom_cmap,interpolation='nearest')
    bottom_tag_matrix.get_xaxis().set_ticks( bottom_x_grid )
    bottom_tag_matrix.get_yaxis().set_ticks( [] )
    bottom_tag_matrix.get_yaxis().set_ticklabels( [] )
    bottom_tag_matrix.get_xaxis().set_ticklabels( [] )
    bottom_tag_matrix.tick_params(colors="white")
    for i in  range( 0,len(bottom_samp) ) :
        bottom_tag_matrix.text( i+0.2,3,'%s' % ( bottom_samp[i]),rotation=270,ha = "right",size=15)
    for sp in bottom_tag_matrix.spines.values():
        sp.set_visible(False)
    
    cbaxes = fig.add_axes(color_pos)
    cbar = plt.colorbar( cax01,cax=cbaxes, orientation='horizontal',ticks=color_tick )
    cbar.ax.yaxis.set_ticks_position('left')
    cbar.set_label( color_label, size=12)
    cbar.ax.tick_params(labelsize=10)
    
    if 'corr' in M_matrix_info:
        corr_matrix = M_matrix_info['corr']['matrix']
        corr_cmap = M_matrix_info['corr']['cmap']
        corr_pos = M_matrix_info['corr']['pos']
        corr_tick = M_matrix_info['corr']['tick']
        corr_label = M_matrix_info['corr']['label'] 
        
        if corr_tick == None:
            corr_tick = [main_mat.min(), main_mat.max()]
        
        if corr_pos is None:
            color_tag_pos_x  = main_mat_pos[0] + main_mat_pos[2] + 0.01
            corr_pos         = [ color_tag_pos_x,main_mat_pos[1],0.01,main_mat_pos[3] ]
        cor_color_tag_pos_x  = corr_pos[0] + corr_pos[2]/3 
        cor_upper_space      = 1 - corr_pos[3] - corr_pos[1]
        cor_color_tag_pos_y  = corr_pos[1] + corr_pos[3] + cor_upper_space/2
        cor_color_pos        = [ cor_color_tag_pos_x,cor_color_tag_pos_y, 0.1,      0.02  ]
        
        cor_info_matrix = fig.add_axes( corr_pos )
        cor_info_matrix.grid(True,color='white',linestyle='-',linewidth=2)
        my_norm = matplotlib.colors.Normalize(-1, 1)
        cax3 = cor_info_matrix.imshow( np.array([corr_matrix]).T, aspect='auto', cmap=corr_cmap,interpolation='nearest',norm=my_norm)
        cor_info_matrix.get_xaxis().set_ticks([])
        cor_info_matrix.get_yaxis().set_ticks(main_y_grid)
        cor_info_matrix.get_yaxis().set_ticklabels([])
        cor_info_matrix.get_xaxis().set_ticklabels([])
        cor_info_matrix.tick_params(colors="white")
        for sp in cor_info_matrix.spines.values():
            sp.set_visible(False)
       
        cbaxes = fig.add_axes(cor_color_pos)
        cbar = plt.colorbar( cax3,cax=cbaxes, orientation='horizontal',ticks=corr_tick )
        cbar.ax.yaxis.set_ticks_position('left')
        cbar.set_label( corr_label, size=12)
        cbar.ax.tick_params(labelsize=10)



def plot_boxPlot(df_mat, outprefix, mat_specific="specific.rev.mat", show_log10=0):
    pd_sam = pd.read_csv(mat_specific, sep="\t", index_col=[0])

    df_mat_melt = pd.melt(df_mat, id_vars=['tags'])
    l_samp = df_mat.columns
    df_mat_melt = df_mat_melt[ df_mat_melt['value'].notnull() ]
    
    if show_log10 == 1:
        df_mat_melt['value'] = np.log10(df_mat_melt['value'].values + 1)
    
    l_state = []

    for tags, samp in zip(df_mat_melt['tags'].values, df_mat_melt['variable'].values):
        if samp in l_samp[pd_sam.loc[tags].values == 1]:
            l_state.append(1)
        else:
            l_state.append(0)

    df_mat_melt['state'] = l_state
    df_mat_melt['sample'] = [ "%s_%d" % (samp, state) for samp, state in zip(df_mat_melt['variable'].values, df_mat_melt['state'].values) ]
    pal = {sample: "r" if sample.split("_")[-1] == "1" else "b" for sample in df_mat_melt['sample'].unique()}

    g = sns.FacetGrid(df_mat_melt, row="tags", margin_titles=True, aspect=1.5)
    g.map(sns.boxplot, "sample", 'value', palette=pal, showfliers=False)
    #g.map(sns.boxplot, "state", 'value', order=[0, 1], showfliers=False)

    g.set_axis_labels("Stages", "GpC Methylation level in NDR region");
    g.fig.subplots_adjust(wspace=0.2, hspace=0.2)

    for ax in g.axes.flat:
        ax.set_xticklabels(l_samp, rotation=90, ha="right")
        if show_log10 != 1:
            ax.set_ylim(0, 100)
    #    ax.set_xticklabels([0, 1], rotation=90, ha="right")


    out_pdf = "%s.boxplot.pdf" % (outprefix)
    g.savefig(out_pdf, format="pdf", figsize=(15, 30))


def show_help():
    print >>sys.stderr,"\n\tpython",sys.argv[0],"MatrixTSS_Exp.Stage.corr.allPositive.pair MatrixTSS_Exp.mat MatrixTSS_GCH.mat MatrixTSS_WCG.mat specific.rev.v2.mat"


def main():
    try:
        in_pair = sys.argv[1]
        in_mat_exp = sys.argv[2]
        in_mat_GCH = sys.argv[3]
        in_mat_WCG = sys.argv[4]
        in_mat_specific = sys.argv[5]
    except IndexError:
        show_help()
        sys.exit(1)
         
    get_highest_stage(in_pair, in_mat_specific)

    pair_sort = "MatrixTSS_Exp.Stage.corr.allPositive.rev.highest_corr.sort.pair"
    get_mat_plot(in_mat_exp, in_mat_GCH, in_mat_WCG, pair_sort, in_mat_specific)
    get_box_plot(in_mat_exp, in_mat_GCH, in_mat_WCG, pair_sort, in_mat_specific)
    
        
if __name__ == '__main__':
    main()
