import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def quality_per_gene(reads,on='quality_mean',gene_name='target',format_base_quality=False):
    if format_base_quality==False:
        reads['quality_all_bases']=reads['quality_all_bases'].str.replace('[','')
        reads['quality_all_bases']=reads['quality_all_bases'].str.replace(']','')
        quality_per_base=pd.DataFrame(list(reads['quality_all_bases'].str.split(','))).astype(float)
        for col in quality_per_base.columns:
            reads['qc_base'+str(col+1)]=list(quality_per_base.loc[:,col])
    score=on
    ordervals=reads.groupby(gene_name).mean()[score].sort_values()
    valsdict=dict(zip(ordervals.index,np.round(ordervals,2)))
    reads['meangenequality']=list(reads[gene_name].map(valsdict))
    reads=reads.sort_values(by='meangenequality')
    plt.figure(figsize=(6,len(valsdict)/1.2))
    sns.violinplot(y=gene_name, x=score, data=reads)
    plt.title([score,'for each gene'])
    return ordervals
def quality_per_cycle(reads,cycles=5,format_base_quality=False):
    bases=cycles
    if format_base_quality==False:
        reads['quality_all_bases']=reads['quality_all_bases'].str.replace('[','')
        reads['quality_all_bases']=reads['quality_all_bases'].str.replace(']','')
        quality_per_base=pd.DataFrame(list(reads['quality_all_bases'].str.split(','))).astype(float)
        for col in quality_per_base.columns:
            reads['qc_cycle'+str(col+1)]=list(quality_per_base.loc[:,col])
    ess=[]
    for e in range(1,bases+1):
        ess.append('qc_cycle'+str(e))
    qualities=reads.loc[:,ess]
    qlt=pd.DataFrame(columns=['quality','base'])
    qualist=[]
    baselist=[]
    for qual in qualities.columns:
        qualist=qualist+list(qualities[qual])
        qa=[qual]*len(qualities)
        baselist=baselist+list(qa)
    output=pd.DataFrame(columns=['quality','cycle'])
    output['quality']=qualist
    output['cycle']=baselist
    plt.figure(figsize=(bases*2,7))
    sns.violinplot(y='quality', x='cycle', data=output)
    plt.title('Qualities for each cycle')
def compare_scores(reads,score1='quality_minimum',score2='quality_mean',kind='kde',color='#3266a8',format_base_quality=False,hue=None): # option for kind are “scatter” | “kde” | “hist” | “hex” | “reg” | “resid” 
    if hue=='assigned':
        reads['assigned']=list(~reads['target'].isna())
    if format_base_quality==False:
        reads['quality_all_bases']=reads['quality_all_bases'].str.replace('[','')
        reads['quality_all_bases']=reads['quality_all_bases'].str.replace(']','')
        quality_per_base=pd.DataFrame(list(reads['quality_all_bases'].str.split(','))).astype(float)
        for col in quality_per_base.columns:
            reads['qc_base'+str(col+1)]=list(quality_per_base.loc[:,col])
    sns.jointplot(x=reads.loc[:,score1], y=reads.loc[:,score2], kind=kind, color=color,hue=reads[hue])

def plot_scores(reads,on='quality_mean',hue='None',log_scale=False,format_base_quality=False,palette='ch:rot=-.25,hue=1,light=.75'):
    if hue=='assigned':
        reads['assigned']=list(~reads['target'].isna())
    if format_base_quality==False:
        reads['quality_all_bases']=reads['quality_all_bases'].str.replace('[','')
        reads['quality_all_bases']=reads['quality_all_bases'].str.replace(']','')
        quality_per_base=pd.DataFrame(list(reads['quality_all_bases'].str.split(','))).astype(float)
        for col in quality_per_base.columns:
            reads['qc_base '+str(col+1)]=list(quality_per_base.loc[:,col])
    sns.histplot(reads,x=on, hue=hue,multiple="stack",palette=palette,edgecolor=".3",linewidth=.5,log_scale=False)
    if hue=='assigned':
        sns.displot(
        data=reads,
        x=on, hue=hue,
        kind="kde", height=6,
        multiple="fill", clip=(0, None),
        palette=palette,
        )
def plot_frequencies(reads,on='targets'):
    readssum=reads.groupby(on).count()
    readssum[on]=list(readssum.index)
    readssum=readssum.sort_values(by='fov')
    plt.figure(figsize=(10,len(readssum)/4))
    plt.title('Number of each '+on)
    ax=sns.barplot(x="fov", y="target", data=readssum)
    ax.set(xlabel='counts', ylabel=on)
    subset=readssum.iloc[:,0:1]
    subset.columns=['counts']
    return subset
def plot_expression(reads,key='target',colorcode="colorblind",xcolumn='xc',ycolumn='yc',genes='all',size=8,background='white',figuresize=(10,7),save=None,format='pdf',title_color='black'): 
    adataobs=reads
    sizecols=len(adataobs[key].unique())
    cls=sns.color_palette(colorcode,sizecols)
    cls2=cls.as_hex()
    colors=dict(zip(adataobs[key].unique(),cls2))
    #cl.apply(lambda x: colors[x])
    plt.rcParams['figure.facecolor'] = background
    if genes=='all':
        cl=adataobs[key]
        plt.figure(figsize=figuresize)
        figa=plt.scatter(x=adataobs[xcolumn],y=adataobs[ycolumn],c=cl.apply(lambda x: colors[x]),s=size,linewidths=0, edgecolors=None)
        plt.axis('off')
        if not save==None:
            plt.savefig(save +'/map_all_genes_'+str(size)+'_'+background+'_'+key+'.'+format)
    elif genes=='individual':
        cl=adataobs[key]
        for each in adataobs[key].unique():
            adatasubobs=adataobs.loc[adataobs.loc[:,key]==each,:]
            plt.figure(figsize=figuresize)
            plt.scatter(x=adataobs[xcolumn],y=adataobs[ycolumn],c='grey',s=size/5,linewidths=0, edgecolors=None)
            cl=adatasubobs[key]
            plt.scatter(x=adatasubobs[xcolumn],y=adatasubobs[ycolumn],c=cl.apply(lambda x: colors[x]),s=size,linewidths=0, edgecolors=None)
            plt.axis('off')
            plt.title(str(each)+':' + str(len(adatasubobs))+' reads',color=title_color)
            if not save==None:
                plt.savefig(save +'/map_individual_cluster_'+str(each)+'_'+str(size)+background+'_'+key+'.'+format)
    else:
        adatasubobs=adataobs.loc[adataobs[key].isin(genes),:]
        plt.figure(figsize=figuresize)
        plt.scatter(x=adataobs.X,y=adataobs.Y,c='grey',s=size/5,linewidths=0, edgecolors=None)
        cl=adatasubobs[key]
        plt.scatter(x=adatasubobs.X,y=adatasubobs.Y,c=cl.apply(lambda x: colors[x]),s=size,linewidths=0, edgecolors=None)
        plt.axis('off')
        plt.legend()
        if not save==None:
                s=''
                for element in genes:
                    s=s+str(element)
                print(s)
                plt.savefig(save +'/map_group_of_clusters_'+str(s)+'_'+str(size)+background+'_'+key+'.'+format)
    #        plt.title('Group: '+ paste(clusters))
    plt.rcParams['figure.facecolor'] = 'white'
def filter_reads(reads,min_quality_mean=False,min_quality_minimum=False,max_distance=False,max_radius=False,min_radius=False,min_intensity=False,max_intensity=False):
    readsfilt=reads
    if not max_distance==False:
        readsfilt=readsfilt.loc[readsfilt['distance']<max_distance,:]
    if not min_quality_mean==False:
        readsfilt=readsfilt.loc[readsfilt['quality_mean']>min_quality_mean,:]
    if not min_quality_minimum==False:
        readsfilt=readsfilt.loc[readsfilt['quality_minimum']>min_quality_minimum,:]
    if not max_radius==False:
        readsfilt=readsfilt.loc[readsfilt['radius']<max_radius,:]
    if not min_radius==False:
        readsfilt=readsfilt.loc[readsfilt['radius']>min_radius,:]
    if not min_intensity==False:
        readsfilt=readsfilt.loc[readsfilt['intensity']>min_intensity,:]
    if not max_intensity==False:
        readsfilt=readsfilt.loc[readsfilt['intensity']<max_intensity,:]
    return readsfilt
