#----------------------------------------                                                                                                                                           
# Outside imports                                                                                                                                                                   
#----------------------------------------                                                                                                                                           
import os
import sys
import pywt
import pandas as pd
import numpy as np
from scipy.stats import entropy
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#----------------------------------------                                                                                                                                           
# Local imports                                                                                                                                                                     
#----------------------------------------                

from region_extraction import *
from wavelet_transforms import *

#load files
def read_annotations(annotations:str):
    
    '''import annotation file
    
    Parameters
    ----------
    annotations : str
        path to annotation in gtf/gff format
            
    Returns
    -------
    gene_coordinates : list of lists
        coordinates contained in gtf/gff files
            
    Raises
    ------
    Exception 
        If file format is not either gff or gtf
    
    '''
    
    gene_coordinates = []

    if annotations.endswith('.gtf'):
        with open(annotations) as annotation_file:
            for line in annotation_file:
                line = line.strip()
                if not line:  
                    continue
                if line.startswith("#"):  
                    continue
                line = line.split('\t')
                gene_id = line[8].split(';')[2].split(' ')[2].strip('"')
                gene_coordinates.append([line[0], line[3], line[4], line[6], gene_id])

    elif annotations.endswith(('gff3', 'gff')):
        with open(annotations) as annotation_file:
            for line in annotation_file:
                line = line.strip()
                if not line:  
                    continue
                if line.startswith("#"):  
                    continue
                line = line.split('\t')
                gene_id = line[8].split(';')[3].split('=')[1]
                gene_coordinates.append([line[0], line[3], line[4], line[6], gene_id])

    elif annotations.endswith('.bed'):
        with open(annotations) as annotation_file:
            for line in annotation_file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    continue
                line = line.split('\t')
                gene_coordinates.append([line[0], line[1], line[2], line[3], line[4]])
    else:
                
        raise Exception('Check annotation file. File format should be .gtf or .gff or .bed with format "chrm, start, stop, strand, id"')
        
    return gene_coordinates

def load_bedgraph(bedgraphs):
    '''Load paths to bedgraph samples to be run                                                                                                                                                                                                                              
    '''
    samples = []
    paths = []
    protocols = []
    lib_preps = []
    
    with open(bedgraphs) as beds:
        for line in beds:
            line = line.strip()
            line = line.split('\t')
            samples.append(line)
            paths.append(line[0])
            protocols.append(line[1])
            lib_preps.append(line[2])
            
    return samples, paths, protocols, lib_preps



def get_pca_summary(dwt_ad:list, ordered_protocol:list,ordered_library:list,pcs=2):

    def color_row(row):
        
        if row['protocol'] == 'GRO':
            val = 0
        elif row['protocol'] == 'PRO':
            val = 1
        else:
            val = 'NA'
        return val

    def color_library(row):
        if row['library'] == 'LIG':
            val = 0
        elif row['library'] == 'CIRC':
            val = 1
        elif row['library'] == 'TSRT':
            val = 2
        elif row['library'] == 'RPR':
            val = 3
        else:
            val = 'NA'
        return val


    #covert approximation or detail coeff lists to np.array
    dwt_array = np.array(dwt_ad)
    
    #perform PCA
    pca = PCA(n_components=pcs) 

    pca.fit(dwt_array)  
    PCA(copy=True, iterated_power='auto', n_components=pcs, random_state=None,
        svd_solver='auto', tol=0.0, whiten=False)
    
    pca_results = pca.fit_transform(dwt_array)
    
    ##convertiong pca results to a pandas dataframe
    pca_df = pd.DataFrame(dict(pc1=pca_results[:, 0], 
                               pc2=pca_results[:, 1],
                               pc3=pca_results[:, 2],
                               pc4=pca_results[:, 3],
                              protocol=ordered_protocol,
                              library=ordered_library))
    pca_df['color'] = pca_df.apply(color_row, axis=1)
    pca_df['lib_color'] = pca_df.apply(color_library, axis=1)
    
    return pca_df, pca


def plot_pca(pca_df,pca,gene,wavefxn,outdir,protocol_color=True):

    #specifying colour variables in hex 
    colors = {'GRO':'#e66101',
            'PRO':'#5e3c99'}

    lib_cols = {'LIG':'#e66101',
                'CIRC':'#fdb863',
                'TSRT':'#b2abd2',
                'RPR':'#5e3c99',}

    lib_cols_set1 = {'LIG':'#E41A1C',
                     'CIRC':'#377EB8',
                     'TSRT':'#4DAF4A',
                     'RPR':'#984EA3'}
    
    if protocol_color == True:
        ##setting up to plot the figures
        fig, ax1 = plt.subplots(figsize=(8,8))

        ax1.scatter(pca_df['pc1'], pca_df['pc2'], s=100, alpha=0.75,
                    c=pca_df['protocol'].apply(lambda x: colors[x]),
                    edgecolors='gray')
        ax1.set_ylabel('PC2 ({}% Variance Explained)'.format(round(pca.explained_variance_ratio_[1]*100,2)),
                    fontsize=18) 
        ax1.set_xlabel('PC1 ({}% Variance Explained)'.format(round(pca.explained_variance_ratio_[0]*100,2)),
                    fontsize=18) 
        ax1.set_title("Detail Coefficients for {} - {}".format(gene, wavefxn),
                    fontsize=20,fontweight='bold')
        ax1.axvline(x=0, color='blue', linestyle='--', alpha=0.5)

        purple_patch = mpatches.Patch(color='#e66101', label='GRO-seq')
        orange_patch = mpatches.Patch(color='#5e3c99', label='PRO-seq')

        plt.legend(handles=[purple_patch, orange_patch],fontsize=14)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.savefig(outdir +'{}_{}_pc1pc2.png'.format(gene,wavefxn),bbox_inches='tight')
        
        plt.close()
    else:
        
        ##setting up to plot the figures
        fig, ax2 = plt.subplots(figsize=(8,8))

        ax2.scatter(pca_df['pc1'], pca_df['pc2'], s=100, alpha=0.75,
                    c=pca_df['library'].apply(lambda x: lib_cols_set1[x]), edgecolors='gray')
        ax2.set_ylabel('PC2 ({}% Variance Explained)'.format(round(pca.explained_variance_ratio_[1]*100,2)),
                    fontsize=18)
        ax2.set_xlabel('PC1 ({}% Variance Explained)'.format(round(pca.explained_variance_ratio_[0]*100,2)),
                    fontsize=18)
        ax2.set_title("Detail Coefficients for {} - {}".format(gene, wavefxn),
                    fontsize=20,fontweight='bold')

        ax2.axvline(x=0, color='blue', linestyle='--', alpha=0.5)
        
        patch1 = mpatches.Patch(color='#E41A1C', label='Ligation')
        patch2 = mpatches.Patch(color='#377EB8', label='Circularization')
        patch3 = mpatches.Patch(color='#4DAF4A', label='TSRT')


        plt.legend(handles=[patch1, patch2, patch3],fontsize=14)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.savefig(outdir +'{}_{}_pc1pc2_libprep.png'.format(gene,wavefxn),bbox_inches='tight')
        
        plt.close()


def run(sample_paths:str, outdir:str, gene_annotations_path:str, wavelet:str, plot=True):
    
    starting_time0 = time.time()
    
    gene_annotations = read_annotations(gene_annotations_path)

    sample_inf, samples, protocol, library = load_bedgraph(sample_paths)

    pca_summary = []
    pc1_results = []
    pc2_results = []
    genes_run = []

    for coords in gene_annotations:
        approx_output = []
        detail_output = []
        coverage = []
        positions_coverage = []
        
        chrom = coords[0]
        start = int(coords[1])
        stop = int(coords[2])
        strand = coords[3]
        gene_name = coords[4]
        
        for bedgraph in samples:
            print("PROCESSING ...... ", bedgraph)
            starting_time = time.time()

            new_name = bedgraph.split('/')[-1].split('.')[0]
            print("Sample ID: ",new_name)
            
            print("Processing transcript coordinates for : {}".format(coords[4]))
            print(coords[0], coords[1], coords[2], coords[3])
                        
            if str(coords[3]) == '+':
                _filtering = bedgraph_awk_strand(bedgraph, chrom, start, stop, outdir, str(new_name),positive=True)
            elif str(coords[3]) == '-':
                _filtering = bedgraph_awk_strand(bedgraph, chrom, start, stop, outdir, str(new_name),positive=False)
            else:
                _filtering = bedgraph_awk(bedgraph, chrom, start, stop, outdir, str(new_name))
                raise Exception('Strand should be either "+"or "-". Check annotation stand column and normalized coverage output.')

            expand_starting_list, expand_coverage_list = expand_bedgraph(_filtering,
                                                                start, stop)
            positions_coverage.append(expand_starting_list)
            coverage.append(expand_coverage_list)

            print('Discrete wavelet transformation....')
            approx, details = dwt_coefficients(_filtering, start, stop, wavelet)
            approx_output.append(approx)
            detail_output.append(details)    
                     
            stoping_time = time.time()
            print("Sample run time: {} seconds \n".format(stoping_time-starting_time))

        print('SAVING DWT FILES .......')
        np.savetxt(outdir+str(coords[4])+'_detail.tsv', np.array(detail_output), delimiter="\t")
        np.savetxt(outdir+str(coords[4])+'_approx.tsv', np.array(approx_output), delimiter="\t")
        np.savetxt(outdir+str(coords[4])+'_coverage.tsv', np.array(coverage), delimiter="\t")
            
        print("SAMPLE RUN END .........\n")
        
        print('Runing PCA on detail coefficients for {}'.format(coords[4]))
        pca_df, pca = get_pca_summary(detail_output, protocol, library, len(detail_output))
        
        pc1_R_gro = pca_df[(pca_df['protocol'] == 'GRO') & (pca_df['pc1'] > 0)].shape[0]
        pc1_L_gro = pca_df[(pca_df['protocol'] == 'GRO') & (pca_df['pc1'] < 0)].shape[0]
        pc1_R_pro = pca_df[(pca_df['protocol'] == 'PRO') & (pca_df['pc1'] > 0)].shape[0]
        pc1_L_pro = pca_df[(pca_df['protocol'] == 'PRO') & (pca_df['pc1'] < 0)].shape[0]
        
        pc1 = round(pca.explained_variance_ratio_[0]*100,2)
        pc2 = round(pca.explained_variance_ratio_[1]*100,2)
        
        pca_summary.append([coords[0], coords[1], coords[2], coords[3], coords[4],
                            pc1, pc2, pc1_R_gro, pc1_L_gro, pc1_R_pro, pc1_L_pro])

        pc1_results.append(list(pca_df['pc1']))
        pc2_results.append(list(pca_df['pc2']))
        genes_run.append(coords[4])
        
        if plot == True:
            print('Plotting PCA plots ')
            plot_pca(pca_df, pca, coords[4], wavelet, outdir)
            plot_pca(pca_df, pca, coords[4], wavelet, outdir, protocol_color=False)
        else:
            pass
        
    
    print('SAVING PC FILES .......')
    np.savetxt(outdir+'PCA_summary_{}.tsv'.format(wavelet),
               np.array(pca_summary), delimiter="\t",fmt='%s')
    
    np.savetxt(outdir+'PC1_results_{}.tsv'.format(wavelet),
               np.array(pc1_results), delimiter="\t", fmt='%s')

    np.savetxt(outdir+'PC2_results_{}.tsv'.format(wavelet),
               np.array(pc2_results), delimiter="\t", fmt='%s')

    np.savetxt(outdir+'gene_names.tsv',
               np.array(genes_run), delimiter="\t", fmt='%s')

    stoping_time0 = time.time()
    print("TOTAL RUN TIME ........ ", stoping_time0-starting_time0, " seconds")
    
    return pca_summary
