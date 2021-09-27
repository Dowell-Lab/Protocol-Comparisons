import os
import sys
import time
import pywt
import numpy as np
import matplotlib.pyplot as plt


from region_extraction import bedgraph_awk,bedgraph_awk_strand,coverage_norm,extract_regions,expand_bedgraph

####################################################################
##perform discrete wavelet transformations on genomic coordinates ##
####################################################################
def dwt_coefficients(annotations:list,
                     start_coor:int,stop_coor:int, 
                     dwtwavelet:str): 
    '''give genomic coordinates get the detail and approximation coefficients

    Parameters
    ----------
    annotations : list
        list of coordinates [chr, start, stop, coverage]

    start_coor : int
        start coordinate

    stop_coor : int
        stop coordinate

    dwtwavelet : str
        discrete wavelet to used (e.g. sym5, db4, coif5, haar)

    Returns
    -------
    aD : list
        Approximation coefficients from DWT

    cD : list
        Detail coeffiecient from DWT
    '''

    ##first expand genomic coordinates
    expand_starting_list, expand_coverage_list = expand_bedgraph(annotations,
                                                                start_coor, stop_coor)
    wavelet = pywt.Wavelet(dwtwavelet)
    no_moments = str(wavelet.vanishing_moments_psi)
    
    ##returns approximation and detail coefficients.
    aD, cD = pywt.dwt(expand_coverage_list, dwtwavelet)
    
    return aD, cD

##############################################################
###Plotting the raw coverage and different transformations ###
##############################################################
def dwt_coverage_figure(annotations:list, chrm:str, 
                        start_coor:int, stop_coor:int,
                        gene_id:str, dwtwavelet:str,
                       ntransforms=2):
    '''Plot and visualize regions being surveyed

    Parameters
    ----------
    annotations : list
        region coordinates [chr, start, stop, coverage]

    chrm : str
        chromosomes ID of filtered region
    
    start_coor : int
        start coordinate

    stop_coor : int
        stop coordinate

    gene_id : str
        name of gene 

    dwtwavelet : str
        wavelet to use for transformation

    ntransforms : int
        number of transformations (default = 2)
    
    Returns
    -------
    plot of normalized coverage, detail and approximation coefficients
    
    '''

    ##expand region to analyze
    expand_starting_list, expand_coverage_list = expand_bedgraph(annotations,
                                                                 start_coor,stop_coor)
    
    wavelet = pywt.Wavelet(dwtwavelet)

    print("--------       Wavelet Transformation       --------")
    print("-----         Using the {} wavelet         -----".format(dwtwavelet))
    print("----- The {} wavelet has {} vanishing moments -----".format(dwtwavelet, 
                                                                                 str(wavelet.vanishing_moments_psi)))
    ##plot region called
    fig, ax = plt.subplots(figsize=(14,2))
    ax.set_title(''.join(["Original coverage: ",
                 gene_id, " ", str(chrm),":",
                 str(start_coor),"-" ,str(stop_coor)]), fontweight='bold')
    ax.plot(expand_starting_list, expand_coverage_list, color='purple', alpha=0.95)
    ax.axhline(y=0, color='blue', linestyle='--', alpha=0.5)
    plt.show()

    data = expand_coverage_list
    waveletname = dwtwavelet

    ##plot coefficients
    fig, axarr = plt.subplots(nrows=ntransforms, ncols=2, figsize=(12,3))
    for ii in range(ntransforms):
        (data, coeff_d) = pywt.dwt(data, waveletname)
        axarr[ii, 0].plot(data, 'r')
        axarr[ii, 0].axhline(y=0, color='blue', linestyle='--', alpha=0.5)
        axarr[ii, 1].plot(coeff_d, 'g')
        axarr[ii, 0].set_ylabel("Level {}".format(ii + 1), 
                                fontsize=14, rotation=90, fontweight='bold')
        axarr[ii, 0].set_yticklabels([])
        if ii == 0:
            axarr[ii, 0].set_title("Approximation coefficients", 
                                   fontsize=14,fontweight='bold')
            axarr[ii, 1].set_title("Detail coefficients", 
                                   fontsize=14,fontweight='bold')
        axarr[ii, 1].set_yticklabels([])
    plt.tight_layout()
    plt.show()


##########################################################
##run the DWT on multiple bedgraphs for a single region###
##########################################################

def multiple_run(samples:list, chrom:str, 
                 start:int, stop:int, outdir:str, 
                 gene_id:str, dwtwavelet:str, positive_strand=True):
    
    '''run the discrete wavelet transformation on multiple samples
    
    Parameters
    ----------
    samples : list
        list with paths to bedgraphs

    chrm : str
        chromosome for region to be extracted

    start : int
        start coordinate

    stop : int
        stop coordinate
    
    outdir : str
        path to save output

    gene_id : str
        gene/transcript name

    dwtwavelet : str
        discrete wavelet to use

    positive : bool
        True if forward strand is being extracted
        False if reverse stranf is being extracted


    Returns
    -------
    approx_output : list of lists
        approximation coefficients for all processed samples

    detail_output : list of lists
        detail coefficients for all processed samples

    '''

    starting_time0 = time.time()

    approx_output = []
    detail_output = []

    for sample in samples:

        print("PROCESSING ...... ", sample)
        starting_time = time.time()

        new_name = sample.split('/')[-1].split('.')[0]
        print("Sample ID: ",new_name)

        if positive_strand == True:
            _filtering = bedgraph_awk_strand(sample, chrom, start, stop, outdir, str(new_name),positive=True)
        else:
            _filtering = bedgraph_awk_strand(sample, chrom, start, stop, outdir, str(new_name),positive=False)
            
        approx, details = dwt_coefficients(_filtering, start, stop, dwtwavelet)

        approx_output.append(approx)
        detail_output.append(details)
        
        dwt_coverage_figure(_filtering, chrom, start, stop, gene_id, dwtwavelet)

        stoping_time = time.time()
        print("Sample run time: ",stoping_time-starting_time, " seconds")
        print("                                                         ") 

    stoping_time0 = time.time()

    print("OVERALL RUN TIME ........ ", stoping_time0-starting_time0, " seconds")
    
    return approx_output, detail_output
