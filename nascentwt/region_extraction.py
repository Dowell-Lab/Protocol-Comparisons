'''This module extracts regions of interest in a bedgraph file and then proceeds 
to expand regions for each base position.
'''

#==============================================================================
__author__ = 'Rutendo F. Sigauke'
__credits__ = ['Rutendo F. Sigauke', 'Jacob T. Stanley',
                'Robin D. Dowell']
__maintainer__ = 'Rutendo F. Sigauke'
__email__ = ['Rutendo.Sigauke@colorado.edu', 'Rutendo.Sigauke@CUAnschutz.edu']
#==============================================================================

import os
import sys
import time
import numpy as np


#############################
##get regions of interest  ##
#############################

def bedgraph_awk(bedgraph, chrm, 
                start, stop,
                outdir, tempbedgraph):

    '''filter bedgraph coordinates (ie. by chromosome, start, and sttop positions.)                                                                     
    Parameters                                                                                                                                          
    ----------                                                                                                                                          
    bedgraph : str                                                                                                                                      
        path to bedgraph                                                                                                                                

    chrm : st                                                                                                                                           
        chromosome to filter                                                                                                                            
                                                                                                                                                        
    start : int                                                                                                                                         
        start coordinate                                                                                                                                
                                                                                                                                                        
    stop : int                                                                                                                                          
        stop coordinate                                                                                                                                 

    outdir : str                                                                                                                           
        path to outout directory   

    tempbedgraph : str                                                                                                                                  
        name of temporary bedgraph file                                                                                                                 
                                       
    Returns                            
    -------                            
    final_bedgraph : list              
        filtered bedgraph coordinates in list format                                                                                                    
    '''

    
    os.system("cat " + bedgraph + " | grep -w " + chrm + " | awk '{ if ($2 >= " + str(start) + " && $3 <= " + str(stop) + " ) print $0 }' > " + outdir + tempbedgraph + ".TEMP.bedgraph")


    infile = str(outdir) + str(tempbedgraph) + ".TEMP.bedgraph"
    
    final_bedgraph = []
    
    with open(infile) as bed:
        for row in bed:
            row in row.strip("\n")
            chrom, start, stop, coverage = row.split("\t")
            final_bedgraph.append([chrom, start, stop, coverage])
    
    os.system("rm " + infile)
    
    return final_bedgraph


##########################################
####Save region counts based on strand ###
##########################################
def bedgraph_awk_strand(bedgraph, chrm, 
                       start, stop,
                       outdir, tempbedgraph,
                      positive=True):
    '''Extranct given a bedgraph region and save as list of coordinates 

    Parameters
    ----------
    bedgraph : str
        path to bedgraph

    chrm : str
        chromosome id

    start : int
        start coordinates

    stop : int
        stop coordinates

    outdir : str
        path to output temporary files
    
    tempbedgraph : str
        name for temporary file
    
    positive : bool
        True if positive strand

    Returns
    -------
    final_bedgraph : list 
        list of bedgraph regions with read counts 
        for each region

    '''

 
    def filter_strand(infile):
        final_bedgraph = []
        with open(infile) as bed:
            for row in bed:
                row in row.strip("\n")
                chrom, start, stop, coverage = row.split("\t")
                final_bedgraph.append([chrom, start, stop, coverage])

        os.system("rm " + infile)
        return final_bedgraph
    
    if positive is True:
        
        os.system("cat " + bedgraph + " | grep -w " + chrm + " | awk '{ if ($2 >= " + str(start-10) + " && $3 <= " + str(stop+10) + " ) print $0 }' | awk '{ if ($4 >= 0) print $0 }' > " + outdir + tempbedgraph + ".TEMP.bedgraph")
        infile_positive = str(outdir) + str(tempbedgraph) + ".TEMP.bedgraph"
        
        return(filter_strand(infile_positive))
          
    else:
        os.system("cat " + bedgraph + " | grep -w " + chrm + " | awk '{ if ($2 >= " + str(start-10) + " && $3 <= " + str(stop+10) + " ) print $0 }' | awk '{ if ($4 <= 0) print $0 }' > " + outdir + tempbedgraph + ".TEMP.bedgraph")
        infile_negative = str(outdir) + str(tempbedgraph) + ".TEMP.bedgraph"
        
        return(filter_strand(infile_negative)) 


########################
####Normalize counts ###
########################
def coverage_norm(cov):
    '''Normalize counts for each region

    Parameters
    ----------
    cov : list
       list with coverage information from coordinates

    Returns
    -------
    cov_new : list
        max normalized coverage 
    '''
    
    max_cov = max(cov)
    min_cov = min(cov)

    cov_new = []
   
    for i in range(len(cov)):
        #try:                                                                                                                                                                                                                                  
        #_cov_new = (cov[i] - min_cov)/(max_cov-min_cov)                                                                                                                                                                                       
        #cov_new.append(_cov_new)   
        try:
            _cov_new = (cov[i] - min_cov)/(max_cov-min_cov)
            cov_new.append(_cov_new)
        except ZeroDivisionError:
            _cov_new = 0
            cov_new.append(_cov_new)
            print('Warning! Coverage may be low. Check bedGraphs.')          

    return cov_new

#################################################################
###make sure coordinates represent each based for called region##
#################################################################
##generate lists with all start and stop regions, since some bedgraphs
##have summarized coordinates and coverage (i.e. will not contain 0 coverage regions)
def extract_regions(annotations:list):
    '''first get all coordates and their respective coverage

    Parameters 
    ----------
    annotations : list
        Cooridnate lists with [chrm, start, stop, coverage]

    Returns
    -------
    sorted_start : list
        sorted start coordinates

    sorted_stop : list
        stop regions sorted by start coordinates
    
    sorted_coverage : list
        coverage sorted by start coordinates
    
    '''

    start_all = []
    stop_all = []
    cover_all = []

    ##storing bedgraph regions to lists
    for rows in annotations:
        start_all.append(rows[1])
        stop_all.append(rows[2])
        cover_all.append(rows[3])
    
    ##below the coverage and stop are sorted by the start coordinate
    sorted_coverage = [x for _,x in sorted(zip(start_all,cover_all))]
    sorted_stop = [x for _,x in sorted(zip(start_all, stop_all))]
    sorted_start = sorted(start_all)
    
    return sorted_start, sorted_stop, sorted_coverage

def expand_bedgraph(annotations:list, 
                    start_coor:int, 
                    stop_coor:int):
    
    '''get all regions including those without coverage

    Parameters
    ----------
    annotations : list
        region coordinates [chrm, start, stop, coverage]

    start_coor : int
        start coordinate

    stop_coor : int
        stop coordinate

    Returns
    -------
    expand_starting_final : list
        all region start coordinates

    expand_coverage_norm : list
        all region normalized coverage includes 0s
    '''

    ##get sorted start, stop and coverage information
    sorted_start, sorted_stop, sorted_coverage = extract_regions(annotations)
    
    ##expanding bedgraph coverage
    expand_coverage = []
    expand_starting = []
    expand_starting_final = []
    expand_coverage_final = []
    
    ##getting expanded coverage 
    ##since bedgraph file format only gives ranges
    for st, sp, cv in zip(sorted_start, sorted_stop, sorted_coverage):
        expand_starting.append(np.arange(int(st), int(sp)+1,1))
        expand_coverage.append([int(cv)] * len(np.arange(int(st), int(sp)+1,1)))
    
    expand_coverage_list = [j for i in expand_coverage for j in i]
    expand_starting_list = [j for i in expand_starting for j in i]
    
    
    ##create a dictionary with all the base positions as keys
    ##and the values be the coverage
    cov_dict = {key: [0] for key in range(int(start_coor), int(stop_coor)+1, 1)}

    for starts, covers in zip(expand_starting_list,expand_coverage_list):
        if starts in cov_dict:
            cov_dict[starts].append(covers)
        else:
            pass

    for k, v in cov_dict.items():
        expand_starting_final.append(k)
        expand_coverage_final.append(v[-1])
    
    expand_coverage_norm = coverage_norm(expand_coverage_final)
    
    return expand_starting_final, expand_coverage_norm
