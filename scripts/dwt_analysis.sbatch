#!/bin/bash                                                                                                                                                                                 

#SBATCH --job-name=04_dwt_low_cv_genes # Job name                                                                                      

#SBATCH --mail-type=FAIL,END # Mail events (NONE, BEGIN, END, FAIL, ALL)                                                                                                                    

#SBATCH --nodes=1                                                                                                                                                                           
#SBATCH --ntasks=1 # Num of CPU                                                                                                                                                            

#SBATCH --time=24:00:00 # Time limit                                                                                                                                                         
                                                                                                                                     
#SBATCH --partition short                                                                                                                                                                   
#SBATCH --mem=50gb # Memory ##whats the limit in memory required                                                                                                                             
 

module load python/3.6.3

outdir=/path/to/temp
samples=${outdir}/samples_table_20210108.tsv
annotations=${outdir}/filtered_cv0.55_tpm150.bed

echo ${outdir}
echo ${samples}
echo ${annotations}

python3 ~/bin/nascentwt/nascentwt/ -s ${samples} -o ${outdir}/ -a ${annotations} -w sym5

