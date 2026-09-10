#!/home/rk2643/miniforge3/envs/loter/bin/python

### This script runs Loter, a tool to estimate ancestry hybridization with 
# haplotype variants
# See this tutorial: https://github.com/bcm-uga/Loter/blob/master/python-package/Local_Ancestry_Example.md

import yaml
import allel
import numpy as np
import matplotlib.pyplot as plt



#sys.path.append('./Loter/python-package')
import loter.locanc.local_ancestry as lc


#Load config settings
with open("../../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

genome = config["GENOME"]

code_path_prefix = config["CODE_PATH"]
code_path = f'{code_path_prefix}/popgen/scripts/loter'

original_vcf_path = config["CM_VCF_PATH"]
vcf_file_path = f'{original_vcf_path}'

loter_script_path=f'{code_path}/popgen/Loter'

output_path_prefix = config["DATA_PATH"]
output_path = f'{output_path_prefix}/popgen/loter'


#this makes a numpy matrix from your vcfs
def vcf2npy(vcfpath):
    callset = allel.read_vcf(vcfpath)
    haplotypes_1 = callset['calldata/GT'][:,:,0]
    haplotypes_2 = callset['calldata/GT'][:,:,1]
    
    m, n = haplotypes_1.shape
    mat_haplo = np.empty((2*n, m))
    mat_haplo[::2] = haplotypes_1.T
    mat_haplo[1::2] = haplotypes_2.T
    
    return mat_haplo.astype(np.uint8)





if __name__ == '__main__':
    populations = ['eyed', 'eyeless', 'surface']

    npy_dict ={}
    for population in populations:


        vcf_file = f'{vcf_file_path}/{genome}.phased.{population}.vcf.gz.vcf.gz'
        
        matrix_output = vcf2npy(vcf_file)
        #make the numpy matrix and collect the outputs
        npy_dict[population] = matrix_output
        print(f'for {population}: {matrix_output.shape} ; dictionary:{len(npy_dict)}, {len(npy_dict[population])}')

    #heres where you call on the loter thing
    ancestry_group = [npy_dict['eyeless'],npy_dict['surface']]
    admixed_group = npy_dict['eyed']
    res_loter = lc.loter_smooth(l_H=ancestry_group, h_adm=admixed_group, num_threads=8)


    np.savetxt(f"{output_path}/loter_output_phased.txt", res_loter, fmt="%i")


    ### SEE loter_analysis.py for figure making
    #plt.imshow(res_loter, interpolation='nearest', aspect='auto')
    #plt.colorbar()

    #plt.show()
    #np.savetext(f"{output_path}/loter_output.txt", new_plot, fmt="%i")





#from subprocess import call
#call(["python", "test.py"])

"""
#pulled straight from tutorial
#makes a numpy matrix with 0,1,2
def vcf2npy(vcfpath):

l_H: a list of "ancestral" or reference haplotypes matrices. Its length is equal to the number of ancestral populations.
h_adm: a matrix of admixed haplotypes

import loter.locanc.local_ancestry as lc

res_loter = lc.loter_smooth(l_H=[H_ceu, H_yri], h_adm=H_mex, num_threads=8) ## set the number of threads
"""
