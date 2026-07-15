#!/home/rk2643/miniforge3/envs/loter/bin/python

import yaml
import allel
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.image
from PIL import Image
import glob
import pandas as pd
import statistics
from scipy.stats import sem

with open("../../config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

code_path_prefix = config["CODE_PATH"]
code_path = f'{code_path_prefix}/popgen/scripts/loter'

data_path_prefix = config["DATA_PATH"]

input_file=f'{data_path_prefix}/popgen/loter/loter_output_phased.txt'
output_path=f'{data_path_prefix}/popgen/loter'

#sys.path.append('./Loter/python-package')
loter_script_path=f'./Loter'


    
def eyeless_parse_matrix(input_file):

    #loter_np (from loter) has rows of haplotypes (2 per individual)
    #each column is a score for eyeless/surface (0/1) at that SNP 

    loter_np=np.loadtxt(input_file)
    #print(loter_np.shape[1])

    #make a df to store all of your stuff
    by_row_0_positions=pd.DataFrame(0, np.arange(loter_np.shape[1]-1),np.arange(0,33)).astype(object)

    index_key=0
    for row in loter_np:

        #take this explanation with a grain of salt... idk I got it from stack overflow
        #so check your switches over: [0 0 0 1 1 1] in loter_df by position in by_row

        #concatenate your zeros 
        iszero=np.concatenate(([0],np.equal(row,0).view(np.int8),[0]))
        #figure out where you not-zeros are
        absdiff=np.abs(np.diff(iszero))
        #find where the switch is
        ranges=np.where(absdiff==1)[0].reshape(-1,2).tolist()

        #print(ranges[0:5], type(ranges[0:5]))
        by_row_0_positions.loc[:,index_key]=pd.Series(ranges)
        index_key+=1

    by_row_0_positions.dropna(axis=0,how='all',inplace=True)
    by_row_0_positions.to_csv('eyeless_tracts.txt', index=False)
    #print(by_row_0_positions)

    #for each column, calculate eyeless tract length
    #then average across each column
    #print(by_row_0_positions.loc[1,0], by_row_0_positions.loc[1,0][1]-by_row_0_positions.loc[1,0][0])
    
    #probs terrible way to do this but lets try
    eyeless_mn_sd_df=pd.DataFrame(columns=['Mean_eyeless_length','SD_eyeless_length','Eyeless_proportion'])
    mean_list=[]
    sd_list=[]
    eyeless_proportion_list=[]
    for column_name,column_list in by_row_0_positions.items():
        #print(column_name)
        #print(len(column_list))
        column_list=column_list.dropna()
        list_of_eyeless_tract_lengths=[]
        for cell in column_list:
            #print(cell)
            length_value = cell[1]-cell[0]

            list_of_eyeless_tract_lengths.append(length_value)
        mean_list.append(statistics.mean(list_of_eyeless_tract_lengths))
        sd_list.append(sem(list_of_eyeless_tract_lengths).item())
        eyeless_proportion_list.append(sum(list_of_eyeless_tract_lengths)/loter_np.shape[1])
        #print(eyeless_proportion_list)
        #print(f'mean and sd list:{mean_list, sd_list}')

    eyeless_mn_sd_df['Mean_eyeless_length']=mean_list
    eyeless_mn_sd_df['SD_eyeless_length']=sd_list
    eyeless_mn_sd_df['Eyeless_proportion']=eyeless_proportion_list
    #print_this=f'total length of eyeless tracts: {sum(list_of_eyeless_tract_lengths)}\n\
    #      Mean eyeless tract length: {statistics.mean(list_of_eyeless_tract_lengths)}\n\
    #      SE of eyeless tract length:{sem(list_of_eyeless_tract_lengths)}\n\
    #      Proportion of eyeless alleles: {(sum(list_of_eyeless_tract_lengths)/(loter_np.shape[1]*34))*100}'
    #with open('eyeless_output.txt','w') as file:
    #    file.write(print_this)
    #print(list_of_eyeless_tract_lengths)
    #print(eyeless_mn_sd_df)


    return(eyeless_mn_sd_df)



def sf_parse_matrix(input_file):
    loter_np=np.loadtxt(input_file)
    #print(loter_np.shape[1])

    by_row_0_positions=pd.DataFrame(0, np.arange(loter_np.shape[1]-1),np.arange(0,33)).astype(object)
    
    #this calculates the eyeless tract lengths same as above
    index_key=0
    for row in loter_np:
        iszero=np.concatenate(([0],np.equal(row,0).view(np.int8),[0]))
        absdiff=np.abs(np.diff(iszero))
        ranges=np.where(absdiff==1)[0].reshape(-1,2).tolist()


        by_row_0_positions.loc[:,index_key]=pd.Series(ranges)
        index_key+=1
    by_row_0_positions.dropna(axis=0,how='all',inplace=True)
    by_row_0_positions.to_csv('surface_tracts.txt', index=False)
    #print(by_row_0_positions)

    #for each column, calculate surface tract length
    #then average across each column?

    #probs terrible way to do this but lets try
    list_of_surface_tract_lengths=[]
    sf_mean_sd_df=pd.DataFrame(columns=['Mean_SF_length','SD_SF_length'])
    mean_list = []
    sd_list = []
    for column_name,column_content in by_row_0_positions.items():
        #print(column_name)
        new_column=column_content.dropna()

        #subtract start from previous stop (calculate the length of 1's instead of 0's)
        #this is how you subtract the row[x] start from row[x-1] stop
        sf_length=[new_column[x][0]-new_column[x-1][1] for x in range(1,len(new_column))] 
        #print(f'sf_length: {sf_length[0:5], len(sf_length), statistics.mean(sf_length)}')

        list_of_surface_tract_lengths.append(sf_length)

        #^ outputs a list of lists, so flatten it for calculation
        output_list=[x for xs in list_of_surface_tract_lengths for x in xs]
        #print(f'mean:{statistics.mean(output_list)}\n\
        #    sem:{sem(output_list)}')
        
        mean_list.append(statistics.mean(sf_length))
        sd_list.append(sem(sf_length))


    sf_mean_sd_df['Mean_SF_length']=mean_list
    sf_mean_sd_df['SD_SF_length']=sd_list

    #print(f'mean: {statistics.mean(mean_list)}')
    one_big_list=[x for xs in list_of_surface_tract_lengths for x in xs]
    #print(len(list_of_surface_tract_lengths), len(one_big_list))
    print_this=f'total length of surface tracts: {sum(one_big_list)}\n\
          Mean surface tract length: {statistics.mean(one_big_list)}\n\
          SE of surface tract length:{sem(one_big_list)}\n\
          Proportion of surface alleles: {(sum(one_big_list)/(loter_np.shape[1]*34))*100}'
    #with open ('surface_output.txt','w') as file:
    #    file.write(print_this)
    #print(sf_mean_sd_df)
    #print(statistics.mean(sf_mean_sd_df['Mean_SF_length']))
    #print(print_this)

    return(sf_mean_sd_df)





#The number of generations since the onset of admixture (T admix ) 
# can be estimated using the following equation:
# where LM is the average ancestry tract length from the minor
# parent in Morgans and pB is the proportion of the genome derived 
# from the major parent (the probability of recombining) (31–33).
#Tadmix = 1/LM*PB
# PB(eyeless) = 0.88
# bp to Morgans 1.16 cM/Mb (0.0000000116 Morgan/bp) 
# LM = surface average bp * 1.13e-8

def calculate_admix_time():
    print('reading in eyeless and surface tract info')
    eyeless_tracts =eyeless_parse_matrix(input_file)
    #print(eyeless_tracts)

    sf_ave_tracts = sf_parse_matrix(input_file)
    #print(sf_ave_tracts)

    #columns: mean SF length, SF SD, SF mean in morgans, mean cave percentage
    for_calc = pd.DataFrame()
    for_calc[['mean_SF_length','SD_SF_length']]=sf_ave_tracts
    for_calc['mean_SF_Morgans']=for_calc['mean_SF_length']*0.0000000116

    for_calc['proportion_eyeless']=eyeless_tracts['Eyeless_proportion']
    print(for_calc)

    for_calc['Tadmix']=(1/(for_calc['mean_SF_Morgans']*for_calc['proportion_eyeless']))

    ave_Tadmix=statistics.mean(for_calc['Tadmix'])
    sem_Tadmix=sem(for_calc['Tadmix'])
    print(ave_Tadmix,sem_Tadmix)
    SF_tract_ave=statistics.mean(for_calc['mean_SF_length'])
    SF_tract_sem=statistics.mean(for_calc['SD_SF_length'])
    eyeless_mean_proportion=statistics.mean(for_calc['proportion_eyeless'])

    for_calc.to_csv('Admixture_calculation_out.csv',index=False)
    with open('Admixture_calculation_out.csv','a') as file:
        file.write(f"Surface tracts average: {SF_tract_ave}\n\
                    Surface tracts sem: {SF_tract_sem}\n\
                    Eyeless proportion average: {eyeless_mean_proportion}\n\
                    Mean generations since admixture: {ave_Tadmix}\n\
                   SEM:{sem_Tadmix}\n\
                    NOTE:Rows are haplotypes of individuals")

    """
    LM_sf_ave= (7159*1e-8)
    #LM_sf_ave=(381445*1e-8)
    PB_eyeless_proportion=0.8851124802885144
    Tadmix=(1/(LM_sf_ave*PB_eyeless_proportion))
    print(LM_sf_ave,PB_eyeless_proportion,Tadmix, sem(Tadmix))
    """

if __name__=='__main__':
    #print('making a basic plot')
    #basic_plot(input_file)
    
    #print('reformatting np matrix')
    #reformat_np_matrix(input_file)

    #eyeless_parse_matrix(input_file)

    #parse_matrix(input_file)
    #sf_parse_matrix(input_file)
    
    calculate_admix_time()


"""
def basic_plot(input_file):
    new_np=np.loadtxt(input_file)
    print(new_np.shape)


    #and returns a matrix indicating the origin (regarding the ancetry) of each SNP in each haplotype: 
    # each entry (i,j)
    #corresponds to the index of the ancestral population in the list of ancestral 
    # haplotypes matrices from which the SNPs j in haplotype i originates 
    # (0 for the first population in l_H and for the second population in l_H).
    #0=eyeless
    #1=surface

    #Read across j?? 1s = surface, 0s = eyeless
    #i=34 cells = number of individuals
    #j=number of snps
    #read down v 
    # i1 i2 i3 i4 i5 ...
    # 0  0  0  1  1
    # 1
    # 1
    # 1

    #how to analyze this?
    # count down, 0s = 57 ; 1=210 ; 0s = 523 ; 1s = 15 ; 
    #new_np
    unique, counts = np.unique(new_np, return_counts=True)

    new_dict= dict(zip(unique, counts))
    print(new_dict)

    py_dict ={}
    for index in new_dict:
        new_index = index.item()
        new_value=new_dict[index].item()
        py_dict[new_index]=new_value

    print(f'Proportion of eyeless/total: {(py_dict[0.0]/(py_dict[0.0]+py_dict[1.0]))*100}')




    #new_np
    #this is actually good, just too too big to render
    sub_mask=new_np[:,1:200]
    #print(sub_mask, sub_mask.shape)
    ind_axis, snp_axis = np.where(sub_mask.astype(bool))

    plt.imshow(sub_mask)
    plt.scatter(ind_axis[::2], snp_axis[::2])
    plt.show()

    plt.imsave('idkbro.png', sub_mask)


def reformat_np_matrix(input_file):
    loter_np=np.loadtxt(input_file)
    print(loter_np.shape)

    reformatted_np=loter_np.transpose()

    better_np_file = f'{output_path}/loter_transposed.txt'
    np.savetxt(better_np_file, reformatted_np, fmt="%i")
    print(reformatted_np.shape)

    return(better_np_file)

#See `collapse_to_bed.sh` for collapse_ancestry.py conversion of numpy matrix to bed file

#generate bed files with coordinates for ancestry tracts in each hybrid chromosome in bp. 
# We converted tract lengths in bp to Morgans using the median recombination rate of 
# 1.16 cM/Mb (0.0000000116 Morgan/bp) calculated using a previously published genetic 
# map for A. mexicanus (5). Using the median recombination rate to convert bp to Morgan 
# does not account for the fact that local recombination rate varies across the genome, 
# but is used here to provide a rough estimate for the number of generations since the 
# onset of hybridization between cave and surface fish. We calculated the average minor parent 
# (i.e. surface) tract length for each hybrid individual. We used the script lai_global.py 
# (https://github.com/armartin/ancestry_pipeline) to calculate genome-wide ancestry proportions 
# for each hybrid. We input these parameters into the formula described above to 
# calculate an approximation for time since admixture for each of the 17 hybrid genomes.



def calc_average():
    print('Reading in bed file')
    input_path = '/n/projects/rk2643/caballo_moro_genomics/streamlined/data/popgen/loter'
    input_files=glob.glob(f'{input_path}/loter_bed.*_A.bed')
    #print(input_files)
    with open('individuals.txt', 'r') as line:
        input= [s.rstrip() for s in line.readlines()] #[s.rstrip() for s in f.readlines()]
    print(type(input), input)
    
    x=0
    dict_samp_file={}
    for file in input_files:
        dict_samp_file[input[x]]=file
        x+=1
        if x >=17:
            break
    
    with open('tract_length_totals.txt', 'w') as file:
        #print(dict_samp_file)
        total_df=pd.DataFrame()
        for input_file in input_files:
            print(input_file)
            col_names = ['chr','start','stop','pop','ignore','ignore_2']
            first_df = pd.read_csv(input_file, sep='\t', names=col_names)

            next_df = first_df[['start','stop','pop']].astype(str)

            next_df[['chr','start']]=next_df['start'].str.split('.',n=1,expand=True)
            next_df[['other_chr','stop']]=next_df['stop'].str.split('.',n=1,expand=True)

            ready_df= next_df[['chr','start','stop','pop']]
            ready_df['tract_length']=ready_df['stop'].astype(int)-ready_df['start'].astype(int)
            
            total_tract_sum=sum(abs(ready_df['tract_length']))

            average_tract_length=abs(ready_df['tract_length'].mean())

            df_eyeless = ready_df[ready_df['pop']=='Eyeless']
            eyeless_sum=sum(abs(df_eyeless['tract_length']))

            df_surface = ready_df[ready_df['pop']=='Surface']
            sf_sum=sum(abs(df_surface['tract_length']))

            ind_dict={'eyeless_sum':eyeless_sum, 'surface_sum':sf_sum,
                       'average_tract_length': average_tract_length,'percent_surface':sf_sum/total_tract_sum}
            print(ind_dict)
            break

            #print(f'Average tract length in eyeless {input_file}: {abs(df_eyeless["tract_length"]).mean()}')
            #file.write(f'{input_file}\t{abs(df_eyeless["tract_length"]).mean()}\t{abs(df_surface["tract_length"]).mean()}')
    

    #    for tract in bed_a:
    #  tract = tract.strip().split()
    #  if tract[3] in pops: #this excludes stuff not listed in pops
    #    lai_props[pops.index(tract[3])] += (float(tract[5]) - float(tract[4]))
    #out.write(ind + '\t' + '\t'.join(map(str, [round(i/sum(lai_props), 4) for i in lai_props])) + '\n')
"""