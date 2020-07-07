############################### This script is to recover the end parameters of coalHMM on simulated dataset and estimate a rough bias for each branch.

import pandas as pd
import pickle
import os.path

### Step 1: recover parameters frome all the 'estimates' directories for each branch

lst = ['NYCBEN_OTOGAR_LEMCAT_HOMSAP', 'LORTAR_NYCBEN_OTOGAR_LEMCAT', 'NYCBEN_NYCPYG_LORTAR_OTOGAR',
    'GALMOH_OTOGAR_NYCBEN_LEMCAT', 'DAUMAD_LEMCAT_NYCBEN_HOMSAP', 'LEMCAT_MICMUR_DAUMAD_HOMSAP',
    'LEMCAT_PROSIM_MICMUR_DAUMAD', 'HYLLAR_MACMUL_SAGMID_LEMCAT', 'HOMSAP_HYLLAR_MACMUL_SAGMID',
    'HOMSAP_PONABE_HYLLAR_MACMUL', 'GORGOR_HOMSAP_PONABE_MACMUL', 'HOMSAP_PANTRO_GORGOR_PONABE',
    'PANPAN_PANTRO_HOMSAP_PONABE', 'PONABE_PONPYG_HOMSAP_MACMUL', 'HYLLAR_NOMLEU_HOMSAP_MACMUL',
    'HYLLAR_SYMSYN_NOMLEU_HOMSAP', 'HOOHOO_SYMSYN_HYLLAR_HOMSAP', 'HYLLAR_HYLPIL_SYMSYN_HOMSAP',
    'COLGUE_MACMUL_HOMSAP_SAGMID', 'CERMON_MACMUL_COLGUE_HOMSAP', 'MACMUL_MANSPH_CERMON_HOMSAP',
    'MACMUL_MACSIL_MANSPH_COLGUE', 'MACNEM_MACSIL_MACMUL_MANSPH', 'MACASS_MACMUL_MACSIL_MANSPH',
    'MANSPH_PAPHAM_MACMUL_COLGUE', 'CERATY_MANSPH_PAPHAM_COLGUE', 'MANLEU_MANSPH_CERATY_COLGUE',
    'PAPHAM_THEGEL_MANSPH_COLGUE', 'LOPATE_PAPHAM_THEGEL_COLGUE', 'PAPANU_PAPHAM_LOPATE_COLGUE',
    'CERMON_CHLSAB_MACMUL_COLGUE', 'CHLSAB_ERYPAT_CERMON_COLGUE', 'CHLAET_CHLSAB_ERYPAT_CERMON',
    'CERMIT_CERMON_CHLSAB_COLGUE', 'COLGUE_TRAPHA_MACMUL_HOMSAP', 'RHIROX_TRAPHA_COLGUE_MACMUL',
    'PYGNEM_RHIROX_TRAPHA_COLGUE', 'RHIROX_RHISTR_PYGNEM_COLGUE', 'COLGUE_PILTEP_TRAPHA_MACMUL',
    'COLANG_COLGUE_PILTEP_MACMUL', 'PITPIT_SAGMID_HOMSAP_LEMCAT', 'ATEFUS_SAGMID_PITPIT_HOMSAP',
    'CEBALB_SAGMID_ATEFUS_HOMSAP', 'CEBALB_SAPAPE_SAGMID_HOMSAP', 'AOTNAN_SAGMID_CEBALB_HOMSAP',
    'CALJAC_SAGMID_AOTNAN_HOMSAP']

dct = {'species':[], 'chr': [], 'start':[], 'stop':[], 'count':[], 'tau1':[], 'tau2':[], 'tau3':[], 'theta1':[], 'theta2':[], 'c2':[], 'rho':[]}
for species in lst:
    print(species)
    for chrom in range(1,6):
        print(chrom, end='\r')
        if chrom == 23:
            chrom_name = 'X'
        else:
            chrom_name = chrom
        if os.path.isfile('{}/chr_{}/tmp/slice_lst.pickle'.format(species, chrom)) == True:
            with open('{}/chr_{}/tmp/slice_lst.pickle'.format(species, chrom), 'rb') as f:
                slice_lst = pickle.load(f)
                for i in range(len(slice_lst)):
                    if os.path.isfile('{}/chr_{}/tmp/outputs/run_{}/estimates'.format(species, chrom, i)) == True:
                        dct['species'].append(species)
                        dct['chr'].append(chrom_name)
                        dct['start'].append(slice_lst[i][0])
                        dct['stop'].append(slice_lst[i][1])
                        dct['count'].append(slice_lst[i][2])
                        with open('{}/chr_{}/tmp/outputs/run_{}/estimates'.format(species, chrom, i), 'r') as new:
                            for tree in new:
                                if 'tau1' in tree:
                                    dct['tau1'].append(float(tree[7:]))
                                if 'tau2' in tree:
                                    dct['tau2'].append(float(tree[7:]))
                                if 'tau3' in tree:
                                    dct['tau3'].append(float(tree[7:]))
                                if 'theta1' in tree:
                                    dct['theta1'].append(float(tree[9:]))
                                if 'theta2' in tree:
                                    dct['theta2'].append(float(tree[9:]))
                                if 'c2' in tree:
                                    dct['c2'].append(float(tree[5:]))
                                if 'rho' in tree:
                                    dct['rho'].append(float(tree[6:]))
                                    break

df = pd.DataFrame.from_dict(dct)
df.to_csv('./simulation_params_per_chromosome.csv', index=False)

### Step 2: avergae parmaters across runs and chromosomes to obtain one estimate per branch

#df = pd.read_csv('simulation_params_total.csv', sep=',', header = 0, encoding = "ISO-8859-1")
dct = {'sp_short':[], 'estimated_mean_tau1':[], 'estimated_median_tau1':[], 'estimated_var_tau1':[], 'estimated_mean_tau2':[], 'estimated_median_tau2':[], 'estimated_var_tau2':[],  'estimated_mean_tau3':[], 'estimated_median_tau3':[], 'estimated_var_tau3':[], 'estimated_mean_theta1':[], 'estimated_median_theta1':[], 'estimated_var_theta1':[], 'estimated_mean_theta2':[], 'estimated_median_theta2':[], 'estimated_var_theta2':[], 'estimated_mean_c2':[], 'estimated_median_c2':[], 'estimated_var_c2':[], 'estimated_mean_rho':[], 'estimated_median_rho':[], 'estimated_var_rho':[]}
for species in lst:
    dct['sp_short'].append(species)
    dct['estimated_mean_tau1'].append(df[(df['species'] == species)]["tau1"].mean())
    dct['estimated_median_tau1'].append(df[(df['species'] == species)]["tau1"].median())
    dct['estimated_var_tau1'].append(df[(df['species'] == species)]["tau1"].var())
    dct['estimated_mean_tau2'].append(df[(df['species'] == species)]["tau2"].mean())
    dct['estimated_median_tau2'].append(df[(df['species'] == species)]["tau2"].median())
    dct['estimated_var_tau2'].append(df[(df['species'] == species)]["tau2"].var())
    dct['estimated_mean_tau3'].append(df[(df['species'] == species)]["tau3"].mean())
    dct['estimated_median_tau3'].append(df[(df['species'] == species)]["tau3"].median())
    dct['estimated_var_tau3'].append(df[(df['species'] == species)]["tau3"].var())
    dct['estimated_mean_theta1'].append(df[(df['species'] == species)]["theta1"].mean())
    dct['estimated_median_theta1'].append(df[(df['species'] == species)]["theta1"].median())
    dct['estimated_var_theta1'].append(df[(df['species'] == species)]["theta1"].var())
    dct['estimated_mean_theta2'].append(df[(df['species'] == species)]["theta2"].mean())    
    dct['estimated_median_theta2'].append(df[(df['species'] == species)]["theta2"].median())
    dct['estimated_var_theta2'].append(df[(df['species'] == species)]["theta2"].var())
    dct['estimated_mean_c2'].append(df[(df['species'] == species)]["c2"].mean())
    dct['estimated_median_c2'].append(df[(df['species'] == species)]["c2"].median())    
    dct['estimated_var_c2'].append(df[(df['species'] == species)]["c2"].var())
    dct['estimated_mean_rho'].append(df[(df['species'] == species)]["rho"].mean())
    dct['estimated_median_rho'].append(df[(df['species'] == species)]["rho"].median())
    dct['estimated_var_rho'].append(df[(df['species'] == species)]["rho"].var())
    
df_estimated = pd.DataFrame.from_dict(dct)
#dfb.to_csv('./simulation_params_summary.csv', index=False)    

### Step 3: group the tables with initial (=simulated=true) parameters vs recovered (=estimated) parameters

df_true = pd.read_csv('parameters.tsv', sep='\t', header = 0, encoding = "ISO-8859-1")
#df_estimated = pd.read_csv('simulation_params_summary.csv', sep=',', header = 0, encoding = "ISO-8859-1")

df=pd.merge(df_true, df_estimated, how='inner', on=['sp_short'])

# estimate the bias:

df['bias_tau1']= (df['estimated_mean_tau1']-df['tau1_mean'])/df['tau1_mean']
df['bias_tau2']= (df['estimated_mean_tau2']-df['tau2_mean'])/df['tau2_mean']
df['bias_tau3']= (df['estimated_mean_tau3']-df['tau3_mean'])/df['tau3_mean']
df['bias_theta1']= (df['estimated_mean_theta1']-df['theta1_mean'])/df['theta1_mean']
df['bias_theta2']= (df['estimated_mean_theta2']-df['theta2_mean'])/df['theta2_mean']
df['bias_c2']= (df['estimated_mean_c2']-df['c2_mean'])/df['c2_mean']
df['bias_rho']= (df['estimated_mean_rho']-df['rho_mean'])/df['rho_mean']

df['bias_median_tau1']= (df['estimated_median_tau1']-df['tau1_median'])/df['tau1_median']
df['bias_median_tau2']= (df['estimated_median_tau2']-df['tau2_median'])/df['tau2_median']
df['bias_median_tau3']= (df['estimated_median_tau3']-df['tau3_median'])/df['tau3_median']
df['bias_median_theta1']= (df['estimated_median_theta1']-df['theta1_median'])/df['theta1_median']
df['bias_median_theta2']= (df['estimated_median_theta2']-df['theta2_median'])/df['theta2_median']
df['bias_median_c2']= (df['estimated_median_c2']-df['c2_median'])/df['c2_median']
df['bias_median_rho']= (df['estimated_median_rho']-df['rho_median'])/df['rho_median']

### Step 4: Obtain real true parameters in good units: 
# T1, T2 are the intervals between speciation times, T3 here is third interval + theta2 (not exactly speciation time,divergence time.
df['true_T1']=(df['estimated_mean_tau1']/(df['bias_tau1']+1))/df['mutation_rate']
df['true_T2']=(df['estimated_mean_tau2']/(df['bias_tau2']+1))/df['mutation_rate']
df['true_T3']=(df['estimated_mean_tau3']/(df['bias_tau3']+1))/df['mutation_rate']
df['true_Ne1']=(df['estimated_mean_theta1']/(df['bias_theta1']+1))/(df['mutation_rate']*2*df['generation_time'])
df['true_Ne2']=(df['estimated_mean_theta2']/(df['bias_theta2']+1))/(df['mutation_rate']*2*df['generation_time'])
df['true_r']=(df['estimated_mean_rho']/(df['bias_rho']+1))*df['mutation_rate']*df['generation_time']*1e8


### Step 5: Obtain estimated parameters in good units: 

df['estimated_T1']=df['estimated_mean_tau1']/df['mutation_rate']
df['estimated_T2']=df['estimated_mean_tau2']/df['mutation_rate']
df['estimated_T3']=df['estimated_mean_tau3']/df['mutation_rate']
df['estimated_Ne1']=df['estimated_mean_theta1']/(df['mutation_rate']*2*df['generation_time'])
df['estimated_Ne2']=df['estimated_mean_theta2']/(df['mutation_rate']*2*df['generation_time'])
df['estimated_r']=df['estimated_mean_rho']*df['mutation_rate']*df['generation_time']*1e8

df['TRUE_TRUE_tau1']=df['tau1_mean']/(df['bias_tau1']+1)
df['TRUE_TRUE_tau2']=df['tau2_mean']/(df['bias_tau2']+1)
df['TRUE_TRUE_tau3']=df['tau3_mean']/(df['bias_tau3']+1)
df['TRUE_TRUE_theta1']=df['theta1_mean']/(df['bias_theta1']+1)
df['TRUE_TRUE_theta2']=df['theta2_mean']/(df['bias_theta2']+1)
df['TRUE_TRUE_r']=df['rho_mean']/(df['bias_rho']+1)

df.to_csv('./rough_bias_estimation.csv', index=False, sep="\t")  








