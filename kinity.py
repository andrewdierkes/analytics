#!/usr/bin/env python
# coding: utf-8

# In[4]:


def kinetic_parser(fp_list,residual_params,display_df=True):
    '''kinetic assay parser, use kinity & kinetic_parser to parse data for tsvs labeled with different enzyme concentrations in same folder
    Parameters:
    -----------
    fp_list: is a list of posixpath's of the same enzyme concentration
    residual_params: +/- int value for deviations between predicted and actual values
    display_df: default True, show what is happening in parser
    
    Displays data'''
    
    import re
    import os
    import sys
    import pathlib as path

    import numpy as np
    import pandas as pd

    import seaborn as sns
    import matplotlib as mpl

    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    from scipy.signal import find_peaks
    from scipy import stats

    from datetime import date
    #import pingouin as pg
    from itertools import combinations

    
    #add custom modules
    sys.path.append(r'C:\Users\AndrewDierkes\OneDrive - LumiraDx Group Ltd\Desktop\Experiments\modules')
    from matchers import matchers
    from spatial import spatial
    import adstat as ad
    import kinetics
    
    #parse out kinetic data from tsv's
    dataset_dict, wl_index, spatial_idx_list = kinetics.enzyme_kinetic_analysis(fp_list,True)
    
    #dig deep into dictionaries to extract data and generate dataframe for each it for each assay
    summary_list = []
    for assay, it_dict in dataset_dict.items():

         for it, kinetic_dict in it_dict.items():

                for kinetic_pos, spatial_dict in kinetic_dict.items():
                #remove last 3 spatial positions to align the length of all kinetic_pos to 26
                    key_list = list(spatial_dict.keys())[:-3]
                    short_spatial_dict = {key: spatial_dict[key] for key in key_list}

                    df = pd.DataFrame(short_spatial_dict, index=wl_index)
                    df.index.name = assay + f'_kinetic_position_{kinetic_pos}_' + f'_{it}_IT'
                    #display(df)

                #label True if you want graphs, spatial_idx_list[0][:-3] takes the first 26 values of the first & only item in list
                    data_dict = kinetics.spd_plotter(df, wl_index, spatial_idx_list[0][:-3], assay, it, int(kinetic_pos),False)

                    df_summary = pd.DataFrame(data_dict,index=[0])

                    summary_list.append(df_summary)
    
    #combine data from each assay
    df_dataset = pd.concat(summary_list).sort_index().reset_index(drop=True)
    
    
    #change type of objects in dataframe
    df_dataset = df_dataset.astype({'kinetic_position':'float'})
    df_dataset.sort_values(by=['kinetic_position'],inplace=True)
    
    df_it = kinetics.kinetic_it_seperator(df_dataset,1,0)

    #plot data from ALL assays together for seperate IT dataframes
    it_dict = {}
    for key, value in df_it.items():
        #reset indices for iterrows
        value.reset_index(drop=True,inplace=True)
        display(value)
    #drop triplicate values
        filtered = []
        prev_val = None

        for idx, row in value.iterrows():

            current_val = row['max_rfu']

            if current_val is None:
                break
            elif current_val != prev_val:
                    filtered.append(value.iloc[idx])

            prev_val = current_val

        df_filtered = pd.DataFrame(filtered)
        
#ADJUSTED BELOW, uncomment for original - basically wanted more comparisons between assays
    #adjust kinetic read number 
        #RIGHT BELOW IS PROPER METHOD:
        # df_filtered['kinetic_position'] = df_filtered['kinetic_position'].apply(lambda x: x/3)
        
    #UPDATED method to apply kinetic positions of arbitary value (remove once we remove triplicate)
        df_filtered = pd.DataFrame(filtered).reset_index(drop=True)

        unique_assay = df_filtered.iloc[:,0].unique()
        idx_list = []
        for var in unique_assay:
            idx = df_filtered[df_filtered['assay_name'] == var].index.to_list()
            idx_list.append(idx)

        kinetic_list = []
        for var in idx_list:
            _min = min(var)
            _max = max(var)+1
            #get num of each 
            kinetic_num = [var for var in range(len(df_filtered.iloc[_min:_max,0]))]
            kinetic_list.append(kinetic_num)
        #delist
        kinetic_list = [var for small_list in kinetic_list for var in small_list]

        df_filtered['kinetic_position'] = kinetic_list
        
        #shape returns tuple: (length in rows, length in columns)
            #kinetic_pos = [var for var in range(df_filtered.shape[0])]
            #print(kinetic_pos)
            #df_filtered['kinetic_position'] = kinetic_pos
                                                 
        if display_df:
            print('DF FILTERED BELOW')
            display(df_filtered)

    #plot data and add the dictionary of desired timepoints and associated RFU values to it dict    
        assay_dict = kinetics.kinetic_plotter(df_filtered,0,2,3,key)
        it_dict[f'{key}']=assay_dict
        
        
    #MICHAELIS MENTEN Prep
    mm_dataset = []
    rxn_rate_dataset = []

    #iterate thru it_dict assay and extract concentration from name 
    for it, assay in it_dict.items():

        if float(it) == 10000.0:

            for tsv, xy in assay.items():

                finder = re.compile(r'\d+(?=TC)')
                matcher = finder.findall(tsv)
            #convert to single list and float
                try:
                    concentration = matcher[0]
                except IndexError:
                    print('Did not find conc number in TSV name')
                try:
                    floater = float(concentration)
                    mm_dataset.append(floater)
                except ValueError:
                    print('could not convert to float')

            #find rep number
                rep_finder = re.compile(r'_\d_')
                rep_matcher = rep_finder.findall(tsv)

                try:
                    rep = rep_matcher[0][1]
                except IndexError:
                    print('Did not find the rep in TSV name')

            #combine conc and rep # 
                conc_rep = str(concentration) + '_' + rep


            #create dataframes for each assay at specific it
                combo = []
                for var2, data in xy.items():


                    df = pd.DataFrame(data)

                    combo.append(df)

                df_combo = pd.concat(combo,axis=1).reset_index(drop=True)
                df_combo.index.set_names(f'{conc_rep}',inplace=True)
                df_combo.columns.set_names(f'{it}_IT',inplace=True)

                #df_combo['kinetic_position'] = df_combo['kinetic_position'].apply(lambda x: x/3)

                df_combo.insert(0,'grouping',np.repeat(0,len(df_combo.index)))
                if display_df:
                    display(df_combo)

                #combo algo
                df_rxn_rate = kinetics.rxn_combinations(df_combo,0,2,1,5)
                if display_df:
                    display(df_rxn_rate)
                rxn_rate_dataset.append(df_rxn_rate)

        else:
             print('other it data passed up')
                
        
    #merge dataframes on time_pair... only compare similar timing intervals

    #start with first dataframe
    df_rxn_merge = rxn_rate_dataset[0]

    for df in rxn_rate_dataset[1:]:
        #print(df.index.name,len(df.index))
        df_rxn_merge = pd.merge(df_rxn_merge,df,on='time_pair')

    if display_df:
        display(df_rxn_merge)

        
    #find the reaction rates
    df_rxn_slope = pd.DataFrame()

    i_slope = 0
    offset_slope = 4

    for var in df_rxn_merge.columns[4:len(df_rxn_merge.columns):4]:    
        df_rxn_slope[f'{var}'] = df_rxn_merge.iloc[:,4+i_slope]
        i_slope += offset_slope

    df_rxn_slope.insert(0,'time_pair', df_rxn_merge['time_pair'])
    if display_df:
        display(df_rxn_slope)

    #transpose dataframe and add in concentration for plotting
    df_rxn_rate = df_rxn_slope.transpose().set_axis(df_rxn_slope.iloc[:,0],axis=1)

    #drop what was added into column
    df_rxn_rate.drop('time_pair',inplace=True)

    #add concentration
    conc_list = []
    for var in df_rxn_rate.index:
        conc_finder = re.compile(r'(?<=rate_)\d+')
        conc_matcher = conc_finder.findall(var)
        try:
            conc = conc_matcher[0]
            conc_list.append(conc)
        except IndexError:
            print('could not find concentration in index value')

    df_rxn_rate.insert(0,'[s]',conc_list)
    df_rxn_rate['[s]'] = df_rxn_rate['[s]'].apply(lambda x: int(x))


    #df_rxn_rate.insert(0,'[s]',sorted(mm_dataset))
    df_rxn_rate.sort_values(by=['[s]'])
    
    
    #filter data taking only kinetic time point rxn rates which have an ascending pattern
    try:
        df_avg_rxn_rate = kinetics.kinetic_rate_funnel(df_rxn_rate,True)
        print('above takes the average rr for each concentration')
    except ValueError as e:
        print('error',e)
        pass
    
    #find the best fit line based on filtered data
    try:
        kinetics.michaelis_menten_optimizer(df_avg_rxn_rate,residual_params,True)
    except Exception as e:
        print('error',e)
        pass
    
    #find the best fit line based on the filtered data 
    kinetics.michaelis_menten_optimizer(df_rxn_rate,residual_params,True)




# In[ ]:




