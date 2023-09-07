#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

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
import pingouin as pg
from itertools import combinations




# In[1]:

def kinetic_analysis(_path,file_type,print_enabled=True):
    '''takes files and parses out kinetic information
    Parameters:
    path = object = path.Path.cwd()
    file_type = string literal '.tsv' '.txt'
    print_enabled = show internal data
    
    Returns dictionary in this order assay:it:kinetic:spatial:pixel'''
    
    
    sys.path.append('/home/jovyan/share/users/Andy D/2023/custom_modules/')
    from matchers import matchers
    from spatial import spatial
    import adstat as ad
    import kinetics
    
    dataset_dict = {}

    it_list = []


    for _filepath in _path.iterdir():
        if _filepath.suffix != r'.tsv':
            continue
        elif _filepath.suffix == r'.tsv':

            with open(_filepath, 'r') as file:

                data = file.read()

            #find specific names of assay
                head, tail = os.path.split(_filepath)
            #extract values
                assay = matchers.name(data)
                spdg = matchers.spdg(data)
                sscparam, wl_start, wl_end = matchers.sscparam(data,3)
                spatial_step_number = matchers.numberofsteps(data)
                spatial_step_size = matchers.stepsize(data)
                num_pixel = matchers.numberpixel(data)
                spatial_start = matchers.occonfig(data)
                it_number = matchers.it_number(data)
                it = matchers.it(data)
                if print_enabled:
                    print(it)
                    print(tail)
                    print('Number of ITs: ',len(it))

                pixel_dataset = {}

            #dictionary for each assay with all ITs and data
                it_dict = {}   
                for _it in it:
                    if _it < 10000:
                        pass
                    else:
                        if print_enabled:
                            print(_it)
                            print('IT=', _it)
                        it_matcher = re.compile(f'dataIndex=[0-9]+, integrationTimeUs={_it}.*[0-9]+\]')


                    #find all data indexes which have it_matcher, some lists will be blank (IE those which have different IT)
                        it_data = []
                        for data_index in spdg:

                            it_ =it_matcher.findall(data_index)

                            it_data.append(it_)


                        it_pixel_data = []

                    #process pixel data information which have same IT
                        for pixel_data in it_data:

                            if len(pixel_data) == 0:
                                continue

                            else:
                            #find pixel data in the strings that are associated with first IT, transform list of one string into a list of float values for each spatial pos 
                                pixel_data_finder = re.compile(r'(?<=pixelData: \[)[-?\d+\.?\d+,\s]+')

                            #find the pixel data in the list of a single string (var[0]),
                                pixel_match = pixel_data_finder.findall(pixel_data[0])


                                try:
                                    px_match = pixel_match[0]

                                except IndexError:
                                    print(f'no data in pixel_match for {_it} in {file}')


                            #split this string
                                pixel_match_list = px_match.split(', ')

                            #convert strings to floats
                                px_match_data = []
                                for var in pixel_match_list:
                                    try:
                                        px_match_data.append(float(var))

                                    except ValueError:
                                        print(f'Couldnt convert:{var} for {_it} in {file}')
                                        continue

                                it_pixel_data.append(px_match_data)

                     #create wavelength profiles
                        wl_index = spatial.wl_index(wl_start,wl_end,num_pixel)


                     #create spatial profiles for each channel
                        spatial_idx_list = []
                        for var in spatial_start:
                            spatial_idx = spatial.spatial_index(var,spatial_step_size, spatial_step_number)
                            spatial_idx_list.append(spatial_idx)


                     #unpack the list of pixel values and seperate into list = to number of steps (number of spatial positions or kinetic reads)
                        pixel_chunk = [*spatial.data_divider(it_pixel_data,spatial_step_number)]

                     #create list for associated kinetic reads
                        kinetic_number = [var for var in range(len(pixel_chunk))]



                    #for every kinetic read, there are n=number of steps (spatial scans) with pixel data within them
                        kinetic_dict = {}

                        i_kinetic = 0
                        offset_kinetic = 1

                        for kinetic_read in pixel_chunk:

                            i_spatial = 0
                            offset_spatial = 1

                            spatial_dict = {}
                            for spatial_position in kinetic_read:

                                #ensure proper number of pixels read
                                if len(spatial_position) != num_pixel:
                                    continue

                                else:

                                #create a dictionary with spatial positions and corresponding pixeldata (0 because we have only 1 channel)
                                    spatial_dict[f'{spatial_idx_list[0][i_spatial]}'] = spatial_position
                                    i_spatial += offset_spatial

                            kinetic_dict[f'{kinetic_number[i_kinetic]}'] = spatial_dict
                            i_kinetic += offset_kinetic



                    #add to dataset for assay
                    it_dict[f'{_it}'] = kinetic_dict

                #add to dataset for whole study
                dataset_dict[f'{tail}'] = it_dict    
        
    return dataset_dict, wl_index, spatial_idx_list

def enzyme_kinetic_analysis(fp_list,print_enabled=True):
    '''this function is meant to be used when you have a dictionary with keys(enzyme conc) and values(list with posixpath), it takes a filepath and parses the kinetic information within the file
    
    Parameters:
    ------------
    _filepath = PosixPath, ideally sorted into lists which hold the same enzyme concentration
    print_enabled = prints out data on the _filepath
    
    Returns dictionary in this order assay:it:kinetic:spatial:pixel'''

    sys.path.append('/home/jovyan/share/users/Andy D/2023/custom_modules/')
    from matchers import matchers
    from spatial import spatial
    import adstat as ad
    import kinetics
    
    dataset_dict = {}

    it_list = []
    for _filepath in fp_list:
        
        with open(_filepath, 'r') as file:

            data = file.read()

        #find specific names of assay
            head, tail = os.path.split(_filepath)
        #extract values
            assay = matchers.name(data)
            spdg = matchers.spdg(data)
            sscparam, wl_start, wl_end = matchers.sscparam(data,3)
            spatial_step_number = matchers.numberofsteps(data)
            spatial_step_size = matchers.stepsize(data)
            num_pixel = matchers.numberpixel(data)
            spatial_start = matchers.occonfig(data)
            it_number = matchers.it_number(data)
            it = matchers.it(data)
            if print_enabled:
                print(it)
                print(tail)
                print('Number of ITs: ',len(it))

            pixel_dataset = {}

        #dictionary for each assay with all ITs and data
            it_dict = {}   
            for _it in it:
                if _it < 10000: #only looking at 100it
                    continue
                else:
                    if print_enabled:
                        print(_it)
                        print('IT=', _it)
                    it_matcher = re.compile(f'dataIndex=[0-9]+, integrationTimeUs={_it}.*[0-9]+\]')


                #find all data indexes which have it_matcher, some lists will be blank (IE those which have different IT)
                    it_data = []
                    for data_index in spdg:

                        it_ =it_matcher.findall(data_index)

                        it_data.append(it_)


                    it_pixel_data = []

                #process pixel data information which have same IT
                    for pixel_data in it_data:

                        if len(pixel_data) == 0:
                            continue

                        else:
                        #find pixel data in the strings that are associated with first IT, transform list of one string into a list of float values for each spatial pos 
                            pixel_data_finder = re.compile(r'(?<=pixelData: \[)[-?\d+\.?\d+,\s]+')

                        #find the pixel data in the list of a single string (var[0]),
                            pixel_match = pixel_data_finder.findall(pixel_data[0])


                            try:
                                px_match = pixel_match[0]

                            except IndexError:
                                if print_enabled: #NEW MAY NEED TO REMOVE
                                    print(f'no data in pixel_match for {_it} in {file}')


                        #split this string
                            pixel_match_list = px_match.split(', ')

                        #convert strings to floats
                            px_match_data = []
                            for var in pixel_match_list:
                                try:
                                    px_match_data.append(float(var))

                                except ValueError:
                                    if print_enabled: #NEW MAY NEED TO REMOVE
                                        print(f'Couldnt convert:{var} for {_it} in {file}')
                                    continue

                            it_pixel_data.append(px_match_data)

                 #create wavelength profiles
                    wl_index = spatial.wl_index(wl_start,wl_end,num_pixel)


                 #create spatial profiles for each channel
                    spatial_idx_list = []
                    for var in spatial_start:
                        spatial_idx = spatial.spatial_index(var,spatial_step_size, spatial_step_number)
                        spatial_idx_list.append(spatial_idx)


                 #unpack the list of pixel values and seperate into list = to number of steps (number of spatial positions or kinetic reads)
                    pixel_chunk = [*spatial.data_divider(it_pixel_data,spatial_step_number)]

                 #create list for associated kinetic reads
                    kinetic_number = [var for var in range(len(pixel_chunk))]



                #for every kinetic read, there are n=number of steps (spatial scans) with pixel data within them
                    kinetic_dict = {}

                    i_kinetic = 0
                    offset_kinetic = 1

                    for kinetic_read in pixel_chunk:

                        i_spatial = 0
                        offset_spatial = 1

                        spatial_dict = {}
                        for spatial_position in kinetic_read:

                            #ensure proper number of pixels read
                            if len(spatial_position) != num_pixel:
                                continue

                            else:

                            #create a dictionary with spatial positions and corresponding pixeldata (0 because we have only 1 channel)
                                spatial_dict[f'{spatial_idx_list[0][i_spatial]}'] = spatial_position
                                i_spatial += offset_spatial

                        kinetic_dict[f'{kinetic_number[i_kinetic]}'] = spatial_dict
                        i_kinetic += offset_kinetic



                #add to dataset for assay
                it_dict[f'{_it}'] = kinetic_dict

            #add to dataset for whole study
            dataset_dict[f'{tail}'] = it_dict    

    return dataset_dict, wl_index, spatial_idx_list


def spd_plotter(df, wl_index, spatial_idx_list, graph_channel,it,kinetic_pos,show):
    '''find max RFU output, peak emission wavelength. From these we will graph and create a dictionary
    ------------------------
    Parameters
    df = spatial dataframe for channel
    wl_index = list of wavelengths scanned in assay
    spatial_idx_list = list of spatial positions scanned in assay 
    wl_index and spatial_idx_list should have same length as df
    graph_channel = list of assay names
    it = integration time
    kinetic_pos = kinetic timepoint
    show= = True or False, if you want to display plots (True)
    all lists are lists of lists and require an offset and chunk method'''

    ch_max = [df.max(axis=1)]
    ch_spatial = spatial_idx_list
    
    
    #iterate thru max list & find max value (max_rfu)
    for var in ch_max:
        max_rfu = (max(var))
        #print(max_rfu)
        
    #find max rfu value and index to find the spatial position in the column names
    try:
        search_spatial = df.isin([max_rfu]).any()
        max_spatial_position = search_spatial[search_spatial].index.values[0]
    except IndexError as e:
        print('Issue finding spatial position in column name',e)
        pass
    

    lmax_index = []

    #iterate thru df columns looking for the lamda max value
    for col in df.columns:
        index_lookup = df[df[col] == max_rfu].index.tolist()
        try:
        #if a variable exists... appended
            lmax_index.append(index_lookup[0])
        except:
            pass

    #wavelength at max rfu
    max_emission_wl = float(lmax_index[0])
    #print(max_emission_wl)

    max_emission_wl_row = wl_index.index(max_emission_wl)


    #avg function
    def average(data, length):
        avg = sum(data)/len(data)
        return avg

    #create an average list similar to plateau if it is within 10% of the RFU
    RFU_max_list = []
    for z in df.iloc[max_emission_wl_row]:
        if z >= max_rfu*.92:
            RFU_max_list.append(z)

    num_spatial_pos = len(RFU_max_list)

    average = average(RFU_max_list, num_spatial_pos)


    #take the spatial data across the max wavelength (df.iloc) and assign variables as arrays for find_peaks
    x = np.array(spatial_idx_list)
    y = np.array(df.iloc[max_emission_wl_row])

                #         #scipy peak finder; use 1d np.array and height = min value for peak 
                #         peaks = find_peaks(y, height = 2000, distance=10)
                #         peak_pos = x[peaks[0]]
                #         height = peaks[1]['peak_heights']


    
    if show is True:
        fig, axs = plt.subplots(1,2)

        fig.suptitle(f'Assay {graph_channel} with an it of {it}; Spatial RFU & Max Emission WL at {round(max_emission_wl, 3)}')


        #plot spatial vs RFU at max wl & maxima via scipy
        axs[0].scatter(x, y)
            #axs[0].scatter(max_spatial_position, max_rfu, color = 'r', s = 30, marker = 'D', label = 'maxima')
        axs[0].set(xlabel='Spatial Position', ylabel='RFU')
            #axs[0].legend(bbox_to_anchor=(1.05,1.05))
        axs[0].grid()

        #plot wavelength vs RFU
        for var in range(len(df.columns)):
            spatial = df.iloc[:,var]

            axs[1].scatter(df.index,spatial,label=df.columns[var])
            axs[1].legend(bbox_to_anchor=(1.05,1.05))
            axs[1].set(xlabel='wavelength',ylabel='RFU')


        fig.tight_layout()
        plt.show()
    else:
        pass


    


    data_dict = {}
    
    
    data_dict['assay_name'] = graph_channel
    data_dict['it'] = it
    data_dict['kinetic_position'] = kinetic_pos
    data_dict['max_rfu'] = max_rfu
    data_dict['average_8'] = average 
    data_dict['number_spatial_positions_for_average_8e'] = num_spatial_pos 
    data_dict['max_emission_wl'] = max_emission_wl
    data_dict['max_emission_spatial'] = max_spatial_position

    return data_dict

def kinetic_plotter(df,groupby,iv,feature,it,df_show=False):
    '''plot data for different groupby(assays) categories, using the feature (dv) and iv from dataframe
    
    Parameters
    ----------
    df = dataframe
    name = col # with assay name
    groupby = col# of category seperating data
    feature = col # for dv
    iv = col # for independent var
    it = integration time
    df_show = bool to show or hide the dataframe(default False)'''
    
    #create variables for plotting
    regex = df.iloc[:,groupby].unique()
    colors = cm.rainbow(np.linspace(0,1,len(regex)))
    
    #extract string values
    iv_name = df.columns[iv]
    groupby_name = df.columns[groupby]
    
    #sort data based on the groupby to align proper min/max 
    df.sort_values(by=[groupby_name,iv_name],inplace=True)
    df.reset_index(drop=True,inplace=True)
    
    #display df
    if df_show:
        display(df)
    
    #find all index values for unique groupby
    idx_list = []
    for var in regex:
        idx = df[df.iloc[:,groupby] == var].index.to_list()
        idx_list.append(idx)
                      
    
    dataset_dict = {}
    
    #zip data for plotting
    plotter_data = zip(regex,colors,idx_list)

    fig, ax = plt.subplots()
    for r,c,i in plotter_data:
        
        assay_dict = {}
        
        #find min and max value for each regex and plot
        _min = min(i)
        _max = max(i)+1
        
        x = df.iloc[_min:_max,iv]
        y = df.iloc[_min:_max,feature]
        
        assay_dict['x'] = x
        assay_dict['y'] = y
        
        dataset_dict[f'{r}'] = assay_dict
        
        print(f'Data for IT of {r} has a min of {_min} and max of {_max-1}')
        
        
        ax.scatter(x,y,marker='.', color=c,label=r)
     
    ax.set(xlabel='kinetic_timepoint',ylabel='rfu',title=f'Kinetics with an IT of {it}')    
    ax.grid()
    ax.legend(bbox_to_anchor=(1.05,1.05))   
    
    plt.show()
        
    return dataset_dict 
    
    
def kinetic_it_seperator(df,groupby,assay):
    '''seperates dataset from 0-x assays into unique dataframes based on IT
    Parameters
    ----------
    df = df_dataset
    groupby = col # with IT values 
    assay = col # with assay name
    
    Returns a dictionary with key(it):value(dataframe with unique it)'''
    
    regex = df.iloc[:,groupby].unique()
    
    groupby_name = df.columns[groupby]
    assay_name = df.columns[assay]
    
    df.sort_values(by= [groupby_name,assay_name],inplace=True)
    df.reset_index(drop=True,inplace=True)
    
    #display(df)
    idx_list = []
    
    for var in regex:
        print(var)
        idx= df[df.iloc[:,groupby] == var].index.to_list()
        idx_list.append(idx)
    
    it_dict = {}
    
    i = 0
    
    for var in idx_list:
        _min = min(var)
        _max = max(var)+1
        print(_min,_max)
        
        df_it = pd.DataFrame(df.iloc[_min:_max])
        display(df_it.head())
        it_dict[f'{regex[i]}'] = df_it
        
        i += 1
    
    return it_dict    

def rxn_combinations(df, groupby, feature, orderby, start):
    '''Function is meant to find the different combinations possible among specific groups (usually IV; groupby). However, to align the
    partitions(groups) we use an order_by before the combination is made, this is important for paired t tests & ANOVAs (as we'd like
    the comparisons the align). We want to do combinations for different groups on the same index

    Parameters
    ----------
    df = dataframe, 
    groupby=col# of the IV, if no IV/all the same, add a column with the same value across all rows
    feature=col# of DV. 
    orderby = col# of indexer for each groupby category, also the time component (kinetic)
    start = index value to start the cominbation algorithim at in the dataframe

    Returns
    -------
    This function returns 1 item; a dataframe of all combinations.'''

    from itertools import combinations
    from statistics import mean

    #get index name before resetting
    df_index = df.index.name

    regex = df.iloc[:, groupby].unique()
    groupby_name = df.columns[groupby]
    orderby_name = df.columns[orderby]

    df.sort_values(by=[groupby_name, orderby_name], inplace=True)
    df.reset_index(drop=True, inplace=True)

    print('check these df insights to ensure the frame is sorted by groupby & the orderby which is repeated across the features')
    #display(df)
    display(df.head())
    display(df.tail())

    #search for values with the same groupby (iv)
    idx_list = []
    for var in regex:
        print(var)
        idx = df[df.iloc[:, groupby] == var].index.to_list()
        idx_list.append(idx)

    combo_rfu = []
    combo_time = []

    #find min and max value from the feature col (RFU usually) using the min/max groupby/regex 
    for var in idx_list:
        _min = start
        _max = max(var)
        print('min', _min, 'max',_max + 1, type(_min)), 'length of data for group:', len(df.iloc[_min:(_max+1),feature])
        print('data taken for group',df.iloc[_min:(_max+1),feature])

        #run combinations on all datapoints
        combo_r = combinations(df.iloc[_min:(_max+1), feature], 2)
        combo_t = combinations(df.iloc[_min:(_max+1), orderby], 2)

        #add to lists
        combo_rfu.append(combo_r)
        combo_time.append(combo_t)


    rfu_diff_list = []
    rfu_pair_list = []
    time_diff_list = []
    time_pair_list = []

    #find difference between two combination values in tuples
    for var in combo_rfu:
        for pair in var:
            diff = pair[1]-pair[0]
            rfu_diff_list.append(diff)
            rfu_pair_list.append(pair)
    #ditto
    for var in combo_time:
        for pair in var:
            diff = pair[1]-pair[0]
            time_diff_list.append(diff)
            time_pair_list.append(pair)

    #create dataframe with combinations and associated substrate concentration (this is to differentiate between different reps of the same substrate concentration)
    #time pair does not have substrate conc bc we will be merging on this column
    df_combo = pd.DataFrame()
    df_combo.index.set_names(df_index,inplace=True)
    df_combo.insert(0,f'rfu_pair_{df_index}',rfu_pair_list)
    df_combo.insert(1,f'rfu_difference_{df_index}',rfu_diff_list)
    df_combo.insert(2,f'time_pair',time_pair_list)
    df_combo.insert(3,f'time_difference_{df_index}',time_diff_list)
    #calculate reaction rate
    df_combo.insert(4,f'reaction_rate_{df_index}', df_combo[f'rfu_difference_{df_index}']/df_combo[f'time_difference_{df_index}'])

    #display(df_combo)

    return df_combo

#VERSION 1
#             def kinetic_rate_funnel(df,show_internal_data=True):
#                 '''iterates thru all columns which are sorted first by [s], then checks to see if the reaction rate ascends with the substrate conc
#                 Parameters
#                 ----------
#                 df = dataframe (df_rxn_rate)
#                 show_internal_data = default True (show plots and dataframes)

#                 Returns timepoint datasets which have ascending average reaction rates'''

#                 #sort by [s] 
#                 df.sort_values(by=['[s]'],inplace=True)

#                 #search for columns with ascending pattern order
#                 asc_reaction_rate = []
#                 list_df_mean_reaction_rate = []

#                 #create a num for each row
#                 col = range(len(df.columns[1:]))

#                 #label the substrate concentrations
#                 sub_unique = df.iloc[:,0].unique()
#                 sub = df.iloc[:,0]

#                 def asc_pattern(col):
#                     '''search col for pattern of low to hi (asc)
#                     all function checks if the condition holds true for all pairs of consecutive elements'''
#                     return all(col.iloc[i] < col.iloc[i+1] for i in range(len(col)-1))

#                 def sub_mean(df, groupby, feature, print_enabled=True):

#                     '''find mean value for different substrate concentrations''

#                     Parameters
#                     ----------
#                     df = dataframe, 
#                     groupby = col # for where to iloc and search for the regex,
#                     feature = col # which you want to aggregate mean, stdev & cv
#                     print_enabled = True default, prints out operations

#                     Returns
#                     -------
#                     dataframe with mean'''

#                     import pandas as pd

#                     regex = df.iloc[:, groupby].unique()

#                     groupby_name = df.columns[groupby]

#                     # sort dataframe to align with regex iterator below(min,max... basically need groupings in one space)
#                     df.sort_values(by=[groupby_name], inplace=True)
#                     df.reset_index(drop=True, inplace=True)

#                     if print_enabled:
#                         display(df.head(), df.tail())
#                         print(f'grouped by {groupby_name}')

#                     # find indexes of regex patterns
#                     idx_list = []
#                     for var in regex:
#                         idx = df[df.iloc[:, groupby] == var].index.to_list()
#                         idx_list.append(idx)

#                     # iterate thru df using idx_list
#                     mean_list = []

#                     for var in idx_list:
#                         # print(var)
#                         if len(var) == 0:
#                             pass
#                         else:
#                             _min = min(var)
#                             _max = max(var)
#                             if print_enabled:
#                                 print('min',_min,'max', _max, type(_min), 'length of data for group:',len(df.iloc[_min:_max+1, feature]),'all data included:',df.iloc[_min:_max+1, feature])
#                             chunk = len(var)
#                             offset = 0

#                             mean = df.iloc[_min:(_max+1), feature].mean()

#                             mean_list.append(mean)

#                     df_final = pd.DataFrame(mean_list, index=[regex], columns=['mean'])

#                     if print_enabled:
#                         display(df_final)

#                     return df_final



#                 #iterate thru all rows in df (omitting [s])
#                 for var in col[1:]:

#                 #create object holding timepoints for each column
#                     rr_timepoint = str(df.columns[var])

#                 #find mean of all reps for each [s]
#                     df_stat = sub_mean(df,0,var,False)

#                     df_stat.columns.set_names(rr_timepoint,inplace=True)
#                     mean_rr = pd.DataFrame(df_stat.iloc[:,0])

#                 #find col with asc value in reaction rates
#                     asc_col = mean_rr.columns[mean_rr.apply(asc_pattern)]

#                 #find asc values amongst the different reaction rates if > 0; it is here
#                     if asc_col.size > 0:

#                     #append col index values to search in df if is True (ASC ORDER)
#                         asc_reaction_rate.append(var)

#                         mean_reaction_rate = df_stat.iloc[:,0]

#                         if show_internal_data:
#                             print('DF_MEAN_RXN_RATE')
#                             display(mean_reaction_rate)

#                         df_mean = pd.DataFrame({f'{rr_timepoint}':mean_reaction_rate})
#                         df_mean.index.set_names('mean_reaction_rate',inplace=True)

#                         if show_internal_data:
#                             print('BELOW DF_MEAN')
#                             display(df_mean)

#                         list_df_mean_reaction_rate.append(df_mean)

#                         if show_internal_data:

#                             fig, ax = plt.subplots()
#                             ax.scatter(sub_unique,df_mean.iloc[:,0])
#                             ax.grid()
#                             ax.set(xlabel='[substrate]',ylabel='reaction_rate',title=f'Slope of reaction rate between timepoints: {rr_timepoint}')
#                             plt.show()

#                 try:
#                     df_avg_rxn_rate = pd.concat(list_df_mean_reaction_rate,axis=1).reset_index(drop=True)
#                     df_avg_rxn_rate.insert(0,'[s]',sub_unique)

#                     return df_avg_rxn_rate

#                 except Exception as e:
#                     print(f'Error: {e}, there were no objects found in the list_df_mean_reaction_rate')
#                     pass

#VERSION 2: 
    #added the raw subplot with the average
def kinetic_rate_funnel(df,show_internal_data=False):
    '''iterates thru all columns which are sorted first by [s], then checks to see if the reaction rate ascends with the substrate conc
    Parameters
    ----------
    df = dataframe (df_rxn_rate)
    show_internal_data = default True (show plots and dataframes)
    
    Returns timepoint datasets which have ascending average reaction rates'''
    
    #sort by [s] 
    df.sort_values(by=['[s]'],inplace=True)
    print('operating on dataframe below:')
    display(df)
    #search for columns with ascending pattern order
    asc_reaction_rate = []
    list_df_mean_reaction_rate = []

    #create a num for each row
    col = range(len(df.columns[1:]))

    #label the substrate concentrations
    sub_unique = df.iloc[:,0].unique()
    sub = df.iloc[:,0]

    def asc_pattern(col):
        '''search col for pattern of low to hi (asc)
        all function checks if the condition holds true for all pairs of consecutive elements'''
        return all(col.iloc[i] < col.iloc[i+1] for i in range(len(col)-1))
    
    def sub_mean(df, groupby, feature, print_enabled=True):
       
        '''find mean value for different substrate concentrations''
        
        Parameters
        ----------
        df = dataframe from df_rxn_rate, 
        groupby = col # for where to iloc and search for the regex,
        feature = col # which you want to aggregate mean, stdev & cv
        print_enabled = True default, prints out operations
        
        Returns
        -------
        dataframe with mean'''
        
        import pandas as pd

        regex = df.iloc[:, groupby].unique()

        groupby_name = df.columns[groupby]

        # sort dataframe to align with regex iterator below(min,max... basically need groupings in one space)
        df.sort_values(by=[groupby_name], inplace=True)
        df.reset_index(drop=True, inplace=True)

        if print_enabled:
            display(df.head(), df.tail())
            print(f'grouped by {groupby_name}')
            
        # find indexes of regex patterns
        idx_list = []
        for var in regex:
            idx = df[df.iloc[:, groupby] == var].index.to_list()
            idx_list.append(idx)

        # iterate thru df using idx_list
        mean_list = []

        for var in idx_list:
            # print(var)
            if len(var) == 0:
                pass
            else:
                _min = min(var)
                _max = max(var)
                if print_enabled:
                    print('min',_min,'max', _max, type(_min), 'length of data for group:',len(df.iloc[_min:_max+1, feature]),'all data included:',df.iloc[_min:_max+1, feature])
                chunk = len(var)
                offset = 0

                mean = df.iloc[_min:(_max+1), feature].mean()

                mean_list.append(mean)

        df_final = pd.DataFrame(mean_list, index=[regex], columns=['mean'])
        
        if print_enabled:
            display(df_final)

        return df_final
    
    
    
    #iterate thru all rows in df (omitting [s])
    for var in col[1:]:

    #create object holding timepoints for each column
        rr_timepoint = str(df.columns[var])

    #find mean of all reps for each [s]
        df_stat = sub_mean(df,0,var,False)
        
        df_stat.columns.set_names(rr_timepoint,inplace=True)
        mean_rr = pd.DataFrame(df_stat.iloc[:,0])

    #find col with asc value in reaction rates
        asc_col = mean_rr.columns[mean_rr.apply(asc_pattern)]

    #find asc values amongst the different reaction rates if > 0; it is here
        if asc_col.size > 0:

        #append col index values to search in df if is True (ASC ORDER)
            asc_reaction_rate.append(var)
            mean_reaction_rate = df_stat.iloc[:,0]
            #dataframe for mean reaction rates of there is ascending order
            df_mean = pd.DataFrame({f'{rr_timepoint}':mean_reaction_rate})
            df_mean.index.set_names('mean_reaction_rate',inplace=True)
            
            list_df_mean_reaction_rate.append(df_mean)
            
            
            #find raw data associated with an ascending order of averages
            for var in df.columns:
                if str(var) == rr_timepoint:
                    location = df.columns.get_loc(var)
            
            #generate df for all reps if df_mean for specific timepoint has asc order
            df_reps = pd.DataFrame({df.columns[0]:df.iloc[:,0],df.columns[location]:df.iloc[:,location]})
            
            #show df_mean and df_reps
            if show_internal_data:
                print('MEAN')
                display(df_mean)
                print('ALL REPS')
                display(df_reps)
                    

#             if show_internal_data:
    
#                 fig, ax = plt.subplots()
#                 ax.scatter(sub_unique,df_mean.iloc[:,0])
#                 ax.grid()
#                 ax.set(xlabel='[substrate]',ylabel='reaction_rate',title=f'Slope of reaction rate between timepoints: {rr_timepoint}')
#                 plt.show()
                
        
            if show_internal_data:
                
                fig, axs = plt.subplots(1,2)
                fig.suptitle(f'Slope of reaction rate between timepoints: {rr_timepoint}')
                fig.subplots_adjust(wspace=0.3)
                #plot averages
                axs[0].scatter(sub_unique,df_mean.iloc[:,0])
                axs[0].grid()
                axs[0].set(xlabel='[substrate]',ylabel='average_reaction_rate')
                
                #plot raw data using location of kinetic timepoint from ascending averages
                axs[1].scatter(df['[s]'],df.iloc[:,location])
                axs[1].grid()
                axs[1].set(xlabel='[substrate]',ylabel='reaction_rate')
                plt.show()
            else:
                pass
     
    try:
        df_avg_rxn_rate = pd.concat(list_df_mean_reaction_rate,axis=1).reset_index(drop=True)
        df_avg_rxn_rate.insert(0,'[s]',sub_unique)
        
        return df_avg_rxn_rate
    
    except Exception as e:
        print(f'Error: {e}, there were no objects found in the list_df_mean_reaction_rate')
        pass
    
def michaelis_menten_optimizer(df,residual_limiter_num,residual_limiter=True):
    '''takes rxn_rate dataframe and finds most optimally fitted line based on data from different kinetic timepoints and performs residual analysis
    Parameters:
    -----------
    df = dataframe (df_avg_rxn_rate or df_rxn_rate)
    residual_limiter_num = min/max residual num
    residual_limiter = default True (limits graphs outputted based on residual_limiter_num)
    
    Returns:
    --------
    equation
    '''
    from scipy.optimize import curve_fit

    #define equation
    def michaelis_menten(substrate, Vmax, Km):
        return Vmax * substrate / (Km + substrate)


    avg_substrate_concentration = np.array(df.iloc[:,0])

    #fit data to michaelis menten that has ASC PATTERN
    for var in range(len(df.columns)):

        if int(var) == 0:
            continue
        else:

        #find col# and take reaction rates
            reaction_rate = np.array(df.iloc[:,var])


            #fit to michealis menten
            params, covar = curve_fit(michaelis_menten, avg_substrate_concentration, reaction_rate)
            fitted_Vmax, fitted_Km = params

            #generate equation
            eq = f'v = ({round(fitted_Vmax,1)} * [s]) / ({round(fitted_Km,1)} + [s])'

            #residual analysis
            predicted_rate = michaelis_menten(avg_substrate_concentration, fitted_Vmax, fitted_Km)


            residuals = reaction_rate - predicted_rate
            
            if residual_limiter:
                if np.min(residuals) > -int(residual_limiter_num) and np.max(residuals) < int(residual_limiter_num):

                #plot model and residual analysis (observed-predicted)
                    fig, axs = plt.subplots(1,2)
                    fig.subplots_adjust(wspace=0.5)
                    fig.suptitle(f'Reaction Rates {eq} between timepoint values: {df.columns[var]} fitted using Michealis Menten & Residual Analysis')
                    
                    #model plot
                    axs[0].scatter(avg_substrate_concentration, reaction_rate)
                    axs[0].plot(avg_substrate_concentration,michaelis_menten(avg_substrate_concentration, fitted_Vmax, fitted_Km), color='r',label='fitted curve')
                    axs[0].grid()
                    axs[0].set(xlabel= '[substrate]',ylabel='reaction_rate (rfu/time)') 

                    #residual plot
                    axs[1].scatter(avg_substrate_concentration, residuals)
                    axs[1].grid()
                    axs[1].set(xlabel='[substrate]',ylabel='residuals')
                    axs[1].yaxis.set_ticks(np.arange(min(residuals)-10,max(residuals)+10,10))
                    
                    plt.show()
            else:
                                                                 
                    fig, axs = plt.subplots(1,2)
                    fig.subplots_adjust(wspace=0.5)
                    fig.suptitle(f'Reaction Rates {eq} between timepoint values: {df.columns[var]} fitted using Michealis Menten & Residual Analysis')
                    #model plot

                    axs[0].scatter(avg_substrate_concentration, reaction_rate)
                    axs[0].plot(avg_substrate_concentration,michaelis_menten(avg_substrate_concentration, fitted_Vmax, fitted_Km), color='r',label='fitted curve')
                    axs[0].grid()
                    axs[0].set(xlabel= '[substrate]',ylabel='reaction_rate (rfu/time)') 

                    #residual plot
                    axs[1].scatter(avg_substrate_concentration, residuals)
                    axs[1].grid()
                    axs[1].set(xlabel='[substrate]',ylabel='residuals')
                    axs[1].yaxis.set_ticks(np.arange(min(residuals)-10,max(residuals)+10,10))
                    
                    plt.show()



# In[ ]:




