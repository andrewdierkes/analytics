#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd


# In[1]:


def hct_dataframe(data):
    '''takes dataset_dict and returns a pd.dataframe
    ---------------
    Parameters:
    data = dataset_dict dictionary in format of {assay_name:[hct_values]}
    
    Returns:
    dataframe'''
    
    df_dataset_list = []
    for key, value in data.items():
       
        df = pd.DataFrame(value)

        #format dataframe
        df_tilt = df.transpose()
        df_tilt.rename(index={0:f'{key}'},inplace=True)
        df_tilt.columns.set_names('hct_read_#',inplace=True)

        df_dataset_list.append(df_tilt)

    #combine assays 
    df_dataset = pd.concat(df_dataset_list,axis=0)
    df_dataset.insert(0,'assay_name',df_dataset.index)
    
    return df_dataset

def hct_stat(data, start, stop):
    '''find the average, stdev & cv for rows with a start value and stop value
    Parameters
    ----------
    data = dataframe from hct.hct_dataframe
    start = col# where the first hct read is 
    stop = col# where the last hct head is'''
    
    hct_avg = []
    hct_std = []

    
    for row in range(len(data.index)):
    
        avg = data.iloc[row,start:stop].mean()
        hct_avg.append(avg)
    
        std = data.iloc[row,start:stop].std()
        hct_std.append(std)
    
    data['avg'] = hct_avg
    data['std'] = hct_std
    data['cv'] = round(((data['std']/data['avg'])*100),3)
    
    return data

def positive_value(df, column):
    '''iterate thru a dataframe's column and change the values to positive if negative
    Parameters
    ----------
    dataframe = dataframe
    column = string, column name'''
    
    df[column] = df[column].apply(lambda cv: cv * -1 if cv < 0 else cv)
    return df

# In[ ]:




