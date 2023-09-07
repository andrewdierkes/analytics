#!/usr/bin/env python
# coding: utf-8

# In[212]:


import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.cm as cm


class plate_reader():
    
    
    def wl_sample_plotter(excel_string, molecule, concentration_unit):
        '''Clean, transpose and plot dataframe
        --------------
        Parameters:
        excel_string = string literal of excel file including .xlsx
        molecule = string literal describing the molecule, usually one word
        concentration_unit = string literal IE mM or uM
        
        Returns a two dataframes, a) transposed and b) normal dataframe from excel, just processed'''
        
        #read excel file
        df_excel = pd.read_excel(excel_string)
        
        #clean data
        df_excel.index = [var for var in df_excel.iloc[:,0]]
        df_excel.drop(columns=[0],inplace=True)
        df_excel.rename(index={df_excel.index[0]:'grouping'},inplace=True)
        
        
        #transpose data and gain index for different titrations
        df = df_excel.transpose()
        regex = df_excel.iloc[0].unique()
        df.index.set_names(molecule,inplace=True)

        i = 0 

        idx_list = []
        for var in regex:
            idx = df[df.iloc[:,0] == var].index.to_list()
            idx_list.append(idx)

            i += 1

        #plot graphs
        for var in idx_list:

            _min = min(var)
            _max = max(var)

            idx_min = _min-1


            fig, ax = plt.subplots()

            while idx_min < _max:
                name = str(df.iloc[idx_min,0]) + '_IT_' + str(df.iloc[idx_min,1]) + f'{concentration_unit}_' + str(df.iloc[idx_min,2])
                ax.plot(df.columns[4:(len(df.columns)-1)],df.iloc[idx_min,4:(len(df.columns)-1)],label=name)
                ax.set(xlabel='wv',ylabel='RFU',title=f'{df.iloc[idx_min,0]} Gain for {df.iloc[0,1]} {concentration_unit} {df.index.name}')
                ax.legend()
                ax.grid()

                idx_min += 1
        return df, df_excel
    
    def wl_dataset_plotter(df):
        '''plots all samples on one graph
        ---------------
        Parameters:
        df = df_b, dataframe b) from wl_sample_plotter (2nd dataframe)
        '''

        unique = []
        for var in range(len(df.columns)):
            unique_name = str(df.iloc[0,var]) + '_IT_' + str(df.iloc[1,var]) + '_mm_' + str(df.iloc[2,var]) + '_rep'
            unique.append(unique_name)

        colors = cm.rainbow(np.linspace(0,1,len(unique)))
        col_i = [i+1 for i in range(len(df.columns))]
        plot = zip(col_i,colors,unique)


        fig, ax= plt.subplots()

        for i,c,u in plot: 
            ax.plot(df.index[4:len(df.index)],df.iloc[4:len(df.index),i-1],label=unique[i-1],color=c)
            ax.legend(bbox_to_anchor=(1.01,1.01))
            ax.grid()
            ax.set(xlabel='wv',ylabel='RFU',title=f'{df.columns.name}')
        
        return unique
        
    
    def max_emission(df,unique):
        '''find the max emission and it's associated wavelength and plot them
        -------------------
        Parameters:
        df = df_b, dataframe b) from wl_sample_plotter
        unique = list of unique names from wl_dataset_plotter'''
        
        #find max_emission values
        max_list = []

        for var in range(len(df.columns)):
            _max = df.iloc[3:(len(df.index)-1),var].max()
            max_list.append(_max)

        idx_list = []
        i = 0
        while i < len(df.columns):
            for var in max_list:
                idx = df[df.iloc[:,i] == var].index.to_list()
                idx_list.append(idx[0])
                i += 1
        color = cm.rainbow(np.linspace(0,1,len(unique)))
        max_emission = zip(unique,idx_list,max_list,color)

        df_emission = pd.DataFrame()
        df_emission.insert(0,'assay_name',unique)
        df_emission.insert(1,'max_emission',max_list)
        df_emission.insert(2,'max_emission_wv',idx_list)
        display(df_emission)

        fig, ax = plt.subplots()
        for u,i,m,c in max_emission:
            if i > 545:
                print('irrelevant',u)
                pass
            else:
                ax.scatter(int(i),int(m),label=u,color=c)
                ax.grid()
                ax.legend(bbox_to_anchor=(1.05,1.05))
                ax.set(xlabel='max_rfu',ylabel='wv_max_rfu',title=f'peak emission and it\'s wavelength for {df.columns.name}')

    

