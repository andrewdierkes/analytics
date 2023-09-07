#!/usr/bin/env python
# coding: utf-8

# In[2]:
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


class spatial():
    
    def wl_index(wl_start,wl_stop,num_pixel):
        '''create wavelength index using wl_start, wl_stop and num_pixel
        ---------------------------
        Parameters:
        wl_start = integer from sscparam.matchers()
        wl_end = integer from sscparam.matchers()
        num_pixel = integer from matchers.numberpixel'''

        #print(wl_stop,wl_start)
        wl_step = (wl_stop-wl_start)/num_pixel
        offset = wl_start-wl_step 

        #create index starting at the wl_start and ending one before end
        wl_index = []
        while offset < (wl_stop-wl_step):

            i = round(float(offset + wl_step),2)
            wl_index.append(i)
            #wl_index_within.append(i)
            offset += wl_step

        print('length of index=',len(wl_index))
        return wl_index  

    def delist(args):
        '''delists a list of lists into 1 list'''
        delist = [var for small_list in args for var in small_list]
        return(delist) 


    def spatial_index(ch_pos, step_size, num_steps):
        '''iterate thru spatial positions using ch_pos we start with, step size and num_steps from assay.. range(start,stop,interval)
        ----------------------
        Parameters:
        ch_pos = integer start of optical channel
        step_size = integer stepsize from matchers.stepsize
        num_steps = number of steps from matchers.numberofsteps'''
        spatial_lister = []
        for i in range(ch_pos, ch_pos+(step_size*num_steps), step_size):
                        spatial_lister.append(i)
        return spatial_lister

    def data_divider(data, number_of_steps):
        '''generator function we divide the data from pixel_data_list into their appropriate channels, 
        using the number of steps in each channel
         ---------------------------
         Parameters:
         data = list of all spatial scans len = total_spatial_scans
         number_of_steps = integer matchers.numberofsteps()
         '''

        for i in range(0, len(data), number_of_steps):
            yield data[i:i+number_of_steps]

    def spd_plotter(df, wl_index, spatial_idx_list, graph_channel,show):
        '''find max RFU output, peak emission wavelength. From these we will graph and create a dictionary
        ------------------------
        Parameters
        df = spatial dataframe for channel
        wl_index = list of wavelengths scanned in assay
        spatial_idx_list = list of spatial positions scanned in assay 
        wl_index and spatial_idx_list should have same length as df
        graph_channel = list of assay names
        show= = True or False, if you want to display plots (True)
        all lists are lists of lists and require an offset and chunk method'''

        ch_max = [df.max(axis=1)]
        ch_spatial = spatial_idx_list
        
        
        #iterate thru max list & find max value (max_rfu)
        for var in ch_max:
            max_rfu = (max(var))
            #print(max_rfu)
            
        #find max rfu value and index to find the spatial position in the column names
        search_spatial = df.isin([max_rfu]).any()
        max_spatial_position = search_spatial[search_spatial].index.values[0]
        

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

            fig.suptitle(f'Assay {graph_channel} Spatial RFU & Max Emission WL at {round(max_emission_wl, 3)}')


            #plot spatial vs RFU at max wl & maxima via scipy
            axs[0].scatter(ch_spatial, y)
            axs[0].scatter(max_spatial_position, max_rfu, color = 'r', s = 30, marker = 'D', label = 'maxima')
            axs[0].set(xlabel='Spatial Position', ylabel='RFU')
            axs[0].legend(bbox_to_anchor=(1.05,1.05))
            axs[0].grid()

            #plot wavelength vs RFU
            for var in range(len(df.columns)):
                spatial = df.iloc[:,var]

                axs[1].scatter(df.index,spatial,label=df.columns[var])
                axs[1].legend(bbox_to_anchor=(1.05,1.05))
                axs[1].set(xlabel='wavelength',ylabel='RFU')


            fig.tight_layout()
        else:
            pass


        


        data_dict = {}

        data_dict['assay_name'] = graph_channel
        data_dict['max_rfu'] = max_rfu
        data_dict['average_8'] = average 
        data_dict['number_spatial_positions_for_average_8e'] = num_spatial_pos 
        data_dict['max_emission_wl'] = max_emission_wl
        data_dict['max_emission_spatial'] = max_spatial_position

        return data_dict
    
    
    def spatial_plotter_3d(df):
        '''takes a dataframe with an index(wavelength) and column(spatial position) and datapoints (RFU) and creates a 3d plot
        Parameters
        ----------
        df = dataframe with correct set up index = wl, columns = spatial, data = rfu
        Returns 3d plot'''

        #create a list which matches a repeated number of columns with the number of wavelengths taken at that position (numpixel) 
        spatial_idx_match = []
        for var in df.columns:
            repeat = np.repeat(var,len(df.index)).tolist()
            #print(repeat)
            spatial_idx_match.append(repeat)

        #create a list which multipies the waveband range by the number of spatial positions
        wl_repeat = [var for var in df.index]*len(df.columns)




        #gain access to all data points organized in order of spatial pos: wlstart to wl end, and so on
        raw_data = []
        offset_col = 0
        chunk = 1
        while offset_col < len(df.columns):

            offset_row = 0

            while offset_row < len(df.index):

                data = df.iloc[offset_row,offset_col]
                raw_data.append(data)


                offset_row += chunk

            offset_col += chunk

        #zip the data and delist the spatial index
        _3d = zip(wl_repeat,spatial.delist(spatial_idx_match),raw_data)

        print(len(raw_data),len(repeat),len(spatial_idx_match))

        #plot data
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(projection='3d')

        for x,y,z in _3d:
            ax.scatter(x,y,z)
            ax.grid()
            ax.set(xlabel='wavelength',ylabel='spatial_position',zlabel='RFU',title=df.index.name)

        plt.show()


