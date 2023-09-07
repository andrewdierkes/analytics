#!/usr/bin/env python
# coding: utf-8




# In[1]:


import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from itertools import combinations
from statistics import mean
import pingouin
import seaborn as sns
from IPython.display import display
from scipy import stats

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn import preprocessing


class visualize():

    def highlight_cell(val,color='blue'):
        
        '''Visualize in posthoc_comparisons which differences are statistically significant
        color='string', default is blue'''
        
        color = color if val < 0.05 else ''

        return 'background-color: {}'.format(color)
    
    def plotter_3d(df,label_loc,x_loc,y_loc,z_loc):
        
        '''3-d plot of dataframe using a dataframe with 4 columns minimum; 1x for each label, x,y,z.
        Parameters
        ----------
        df = dataframe
        label_loc = col # for label of datapoint
        x_loc,y_loc,z_loc = col # for each variable

        Returns
        -------
        NA'''

        label = df.iloc[:,label_loc].unique()

        colors = cm.rainbow(np.linspace(0,1, len(label)))
        _3d = zip(df.iloc[:,x_loc],df.iloc[:,y_loc],df.iloc[:,z_loc],label,colors)

        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(projection='3d')

        for x,y,z,l,c in _3d:
            ax.scatter(x,y,z,color=c,label=l+f': {round(x,2),round(y,2),round(z,2)}' ,s=50)
            ax.legend()
            ax.set(xlabel=f'{df.columns[x_loc]}',ylabel=f'{df.columns[y_loc]}',zlabel=f'{df.columns[z_loc]}')
            ax.text(x+.01,y+.05,z+.05,f'{round(x,2),round(y,2),round(z,2)}',fontsize=6)
            ax.set(title = f'Plotting {df.columns[x_loc]} v {df.columns[y_loc]} v {df.columns[z_loc]}') 


    


class statistics():

    def paired_t_normality(data, str_comparison):
        
        '''Use this function when visualizing a paired T test, share 'data' 
        with the function in format of df['col']
        'str_comparison' is a string describing the two groups compared
        
        Parameters
        ----------
        data = dataframe['col_name']
        str_comparison = string describing the two groups compared
        
        Returns
        -------
        NA'''

        plt.hist(data)
        df_norm = pingouin.normality(data, method='shapiro')

        display(df_norm)
        for var in df_norm.iloc[:, 2]:

            if var is False:

                plt.axvline(np.median(data), color='r', linewidth=4)
                plt.text(np.median(data) + .05, (data.count() / 5),
                         f'median value of differences between {str_comparison}:')
                plt.text(np.median(data) + 0.05, (data.count() / 5.5), round(np.median(data), 3))
                plt.title(f'Non-Gaussian Distribution of the Differences between {str_comparison}')
                plt.xlabel(f'difference in {str_comparison}')
                plt.ylabel('count')

            else:

                plt.axvline(np.mean(data), color='r', linewidth=4)
                plt.text(np.mean(data) + .05, (data.count() / 4.25),
                         f'mean value of differences between {str_comparison}:')
                plt.text(np.mean(data) + 0.05, (data.count() / 4.5), round(np.mean(data), 3))
                plt.title(f'Gaussian Distribution of the Differences between {str_comparison}')
                plt.xlabel(f'difference in {str_comparison}')
                plt.ylabel('count')

    def kruskal_array(df, groupby, feature):
            
        '''Convert a dataframe into a list of arrays divided by the groupby
        function (IV). A non-parametric statistical method which we use when 
        the data is not gaussian distributed. 
        
        Parameters
        ----------
        df = dataframe, groupby=col # of iv,
        feature = col # of dv
        
        Returns
        -------
        1) list of arrays seperated by groupby
        2) list of unique groups'''
        
        regex = df.iloc[:, groupby].unique()

        groupby_name = df.columns[groupby]

        # sort dataframe to align with regex iterator below(min,max... basically need groupings in one space)
        df.sort_values(by=[groupby_name], inplace=True)
        df.reset_index(drop=True, inplace=True)

        # display(df)
        display(df.head(), df.tail())

        print(f'grouped by {groupby_name}')

        # find indexes of regex patterns
        idx_list = []
        for var in regex:
            idx = df[df.iloc[:, groupby] == var].index.to_list()
            idx_list.append(idx)

        for var in idx_list:
            print(f'length of group = {len(var)}')

        # iterate thru df using idx_list
        kw_list = []
        kw_dict = {}
        for var in idx_list:
            if len(var) == 0:
                pass
            else:
                _min = min(var)
                _max = max(var)
                chunk = len(var)
                offset = 0
                print('min',_min,'max', _max+1, type(_min), 'length of data for group:',len(df.iloc[_min:_max+1]),'chunk',chunk)
                print('data used in a group:',df.iloc[_min:_max+1])
                kw_list.append(np.array(df.iloc[_min:(_max+1), feature]))
                kw_dict[f'{df.iloc[offset, groupby]}'] = np.array(df.iloc[_min:(_max+1), feature])

        return kw_list, regex

    def normality(df, groupby, feature):
        
        '''Returns two plots KDE & QQ plots for each group, used to analyze the distribution 
        
        Parameters
        ----------
        df=dataframe name, 
        groupby=col # different categories were comparing(IV),
        feature=DV col #,
        
        Returns
        -------
        NA'''

        import pandas as pd
        import pylab

        regex = df.iloc[:, groupby].unique()

        groupby_name = df.columns[groupby]

        # sort dataframe to align with regex iterator below(min,max... basically need groupings in one space)
        df.sort_values(by=[groupby_name], inplace=True)
        df.reset_index(drop=True, inplace=True)

        display(df.head(), df.tail())

        print(f'grouped by {groupby_name}')
        # find indexes of regex patterns
        idx_list = []
        for var in regex:
            idx = df[df.iloc[:, groupby] == var].index.to_list()
            idx_list.append(idx)

        # iterate thru df using idx_list
        cv_list = []
        i = 0
        for var in idx_list:
            if len(var) == 0:
                pass
            else:
                _min = min(var)
                _max = max(var)
                chunk = len(var)
                offset = 0
                print('min',_min,'max', _max+1, type(_min), 'length of data for group:',len(df.iloc[_min:_max+1]))
                print('data grouping:',df.iloc[_min:_max+1])
                # KDE & QQ plots
                plt.figure(figsize=(10, 5))
                plt.subplot(1, 2, 1, title=f"Grouping {''.join(str(regex[i]))}")
                sns.kdeplot(df.iloc[_min:(_max+1), feature])
                plt.subplot(1, 2, 2)
                stats.probplot(df.iloc[_min:(_max+1), feature], plot=pylab)
                plt.show()
                i += 1

    def regex_cv(df, groupby, feature, print_enabled=True):
       
        '''Regex_cv finds mean,meedian,stdev, and CV for a feature in a df grouped by groupby.
        
        Parameters
        ----------
        df = dataframe, 
        groupby = col # for where to iloc and search for the regex,
        feature = col # which you want to aggregate mean, stdev & cv
        print_enabled = True default, prints out operations
        
        Returns
        -------
        dataframe with statistics'''
        
        

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
        cv_list = []
        mean_list = []
        median_list = []
        stdev_list = []
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
                median = df.iloc[_min:(_max+1), feature].median()
                stdev = df.iloc[_min:(_max+1), feature].std()

                cv = (stdev / mean) * 100

                mean_list.append(mean)
                median_list.append(median)
                stdev_list.append(stdev)
                cv_list.append(cv)

        final = list(zip(cv_list, median_list, mean_list, stdev_list))
        df_final = pd.DataFrame(final, index=[regex], columns=['cv', 'mean', 'median', 'stdev'])
        
        if print_enabled:
            display(df_final)

        return df_final
    
    def ph_nemenyi(df,y_col,group_col,block_col):
        
        '''This function runs a nemenyi posthoc if the friedman test rejects it's null hypothesis, to determine the statistically different pairwise comparisons,
        
        Parameters
        ----------
        df = dataframe with data containing atleast 3 mandatory columns (use dataframe returned by sts_stat; y_col= col# with data (dependent_var ie CV), 
        group_col = col# for repeated samples across different IVs/group_col, 
        block_col = overarching groups(iv) which hold the repeated group_col
        
        Returns
        -------
        dataframe with post hoc comparisons of different blocks'''
        
        
        import scikit_posthocs as sp

        print('data operator on:')
        display(df.head())

        y_col = df.columns[y_col]
        group_col = df.columns[group_col]
        block_col = df.columns[block_col]


        print('DV:',y_col,'IV:',block_col,'subject_repeat:',group_col)

        df_nemenyi = sp.posthoc_nemenyi_friedman(df,y_col,group_col,block_col,melted=True)

        iv = df[block_col].unique()

        df_color = df_nemenyi.style.applymap(visualize.highlight_cell)
        display(df_color)
        
        return df_color
    
    def ph_dunn(kw_sts, regex_sts):

        '''Runs a post hoc test for a kruskal wallis test when you reject it's the null;
        takes 2 objects returned from kruskal array (df aka kw_sts & list aka regex_sts)...plug in and rip
        
        Parameters
        ----------
        kw_sts = list of arrays from kruskal_array
        regex_sts = list of regex from kruskal_array
        
        Returns
        -------
        dataframe with post hoc analysis'''

        import scikit_posthocs as sp

        sts_dunn = sp.posthoc_dunn(kw_sts, p_adjust='bonferroni')

        sts_dunn.set_index(regex_sts, inplace=True)
        sts_dunn.columns = regex_sts
        sts_dunn.columns.set_names('dunn_post_hoc', inplace=True)

        r = 0
        chunk_r = 1

        while r < len(sts_dunn.columns):

            i = 0
            chunk_i = 1

            while i < len(sts_dunn.index):
                p_val = sts_dunn.iloc[i, r]
                # print(r,i)

                if p_val < 0.05:
                    print(
                        f'there IS a statistically significant difference between {sts_dunn.index[i]} and {sts_dunn.columns[r]}')
                else:
                    print(
                        f'there is NOT a statistically significant difference between {sts_dunn.index[i]} and {sts_dunn.columns[r]}')

                i += chunk_i

            r += chunk_r
        
        

        return sts_dunn



# In[7]:


class chunk():

    def delist(args):
        '''delists a list of lists into 1 list'''
        delist = [var for small_list in args for var in small_list]
        return (delist)

    def slice_name(df, row, name_len, end_offset):
        
        '''Function slices strings (usually assay name) that are similar in nature (df row). It works to
        remove the ends (end_offset) which are different (usually rep number) to allow the UNIQUE function to work
        
        Parameters
        ----------
        df = dataframe to use
        row = the index of the column you want to slice
        name_len = the length of each assay name, they all should be the same
        end_offset = what position you want to end with
        
        Returns
        -------
        dataframe with unique names'''

        assay__name = [var for var in df.iloc[:, row]]

        unique_name = []

        for var in assay__name:
            if len(var) == name_len:
                unique_name.append(var[0:end_offset])

            else:
                print(file)
                print(f'labeling error for assay:{var}, length {len(var)}')
                for count, var2 in enumerate(var):
                    print(count, var2)

        return pd.DataFrame(unique_name)

    def unique(df_col, chunk=4):

        '''Function will iterate over a df_col and return only unique values using chunk_iterator.. so if you'd like the unique value of something occuring every 3 rows... use 3 as your chunk iterator
        
        Parameters
        ----------
        df_col = df col # to operator on'''

        def delist(args):
            '''delists a list of lists into 1 list'''
            delist = [var for small_list in args for var in small_list]
            return (delist)

        offset = 0

        number_list = [var for var in range(len(df_col))]

        dataset_array = []
        dataset_list = []

        while offset < len(number_list):
            i = number_list[offset:chunk + offset]
            _array = df_col.iloc[i].unique()

            dataset_array.append(_array)

            offset += chunk

        for var in dataset_array:
            dataset_list.append(var.tolist())

        unique = delist(dataset_list)

        return unique

    def chunk_cv(df_col, chunk=4):
        '''This function will iterate over a df_col and return the average and stdev, using chunk_iterator.. so if you'd like the average & stdev values occuring every 3 rows... use 3 as your chunk iterator'''

        offset_mean = 0
        offset_stdev = 0

        number_list = [var for var in range(len(df_col))]

        dataset_average = []
        dataset_stdev = []

        while offset_mean < len(number_list):
            i_mean = number_list[offset_mean:chunk + offset_mean]
            average = df_col.iloc[i_mean].mean(axis=0)

            dataset_average.append(average)
            # dataset_array.append(_array)

            offset_mean += chunk

        while offset_stdev < len(number_list):
            i_stdev = number_list[offset_stdev:chunk + offset_stdev]
            stdev = df_col.iloc[i_stdev].std(ddof=1)

            dataset_stdev.append(stdev)

            offset_stdev += chunk

        return dataset_average, dataset_stdev
    
        
    
    
    


# In[6]:


def ch_comparison(df, sort):
    
    '''A channel comparison
    
    Parameters
    ----------
    df = dataframe
    sort = column to sort values on
    
    Returns
    -------
    dataframe and list'''
    
    from adstat import statistics

    df.sort_values(by=[sort])

    ch_1400 = [var for var in range(0, len(df), 4)]
    ch_2200 = [var for var in range(1, len(df), 4)]
    ch_3000 = [var for var in range(2, len(df), 4)]
    ch_3800 = [var for var in range(3, len(df), 4)]

    ch_1400_df = []
    for var in ch_1400:
        _1400 = df.iloc[var]
        ch_1400_df.append(_1400)
    df_ch1400 = pd.DataFrame(ch_1400_df)
    # display(df_ch1400)

    ch_2200_df = []
    for var in ch_2200:
        _2200 = df.iloc[var]
        ch_2200_df.append(_2200)
    df_ch2200 = pd.DataFrame(ch_2200_df)
    # display(df_ch2200)

    ch_3000_df = []
    for var in ch_3000:
        _3000 = df.iloc[var]
        ch_3000_df.append(_3000)
    df_ch3000 = pd.DataFrame(ch_3000_df)
    # display(df_ch3000)

    ch_3800_df = []
    for var in ch_3800:
        _3800 = df.iloc[var]
        ch_3800_df.append(_3800)
    df_ch3800 = pd.DataFrame(ch_3800_df)

    df_channel_dataset = pd.concat([df_ch1400.reset_index(drop=True),
                                    df_ch2200.reset_index(drop=True),
                                    df_ch3000.reset_index(drop=True),
                                    df_ch3800.reset_index(drop=True)])

    spatial_stos = {}

    spatial_stos['cv_ch1400'] = regex_cv(df_ch1400, 2, 4)
    spatial_stos['cv_ch2200'] = regex_cv(df_ch2200, 2, 4)
    spatial_stos['cv_ch3000'] = regex_cv(df_ch3000, 2, 4)
    spatial_stos['cv_ch3800'] = regex_cv(df_ch3800, 2, 4)

    df_cv = []

    for key, value in spatial_stos.items():
        value.sort_values(by=['cv'], inplace=True)
        value.index.set_names(key, inplace=True)
        display(value)
        df_cv.append(value)

    return df_cv, spatial_stos

    

class strip_to_strip():

    def sts_combinations(df, groupby, feature, orderby):
        '''Function is meant to find the different combinations possible among specific groups (usually IV; groupby). However, to align the partitions(groups) we use an order_by before the combination is made, this is important for paired t tests & ANOVAs (as we'd like the comparisons the align). We want to do combinations for different groups on the same index
        
        Parameters
        ----------
        df = dataframe, 
        groupby=col# of the IV, 
        feature=col# of DV. 
        orderby = col# of indexer for each groupby category.
        
        Returns
        -------
        This function returns 2 items, you can use the returned dictionary of df (cv_final) in sts_stat
        to make assumptions on the data:
        1) dictionary with groupby: dataframe of combinations & stats
        2) dictionary with statistics on the dataframes.'''

        from itertools import combinations
        from statistics import mean

        regex = df.iloc[:, groupby].unique()
        groupby_name = df.columns[groupby]
        orderby_name = df.columns[orderby]

        df.sort_values(by=[groupby_name, orderby_name], inplace=True)
        df.reset_index(drop=True, inplace=True)

        print('check these df insights to ensure the frame is sorted by groupby & the orderby which is repeated across the features')
        display(df)
        #display(df.head())
        #display(df.tail())

        idx_list = []
        for var in regex:
            idx = df[df.iloc[:, groupby] == var].index.to_list()
            idx_list.append(idx)

        combo_list = []


        for var in idx_list:
            print(var)
            _min = min(var)
            _max = max(var)
            print('min', _min, 'max',_max + 1, type(_min)), 'length of data for group:', len(df.iloc[_min:(_max+1),feature])
            print('data taken for group',df.iloc[_min:(_max+1),feature])
            combo = combinations(df.iloc[_min:(_max+1), feature], 2)
            combo_list.append(combo)

                # idxer =df[df.iloc[_min:(_max+1),feature] == var].index.to_list()
                # print(idxer)
        # combinations stored in a list of combo_object
        cv_final = {}
        dict_final = {}

        # for index labeling
        i = 0
        chunk = 1

        for var in combo_list:

            mean_list = []
            std_list = []
            cv_list = []
            pair_list = []

            for pair in var:
                mean = np.mean(pair)
                stdev = np.std(pair)
                cv = round((stdev / mean) * 100, 3)

                pair_list.append(pair)
                mean_list.append(mean)
                std_list.append(stdev)
                cv_list.append(cv)

            df_group = pd.DataFrame()
            df_group.insert(0, 'sts_combination', pair_list)
            df_group.insert(1, 'sts_average', mean_list)
            df_group.insert(2, 'sts_stdev', std_list)
            df_group.insert(3, 'sts_cv', cv_list)
            df_group.index.set_names(f'{regex[i]}', inplace=True)

            cv_final[f'{regex[i]}'] = df_group

            d_stats = {}

            d_stats['mean_cv'] = df_group['sts_cv'].mean()
            d_stats['median_cv'] = df_group['sts_cv'].median()

            df_stats = pd.DataFrame(d_stats, index=[f'{regex[i]}'])
            df_stats.index.set_names('iv', inplace=True)

            dict_final[f'{regex[i]}'] = d_stats

            print('below is the mean and median cv values for the entire combination df')
            display(df_stats)

            i += chunk

            print('total number of combinations executed:', len(df_group))

            print('below is the head portion of the combination df')
            #display(df_group.head())
            display(df_group)
        return cv_final, dict_final

    def sts_stat(sts_cv_dict, column, frac=0.1, random_state=123):
        
        '''Reveals information about the distribution, homoscedasticity, & plots KDE & qq plots.
        It also combines the dictionary returned by the function sts_combinations and adds a new column for IV.
        
        Parameters
        ----------
        sts_cv_dict = dictionary; 1) from sts_combinations (dictionary of dataframes)
        column= string; column name to act on
        frac= float <= 1; % of population to use as sample size, default @ 0.1
        random_state= seed to ensure the function takes the same sample every time
        
        Returns
        -------
        df_final, dataframe '''

        # take dictionary of dataframes, sample the population, add an additional column with IV and concat into a workable df
        # perform statistical analysis on df

        sample_collection = []

        for key, value in sts_cv_dict.items():
            sample = value[column].sample(frac=frac, random_state=random_state)
            repeat_sample = len(sample)
            

            print('number of samples taken:', repeat_sample)

            df = pd.DataFrame(sample)
            df.insert(0, 'iv', [key for var in range(repeat_sample)])
            df.reset_index(drop=True, inplace=True)

            sample_collection.append(df)

        df_final = pd.concat([var for var in sample_collection])
        df_final.reset_index(drop=True,inplace=True)
        display(df_final)
        
        df_norm = pingouin.normality(df_final, dv=column, group='iv', method='shapiro')
        df_norm.index.set_names('shapiro_normality', inplace=True)

        df_homosc = pingouin.homoscedasticity(df_final, dv=column, group='iv', method='levene')
        df_homosc.index.set_names('levene_homoscedasticity', inplace=True)

        display(df_norm)
        display(df_homosc)
        display(statistics.normality(df_final, 0, 1))

        return (df_final)


def cartesian_pdt(df, groupby, feature1, feature2):
    '''Use this function to return the cartesian product of two sets(feature1 and 2) for each of the unique categories in the groupby column
    df = target df, groupby = IV/column with labels, feature1/2 = data to perform the cartesian product on'''

    regex = df.iloc[:, groupby].unique()
    groupby_name = df.columns[groupby]

    df.sort_values(by=[groupby_name], inplace=True)
    df.reset_index(drop=True, inplace=True)

    idx_list = []

    for var in regex:
        idx = df[df.iloc[:, groupby] == var].index.to_list()
        idx_list.append(idx)

    # print(regex)
    # print(idx_list)
    cartesian_list = []

    # take the cartesian product for in rows defined by groupby min/max and using feature1/2 for columns
    for var in idx_list:
        _min = min(var)
        _max = max(var)
        print('min',_min,'max', _max+1, type(_min), 'length of data for group:',len(df.iloc[_min:_max+1]))
        print('data taken for feature1',df.iloc[_min:(_max+1),feature1])
        print('data taken for feature1',df.iloc[_min:(_max+1),feature2])
        _product = list(product(df.iloc[_min:(_max+1), feature1], df.iloc[_min:(_max+1), feature2]))
        # print(combo)

        cartesian_list.append(_product)

    # store individual dataframes for each groupby
    dict_final = {}
    dict_stat = {}

    # for index labeling
    i = 0
    chunk = 1

    # apply correction (division) for all cartesian products & create dataframe for each individual groupby category
    for var in cartesian_list:
        # print(var)

        correction_list = []

        for pair in var:
            print(pair)
            correction = pair[0] / pair[1]
            correction_list.append(correction)

        # print(correction_list)
        df_var = pd.DataFrame(correction_list, columns=['correction'])
        df_var.index.set_names(f'{regex[i]}', inplace=True)

        print('below is the head portion of the combination df')
        display(df_var)

        dict_final[f'{regex[i]}'] = df_var

        # create dictionary & df with stats on correction values
        d_stats = {}
        d_stats['mean_correction_val'] = df_var['correction'].mean()
        d_stats['median_correction_val'] = df_var['correction'].median()
        d_stats['stdev_correction_val'] = df_var['correction'].std()
        d_stats['cv_correction_val'] = round((d_stats['stdev_correction_val']/d_stats['mean_correction_val'])*100,3)

        df_stat = pd.DataFrame(d_stats, index=[f'{regex[i]}'])
        display(df_stat)
        # add to outside dictionary
        dict_stat[f'{regex[i]}'] = d_stats

        i += chunk

        print('total number of combinations executed:', len(df_var))

    # df = pd.DataFrame(dict_stat)
    print(dict_stat)
    return dict_final, dict_stat

#least squares regression:
class models():
    

    #regression function:
    def lr_train(df, x, y,train_size=0.7):
        '''linear regression for a single x variable
        df = dataframe with x & y, x = col # of IV, y = col # of DV
        train_size = decimal number for % of data to use to train'''
        
        print('DV:',df.columns[x],'IV:',df.columns[y])

        #create array and reshape for x value
        regress_x = np.array(df.iloc[:,x]).reshape(-1,1)
        
        X_train,X_test,y_train,y_test = train_test_split(regress_x,df.iloc[:,y],test_size=1-train_size,random_state=32)

 
        regr = LinearRegression()
        regr.fit(X_train,y_train)
        
            
        r2_val = round(regr.score(X_train,y_train), 4)
        intercept = round(regr.intercept_, 4)
        slope = regr.coef_
        
        
        
        y_pred = regr.predict(X_test)
        
        fig, ax = plt.subplots()
        ax.scatter(X_test,y_test,c='b',label='true_data')
        ax.scatter(X_test,y_pred,color='k',label='predicted_data')
        ax.grid()
        ax.legend()
        ax.set(xlabel=f'{df.columns[x]}',ylabel=f'{df.columns[y]}',title='Actual test data v. Regression model predicted values')
        equation = f'y = {slope}x + {intercept}'
        
        #MAE - measure of errors bt paired observations (predicted v observed)
        print('evaluate model:')
        mae = mean_absolute_error(y_true=y_test,y_pred=y_pred)
        #squared True returns MSE value, False returns RMSE value.
        mse = mean_squared_error(y_true=y_test,y_pred=y_pred) #default=True
        rmse = mean_squared_error(y_true=y_test,y_pred=y_pred,squared=False)

        print("MAE:",mae)
        print("MSE:",mse)
        print("RMSE:",rmse)
        
        return f'Regression Model: {equation} with an R2 value of {r2_val}'

    def lr_simple(df, x, y, title, x_label, y_label, ytick_adjust=False, ytick_spacer=0, ytick_interval=0):
        '''Regression function using all data, no predictions are made. Returns the linear equation
        ------------------------
        Parameters:
        df = dataframe
        x = col # of IV
        y = col # of DV
        title = string title
        x_label = string label
        y_label = string label
        ytick_adjust = if you want to adjust the y ticks change to True (default False)
        y_spacer = integer you subtract from min and add to max values of the y dataset. If you're trying to show a zoomed out view of the data (default 0)
        ytick_interval = spacing between the y ticks on the y axis (default 0)'''


        #create array and reshape for x value
        x=df.iloc[:,x]
        y=df.iloc[:,y]
        regress_x = np.array(x).reshape(-1,1)


        #allow model to be used for predictions
        global model 
        model = LinearRegression().fit(regress_x,y)
        r2_val = round(model.score(regress_x,y), 4)
        global intercept
        intercept = round(model.intercept_, 4)
        global slope
        slope = model.coef_

        equation = f'y = {slope}x + {intercept}'

        #plot data
        fig, ax = plt.subplots()
        plt.scatter(x, y, color='r')
        ax.grid()
        #ax.set(xlabel=x_label, ylabel=y_label,title=f'Linear Regression {equation} for {title}')
        ax.set(xlabel=x_label, ylabel=y_label,title=f'{title}')
        
        #change distance between y ticks 
        if ytick_adjust: 
            y_min = min(y)
            y_max = max(y)
            ax.yaxis.set_ticks(np.arange(y_min-ytick_spacer,y_max+ytick_spacer,ytick_interval))


        print(f'Regression Model: {equation} with an R2 value of {r2_val}')

        return equation



# In[ ]:




