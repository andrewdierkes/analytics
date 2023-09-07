#!/usr/bin/env python
# coding: utf-8

# In[2]:


import re
import pandas as pd


# In[1]:


def assay_name_divider(df,feature,regex,new_feature_name):
    '''meant to extract relevant information from a column and create a new column with information from regex
    Parameters
    -----------
    df = dataframe 
    feature = column # to search the regex in
    regex = string literal in format: r'(?<=)strip_design' to search with
    new_feature_name = string literal which names the new column
    
    Returns an updated dataframe'''
    pattern = re.compile(regex)
    
    finder_list = []
    for var in df.iloc[:,feature]:
        finder = pattern.findall(var)
        finder_list.append(finder[0])
        
    df.insert(feature+1, new_feature_name, finder_list)
    return df


# In[ ]:


class name_sort():
    
    def esterase_ul(ul=r'\d+(?=UL)'):
        ''' ul regex is defaulted, adjust if needed
        
        Parameters:
        -----------
        ul = string regex'''
        return re.compile(ul)
    
    def ph(ph=r'(?<=_)\d.\d'):
        '''ph regex is defaulted, adjust if needed
        
        Parameters:
        -----------
        ph = string regex'''
        return re.compile(ph)
    
    def sort(_path,func_matcher):
        '''iterate thru tsv filenames to find different esterate concentrations

        Parameters:
        _path = path.Path.cwd() object
        func_matcher = function object with re.compile('regex') returned

        Returns a dictonary of conc:[posixpath]'''

        filepath_dict = {}

        for _filepath in _path.iterdir():
            if _filepath.suffix != '.tsv':
                continue
            elif _filepath.suffix == '.tsv':

                finder = func_matcher.findall(str(_filepath))
                # esterase_finder = re.compile(r'\d+(?=UL)')
                # esterase_matcher = esterase_finder.findall(str(_filepath))


                try:
                    found = finder[0]
                    #esterase_conc = esterase_matcher[0]

                     #if conc key exists, append to list..

                    if filepath_dict.get(found):
                        #print(f'already a key named {esterase_conc}')
                        filepath_dict[f'{found}'].append(_filepath)

                    #if none create a new key with a value list
                    else:
                        filepath_dict[f'{found}'] = [_filepath]

    #                 if tsv_conc.get(esterase_conc):
    #                     #print(f'already a key named {esterase_conc}')
    #                     tsv_conc[f'{esterase_conc}'].append(_filepath)

    #                 #if none create a new key with a value list
    #                 else:
    #                     tsv_conc[f'{esterase_conc}'] = [_filepath]


                except IndexError or ValueError:
                    print('error in Index or Value type')
        return filepath_dict

