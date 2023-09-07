#!/usr/bin/env python
# coding: utf-8

# In[1]:
import re
import pandas as pd

def name_sort(_path,func):
    '''iterate thru tsv filenames to find different esterate concentrations
    
    Parameters:
    _path = path.Path.cwd() object
    
    Returns a dictonary of conc:[posixpath]'''

    esterase_list = []
    tsv_conc = {}

    for _filepath in _path.iterdir():
        if _filepath.suffix != '.tsv':
            continue
        elif _filepath.suffix == '.tsv':
            esterase_finder = re.compile(r'\d+(?=UL)')
            esterase_matcher = esterase_finder.findall(str(_filepath))


            try:
                esterase_conc = esterase_matcher[0]

                 #if conc key exists, append to list..
                if tsv_conc.get(esterase_conc):
                    #print(f'already a key named {esterase_conc}')
                    tsv_conc[f'{esterase_conc}'].append(_filepath)

                #if none create a new key with a value list
                else:
                    tsv_conc[f'{esterase_conc}'] = [_filepath]


            except IndexError or ValueError:
                print('error in Index or Value type')
    return tsv_conc

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


def add_sub(concentration,substitute,tsv):
    '''find and substitue values using function in a for loop iterating over a zip of two lists and return the values which changed

    IE: 
    fcon = ['con1','con2']
    rcon = [1,2]

    for f,r in zip(fcon,rcon):
        add_sub(f,r,tsv)

    Parameters
    ----------
    concentration = raw string literal with concentration placeholder (IE con3)
    substitute =  holding string literal with concentration you'd like to substitute
    tsv = string literal

    Returns tsv with updated value'''

    con_finder = re.compile(concentration)
    con_match = con_finder.findall(tsv)

    try:
        con_matcher = con_match[0]
        con_update = re.sub(con_matcher,substitute,tsv)

        return con_update

    except IndexError:
        print(f'did not match {concentration} to {substitute}')
        pass
    
            #EXAMPLE FOR ABOVE
#             fcon = ['con1','con2','con3','con4']
#             rcon = ['50mgdl','118mgdl','208mgdl','277mgdl']
#             con_data = zip(fcon,rcon)

#             updated_name = []
#             for f,r in con_data:

#                 update = add_sub(f,r,tsv)
#                 #only add variables that matched
#                 if isinstance(update,str):
#                     updated_name.append(update)

#             #find conc in updated name to add to index
#             concentration = []
#             for var in rcon:
                
#                 conc_finder = re.compile(var)
#                 conc_matcher = conc_finder.findall(updated_name[0])
#                 try:
#                     conc_matched = conc_matcher[0]
#                     concentration.append(conc_matched)
#                 except IndexError:
#                     print('no match for var in rcon when converting tsv to concentration')
#                     pass
                    



# In[ ]:




