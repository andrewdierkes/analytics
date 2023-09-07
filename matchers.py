    #!/usr/bin/env python
# coding: utf-8

# In[1]:


import re


####################################NOT FINISHED FUNCTIONS

#                 def check_value(val1,val2):
#                     '''raises value error if val1==val2'''
#                     if val1 == val2:
#                         raise ValueError('Values are equal')


    #             def add_sub(concentration,substitute,tsv):
    #                 '''find and substitue values using function in a for loop iterating over a zip of two lists and return the values which changed

    #                 IE: 
    #                 fcon = ['con1','con2']
    #                 rcon = [1,2]

    #                 for f,r in zip(fcon,rcon):
    #                     add_sub(f,r,tsv)

    #                 Parameters
    #                 ----------
    #                 concentration = raw string literal with concentration placeholder (IE con3)
    #                 substitute =  holding string literal with concentration you'd like to substitute
    #                 tsv = string literal

    #                 Returns tsv with updated value'''

    #                 con_finder = re.compile(concentration)
    #                 con_match = con_finder.findall(tsv)

    #                 try:
    #                     con_matcher = con_match[0]
    #                     con_update = re.sub(con_matcher,substitute,tsv)

    #                     return con_update

    #                 except IndexError:
    #                     print(f'did not match {concentration} to {substitute}')
    #                     pass


#FOR TIME COMPONENT
    #     def spdg(data,time=True):
    #         '''find spdg data                    
    #         ---------------------------
    #         Parameters:
    #         data = string of data (file.read())
    #         time = whether or not to include time component, default = True

    #         Returns spdg and time... if there is processing for time, use 'if time' in script to prevent errors if time set to false
    #         '''


    #         #spdg_finder = re.compile(r'SOURCE:INS,ID:[0-9]+,NAME:XMSG_SPECTRO_PRINT, PAYLOAD:{dataIndex=[0-9]+, integrationTimeUs=[0-9]+.[0-9]+, delayTimeUs=[0-9]+, excitationTimeUs=[0-9]+, excitationMode=[0-9]+, averageCount=[0-9]+, numberPixels=[0-9]+, windowStartWavelength=[0-9]+.[0-9]+ ,pixelData:.*')
    #         if time:
    #             spdg_finder = re.compile(r'\d+_\d+.\d+	EVENT:\{SOURCE:INS,ID:[0-9]+,NAME:XMSG_SPECTRO_PRINT.*[0-9]\]}}')
    #         else:
    #             spdg_finder = re.compile(r'SOURCE:INS,ID:[0-9]+,NAME:XMSG_SPECTRO_PRINT.*[0-9]\]}}')                             
    #         spdg = spdg_finder.findall(data)
    #         #spdg_data = ' '.join(spdg)

    #         return spdg, time



# In[2]:


class matchers():
    
    def unique(list1):
        '''used to reduce repitions in lists'''
        list_set = set(list1)

        # convert the set to the list
        unique_list = (list(list_set))                    
        return unique_list

    
    def des(data):
        '''find des from @ codes
        ---------------------------
        Parameters:
        data = string of data (file.read())'''
        
        finder = re.compile('(?<=@DES: ).*')
        matcher = finder.findall(data)
        return matcher[0]
    
    def name(data):
        '''finds name using @name
        Parameters
        ----------
        data = string of data (file.read())'''
        
        #\d?.*[a-zA-Z]?.*_\d?.*[a-zA-Z]?.*'
        finder = re.compile(r'(?<=@name ).*')
        matcher = finder.findall(data)
        return matcher[0]
    
    def numberpixel(data):
        '''finds the number of wavelengths scanned
        ---------------------------
        Parameters:
        data = string of data (file.read())'''

        finder = re.compile('(?<=numberPixels=)[0-9]+')
        matcher = finder.findall(data)
        #keep one numberpixel and convert to int
        match = int(matcher[0])
        return match
    def stepsize(data):
        '''finds the number of steps/spatial scans
        ---------------------------
        Parameters:
        data = string of data (file.read())'''
        finder = re.compile('(?<=stepSize).*[0-9]+')
        matcher = finder.findall(data)
        #find stepsize, remove spaces and keep number
        matcher_final = matcher[0].split(' ')
        match = []
        for var in matcher_final:
            try:
                match.append(int(var))
            except:
                pass
        return match[0]


    def numberofsteps(data):
        '''finds the number of steps/spatial scans
        ---------------------------
        Parameters:
        data = string of data (file.read())'''

        finder = re.compile('(?<=numberOfSteps).*[0-9]+')
        matcher = finder.findall(data)
        #find stepsize, remove spaces and keep number
        matcher_final = matcher[0].split(' ')
        match = []
        for var in matcher_final:
            try:
                match.append(int(var))
            except:
                pass
        return match[0]

    def sscparam(data,num_it):
        '''find wl_start, wl_end in sscparam using a string
        ---------------------------
        Parameters:
        data = string of data (file.read())
        num_it = number of ITs
        '''

        finder = re.compile(r'(?<=SSCPARAM )[0-9]+.*')
        matcher = finder.findall(data)
        matcher_list = matcher[0].split(' ')

        wl_start = int(matcher_list[5])
        wl_end = int(matcher_list[6])

        print(f'This assay had wavebands starting at {wl_start} and ending at {wl_end}')
        print(f'This assay used {num_it} ITs') 
        return matcher, wl_start, wl_end

    def occonfig(data):
        '''collect start information on the optical channel edges within strip
        ---------------------------
        Parameters:
        data = string of data (file.read())
        '''
        finder = re.compile(r'(?<=OCCONFIG [0-9] )[0-9]+')
        #find unique channels in data
        matcher = sorted([int(var) for var in matchers.unique(finder.findall(data))])

        return matcher   

    def it(data):
        '''find unique IT                    
        ---------------------------
        Parameters:
        data = string of data (file.read())
        '''

        finder = re.compile(r'(?<=integrationTimeUs=)[0-9]+')
        matcher = matchers.unique(finder.findall(''.join(data)))
        it_float = sorted([float(var) for var in matcher],reverse=True)
        return it_float
    
    def it_number(data):
        '''find number of its per assay
        ---------------------------
        Parameters:
        data = string of data (file.read())'''
        
        finder = re.compile(r'(?<=numberOfSpectraPerPosition).*[0-9]+')
        matcher = finder.findall(data)
        matcher_final = matcher[0].split(' ')
        match = []
        for var in matcher_final:
            try:
                match.append(int(var))
            except:
                pass
        return match[0]
    
    def spdg(data):
        '''find spdg data                    
        ---------------------------
        Parameters:
        data = string of data (file.read())
        '''
        
        #spdg_finder = re.compile(r'SOURCE:INS,ID:[0-9]+,NAME:XMSG_SPECTRO_PRINT, PAYLOAD:{dataIndex=[0-9]+, integrationTimeUs=[0-9]+.[0-9]+, delayTimeUs=[0-9]+, excitationTimeUs=[0-9]+, excitationMode=[0-9]+, averageCount=[0-9]+, numberPixels=[0-9]+, windowStartWavelength=[0-9]+.[0-9]+ ,pixelData:.*')
        spdg_finder = re.compile(r'SOURCE:INS,ID:[0-9]+,NAME:XMSG_SPECTRO_PRINT.*[0-9]\]}}')                             
        spdg = spdg_finder.findall(data)
        #spdg_data = ' '.join(spdg)
        
        return spdg
            
    def spd(data):
        '''find pixel data                    
        ---------------------------
        Parameters:
        data = string of data (file.read())
        '''
        
        finder = re.compile(r'(?<=TYPE:SPD, STATUS:OK_STATUS, PAYLOAD).*[0-9]+]')
        matcher = finder.findall(data)
        
        return matcher
        
    def pixeldata(data):
        '''find pixel data                    
        ---------------------------
        Parameters:
        data = string of data (file.read())
        '''
        
        finder = re.compile(r'(?<=pixelData: \[).*[0-9]')
        matcher = finder.findall(data)
        
        return matcher
    
    def qos(data):
        '''find qos data                    
        ---------------------------
        Parameters:
        data = string of data (file.read())
        '''
        finder = re.compile(r'(?<=QOS Optical Signal: )[0-9]+.[0-9]+') 
        matcher = finder.findall(data)
        qos_float = [float(var) for var in matcher]
        
        return qos_float 
    
    def ir(data):
        '''find ir values from tsv file; only searching for those associated with HCT Value block
        ---------------------------
        Parameters:
        data = string of data (file.read())'''

        #? = optional 0 or 1
        finder = re.compile(r'(?<= \| IR = )-?[0-9].[0-9]+')
        matcher = finder.findall(data)

        ir_float = [float(var) for var in matcher]

        return ir_float
    
    def cyan(data):
        '''find cyan values from tsv file; only searching for those associated with HCT Value block
        ---------------------------
        Parameters:
        data = string of data (file.read())'''

        #? = optional 0 or 1
        finder = re.compile(r'(?<= \| CYAN = )-?[0-9].[0-9]+')
        matcher = finder.findall(data)

        cyan_float = [float(var) for var in matcher]

        return cyan_float  
    
    def hct_percent(data):
        '''find hct % from tsv file'''

        finder = re.compile('(?<=HCT Value  \|      = )\d+.\d+')
        matcher = finder.findall(data)

        hct_float = [float(var) for var in matcher]

        return hct_float
##UPDATE CLASS GOES HERE 

        

# In[ ]:






# In[ ]:




