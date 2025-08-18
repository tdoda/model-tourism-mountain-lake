"""
Functions used for the two-box lake model.

Author: T. Doda, Faculty of Geoscience and Environment, Institute of Earth Surface Dynamics, University of Lausanne
Contact: tomy.doda@unil.ch
Date: xx.dd.2025
"""

import os
import math
import numpy as np
import pandas as pd
from datetime import datetime, timedelta, timezone, UTC


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def load_parameters(file):
    """Function load_parameters

    Load model parameters from a csv file.

    Inputs:
        file (string): file name (including path).
        
    Outputs:
        values_dict (dictionary): values of the parameters from the csv file
        units_dict (dictionary): units of the parameters from the csv file
        ref_dict (dictionary): reference of the parameter values from the csv file
    """
    
    df=pd.read_csv(file, sep=",",index_col=0,header=0,dtype=str) 
    data_dict=df.to_dict()
    if "Values" in data_dict.keys():
        values_dict=data_dict["Values"]
        for key in values_dict.keys(): # Convert values to float if possible
            if "Units" in data_dict.keys():
                if data_dict["Units"][key]=="yyyymmdd":
                    values_dict[key]=datetime.strptime(values_dict[key],"%Y%m%d") 
                elif data_dict["Units"][key]=="yyyymmddHHMM":
                    values_dict[key]=datetime.strptime(values_dict[key],"%Y%m%d%H%M")
                elif data_dict["Units"][key]=="HHMM":
                    values_dict[key]=datetime.strptime(values_dict[key],"%H%M")
            if type(values_dict[key])==str and "--" in values_dict[key]: # Range
                range_str=values_dict[key]
                values_dict[key]=np.array([float(range_str[:range_str.find("--")]),float(range_str[range_str.find("--")+2:])])            
            else:
                try:
                    values_dict[key]=float(values_dict[key])
                except:
                    continue
    else:
        raise Exception("Parameter values are missing")
        
    if "Units" in data_dict.keys():
        units_dict=data_dict["Units"]
    else:
        print("*WARNING*: Parameter units are missing")
        units_dict=dict()
        
    if "Reference" in data_dict.keys():
        ref_dict=data_dict["Reference"]
    else:
        ref_dict=dict()
    
    
    return values_dict, units_dict, ref_dict



