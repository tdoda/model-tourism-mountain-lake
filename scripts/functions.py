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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def test_param(param_dict,param_names):
    """Function test_param

    Returns an error if one of the specified parameters is not in the dictionary.

    Inputs:
        
    Outputs: None
    """
    for p in param_names:
        if p not in param_dict.keys():
            raise Exception("Parameter {} is missing in the file",p)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def create_time(param_sim):
    """Function load_parameters

    Create time variables from simulation parameters.

    Inputs:
        param_sim (dict): parameters of the simulation (from csv file)
        
    Outputs:
        tnum (numpy array of float): time values as number of seconds since 01.01.1970
        tdate (numpy array of datetime): datetime values corresponding to tnum
        
    """

    tdate=np.arange(param_sim["Start_date"],param_sim["End_date"],timedelta(hours=param_sim["Time_step"])).astype(datetime)
    tnum=tdate.astype("datetime64[ns]").astype(np.float64)*1e-9
    
    return tnum, tdate

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_outflow(outflow_type,param_lake,tnum):
    """Function get_outflow

    Create time series of outflow discharges based on the specific data type (constant, function or file).

    Inputs:
       outflow_type (str): type of outflow data, options: 
           "constant" = use constant discharge from param_lake,
           "function_name" use the temporal function "function_name" defined in functions.py, 
           "file" = use data file lakename_outflow.csv in lakes folder (with lakename the name of the specific lake)
       param_lake (dict): lake parameters from csv file 
       tnum (numpy array of floats, shape (n,)): time values [seconds since 01.01.1970] 
        
    Outputs:
        Qout (numpy array of floats, shape (n,)): outflow discharge values [m3/s]
    """

    if outflow_type=="constant":
        test_param(param_lake,["Outflow_constant"])
        Qout=np.full(tnum.shape,param_lake["Outflow_constant"])
        
    elif outflow_type=="file":
        df=df=pd.read_csv(os.path.join("lakes",param_lake["Lake_name"],param_lake["Lake_name"]+"_outflow.csv"), sep=",",header=0,names=["Time","Qout"],parse_dates=["Time"]) 
        tnum_Q=df["Time"].values.astype("datetime64[ns]").astype(float)/1e9
        Qout=np.interp(tnum, tnum_Q, df["Qout"].values) # Linear interpolation with constant values outside data period
        
    elif outflow_type=="Qout_peak":
        test_param(param_lake,["Outflow_maximum","Outflow_start_date","Outflow_end_date","Outflow_max_date"])      
        t_start=param_lake["Outflow_start_date"].replace(tzinfo=timezone.utc).timestamp()
        t_end=param_lake["Outflow_end_date"].replace(tzinfo=timezone.utc).timestamp()
        t_max=param_lake["Outflow_max_date"].replace(tzinfo=timezone.utc).timestamp()  
        Qout=Qout_peak(tnum,param_lake["Outflow_maximum"],t_start,t_max,t_end)
        
    return Qout

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def Qout_peak(tnum,Qmax,t_start,t_max,t_end):
    """Function Qout_parabolic

    Create time series of outflow discharges as an increase from t_start to t_max and a decrease from t_max to t_end.

    Inputs:
        
    Outputs:
        Qout (numpy array of floats, shape (n,)): outflow discharge values [m3/s]
    """
    Qout=np.zeros(tnum.shape)
    ratio1=(tnum[(tnum>=t_start)&(tnum<t_max)] - t_start) / (t_max - t_start)
    ratio2=(tnum[(tnum>=t_max)&(tnum<t_end)] - t_max) / (t_end - t_max)
    Qout[(tnum>=t_start)&(tnum<t_max)]=Qmax*(3*ratio1**2-2*ratio1**3)
    Qout[(tnum>=t_max)&(tnum<t_end)]=Qmax*(1-(3*ratio2**2-2*ratio2**3))
    
    
    return Qout
