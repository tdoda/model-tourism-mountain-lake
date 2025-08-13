# -*- coding: utf-8 -*-
import os
import sys
import netCDF4
import numpy as np
import xarray as xr
from datetime import datetime, timedelta


class COMP:
    def __init__(self,name_comp,conc_units):
        self.general_attributes = {
            "name": name_comp
        }
        
        self.dimensions = {
            'time': {'dim_name': 'time', 'dim_size': None}
        }

        self.variables = {
            'time': {'var_name': 'time', 'dim': ('time',), 'unit': 'seconds since 1970-01-01 00:00:00', 'long_name': 'time'},
            'Conc': {'var_name': 'Conc', 'dim': ('time',), 'unit': conc_units, 'long_name': 'concentration'},
        }
        
        self.data={}
        
    def set_time_series(self,param_sim,param_process):
        
        tdate=np.arange(param_sim["Start_date"],param_sim["End_date"],timedelta(hours=param_sim["Time_step"]))
        tnum=tdate.astype("datetime64[ns]").astype(np.float64)*1e-9
        
        self.data["time"]=tnum
        self.data["tdate"]=tdate
        self.data["conc"]=np.full(tnum.shape,np.nan)
        self.data["conc"][0]=param_process["Initial_concentration"]
        


                                
                        
        