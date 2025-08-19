# -*- coding: utf-8 -*-
"""
Compound class containing functions to setup, compute and export concentration values of a specific compound.

Author: T. Doda, Faculty of Geoscience and Environment, Institute of Earth Surface Dynamics, University of Lausanne
Contact: tomy.doda@unil.ch
Date: xx.dd.2025
"""


import os
import sys
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta, timezone


class COMP:
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    def __init__(self,name_comp,act):
        self.general_attributes = {
            "name": name_comp,
            "activity":act
        }
        
        self.dimensions = {
            'time': {'dim_name': 'time', 'dim_size': None}
        }

        self.variables = { 
            'tdate': {'var_name': 'tdate', 'dim': ('time',), 'unit': 'yyyy-mm-dd HH:MM:SS', 'long_name': 'datetime'},
            'time': {'var_name': 'time', 'dim': ('time',), 'unit': 'seconds since 1970-01-01 00:00:00', 'long_name': 'time'},
            'conc_epi_min': {'var_name': 'conc_epi_min', 'dim': ('time',), 'unit': 'g.m-3', 'long_name': 'minimum concentration in the epilimnion'},
            'conc_epi_max': {'var_name': 'conc_epi_max', 'dim': ('time',), 'unit': 'g.m-3', 'long_name': 'maximum concentration in the epilimnion'},
            'conc_hypo_min': {'var_name': 'conc_hypo_min', 'dim': ('time',), 'unit': 'g.m-3', 'long_name': 'minimum concentration in the hypolimnion'},
            'conc_hypo_max': {'var_name': 'conc_hypo_max', 'dim': ('time',), 'unit': 'g.m-3', 'long_name': 'maximum concentration in the hypolimnion'}
        }
        
        self.data={}
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    def set_time_series(self,tnum,tdate,param_process,Qout):
        """Function set_time_series

        Initiate the variables as time series.

        Inputs:
            tnum (numpy array of float): time values as number of seconds since 01.01.1970
            tdate (numpy array of datetime): datetime values corresponding to tnum
            param_process (dict): parameters related to the specific compound (from csv file)
            Qout (numpy array): time series of outflow discharge [m3/s]
            
        Outputs:
            None (addition of data in COMP object only)
        """
        
        
        self.data["time"]=tnum
        self.data["tdate"]=tdate
        self.data["conc_epi"]=np.full((2,len(tnum)),np.nan)
        self.data["conc_epi"][:,0]=param_process["Initial_concentration"]
        self.data["conc_hypo"]=np.copy(self.data["conc_epi"])
        self.data["Fact"]=np.zeros((2,len(tnum)))
        self.data["Fz"]=np.full((2,len(tnum)),np.nan)
        self.data["Fout"]=np.full((2,len(tnum)),np.nan)
        self.data["Qout"]=Qout
  
        
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    def compute_hepi_constant(self,datesval,tdate_start,tdate_end,hepi,zmax,filter_stratif=np.array([])):
        # To fill
        return 0
            
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    def compute_hepi_varying(self,tnum_TP,tnum_T,depth_T,tempval,zmax,mingrad=0.05,windowsize=4):
        # To fill
        return 0
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            
    def predict_concentration(self,param_lake,param_act,param_process,units_process):
        """Function predict_concentration
    
        Predicts concentrations in the epilimnion and hypolimnion based on a two-box model. 
    
        Inputs:
        
        Outputs:
            
        """  
        
        # Conversion to same concentration units (g.m-3):
        if units_process["Inflow_concentration"]=="g.m-3" or units_process["Inflow_concentration"]=="mg.l-1":
            fact_conc=1
        elif units_process["Inflow_concentration"]=="kg.m-3" or units_process["Inflow_concentration"]=="g.l-1":
            fact_conc=1000
        elif units_process["Inflow_concentration"]=="t.m-3" or units_process["Inflow_concentration"]=="kg.l-1":
            fact_conc=1e6   
        elif units_process["Inflow_concentration"]=="mg.m-3" or units_process["Inflow_concentration"]=="ug.l-1":
            fact_conc=1e-3
        elif units_process["Inflow_concentration"]=="ug.m-3" or units_process["Inflow_concentration"]=="ng.l-1":
            fact_conc=1e-6
        elif units_process["Inflow_concentration"]=="ng.m-3":
            fact_conc=1e-9
        else: 
            raise Exception("Unknown concentration units")
            
        # 1-Input from activity:
        if self.general_attributes["name"]=="UV_filter":
            self.input_UV_filter(param_act,param_process)
        else:
            raise Exception("Only UV filter can be simulated for now")   
            
        # 2-Input from inflows:
        self.data["Fin"]=np.full(self.data["time"].shape,
                    self.data["Qout"]*param_process["Inflow_concentration"]*fact_conc) # [g.s-1]
        
        # 3-Reaction (to add)
        self.data["Repi"]=np.zeros(self.data["time"].shape)
        self.data["Rhypo"]=np.zeros(self.data["time"].shape)
        
        # 4-Other fluxes and mass budget (iterative):
        h_therm=np.full(self.data["time"].shape,param_lake["Thermocline_depth"]) # Assumed constant for now
        h_epi,h_hypo,p=self.mass_budget(param_lake["Volume"],param_lake["Surface_area"],h_therm,param_lake["Maximum_depth"],self.data["Qout"],param_lake["Kz"])
        
        return h_epi,h_hypo,p

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            
    def input_UV_filter(self,param_act,param_process):
        """Function input_UV_filter
    
        Creates the time series of the compound input from the tourism activity.
    
        Inputs:
        
        Outputs:
            
        """  
        
        
        tnum_start_act=param_act["Start_swimming_time"].replace(tzinfo=timezone.utc).timestamp()
        tnum_end_act=param_act["End_swimming_time"].replace(tzinfo=timezone.utc).timestamp()
        
        tnum_start_season=param_act["Start_swimming_day"].replace(tzinfo=timezone.utc).timestamp()
        tnum_end_season=param_act["End_swimming_day"].replace(tzinfo=timezone.utc).timestamp()
        
        avg_input=param_act["Number_swimmers"]*param_process["Amount_per_person"]*\
        param_process["Percentage_UV_filter"]/100*param_process["Percentage_release"]/100*param_act["Swimming_duration"]\
         /(tnum_end_act-tnum_start_act)   # average input during the activity [g/s]
         
        
        if isinstance(avg_input, np.ndarray) and len(avg_input)==2: # Range type
            avg_input=avg_input[:,None] # Transpose
        else:
            avg_input=np.full((2,),avg_input) # Repeat the value
            
        dt_sim=self.data["time"][1]-self.data["time"][0]
        for kt in range(len(self.data["time"])):
            if self.data["time"][kt]>=tnum_start_season and self.data["time"][kt]<tnum_end_season:
                if (tnum_end_act-tnum_start_act)>dt_sim: # The activity is longer than the simulation time step
                    nsec_day=self.data["time"][kt]-np.floor(self.data["time"][kt]/86400)*86400
                    if nsec_day>=(tnum_start_act-np.floor(tnum_start_act/86400)*86400) and nsec_day<(tnum_end_act-np.floor(tnum_end_act/86400)*86400): # During activity period  
                        self.data["Fact"][:,kt]=avg_input[:,0]
                else: # The activity is shorter than the simulation time step: use daily average
                    self.data["Fact"][:,kt]=avg_input[:,0]*(tnum_end_act-tnum_start_act)/86400 # Average over 24h [g/s]

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            
    def mass_budget(self,Vlake,A0,h_therm,hmax,Qout,Kz):
        """Function mass_budget
    
        Computes the changes in epilimnetic and hypolimnetic concentrations over time according to the mass budget.
    
        Inputs:
        
        Outputs:
            
        """  
        
        # Assume power-law depth-dependent volume to compute epilimnion volume such that (i) hepi+hypo=hmean=Vlake/A0 and (ii) Vintegrated from surface to hmax=Vlake.
        # General hypsometry: A(z)=A0​(1−z/zmax​​)^p, with p>2 for steep walls, p=2 for a cone and p<2 for shallow lakes (e.g., p=1 for a parabole). 
        # Corresponding volume: Vlake=A0*zmax/(1+p) --> calculate p such that this volume matches the actual lake volume
        
        p=(A0*hmax/Vlake)-1
        
        #Vepi=A0*hmax/3*(1-(1-h_therm/hmax)**3) # Volume of the truncated cone [m3]
        Vepi=A0*hmax*(1-(1-h_therm/hmax)**(p+1))/(p+1) # Volume of the truncated cone [m3]
        h_epi=Vepi/A0 # [m] Thermocline depth of a rectangular epilimnion 
        Vhypo=Vlake-Vepi # [m3]
        h_hypo=Vhypo/A0 # [m]
        Atherm=A0 # [m2]
        h_mean=Vlake/A0
        
        if np.sum(h_mean==(h_epi+h_hypo))<len(h_epi):
            print(h_epi+h_hypo)
            raise Exception("Incorrect calculation of mean lake depth for {} values".format(np.sum((h_epi+h_hypo)!=h_mean)))
        
        
        for kt in range(len(self.data["time"])-1):
            self.data["Fz"][:,kt]=2*Kz*Atherm*(self.data["conc_hypo"][:,kt]-self.data["conc_epi"][:,kt])/h_mean # [g/s]
            self.data["Fout"][:,kt]=self.data["conc_epi"][:,kt]*Qout[kt] # [g/s]
            
            # Epilimnion mass balance:
            F_balance_epi=self.data["Fin"][kt]+self.data["Fz"][:,kt]+self.data["Fact"][:,kt]-self.data["Fout"][:,kt] # [g/s]
            if h_epi[kt+1]>h_epi[kt]:
                Cstar=self.data["conc_hypo"][:,kt]
            else:
                Cstar=self.data["conc_epi"][:,kt]
            new_conc_epi=self.data["conc_epi"][:,kt]+F_balance_epi*(self.data["time"][kt+1]-self.data["time"][kt])/Vepi[kt]\
                -self.data["Repi"][kt]*(self.data["time"][kt+1]-self.data["time"][kt])\
                    +Cstar/h_epi[kt]*(h_epi[kt+1]-h_epi[kt]) # [g.m-3]
            new_conc_epi[new_conc_epi<0]=0
            self.data["conc_epi"][:,kt+1]=new_conc_epi
                    
            # Hypolimnion mass balance:
            F_balance_hypo=-self.data["Fz"][:,kt]
            
                
            new_conc_hypo=self.data["conc_hypo"][:,kt]+F_balance_hypo*(self.data["time"][kt+1]-self.data["time"][kt])/Vhypo[kt]\
                -self.data["Rhypo"][kt]*(self.data["time"][kt+1]-self.data["time"][kt])\
                    -Cstar/h_hypo[kt]*(h_epi[kt+1]-h_epi[kt]) # [g.m-3]
            new_conc_hypo[new_conc_hypo<0]=0
            self.data["conc_hypo"][:,kt+1]=new_conc_hypo
        
        return h_epi,h_hypo,p
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            
    def results_to_csv(self,file):
        """Function results_to_csv
    
        Export the compound data to a csv file.
    
        Inputs:
        
        Outputs:
            
        """ 
        
        df = pd.DataFrame()
        
        for key, values in self.variables.items():
            csv_varname=key+"["+values["unit"]+"]"
            
            if key=="tdate":
                df[csv_varname]=pd.to_datetime(self.data["tdate"],format="%Y-%m-%d %H:%M:%S")
            else: 
                if "_min" in key:
                    varname=key[:key.find("_min")]
                    df[csv_varname]=self.data[varname][0,:]
                elif "_max" in key:
                    varname=key[:key.find("_max")]
                    df[csv_varname]=self.data[varname][1,:]
                else:
                    df[csv_varname]=self.data[key]
                    
            if "conc" in csv_varname:
                df[csv_varname] = df[csv_varname].map(lambda x: '{0:.3e}'.format(x))
                    
        df.to_csv(file, sep=",",header=True,index=False)
        print("Results exported to {}".format(file))    
        
        
        