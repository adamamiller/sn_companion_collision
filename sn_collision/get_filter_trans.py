import pandas as pd
from sn_companion_collision.sn_collision.kasen import Filter

import os
dirname = os.path.dirname(__file__)

def get_uvw2_tc():
    filt_df = pd.read_csv(dirname + '/filter_data/UVOT/UVW2_synphot.txt', 
                          delim_whitespace=True, skiprows=[0], 
                          header=None, names = ['wavelength', 'transmission']
                          )
    return Filter("uvw2", filt_df.wavelength.values, filt_df.transmission.values)
    
def get_uvm2_tc():
    filt_df = pd.read_csv(dirname + '/filter_data/UVOT/UVM2_synphot.txt', 
                          delim_whitespace=True, skiprows=[0], 
                          header=None, names = ['wavelength', 'transmission']
                          )
    return Filter("uvm2", filt_df.wavelength.values, filt_df.transmission.values)

def get_uvw1_tc():
    filt_df = pd.read_csv(dirname + '/filter_data/UVOT/UVW1_synphot.txt', 
                          delim_whitespace=True, skiprows=[0], 
                          header=None, names = ['wavelength', 'transmission']
                          )
    return Filter("uvw1", filt_df.wavelength.values, filt_df.transmission.values)
    
def get_u_uvot_tc():
    filt_df = pd.read_csv(dirname + '/filter_data/UVOT/U_UVOT_synphot.txt', 
                          delim_whitespace=True, skiprows=[0], 
                          header=None, names = ['wavelength', 'transmission']
                          )
    return Filter("u_uvot", filt_df.wavelength.values, filt_df.transmission.values)

def get_b_uvot_tc():
    filt_df = pd.read_csv(dirname + '/filter_data/UVOT/B_UVOT_synphot.txt', 
                          delim_whitespace=True, skiprows=[0], 
                          header=None, names = ['wavelength', 'transmission']
                          )
    return Filter("b_uvot", filt_df.wavelength.values, filt_df.transmission.values)

def get_v_uvot_tc():
    filt_df = pd.read_csv(dirname + '/filter_data/UVOT/V_UVOT_synphot.txt', 
                          delim_whitespace=True, skiprows=[0], 
                          header=None, names = ['wavelength', 'transmission']
                          )
    return Filter("v_uvot", filt_df.wavelength.values, filt_df.transmission.values)
