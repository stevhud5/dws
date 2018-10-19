# -*- coding: utf-8 -*-
"""
Template
Created on Tue Apr 24 13:59:42 2018
@author: sdh
"""
#
# Copyright 2008 - 2013 Brian R. D'Urso
#
# This file is part of Python Instrument Control System, also known as Pythics.
#
# Pythics is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Pythics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Pythics.  If not, see <http://www.gnu.org/licenses/>.
#
"""
template is the namespace from the outside, and we define methods and attributes under it. 
Q: We can define a class or not to manage the namespace inside; encapsulation etc? 
"""

#load libraries
import logging
import os.path
import numpy as np
import DWS_modules.DWS_DataTreatment_1 as DWS
import mpl
#import task
#import json

# private data shared among methods
class Private(object):
    pass
private = Private()

#methods
def initialize(log_text_box, paramload, **kwargs):
    # setup the logger and log display box
    # examples of how to use  (I need explanation):
#    print(logging.getLogger('test'))
    private.logger = logging.getLogger('log')
    #private.logger.setLevel(logging.DEBUG)
    private.logger.setLevel(logging.INFO)
    sh = logging.StreamHandler(log_text_box)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    sh.setFormatter(formatter)
    private.logger.addHandler(sh)
    private.logger.info('initialize complete')
    clear(**kwargs)

def run(folder, **kwargs):
    #
    #open and read next mra file. 
    pass

def process(fn, **kwargs)

def clear_all(data_plot, **kwargs):
    clear(data_plot)

def clear(data_plot, **kwargs):
    data_plot.clear()
    data_plot.set_plot_properties(
        x_label='t(s)',
        x_scale='log',
        x_limits=(1e-8,100),
        y_label='g2_1',
        dpi=120,
        aspect_ratio='auto')
        
#    data_plot.new_curve('data', line_color='blue', line_style='',
#                     marker_style='o', marker_color='blue', marker_width=3)
    data_plot.new_curve('mean', line_color='green', 
                     marker_style='+', marker_color='green', marker_width=3)
    
def array_to_string(f_num_array):
    text = ''
    for i in range(0,len(f_num_array)):
        text = text+'_'+str(f_num_array[i])
    return text
    
def build_mean_fn():
    #first draft will assume that all in the same folder, just different numbers at the end
    Na = len(private.filename_array)
    f_root = []
    f_ext = []
    for i in range(0,Na):
        _f_root,_f_ext = os.path.splitext(private.filename_array[i])
        f_root.append(_f_root)
        f_ext.append(_f_ext)
    f_root_common = os.path.commonprefix(f_root[0:2])
    Ncommon = len(f_root_common)
    while f_root_common[Ncommon-1]!='_':
        Ncommon = Ncommon-1
    f_root_common = f_root_common[:Ncommon-1]
    #compute f_num_array
    f_num_array = []
    for i in range (0,Na):
        temp_str=f_root[i]
        f_num_array.append(temp_str[Ncommon:])
    fn_out = f_root_common+'_mean'+array_to_string(f_num_array)+f_ext[0]
    return fn_out
    
#do this when file is selected
def add_to_array(filename_in, filename_out, **kwargs):  
    private.filename_array.append(filename_in.value) 
    filename_out.value = build_mean_fn()
    
def del_from_array(filename_in, filename_out, index, **kwargs):  
    private.filename_array.pop(int(index.value)) 
    filename_out.value = build_mean_fn()
    
#do this when trigger_B button is pressed
def plot_all(data_plot, **kwargs): #Bstatus is the indicator that displays the result. 
    maxCRtolerated=500
    clear(data_plot)
    
    ICF_array = []
    LagTime_array = []
    for j in range(0,len(private.filename_array)):
        fname=private.filename_array[j]
        MTDuration,EchoDuration,wavelength,T,L,labs,n,R,lstar,CRmean,g2,CRhistory=DWS.loadDWSdata(fname) 
        LagTime=g2[:,0]
        ICF=g2[:,1]
        LagTime_array.append(LagTime)
        ICF_array.append(ICF)
        private.logger.info(str(j)+', '+str(T))
        data_plot.new_curve('data'+str(j), line_color='blue', line_style='',
                     marker_style='o', marker_color='blue', marker_width=3)
        data_plot.set_data('data'+str(j),g2)
    
    if len(private.filename_array)>1:
        private.ICF_mean = np.mean(ICF_array, axis=0)
        ICF_sd = np.std(ICF_array, axis=0)
        ICF_grad = np.gradient(ICF_array, axis=0)
        data_plot.set_data('mean',np.column_stack([LagTime, private.ICF_mean]))
    else:
        private.ICF_mean = ICF
    pass
    
def write(filename_out, **kwargs):
    #get header to first file
    fname=private.filename_array[0]
    with open(fname, 'r') as f_in:
        header = ''
        MTDuration,EchoDuration,wavelength,T,L,labs,n,R,lstar,CRmean,g2,CRhistory=DWS.loadDWSdata(fname) 
        LagTime=g2[:,0]
        for index in range (21):
            header = header+f_in.readline()
        
    with open(filename_out.value, 'w') as f_out:
        f_out.write(header)
        for i in range(len(LagTime)):
            #f_out.write(LagTime[i]+'\t'+ICF_mean[i]+'\n')
            f_out.write('%8.6G\t%8.6G\n' % (LagTime[i],private.ICF_mean[i]))
            
        f_out.write('\n\nCount Rate  History (kHz) \nElapsed Time (s)\tCount Rate (kHz)\n\n')
        for i in range(len(CRhistory)):
            f_out.write('%10.6f\t%10.6f\n' % (CRhistory[i,0],CRhistory[i,1]))
            pass
        
        f_out.write('\nMeasurement notes:\n\n')
    
    
def read_config(config_filename_in, paramload, **kwargs):
    try:
        paramload.load_parameters(config_filename_in.value)
        private.logger.info('Should add a check to see if parameter file is correct format.')
    except:
        private.logger.exception('parameters did not save') 
