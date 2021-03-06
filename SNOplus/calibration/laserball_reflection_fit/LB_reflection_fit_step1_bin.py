
# coding: utf-8

# In[1]:

import ROOT, rat
import os, sys, pickle
import matplotlib
matplotlib.use('Agg') # Avoid display

import matplotlib.pyplot as plt
import numpy as np
import collections


import jp_mpl as jplot
import geo_studies, rat_misc

from scipy import optimize
reload(rat_misc)
from matplotlib.colors import LogNorm
import pickle

# My tools
import jp_analysis as analysis


# Airplane mode - just histogramming, not more
db = rat.RAT.DB.Get()
db.SetAirplaneModeStatus(True)
db.SetDefaultPlaneLockStatus(False)
print 'This is AIRPLANE MODE - be careful!'


socdir = '/sb/project/qbs-015-ac/jpyanez/data/SOC_files'
histdir = '/sb/project/qbs-015-ac/jpyanez/data/SOC_histograms_python'


npmts = 9728

# First, bin the data

bin_width   = 0.5 # in ns
time_window = 100. # 100 ns is enough for all runs (at least central ones)
start_time  = -10
residual_axis = np.arange(start_time, time_window + start_time+ bin_width/2., bin_width)
hist_size  = [residual_axis.size-1, npmts]


print residual_axis.size, hist_size
test_pmts = [304, 1002, 4035, 7932]


# In[6]:

def makeHistogram(local_time_array):
        
    #### Bin the TOA, but *first* subtract the time delay of maximum
    
    # Large ToA binning
    
    # Bin the ToA
    ybins  = np.linspace(local_time_array.mean()-100., 
                         local_time_array.mean()+200, 
                         (300/(bin_width/4.))+1)
    
    # Y-centers
    ycenters = (ybins[1:] + ybins[:-1])/2.
    
    # First histogram
    n, x = np.histogram(local_time_array, ybins)
    
    # Find the highest peak time
    prompt_peak_index = n.argmax()
    prompt_peak_time  = ycenters[prompt_peak_index] 
  
    # Now re-do the histogram minus this time
    
    h, x = np.histogram( local_time_array - prompt_peak_time, residual_axis )
    
    return h, prompt_peak_time


# In[7]:

infile_list = os.listdir(socdir)
#for i, x in enumerate(infile_list): print i, x 


# ## Running over all the files

# In[ ]:

# Complicated because of subruns

while len(infile_list) > 0:
    one_file = infile_list[0]
    data = None
    
    
    # Need to establish if the file has subfiles
    run_number = one_file[:14]
    new_list   = []
    for one_name in infile_list:
        if run_number in one_name:
            new_list.append(one_name)
    for one_name in new_list:
        infile_list.remove(one_name)
        
    
    outfile_name =   os.path.join(histdir, run_number + '.pckl')
    if os.path.isfile(outfile_name):
        print run_number, 'Run already processed. Skipping it ...'

        # Open it and add the residual axis
        #data = pickle.load(open(outfile_name))
        #data['residual_axis'] = residual_axis
        #data['bin_width'] = bin_width
        #pickle.dump(data, open(outfile_name, 'w'))

        continue

    # Going over the new list that includes subruns, open all SOC files
    reader_list = []
    soc_list = []
    zero_counter = 0
    skip_counter = 0
    soc_pmts = np.array([0])
    already_processed = False
    for i, this_name in enumerate(new_list):
        infile_name = os.path.join(socdir, this_name)        
        if i == 0:
            print '\nRun', this_name, ' - Subruns', len(new_list)
            outfile_name = os.path.join(histdir, run_number + '.pckl')
            print outfile_name
                
        
        # Need to open all the readers
        reader_list.append(rat.socreader(infile_name))
        this_soc, this_run = reader_list[-1].next()
        soc_list.append(this_soc)
        soc_pmts = np.concatenate((soc_pmts, np.array(this_soc.GetSOCPMTIDs())))

        if i == 0:
            # Get the basic info from the first file
            manip_pos = np.array(this_soc.calib.GetPos())
            fit_pos  = np.array(this_soc.GetFitResult(this_soc.GetFitNames()[0]).GetVertex(0).GetPosition())

            data = {'manip_pos':manip_pos,
                    'fit_pos':fit_pos,
                    'soc_pmts':soc_pmts,
                    'wavelength':this_soc.calib.GetMode(),
                    'residual_axis': residual_axis,
                    'bin_width':bin_width}
    
            time_residuals = np.zeros(hist_size)
        
            print 'Going over the ', len(soc_pmts), ' PMTs'
                              
    if already_processed:
        continue     

        
    # Now that SOC files have been opened, fill the histogram
    soc_pmts    = np.unique(soc_pmts)
    pmt_delays  = np.zeros(npmts)

        
    for ipmt, one_pmt in enumerate(soc_pmts):
        time_array = np.array([0])
        for isoc, soc in enumerate(soc_list):
            try:
                this_array = np.array(soc.GetSOCPMT(one_pmt).GetTimes())
                time_array = np.concatenate((time_array,
                                             this_array[this_array > 0]))
            except:
                print 'PMT problematic', one_pmt, new_list[isoc]
                continue
                
            if ipmt % 500 == 0:
                print ipmt
                
        # PMT ready. Make histogram and store time delay
        time_residuals[:, one_pmt], pmt_delays[one_pmt] =  makeHistogram(time_array)
        
    data['time_residuals'] = time_residuals
    data['pmt_delays']     = pmt_delays
    pickle.dump(data, open(outfile_name,'w'))
    
    for reader in reader_list:
        reader.close()
    
    # Verification for one file
    fig = plt.figure()
    colors = ['b','r','g','m']
    for iTest, pmtid in enumerate(test_pmts):
        jplot.unfilledBar(residual_axis + data['pmt_delays'][pmtid],
                          data['time_residuals'][:,pmtid], color = colors[iTest])
    plt.yscale('log')
    fig.savefig(outfile_name.rstrip('.pckl') + '.png')



