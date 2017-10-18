
# coding: utf-8

# In[1]:

import ROOT, rat
import os, sys, pickle
import matplotlib
matplotlib.use('Agg') # Avoid display

import matplotlib.pyplot as plt
import numpy as np
import collections
from scipy import interpolate


import jp_mpl as jplot


from scipy import optimize
from matplotlib.colors import LogNorm
import pickle

# My tools
import jp_analysis as analysis

# Define the run number, don't do all files at once
run_nr = sys.argv[1]
wavelength = 420.
manip_pos = np.array([0, -254., 25.])
max_events = 500000

# These are the conditions for selecting events
# For run 101428
nhit_range = [40, 1000]
median_range = [300, 310]
global_offset = np.mean(median_range)


# Correction of tof
# Optics DB water
dbwl = np.array([200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0, 390.0, 395.0, 400.0, 405.0, 410.0, 415.0, 420.0, 425.0, 430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0, 465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0, 500.0, 505.0, 510.0, 515.0, 520.0, 525.0, 530.0, 535.0, 540.0, 545.0, 550.0, 555.0, 560.0, 565.0, 570.0, 575.0, 580.0, 585.0, 590.0, 595.0, 600.0, 605.0, 610.0, 615.0, 620.0, 625.0, 630.0, 635.0, 640.0, 645.0, 650.0, 655.0, 660.0, 665.0, 670.0, 675.0, 680.0, 685.0, 690.0, 695.0, 700.0, 705.0, 710.0, 715.0, 720.0, 725.0, 730.0, 735.0, 740.0, 745.0, 750.0, 755.0, 760.0, 765.0, 770.0, 775.0, 780.0, 785.0, 790.0, 795.0, 800.0, ])
gvel = np.array([182.709, 185.104, 187.309, 189.341, 191.217, 192.951, 194.558, 196.047, 197.432, 198.72, 199.92, 201.04, 202.088, 203.068, 203.986, 204.849, 205.659, 206.422, 207.14, 207.818, 208.458, 209.063, 209.636, 210.179, 210.694, 211.183, 211.648, 212.09, 212.51, 212.911, 213.294, 213.659, 214.008, 214.342, 214.661, 214.967, 215.26, 215.541, 215.81, 216.069, 216.318, 216.557, 216.787, 217.008, 217.222, 217.427, 217.626, 217.817, 218.001, 218.18, 218.352, 218.518, 218.679, 218.835, 218.986, 219.132, 219.274, 219.411, 219.544, 219.673, 219.798, 219.92, 220.038, 220.153, 220.264, 220.372, 220.478, 220.58, 220.68, 220.777, 220.872, 220.964, 221.054, 221.141, 221.226, 221.31, 221.391, 221.47, 221.547, 221.622, 221.696, 221.768, 221.838, 221.907, 221.974, 222.039, 222.104, 222.166, 222.228, 222.288, 222.346, 222.404, 222.46, 222.515, 222.569, 222.622, 222.674, 222.725, 222.775, 222.824, 222.871, 222.918, 222.964, 223.01, 223.054, 223.097, 223.14, 223.182, 223.223, 223.264, 223.303, 223.342, 223.381, 223.418, 223.455, 223.491, 223.527, 223.562, 223.597, 223.631, 223.664,])
group_vel_fcn = interpolate.InterpolatedUnivariateSpline(dbwl, gvel)
light_c = group_vel_fcn(wavelength)
pmt_info = pickle.load(open('/home/jpyanez/snoplus/snoplus_python/pmt_positions.pckl'))





# Airplane mode - just histogramming, not more
db = rat.RAT.DB.Get()
db.SetAirplaneModeStatus(True)
db.SetDefaultPlaneLockStatus(False)
print 'This is AIRPLANE MODE - be careful!'


# Directory where the ROOT files (treated as normal data in ds format) are
datadir  = '/home/jpyanez/snoplus/data/laserball_runs/root'
histdir  = '/home/jpyanez/snoplus/data/laserball_runs/oca_mine'

npmts = 9728

# First, bin the data

bin_width   = 0.5 # in ns
time_window = 100. # 100 ns is enough for all runs (at least central ones)
start_time  = -10
residual_axis = np.arange(start_time, time_window + start_time+ bin_width/2., bin_width)
hist_size  = [residual_axis.size-1, npmts]


print residual_axis.size, hist_size
test_pmts = [304, 1002, 4035, 7932]




infile  = os.path.join(datadir, 'SNOP_0000'+run_nr+'.root')
outfile = os.path.join(histdir, 'SNOP_0000'+run_nr+'.pckl')
    
if not os.path.isfile(infile):
    print 'File doesnt exist!'
    print infile
    sys.exit()

data = {'residual_axis': residual_axis,
        'bin_width':bin_width}
    
time_residuals = np.zeros(hist_size)
        

reader = rat.dsreader(infile)

counter = 0
end_file = False
for ds, run in reader:
    for iEV in range(ds.GetEVCount()):
        event = ds.GetEV(iEV)
        nhits = int(event.GetNhits())
        if (nhits <  nhit_range[0]) or (nhits > nhit_range[1]):
            continue
        toa_series_manip = np.zeros(nhits)
        pmt_ids = np.zeros(nhits, dtype=int)
        pmts = event.GetCalPMTs()

        for iPMT in range(nhits):
            one_pmt = pmts.GetPMT(iPMT)
            pmt_ids[iPMT] = int(one_pmt.GetID())
            light_path_manip  =  np.linalg.norm(manip_pos - pmt_info['xyz'][one_pmt.GetID()])
            
            travel_time_manip =  light_path_manip/light_c
            #one_pmt.GetID() # In case I need the position
            this_time = one_pmt.GetTime()
            if this_time < 0:
                continue

            toa_series_manip[iPMT] = one_pmt.GetTime() - travel_time_manip

        toa_series_manip = toa_series_manip[toa_series_manip>0]
        event_median = np.median(toa_series_manip)

        if (event_median < median_range[0]) or (event_median > median_range[1]) :
            continue

        pmt_ids = pmt_ids[toa_series_manip > 0]

        for i, one_pmt in enumerate(pmt_ids):
            # Find where the hits go in the PMT histogram
            bin_index = np.digitize([toa_series_manip[i]]-global_offset, residual_axis)[0]
            if bin_index >= time_residuals.shape[0]:
                continue
            time_residuals[bin_index, one_pmt] += 1
            #print one_pmt, bin_index
        
        counter += 1

        if counter % 1000 == 0:
            print 'Event', counter

        if counter == max_events:
            end_file = True
            break
    if end_file:
        print 'Reached maximum number of events'
        break

        
data['time_residuals'] = time_residuals

pickle.dump(data, open(outfile, 'w'))

# Verification for one file
for iTest, pmtid in enumerate(test_pmts):
    

    if np.sum(data['time_residuals'][:,pmtid]) == 0:
        print 'No hits in ', pmtid
        continue
    fig = plt.figure()

    jplot.unfilledBar(residual_axis,
                      data['time_residuals'][:,pmtid])
    plt.yscale('log')
    plt.title('PMT ' + "%i" % pmtid)
    fig.savefig(outfile.rstrip('.pckl') + '_'+ "%i" % pmtid + '.png')

# Doing the full plot
fig = plt.figure(figsize=(12,7))

plt.pcolor(np.array(range(hist_size[1])), residual_axis, np.log10(data['time_residuals']))
plt.colorbar()
fig.savefig(outfile.rstrip('.pckl') + '_full.png')


