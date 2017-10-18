
import ROOT, rat
import os, sys, pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import collections

import jp_mpl as jplot

from scipy import interpolate

import pickle



#### This is the part that changes depending on the run, although it does't seem to make much of a difference

run_nr = sys.argv[1]

if run_nr == '102529':
    wavelength = 420.
    manip_pos = np.array([0,-254.,25.])
    fit_pos = np.array([-10., -268., -3.])
elif run_nr == '100556':
    wavelength = 420.
    manip_pos = np.array([0, -254., 25.])
    fit_pos = np.array([3.83, -207.08, 24.79])
elif run_nr == '101161':
    wavelength = 420.
    manip_pos = np.array([0, -254., 25.])
    fit_pos = np.array([-10., -268., -3.]) # This is the wrong position
elif run_nr == '101428':
    wavelength = 420.
    manip_pos = np.array([0, -254., 25.])
    fit_pos = np.array([-10., -268., -3.]) # This is the wrong position
else:
    print 'Select a run that has been defined!'
    sys.exit()


max_events = 500000






# In[2]:

# Airplane mode
db = rat.RAT.DB.Get()
db.SetAirplaneModeStatus(True)
db.SetDefaultPlaneLockStatus(False)
print 'This is AIRPLANE MODE - be careful!'


# In[3]:

# Optics DB water
dbwl = np.array([200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0, 390.0, 395.0, 400.0, 405.0, 410.0, 415.0, 420.0, 425.0, 430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0, 465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0, 500.0, 505.0, 510.0, 515.0, 520.0, 525.0, 530.0, 535.0, 540.0, 545.0, 550.0, 555.0, 560.0, 565.0, 570.0, 575.0, 580.0, 585.0, 590.0, 595.0, 600.0, 605.0, 610.0, 615.0, 620.0, 625.0, 630.0, 635.0, 640.0, 645.0, 650.0, 655.0, 660.0, 665.0, 670.0, 675.0, 680.0, 685.0, 690.0, 695.0, 700.0, 705.0, 710.0, 715.0, 720.0, 725.0, 730.0, 735.0, 740.0, 745.0, 750.0, 755.0, 760.0, 765.0, 770.0, 775.0, 780.0, 785.0, 790.0, 795.0, 800.0, ])
gvel = np.array([182.709, 185.104, 187.309, 189.341, 191.217, 192.951, 194.558, 196.047, 197.432, 198.72, 199.92, 201.04, 202.088, 203.068, 203.986, 204.849, 205.659, 206.422, 207.14, 207.818, 208.458, 209.063, 209.636, 210.179, 210.694, 211.183, 211.648, 212.09, 212.51, 212.911, 213.294, 213.659, 214.008, 214.342, 214.661, 214.967, 215.26, 215.541, 215.81, 216.069, 216.318, 216.557, 216.787, 217.008, 217.222, 217.427, 217.626, 217.817, 218.001, 218.18, 218.352, 218.518, 218.679, 218.835, 218.986, 219.132, 219.274, 219.411, 219.544, 219.673, 219.798, 219.92, 220.038, 220.153, 220.264, 220.372, 220.478, 220.58, 220.68, 220.777, 220.872, 220.964, 221.054, 221.141, 221.226, 221.31, 221.391, 221.47, 221.547, 221.622, 221.696, 221.768, 221.838, 221.907, 221.974, 222.039, 222.104, 222.166, 222.228, 222.288, 222.346, 222.404, 222.46, 222.515, 222.569, 222.622, 222.674, 222.725, 222.775, 222.824, 222.871, 222.918, 222.964, 223.01, 223.054, 223.097, 223.14, 223.182, 223.223, 223.264, 223.303, 223.342, 223.381, 223.418, 223.455, 223.491, 223.527, 223.562, 223.597, 223.631, 223.664,])
group_vel_fcn = interpolate.InterpolatedUnivariateSpline(dbwl, gvel)
light_c = group_vel_fcn(wavelength)


# In[9]:

pmt_info = pickle.load(open('/home/jpyanez/snoplus/snoplus_python/pmt_positions.pckl'))


# In[84]:

outdir = '/sb/project/qbs-015-aa/jpyanez/data/laserball_runs/trigger_verification_figures'


# In[73]:




# In[74]:




# In[75]:

infile = '/home/jpyanez/snoplus/data/laserball_runs/root/SNOP_0000' + run_nr + '.root'


# In[76]:

median_tolerance = 1000.


# In[82]:

reader = rat.dsreader(infile)

all_median_manip = np.zeros(max_events)
all_median_fit = np.zeros(max_events)
event_time = np.zeros(max_events)

all_qhs = np.zeros(max_events)
all_qhl = np.zeros(max_events)
all_nhit = np.zeros(max_events)

counter = 0
exit_loop = False
for ds, run in reader:
    for iEV in range(ds.GetEVCount()):
        event = ds.GetEV(iEV)
        toa_series_manip = np.zeros(event.GetNhits())
        toa_series_fit = np.zeros(event.GetNhits())
        pmts = event.GetCalPMTs()
        tot_qhl = 0
        tot_qhs = 0
        all_nhit[counter] = pmts.GetCount()
        event_time[counter] = event.GetUniversalTime().GetNanoSeconds()
        for iPMT in range(int(all_nhit[counter])):
            one_pmt = pmts.GetPMT(iPMT)

            light_path_manip  =  np.linalg.norm(manip_pos - pmt_info['xyz'][one_pmt.GetID()])
            light_path_fit  =  np.linalg.norm(fit_pos - pmt_info['xyz'][one_pmt.GetID()])
            travel_time_manip =  light_path_manip/light_c
            travel_time_fit =  light_path_fit/light_c
            #one_pmt.GetID() # In case I need the position
            this_time = one_pmt.GetTime()
            if this_time < 0:
                continue

            toa_series_manip[iPMT] = one_pmt.GetTime() - travel_time_manip
            toa_series_fit[iPMT] = one_pmt.GetTime() - travel_time_fit
            all_qhl[counter] += one_pmt.GetQHL()
            all_qhs[counter] += one_pmt.GetQHS()

        # First easy check - do the median
        all_median_manip[counter] = np.median(toa_series_manip[toa_series_manip>0])
        all_median_fit[counter] = np.median(toa_series_fit[toa_series_fit>0])

        if np.abs(all_median_manip[counter]-310) > median_tolerance:
            xaxis = np.arange(100, 325, 0.25)
            counts, x = np.histogram(toa_series_manip, xaxis)
            jplot.unfilledBar(xaxis, counts)
            plt.show()
            raw_input()
            
        if counter %1000 == 0:
            print counter
        
        
        
        counter += 1
        if counter == max_events:
            print 'Max number of events reached. Exiting!'
            exit_loop = True
        
        if exit_loop:
            break
    if exit_loop:
        break
        
all_median_manip = all_median_manip[:counter]
all_median_fit = all_median_fit[:counter]

all_qhl = all_qhl[:counter]
all_qhs = all_qhs[:counter]
all_nhit = all_nhit[:counter]
event_time = event_time[:counter]
reader.close()


data = {'all_median_manip':all_median_manip,
        'all_median_fit':all_median_fit,
        'all_qhl':all_qhl,
        'all_qhs':all_qhs,
        'all_nhit':all_nhit,
        'event_time':event_time}

import pickle
pickle.dump(data, open(os.path.join(outdir, run_nr + '_data.pckl'), 'w'))


