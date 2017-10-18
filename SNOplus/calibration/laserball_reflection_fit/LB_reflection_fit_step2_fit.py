
import os, sys, pickle

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import collections


import jp_mpl as jplot

from scipy import optimize
from matplotlib.colors import LogNorm
import pickle

# My tools
import jp_analysis as analysis
import reflection_fit_v2 as reflection_fit
reload(reflection_fit)
from scipy import interpolate
from copy import deepcopy


# ## This step uses the histograms and calculates the peak positions of the first peak when the fit is initialized. The second peak is calculated as part of the minimization.

###
### Global variables
###
c            = 0.299792458*1000 # mm/ns
# This is the wl function from the optics DB
dbwl = np.array([200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 355.0, 360.0, 365.0, 370.0, 375.0, 380.0, 385.0, 390.0, 395.0, 400.0, 405.0, 410.0, 415.0, 420.0, 425.0, 430.0, 435.0, 440.0, 445.0, 450.0, 455.0, 460.0, 465.0, 470.0, 475.0, 480.0, 485.0, 490.0, 495.0, 500.0, 505.0, 510.0, 515.0, 520.0, 525.0, 530.0, 535.0, 540.0, 545.0, 550.0, 555.0, 560.0, 565.0, 570.0, 575.0, 580.0, 585.0, 590.0, 595.0, 600.0, 605.0, 610.0, 615.0, 620.0, 625.0, 630.0, 635.0, 640.0, 645.0, 650.0, 655.0, 660.0, 665.0, 670.0, 675.0, 680.0, 685.0, 690.0, 695.0, 700.0, 705.0, 710.0, 715.0, 720.0, 725.0, 730.0, 735.0, 740.0, 745.0, 750.0, 755.0, 760.0, 765.0, 770.0, 775.0, 780.0, 785.0, 790.0, 795.0, 800.0, ])
gvel = np.array([182.709, 185.104, 187.309, 189.341, 191.217, 192.951, 194.558, 196.047, 197.432, 198.72, 199.92, 201.04, 202.088, 203.068, 203.986, 204.849, 205.659, 206.422, 207.14, 207.818, 208.458, 209.063, 209.636, 210.179, 210.694, 211.183, 211.648, 212.09, 212.51, 212.911, 213.294, 213.659, 214.008, 214.342, 214.661, 214.967, 215.26, 215.541, 215.81, 216.069, 216.318, 216.557, 216.787, 217.008, 217.222, 217.427, 217.626, 217.817, 218.001, 218.18, 218.352, 218.518, 218.679, 218.835, 218.986, 219.132, 219.274, 219.411, 219.544, 219.673, 219.798, 219.92, 220.038, 220.153, 220.264, 220.372, 220.478, 220.58, 220.68, 220.777, 220.872, 220.964, 221.054, 221.141, 221.226, 221.31, 221.391, 221.47, 221.547, 221.622, 221.696, 221.768, 221.838, 221.907, 221.974, 222.039, 222.104, 222.166, 222.228, 222.288, 222.346, 222.404, 222.46, 222.515, 222.569, 222.622, 222.674, 222.725, 222.775, 222.824, 222.871, 222.918, 222.964, 223.01, 223.054, 223.097, 223.14, 223.182, 223.223, 223.264, 223.303, 223.342, 223.381, 223.418, 223.455, 223.491, 223.527, 223.562, 223.597, 223.631, 223.664,])
group_vel_fcn = interpolate.InterpolatedUnivariateSpline(dbwl, gvel)


# Load PMT positions
pmt_info = pickle.load(open('/home/jpyanez/snoplus/snoplus_python/pmt_positions.pckl'))

mode = sys.argv[1] # Can be simple or advanced
if len(sys.argv) > 2:
    print 'Fixing XY position'
    fix_xy = True
else:
    fix_xy = False

min_pmts     = 4000
fits_desired = 10 # Number of successful fits that I would like to have
max_fits     = 30 # Maximum number of fits that will be attempted

test   = False
if test:
    print '\n****Running in test mode! Only one fit for a single run******'
if fix_xy:
    print '**** FIXING X-Y POSITIONS TO MANIP'

# Settings
#myinput = '/home/jpyanez/scratch/laserball/337/p300_n300_0'
histdir = '/sb/project/qbs-015-ac/jpyanez/data/SOC_histograms_python'


if mode == 'simple':
    print '\n\n**********Doing the SIMPLE fits***************'
    fitdir  = '/sb/project/qbs-015-ac/jpyanez/data/SOC_reflection_fits_simple'
else:
    print '\n\n**********Doing the ADVANCED fits***************'
    fitdir = '/sb/project/qbs-015-ac/jpyanez/data/SOC_reflection_fits_advanced'

if fix_xy:
    fitdir += '_fix_xy'


# In[12]:

usable_run_list =  [17375,
100556,
100558,
101427,
101428,
101432,
101433,
102518,
102529,
102552,
102554,
102570,
102572,
102574]


# In[13]:

original_infile_list = os.listdir(histdir)
socfiles = []
for one_file in original_infile_list:
    for usable_run in usable_run_list:
        if ("%i" % usable_run in one_file) and ('pckl' in one_file):
            socfiles.append(one_file)
            break
            
print 'List of good run files'
socfiles


# In[14]:

inward_looking = (pmt_info['type'] == 1)+(pmt_info['type']==7)
non_bottom = (pmt_info['xyz'][:,2]>-6000)#*(pmt_info['xyz'][:,2]<6000)


# In[ ]:

runs_with_issues = []
for one_file in socfiles:

    print '\n\n ********** Using file', one_file
    histfile = os.path.join(histdir, one_file)
    out_file = os.path.join(fitdir, one_file)
    
    if not test:
        if os.path.isfile(out_file):
            print 'Run already fit - skipping'
            continue
    
    print out_file
    
    data_file = open(histfile)
    data = pickle.load(data_file)
    
    
    peakfit = reflection_fit.FitLBpos( data = data,
                                       pmt_xyz = pmt_info['xyz'],
                                       pmtbool = non_bottom*inward_looking,
                                       fit_mode = mode,
                                       fix_xy = fix_xy,
                                       print_call = True)
    my_x0 = np.concatenate((data['manip_pos']+np.array([10, -10, 10]), 
                            [1.03*c/group_vel_fcn(data['wavelength'])]))
    if test:
        my_x0 = np.concatenate((data['manip_pos']+ np.array([10.]*3), [1.36]))
    print 'Difference wrt MANIP', data['manip_pos'] - my_x0[:3]
    
    bestfit_value = 9E9
    results = errors = None
    all_results = {}
    peakfit.second_peak_error = []
    
    keep_fitting = True
    good_fit_counter = 0
    fits_attempted   = 0

    if fix_xy:
        my_x0 = my_x0[2:]
        bounds = ((-6500,6500), (0.5, 2.))
        peakfit.x = data['manip_pos'][0]
        peakfit.y = data['manip_pos'][1]
    else:
        bounds = ((-6500, 6500),
                  (-6500, 6500),
                  (-6500, 6500),
                  (0.5, 2.))
        print bounds

    while keep_fitting:

        print '\n****************Iteration ', good_fit_counter, '/', fits_attempted, '********************'
        wrapfcn = lambda p: peakfit(*p)

        peakr = optimize.minimize(wrapfcn,
                                  x0 = my_x0,
                                  method = 'SLSQP',
                                  bounds=bounds,
                                  options={'ftol':1E-7, 'maxiter':1000})            




        if peakr['success'] and np.sum(peakfit.thisbool)> min_pmts:
            good_fit_counter += 1
            all_results[good_fit_counter] = {'fit':deepcopy(peakr),
                                 'npmts':np.sum(peakfit.thisbool)}
            if fix_xy:
                my_x0 = np.concatenate(( [peakr.x[0]+np.random.normal(0., 5.)], 
                                         [peakr.x[1] + np.random.normal(0., 0.015)]))
            else:
                my_x0 = np.concatenate(( peakr.x[:3]+np.random.normal(0., 5., size=3), 
                                         [peakr.x[3] + np.random.normal(0., 0.015)]))
            
            if peakr['fun'] < bestfit_value:
                bestfit_value = peakr['fun']
                results = deepcopy(peakr)
        else:
            if fix_xy:
                my_x0 = np.concatenate(([data['manip_pos'][2]+np.random.normal(0., 40.)], 
                                        [c/group_vel_fcn(data['wavelength']) + np.random.normal(0., 0.04)]))
            else:
                my_x0 = np.concatenate((data['manip_pos']+np.random.normal(0., 40., size=3), 
                                        [c/group_vel_fcn(data['wavelength']) + np.random.normal(0., 0.04)]))                
        if test:
            print 'Final result'
            print peakr.x
            sys.exit()

        if good_fit_counter>= fits_desired:
            break


        # Making sure this gets reinitialized
        peakfit.second_peak_error = []
        fits_attempted += 1

        if fits_attempted >= max_fits:
            runs_with_issues.append(one_file)
            break


    data_file.close()

    if fits_attempted >= max_fits:
            print 'Reached the maximum number of fits allowed! not saving the result'
    else:
        pickle.dump(all_results, open(out_file, 'w'))

print '\n**** Files that didnt work ****'
print runs_with_issues






