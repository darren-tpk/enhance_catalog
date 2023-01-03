# Import all dependencies
import numpy as np
import matplotlib.pyplot as plt
from eqcorrscan.utils.plotting import detection_multiplot
from obspy import UTCDateTime
from toolbox import get_detection, reader
from phase_processing.read_hypoddpha import read_hypoddpha

#%% Define variables

# Define variables
main_dir = '/home/ptan/enhance_catalog/'
data_dir = '/home/data/redoubt/'
party_dir = main_dir + 'output/redoubt2/scan_data/'
party_filename = 'party_special.tgz'
PEC_dir = main_dir + 'data/avo/'
hypoi_file = 'redoubt_20080501_20090901_hypoi.txt'
hypoddpha_file = 'redoubt_20080501_20090901_hypoddpha.txt'
thres_min = 16 # MAD
thres_max = 22 # MAD
thres_vec = np.linspace(thres_min,thres_max,4)
num_days = 14
base_time = UTCDateTime(2009, 2, 17)

PEC = read_hypoddpha(PEC_dir + hypoi_file,PEC_dir + hypoddpha_file,True)
for event in PEC:
    if event.origins[0].time < base_time or event.origins[0].time > (base_time+num_days*86400):
        PEC.events.remove(event)

party_raw = reader(party_dir + party_filename)
party = party_raw.copy()

tick_times = [UTCDateTime(2009, 2, 17),
              UTCDateTime(2009, 2, 19),
              UTCDateTime(2009, 2, 21),
              UTCDateTime(2009, 2, 23),
              UTCDateTime(2009, 2, 25),
              UTCDateTime(2009, 2, 27),
              UTCDateTime(2009, 3, 1),
              UTCDateTime(2009, 3, 3)]

# Extract detection value, number of channels and detection times
all_detect_time_lists = []
all_detect_day_lists = []

for thres in thres_vec:

    if thres != 9:
        party.rethreshold(thres,'MAD',abs_values=True)

    # initialize detect time list
    detect_time_list = []
    detect_day_list = []

    # Loop through families in party
    for i,family in enumerate(party):

        # Loop through detections in each family
        for detection in family:

            # Append detect time
            detect_time_list.append(detection.detect_time)
            detect_day_list.append((detection.detect_time - base_time) / 86400)

    all_detect_time_lists.append(detect_time_list)
    all_detect_day_lists.append(detect_day_list)

# Get AVO catalog information for plotting later
PEC_detect_times = [event.origins[0].time for event in PEC]
PEC_detect_days = [(t - base_time)/86400 for t in PEC_detect_times]

# Get cumulative numbers for plotting PEC
PEC_detect_num = []
PEC_days = []
num_points = 1000
dt = num_days/1000
for i in np.linspace(0,num_days,num_points+1):
    detect_num = sum(PEC_detect_days <= i)
    PEC_detect_num.append(detect_num)
    PEC_days.append(i)

# Start plot
fig, ax = plt.subplots(figsize=(12,9))

# Loop through thresholds
for i in range(len(thres_vec)):

    detect_day_list = all_detect_day_lists[i]

    # Retrieve detection times to plot, and convert to hours
    MAD_detect_num = []
    MAD_days = []
    for j in np.linspace(0, num_days, num_points + 1):
        detect_num = sum(detect_day_list <= j)
        MAD_detect_num.append(detect_num)
        MAD_days.append(j)

    # Plot cumulative trend in a stepwise fashion
    ax.step(MAD_days,MAD_detect_num,linestyle='-',linewidth=2,label='Threshold = %d*MAD' % thres_vec[i])

# plot redpy
filename = '/home/ptan/enhance_catalog/redpy/runs/redoubt2/catalog.txt'
all_redpy_times = []
with open(filename) as file:
    for line in file:
        all_redpy_times.append(UTCDateTime(line.rstrip().split(' ')[1]))
all_redpy_times = [t for t in all_redpy_times if t > base_time and t < (base_time+num_days*86400)]
redpy_day_list = [(t-base_time)/86400 for t in all_redpy_times]
# Retrieve detection times to plot, and convert to hours
redpy_detect_num = []
redpy_days = []
for j in np.linspace(0, num_days, num_points + 1):
    detect_num = sum(redpy_day_list <= j)
    redpy_detect_num.append(detect_num)
    redpy_days.append(j)
ax.step(redpy_days,redpy_detect_num,color='darkred',linestyle='-',linewidth=2,label='REDPy repeater catalog')

# Add PEC's cumulative trend in black
ax.step(PEC_days,PEC_detect_num,color='k',linestyle='-',linewidth=3,label='Original AVO catalog')

# Tidy up plot
ax.legend(fontsize=18)
ax.set_ylim(0,6000)
ax.grid()
ax.tick_params(axis='x',labelsize=18)
ax.tick_params(axis='y',labelsize=18)
ax.set_xlim([0,num_days])
ax.set_xticks((np.array(tick_times) - base_time) / 86400)
ax.set_xticklabels([tick_time.strftime('%m/%d') for tick_time in tick_times],fontsize=18,rotation=30,ha='right')
ax.set_ylabel('Number of detections',fontsize=20)
ax.set_xlabel('UTC Datetime',fontsize=20)
ax.set_title('Cumulative detection plot for different MAD thresholds',fontsize=20)
fig.show()
fig.savefig('/home/ptan/enhance_catalog/fig_Sx_MAD_comparison.pdf')
