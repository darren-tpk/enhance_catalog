# Scan data to make detections
# Data is from 20km radius around Augustine

# Import packages and functions that we need
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from eqcorrscan import Party
from eqcorrscan import Family
from eqcorrscan import Template
from eqcorrscan import Detection
from obspy import UTCDateTime
from phase_processing.ncsn2pha import ncsn2pha
from phase_processing.read_hypoddpha import read_hypoddpha

# Load one party's data
save_dir = '/home/ptan/project/output/'
party_file = save_dir + 'sample_party.pickle'
party_in = open(party_file,'rb')
party = pickle.load(party_in)
st = pickle.load(party_in)
print(party)
print(st)

# Load all party data and append
party_all = Party()
#st_all = Stream()
main_dir = '/home/ptan/output/'
save_dir =  main_dir + 'pickles/'
all_files = os.listdir(save_dir)
for i in range(len(all_files)):
#for i in range(5):
    file_string = all_files[i]
    if '_party.pickle' in file_string:
        party_file = save_dir + file_string
        party_in = open(party_file,'rb')
        party = pickle.load(party_in)
        #st = pickle.load(party_in)
        party_all = party_all + party
        #st_all = st_all + st

# Save combined party and st data in pickle
party_all_file = save_dir + 'party_all.pickle'
party_all_out = open(party_all_file,'wb')
pickle.dump(party_all,party_all_out)
#pickle.dump(st_all,party_all_out)
party_all_out.close()

### Clean repeating events
# Craft an inclusion list, list of indices i&j
# and arrays of detection times & values
include_list = []
all_i = []
all_j = []
all_detection_time = []
all_detection_val = []
all_detection_threshold = []
for i in range(len(party_all.families)):
    include_list.append(np.ones(len(party_all[i])))
    family = party_all[i]
    for j in range(len(family)):
        all_i.append(i)
        all_j.append(j)
        all_detection_time.append(family[j].detect_time)
        all_detection_val.append(family[j].detect_val)
        all_detection_threshold.append(family[j].threshold)
all_detection_time = np.array(all_detection_time)
all_detection_val = np.array(all_detection_val)
all_detection_threshold = np.array(all_detection_threshold)
# Now loop through all detections again, making comparisons
for i in range(len(party_all.families)):
    family = party_all[i]
    for j in range(len(family)):
        # If the event is already excluded, continue
        if include_list[i][j] == 0:
            continue
        # Calculate time difference for all other detections
        time_difference = abs(all_detection_time-family[j].detect_time)
        matching_index = np.where(time_difference<30)[0]
        # If the only matching event is itself, continue
        if len(matching_index) == 1:
            continue
        # Otherwise, find highest detection value, and set exclude other events
        else:
            matching_detection_diff = all_detection_val[matching_index] - all_detection_threshold[matching_index]
            max_detection_index = matching_index[np.argmax(abs(matching_detection_diff))]
            exclude_index = matching_index[matching_index!=max_detection_index]
            for k in exclude_index:
                include_list[all_i[k]][all_j[k]] = 0
# Lastly, populate new party with only unrepeated detections
party_clean = Party()
for i in range(len(party_all.families)):
    family_clean = Family(template=party_all[i].template)
    for j in range(len(party_all[i])):
        if include_list[i][j] == 1:
            family_clean = family_clean + party_all[i][j]
    party_clean = party_clean + family_clean

# Save cleaned party in pickle
party_clean_file = save_dir + 'party_clean.pickle'
party_clean_out = open(party_clean_file,'wb')
pickle.dump(party_clean,party_clean_out)
party_clean_out.close()

### create detection-threshold histogram
# Initialize lists for detection and threshold values
detection_list = []
threshold_list = []
for i in range(len(party_clean.families)):
    family = party_clean[i]
    for j in range(len(family)):
        # append detection & threshold values for every detection in the current family
        detection_value = abs(family[j].detect_val)
        detection_list.append(detection_value)
        threshold_value = family[j].threshold
        threshold_list.append(threshold_value)
# plot distribution of detections as histogram and save
detection_array = np.array(detection_list)-np.array(threshold_list)
detection_floor = np.floor(min(detection_array))
detection_ceil = np.ceil(max(detection_array))
threshold_array = np.unique(threshold_list)
fig,ax = plt.subplots()
ax.grid(True)
ax.hist(detection_array, bins=np.arange(detection_floor,detection_ceil,0.1),color='teal',edgecolor='black')
ax.set_xlim([detection_floor,max(detection_array)])
ax.set_xlabel('Detection-Threshold Gap');
ax.set_ylabel('Frequency');
ax.set_title('Histogram of Detection Values');
outfile = main_dir + 'figures/det_hist.png'
plt.savefig(outfile)

### Create plot cumulative detections in time (inbuilt for EQcorrscan)
party_clean.plot(save='True',savefile=main_dir + 'figures/cumulative_det.png')

### Create events per day plot (FOR DETECTIONS)
# Get all detection times from cleaned party object
all_detection_time = []
for family in party_clean:
    for detection in family:
        all_detection_time.append(detection.detect_time)
all_detection_time = np.array(all_detection_time)
# Bin by day
first_time = min(all_detection_time)
start_day = UTCDateTime(first_time.year,first_time.month,first_time.day)
num_day = np.ceil((max(all_detection_time) - min(all_detection_time))/86400)
bin_counts_det = []
for i in range(int(num_day)):
    bin_start = start_day + i*86400
    bin_end = start_day + (i+1)*86400
    bin_count = sum((all_detection_time>bin_start) & (all_detection_time <= bin_end))
    bin_counts_det.append(bin_count)
# Plot histogram
fig,ax = plt.subplots()
ax.grid(True)
ax.hist(range(int(num_day)),range(int(num_day+1)),weights=bin_counts_det,color='royalblue',edgecolor='black')
ax.set_xlim([0,num_day])
#ax.set_ylim([0,170])
ax.set_xlabel('Days since ' + str(start_day)[:10]);
ax.set_ylabel('Daily frequency');
ax.set_title('Plot of DETECTIONS per day, from ' + str(start_day)[:10])
outfile = main_dir + 'figures/det_perday.png'
plt.savefig(outfile)

### Create events per day plot (FOR RAW CATALOG)
# Load raw catalog and extract all event times
hypoi_file = '/home/ptan/anaconda3/envs/eqcorrscan/attempt_eqcorrscan/avo_data/redoubt_20080101_20100101_hypoi.txt'
hypoddpha_file = '/home/ptan/anaconda3/envs/eqcorrscan/attempt_eqcorrscan/avo_data/redoubt_20080101_20100101_hypoddpha.txt'
ncsn2pha(hypoi_file, hypoddpha_file)
catalog_raw = read_hypoddpha(hypoi_file, hypoddpha_file)
all_event_time = []
fetch_time_start = UTCDateTime(2009,2,15,21,0,0)
fetch_time_end = UTCDateTime(2009,3,15,21,0,0)
for event in catalog_raw:
    event_time = event.origins[0].time
    if event_time >= fetch_time_start and event_time <= fetch_time_end:
        all_event_time.append(event.origins[0].time)
all_event_time = np.array(all_event_time)
# Bin by day
bin_counts_cat = []
for i in range(int(num_day)):
    bin_start = start_day + i*86400
    bin_end = start_day + (i+1)*86400
    bin_count = sum((all_event_time>bin_start) & (all_event_time <= bin_end))
    bin_counts_cat.append(bin_count)
# Plot histogram
fig,ax = plt.subplots()
ax.grid(True)
ax.hist(range(int(num_day)),range(int(num_day+1)),weights=bin_counts_cat,color='indianred',edgecolor='black')
ax.set_xlim([0,num_day])
#ax.set_ylim([0,170])
ax.set_xlabel('Days since ' + str(start_day)[:10]);
ax.set_ylabel('Daily frequency');
ax.set_title('Plot of EVENTS per day, from ' + str(start_day)[:10]);
outfile = main_dir + 'figures/event_perday.png'
plt.savefig(outfile)

### Plot stacked event vs detections per day histogram for comparison
fig,ax = plt.subplots()
ax.grid(True)
ax.hist(range(int(num_day)),range(int(num_day+1)),weights=bin_counts_det,color='royalblue',edgecolor='black',label = 'EQcorrscan Detections')
ax.hist(range(int(num_day)),range(int(num_day+1)),weights=bin_counts_cat,color='indianred',edgecolor='black',label = 'AVO Catalog',rwidth=0.8)
ax.legend()
ax.set_xlim([0,num_day])
#ax.set_ylim([0,170])
ax.set_xlabel('Days since ' + str(start_day)[:10]);
ax.set_ylabel('Daily frequency');
ax.set_title('Plot of EVENTS vs DETECTIONS per day, from ' + str(start_day)[:10]);
outfile = main_dir + 'figures/comparison_perday.png'
plt.savefig(outfile)

# Create family using raw catalog events
aligned_event_time = all_event_time[all_event_time>start_day]
rawcat_family = Family(template=Template(name='original_catalog'))
for detect_time in aligned_event_time:
    detection = Detection('original_catalog',detect_time,no_chans=99,detect_val=99,
            threshold=0,typeofdet='corr',threshold_type='MAD',threshold_input=20)
    rawcat_family = rawcat_family + detection
# Create clean family, combining all clean families
all_detection_time = []
for family in party_clean:
    for detection in family:
        all_detection_time.append(detection.detect_time)
cleandet_family = Family(template=Template(name='eqcorrscan_detections'))
for detect_time in all_detection_time:
    detection = Detection('eqcorrscan_detections',detect_time,no_chans=99,detect_val=99,
            threshold=0,typeofdet='corr',threshold_type='MAD',threshold_input=20)
    cleandet_family = cleandet_family + detection
# Combine both to craft comparison cumulative plot
party_combine = Party()
party_combine = party_combine + cleandet_family + rawcat_family
party_combine.plot(save='True',savefile=main_dir + 'figures/cumulative_comparison.png')

# ### Merge streams and create repicked catalog
# st_merged = st_all.merge(method=1)
# repicked_catalog = party.lag_calc(st_merged, pre_processed=False, shift_len=0.5, min_cc=0.4)
# # print UTCDatetime of repicked catalog events
# all_time = []
# for i in range(len(repicked_catalog)):
#     time_str = str(repicked_catalog[i].resource_id)[20:]
#     time = UTCDateTime(int(time_str[0:4]),int(time_str[4:6]),int(time_str[6:8]),int(time_str[9:11]),int(time_str[11:13]),int(time_str[13:15]),int(time_str[15:]))
#     all_time.append(time)
#     print('Event ' + str(i+1) + ':  ' + str(time))
#     for j in range(len(repicked_catalog[i].picks)):
#         print(repicked_catalog[i].picks[j].comments)
# all_time = list(np.sort(all_time))

