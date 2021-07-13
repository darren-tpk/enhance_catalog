#%% SENSITIVITY TEST

# This script loads in a party and rethresholds it over a range of thresholds to be explored.
# (This range of thresholds must be lower than the original party threshold.)
# The user can then choose to plot a detection value histogram, a cumulative event count plot,
# and template-detection waveform comparisons either separately or together.

# Import all dependencies
import numpy as np
import matplotlib.pyplot as plt
from eqcorrscan.utils.plotting import detection_multiplot
from obspy import UTCDateTime
from toolbox import get_detection, reader

#%% Define variables

# Define variables
main_dir = '/Users/darrentpk/Desktop/Github/enhance_catalog/'
party_dir = main_dir + 'output/scan_data/'
party_filename = 'party_redoubt.tgz'
plot_hist = False
plot_cumu = False
plot_wave = True
separate_wave = False
thres_min = 0.60
thres_max = 0.76
thres_vec = np.linspace(thres_min,thres_max,9)
num_days = 120

# Define a base time for x-axis
base_time = UTCDateTime(2009, 1, 1, 0, 0, 0)

# Define swarm start and swarm end hours for vertical dashed lines
swarm_times = [(UTCDateTime(2009, 2, 26, 6, 0, 0), UTCDateTime(2009, 2, 27, 13, 0, 0)),
               (UTCDateTime(2009, 3, 20, 12, 0, 0), UTCDateTime(2009, 3, 23, 6, 34, 0)),
               (UTCDateTime(2009, 3, 27, 0, 0, 0), UTCDateTime(2009, 3, 27, 8, 28, 0)),
               (UTCDateTime(2009, 3, 29, 7, 50, 0), UTCDateTime(2009, 3, 29, 9, 0, 0)),
               (UTCDateTime(2009, 4, 2, 19, 0, 0), UTCDateTime(2009, 4, 4, 13, 58, 0))]
               # (UTCDateTime(2009,5,2,21,0,0),UTCDateTime(2009,5,8,1,0,0)) May's swarm

# Define tick marks over span of party
tick_times = [UTCDateTime(2009, 1, 1), UTCDateTime(2009, 1, 15),
              UTCDateTime(2009, 2, 1), UTCDateTime(2009, 2, 15),
              UTCDateTime(2009, 3, 1), UTCDateTime(2009, 3, 15),
              UTCDateTime(2009, 4, 1), UTCDateTime(2009, 4, 15),
              UTCDateTime(2009, 5, 1)]

#%% Define functions

# NIL

#%% Load party and retrieve all detect_val, detect_time and no_chans

# Load party
party = reader(party_dir + party_filename)

# Extract detection value, number of channels and detection times
all_detect_val = []
all_no_chans = []
all_detect_time = []

# Loop through families in party
for i,family in enumerate(party):

    # Loop through detections in each family
    for detection in family:

        # Append detect value, number of channels and detection time
        all_detect_val.append(detection.detect_val)
        all_no_chans.append(detection.no_chans)
        all_detect_time.append(detection.detect_time)

# Convert lists to arrays for convenience
all_detect_val = np.abs(np.array(all_detect_val))
all_no_chans = np.array(all_no_chans)
all_detect_time = np.array(all_detect_time)

# If a detection value histogram plot is desired
if plot_hist:

    # Get the average channel correlation of all detections
    all_av_chan_corr = all_detect_val / all_no_chans

    # Plot histogram
    fig, ax = plt.subplots()
    ax.hist(all_av_chan_corr, bins=40, rwidth=0.8, color='brown')  # density=False would make counts
    ax.set_ylim([0,2000])
    ax.set_ylabel('Frequency')
    ax.set_xlabel('Average Channel Correlation')
    ax.set_title('Located Detections (N=%d)' % len(all_av_chan_corr))
    ax.grid()
    fig.show()

# If a cumulative event count plot for different thresholds is desired
if plot_cumu:

    # Start plot
    fig, ax = plt.subplots(figsize=(9,6))

    # Loop through thresholds
    for thres in thres_vec:

        # Calculate threshold sum
        all_thresholds = thres*all_no_chans

        # Check if detections pass this adjusted threshold, and print number of detections
        valid_bool = (all_detect_val >= all_thresholds)
        print('For thres =',thres,', number of detections =',sum(valid_bool))

        # Retrieve detection times to plot, and convert to hours
        plot_detect_time = all_detect_time[valid_bool]
        plot_detect_hours = (plot_detect_time - base_time) / 3600
        plot_detect_num = []
        plot_hours = []

        # Get cumulative numbers for plotting
        for i in range(num_days*24+1):
            detect_num = sum(plot_detect_hours <= i)
            plot_detect_num.append(detect_num)
            plot_hours.append(i)

        # Plot cumulative trend in a stepwise fashion
        ax.step(plot_hours,plot_detect_num,label=str(thres))

    # Add vertical spans indicating swarms
    for swarm_time in swarm_times:
        swarm_start_hour = (swarm_time[0] - base_time) / 3600
        swarm_end_hour = (swarm_time[1] - base_time) / 3600
        ax.axvspan(swarm_start_hour,swarm_end_hour,color='grey',alpha=0.5)

    # Tidy up plot
    ax.legend(title='Threshold',fontsize=12)
    ax.grid()
    ax.set_xlim([0,24*num_days])
    ax.set_xticks((np.array(tick_times) - base_time) / 3600)
    ax.set_xticklabels([tick_time.strftime('%b %d') for tick_time in tick_times],rotation=30,ha='right')
    ax.set_ylabel('Number of detections',fontsize=15)
    ax.set_xlabel('UTC Time',fontsize=15)
    ax.set_title('Cumulative detection plot for different thresholds',fontsize=15)
    fig.show()

# Conduct waveform review if desired
if plot_wave:

    # Craft matrices of appropriate size
    max_len = 0
    for family in party:
        if len(family) > max_len:
            max_len = len(family)
    detect_val_matrix = np.zeros((len(party.families),max_len))
    detect_val_matrix[:] = np.nan
    no_chans_matrix = np.zeros((len(party.families),max_len))
    detect_time_matrix = np.zeros((len(party.families),max_len))

    # Fill up matrix with detection values, no. of chans and detection time
    for i,family in enumerate(party):
        for j,detection in enumerate(family):
            detect_val_matrix[i][j] = abs(detection.detect_val)
            no_chans_matrix[i][j] = detection.no_chans
            detect_time_matrix[i][j] = detection.detect_time

    # Sequentially compare with the different thresholds in the threshold vector
    for thres in thres_vec:

        # Get a matrix of differences
        diff_matrix = abs(detect_val_matrix - (thres*no_chans_matrix))

        # Get the index of the detection that marginally passes the threshold, and pull the detection out of the party
        marginal_index = np.nanargmin(diff_matrix[np.nonzero(diff_matrix)])
        marginal_position = (int(np.floor(marginal_index/max_len)),marginal_index%max_len)
        marginal_detection = party[marginal_position[0]][marginal_position[1]]

        # Pull out the family from which the detection belongs
        family = party[marginal_position[0]]

        # If plotting the template and detection separately
        if separate_wave:
            family.template.st.plot(equal_scale=False, size=(800, 600))
            _ = get_detection(marginal_detection, data_dir='/home/data/redoubt/', length=10, plot=True)

        # If plotting the template and detection together
        else:
            detect_stream = get_detection(marginal_detection)
            detection_multiplot(detect_stream,family.template.st,[marginal_detection.detect_time])
