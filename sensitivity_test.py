# conduct sensitivity test for different threshold values (av_chan_corr)

# import packages that we need
import numpy as np
import matplotlib.pyplot as plt
from eqcorrscan.utils.plotting import detection_multiplot
from obspy import UTCDateTime
from toolbox import get_detection, reader

# define some inputs
party_input = party_all
#party_input = reader('/home/ptan/attempt_eqcorrscan/output/party_marchswarm.tgz')
plot_cumu = True
plot_wave = True
separate_wave = False
thres_min = 0.65
thres_max = 0.85
thres_vec = np.linspace(thres_min,thres_max,5)
num_days = 9

# retrieve all detect_val, detect_time and no_chans
all_detect_val = []
all_no_chans = []
all_detect_time = []
# loop through each family
for i,family in enumerate(party_input):
    # loop through each detection
    for detection in family:
        # append detect value, number of channels and detection time
        all_detect_val.append(detection.detect_val)
        all_no_chans.append(detection.no_chans)
        all_detect_time.append(detection.detect_time)
# convert all to arrays for convenience
all_detect_val = np.abs(np.array(all_detect_val))
all_no_chans = np.array(all_no_chans)
all_detect_time = np.array(all_detect_time)

# now create plot to see the influence of the threshold
if plot_cumu:
    # define a base time for x-axis
    base_time = UTCDateTime(2009, 3, 18, 0, 0, 0)
    # define swarm start and swarm end hours for vertical dashed lines
    swarm_start_hour = (UTCDateTime(2009,3,20,12,0) - base_time) / 3600
    swarm_end_hour = (UTCDateTime(2009,3,23,6,34) - base_time) / 3600
    # start plot
    fig, ax = plt.subplots(figsize=(9,6))
    # loop through thresholds
    for thres in thres_vec:
        # calculate threshold sum
        all_thresholds = thres*all_no_chans
        # now check if detections pass this adjusted threshold, and print number of detections
        valid_bool = (all_detect_val >= all_thresholds)
        print('for thres =',thres,', number of detections =',sum(valid_bool))
        # now retrieve detection times to plot, and convert to hours
        plot_detect_time = all_detect_time[valid_bool]
        plot_detect_hours = (plot_detect_time - base_time) / 3600
        plot_detect_num = []
        plot_hours = []
        # get cumulative numbers for plotting
        for i in range(num_days*24+1):
            detect_num = sum(plot_detect_hours <= i)
            plot_detect_num.append(detect_num)
            plot_hours.append(i)
        # plot cumulative trend in a stepwise fashion
        ax.step(plot_hours,plot_detect_num,label=str(thres))
    # add vertical lines indicating swarm
    ax.axvline(x=swarm_start_hour,color='black',ls='--')
    ax.axvline(x=swarm_end_hour,color='black',ls='--')
    # tidy up plot
    ax.legend(title='Threshold',fontsize=12)
    ax.grid()
    ax.set_xlim([0,24*num_days])
    ax.set_xticks(range(0,24*num_days+1,24))
    ax.set_xticklabels(['2009_03_' + str(int(day)) for day in range(18,28)],rotation=30,ha='right')
    ax.set_ylabel('Number of detections',fontsize=15)
    ax.set_xlabel('UTC Time',fontsize=15)
    ax.set_title('Cumulative detection plot for different thresholds',fontsize=15)
    fig.show()

# now conduct waveform review
if plot_wave:
    # we start by crafting a matrix of appropriate size
    max_len = 0
    for family in party_input:
        if len(family) > max_len:
            max_len = len(family)
    detect_val_matrix = np.zeros((len(party_input.families),max_len))
    detect_val_matrix[:] = np.nan
    no_chans_matrix = np.zeros((len(party_input.families),max_len))
    detect_time_matrix = np.zeros((len(party_input.families),max_len))
    # fill up matrix with detection values
    for i,family in enumerate(party_input):
        for j,detection in enumerate(family):
            detect_val_matrix[i][j] = abs(detection.detect_val)
            no_chans_matrix[i][j] = detection.no_chans
            detect_time_matrix[i][j] = detection.detect_time
    for thres in thres_vec:
        diff_matrix = abs(detect_val_matrix - (thres*no_chans_matrix))
        marginal_index = np.nanargmin(diff_matrix[np.nonzero(diff_matrix)])
        marginal_position = (int(np.floor(marginal_index/max_len)),marginal_index%max_len)
        marginal_detection = party_input[marginal_position[0]][marginal_position[1]]
        # plot template and detection
        family = party_input[marginal_position[0]]
        # plot separately
        if separate_wave:
            family.template.st.plot(equal_scale=False, size=(800, 600))
            _ = get_detection(marginal_detection, data_dir='/home/data/redoubt/', length=10, plot=True)
        # or plot together
        else:
            detect_stream = get_detection(marginal_detection, data_dir='/home/data/redoubt/', length=10, plot=False)
            detection_multiplot(detect_stream,family.template.st,[marginal_detection.detect_time])

