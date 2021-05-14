# repository for all functions

# function to remove artificial boxcar-like signals from stream
def remove_boxcars(st,tolerance):
    # import packages
    import numpy as np
    from obspy import Stream
    # initialize output stream
    st_out = Stream()
    # loop through traces in stream
    for tr in st:
        # detrend data and find spikes above tolerance (median(abs(data)) multiplier)
        data = tr.data
        tr_detrended = tr.copy().detrend()
        data_detrended = tr_detrended.data
        spike_indices = np.where(data_detrended > tolerance * (np.median(abs(data_detrended))+1))[0]
        # extract isolated spikes
        isolated_spike_index = []
        for i, spike_index in enumerate(spike_indices):
            # if spike indices is 1x1
            if len(spike_indices) == 1:
                isolated_spike_index.append(i)
            # if isolated spike is the first element
            elif i == 0 and (spike_index != spike_indices[i+1]-1):
                isolated_spike_index.append(i)
            # if isolated spike is the last element
            elif i == len(spike_indices)-1 and (spike_index != spike_indices[i-1]+1):
                isolated_spike_index.append(i)
            # if the isolated spike is in the middle of the list
            elif (spike_index != spike_indices[i-1]+1) and (spike_index != spike_indices[i+1]-1):
                isolated_spike_index.append(i)
        # now determine where the isolated spikes and boxcars are by index
        isolated_spikes = [spike_indices[i] for i in isolated_spike_index]
        boxcar_indices = [spike_index for spike_index in spike_indices if spike_index not in spike_indices[isolated_spike_index]]
        # knowing our spike limits, we execute interpolation
        for isolated_spike in isolated_spikes:
            tr.data[isolated_spike] = (tr.data[isolated_spike-1] + tr.data[isolated_spike+1])/2
        # now extract boxcar limits
        boxcar_limits = []
        for i in range(len(boxcar_indices)):
            index = boxcar_indices[i]
            if index == boxcar_indices[0]:
                boxcar_limits.append(index)
            elif index == boxcar_indices[-1]:
                boxcar_limits.append(index)
            elif index != boxcar_indices[i - 1] + 1 or index != boxcar_indices[i + 1] - 1:
                boxcar_limits.append(index)
        # if len(boxcar_limits) % 2 != 0:
        #     print('WARNING: Data has a non-even number of spikes, skipping trace.')
        #     continue
        # knowing our spike limits, we get data limits for extraction
        boxcar_limits = np.array(boxcar_limits)
        extract_limits = [*sum(zip(boxcar_limits[0::2] - 1, boxcar_limits[1::2] + 1), ())]
        # add data's head and tail index to extraction limits
        extract_limits = np.append(extract_limits, len(data) - 1)
        extract_limits = np.insert(extract_limits, 0, 0)
        # using the sampling rate, we determine start time and end time of extraction
        tr_starttime = tr.stats.starttime
        dt = 1 / tr.stats.sampling_rate
        # loop through limit pairs
        for j in range(0, len(extract_limits), 2):
            extract_starttime = tr_starttime + (extract_limits[j] * dt)
            extract_endtime = tr_starttime + (extract_limits[j + 1] * dt)
            # extract trace and store in output stream
            tr_out = tr.slice(extract_starttime, extract_endtime)
            st_out = st_out + tr_out
    return st_out

# function to compare trace with data on IRIS
def compare_tr(tr):
    # import packages
    from obspy import Stream
    from waveform_collection import gather_waveforms
    # define net-sta-chan, starttime, endtime
    network = tr.stats.network
    station = tr.stats.station
    channel = tr.stats.channel
    starttime = tr.stats.starttime
    endtime = tr.stats.endtime
    # gather waveforms
    tr_client = gather_waveforms(source='IRIS', network=network, station=station,
                          location='*', channel=channel, starttime=starttime,
                          endtime=endtime)
    # combine and plot on same figure
    combine = Stream()
    combine = combine + tr
    combine = combine + tr_client
    combine.plot()

# function to clean repeating events in a party
def remove_repeats(party,time_interval):
    # import packages
    import numpy as np
    from eqcorrscan import Party, Family
    # initialize an inclusion list, list of all i&j, and arrays of detection times, values, and thresholds
    include_list = []
    all_i = []
    all_j = []
    all_detection_time = []
    all_detection_val = []
    all_detection_threshold = []
    # loop through each family in the party
    for i in range(len(party.families)):
        # populate include list with ones (we start by including all detections in each family)
        include_list.append(np.ones(len(party[i])))
        family = party[i]
        # loop through each detection in the family, storing key info
        for j in range(len(family)):
            all_i.append(i)
            all_j.append(j)
            all_detection_time.append(family[j].detect_time)
            all_detection_val.append(family[j].detect_val)
            all_detection_threshold.append(family[j].threshold)
    # convert to arrays
    all_detection_time = np.array(all_detection_time)
    all_detection_val = np.array(all_detection_val)
    all_detection_threshold = np.array(all_detection_threshold)
    # now loop through all families' detections again, making comparisons
    for i in range(len(party.families)):
        # extract family
        family = party[i]
        # loop through family's detections
        for j in range(len(family)):
            # if the event is already excluded, continue
            if include_list[i][j] == 0:
                continue
            # calculate time difference for all other detections
            time_difference = abs(all_detection_time-family[j].detect_time)
            matching_index = np.where(time_difference<time_interval)[0]
            # if the only matching event is itself, continue
            if len(matching_index) == 1:
                continue
            # otherwise, find highest detection value, and exclude other events
            else:
                matching_detection_diff = all_detection_val[matching_index] - all_detection_threshold[matching_index]
                max_detection_index = matching_index[np.argmax(abs(matching_detection_diff))]
                exclude_index = matching_index[matching_index!=max_detection_index]
                for k in exclude_index:
                    include_list[all_i[k]][all_j[k]] = 0
    # lastly, populate new party with only unrepeated detections
    party_clean = Party()
    for i in range(len(party.families)):
        family_clean = Family(template=party[i].template)
        for j in range(len(party[i])):
            if include_list[i][j] == 1:
                family_clean = family_clean + party[i][j]
        party_clean = party_clean + family_clean
    return party_clean

# function to get all stations within a radius of a volcano
def get_local_stations(volcano_name,radius):
    # import packages
    import numpy as np
    import pandas as pd
    import geopy.distance
    # get volcano lat and lon
    volcano_list = pd.read_csv('/home/ptan/attempt_eqcorrscan/avo_data/volcano_list.csv')
    station_list = pd.read_csv('/home/ptan/attempt_eqcorrscan/avo_data/station_list.csv')
    volcano_names = [name.lower() for name in volcano_list.volcano]
    volcano_index = volcano_names.index(volcano_name.lower())
    volcano_lat = volcano_list.latitude[volcano_index]
    volcano_lon = volcano_list.longitude[volcano_index]
    volcano_coord = (volcano_lat,volcano_lon)
    # now get distances with all stations
    station_distances = [];
    for i in range(len(station_list)):
        station_lat = station_list.latitude[i]
        station_lon = station_list.longitude[i]
        station_coord = (station_lat,station_lon)
        station_distance = geopy.distance.GeodesicDistance(volcano_coord,station_coord).km
        station_distances.append(station_distance)
    # use radius condition to derive list
    local_stations = list(station_list.station[(np.array(station_distances) < radius)])
    return local_stations

# function to create detection-threshold histogram
def gen_detect_hist(party):
    # import pacakges
    import numpy as np
    import matplotlib.pyplot as plt
    # Initialize lists for detection and threshold values
    detection_list = []
    threshold_list = []
    # loop through each family
    for i in range(len(party.families)):
        family = party[i]
        # loop through detections in each family
        for j in range(len(family)):
            # append detection & threshold values for every detection
            detection_value = abs(family[j].detect_val)
            detection_list.append(detection_value)
            threshold_value = family[j].threshold
            threshold_list.append(threshold_value)
    # plot distribution of detections as histogram and save
    detection_array = np.array(detection_list) - np.array(threshold_list)
    detection_floor = np.floor(min(detection_array))
    detection_ceil = np.ceil(max(detection_array))
    fig, ax = plt.subplots()
    ax.grid(True)
    ax.hist(detection_array, bins=np.arange(detection_floor, detection_ceil, 0.1), color='teal', edgecolor='black')
    ax.set_xlim([detection_floor, max(detection_array)])
    ax.set_xlabel('Detection-Threshold Gap')
    ax.set_ylabel('Frequency')
    ax.set_title('Histogram of Detection Values')
    fig.show()

def writer(outpath,object):
    # check for file extension
    if (outpath[-4:] != '.tgz') and (outpath[-4:] != '.xml'):
        raise ValueError('Invalid file extension, please check path input.')
    # import pacakges and object classes
    import os
    from obspy import Catalog
    from eqcorrscan import Party
    from eqcorrscan.core.match_filter.tribe import Tribe
    # if the file already exists, remove it
    if os.path.exists(outpath):
        os.remove(outpath)
    # initialize class samplers
    catalog_type = Catalog()
    party_type = Party()
    tribe_type = Tribe()
    # use if condition to choose saving style
    if type(object) == type(catalog_type):
        object.write(outpath, format='QUAKEML')
    elif type(object) == type(party_type):
        object.write(outpath, format='tar')
    elif type(object) == type(tribe_type):
        object.write(outpath)
    else:
        raise ValueError('Review writer(): input has invalid object type!')

def reader(inpath):
    # check for file extension
    if (inpath[-4:] != '.tgz') and (inpath[-4:] != '.xml'):
        raise ValueError('Invalid file extension, please check path input.')
    # import packages
    from obspy import read_events
    from eqcorrscan import Party
    from eqcorrscan.core.match_filter.tribe import read_tribe
    # use if condition to choose reading style
    if ('catalog' in inpath) and (inpath[-4:] == '.xml'):
        object = read_events(inpath, format='QUAKEML')
    elif ('party' in inpath) and (inpath[-4:] == '.tgz'):
        object = Party().read(inpath)
    elif ('tribe' in inpath) and (inpath[-4:] == '.tgz'):
        object = read_tribe(inpath)
    else:
        raise ValueError('Review reader(): input path does not follow naming convention.')
    return object

# read traces from local data files (NOTE: limited to less than 24h in duration)
def read_trace(data_dir,station,channel,starttime,endtime,tolerance=4e4):
    # import packages
    import glob
    from obspy import Stream, Trace, read
    from toolbox import remove_boxcars
    # initialize list of data filenames for reading
    data_filenames = []
    # extract year and day from starttime
    start_year = starttime.year
    start_julday = starttime.julday
    end_year = endtime.year
    end_julday = endtime.julday
    # craft file string and append
    data_filename = data_dir + station + '.' + channel + '.' + str(start_year) + ':' + f'{start_julday:03}' + ':*'
    data_filenames.append(data_filename)
    # add another day if starttime and endtime are on different data files
    if start_year != end_year or start_julday != end_julday:
        data_filename = data_dir + station + '.' + channel + '.' + str(end_year) + ':' + f'{end_julday:03}' + ':*'
        data_filenames.append(data_filename)
    # read in all traces we need
    tr_unmerged = Stream()
    for data_filename in data_filenames:
        matching_filenames = (glob.glob(data_filename))
        for matching_filename in matching_filenames:
            try:
                tr_contribution = read(matching_filename)
                tr_unmerged = tr_unmerged + tr_contribution
            except:
                continue
    # remove boxcar spikes, detrend, merge, and trim to desired starttime and endtime
    tr_merged = tr_unmerged.copy()
    tr_merged.trim(starttime=starttime, endtime=endtime)
    tr_merged = remove_boxcars(tr_merged, tolerance)
    tr_merged.detrend("simple").merge()
    # since a single station and channel is provided, the traces should merge into 1
    if len(tr_merged) == 0:
        output = Stream()
    elif len(tr_merged) == 1:
        output = tr_merged[0]
    else:
        raise ValueError('Error in read_trace(), function returns stream rather than trace.')
    return output

def prepare_catalog_stream(data_dir,catalog,resampling_frequency,tolerance):
    # import packages
    import glob
    from obspy import Stream, read
    from toolbox import remove_boxcars
    # get unique list of all station channel data needed for the day's events
    data_filenames = []
    for event in catalog:
        for pick in event.picks:
            # extract key information from pick
            sta = pick.waveform_id.station_code
            chan = pick.waveform_id.channel_code
            pick_year = pick.time.year
            pick_julday = pick.time.julday
            # craft file string and append
            data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year) + ':' + f'{pick_julday:03}' + ':*'
            data_filenames.append(data_filename)
            # add next day if pick occurs in the first 15 minutes of the day
            if pick.time.hour == 0 and pick.time.minute < 15:
                if pick.time.julday == 1:  # special case for first day of the year
                    data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year - 1) + ':365:*'
                    data_filenames.append(data_filename)
                else:
                    data_filename = data_dir + sta + '.' + chan + '.' + str(
                        pick_year) + ':' + f'{(pick_julday - 1):03}' + ':*'
                    data_filenames.append(data_filename)
            # add previous day if pick occurs in the last 15 minutes of the day
            if pick.time.hour == 23 and pick.time.minute > 45:
                if pick.time.julday == 365:  # special case for last day of the year
                    data_filename = data_dir + sta + '.' + chan + '.' + str(pick_year + 1) + ':001:*'
                    data_filenames.append(data_filename)
                else:
                    data_filename = data_dir + sta + '.' + chan + '.' + str(
                        pick_year) + ':' + f'{(pick_julday + 1):03}' + ':*'
                    data_filenames.append(data_filename)
    # now compile unique and sort
    data_filenames = list(set(data_filenames))
    data_filenames.sort()
    # read in all streams we need
    stream = Stream()
    for data_filename in data_filenames:
        matching_filenames = (glob.glob(data_filename))
        for matching_filename in matching_filenames:
            try:
                stream_contribution = read(matching_filename)
                stream = stream + stream_contribution
            except:
                continue
    # remove boxcar spikes, resample and merge
    stream = remove_boxcars(stream, tolerance)
    stream = stream.resample(resampling_frequency)
    stream = stream.merge()
    return stream

# function to use IRIS to download same streams as input
def client_download(stream_in,source='IRIS'):
    # import packages
    from obspy import Stream
    from waveform_collection import gather_waveforms
    # initialize stream_out
    stream_out = Stream()
    # loop through input stream
    for trace_in in stream_in:
        network = trace_in.stats.network
        station = trace_in.stats.station
        location = trace_in.stats.location
        channel = trace_in.stats.channel
        starttime = trace_in.stats.starttime
        endtime = trace_in.stats.endtime
        trace_out = gather_waveforms(source=source, network=network, station=station,
                      location=location, channel=channel, starttime=starttime,
                      endtime=endtime)
        stream_out = stream_out + trace_out
    return(stream_out)

# get detection stream by downloading necessary streams
def get_detection(detection,data_dir='/home/data/redoubt/',length=10,resampling_frequency=50,tolerance=4e4,lowcut=1,highcut=10,plot=False):
    # import packages
    import glob
    from obspy import Stream, read
    from toolbox import remove_boxcars
    # extract detection time and stachan combinations
    detect_time = detection.detect_time
    year = detect_time.year
    julday = detect_time.julday
    stachans = detection.chans
    # initialize data filename list
    data_filenames = []
    # loop through stachan combinations
    for sta, chan in stachans:
        # stitch data filename and append
        data_filename = data_dir + sta + '.' + chan + '.' + str(year) + ':' + f'{julday:03}' + ':*'
        data_filenames.append(data_filename)
    # read in all streams we need
    stream = Stream()
    for data_filename in data_filenames:
        matching_filenames = (glob.glob(data_filename))
        for matching_filename in matching_filenames:
            try:
                stream_contribution = read(matching_filename)
                stream = stream + stream_contribution
            except:
                continue
    # trim streams to what we need
    stream.trim(starttime=detect_time,endtime=detect_time+length)
    # resample all downloaded streams to prevent inconsistencies
    stream.resample(resampling_frequency)
    # remove boxcar spikes, detrend and merge
    stream = remove_boxcars(stream, tolerance)
    stream = stream.detrend("simple")
    # filter and taper
    stream = stream.filter('bandpass',freqmax=highcut,freqmin=lowcut)
    stream = stream.taper(0.05, type='hann', max_length=(0.75 * 1024 / 50))
    stream = stream.merge()
    if plot:
        stream.plot(color='b',equal_scale=False, size=(800, 600))
    return stream

# # plot detection-threshold gap by UTCDateTime
# def gen_detect_scatter(party,option):
#     # import packages
#     import numpy as np
#     import matplotlib.pyplot as plt
#     # initialize lists
#     detect_times = []
#     channel_numbers = []
#     detection_values = []
#     threshold_values = []
#     # loop through families
#     for family in party.families:
#         # loop through detections in the family
#         for detection in family:
#             # append essential information
#             detect_times.append(detection.detect_time)
#             channel_numbers.append(len(detection.chans))
#             detection_values.append(abs(detection.detect_val))
#             threshold_values.append(detection.threshold)
#     # calculate surplus and convert lists to arrays for plotting
#     surplus = np.array(detection_values) - np.array(threshold_values)
#     detect_times = np.array(detect_times)
#     channel_numbers = np.array(channel_numbers)
#     # plot scatters
#     fig, ax = plt.subplots(figsize=(8, 7))
#     ax.plot(detect_times,surplus,'.',color='teal')
#     ax.axhline(y=0,color='red')
#     ax.set_xticks(xticks)
#     ax.set_xticklabels(xticklabels,rotation=20,horizontalalignment='right')
#     ax.set_xlabel('UTCDate')
#     ax.set_ylabel('Detection-Threshold Gap')
#     ax.set_title('Detection Quality vs Time')
#     ax.grid()
#     fig.show()